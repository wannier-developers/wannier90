!-*- mode: F90; mode: font-lock -*-!

! =================================
!
! CHANGE INPUT VARIABLES AS FOLLOWS
!
! optics_task = ahc, morb, spectrum
! spectrum_type = sigma_ab, sigma_abc 
!                         (currently mcd, ord, ahe, orb, gyro, noa, mespn)
! time_parity = even, odd
! sigma_abc_onlyorb = F,T (default F)
!
! =================================

module w90_berry_wanint

  use w90_constants, only : dp

  implicit none

  private

  public :: berry,get_imf_ab_k,get_img_ab_k,get_imh_ab_k
  real(kind=dp), parameter :: eps=1.0e-7

  integer       :: nfreq,nfreq_cut,aa(10),bb(10),cc(10)
  real(kind=dp) :: d_freq,d_freq_cut
  logical       :: T_odd

  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine berry
  !=============================================================!
  !                                                             !
  ! Computes the following quantities:                          !
  !                                                             !
  !   (i) Intrinsic anomalous Hall conductivity                 ! 
  !  (ii) Interband optical conductivity in the electric-dipole ! 
  !       approximation (order q^0)                             !
  ! (iii) Joint density of states                               !
  ! (iii) Orbital magnetization                                 !
  !  (iv) Interband optical conductivity in the magnetic-dipole !
  !       electric-quadrupole approximation (order q^1)         !
  !                                                             !
  !=============================================================!

    use w90_constants, only     : dp,pi,cmplx_i,elem_charge_SI,hbar_SI,&
                                  elec_mass_SI,eV_au,bohr,bohr_magn_SI,&
                                  eps0_SI,speedlight_SI
    use w90_comms, only         : on_root,num_nodes,my_node_id,comms_reduce
    use w90_io, only            : io_error,stdout,io_file_unit,seedname,&
                                  io_stopwatch
    use w90_wanint_common, only : nrpts,irvec,num_int_kpts_on_node,int_kpts,&
                                  adkpt,weight,rpt_origin
    use w90_parameters, only    : timing_level,alpha,beta,gamma,num_wann,&
                                  optics_num_points,&
                                  optics_adaptive_pts,optics_adaptive_thresh,&
                                  wanint_kpoint_file,cell_volume,transl_inv,&
                                  optics_task,optics_min_energy,&
                                  optics_max_energy,optics_energy_step,&
                                  adpt_smr_steps,adpt_smr_width,band_by_band,&
                                  spn_decomp,found_fermi_energy,fermi_energy,&
                                  omega_from_FF,sigma_abc_onlyorb,&
                                  ecut_spectralsum
    use w90_get_oper, only      : get_HH_R,get_AA_R,get_BB_R,get_CC_R,get_FF_R,&
                                  get_SS_R


    real(kind=dp) :: imf_ab_k(3),imf_ab(3)
    real(kind=dp) :: img_ab_k(3),img_ab(3)
    real(kind=dp) :: imh_ab_k(3),imh_ab(3)
    real(kind=dp) :: ahc(3),LCtil(3),ICtil(3)

    ! Repeat for optical conductivity
    !
    ! 'sig_ab_k' is the contrib to the optical conductivity from one k-point,
    ! 'sig_ab_node' from k-points in one node, 'sig_ab' is the BZ average.
    ! 'ahc_kk' is the cumulative Kramers-Kronig transform of the MCD spectrum
    ! 'm_orb_sr' is the cumulative dichroic f-sum rule
    !
    !   sig_ab_(1,:,:,:) is the 'DD' term in the notation of YWSV07 
    !   sig_ab_(2,:,:,:) is the 'DA' term
    !   sig_ab_(3,:,:,:) is the 'AA' term

    real(kind=dp), allocatable :: sig_ab_k(:,:,:,:)
    real(kind=dp), allocatable :: sig_ab(:,:,:,:)
    real(kind=dp), allocatable :: ahc_kk(:,:,:)
    real(kind=dp), allocatable :: m_orb_sr(:,:)

    ! Optical conductivity divided by the frequency
    ! (used to compute the anomalous Hall conductivity via Kramers-Kronig)
    !
    real(kind=dp), allocatable :: sig_ab_over_freq_k(:,:,:,:)
    real(kind=dp), allocatable :: sig_ab_over_freq(:,:,:,:)
    
    ! Joint density of states
    !
    real(kind=dp), allocatable :: jdos_k(:,:,:)
    real(kind=dp), allocatable :: jdos(:,:,:)

    ! Spatially-dispersive optical conductivity
    ! sig_abc_(:,1) is the "matrix element term" part of the orbital contrib
    ! sig_abc_(:,2) is the "energy term" part of the orbital contribution
    ! sig_abc_(:,3) is the spin contribution ("matrix element term"-only)
    !
    real(kind=dp), allocatable :: sig_abc_k(:,:)
    real(kind=dp), allocatable :: sig_abc(:,:)

    ! Same as above, but contrib to the static value from optical transitions 
    ! below a given frequency, scanned between zero and ecut_spectralsum
    !
    real(kind=dp), allocatable :: sig_abc_cut_k(:,:)
    real(kind=dp), allocatable :: sig_abc_cut(:,:)

    ! Traceless optical magnetoelectric (ME) tensor alpha_ij
    ! and totally symmetric electric quadrupole (EQ) response tensor gamma_ijl 
    !
    real(kind=dp), allocatable :: alpha_k(:,:,:,:)
    real(kind=dp), allocatable :: alpha_me(:,:,:,:)
    real(kind=dp), allocatable :: imgamma_k(:,:,:)
    real(kind=dp), allocatable :: imgamma(:,:,:)
    
    ! Static spin-electronic magnetoelectric tensor alphaspn_ij, and 
    ! contribution from optical transitions below a given frequency, scanned
    ! between zero and ecut_spectralsum
    !
    real(kind=dp) :: alphaspn_k(3,3)
    real(kind=dp) :: alphaspn(3,3)
    !
    real(kind=dp), allocatable :: alphaspn_cut_k(:,:,:)
    real(kind=dp), allocatable :: alphaspn_cut(:,:,:)

    real(kind=dp)     :: kweight,kweight_adpt,kpt(3),freq,adpt_trigger
    integer           :: i,j,p,loop_x,loop_y,loop_z,loop_kpt,loop_adpt,&
                         loop_f,adpt_counter,ndim,&
                         ierr,jdos_unit,AA_unit,DA_unit,DD_unit,&
                         tot_unit,alpha_unit(3,3),gamma_unit(10)
    character(len=20) :: file_name
    logical           :: eval_ahe,eval_orb,eval_sig_ab,eval_sig_abc,eval_ME_EQ,eval_MEspn

    ! Units conversion factors
    !
    real(kind=dp) :: ahc_conv,orb_conv,sig_ab_conv,Morb_conv,fac

    if(.not.found_fermi_energy) call io_error&
         (&
    'Must set either "fermi_energy" or "num_elec_cell" for optical properties'&
         )
    
    ! Must initialize to .false. all eval_ flags
    !
    eval_ahe=.false.
    eval_orb=.false.
    eval_sig_ab=.false.
    eval_sig_abc=.false.
    eval_ME_EQ=.false.
    eval_MEspn=.false.
    T_odd=.false.
    if(index(optics_task,'ahe')>0) then
       eval_ahe=.true.
    elseif(index(optics_task,'orb')>0) then
       eval_orb=.true.
    end if
    if(index(optics_task,'mcd')>0) then
       T_odd=.true.
       eval_sig_ab=.true.
    elseif(index(optics_task,'ord')>0) then
       T_odd=.false.
       eval_sig_ab=.true.
    elseif(index(optics_task,'gyro')>0) then
       T_odd=.true.
       eval_sig_abc=.true.
       eval_ME_EQ=.true.
    elseif(index(optics_task,'noa')>0) then
       T_odd=.false.
       eval_sig_abc=.true.
    end if
    if(index(optics_task,'mespn')>0) then
       T_odd=.true.
       eval_MEspn=.true.
    endif
    if(alpha==0.or.beta==0) call io_error&
     ('Must specify cartesian directions alpha and beta for optical properties')
    if(gamma==0.and.eval_sig_abc) call io_error&
         ('Must specify cartesian direction gamma for optics_task=gyro,noa')

    ! Set up required Wannier matrix elements
    !
    call get_HH_R
    call get_AA_R
    if(eval_orb.or.eval_sig_abc.or.eval_ME_EQ) then
       call get_BB_R
       call get_CC_R
    endif
    if(omega_from_FF.or.eval_sig_abc.or.eval_ME_EQ) call get_FF_R
    if((eval_sig_ab.and.spn_decomp).or.&
       (eval_sig_abc .and. .not.(sigma_abc_onlyorb)) .or.&
       (eval_ME_EQ .and. .not.(sigma_abc_onlyorb)) .or. eval_MEspn) call get_SS_R

    ! Start timing this routine after setting up the real-space matrix elements
    if (timing_level>1.and.on_root) call io_stopwatch('berry_wanint: berry',1)

    if(eval_sig_ab.or.eval_sig_abc.or.eval_ME_EQ) then
       nfreq=nint((optics_max_energy-optics_min_energy)/optics_energy_step)
       d_freq=(optics_max_energy-optics_min_energy)/(nfreq-1)
    endif
    if(eval_sig_ab) then
       if(spn_decomp) then
          !
          ! The extra three entries contain the decomposition
          ! into up-->up, down-->down, and spin-flip transitions
          !
          ndim=4
       else
          ndim=1
       end if
       allocate(ahc_kk(3,adpt_smr_steps,ndim))
       allocate(m_orb_sr(3,adpt_smr_steps))
       allocate(sig_ab_k(3,nfreq,adpt_smr_steps,ndim))
       allocate(sig_ab(3,nfreq,adpt_smr_steps,ndim))
       allocate(sig_ab_over_freq_k(3,nfreq,adpt_smr_steps,ndim))
       allocate(sig_ab_over_freq(3,nfreq,adpt_smr_steps,ndim))
       allocate(jdos_k(nfreq,adpt_smr_steps,ndim))
       allocate(jdos(nfreq,adpt_smr_steps,ndim))
    end if
    if(eval_sig_abc) then
       allocate(sig_abc_k(nfreq,3))
       allocate(sig_abc(nfreq,3))
       nfreq_cut=nint(ecut_spectralsum/optics_energy_step)
       if(nfreq_cut<2) call io_error('Increase value of ecut_spectralsum')
       d_freq_cut=ecut_spectralsum/(nfreq_cut-1)
       allocate(sig_abc_cut_k(nfreq_cut,3))
       allocate(sig_abc_cut(nfreq_cut,3))
    endif
    if(eval_ME_EQ) then
       allocate(alpha_k(3,nfreq,3,3))
       allocate(alpha_me(3,nfreq,3,3))
       !
       ! N.B.: There should be no spin contribution to gamma, so ultimately want
       !       to use the dimensions (2,nfreq,10). For the moment include spin
       !       to check whether its contribution vanishs numerically
       !
       allocate(imgamma_k(3,nfreq,10))
       allocate(imgamma(3,nfreq,10))
       !
       ! Convention for totally symmetric gamma_abc --> gamma_p
       !
       ! a  b  c  |  p
       ! -------------
       ! 1  1  1  |  1
       ! 2  2  2  |  2
       ! 3  3  3  |  3
       ! 1  1  2  |  4
       ! 1  1  3  |  5
       ! 2  2  1  |  6
       ! 3  3  1  |  7
       ! 3  3  2  |  8
       ! 2  2  3  |  9
       ! 1  2  3  | 10
       !
       aa(1)= 1; bb(1)= 1; cc(1)= 1
       aa(2)= 2; bb(2)= 2; cc(2)= 2
       aa(3)= 3; bb(3)= 3; cc(3)= 3
       aa(4)= 1; bb(4)= 1; cc(4)= 2
       aa(5)= 1; bb(5)= 1; cc(5)= 3
       aa(6)= 2; bb(6)= 2; cc(6)= 1
       aa(7)= 3; bb(7)= 3; cc(7)= 1
       aa(8)= 3; bb(8)= 3; cc(8)= 2
       aa(9)= 2; bb(9)= 2; cc(9)= 3
       aa(10)=1; bb(10)=2; cc(10)=3
    endif
    if(eval_MEspn) then
       allocate(alphaspn_cut_k(nfreq_cut,3,3))
       allocate(alphaspn_cut(nfreq_cut,3,3))
    endif

    if(on_root) then

       write(stdout,'(/,/,1x,a)') '============================='
       write(stdout,'(1x,a)')     'Optical/transport properties:'
       write(stdout,'(1x,a)')     '============================='

       if(eval_ahe) write(stdout,'(/,3x,a)') '* Re[sigma_{A,'//&
            achar(119+alpha)//achar(119+beta)//&
      '}(0)]: dc anomalous Hall (antisymm) conductivity'

       if(eval_orb) write(stdout,'(/,3x,a)') '* M^orb_{'//&
            achar(119+alpha)//achar(119+beta)//&
      '}: Orbital magnetization'

       if(eval_sig_ab) then
          write(stdout,'(/,3x,a)') '* Joint density of states '//&
             '(file '//trim(trim(seedname))//'-jdos)'
          if(T_odd) then
             write(stdout,'(/,3x,a)') '* Im[sigma_{A,'//&
                  achar(119+alpha)//achar(119+beta)//&
                  '}(omega)]: Dichroic (antisymm) interband'
             write(stdout,'(5x,a)')&
                  'absorptive conductivity in units of S/cm '//&
             '(file '//trim(trim(seedname))//'-sig_ab^A)'
             write(stdout,'(/,3x,a)')&
                  '* Cumulative anomalous Hall conductivity from Kramers-Kronig'
             write(stdout,'(5x,a)') 'in units of S/cm '//&
             '(file '//trim(trim(seedname))//'-ahc_kk)'
             write(stdout,'(/,3x,a)') '* Dichroic f-sum rule'
             write(stdout,'(5x,a)') 'in units of (Bohr magn)/cell'//&
             '(file '//trim(trim(seedname))//'-m_orb_sr)' 
          else
             write(stdout,'(/,3x,a)') '* Re[sigma_{S,'//&
                  achar(119+alpha)//achar(119+beta)//&
                  '}(omega)]: Ordinary (symm) interband' 
             write(stdout,'(5x,a)')&
                  'absorptive conductivity in units of S/cm '//&
             '(file '//trim(trim(seedname))//'-sig_ab^S)'
          end if
          write(stdout,'(/,1x,a,6(f6.3,1x))')&
            'Adaptive smearing width prefactors for optical properties: ',&
            (adpt_smr_width(i),i=1,adpt_smr_steps)
      
          write(stdout,'(/,1x,a,f8.4)')&
               'Frequency step for optical properties (eV): ',d_freq
       end if

       if(eval_sig_abc) then
          if(T_odd) then
             write(stdout,'(/,3x,a)')&
                  '* Interband spatially dispersive conductance Im[sigma_{S,'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//&
                  '}(omega)]'
             write(stdout,'(5x,a)')
             write(stdout,'(5x,a)')&
                  'Time-odd (symm) reactive part in e^2/h units '
              write(stdout,'(5x,a)') 'File: '//&
                   trim(trim(seedname))//'-sig^S_'//&
             achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat'
          else
             write(stdout,'(/,3x,a)')&
                  '* Optical rotatory power rho divided by (hbar.omega)^2'
             write(stdout,'(5x,a)') 'Units: deg/[mm.(eV)^2]'
             write(stdout,'(5x,a)') 'File: '//trim(trim(seedname))//&
                  '-rhobw2_'//&
             achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat'
             write(stdout,'(/,5x,a)') 'rho/(hbar omega)^2=Re[sigma_{A,'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//&
                  '}(omega)]/(hbar omega)/(2.hbar.c^2.eps_0)'
          endif
       endif

       if(eval_ME_EQ) write(stdout,'(/,3x,a)')&
            '* Optical magnetoelectric and electric-quadrupole responses'

       if(eval_MEspn) write(stdout,'(/,3x,a)')&
            '* Static spin-electronic magnetoelectric tensor'

       write(stdout,'(/,1x,a,f10.6,1x,a)') 'Fermi level: ',&
            fermi_energy,'eV'
       if(transl_inv) then
          if(eval_orb.or.eval_sig_abc) then
             write(stdout,*) 'transl_inv=T disabled for Morb'
             stop
          endif
          write(stdout,'(/,1x,a)')&
               'Using a translationally-invariant discretization for the'
          write(stdout,'(1x,a)')&
               'band-diagonal Wannier matrix elements of r, etc.'
!       else
!          write(stdout,'(/,1x,a)')&
!               'Using the discretization in Appendix B of MV97 for the'
!          write(stdout,'(1x,a)')&
!               'Wannier matrix elements containing the position operator'
       endif

    end if !on_root

    ! ---------------------------------------
    ! TO DO: Calculate here the adaptive grid
    ! ---------------------------------------
       
    ! Must initialize to zero all arrays
    !
    adpt_counter=0
    imf_ab=0.0_dp
    img_ab=0.0_dp
    imh_ab=0.0_dp
    sig_ab=0.0_dp
    sig_ab_over_freq=0.0_dp
    jdos=0.0_dp
    sig_abc=0.0_dp
    sig_abc_cut=0.0_dp
    alpha_me=0.0_dp
    imgamma=0.0_dp
    alphaspn=0.0_dp
    alphaspn_cut=0.0_dp

    if(wanint_kpoint_file) then
       
       if(on_root) then
          !
          ! As implemented, this does not always work. For example,
          ! it gives the wrong ordinary absorption sigma_{S,xx}(omega) for bcc
          ! iron with magnetization along z (but it does give the correct
          ! sigma_{S,zz}(omega), as well as the dichroic sigma_{A,xy}. The 
          ! correct general implementation would be as outlines in, e.g., the 
          ! 3rd paragraph of Sec. III.C of PRB 60, 14105 (1999).
          !
          write(stdout,'(/,1x,a)') 'Sampling the irreducible BZ only'
          write(stdout,'(3x,a)')&
          'WARNING: - IBZ implementation is currently limited to simple cases:'
          write(stdout,'(3x,a)')&
           '           Check results agains a full BZ calculation!'
       end if

       ! Loop over k-points on the irreducible wedge of the Brillouin zone,
       ! read from file 'kpoint.dat'
       !
       do loop_kpt=1,num_int_kpts_on_node(my_node_id)
          kpt(:)=int_kpts(:,loop_kpt)
          kweight=weight(loop_kpt)
          kweight_adpt=kweight/optics_adaptive_pts**3
          if(eval_ahe) then 
             call get_imf_ab_k(kpt,imf_ab_k)
            ! 
             ! Use the same adaptive grid for the entire set of Fermi levels.
             ! Apply the threshold criterion to the *true* scf Fermi level
             !
             adpt_trigger=abs(sum(imf_ab_k))
          elseif(eval_orb) then
             call get_imf_ab_k(kpt,imf_ab_k)
             call get_img_ab_k(kpt,img_ab_k)
             call get_imh_ab_k(kpt,imh_ab_k)
             adpt_trigger=abs(sum(img_ab_k)+sum(imh_ab_k)&
                  -2.0_dp*fermi_energy*sum(imf_ab_k))
          else
             !
             ! Ensure that adaptive refinement is not triggered
             !
             adpt_trigger=optics_adaptive_thresh-1.0_dp
          end if
          if(eval_sig_ab)&
               call get_sig_ab_k(kpt,sig_ab_k,sig_ab_over_freq_k,jdos_k)
          if(eval_sig_abc) call get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
          if(eval_ME_EQ) call get_alpha_k(kpt,alpha_k,imgamma_k)
          if(eval_MEspn) call get_alphaspn_k(kpt,alphaspn_k,alphaspn_cut_k)
          !
          ! Decide whether or not to trigger adaptive refinement of 
          ! integration k-mesh, when computing AHC (and possibly MCD)
          !
          if(adpt_trigger>optics_adaptive_thresh) then
             adpt_counter=adpt_counter+1
             do loop_adpt=1,optics_adaptive_pts**3
                if(eval_ahe) then
                   call get_imf_ab_k(kpt(:)+adkpt(:,loop_adpt),imf_ab_k)
                   imf_ab=imf_ab+imf_ab_k*kweight_adpt
                elseif(eval_orb) then
                   call get_imf_ab_k(kpt(:)+adkpt(:,loop_adpt),imf_ab_k)
                   imf_ab=imf_ab+imf_ab_k*kweight_adpt
                   call get_img_ab_k(kpt(:)+adkpt(:,loop_adpt),img_ab_k)
                   img_ab=img_ab+img_ab_k*kweight_adpt
                   call get_imh_ab_k(kpt(:)+adkpt(:,loop_adpt),imh_ab_k)
                   imh_ab=imh_ab+imh_ab_k*kweight_adpt
                end if
                if(eval_sig_ab) then
                   call get_sig_ab_k(kpt(:)+adkpt(:,loop_adpt),&
                        sig_ab_k,sig_ab_over_freq_k,jdos_k)
                   sig_ab=sig_ab+sig_ab_k*kweight_adpt
                   sig_ab_over_freq=sig_ab_over_freq&
                        +sig_ab_over_freq_k*kweight_adpt
                   jdos=jdos+jdos_k*kweight_adpt
                end if
                if(eval_sig_abc) then
                   call get_sig_abc_k(kpt(:)+adkpt(:,loop_adpt),sig_abc_k,&
                        sig_abc_cut_k)
                   sig_abc=sig_abc+sig_abc_k*kweight_adpt
                   sig_abc_cut=sig_abc_cut&
                                     +sig_abc_cut_k*kweight_adpt
                end if
                if(eval_ME_EQ) then
                   call get_alpha_k(kpt(:)+adkpt(:,loop_adpt),alpha_k,imgamma_k)
                   alpha_me=alpha_me+alpha_k*kweight_adpt
                   imgamma=imgamma+imgamma_k*kweight_adpt
                end if
                if(eval_MEspn) then
                   call get_alphaspn_k(kpt(:)&
                        +adkpt(:,loop_adpt),alphaspn_k,alphaspn_cut_k)
                   alphaspn=alphaspn+alphaspn_k*kweight_adpt
                   alphaspn_cut=alphaspn_cut&
                        +alphaspn_cut_k*kweight_adpt
                end if
             end do
          else
             if(eval_ahe) then
                imf_ab=imf_ab+imf_ab_k*kweight
             elseif(eval_orb) then 
                imf_ab=imf_ab+imf_ab_k*kweight
                img_ab=img_ab+img_ab_k*kweight
                imh_ab=imh_ab+imh_ab_k*kweight
             end if
             if(eval_sig_ab) then
                sig_ab=sig_ab+sig_ab_k*kweight
                sig_ab_over_freq=sig_ab_over_freq&
                     +sig_ab_over_freq_k*kweight
                jdos=jdos+jdos_k*kweight
             end if
             if(eval_sig_abc) then
                sig_abc=sig_abc+sig_abc_k*kweight
                sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight
             endif
             if(eval_ME_EQ) then
                alpha_me=alpha_me+alpha_k*kweight
                imgamma=imgamma+imgamma_k*kweight
             endif
             if(eval_MEspn) then
                alphaspn=alphaspn+alphaspn_k*kweight
                alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight
             endif
          end if
       end do

    else
       !
       ! Do not read 'kpoint.dat'. Loop over a uniform grid in the full BZ
       !
       if (on_root) write(stdout,'(/,1x,a)') 'Sampling the full BZ'
       
       kweight=1.0_dp/optics_num_points**3
       kweight_adpt=kweight/optics_adaptive_pts**3
       do loop_kpt=my_node_id,optics_num_points**3-1,num_nodes
          loop_x=loop_kpt/optics_num_points**2
          loop_y=(loop_kpt-loop_x*optics_num_points**2)/optics_num_points
          loop_z=loop_kpt-loop_x*optics_num_points**2-loop_y*optics_num_points
          kpt(1)=real(loop_x,dp)/optics_num_points
          kpt(2)=real(loop_y,dp)/optics_num_points
          kpt(3)=real(loop_z,dp)/optics_num_points
          
          if(eval_ahe) then
             call get_imf_ab_k(kpt,imf_ab_k)
             adpt_trigger=abs(sum(imf_ab_k))
          elseif(eval_orb) then
             call get_imf_ab_k(kpt,imf_ab_k)
             call get_img_ab_k(kpt,img_ab_k)
             call get_imh_ab_k(kpt,imh_ab_k)
             adpt_trigger=abs(sum(img_ab_k)&
                  +sum(imh_ab_k)-2.0_dp*fermi_energy*sum(imf_ab_k))
          else
             !
             ! Ensure that adaptive refinement is not triggered
             !
             adpt_trigger=optics_adaptive_thresh-1.0_dp
          end if
          if(eval_sig_ab)&
               call get_sig_ab_k(kpt,sig_ab_k,sig_ab_over_freq_k,jdos_k)
          if(eval_sig_abc) call get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
          if(eval_ME_EQ) call get_alpha_k(kpt,alpha_k,imgamma_k)
          if(eval_MEspn) call get_alphaspn_k(kpt,alphaspn_k,alphaspn_cut_k)

          if(adpt_trigger>optics_adaptive_thresh) then
             adpt_counter=adpt_counter+1
             do loop_adpt=1,optics_adaptive_pts**3
                if(eval_ahe) then
                   call get_imf_ab_k(kpt(:)+adkpt(:,loop_adpt),imf_ab_k)
                   imf_ab=imf_ab+imf_ab_k*kweight_adpt
                elseif(eval_orb) then
                   call get_imf_ab_k(kpt(:)+adkpt(:,loop_adpt),imf_ab_k)
                   imf_ab=imf_ab+imf_ab_k*kweight_adpt
                   call get_img_ab_k(kpt(:)+adkpt(:,loop_adpt),img_ab_k)
                   img_ab=img_ab+img_ab_k*kweight_adpt
                   call get_imh_ab_k(kpt(:)+adkpt(:,loop_adpt),imh_ab_k)
                   imh_ab=imh_ab+imh_ab_k*kweight_adpt
                end if
                if(eval_sig_ab) then
                   call get_sig_ab_k(kpt(:)+adkpt(:,loop_adpt),&
                        sig_ab_k,sig_ab_over_freq_k,jdos_k)
                   sig_ab=sig_ab+sig_ab_k*kweight_adpt
                   sig_ab_over_freq=sig_ab_over_freq&
                        +sig_ab_over_freq_k*kweight_adpt
                   jdos=jdos+jdos_k*kweight_adpt
                end if
                if(eval_sig_abc) then
                   call get_sig_abc_k(kpt(:)+adkpt(:,loop_adpt),sig_abc_k,&
                        sig_abc_cut_k)
                   sig_abc=sig_abc+sig_abc_k*kweight_adpt
                   sig_abc_cut=sig_abc_cut&
                                     +sig_abc_cut_k*kweight_adpt
                end if
                if(eval_ME_EQ) then
                   call get_alpha_k(kpt(:)+adkpt(:,loop_adpt),alpha_k,imgamma_k)
                   alpha_me=alpha_me+alpha_k*kweight_adpt
                   imgamma=imgamma+imgamma_k*kweight_adpt
                end if
                if(eval_MEspn) then
                   call get_alphaspn_k(kpt(:)+adkpt(:,loop_adpt),&
                        alphaspn_k,alphaspn_cut_k)
                   alphaspn=alphaspn+alphaspn_k*kweight_adpt
                   alphaspn_cut=alphaspn_cut&
                        +alphaspn_cut_k*kweight_adpt
                end if
             end do
          else
             if(eval_ahe) then
                imf_ab=imf_ab+imf_ab_k*kweight
             elseif(eval_orb) then 
                imf_ab=imf_ab+imf_ab_k*kweight
                img_ab=img_ab+img_ab_k*kweight
                imh_ab=imh_ab+imh_ab_k*kweight
             end if
             if(eval_sig_ab) then
                sig_ab=sig_ab+sig_ab_k*kweight
                sig_ab_over_freq=sig_ab_over_freq&
                     +sig_ab_over_freq_k*kweight
                jdos=jdos+jdos_k*kweight
             end if
             if(eval_sig_abc) then
                sig_abc=sig_abc+sig_abc_k*kweight
                sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight
             endif
             if(eval_ME_EQ) then
                alpha_me=alpha_me+alpha_k*kweight
                imgamma=imgamma+imgamma_k*kweight
             endif
             if(eval_MEspn) then
                alphaspn=alphaspn+alphaspn_k*kweight
                alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight
             endif
          end if
       end do !loop_kpt
       
    end if !wanint_kpoint_file

! Collect contributions from all nodes    
!
    if(eval_ahe) then
       call comms_reduce(imf_ab(1),3,'SUM')
       call comms_reduce(adpt_counter,1,'SUM')
    elseif(eval_orb) then
       call comms_reduce(imf_ab(1),3,'SUM')
       call comms_reduce(img_ab(1),3,'SUM')
       call comms_reduce(imh_ab(1),3,'SUM')
       call comms_reduce(adpt_counter,1,'SUM')
    end if
    if(eval_sig_ab) then
       call comms_reduce(sig_ab(1,1,1,1),3*nfreq*adpt_smr_steps*ndim,'SUM')
       call comms_reduce(sig_ab_over_freq(1,1,1,1),3*nfreq*adpt_smr_steps*ndim,'SUM')
       call comms_reduce(jdos(1,1,1),nfreq*adpt_smr_steps*ndim,'SUM')
    endif
    if(eval_sig_abc) then
       call comms_reduce(sig_abc(1,1),3*nfreq,'SUM')
       call comms_reduce(sig_abc_cut(1,1),3*nfreq_cut,'SUM')
    endif
    if(eval_ME_EQ) then
       call comms_reduce(alpha_me(1,1,1,1),3*nfreq*3*3,'SUM')
       !
       ! *****************************************
       ! Eventually change 3*nfreq*3 --> 2*nfreq*3
       ! *****************************************
       !
       call comms_reduce(imgamma(1,1,1),3*nfreq*10,'SUM')
    endif
    if(eval_MEspn) then
!jry corrected alphasp_node to alphaspn_node !! check!!
       call comms_reduce(alphaspn(1,1),3*3,'SUM')
       call comms_reduce(alphaspn_cut(1,1,1),nfreq_cut*3*3,'SUM')
    endif
    
    if(on_root) then

       if(eval_ahe) then

          ! --------------------------------------------------------------------
          ! At this point ms_imf_ab contains 
          !
          ! (1/N) sum_k Omega_{alpha beta}(k), 
          !
          ! an approximation to 
          !
          ! V_c.int dk/(2.pi)^3 Omega_{alpha beta}(k) dk 
          !
          ! (V_c is the cell volume). We want
          !
          ! sigma_{alpha beta}=-(e^2/hbar) int dk/(2.pi)^3 Omega(k) dk 
          !
          ! Hence need to multiply by -(e^2/hbar.V_c). 
          ! To get a conductivity in units of (Ohm. cm)^{-1},
          !
          ! (i)   Divide by V_c to obtain (1/N) sum_k omega(k)/V_c, which has 
          !       units of [L]^{-1} (Berry curvature Omega(k) has units of 
          !       [L]^2)
          ! (ii)  I am assuming the working unit of length is Angstroms. 
          !       Multiply by 10^8 to convert above quantity to (cm)^{-1}
          ! (iii) Multiply by -e^2/hbar expressed in SI, which has units of 
          !       conductance, (Ohm)^{-1}, or Siemens (S), to get the final 
          !       result in S/cm
          !
          ! ================================
          ! ahc_conv = -e^2/(hbar.V_c*10^-8) 
          ! ================================
          !
          ! with 'V_c' in Angstroms^3, and 'e' and 'hbar' in SI units.
          ! --------------------------------------------------------------------
          !
          ahc_conv=-1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
          ahc(:)=imf_ab(:)*ahc_conv

          write(stdout,'(/,/,1x,a)') 'Integrated AHC (S/cm)'
          write(stdout,'(1x,a)')     '=============='
          write(stdout,'(/,1x,a)') 'Using new trace formulation'
          write(stdout,'(1x,a,f14.5,2x,f10.5,a)') 'Omega term :',&
               ahc(1),ahc(1)*100.0_dp/sum(ahc),' %'
          write(stdout,'(1x,a,f14.5,2x,f10.5,a)') 'J.AA term  :',&
               ahc(2),ahc(2)*100.0_dp/sum(ahc),' %'
          write(stdout,'(1x,a,f14.5,2x,f10.5,a)') 'J.J term   :',&
               ahc(3),ahc(3)*100.0_dp/sum(ahc),' %'
          write(stdout,'(1x,a)')&
               '--------------------------------------------'
          write(stdout,'(1x,a,f14.5,2x)')         'Total AHC  :',sum(ahc)
          
       elseif(eval_orb) then

          ! --------------------------------------------------------------------
          ! At this point X=img_ab(:)-fermi_energy*imf_ab(:) and
          !               Y=imh_ab(:)-fermi_energy*imf_ab(:) 
          ! contain, eg,
          !
          ! (1/N) sum_k X(k), where X(k)=-2*Im[g(k)-E_F.f(k)] 
          !
          ! This is an approximation to 
          !
          ! V_c.int dk/(2.pi)^3 X(k) dk 
          !
          ! (V_c is the cell volume). We want a magnetic moment per cell,
          ! in units of the Bohr magneton. The magnetization-like quantity is
          !
          ! \tilde{M}^LC=-(e/2.hbar) int dk/(2.pi)^3 X(k) dk 
          !
          ! So we take X and
          !
          !  (i)  The summand is an energy in eV times a Berry curvature in
          !       Ang^2. To convert to a.u., divide by 27.2 and by 0.529^2
          !  (ii) Multiply by -(e/2.hbar)=-1/2 in atomic units
          ! (iii) At this point we have a magnetic moment (per cell) in atomic 
          !       units. 1 Bohr magneton = 1/2 atomic unit, so need to multiply
          !       by 2 to convert it to Bohr magnetons
          ! --------------------------------------------------------------------
          !
          orb_conv=-eV_au/bohr**2
          LCtil(:)=(img_ab(:)-fermi_energy*imf_ab(:))*orb_conv
          ICtil(:)=(imh_ab(:)-fermi_energy*imf_ab(:))*orb_conv

          write(stdout,'(/,/,1x,a)')        '============================'
          write(stdout,'(1x,a,1x,e13.6,a)') 'Orbital magnetization M_orb:',&
               sum(LCtil)+sum(ICtil),   ' bohr magn/cell'
          write(stdout,'(1x,a)')            '============================'
          write(stdout,'(/,1x,a)')&
               'M_orb = LCtil + ICtil'
          write(stdout,'(/,29x,a4,12x,a4,12x,a4)') '(J0)','(J1)','(J2)'
          write(stdout,'(1x,a,1x,e13.6,3(1x,a,1x,e13.6))')&
               'LCtil=',sum(LCtil),'=',LCtil(1),'+',LCtil(2),'+',LCtil(3)
          write(stdout,'(1x,a,1x,e13.6,3(1x,a,1x,e13.6))')&
               'ICtil=',sum(ICtil),'=',ICtil(1),'+',ICtil(2),'+',ICtil(3)

       end if ! eval_ahe or eval_orb

       ! ---------!
       ! sigma_ab !
       ! ---------!
       !
       if(eval_sig_ab) then

          write(stdout,'(/,1x,a)')&
               '--------------------------------------------------------------'
          write(stdout,'(1x,a)')&
               'Output data files related to long-wavelength optical spectrum:'
          write(stdout,'(1x,a)')&
               '--------------------------------------------------------------'

          file_name=trim(seedname)//'_jdos.dat'
          write(stdout,'(/,3x,a)') '* '//file_name
          jdos_unit=io_file_unit()
          open(jdos_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
          do loop_f=1,nfreq
             freq=optics_min_energy+loop_f*d_freq
             write(jdos_unit,'(20E16.8)') freq,jdos(loop_f,:,:)
          enddo
          close(jdos_unit)

          ! --------------------------------------------------------------------
          ! At this point sig_ab contains 
          !
          ! (1/N) sum_k Optical_matrix_elem_{alpha beta}(k), 
          !
          ! an approximation to 
          !
          ! V_c.int dk/(2.pi)^3 Optical_matrix_elem_{alpha beta}(k) dk 
          !
          ! (V_c is the cell volume). We want the absorptive optical 
          ! conductivity
          !
          ! sigma=(pi.e^2/hbar)int dk/(2.pi)^3 Optical_matrix_elem(k)dk,
          !  
          ! where the optical matrix element is (for MCD; replace "Im" by "Re" 
          ! for ordinary absorption)
          !
          ! sum_n^occ sum_m^empty Im(A_{nm,alpha}A_{mn,beta}).omega_{mn}.
          !                       .delta(omega-omega_{mn}),
          !
          ! which has the same units as Im(...), i.e., [L]^2
          ! (Berry connection has units of length, and delta(...) has units of
          ! inverse frequency)
          !
          ! Hence need to multiply by pi.e^2/(hbar.V_c). To get a conductivity,
          !
          ! (i)   Divide by V_c to obtain 
          !                             (1/N) sum_k Optical_matrix_elem(k)/V_c, 
          !       which has units of 1/[L] 
          ! (ii)  I am assuming the working unit of length is Angstroms. 
          !       Multiply by 10^8 to convert to centimeters
          ! (iii) Multiply by pi.e^2/hbar expressed in SI (e^2/hbar is the 
          !       quantum of conductance, which has SI units of Siemens) to 
          !       get the final result for the conductivity in units of 
          !       Siemens/cm
          !
          ! ========================================================
          ! sig_ab_conv = pi*10^8*elem_charge_SI^2/(hbar_SI*V_c_ang) 
          ! ========================================================
          !
          ! --------------------------------------------------------------------
          !
          sig_ab_conv=pi*1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
          sig_ab=sig_ab*sig_ab_conv
          
!          if(T_odd) then
!             file_name=trim(trim(seedname)//'-sig_ab^A_'//&
!                  achar(119+alpha)//achar(119+beta)//'_DD.dat')
!          else
!             file_name= trim(trim(seedname)//'-sig_ab^S_'//&
!                  achar(119+alpha)//achar(119+beta)//'_DD.dat')
!          end if
!          write(stdout,'(/,3x,a)') '* '//file_name
!          DD_unit=io_file_unit()
!          open(DD_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')

!          if(T_odd) then
!             file_name= trim(trim(seedname)//'-sig_ab^A_'//&
!                  achar(119+alpha)//achar(119+beta)//'_DA.dat')
!          else
!             file_name= trim(trim(seedname)//'-sig_ab^S_'//&
!                  achar(119+alpha)//achar(119+beta)//'_DA.dat')
!          end if
!          write(stdout,'(/,3x,a)') '* '//file_name
!          DA_unit=io_file_unit()
!          open(DA_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
          
!          if(T_odd) then
!             file_name= trim(trim(seedname)//'-sig_ab^A_'//&
!                  achar(119+alpha)//achar(119+beta)//'_AA.dat')
!          else
!             file_name= trim(trim(seedname)//'-sig_ab^S_'//&
!                  achar(119+alpha)//achar(119+beta)//'_AA.dat')
!          end if
!          write(stdout,'(/,3x,a)') '* '//file_name
!          AA_unit=io_file_unit()
!          open(AA_unit,FILE=file_name,FORM='FORMATTED')
          
          if(T_odd) then
             file_name= trim(trim(seedname)//'-sigA_'//&
                  achar(119+alpha)//achar(119+beta)//'.dat')
          else
             file_name= trim(trim(seedname)//'-sigS_'//&
                  achar(119+alpha)//achar(119+beta)//'.dat')
          end if
          write(stdout,'(/,3x,a)') '* '//file_name
          tot_unit=io_file_unit()
          open(tot_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   

          do loop_f=1,nfreq
             freq=optics_min_energy+loop_f*d_freq
!             write(DD_unit,'(20E16.8)') freq,sig_ab(1,loop_f,:,:)
!             write(DA_unit,'(20E16.8)') freq,sig_ab(2,loop_f,:,:)
!             write(AA_unit,'(20E16.8)') freq,sig_ab(3,loop_f,:,:)
             write(tot_unit,'(20E16.8)') freq,(sig_ab(1,loop_f,:,:)&
                                              +sig_ab(2,loop_f,:,:)&
                                              +sig_ab(3,loop_f,:,:))
             
          enddo
          
 !         close(DD_unit)
 !         close(DA_unit)
 !         close(AA_unit)
          close(tot_unit)
          
          if(T_odd) then
             
             ! Anomalous Hall conductivity as cumulative Kramers-Kronig 
             ! transform of the magnetic circular dichroism spectrum: 
             ! Fig. 5b & Eq. (43) of YWVS07 
          
             file_name=trim( trim(seedname)//'_ahc_kk_DD.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             DD_unit=io_file_unit()
             open(DD_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
                             
             file_name=trim(trim(seedname)//'_ahc_kk_DA.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             DA_unit=io_file_unit()
             open(DA_unit,FILE=file_name,STATUS='UNKNOWN',&
                  FORM='FORMATTED')
             
             file_name=trim(trim(seedname)//'_ahc_kk_AA.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             AA_unit=io_file_unit()
             open(AA_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
             
             file_name=trim( trim(seedname)//'_ahc_kk_tot.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             tot_unit=io_file_unit()
             open(tot_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   

             ahc_kk=0.0_dp
             do loop_f=nfreq,1,-1
                ahc_kk(1,:,:)=ahc_kk(1,:,:)&
                     +(2.0_dp/pi)*sig_ab_over_freq(1,loop_f,:,:)*d_freq
                ahc_kk(2,:,:)=ahc_kk(2,:,:)&
                     +(2.0_dp/pi)*sig_ab_over_freq(2,loop_f,:,:)*d_freq
                ahc_kk(3,:,:)=ahc_kk(3,:,:)&
                     +(2.0_dp/pi)*sig_ab_over_freq(3,loop_f,:,:)*d_freq
                freq=optics_min_energy+loop_f*d_freq
                !
                ! Since the factor (1/freq)*d_freq in Eq.(43) YWVS07 is 
                ! dimensionless, the conversion factor to obtain the AHC in 
                ! S/cm is the same as above for sig_ab
                !
                write(DD_unit,'(20E16.8)') freq,ahc_kk(1,:,:)*sig_ab_conv
                write(DA_unit,'(20E16.8)') freq,ahc_kk(2,:,:)*sig_ab_conv
                write(AA_unit,'(20E16.8)') freq,ahc_kk(3,:,:)*sig_ab_conv
                write(tot_unit,'(20E16.8)')freq,sig_ab_conv*(ahc_kk(1,:,:)&
                                                               +ahc_kk(2,:,:)&
                                                               +ahc_kk(3,:,:))
             end do
             
             close(DD_unit)
             close(DA_unit)
             close(AA_unit)
             close(tot_unit)
             
             ! Dichroic f-sum rule
             !
             file_name=trim( trim(seedname)//'_m_orb_sr_DD.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             DD_unit=io_file_unit()
             open(DD_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
             
             file_name=trim(trim(seedname)//'_m_orb_sr_DA.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             DA_unit=io_file_unit()
             open(DA_unit,FILE=file_name,STATUS='UNKNOWN',&
                  FORM='FORMATTED')
             
             file_name=trim(trim(seedname)//'_m_orb_sr_AA.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             AA_unit=io_file_unit()
             open(AA_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
             
             file_name=trim( trim(seedname)//'_m_orb_sr_tot.dat')
             write(stdout,'(/,3x,a)') '* '//file_name
             tot_unit=io_file_unit()
             open(tot_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
 
             ! At this point sig_ab contains the conductivity in S/cm (see 
             ! above). Multiply by 100 to convert to S/m (SI). The sum rule is
             !
             ! M_SR^(I)=(hbar/pi.e)int_0^infty sigma(omega).d(omega)
             !         =(1/pi)int_0^infty sigma(E).dE
             !
             ! with the energy in eV (converted from angular freq. to energy
             ! by multiplying by hbar), and converted energy to eV when
             ! dividing by the elementary charge). To get a magnetic moment
             ! per unit cell in units of the bohr magneton, multiply by the
             ! cell volume (in Ang^3), convert to m^3, divide by bohr_magn_SI,
             ! and finally divide by pi:
             ! 
             Morb_conv=100.0_dp*(cell_volume*1.0e-30_dp)/bohr_magn_SI/pi
             
             m_orb_sr=0.0_dp             
             do loop_f=1,nfreq
                m_orb_sr(1,:)=m_orb_sr(1,:)+sig_ab(1,loop_f,:,1)*d_freq
                m_orb_sr(2,:)=m_orb_sr(2,:)+sig_ab(2,loop_f,:,1)*d_freq
                m_orb_sr(3,:)=m_orb_sr(3,:)+sig_ab(3,loop_f,:,1)*d_freq
                freq=optics_min_energy+loop_f*d_freq
                write(DD_unit,'(10E16.8)') freq,m_orb_sr(1,:)*Morb_conv
                write(DA_unit,'(10E16.8)') freq,m_orb_sr(2,:)*Morb_conv
                write(AA_unit,'(10E16.8)') freq,m_orb_sr(3,:)*Morb_conv
                write(tot_unit,'(10E16.8)') freq,Morb_conv*(m_orb_sr(1,:)&
                                                                +m_orb_sr(2,:)&
                                                                +m_orb_sr(3,:))
             end do
             
             close(DD_unit)
             close(DA_unit)
             close(AA_unit)
             close(tot_unit)
             
          end if !T_odd

       end if !eval_sig_ab

       ! ----------!
       ! sigma_abc !
       ! ----------!
       !
       if(eval_sig_abc) then

          write(stdout,'(/,/,1x,a)')&
             '-----------------------------------------------------------------'
          write(stdout,'(1x,a)')&
             'Output data files related to spatial-dispersion optical spectrum:'
          write(stdout,'(1x,a)')&
             '-----------------------------------------------------------------'

          if(T_odd) then
             file_name= trim(trim(seedname)//'-sigS_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat')
          else
             file_name= trim(trim(seedname)//'-rhobw2_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat')
          end if
          write(stdout,'(/,3x,a)') '* '//file_name
          tot_unit=io_file_unit()
          open(tot_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   

          ! After the operation below, sig_abc contains:
          !
          ! If T_odd==T: Im[sigma_{S,abc}] in units of e^2/hbar
          ! If T_odd==F: Re[sigma_{A,abc}]/(hbar.omega) 
          !              in units of (e^2/hbar)/(eV)
          !
          sig_abc=sig_abc/cell_volume
          sig_abc_cut=sig_abc_cut/cell_volume
 
          if(T_odd) then
             ! convert Im[sigma_{S,abc}], which has units of conductance, 
             ! from units of e^2/hbar to units of e^2/h (quantum of conductance)
             fac=2.0_dp*pi
          else
             !
             ! -------------------------------------------------------
             ! Optical rotatory power rho divided by (photon energy)^2
             ! -------------------------------------------------------
             !
             ! rho_ab/(hbar.omega)^2=
             !                 Re[sigma_{A,abc}]/(2.hbar^2.c^2.eps_0.omega)
             ! 
             ! In SI units: rad/(m.J^2)
             fac=elem_charge_SI/(2* hbar_SI**2*speedlight_SI**2*eps0_SI)
             ! Now convert to deg/[mm.(eV)^2]
             fac=fac*(180.0_dp/pi)*elem_charge_SI**2/1000.0_dp
          endif
          sig_abc=fac*sig_abc

         do loop_f=1,nfreq
             freq=optics_min_energy+(loop_f-1)*d_freq
             write(tot_unit,'(5E18.8)') freq,&
                  sig_abc(loop_f,1),&   !orbital, matrix-element-term
                  sig_abc(loop_f,2),&   !orbital, energy-term
                  sig_abc(loop_f,3),&   !spin
                  sum(sig_abc(loop_f,:))
          end do

          close(tot_unit)
          
          ! Now write the cumulative contribution as a function of frequency
          !
          if(T_odd) then
             file_name= trim(trim(seedname)//'-cutsigS_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat')
          else
             file_name= trim(trim(seedname)//'-cutrhobw2_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat')
          end if
          write(stdout,'(/,3x,a)') '* '//file_name
          tot_unit=io_file_unit()
          open(tot_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')  
          sig_abc_cut=fac*sig_abc_cut
          do loop_f=1,nfreq_cut
             write(tot_unit,'(5E18.8)') (loop_f-1)*d_freq_cut,&
                  sig_abc_cut(loop_f,1),sig_abc_cut(loop_f,2),&
                  sig_abc_cut(loop_f,3),sum(sig_abc_cut(loop_f,:))
          end do
          close(tot_unit)

       endif !eval_sig_abc

       ! ------------------------------------------------------------------!
       ! Optical magnetoelectric and electric-quadrupolar response tensors !
       ! ------------------------------------------------------------------!
       !
       if(eval_ME_EQ) then

          write(stdout,'(/,/,1x,a)')&
             '--------------------------------------------------------------'
          write(stdout,'(1x,a)')&
             'Output data files related to optical magnetoelectric spectrum:'
          write(stdout,'(1x,a)')&
             '--------------------------------------------------------------'

          ! Print out the nine elements alpha(i,j)
          ! (only eight are independent, because optical alpha is traceless)
          !
          do j=1,3
             do i=1,3
                file_name= trim(trim(seedname))//'-alpha_'&
                     //achar(119+i)//achar(119+j)//'.dat'
                write(stdout,'(/,3x,a)') '* '//file_name
                alpha_unit(i,j)=io_file_unit()
                open(alpha_unit(i,j),FILE=file_name,&
                     STATUS='UNKNOWN',FORM='FORMATTED')
             enddo
          enddo

          ! After the operation below, alpha contains
          ! alpha_{ab}(omega) in units of e^2/hbar
          !
          alpha_me=alpha_me/cell_volume
          !
          ! Now convert to SI
          !
          alpha_me=alpha_me*elem_charge_SI**2/hbar_SI
          !
          ! At this point we have alpha_ij=del M_j/del E_i, which has units of
          ! admittance (conductance). Following discussion in Coh et al PRB11, 
          ! we now convert to alpha^{EH} (units of inverse velocity) by 
          ! multiplying by the vaccum magnetic permeability mu_0
          !
          alpha_me=alpha_me*4.0_dp*pi*1.0e-7_dp
          !
          ! Finally convert from s/m to ps/m
          !
          alpha_me=alpha_me*1.0e12_dp

         do loop_f=1,nfreq
             freq=optics_min_energy+(loop_f-1)*d_freq
             do j=1,3
                do i=1,3
                   write(alpha_unit(i,j),'(5E18.8)') freq,&
                        alpha_me(1,loop_f,i,j),& !orbital, matrix-element-term
                        alpha_me(2,loop_f,i,j),& !orbital, energy-term
                        alpha_me(3,loop_f,i,j),& !spin
                        sum(alpha_me(:,loop_f,i,j))
                enddo
             enddo
          end do

          ! Print out the ten independent elements of gamma_ijl
          !
          do p=1,10
             file_name= trim(trim(seedname))//'-gamma_'&
                  //achar(119+aa(p))//achar(119+bb(p))//achar(119+cc(p))//'.dat'
             write(stdout,'(/,3x,a)') '* '//file_name
             gamma_unit(p)=io_file_unit()
             open(gamma_unit(p),FILE=file_name,&
                  STATUS='UNKNOWN',FORM='FORMATTED')
          enddo

          ! After the operation below, gamma contains
          ! gamma_{abc}(omega)=(1/3)(sigma_{S,abc}+sigma_{S,cab}+sigma_{S,bca}) 
          ! in units of e^2/hbar (conductance)
          !
          imgamma=imgamma/cell_volume
          !
          ! Now multiply by c.mu_0, the inverse of the vaccum admittance, to 
          ! obtain a dimensionless quantity (see Sec. II.A Coh et al PRB11 for 
          ! discussion of units)
          !
          imgamma=imgamma*speedlight_SI*4.0_dp*pi*1.0e-7

         do loop_f=1,nfreq
             freq=optics_min_energy+(loop_f-1)*d_freq
             do p=1,10
                !*************************
                ! eventually change to 1,2
                !*************************
                write(gamma_unit(p),'(5E18.8)') freq,&
                     imgamma(1,loop_f,p),& !orbital, matrix-element-term
                     imgamma(2,loop_f,p),& !orbital, energy-term
                     imgamma(3,loop_f,p),& !spin (***should vanish***)
                     sum(imgamma(:,loop_f,p))
             enddo
          enddo

          do j=1,3
             do i=1,3
                close(alpha_unit(i,j))
             enddo
          enddo
          do p=1,10
             close(gamma_unit(p))
          enddo

       endif !eval_ME_EQ

       !----------------------------------------------
       ! Static spin-electronic magnetoelectric tensor
       !----------------------------------------------
       !
       if(eval_MEspn) then

          ! alphaspn contains alphaspn_{ab}(omega) in units of e^2/hbar
          ! Now convert to SI
          !
          alphaspn=alphaspn*elem_charge_SI**2/hbar_SI
          alphaspn_cut=alphaspn_cut*elem_charge_SI**2/hbar_SI
          !
          ! At this point we have alpha_ij=del M_j/del E_i, which has units of
          ! admittance (conductance). Following discussion in Coh et al PRB11, 
          ! we now convert to alpha^{EH} (units of inverse velocity) by 
          ! multiplying by the vaccum magnetic permeability mu_0. Finally,
          ! convert from s/m to ps/m
          !
          alphaspn=alphaspn*4.0_dp*pi*1.0e-7_dp*1.0e12_dp
          alphaspn_cut=alphaspn_cut*4.0_dp*pi*1.0e-7_dp*1.0e12_dp

          write(stdout,'(/,/,1x,a)')&
   '---------------------------------------------------------------------------'
          write(stdout,'(1x,a)')&
   'Output data files related to static spin-electronic magnetoelectric tensor:'
          write(stdout,'(1x,a,/)')&
   '---------------------------------------------------------------------------'
          do j=1,3
             do i=1,3
                file_name= trim(trim(seedname))//'-alphaspn_'&
                     //achar(119+i)//achar(119+j)//'.dat'
                write(stdout,'(3x,a)') '* '//file_name
                tot_unit=io_file_unit()
                open(tot_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
                do loop_f=1,nfreq_cut
                   write(tot_unit,'(2(E18.8,1x))') (loop_f-1)*d_freq_cut,&
                        alphaspn_cut(loop_f,i,j)
                enddo
                close(tot_unit)
             enddo
          enddo

          write(stdout,'(/,/,1x,a)')&
               'Spin-electronic static magnetoelectric tensor (ps/m)'
          write(stdout,'(1x,a,/)')&
               '============================================='
          do j=1,3
             do i=1,3
                write(stdout,'(a,f12.6)')&
                     'alpha_'//achar(119+i)//achar(119+j)//'=',&
                     alphaspn(i,j)
             enddo
          enddo
                
       endif !eval_MEspn
       
       write(stdout,'(1x,a)') ' '
       if(wanint_kpoint_file) then

          write(stdout,'(1x,a47,i10,a)')&
               'Nominal interpolation mesh in IBZ: ',&
               sum(num_int_kpts_on_node),' points'

          if(eval_ahe .or. eval_orb) then
             if(eval_ahe) then
                write(stdout,'(1x,a47,3x,f8.4,a)')&
                     'Adaptive refinement triggered when Omega(k) >',&
                     optics_adaptive_thresh,' Ang^2'
             else
                write(stdout,'(1x,a48,f8.4,a)')&
                     'Adaptive refinement triggered when k-integrand >',&
                     optics_adaptive_thresh,' eV.Ang^2'
             end if
             write(stdout,'(1x,a47,i7,a,f8.4,a)')&
                  'How many points triggered adaptive refinement: ',&
                  adpt_counter,' (',&
                  100*real(adpt_counter,dp)/sum(num_int_kpts_on_node),' %)'
             if(optics_adaptive_pts < 10) then
                write(stdout,'(1x,a47,4x,i1,a,i1,a,i1)')&
                     'Adaptive mesh: ',optics_adaptive_pts,'x',&
                     optics_adaptive_pts,'x',optics_adaptive_pts
             elseif(optics_adaptive_pts < 100) then
                write(stdout,'(1x,a47,i2,a,i1,a,i1)')&
                     'Adaptive mesh: ',optics_adaptive_pts,'x',&
                     optics_adaptive_pts,'x',optics_adaptive_pts
             else
                write(stdout,'(1x,a47,i3,a,i1,a,i1)')&
                     'Adaptive mesh: ',optics_adaptive_pts,'x',&
                     optics_adaptive_pts,'x',optics_adaptive_pts
             end if
             write(stdout,'(1x,a47,i10)') 'Total number of points: ',&
                  sum(num_int_kpts_on_node)-adpt_counter+&
                  adpt_counter*optics_adaptive_pts**3
          end if

       else

          if(optics_num_points < 10) then
             write(stdout,'(1x,a47,i1,a,i1,a,i1,a,i12,a)')&
                  'Nominal interpolation mesh in full BZ: ',&
                  optics_num_points,'x',optics_num_points,'x',&
                  optics_num_points,'=',optics_num_points**3,' points'
          elseif(optics_num_points < 100) then
             write(stdout,'(1x,a47,i2,a,i2,a,i2,a,i10,a)')&
                  'Nominal interpolation mesh in full BZ: ',&
                  optics_num_points,'x',optics_num_points,'x',&
                  optics_num_points,'=',optics_num_points**3,' points'
          elseif(optics_num_points < 1000) then
             write(stdout,'(1x,a47,i3,a,i3,a,i3,a,i12,a)')&
                  'Nominal interpolation mesh in full BZ: ',&
                  optics_num_points,'x',optics_num_points,'x',&
                  optics_num_points,'=',optics_num_points**3,' points'
          else
             write(stdout,'(1x,a47,i4,a,i4,a,i4,a,i12,a)')&
                  'Nominal interpolation mesh in full BZ: ',&
                  optics_num_points,'x',optics_num_points,'x',&
                  optics_num_points,'=',optics_num_points**3,' points'
          end if

          if(eval_ahe .or. eval_orb) then
             if(eval_ahe) then
                write(stdout,'(1x,a47,f8.4,a)')&
                     'Adaptive refinement triggered when Omega(k) >',&
                     optics_adaptive_thresh,' Ang^2'
             else
                write(stdout,'(1x,a48,f8.4,a)')&
                     'Adaptive refinement triggered when k-integrand >',&
                     optics_adaptive_thresh,' eV.Ang^2'
             end if
             write(stdout,'(1x,a47,i7,a,f8.4,a)')&
                  'How many points triggered adaptive refinement: ',&
                  adpt_counter,' (',&
                  100*real(adpt_counter,dp)/optics_num_points**3,' %)'
             if(optics_adaptive_pts < 10) then
                write(stdout,'(1x,a47,i1,a,i1,a,i1)')&
                     'Adaptive mesh: ',optics_adaptive_pts,'x',&
                     optics_adaptive_pts,'x',optics_adaptive_pts
             elseif(optics_adaptive_pts < 100) then
                write(stdout,'(1x,a47,i2,a,i1,a,i1)')&
                     'Adaptive mesh: ',optics_adaptive_pts,'x',&
                     optics_adaptive_pts,'x',optics_adaptive_pts
             else
                write(stdout,'(1x,a47,i3,a,i1,a,i1)')&
                     'Adaptive mesh: ',optics_adaptive_pts,'x',&
                     optics_adaptive_pts,'x',optics_adaptive_pts
                write(stdout,'(1x,a47,i12)') 'Total number of points: ',&
                     optics_num_points**3-adpt_counter+&
                     adpt_counter*optics_adaptive_pts**3
             end if
          end if

       end if !wanint_kpoint_file

       if (timing_level>1) call io_stopwatch('berry_wanint: berry',2)

    end if !on_root

  end subroutine berry

  subroutine get_imf_ab_k(kpt,imf)
  !=====================================!
  !                                     !
  ! Calculates -2Im[f_{alpha beta}(k)], !
  ! Eq.33 Ceresoli et al PRB 74, 024408 !
  !                                     !
  !=====================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i
    use w90_utility, only       : utility_diagonalize,utility_re_tr,utility_im_tr
    use w90_parameters, only    : alpha,beta,num_wann,omega_from_FF
    use w90_wanint_common, only : fourier_R_to_k
    use w90_wan_ham, only   : get_JJplus,get_JJminus,get_occ_mat
    use w90_get_oper, only      : HH_R,AA_R,FF_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: imf(3)

    ! Physics
    !
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f(:,:)
    complex(kind=dp), allocatable :: g(:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: FF(:,:)
    complex(kind=dp), allocatable :: OOmega(:,:)
    complex(kind=dp), allocatable :: JJminus_a(:,:)
    complex(kind=dp), allocatable :: JJplus_b(:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)


    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,2))
    allocate(UU(num_wann,num_wann))
    allocate(f(num_wann,num_wann))
    allocate(g(num_wann,num_wann))
    allocate(JJminus_a(num_wann,num_wann))
    allocate(JJplus_b(num_wann,num_wann))
    allocate(AA(num_wann,num_wann,2))
    allocate(OOmega(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    if(omega_from_FF) allocate(FF(num_wann,num_wann))

    ! gather W-gauge matrix objects

    ! hamiltonian
    call fourier_R_to_k(kpt,HH_R,HH,0)
    call utility_diagonalize(HH,num_wann,eig,UU)

    ! occupation f, and g=1-f
    call get_occ_mat(eig,UU,f,g)

    ! JJminus_a and JJplus_b
    call fourier_R_to_k(kpt,HH_R,mdum,alpha) 
    call get_JJminus(mdum,UU,eig,JJminus_a(:,:))
    call fourier_R_to_k(kpt,HH_R,mdum,beta)
    call get_JJplus(mdum,UU,eig,JJplus_b(:,:))

    ! AA_a = i<u|del_a u> (a=alpha,beta)
    call fourier_R_to_k(kpt,AA_R(:,:,:,alpha),AA(:,:,1),0)
    call fourier_R_to_k(kpt,AA_R(:,:,:,beta),AA(:,:,2),0)

    if(omega_from_FF) then
       ! First equality in Eq.(31) LVTS12: OOmega = i(F-F^dag)
       call fourier_R_to_k(kpt,FF_R(:,:,:,alpha,beta),FF(:,:),0)
       OOmega=cmplx_i*(FF-conjg(transpose(FF)))
    else
       ! Second equality in Eq.(31) LVTS12:
       ! OOmega = del_alpha AA_beta - del_beta AA_alpha
       call fourier_R_to_k(kpt,AA_R(:,:,:,beta),mdum,alpha)
       Oomega=mdum
       call fourier_R_to_k(kpt,AA_R(:,:,:,alpha),mdum,beta)
       OOmega=Oomega-mdum
    endif

    ! Berry curvature of occ. states as a trace, Eq.(51) LVTS12

    ! J0 (Omega_bar) term
    mdum=matmul(f,OOmega)
    imf(1)=utility_re_tr(mdum)

    ! J1 (DA) term
    mdum=matmul(AA(:,:,1),JJplus_b)+matmul(JJminus_a,AA(:,:,2))
    imf(2)=-2.0_dp*utility_im_tr(mdum)

    ! J2 (DD) term
    mdum=matmul(JJminus_a,JJplus_b)
    imf(3)=-2.0_dp*utility_im_tr(mdum)

  end subroutine get_imf_ab_k


  subroutine get_imh_ab_k(kpt,imh)
  !=====================================!
  !                                     !
  ! Calculates -2Im[h_{alpha beta}(k)], !
  ! Eq.35 Ceresoli et al PRB 74, 024408 !
  !                                     !
  !=====================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i
    use w90_utility, only       : utility_diagonalize,utility_re_tr,utility_im_tr
    use w90_parameters, only    : alpha,beta,num_wann,omega_from_FF
    use w90_wanint_common, only : fourier_R_to_k
    use w90_wan_ham, only   : get_JJplus,get_JJminus,get_occ_mat
    use w90_get_oper, only      : HH_R,AA_R,FF_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: imh(3)

    ! Physics
    !
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f(:,:)
    complex(kind=dp), allocatable :: g(:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: FF(:,:)
    complex(kind=dp), allocatable :: OOmega(:,:)
    complex(kind=dp), allocatable :: JJminus_a(:,:)
    complex(kind=dp), allocatable :: JJplus_b(:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)


    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
    allocate(f(num_wann,num_wann))
    allocate(g(num_wann,num_wann))
    allocate(JJminus_a(num_wann,num_wann))
    allocate(JJplus_b(num_wann,num_wann))
    allocate(AA(num_wann,num_wann,2))
    allocate(OOmega(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    if(omega_from_FF) allocate(FF(num_wann,num_wann))

    ! hamiltonian
    call fourier_R_to_k(kpt,HH_R,HH,0)
    call utility_diagonalize(HH,num_wann,eig,UU)

    ! occupation f, and g=1-f
    call get_occ_mat(eig,UU,f,g)

    ! JJminus_a and JJplus_b
    call fourier_R_to_k(kpt,HH_R,mdum,alpha)
    call get_JJminus(mdum,UU,eig,JJminus_a(:,:))
    call fourier_R_to_k(kpt,HH_R,mdum,beta)
    call get_JJplus(mdum,UU,eig,JJplus_b(:,:))

    ! AA_a = i<u|del_a u> (a=alpha,beta)
    call fourier_R_to_k(kpt,AA_R(:,:,:,alpha),AA(:,:,1),0)
    call fourier_R_to_k(kpt,AA_R(:,:,:,beta),AA(:,:,2),0)

    if(omega_from_FF) then
       ! First equality in Eq.(31) LVTS12: OOmega = i(F-F^dag)
       call fourier_R_to_k(kpt,FF_R(:,:,:,alpha,beta),FF(:,:),0)
       OOmega=cmplx_i*(FF-conjg(transpose(FF)))
    else
       ! Second equality in Eq.(31) LVTS12:
       ! OOmega = del_alpha AA_beta - del_beta AA_alpha
       call fourier_R_to_k(kpt,AA_R(:,:,:,beta),mdum,alpha)
       Oomega=mdum
       call fourier_R_to_k(kpt,AA_R(:,:,:,alpha),mdum,beta)
       OOmega=Oomega-mdum
    endif

    mdum=matmul(f,matmul(HH,OOmega))
    imh(1)=utility_re_tr(mdum)
    mdum=matmul(f,matmul(HH,matmul(AA(:,:,1),matmul(f,AA(:,:,2)))))
    imh(1)=imh(1)+2.0_dp*utility_im_tr(mdum)

    mdum=matmul(HH,matmul(AA(:,:,1),JJplus_b))&
        +matmul(HH,matmul(JJminus_a,AA(:,:,2)))
    imh(2)=-2.0_dp*utility_im_tr(mdum)

    mdum=matmul(HH,matmul(JJminus_a,JJplus_b))
    imh(3)=-2.0_dp*utility_im_tr(mdum)

  end subroutine get_imh_ab_k


  subroutine get_img_ab_k(kpt,img)
  !=====================================!
  !                                     !
  ! Calculates -2Im[g_{alpha beta}(k)], !
  ! Eq.34 Ceresoli et al PRB 74, 024408 !
  !                                     !
  !=====================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i
    use w90_utility, only   : utility_diagonalize,utility_re_tr,utility_im_tr
    use w90_parameters, only    : alpha,beta,num_wann
    use w90_wanint_common, only : fourier_R_to_k
    use w90_wan_ham, only   : get_JJplus,get_JJminus,get_occ_mat
    use w90_get_oper, only      : HH_R,AA_R,BB_R,CC_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: img(3)

    ! Physics
    !
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f(:,:)
    complex(kind=dp), allocatable :: g(:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: BB(:,:,:)
    complex(kind=dp), allocatable :: CC(:,:)
    complex(kind=dp), allocatable :: LLambda(:,:)
    complex(kind=dp), allocatable :: JJminus_a(:,:)
    complex(kind=dp), allocatable :: JJminus_b(:,:)
    complex(kind=dp), allocatable :: JJplus_b(:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)

    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
    allocate(f(num_wann,num_wann))
    allocate(g(num_wann,num_wann))
    allocate(JJminus_a(num_wann,num_wann))
    allocate(JJminus_b(num_wann,num_wann))
    allocate(JJplus_b(num_wann,num_wann))
    allocate(AA(num_wann,num_wann,2))
    allocate(BB(num_wann,num_wann,2))
    allocate(CC(num_wann,num_wann))
    allocate(LLambda(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    ! hamiltonian
    call fourier_R_to_k(kpt,HH_R,HH,0)
    call utility_diagonalize(HH,num_wann,eig,UU)

    ! occupation f, and g=1-f
    call get_occ_mat(eig,UU,f,g)

    ! JJminus_a, JJminus_b, and JJplus_b
    call fourier_R_to_k(kpt,HH_R,mdum,alpha)
    call get_JJminus(mdum,UU,eig,JJminus_a(:,:))
    call fourier_R_to_k(kpt,HH_R,mdum,beta)
    call get_JJminus(mdum,UU,eig,JJminus_b(:,:))
    call get_JJplus(mdum,UU,eig,JJplus_b(:,:))

    ! AA_a = i<u|del_a u> 
    call fourier_R_to_k(kpt,AA_R(:,:,:,alpha),AA(:,:,1),0)
    call fourier_R_to_k(kpt,AA_R(:,:,:,beta),AA(:,:,2),0)

    ! BB_a = i<u|H|del_a u>
    call fourier_R_to_k(kpt,BB_R(:,:,:,alpha),BB(:,:,1),0)
    call fourier_R_to_k(kpt,BB_R(:,:,:,beta),BB(:,:,2),0)

    ! CC_ab = i<del_a u|H|del_b u>
    call fourier_R_to_k(kpt,CC_R(:,:,:,alpha,beta),CC(:,:),0)

    ! LLambda_ab = [CC - CC^dag]_ab
    LLambda=cmplx_i*(CC-conjg(transpose(CC)))

    mdum=matmul(f,LLambda)
    img(1)=utility_re_tr(mdum)
    mdum=matmul(f,matmul(HH,matmul(AA(:,:,1),matmul(f,AA(:,:,2)))))
    img(1)=img(1)-2.0_dp*utility_im_tr(mdum)

    mdum=matmul(JJminus_a,BB(:,:,2))-matmul(JJminus_b,BB(:,:,1))
    img(2)=-2.0_dp*utility_im_tr(mdum)

    mdum=matmul(JJminus_a,matmul(HH,JJplus_b))
    img(3)=-2.0_dp*utility_im_tr(mdum)

  end subroutine get_img_ab_k


  !===========================================================!
  !                   PRIVATE PROCEDURES                      ! 
  !===========================================================!

  subroutine get_sig_ab_k(kpt,sig_ab_k,sig_ab_over_freq_k,jdos_k)
  !================================================================!
  !                                                                !
  ! Optical conductivity (absorptive) in the long-wavelength limit !
  !                                                                !
  !================================================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i
    use w90_utility, only       : utility_diagonalize,utility_rotate,w0gauss
    use w90_parameters, only    : num_wann,optics_min_energy,&
                                  optics_num_points,alpha,beta,&
                                  adpt_smr_steps,adpt_smr_width,&
                                  spn_decomp,fermi_energy
    use w90_wanint_common, only : get_occ,kmesh_spacing,fourier_R_to_k
    use w90_wan_ham, only   : get_D_h_a,get_deleig_a
    use w90_get_oper, only      : HH_R,AA_R
    use w90_spin_wanint, only   : get_spn_nk

    ! Arguments
    !
    real(kind=dp),                     intent(in)  :: kpt(3)
    real(kind=dp), dimension(:,:,:,:), intent(out) :: sig_ab_k
    real(kind=dp), dimension(:,:,:,:), intent(out) :: sig_ab_over_freq_k
    real(kind=dp), dimension(:,:,:),   intent(out) :: jdos_k

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: AA_bar(:,:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    
    ! Adaptive smearing
    !
    real(kind=dp) :: del_eig(num_wann,3),joint_level_spacing,&
                     smear,Delta_k,arg
 
    ! Misc/Dummy
    !
    integer          :: i,j,loop_f,loop_s,case
    real(kind=dp)    :: rdum,rvdum(3),eig(num_wann),occ(num_wann),occ_prod,&
                        freq,A2_ij(3),spn_nk(num_wann) 

    allocate(HH(num_wann,num_wann))
    allocate(D_h(num_wann,num_wann,2))
    allocate(UU(num_wann,num_wann))
    allocate(AA_bar(num_wann,num_wann,2))
    allocate(mdum(num_wann,num_wann))

    call fourier_R_to_k(kpt,HH_R,HH,0) 
    call utility_diagonalize(HH,num_wann,eig,UU) 
    call fourier_R_to_k(kpt,HH_R,mdum,alpha) 
    call get_D_h_a(mdum,UU,eig,D_h(:,:,1))
    call fourier_R_to_k(kpt,HH_R,mdum,beta) 
    call get_D_h_a(mdum,UU,eig,D_h(:,:,2))

    ! band gradients
    call fourier_R_to_k(kpt,HH_R,mdum,1) 
    call get_deleig_a(del_eig(:,1),eig,mdum,UU)
    call fourier_R_to_k(kpt,HH_R,mdum,2) 
    call get_deleig_a(del_eig(:,2),eig,mdum,UU)
    call fourier_R_to_k(kpt,HH_R,mdum,3) 
    call get_deleig_a(del_eig(:,3),eig,mdum,UU)

    call fourier_R_to_k(kpt,AA_R(:,:,:,alpha),mdum(:,:),0)
    AA_bar(:,:,1)=utility_rotate(mdum,UU,num_wann)
    call fourier_R_to_k(kpt,AA_R(:,:,:,beta),mdum(:,:),0)
    AA_bar(:,:,2)=utility_rotate(mdum,UU,num_wann)

    call get_occ(eig,occ,fermi_energy)

    ! Find spacing of interpolation mesh 
    !
    Delta_k=kmesh_spacing(optics_num_points)
    
    ! Get spin projections along chosen quantization axis for every band
    !
    if(spn_decomp) call get_spn_nk(kpt,spn_nk)

    sig_ab_k=0.0_dp
    sig_ab_over_freq_k=0.0_dp
    jdos_k=0.0_dp
    case=0 ! just so that compilation with 'g95 -Wall' doesn't complain!
           ! the value 0 for 'case' is actually meaningless
    do i=1,num_wann
       do j=1,num_wann
          if(spn_decomp) then
             if(spn_nk(i)>=0 .and. spn_nk(j)>=0) then 
                case=1 ! up --> up transition 
                       ! NOTE: These are *minority* spins if magnetization 
                       !       is positive
             elseif(spn_nk(i)<0 .and. spn_nk(j)<0) then 
                case=2 ! down --> down
             else
                case=3 ! spin-flip transition
             end if
          end if
          occ_prod=occ(i)*(1.0_dp-occ(j))
!          if(occ_prod < eps_occ) cycle
          rvdum(:)=del_eig(j,:)-del_eig(i,:)
          joint_level_spacing=sqrt(dot_product(rvdum(:),rvdum(:)))*Delta_k
          if(T_odd) then 
             ! Dichoric [Eqs.(39),(42) YWVS07]
             call get_A2(D_h,AA_bar,i,j,A2_ij,'im')
          else
             ! Ordinary
             call get_A2(D_h,AA_bar,i,j,A2_ij,'re')
          end if
          do loop_s=1,adpt_smr_steps
             
             ! Except for the factor 1/sqrt(2), this is Eq.(35) YWVS07
             ! !!!UNDERSTAND THAT FACTOR!!!
             !
             smear=joint_level_spacing*adpt_smr_width(loop_s)/sqrt(2.0_dp)
             do loop_f=1,nfreq   
                
                ! NOTE: Skipping optics_min_energy. When it is zero can give
                !       problems in the case of ordinary spectrum 
                !
                freq=optics_min_energy+loop_f*d_freq
                arg=(eig(j)-eig(i)-freq)/smear
                if(abs(arg) > 10.0_dp) then !optimisation
                   cycle
                else
                   !
                   ! Adaptive broadening of the delta-function in Eq.(39) YWVS07
                   !
                   rdum=w0gauss(arg)/smear
                   !
                   ! Fermi occupancy factor
                   !
                   rdum=rdum*occ_prod
                end if
                !
                ! Joint density of states
                !
                jdos_k(loop_f,loop_s,1)=jdos_k(loop_f,loop_s,1)+rdum
                if(spn_decomp) jdos_k(loop_f,loop_s,1+case)&
                     =jdos_k(loop_f,loop_s,1+case)+rdum
                !
                ! Optical conductivity
                !
                ! NOTE: Minus sign in Eq.(39) YWVS07 is *wrong*
                
                ! sigma(freq)/freq: Exclude factor of freq in Eq.(39) YWVS07
                !
                sig_ab_over_freq_k(1,loop_f,loop_s,1)&
                     =sig_ab_over_freq_k(1,loop_f,loop_s,1)+A2_ij(1)*rdum
                if(spn_decomp) sig_ab_over_freq_k(1,loop_f,loop_s,1+case)&
                     =sig_ab_over_freq_k(1+case,loop_f,loop_s,1)+A2_ij(1)*rdum
                sig_ab_over_freq_k(2,loop_f,loop_s,1)&
                     =sig_ab_over_freq_k(2,loop_f,loop_s,1)+A2_ij(2)*rdum
                sig_ab_over_freq_k(3,loop_f,loop_s,1)&
                     =sig_ab_over_freq_k(3,loop_f,loop_s,1)+A2_ij(3)*rdum
                sig_ab_over_freq_k(4,loop_f,loop_s,1)&
                     =sig_ab_over_freq_k(4,loop_f,loop_s,1)+sum(A2_ij(1:3))*rdum
                if(spn_decomp) then
                   sig_ab_over_freq_k(2,loop_f,loop_s,1+case)&
                        =sig_ab_over_freq_k(2,loop_f,loop_s,1+case)+A2_ij(2)*rdum
                   sig_ab_over_freq_k(3,loop_f,loop_s,1+case)&
                        =sig_ab_over_freq_k(3,loop_f,loop_s,1+case)+A2_ij(3)*rdum
                   sig_ab_over_freq_k(4,loop_f,loop_s,1+case)&
                        =sig_ab_over_freq_k(4,loop_f,loop_s,1+case)&
                        +sum(A2_ij(1:3))*rdum
                end if
                !
                ! sigma(freq): Include factor of freq in Eq.(39) YWVS07
                !
                rdum=rdum*(eig(j)-eig(i))
                sig_ab_k(1,loop_f,loop_s,1)=sig_ab_k(1,loop_f,loop_s,1)&
                     +A2_ij(1)*rdum
                if(spn_decomp) sig_ab_k(1,loop_f,loop_s,1+case)&
                     =sig_ab_k(1,loop_f,loop_s,1+case)+A2_ij(1)*rdum
                sig_ab_k(2,loop_f,loop_s,1)=&
                     sig_ab_k(2,loop_f,loop_s,1)+A2_ij(2)*rdum
                sig_ab_k(3,loop_f,loop_s,1)=sig_ab_k(3,loop_f,loop_s,1)&
                     +A2_ij(3)*rdum
                sig_ab_k(4,loop_f,loop_s,1)=sig_ab_k(4,loop_f,loop_s,1)&
                     +sum(A2_ij(1:3))*rdum
                if(spn_decomp) then
                   sig_ab_k(2,loop_f,loop_s,1+case)&
                        =sig_ab_k(2,loop_f,loop_s,1+case)+A2_ij(2)*rdum
                   sig_ab_k(3,loop_f,loop_s,1+case)&
                        =sig_ab_k(3,loop_f,loop_s,1+case)+A2_ij(3)*rdum
                   sig_ab_k(4,loop_f,loop_s,1+case)&
                        =sig_ab_k(4,loop_f,loop_s,1+case)+sum(A2_ij(1:3))*rdum
                end if
             end do !loop_f
          end do !loop_s
       end do !j
    end do !i

  end subroutine get_sig_ab_k


  subroutine get_A2(D_h,AA_bar,i,j,A2_ij,imre)
  !=================================================================!
  !                                                                 !
  ! Evaluate the real or imaginary part of A_{ij,alpha}A_{ji,beta}, !
  ! Eq.(42) YWVS07                                                  !
  !                                                                 !
  ! A2_ij(1) is the DD term                                         !
  ! A2_ij(2) is the DA term                                         !
  ! A2_ij(3) is the AA term                                         !
  !=================================================================!

    use w90_constants, only     : dp,cmplx_i,cmplx_0
    use w90_io, only            : io_error
!    use w90_utility, only       : utility_rotate,utility_rotate_diag,utility_commutator_diag

    ! Arguments
    !
    complex(kind=dp), dimension(:,:,:), intent(in)  :: D_h
    complex(kind=dp), dimension(:,:,:), intent(in)  :: AA_bar
    integer,                            intent(in)  :: i,j
    real(kind=dp),    dimension(:,:,:), intent(out) :: A2_ij(3)
    character(len=2),                   intent(in)  :: imre
 
    ! Misc/Dummy
    !
    complex(kind=dp) :: D_ji_a,D_ji_b,A_ji_a,A_ji_b

    D_ji_a=D_h(j,i,1)
    D_ji_b=D_h(j,i,2)
    A_ji_a=AA_bar(j,i,1)
    A_ji_b=AA_bar(j,i,2)
    
    if(imre=='im') then
       !
       ! DD term in Eq.(42) YWVS07
       !
       A2_ij(1)=aimag(conjg(D_ji_a)*D_ji_b) 
       !
       ! DA terms
       !
       A2_ij(2)=aimag(cmplx_i*(conjg(A_ji_a)*D_ji_b&
                              -conjg(D_ji_a)*A_ji_b))
       !
       ! AA term
       !
       A2_ij(3)=aimag(conjg(A_ji_a)*A_ji_b)   
    elseif(imre=='re') then
       A2_ij(1)=aimag(cmplx_i*conjg(D_ji_a)*D_ji_b)
       A2_ij(2)=aimag(conjg(D_ji_a)*A_ji_b&
            -conjg(A_ji_a)*D_ji_b)
       A2_ij(3)=aimag(cmplx_i*conjg(A_ji_a)*A_ji_b) 
    else 
       call io_error('error in value of variable imre in get_A2')
    end if

  end subroutine get_A2


  subroutine get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
  !==================================================================!
  !                                                                  !
  ! Optical conductivity of an insulator at first order in           !
  ! the wavevector of light. The reactive component is               !
  ! evaluated, for frequencies assumed to be in the gap.             !
  !                                                                  !
  ! * If T_odd==F, on output sig_abc_k contains Im[sigma_{S,abc}(k)] !
  !   ("S" stands for symmetric part under a<-->b)                   !
  !                                                                  !
  ! * If T_odd=T, on output sig_abc_k contains                       !
  !   Re[sigma_{A,abc}(k)]/(hbar.omega)                              !
  !   ("A" stands for antisymmetric under a<-->b)                    !
  !                                                                  !
  ! NOTE: Here a=alpha, b=beta, c=gamma                              !
  !                                                                  !
  !==================================================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i,bohr,eV_au
    use w90_utility, only   : utility_diagonalize,utility_rotate
    use w90_parameters, only    : num_wann,optics_min_energy,fermi_energy,&
                                  optics_num_points,alpha,beta,gamma,&
                                  have_disentangled,dis_froz_max,&
                                  sigma_abc_onlyorb
    use w90_wanint_common, only : fourier_R_to_k
    use w90_wan_ham, only   : get_D_h_a,get_deleig_a
    use w90_get_oper, only      : HH_R,AA_R,BB_R,CC_R,FF_R,SS_R

    ! Arguments
    !
    real(kind=dp),                 intent(in)  :: kpt(3)
    real(kind=dp), dimension(:,:), intent(out) :: sig_abc_k
    real(kind=dp), dimension(:,:), intent(out) :: sig_abc_cut_k

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: HH_h(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: AA_bar(:,:,:)
    complex(kind=dp), allocatable :: BB_bar(:,:,:)
    complex(kind=dp), allocatable :: CC_bar(:,:,:,:)
    complex(kind=dp), allocatable :: FF_bar(:,:,:,:)
    complex(kind=dp), allocatable :: AA_h(:,:,:)
    complex(kind=dp), allocatable :: BB_h(:,:,:)
    complex(kind=dp), allocatable :: CC_h(:,:,:,:)
    complex(kind=dp), allocatable :: FF_h(:,:,:,:)
    complex(kind=dp), allocatable :: SS_h(:,:,:)
    complex(kind=dp), allocatable :: KKspn_h(:,:,:,:)
    complex(kind=dp), allocatable :: KKorb_h(:,:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)

    integer                       :: i,j,n,m,loop_f
    real(kind=dp)                 :: eig(num_wann),del_eig(num_wann,3),eig_cut,&
                                     sig_matel(2),sig_en,freq,omega_mn,&
                                     omega_fac_matel,omega_fac_en
    complex(kind=dp)              :: cdum

    allocate(HH(num_wann,num_wann))
    allocate(HH_h(num_wann,num_wann))
    allocate(D_h(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(AA_bar(num_wann,num_wann,3))
    allocate(BB_bar(num_wann,num_wann,3))
    allocate(CC_bar(num_wann,num_wann,3,3))
    allocate(FF_bar(num_wann,num_wann,3,3))
    allocate(AA_h(num_wann,num_wann,3))
    allocate(BB_h(num_wann,num_wann,3))
    allocate(CC_h(num_wann,num_wann,3,3))
    allocate(FF_h(num_wann,num_wann,3,3))
    if(.not.(sigma_abc_onlyorb)) allocate(SS_h(num_wann,num_wann,3))
    allocate(KKspn_h(num_wann,num_wann,3,3))
    allocate(KKorb_h(num_wann,num_wann,3,3))
    allocate(mdum(num_wann,num_wann))

    call fourier_R_to_k(kpt,HH_R,HH,0)
    call utility_diagonalize(HH,num_wann,eig,UU)
    do i=1,3
       ! Matrix objects with one Cartesian index
       call fourier_R_to_k(kpt,HH_R,mdum,i) ! del_i HH
       call get_D_h_a(mdum,UU,eig,D_h(:,:,i))
       call get_deleig_a(del_eig(:,i),eig,mdum,UU)
       call fourier_R_to_k(kpt,AA_R(:,:,:,i),mdum(:,:),0)
       AA_bar(:,:,i)=utility_rotate(mdum,UU,num_wann)
       call fourier_R_to_k(kpt,BB_R(:,:,:,i),mdum(:,:),0)
       BB_bar(:,:,i)=utility_rotate(mdum,UU,num_wann)
       if(.not.(sigma_abc_onlyorb)) then
          call fourier_R_to_k(kpt,SS_R(:,:,:,i),mdum(:,:),0)
          ! SS_h = SS_bar (SS is "gauge-covariant" in the sense of WYSV06)
          SS_h(:,:,i)=utility_rotate(mdum,UU,num_wann)
       endif
       do j=1,3
          ! Matrix objects with two Cartesian indices
          call fourier_R_to_k(kpt,CC_R(:,:,:,j,i),mdum(:,:),0)
          CC_bar(:,:,j,i)=utility_rotate(mdum,UU,num_wann)
          call fourier_R_to_k(kpt,FF_R(:,:,:,j,i),mdum(:,:),0)
          FF_bar(:,:,j,i)=utility_rotate(mdum,UU,num_wann)
       enddo
    enddo
    
    ! Hamiltonian-gauge matrices
    !
    HH_h=cmplx_0
    AA_h=cmplx_0
    CC_h=cmplx_0
    FF_h=cmplx_0
    do i=1,3
       HH_h(i,i)=eig(i)
       ! Eq.(25) WYSV06
       AA_h(:,:,i)=AA_bar(:,:,i)+cmplx_i*D_h(:,:,i)
    enddo
    do j=1,3
       do i=1,j
          ! See 2012-02-06 notes, p.4
          FF_h(:,:,i,j)=FF_bar(:,:,i,j)&
            +cmplx_i*matmul(AA_bar(:,:,i),D_h(:,:,j))&
            +cmplx_i*matmul(D_h(:,:,i),AA_bar(:,:,j))&
            -matmul(D_h(:,:,i),D_h(:,:,j))
          FF_h(:,:,j,i)=conjg(transpose(FF_h(:,:,i,j)))
          CC_h(:,:,i,j)=CC_bar(:,:,i,j)&
            +cmplx_i*matmul(conjg(transpose(BB_bar(:,:,i))),D_h(:,:,j))&
            +cmplx_i*matmul(D_h(:,:,i),BB_bar(:,:,j))&
            -matmul(D_h(:,:,i),matmul(HH_h,D_h(:,:,j)))
          CC_h(:,:,j,i)=conjg(transpose(CC_h(:,:,i,j)))
       enddo
    enddo

    ! Evaluate, in the H-gauge of  the projected subspace, the matrices 
    ! K_ij=hbar.B_ij, where B_ij is defined in Eqs. (35-37) MS10.
    !
    ! ------------
    ! Orbital part
    ! ------------
    !
    do j=1,3
       do i=1,3
          do m=1,num_wann
             do n=1,m
                cdum=-cmplx_i*(del_eig(n,i)+del_eig(m,i))*AA_h(n,m,j)&
                     +CC_h(n,m,j,i)-CC_h(n,m,i,j)&
                     +eig(n)*FF_h(n,m,i,j)-eig(m)*FF_h(n,m,j,i)
                KKorb_h(n,m,i,j)=cdum/2.0_dp
                !
                ! Sanity check: result unchanged upon removing the
                ! next line and extending loop over n from 1 to num_wann
                KKorb_h(m,n,i,j)=-conjg(KKorb_h(n,m,i,j))
             enddo !n
          enddo !m
       enddo  !i
    enddo !j
    !
    ! ---------
    ! Spin part
    ! ---------
    !
    KKspn_h=cmplx_0 
    !Diagonal elements in Cartesian indices will remain zero matrices
    if(.not.sigma_abc_onlyorb) then
       KKspn_h(:,:,1,2)=SS_h(:,:,3)
       KKspn_h(:,:,2,1)=-SS_h(:,:,3)
       KKspn_h(:,:,1,3)=-SS_h(:,:,2)
       KKspn_h(:,:,3,1)=SS_h(:,:,2)
       KKspn_h(:,:,2,3)=SS_h(:,:,1)
       KKspn_h(:,:,3,2)=-SS_h(:,:,1)
       ! Multiply by (i/2)(hbar^2)/m_e. In atomic units hbar=m_e=1, and we
       ! would get i/2. But our unit of energy is eV, not Hartee, and the 
       ! unit of length is Angstrom, not Bohr radius. Since hbar^2/m_e has 
       ! units of [E][L]^2, need to multiply by 
       ! 1 Ha = 27.21138387 eV and by (1 bohr)^2 = 0.529^2 Ang^2
       KKspn_h=KKspn_h*(cmplx_i/2.0_dp)*bohr**2/eV_au
    endif
    
    if(have_disentangled) then
       eig_cut=dis_froz_max+0.001_dp !some small number 
    else
       eig_cut=maxval(eig)+0.001_dp
    endif
    
    sig_abc_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann !conduction, inside frozen window
          if(eig(m)<fermi_energy.or.eig(m)>eig_cut) cycle
          omega_mn=eig(m)-eig(n)
          do loop_f=1,nfreq
             freq=optics_min_energy+(loop_f-1)*d_freq
             if(T_odd) then
                !
                ! ------------------------------
                ! Im[sigma_{S,alpha beta gamma}]
                ! ------------------------------
                !
                omega_fac_matel=2.0_dp*omega_mn/(omega_mn**2-freq**2)
                omega_fac_en=2.0_dp*omega_mn**3/(omega_mn**2-freq**2)**2
                !
                ! "Matrix element term", Eq. (30) MS10. Integrand,
                ! including factor of two in front.
                !
                ! orbital
                sig_matel(1)=omega_fac_matel&
                     *aimag(AA_h(m,n,beta)*KKorb_h(n,m,alpha,gamma)&
                           +AA_h(m,n,alpha)*KKorb_h(n,m,beta,gamma))
                ! spin
                sig_matel(2)=omega_fac_matel&
                     *aimag(AA_h(m,n,beta)*KKspn_h(n,m,alpha,gamma)&
                           +AA_h(m,n,alpha)*KKspn_h(n,m,beta,gamma))
                !
                ! "Energy term", Eq. (31) MS10 (purely orbital). Integrand,
                ! including factor of two in front.
                !
                sig_en=omega_fac_en*(del_eig(m,gamma)+del_eig(n,gamma))&
                     *aimag(cmplx_i*AA_h(m,n,alpha)*AA_h(n,m,beta))
                !
             else
                !
                ! -------------------------------------------
                ! Re[sigma_{A,alpha beta gamma}]/(hbar.omega)
                ! -------------------------------------------
                !
                omega_fac_matel=2.0_dp/(omega_mn**2-freq**2)
                omega_fac_en=(3.0_dp*omega_mn**2-freq**2)&
                     /(omega_mn**2-freq**2)**2
                !
                ! "Matrix element term", integrand of Eq. (32) MS10 
                ! divided by omega and including factor of two in front 
                !
                ! orbital
                sig_matel(1)=omega_fac_matel&
                     *aimag(cmplx_i*(AA_h(m,n,beta)*KKorb_h(n,m,alpha,gamma)&
                                    -AA_h(m,n,alpha)*KKorb_h(n,m,beta,gamma)))
                ! spin
                sig_matel(2)=omega_fac_matel&
                     *aimag(cmplx_i*(AA_h(m,n,beta)*KKspn_h(n,m,alpha,gamma)&
                                    -AA_h(m,n,alpha)*KKspn_h(n,m,beta,gamma)))
                !
                ! "Energy term", Eq. (33) MS10 (purely orbital) divided by
                ! omega. Note: Changed the sign in front from "-" to "+", and
                ! exchanged the order of the band indices inside Im(...)
                !
                sig_en=omega_fac_en*(del_eig(m,gamma)+del_eig(n,gamma))&
                     *aimag(AA_h(m,n,alpha)*AA_h(n,m,beta))
             endif !T_odd
             sig_abc_k(loop_f,1)=sig_abc_k(loop_f,1)+sig_matel(1) !orb (m)-term
             sig_abc_k(loop_f,2)=sig_abc_k(loop_f,2)+sig_en !orb (e)-term
             sig_abc_k(loop_f,3)=sig_abc_k(loop_f,3)+sig_matel(2) !spin
          enddo !loop_f
       enddo !m
    enddo !n

    ! Now recalculate the dc (static) value sig_abc_k(omega=0), 
    ! including only transitions *below* a given energy, scanned 
    ! between zero and ecut_spectralsum
    !
    sig_abc_cut_k=0.0_dp
    do n=1,num_wann
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann
          if(eig(m)<fermi_energy.or.eig(m)>eig_cut) cycle
          omega_mn=eig(m)-eig(n)
          do loop_f=nfreq_cut,1,-1
             !---------------------------------------!
             if(omega_mn>(loop_f-1)*d_freq_cut) exit !
             !---------------------------------------!
             if(T_odd) then
                omega_fac_matel=2.0_dp/omega_mn
                omega_fac_en=2.0_dp/omega_mn
                sig_matel(1)=omega_fac_matel&
                     *aimag(AA_h(m,n,beta)*KKorb_h(n,m,alpha,gamma)&
                           +AA_h(m,n,alpha)*KKorb_h(n,m,beta,gamma))
                sig_matel(2)=omega_fac_matel&
                     *aimag(AA_h(m,n,beta)*KKspn_h(n,m,alpha,gamma)&
                           +AA_h(m,n,alpha)*KKspn_h(n,m,beta,gamma))
                sig_en=omega_fac_en*(del_eig(m,gamma)+del_eig(n,gamma))&
                     *aimag(cmplx_i*AA_h(m,n,alpha)*AA_h(n,m,beta))
             else
                omega_fac_matel=2.0_dp/omega_mn**2
                omega_fac_en=3.0_dp/omega_mn**2
                sig_matel(1)=omega_fac_matel&
                     *aimag(cmplx_i*(AA_h(m,n,beta)*KKorb_h(n,m,alpha,gamma)&
                                    -AA_h(m,n,alpha)*KKorb_h(n,m,beta,gamma)))
                sig_matel(2)=omega_fac_matel&
                     *aimag(cmplx_i*(AA_h(m,n,beta)*KKspn_h(n,m,alpha,gamma)&
                                    -AA_h(m,n,alpha)*KKspn_h(n,m,beta,gamma)))
                sig_en=omega_fac_en*(del_eig(m,gamma)+del_eig(n,gamma))&
                     *aimag(AA_h(m,n,alpha)*AA_h(n,m,beta))
             endif !T_odd
             sig_abc_cut_k(loop_f,1)=sig_abc_cut_k(loop_f,1)+sig_matel(1)
             sig_abc_cut_k(loop_f,2)=sig_abc_cut_k(loop_f,2)+sig_en
             sig_abc_cut_k(loop_f,3)=sig_abc_cut_k(loop_f,3)+sig_matel(2)
          enddo !loop_f
       enddo !m
    enddo !n
    
  end subroutine get_sig_abc_k


  subroutine get_alpha_k(kpt,alpha_k,imgamma_k)
  !====================================================!
  !                                                    !
  ! Traceless optical magnetoelectric tensor alpha_ij  !
  ! and totally symmetric quadrupolar tensor gamma_ijl !
  !                                                    ! 
  !====================================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i,bohr,eV_au
    use w90_utility, only       : utility_diagonalize,utility_rotate
    use w90_parameters, only    : num_wann,fermi_energy,optics_min_energy,&
                                  optics_num_points,have_disentangled,&
                                  dis_froz_max,sigma_abc_onlyorb
    use w90_wanint_common, only : fourier_R_to_k
    use w90_wan_ham, only   : get_D_h_a,get_deleig_a
    use w90_get_oper, only      : HH_R,AA_R,BB_R,CC_R,FF_R,SS_R

    ! Arguments
    !
    real(kind=dp),                     intent(in)  :: kpt(3)
    real(kind=dp), dimension(:,:,:,:), intent(out) :: alpha_k
    real(kind=dp), dimension(:,:,:),   intent(out) :: imgamma_k

    real(kind=dp),    allocatable :: imsig_s_k(:,:,:,:,:)
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: HH_h(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: AA_bar(:,:,:)
    complex(kind=dp), allocatable :: BB_bar(:,:,:)
    complex(kind=dp), allocatable :: CC_bar(:,:,:,:)
    complex(kind=dp), allocatable :: FF_bar(:,:,:,:)
    complex(kind=dp), allocatable :: AA_h(:,:,:)
    complex(kind=dp), allocatable :: BB_h(:,:,:)
    complex(kind=dp), allocatable :: CC_h(:,:,:,:)
    complex(kind=dp), allocatable :: FF_h(:,:,:,:)
    complex(kind=dp), allocatable :: SS_h(:,:,:)
    complex(kind=dp), allocatable :: KKspn_h(:,:,:,:)
    complex(kind=dp), allocatable :: KKorb_h(:,:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)

    integer                       :: i,j,n,m,loop_f,a,b,c,d,p
    real(kind=dp)                 :: eig(num_wann),del_eig(num_wann,3),&
                                     eig_cut,occ(num_wann),sig_matel(2),&
                                     sig_en,freq,omega_mn,eps_abc,&
                                     omega_fac_matel,omega_fac_en
    complex(kind=dp)              :: cdum

    allocate(imsig_s_k(3,nfreq,3,3,3)) ! MOVE TO ARGUMENT LIST?
    allocate(HH(num_wann,num_wann))
    allocate(HH_h(num_wann,num_wann))
    allocate(D_h(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(AA_bar(num_wann,num_wann,3))
    allocate(BB_bar(num_wann,num_wann,3))
    allocate(CC_bar(num_wann,num_wann,3,3))
    allocate(FF_bar(num_wann,num_wann,3,3))
    allocate(AA_h(num_wann,num_wann,3))
    allocate(BB_h(num_wann,num_wann,3))
    allocate(CC_h(num_wann,num_wann,3,3))
    allocate(FF_h(num_wann,num_wann,3,3))
    if(.not.(sigma_abc_onlyorb)) allocate(SS_h(num_wann,num_wann,3))
    allocate(KKspn_h(num_wann,num_wann,3,3))
    allocate(KKorb_h(num_wann,num_wann,3,3))
    allocate(mdum(num_wann,num_wann))

    call fourier_R_to_k(kpt,HH_R,HH,0)
    call utility_diagonalize(HH,num_wann,eig,UU)
    do i=1,3
       ! Matrix objects with one Cartesian index
       call fourier_R_to_k(kpt,HH_R,mdum,i) ! del_i HH
       call get_D_h_a(mdum,UU,eig,D_h(:,:,i))
       call get_deleig_a(del_eig(:,i),eig,mdum,UU)
       call fourier_R_to_k(kpt,AA_R(:,:,:,i),mdum,0)
       AA_bar(:,:,i)=utility_rotate(mdum,UU,num_wann)
       call fourier_R_to_k(kpt,BB_R(:,:,:,i),mdum,0)
       BB_bar(:,:,i)=utility_rotate(mdum,UU,num_wann)
       if(.not.(sigma_abc_onlyorb)) then
          call fourier_R_to_k(kpt,SS_R(:,:,:,i),mdum,0)
          ! SS_h = SS_bar (SS is "gauge-covariant" in the sense of WYSV06)
          SS_h(:,:,i)=utility_rotate(mdum,UU,num_wann)
       endif
       do j=1,3
          ! Matrix objects with two Cartesian indices
          call fourier_R_to_k(kpt,CC_R(:,:,:,j,i),mdum,0)
          CC_bar(:,:,j,i)=utility_rotate(mdum,UU,num_wann)
          call fourier_R_to_k(kpt,FF_R(:,:,:,j,i),mdum,0)
          FF_bar(:,:,j,i)=utility_rotate(mdum,UU,num_wann)
       enddo
    enddo
    
    ! Hamiltonian-gauge matrices
    !
    HH_h=cmplx_0
    AA_h=cmplx_0
    CC_h=cmplx_0
    FF_h=cmplx_0
    do i=1,3
       HH_h(i,i)=eig(i)
       ! Eq.(25) WYSV06
       AA_h(:,:,i)=AA_bar(:,:,i)+cmplx_i*D_h(:,:,i)
    enddo
    do j=1,3
       do i=1,j
          ! See 2012-02-06 notes, p.4
          FF_h(:,:,i,j)=FF_bar(:,:,i,j)&
            +cmplx_i*matmul(AA_bar(:,:,i),D_h(:,:,j))&
            +cmplx_i*matmul(D_h(:,:,i),AA_bar(:,:,j))&
            -matmul(D_h(:,:,i),D_h(:,:,j))
          FF_h(:,:,j,i)=conjg(transpose(FF_h(:,:,i,j)))
          CC_h(:,:,i,j)=CC_bar(:,:,i,j)&
            +cmplx_i*matmul(conjg(transpose(BB_bar(:,:,i))),D_h(:,:,j))&
            +cmplx_i*matmul(D_h(:,:,i),BB_bar(:,:,j))&
            -matmul(D_h(:,:,i),matmul(HH_h,D_h(:,:,j)))
          CC_h(:,:,j,i)=conjg(transpose(CC_h(:,:,i,j)))
       enddo
    enddo

    ! Evaluate in the H-gauge the matrices K_ij=hbar.B_ij, 
    ! where B_ij is defined in Eqs. (35-37) MS10.
    !
    ! ------------
    ! Orbital part
    ! ------------
    !
    do j=1,3
       do i=1,3
          do m=1,num_wann
             do n=1,m
                cdum=-cmplx_i*(del_eig(n,i)+del_eig(m,i))*AA_h(n,m,j)&
                     +CC_h(n,m,j,i)-CC_h(n,m,i,j)&
                     +eig(n)*FF_h(n,m,i,j)-eig(m)*FF_h(n,m,j,i)
                KKorb_h(n,m,i,j)=cdum/2.0_dp
                KKorb_h(m,n,i,j)=-conjg(KKorb_h(n,m,i,j))
             enddo !n
          enddo !m
       enddo  !i
    enddo !j
    !
    ! ---------
    ! Spin part
    ! ---------
    !
    KKspn_h=cmplx_0 
    !Diagonal elements in Cartesian indices will remain zero matrices
    if(.not.sigma_abc_onlyorb) then
       KKspn_h(:,:,1,2)= SS_h(:,:,3)
       KKspn_h(:,:,2,1)=-SS_h(:,:,3)
       KKspn_h(:,:,1,3)=-SS_h(:,:,2)
       KKspn_h(:,:,3,1)= SS_h(:,:,2)
       KKspn_h(:,:,2,3)= SS_h(:,:,1)
       KKspn_h(:,:,3,2)=-SS_h(:,:,1)
       ! Multiply by (i/2)(hbar^2)/m_e. In atomic units hbar=m_e=1, and we
       ! would get i/2. But our unit of energy is eV, not Hartee, and the 
       ! unit of length is Angstrom, not Bohr radius. Since hbar^2/m_e has 
       ! units of [E][L]^2, need to multiply by 
       ! 1 Ha = 27.21138387 eV and by (1 bohr)^2 = 0.529^2 Ang^2
       KKspn_h=KKspn_h*(cmplx_i/2.0_dp)*bohr**2/eV_au
    endif
    
    if(have_disentangled) then
       eig_cut=dis_froz_max+0.001_dp !some small number 
    else
       eig_cut=maxval(eig)+0.001_dp
    endif
    
    ! -----------------
    ! Im[sigma_{S,abc}] (symmetric under a <--> b)
    ! -----------------
    !
    imsig_s_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann !conduction, inside frozen window
          if(eig(m)<fermi_energy.or.eig(m)>eig_cut) cycle
          omega_mn=eig(m)-eig(n)
          do c=1,3
          do b=1,3
          do a=1,b ! Will symmetrize at the end
             do loop_f=1,nfreq
                freq=optics_min_energy+(loop_f-1)*d_freq
                omega_fac_matel=2.0_dp*omega_mn/(omega_mn**2-freq**2)
                omega_fac_en=2.0_dp*omega_mn**3/(omega_mn**2-freq**2)**2
                !
                ! "Matrix element term", Eq. (30) MS10. Integrand,
                ! including factor of two in front
                !
                ! orbital
                sig_matel(1)=omega_fac_matel&
                     *aimag(AA_h(m,n,b)*KKorb_h(n,m,a,c)&
                           +AA_h(m,n,a)*KKorb_h(n,m,b,c))
                ! spin
                sig_matel(2)=omega_fac_matel&
                     *aimag(AA_h(m,n,b)*KKspn_h(n,m,a,c)&
                           +AA_h(m,n,a)*KKspn_h(n,m,b,c))
                !
                ! "Energy term", Eq. (31) MS10 (purely orbital). Integrand,
                ! including factor of two in front
                !
                sig_en=omega_fac_en*(del_eig(m,c)+del_eig(n,c))&
                     *aimag(cmplx_i*AA_h(m,n,a)*AA_h(n,m,b))
                !
                imsig_s_k(1,loop_f,a,b,c)=imsig_s_k(1,loop_f,a,b,c)&
                     +sig_matel(1) !orb (m)-term
                imsig_s_k(2,loop_f,a,b,c)=imsig_s_k(2,loop_f,a,b,c)&
                     +sig_en !orb (e)-term
                imsig_s_k(3,loop_f,a,b,c)=imsig_s_k(3,loop_f,a,b,c)&
                     +sig_matel(2) !spin
             enddo !loop_f
          enddo !a
          enddo !b
          enddo !c
       enddo !m
    enddo !n
    
    ! ***CHECK*** Does this make any difference in rendering
    ! gamma_abc^(spin) *exactly* zero? If not, revert to previous code
    !
    ! Enforce symmetry under a <--> b
    !
    do b=1,3
       do a=1,b
          imsig_s_k(:,:,b,a,:)=imsig_s_k(:,:,a,b,:)
       enddo
    enddo

    ! Eq. (26) MS10 at transparent frequencies:
    ! alpha_da=(1/3).eps_abc.Im[sig_{dbc}],
    !
    alpha_k=0.0_dp
    do a=1,3
       do d=1,3
          do b=1,3
             do c=1,3
                eps_abc=real((b-a)*(c-a)*(c-b)/2,kind=dp) ! Levi-Civita symbol
                do loop_f=1,nfreq
                   do i=1,3 ! orb (m)-term, orb (e)-term, and spin
                      alpha_k(i,loop_f,d,a)=alpha_k(i,loop_f,d,a)&
                           +eps_abc*imsig_s_k(i,loop_f,d,b,c)/3.0_dp
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    ! Eq.(27) MS10 at non-absorbing frequencies, but *without* factor
    ! of 1/omega (more convenient for comparing with experiment)
    !
    imgamma_k=0.0_dp
    do p=1,10
       do loop_f=1,nfreq
          !***************************
          ! eventually change to i=1,2
          !***************************
          do i=1,3
             imgamma_k(i,loop_f,p)=imsig_s_k(i,loop_f,aa(p),bb(p),cc(p))&
                                  +imsig_s_k(i,loop_f,cc(p),aa(p),bb(p))&
                                  +imsig_s_k(i,loop_f,bb(p),cc(p),aa(p))
          enddo
       enddo
    enddo
    imgamma_k=imgamma_k/3.0_dp

  end subroutine get_alpha_k


  subroutine get_alphaspn_k(kpt,alphaspn_k,alphaspn_cut_k)
  !===========================================================!
  !                                                           !
  ! Spin-electronic static magnetoelectric tensor alphaspn_ij !
  !                                                           !
  !===========================================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i,bohr,eV_au
    use w90_utility, only       : utility_diagonalize,utility_rotate
    use w90_parameters, only    : num_wann,fermi_energy,cell_volume,&
                                  have_disentangled,dis_froz_max
    use w90_wanint_common, only : fourier_R_to_k
    use w90_wan_ham, only   : get_D_h_a
    use w90_get_oper, only      : HH_R,AA_R,SS_R

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), intent(out)                   :: alphaspn_k(3,3)
    real(kind=dp), dimension(:,:,:), intent(out) :: alphaspn_cut_k

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h_i(:,:)
    complex(kind=dp), allocatable :: AA_h(:,:,:)
    complex(kind=dp), allocatable :: SS_h(:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)

    integer                       :: i,j,n,m,loop_f
    real(kind=dp)                 :: eig(num_wann),eig_cut,fac
    complex(kind=dp)              :: cdum

    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
    allocate(D_h_i(num_wann,num_wann))
    allocate(AA_h(num_wann,num_wann,3))
    allocate(SS_h(num_wann,num_wann,3))
    allocate(mdum(num_wann,num_wann))

    call fourier_R_to_k(kpt,HH_R,HH,0)
    call utility_diagonalize(HH,num_wann,eig,UU)
    do i=1,3
       call fourier_R_to_k(kpt,HH_R,mdum,i) ! del_i HH
       call get_D_h_a(mdum,UU,eig,D_h_i(:,:))
       call fourier_R_to_k(kpt,AA_R(:,:,:,i),mdum(:,:),0)
       AA_h(:,:,i)=utility_rotate(mdum,UU,num_wann)+cmplx_i*D_h_i(:,:)
       call fourier_R_to_k(kpt,SS_R(:,:,:,i),mdum(:,:),0)
       SS_h(:,:,i)=utility_rotate(mdum,UU,num_wann)
    enddo

    if(have_disentangled) then
       eig_cut=dis_froz_max+0.001_dp !some small number 
    else
       eig_cut=maxval(eig)+0.001_dp
    endif

    ! hbar^2/(m_e.V_cell) after converting hbar^2/m_e from a.u. to eV.Angstrom^2
    !
    fac=bohr**2/eV_au/cell_volume ! units of ev/Angstrom

    ! Compute alphaspn_k (dimensionless due to multiplication by fac)
    !
    alphaspn_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann !conduction, inside frozen window
          if(eig(m)<fermi_energy.or.eig(m)>eig_cut) cycle
          do j=1,3
             do i=1,3
                alphaspn_k(i,j)=alphaspn_k(i,j)&
                     +fac*aimag(cmplx_i*SS_h(n,m,j)*AA_h(m,n,i))/(eig(m)-eig(n))
             enddo
          enddo
       enddo
    enddo
     
    ! Now recalculate alphaspn_k including only transitions below a given
    ! energy, scanned between zero and ecut_spectralsum
    !
    alphaspn_cut_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann !conduction
          ! Exclude transitions to states outside the frozen window
          if(eig(m)<fermi_energy.or.eig(m)>eig_cut) cycle
          do loop_f=nfreq_cut,1,-1
             !----------------------------------------------!
             if((eig(m)-eig(n))>(loop_f-1)*d_freq_cut) exit !
             !----------------------------------------------!
             do j=1,3
                do i=1,3
                   alphaspn_cut_k(loop_f,i,j)=alphaspn_cut_k(loop_f,i,j)&
                     +fac*aimag(cmplx_i*SS_h(n,m,j)*AA_h(m,n,i))/(eig(m)-eig(n))
                enddo
             enddo
          enddo !loop_f
       enddo !m
    enddo !n

  end subroutine get_alphaspn_k

end module w90_berry_wanint
