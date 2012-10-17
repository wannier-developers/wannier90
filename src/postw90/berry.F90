!-*- mode: F90; mode: font-lock -*-!

! ---------------------------------------------------------------
! REFERENCES
!
!   WYSV06 = PRB 74, 195118 (2006)  (Anomalous Hall conductivity)
!   YWVS07 = PRB 75, 195121 (2007)  (Interband optical conductivity)
!   LVTS12 = PRB 85, 014435 (2012)  (Orbital magnetization)
!   CTVR06 = PRB 74, 024408 (2006)  (  "          "       )
!
!   Note: AHC is now coded using the trace formulation of LVTS12
!         instead of the original formulation of WYSV06
! ---------------------------------------------------------------
!
! * Implementation in progress (undocumented, largely untested):
!
!   MV10   = PRB 82, 245118 (2010)  (Spatial dispersion)
!     -- Add logical input keyword 'optics_spatial_dispersion' (default F)
!        If 'T', calculate sig_abc, and the conversion into
!        optical rotatory power (if time-even), or into
!        magnetoelectric + electric-quadrupole response tensors (if time-odd)
!
! * Also undocumented (needs re-testing): reading k-points and weights from file

module w90_berry

  use w90_constants, only : dp

  implicit none

  private

  public :: berry_main,get_imf_k,get_img_k,get_imh_k
  real(kind=dp), parameter :: eps=1.0e-7
  !
  ! Pseudovector <--> antisymmetric tensor
  !
  ! x <--> (y,z)
  ! y <--> (z,x)
  ! z <--> (x,y)
  !
  integer, dimension(3), parameter :: alp=(/ 2,3,1 /)
  integer, dimension(3), parameter :: bet=(/ 3,1,2 /)

  ! Independent components of a symmetric tensor
  !
  ! 1 <--> xx
  ! 2 <--> yy
  ! 3 <--> zz
  ! 4 <--> xy
  ! 5 <--> xz
  ! 6 <--> yz
  !
  integer, dimension(6), parameter :: alps=(/ 1,2,3,1,1,2 /)
  integer, dimension(6), parameter :: bets=(/ 1,2,3,2,3,3 /)

  integer       :: nfreq,nfreq_cut,aa(10),bb(10),cc(10)
  real(kind=dp) :: d_freq,d_freq_cut
  logical       :: T_odd

  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine berry_main
  !=============================================================!
  !                                                             !
  ! Computes the following quantities:                          !
  !                                                             !
  !   (i) Intrinsic anomalous Hall conductivity                 ! 
  !  (ii) Interband optical conductivity in the long-wavelength ! 
  !       limit (order q^0)                                     !
  ! (iii) Joint density of states                               !
  !  (iv) Orbital magnetization                                 !
  !                                                             !
  !   -------------- WORK IN PROGRESS ---------------           !
  !   (v) Interband optical conductivity at order q^1           !
  !                                                             !
  !=============================================================!

    use w90_constants, only     : dp,pi,cmplx_i,elem_charge_SI,hbar_SI,&
                                  elec_mass_SI,eV_au,bohr,bohr_magn_SI,&
                                  eps0_SI,speedlight_SI
    use w90_comms, only         : on_root,num_nodes,my_node_id,comms_reduce
    use w90_io, only            : io_error,stdout,io_file_unit,seedname,&
                                  io_stopwatch
    use w90_postw90_common, only : nrpts,irvec,num_int_kpts_on_node,int_kpts,&
                                  adkpt,weight,rpt_origin
    use w90_parameters, only    : timing_level,num_wann,berry_interp_mesh,&
                                  !---------remove eventually---------
                                  alpha,beta,gamma,sigma_abc_onlyorb,&
                                  !-----------------------------------
                                  berry_adpt_mesh,berry_adpt_thresh,&
                                  wanint_kpoint_file,cell_volume,transl_inv,&
                                  berry_task,optics_time_parity,&
                                  optics_energy_min,optics_energy_max,&
                                  optics_energy_step,spn_decomp,&
                                  found_fermi_energy,fermi_energy
    use w90_get_oper, only      : get_HH_R,get_AA_R,get_BB_R,get_CC_R,&
                                  get_FF_R,get_SS_R

    ! AHC and orbital magnetization
    !
    ! First index labels J0,J1,J2 terms (in the notation of LVTS12) 
    ! Second index is the Cartesian component 
    !
    real(kind=dp) :: imf_k(3,3),imf(3,3)
    real(kind=dp) :: img_k(3,3),img(3,3)
    real(kind=dp) :: imh_k(3,3),imh(3,3)
    real(kind=dp) :: ahc(3,3)
    real(kind=dp) :: LCtil(3,3),ICtil(3,3),Morb(3,3) 

    ! Joint density of states
    !
    real(kind=dp), allocatable :: jdos_k(:,:)
    real(kind=dp), allocatable :: jdos(:,:)
    
    ! Optical conductivity (long wavelength limit)
    !
    integer                    :: ncomp ! number of components (3 or 6)
    real(kind=dp), allocatable :: sig_k(:,:,:,:)
    real(kind=dp), allocatable :: sig(:,:,:,:)
    real(kind=dp), allocatable :: ahc_kk(:,:)

    ! Optical conductivity (spatial dispersion)
    !
    ! sig_abc_(:,1) is the "matrix element term" part of the orbital contrib
    ! sig_abc_(:,2) is the "energy term" part of the orbital contribution
    ! sig_abc_(:,3) is the spin contribution ("matrix element term"-only)
    !
    real(kind=dp), allocatable :: sig_abc_k(:,:)
    real(kind=dp), allocatable :: sig_abc(:,:)

    ! Same as above, but contrib to the static value from optical transitions 
    ! below a given frequency, scanned between zero and optics_energy_max
    !
    real(kind=dp), allocatable :: sig_abc_cut_k(:,:)
    real(kind=dp), allocatable :: sig_abc_cut(:,:)

    ! Traceless optical magnetoelectric (ME) tensor alpha_ij
    ! and totally symmetric electric quadrupole (EQ) response tensor gamma_ijl 
    !
    real(kind=dp), allocatable :: alphaME_k(:,:,:,:)
    real(kind=dp), allocatable :: alphaME(:,:,:,:)
    real(kind=dp), allocatable :: imgamma_k(:,:,:)
    real(kind=dp), allocatable :: imgamma(:,:,:)
    
    ! Static spin-electronic magnetoelectric tensor alphaspn_ij, and 
    ! contribution from optical transitions below a given frequency, scanned
    ! between zero and optics_energy_max
    !
    real(kind=dp) :: alphaspn_k(3,3)
    real(kind=dp) :: alphaspn(3,3)
    !
    real(kind=dp), allocatable :: alphaspn_cut_k(:,:,:)
    real(kind=dp), allocatable :: alphaspn_cut(:,:,:)

    real(kind=dp)     :: kweight,kweight_adpt,kpt(3),kpt_ad(3),freq,&
                         adpt_trigger,fac,vdum(3)
    integer           :: i,j,p,loop_x,loop_y,loop_z,loop_tot,loop_adpt,&
                         ifreq,adpt_counter,ndim,jdos_unit,file_unit,&
                         tot_unit,alpha_unit(3,3),gamma_unit(10)
    character(len=20) :: file_name
    logical           :: eval_ahc,eval_morb,eval_sig_ab,eval_sig_abc,&
                         eval_ME_EQ,eval_MEspn

    if(.not.found_fermi_energy) call io_error&
         (&
 'Must set either "fermi_energy" or "num_valence_bands" for optical properties'&
         )

    if (timing_level>1.and.on_root)&
         call io_stopwatch('berry: prelims',1)

    nfreq=max(nint((optics_energy_max-optics_energy_min)/optics_energy_step),2)
    d_freq=(optics_energy_max-optics_energy_min)/(nfreq-1)

    nfreq_cut=max(nint(optics_energy_max/optics_energy_step),2)
    d_freq_cut=optics_energy_max/(nfreq_cut-1)
    
    ! Must initialize to .false. all eval_ flags
    !
    eval_ahc=.false.
    eval_morb=.false.
    eval_sig_ab=.false.
    eval_sig_abc=.false.
    eval_ME_EQ=.false.
    eval_MEspn=.false.
    T_odd=.false.
    if(index(berry_task,'ahc')>0) then
       eval_ahc=.true.
    elseif(index(berry_task,'morb')>0) then
       eval_morb=.true.
    end if
    if(index(berry_task,'optics')>0) then
       eval_sig_ab=.true.
       if(index(optics_time_parity,'odd')>0) T_odd=.true.
    elseif(index(berry_task,'gyro')>0) then
       T_odd=.true.
       eval_sig_abc=.true.
       eval_ME_EQ=.true.
    elseif(index(berry_task,'noa')>0) then
       T_odd=.false.
       eval_sig_abc=.true.
    end if
    if(index(berry_task,'mespn')>0) eval_MEspn=.true.

    ! NOTE: This may change, if eval_ME_EQ is done separately from eval_sig_abc
!--------------------------------------------
!    if(alpha==0.or.beta==0) call io_error&
!('Must specify cartesian directions alpha and beta for optical properties')
!--------------------------------------------
    if(gamma==0.and.eval_sig_abc) call io_error&
         ('Must specify cartesian direction gamma for berry_task=gyro,noa')

    ! ---------------------------------------
    ! TO DO: Calculate here the adaptive grid (?)
    ! ---------------------------------------

    ! Wannier matrix elements, allocations and initializations
    !
    adpt_counter=0
    if(eval_ahc) then
       call get_HH_R 
       call get_AA_R
       imf=0.0_dp
    endif
    if(eval_morb) then
       call get_HH_R 
       call get_AA_R
       call get_BB_R
       call get_CC_R
       imf=0.0_dp
       img=0.0_dp
       imh=0.0_dp
    endif
    if(eval_sig_ab) then
       call get_HH_R 
       call get_AA_R
       if(T_odd) then
          ncomp=3 ! antisymmetric conductivity, 3 indep. components
       else
          ncomp=6 ! symmetric conductivity, 6 indep. components
       endif
       if(spn_decomp) then
          call get_SS_R
          !
          ! Extra entries contain the decomposition into 
          ! up-->up, down-->down, and spin-flip transitions
          !
          ndim=4
       else
          ndim=1
       end if
       allocate(sig_k(3,ncomp,nfreq,ndim))
       allocate(sig(3,ncomp,nfreq,ndim)) 
       if(T_odd) allocate(ahc_kk(3,ndim))
       allocate(jdos_k(nfreq,ndim))
       allocate(jdos(nfreq,ndim))
       sig=0.0_dp
       jdos=0.0_dp
    elseif(eval_sig_abc) then
       call get_HH_R 
       call get_AA_R
       call get_BB_R
       call get_CC_R
       call get_FF_R
       if(.not.sigma_abc_onlyorb) call get_SS_R
       allocate(sig_abc_k(nfreq,3))
       allocate(sig_abc(nfreq,3))
       allocate(sig_abc_cut_k(nfreq_cut,3))
       allocate(sig_abc_cut(nfreq_cut,3))
       sig_abc=0.0_dp
       sig_abc_cut=0.0_dp
    endif
    if(eval_ME_EQ) then
       call get_HH_R 
       call get_AA_R
       call get_BB_R
       call get_CC_R
       call get_FF_R
       call get_SS_R
       allocate(alphaME_k(3,nfreq,3,3))
       allocate(alphaME(3,nfreq,3,3))
       !
       ! N.B.: There should be no spin contribution to gamma, so ultimately want
       !       to use the dimensions (2,nfreq,10). For the moment include spin
       !       to check whether its contribution vanishes numerically
       !
       allocate(imgamma_k(3,nfreq,10))
       allocate(imgamma(3,nfreq,10))
       alphaME=0.0_dp
       imgamma=0.0_dp
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
       call get_HH_R 
       call get_AA_R
       call get_SS_R
       allocate(alphaspn_cut_k(nfreq_cut,3,3))
       allocate(alphaspn_cut(nfreq_cut,3,3))
       alphaspn=0.0_dp
       alphaspn_cut=0.0_dp
    endif

    if(on_root) then

       write(stdout,'(/,/,1x,a)')&
            'Properties calculated in module  b e r r y'
       write(stdout,'(1x,a)')&
            '------------------------------------------'

       if(eval_ahc) write(stdout,'(/,3x,a)')&
            '* Anomalous Hall conductivity (AHC)'

       if(eval_morb) write(stdout,'(/,3x,a)')&
            '* Orbital magnetization'

       if(eval_sig_ab) then
          write(stdout,'(/,3x,a)') '* Joint density of states'
          write(stdout,'(6x,a)')&
               '-- File '//trim(trim(seedname))//'-jdos.dat'
          if(T_odd) then
             write(stdout,'(/,3x,a)')&
             '* Interband optical conductivity'
             write(stdout,'(6x,a)')&
                  '-- Imaginary, antisymmetric part (dichroic absorption)'
             write(stdout,'(6x,a)')&
                  '-- omega*sigma^A_{ij}(omega) in units of 10^29 sec^{-2}'
             write(stdout,'(6x,a)')&
                  '-- Files '//trim(trim(seedname))//'-sigA_ij.dat'
             write(stdout,'(/,3x,a)')&
        '* Cumulative anomalous Hall conductivity from Kramers-Kronig in S/cm '
             write(stdout,'(6x,a)')&
                  '-- Files '//trim(trim(seedname))//'-kk_ij.dat'
          else
             write(stdout,'(/,3x,a)')&
                  '* Dielectric function epsilon^r_{ij}(omega)'
             write(stdout,'(6x,a)')&
                  '-- Imaginary, symmetric part (ordinary absorption)'
             write(stdout,'(6x,a)')&
                  '-- Files '//trim(trim(seedname))//'-epsS_ij.dat'
             write(stdout,'(/,3x,a)')&
                  '* Interband optical conductivity epsilon'
             write(stdout,'(6x,a)')&
                  '-- Real, symmetric part (ordinary absorption)'
             write(stdout,'(6x,a)')&
                  '-- Sigma^S_{ij}(omega) in units of 10^15 sec^{-1}'
             write(stdout,'(6x,a)')&
                  '-- Files '//trim(trim(seedname))//'-sigS_ij.dat'
          end if

       end if !eval_sig_ab

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
                   trim(trim(seedname))//'-sigS_'//&
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

       if(transl_inv) then
          if(eval_morb.or.eval_sig_abc.or.eval_ME_EQ)&
            call io_error('transl_inv=T disabled for orb, sig_abc, eval_ME_EQ')
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

    if (timing_level>1.and.on_root) then
         call io_stopwatch('berry: prelims',2)
         call io_stopwatch('berry: kpts',1)
      endif
       
    ! Loop over interpolation k-points
    !
    if(wanint_kpoint_file) then
       
!       if(on_root) then
          !
!          write(stdout,'(/,1x,a)') 'Sampling the irreducible BZ only'
!          write(stdout,'(3x,a)')&
!          'WARNING: - IBZ implementation is currently limited to simple cases:'
!          write(stdout,'(3x,a)')&
!           '           Check results agains a full BZ calculation!'
!       end if

       ! Loop over k-points on the irreducible wedge of the Brillouin zone,
       ! read from file 'kpoint.dat'
       !
       do loop_tot=1,num_int_kpts_on_node(my_node_id)
          kpt(:)=int_kpts(:,loop_tot)
          kweight=weight(loop_tot)
          kweight_adpt=kweight/berry_adpt_mesh**3
          if(eval_ahc) then 
             call get_imf_k(kpt,imf_k)
             do i=1,3
                vdum(i)=sum(imf_k(1:3,i))
             enddo
             adpt_trigger=sqrt(dot_product(vdum,vdum))
          elseif(eval_morb) then
             call get_imf_k(kpt,imf_k)
             call get_img_k(kpt,img_k)
             call get_imh_k(kpt,imh_k)
             do i=1,3
                vdum(i)=sum(img_k(1:3,i))+sum(imh_k(1:3,i))&
                     -2.0_dp*fermi_energy*sum(imf_k(1:3,i))
             enddo
             adpt_trigger=sqrt(dot_product(vdum,vdum))
          else ! Ensure that adaptive refinement is not triggered
             adpt_trigger=berry_adpt_thresh-1.0_dp
          end if
          if(eval_sig_ab) then
               call get_sig_k(kpt,sig_k,jdos_k)
          elseif(eval_sig_abc) then
             call get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
          endif
          if(eval_ME_EQ) call get_ME_EQ_k(kpt,alphaME_k,imgamma_k)
          if(eval_MEspn) call get_MEspn_k(kpt,alphaspn_k,alphaspn_cut_k)
          !
          ! Decide whether or not to trigger adaptive refinement of 
          ! integration k-mesh, when computing AHC (and possibly MCD) or Morb
          !
          if(adpt_trigger>berry_adpt_thresh) then
             adpt_counter=adpt_counter+1
             do loop_adpt=1,berry_adpt_mesh**3
                kpt_ad(:)=kpt(:)+adkpt(:,loop_adpt)
                if(eval_ahc) then
                   call get_imf_k(kpt_ad,imf_k)
                   imf=imf+imf_k*kweight_adpt
                elseif(eval_morb) then
                   call get_imf_k(kpt_ad,imf_k)
                   imf=imf+imf_k*kweight_adpt
                   call get_img_k(kpt_ad,img_k)
                   img=img+img_k*kweight_adpt
                   call get_imh_k(kpt_ad,imh_k)
                   imh=imh+imh_k*kweight_adpt
                end if
                if(eval_sig_ab) then
                   call get_sig_k(kpt_ad,sig_k,jdos_k)
                   sig=sig+sig_k*kweight_adpt
                   jdos=jdos+jdos_k*kweight_adpt
                elseif(eval_sig_abc) then
                   call get_sig_abc_k(kpt_ad,sig_abc_k,sig_abc_cut_k)
                   sig_abc=sig_abc+sig_abc_k*kweight_adpt
                   sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight_adpt
                end if
                if(eval_ME_EQ) then
                   call get_ME_EQ_k(kpt_ad,alphaME_k,imgamma_k)
                   alphaME=alphaME+alphaME_k*kweight_adpt
                   imgamma=imgamma+imgamma_k*kweight_adpt
                end if
                if(eval_MEspn) then
                   call get_MEspn_k(kpt_ad,alphaspn_k,alphaspn_cut_k)
                   alphaspn=alphaspn+alphaspn_k*kweight_adpt
                   alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight_adpt
                end if
             end do
          else
             if(eval_ahc) then
                imf=imf+imf_k*kweight
             elseif(eval_morb) then 
                imf=imf+imf_k*kweight
                img=img+img_k*kweight
                imh=imh+imh_k*kweight
             end if
             if(eval_sig_ab) then
                sig=sig+sig_k*kweight
                jdos=jdos+jdos_k*kweight
             elseif(eval_sig_abc) then
                sig_abc=sig_abc+sig_abc_k*kweight
                sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight
             endif
             if(eval_ME_EQ) then
                alphaME=alphaME+alphaME_k*kweight
                imgamma=imgamma+imgamma_k*kweight
             endif
             if(eval_MEspn) then
                alphaspn=alphaspn+alphaspn_k*kweight
                alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight
             endif
          end if
       end do

    else ! Do not read 'kpoint.dat'. Loop over a uniform grid in the full BZ

!       if (on_root) write(stdout,'(/,1x,a)') 'Sampling the full BZ'
!       kweight=1.0_dp/optics_num_points**3
       kweight = 1.0_dp / real(PRODUCT(berry_interp_mesh),kind=dp)
       kweight_adpt=kweight/berry_adpt_mesh**3

       do loop_tot=my_node_id,PRODUCT(berry_interp_mesh)-1,num_nodes
          loop_x= loop_tot/(berry_interp_mesh(2)*berry_interp_mesh(3))
          loop_y=(loop_tot-loop_x*(berry_interp_mesh(2)&
               *berry_interp_mesh(3)))/berry_interp_mesh(3)
          loop_z=loop_tot-loop_x*(berry_interp_mesh(2)*berry_interp_mesh(3))&
                -loop_y*berry_interp_mesh(3)
          kpt(1)=real(loop_x,dp)/real(berry_interp_mesh(1),dp)
          kpt(2)=real(loop_y,dp)/real(berry_interp_mesh(2),dp)
          kpt(3)=real(loop_z,dp)/real(berry_interp_mesh(3),dp)
          if(eval_ahc) then
             call get_imf_k(kpt,imf_k)
             do i=1,3
                vdum(i)=sum(imf_k(1:3,i))
             enddo
             adpt_trigger=sqrt(dot_product(vdum,vdum))
          elseif(eval_morb) then
             call get_imf_k(kpt,imf_k)
             call get_img_k(kpt,img_k)
             call get_imh_k(kpt,imh_k)
             do i=1,3
                vdum(i)=sum(img_k(1:3,i))+sum(imh_k(1:3,i))&
                     -2.0_dp*fermi_energy*sum(imf_k(1:3,i))
             enddo
             adpt_trigger=sqrt(dot_product(vdum,vdum))
          else ! Ensure that adaptive refinement is not triggered
             adpt_trigger=berry_adpt_thresh-1.0_dp
          end if
          if(eval_sig_ab) then
             call get_sig_k(kpt,sig_k,jdos_k)
          elseif(eval_sig_abc) then
             call get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
          endif
          if(eval_ME_EQ) call get_ME_EQ_k(kpt,alphaME_k,imgamma_k)
          if(eval_MEspn) call get_MEspn_k(kpt,alphaspn_k,alphaspn_cut_k)
          if(adpt_trigger>berry_adpt_thresh) then
             adpt_counter=adpt_counter+1
             do loop_adpt=1,berry_adpt_mesh**3
                kpt_ad(:)=kpt(:)+adkpt(:,loop_adpt)
                 if(eval_ahc) then
                   call get_imf_k(kpt_ad,imf_k)
                   imf=imf+imf_k*kweight_adpt
                elseif(eval_morb) then
                   call get_imf_k(kpt_ad,imf_k)
                   imf=imf+imf_k*kweight_adpt
                   call get_img_k(kpt_ad,img_k)
                   img=img+img_k*kweight_adpt
                   call get_imh_k(kpt_ad,imh_k)
                   imh=imh+imh_k*kweight_adpt
                end if
                if(eval_sig_ab) then
                   call get_sig_k(kpt_ad,sig_k,jdos_k)
                   sig=sig+sig_k*kweight_adpt
                   jdos=jdos+jdos_k*kweight_adpt
                elseif(eval_sig_abc) then
                   call get_sig_abc_k(kpt_ad,sig_abc_k,sig_abc_cut_k)
                   sig_abc=sig_abc+sig_abc_k*kweight_adpt
                   sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight_adpt
                end if
                if(eval_ME_EQ) then
                   call get_ME_EQ_k(kpt_ad,alphaME_k,imgamma_k)
                   alphaME=alphaME+alphaME_k*kweight_adpt
                   imgamma=imgamma+imgamma_k*kweight_adpt
                end if
                if(eval_MEspn) then
                   call get_MEspn_k(kpt_ad,alphaspn_k,alphaspn_cut_k)
                   alphaspn=alphaspn+alphaspn_k*kweight_adpt
                   alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight_adpt
                end if
             end do
          else
             if(eval_ahc) then
                imf=imf+imf_k*kweight
             elseif(eval_morb) then 
                imf=imf+imf_k*kweight
                img=img+img_k*kweight
                imh=imh+imh_k*kweight
             end if
             if(eval_sig_ab) then
                sig=sig+sig_k*kweight
                jdos=jdos+jdos_k*kweight
             elseif(eval_sig_abc) then
                sig_abc=sig_abc+sig_abc_k*kweight
                sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight
             endif
             if(eval_ME_EQ) then
                alphaME=alphaME+alphaME_k*kweight
                imgamma=imgamma+imgamma_k*kweight
             endif
             if(eval_MEspn) then
                alphaspn=alphaspn+alphaspn_k*kweight
                alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight
             endif
          end if
       end do !loop_tot
       
    end if !wanint_kpoint_file

! Collect contributions from all nodes    
!
    if(eval_ahc) then
       call comms_reduce(imf(1,1),3*3,'SUM')
       call comms_reduce(adpt_counter,1,'SUM')
    elseif(eval_morb) then
       call comms_reduce(imf(1,1),3*3,'SUM')
       call comms_reduce(img(1,1),3*3,'SUM')
       call comms_reduce(imh(1,1),3*3,'SUM')
       call comms_reduce(adpt_counter,1,'SUM')
    end if
    if(eval_sig_ab) then
       call comms_reduce(sig(1,1,1,1),3*ncomp*nfreq*ndim,'SUM')
       call comms_reduce(jdos(1,1),nfreq*ndim,'SUM')
    endif
    if(eval_sig_abc) then
       call comms_reduce(sig_abc(1,1),3*nfreq,'SUM')
       call comms_reduce(sig_abc_cut(1,1),3*nfreq_cut,'SUM')
    endif
    if(eval_ME_EQ) then
       call comms_reduce(alphaME(1,1,1,1),3*nfreq*3*3,'SUM')
       !
       ! *****************************************
       ! Eventually change 3*nfreq*3 --> 2*nfreq*3
       ! *****************************************
       !
       call comms_reduce(imgamma(1,1,1),3*nfreq*10,'SUM')
    endif
    if(eval_MEspn) then
       call comms_reduce(alphaspn(1,1),3*3,'SUM')
       call comms_reduce(alphaspn_cut(1,1,1),nfreq_cut*3*3,'SUM')
    endif
    
    if(on_root) then

       write(stdout,'(1x,a)') ' '
       if(wanint_kpoint_file) then

          write(stdout,'(1x,a47,i10,a)')&
               'Nominal interpolation mesh in IBZ: ',&
               sum(num_int_kpts_on_node),' points'

          if(eval_ahc .or. eval_morb) then
             if(eval_ahc) then
                write(stdout,'(1x,a47,3x,f8.4,a)')&
                     'Adaptive refinement triggered when Omega(k) >',&
                     berry_adpt_thresh,' Ang^2'
             else
                write(stdout,'(1x,a48,f8.4,a)')&
                     'Adaptive refinement triggered when k-integrand >',&
                     berry_adpt_thresh,' eV.Ang^2'
             end if
             write(stdout,'(1x,a47,i7,a,f8.4,a)')&
                  'How many points triggered adaptive refinement: ',&
                  adpt_counter,' (',&
                  100*real(adpt_counter,dp)/sum(num_int_kpts_on_node),' %)'
             if(berry_adpt_mesh < 10) then
                write(stdout,'(1x,a47,4x,i1,a,i1,a,i1)')&
                     'Adaptive mesh: ',berry_adpt_mesh,'x',&
                     berry_adpt_mesh,'x',berry_adpt_mesh
             elseif(berry_adpt_mesh < 100) then
                write(stdout,'(1x,a47,i2,a,i1,a,i1)')&
                     'Adaptive mesh: ',berry_adpt_mesh,'x',&
                     berry_adpt_mesh,'x',berry_adpt_mesh
             else
                write(stdout,'(1x,a47,i3,a,i1,a,i1)')&
                     'Adaptive mesh: ',berry_adpt_mesh,'x',&
                     berry_adpt_mesh,'x',berry_adpt_mesh
             end if
             write(stdout,'(1x,a47,i10)') 'Total number of points: ',&
                  sum(num_int_kpts_on_node)-adpt_counter+&
                  adpt_counter*berry_adpt_mesh**3
          end if

       else

          write(stdout,*) ' '
          if((eval_ahc .or. eval_morb) .and. berry_adpt_mesh.ne.1) then
             write(stdout, '(1x,a30,3(i0,1x))')&
                  'Interpolation grid (nominal): ',&
                  berry_interp_mesh(1),berry_interp_mesh(2),berry_interp_mesh(3)
             write(stdout, '(1x,a30,i0,a)') 'Adaptive refinement grid: ',&
                  berry_adpt_mesh,'^3'
             write(stdout, '(1x,a30,f6.2,a)')&
                  'Points triggering refinement: ',&
                  100*real(adpt_counter,dp)/product(berry_interp_mesh),'%'
          else
             write(stdout, '(1x,a20,3(i0,1x))') 'Interpolation grid: ',&
                  berry_interp_mesh(1:3)
             !,berry_interp_mesh(2),berry_interp_mesh(3)
          endif

       end if !wanint_kpoint_file

       if(eval_ahc) then

          ! --------------------------------------------------------------------
          ! At this point imf contains 
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
          ! To get a conductivity in units of S/cm,
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
          ! ===========================
          ! fac = -e^2/(hbar.V_c*10^-8) 
          ! ===========================
          !
          ! with 'V_c' in Angstroms^3, and 'e' and 'hbar' in SI units.
          ! --------------------------------------------------------------------
          !
          fac=-1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
          ahc(:,:)=imf(:,:)*fac

          write(stdout,'(/,/,1x,a)')&
               'AHC (S/cm)       x          y          z'
          write(stdout,'(1x,a)')&
               '=========='
          write(stdout,'(1x,a9,2x,3(f10.4,1x))') 'J0 term :',&
               ahc(1,1),ahc(1,2),ahc(1,3)
          write(stdout,'(1x,a9,2x,3(f10.4,1x))') 'J1 term :',&
               ahc(2,1),ahc(2,2),ahc(2,3)
          write(stdout,'(1x,a9,2x,3(f10.4,1x))') 'J2 term :',&
               ahc(3,1),ahc(3,2),ahc(3,3)
          write(stdout,'(1x,a)')&
               '-------------------------------------------'
          write(stdout,'(1x,a9,2x,3(f10.4,1x),/)') 'Total   :',&
               sum(ahc(:,1)),sum(ahc(:,2)),sum(ahc(:,3))
          
       elseif(eval_morb) then

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
          fac=-eV_au/bohr**2
          LCtil(:,:)=(img(:,:)-fermi_energy*imf(:,:))*fac
          ICtil(:,:)=(imh(:,:)-fermi_energy*imf(:,:))*fac
          Morb=LCtil+ICtil

          write(stdout,'(/,/,1x,a)')&
               'M_orb (bohr magn/cell)        x          y          z'
          write(stdout,'(1x,a)')&
               '======================'
          write(stdout,'(1x,a22,2x,3(f10.4,1x))') 'Local circulation :',&
               sum(LCtil(1:3,1)),sum(LCtil(1:3,2)),sum(LCtil(1:3,3))
          write(stdout,'(1x,a22,2x,3(f10.4,1x))') 'Itinerant circulation:',&
               sum(ICtil(1:3,1)),sum(ICtil(1:3,2)),sum(ICtil(1:3,3))
          write(stdout,'(1x,a)')&
               '--------------------------------------------------------'
          write(stdout,'(1x,a22,2x,3(f10.4,1x),/)') 'Total   :',&
               sum(Morb(1:3,1)),sum(Morb(1:3,2)),sum(Morb(1:3,3))

          !DEBUG
          write(stdout,'(/,/,1x,a)') 'DEBUG: Breakdown of z-component'
          write(stdout,'(1x,a)') 'M_orb = LCtil + ICtil'
          write(stdout,'(/,29x,a4,12x,a4,12x,a4)') '(J0)','(J1)','(J2)'
          write(stdout,'(1x,a,1x,e13.6,3(1x,a,1x,e13.6))')&
            'LCtil=',sum(LCtil(:,3)),'=',LCtil(1,3),'+',LCtil(2,3),'+',LCtil(3,3)
          write(stdout,'(1x,a,1x,e13.6,3(1x,a,1x,e13.6))')&
            'ICtil=',sum(ICtil(:,3)),'=',ICtil(1,3),'+',ICtil(2,3),'+',ICtil(3,3)
          !ENDDEBUG

       end if ! eval_ahc or eval_morb

       ! ---------------------!
       ! Optical conductivity !
       ! ---------------------!
       !
       if(eval_sig_ab) then

          write(stdout,'(/,1x,a)')&
               '--------------------------------------------------'
          write(stdout,'(1x,a)')&
               'Output data files related to optical conductivity:'
          write(stdout,'(1x,a)')&
               '--------------------------------------------------'

          file_name=trim(seedname)//'-jdos.dat'
          write(stdout,'(/,3x,a)') '* '//file_name
          jdos_unit=io_file_unit()
          open(jdos_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
          do ifreq=1,nfreq
             freq=optics_energy_min+(ifreq-1)*d_freq
             write(jdos_unit,'(5E16.8)') freq,jdos(ifreq,:)
          enddo
          close(jdos_unit)

          ! --------------------------------------------------------------------
          ! At this point sig contains 
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
          ! sigma=(pi.e^2/hbar)(omega)int dk/(2.pi)^3 Optical_matrix_elem(k)dk,
          !  
          ! where the optical matrix element is (for MCD; replace "Im" by "Re" 
          ! for ordinary absorption)
          !
          ! sum_n^occ sum_m^empty Im(A_{nm,alpha}A_{mn,beta}).
          !                       .delta(omega-omega_{mn}),
          !
          ! which has units of [L]^2[t]
          ! (Berry connection has units of length, and delta(...) has units of
          ! time)
          !
          ! Hence need to multiply by pi*omega*e^2/(hbar.V_c) to get a 
          ! conductivity,
          !
          ! (i)   Divide by V_c to obtain 
          !       (1/N) sum_k Optical_matrix_elem(k)/V_c (units of 1/[L]) 
          ! (ii)  I am assuming the working unit of length is Angstroms. 
          !       Multiply by 10^8 to convert to centimeters
          ! (iii) Multiply by pi.e^2/hbar expressed in SI (e^2/hbar is the 
          !       quantum of conductance, which has SI units of Siemens) to 
          !       get the final result for the conductivity in units of S/cm
          !
          ! ================================================
          ! fac = pi*10^8*elem_charge_SI^2/(hbar_SI*V_c_ang) 
          ! ================================================
          !
          ! (iv)   Multiply by omega 
          ! --------------------------------------------------------------------
          !
          ! Steps (i-iii) combined (do (iv) later, after Kramers-Kronig):
          fac=pi*1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
          sig=sig*fac
       
          if(T_odd) then
             fac=9.0e-18_dp*elem_charge_SI/hbar_SI ! ~0.0137
             do i=1,ncomp
                ! Anomalous Hall conductivity as cumulative Kramers-Kronig 
                ! transform of the magnetic circular dichroism spectrum: 
                ! Fig. 5b & Eq. (43) of YWVS07 
                file_name= trim(seedname)//'-kk_'//&
                     achar(119+alp(i))//achar(119+bet(i))//'.dat'
                file_name=trim(file_name)
                file_unit=io_file_unit()
                write(stdout,'(/,3x,a)') '* '//file_name
                open(file_unit,FILE=file_name,STATUS='UNKNOWN',&
                     FORM='FORMATTED')   
                ahc_kk=0.0_dp
                do ifreq=nfreq,1,-1
                   ahc_kk(1,:)=ahc_kk(1,:)+(2.0_dp/pi)*sig(1,i,ifreq,:)*d_freq
                   ahc_kk(2,:)=ahc_kk(2,:)+(2.0_dp/pi)*sig(2,i,ifreq,:)*d_freq
                   ahc_kk(3,:)=ahc_kk(3,:)+(2.0_dp/pi)*sig(3,i,ifreq,:)*d_freq
                   freq=optics_energy_min+(ifreq-1)*d_freq
                   !
                   ! Since the factor (1/freq)*d_freq in Eq.(43) YWVS07 is 
                   ! dimensionless, the conversion factor to obtain the AHC in 
                   ! S/cm is the same as for sig
                   !
                   write(file_unit,'(5E16.8)') freq,&
                        ahc_kk(1,:)+ahc_kk(2,:)+ahc_kk(3,:)
                end do
                close(file_unit)
                file_name= trim(seedname)//'-sigA_'//&
                     achar(119+alp(i))//achar(119+bet(i))//'.dat'
                file_name=trim(file_name)
                file_unit=io_file_unit()
                write(stdout,'(/,3x,a)') '* '//file_name
                open(file_unit,FILE=file_name,STATUS='UNKNOWN',&
                     FORM='FORMATTED')   
                do ifreq=1,nfreq
                   freq=optics_energy_min+(ifreq-1)*d_freq
                   !
                   ! After Kramers-Kronig can multiply by omega (step (iv) 
                   ! above). To write to file sigma^A_{ij}(omega) in S/cm, use
                   ! -----------------------------------------------------------
                   ! write(sig_unit,'(5E16.8)') freq,&
                   !   freq*(sig(1,i,ifreq,:)+sig(2,i,ifreq,:)+sig(3,i,ifreq,:))
                   ! -----------------------------------------------------------
                   !
                   ! To convert a conductivity from S/cm to 1/sec (Gaussian 
                   ! units), multiply by 9E11. For the MCD spectrum it is 
                   ! customary to plot omega*sigma(omega), in units of 
                   ! 10^29 sec^{-2}. Using omega= E(ev)(e/hbar) sec^{-1}, we 
                   ! find the conversion factor given above, 
                   ! fac=9.0e-18_dp*elem_charge_SI/hbar_SI
                   !
                   write(file_unit,'(5E16.8)') freq,fac*freq**2&
                        *(sig(1,i,ifreq,:)+sig(2,i,ifreq,:)+sig(3,i,ifreq,:))
                enddo
                close(file_unit)
             enddo !ncomp
          else ! time-even
             !
             ! epsilon^r = 1+ 4.pi.i.sigma/omega
             !
             ! ------------------------------------
             ! Im[epsilon^r] = 4.pi.Re[sigma]/omega
             ! ------------------------------------
             !
             ! sig contains Re[sigma] in S/cm divided by energy in eV. Want 
             ! Re[sigma] in 1/sec divided by omega in 1/sec (dimensionless),
             ! then converted to Im[epsilon^r]
             !
!             fac=4.0_dp*pi*9.0e11_dp*hbar_SI/elem_charge_SI ! ~0.0137
!             !
!             do i=1,ncomp
!                file_name= trim(seedname)//'-epsS_'//&
!                     achar(119+alps(i))//achar(119+bets(i))//'.dat'
!                file_name=trim(file_name)
!                file_unit=io_file_unit()
!                write(stdout,'(/,3x,a)') '* '//file_name
!                open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
!                do ifreq=1,nfreq
!                   freq=optics_energy_min+(ifreq-1)*d_freq
!                   write(file_unit,'(5E16.8)') freq,&
!                        fac*(sig(1,i,ifreq,:)+sig(2,i,ifreq,:)+sig(3,i,ifreq,:))
!                enddo
!             enddo
             !
             ! ------------------------------------
             ! Re[sigma] in units of 10^15 sec^{-1}
             ! ------------------------------------
             !
             fac=9.0e-4_dp
             !
             do i=1,ncomp
                file_name= trim(seedname)//'-sigS_'//&
                     achar(119+alps(i))//achar(119+bets(i))//'.dat'
                file_name=trim(file_name)
                file_unit=io_file_unit()
                write(stdout,'(/,3x,a)') '* '//file_name
                open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
                do ifreq=1,nfreq
                   freq=optics_energy_min+(ifreq-1)*d_freq
                   write(file_unit,'(5E16.8)') freq,fac*freq*&
                        (sig(1,i,ifreq,:)+sig(2,i,ifreq,:)+sig(3,i,ifreq,:))
                enddo
             enddo
          endif !Todd

          write(stdout,*) ' '

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
             file_name= trim(seedname)//'-sigS_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat'
             file_name=trim(file_name)
          else
             file_name= trim(seedname)//'-rhobw2_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat'
             file_name=trim(file_name)
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
             fac=elem_charge_SI/(2*hbar_SI**2*speedlight_SI**2*eps0_SI)
             ! Now convert to deg/[mm.(eV)^2]
             fac=fac*(180.0_dp/pi)*elem_charge_SI**2/1000.0_dp
          endif
          sig_abc=fac*sig_abc

         do ifreq=1,nfreq
             freq=optics_energy_min+(ifreq-1)*d_freq
             write(tot_unit,'(5E18.8)') freq,&
                  sig_abc(ifreq,1),&   !orbital, matrix-element-term
                  sig_abc(ifreq,2),&   !orbital, energy-term
                  sig_abc(ifreq,3),&   !spin
                  sum(sig_abc(ifreq,:))
          end do

          close(tot_unit)
          
          ! Now write the cumulative contribution as a function of frequency
          !
          if(T_odd) then
             file_name= trim(seedname)//'-cutsigS_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat'
             file_name=trim(file_name)
          else
             file_name= trim(seedname)//'-cutrhobw2_'//&
                  achar(119+alpha)//achar(119+beta)//achar(119+gamma)//'.dat'
             file_name=trim(file_name)
          end if
          write(stdout,'(/,3x,a)') '* '//file_name
          tot_unit=io_file_unit()
          open(tot_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')  
          sig_abc_cut=fac*sig_abc_cut
          do ifreq=1,nfreq_cut
             write(tot_unit,'(5E18.8)') (ifreq-1)*d_freq_cut,&
                  sig_abc_cut(ifreq,1),sig_abc_cut(ifreq,2),&
                  sig_abc_cut(ifreq,3),sum(sig_abc_cut(ifreq,:))
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
          alphaME=alphaME/cell_volume
          !
          ! Now convert to SI
          !
          alphaME=alphaME*elem_charge_SI**2/hbar_SI
          !
          ! At this point we have alpha_ij=del M_j/del E_i, which has units of
          ! admittance (conductance). Following discussion in Coh et al PRB11, 
          ! we now convert to alpha^{EH} (units of inverse velocity) by 
          ! multiplying by the vaccum magnetic permeability mu_0. Finally,
          ! convert from s/m to ps/m
          !
          alphaME=alphaME*4.0_dp*pi*1.0e-7_dp*1.0e12_dp

          do ifreq=1,nfreq
             freq=optics_energy_min+(ifreq-1)*d_freq
             do j=1,3
                do i=1,3
                   write(alpha_unit(i,j),'(5E18.8)') freq,&
                        alphaME(1,ifreq,i,j),& !orbital, matrix-element-term
                        alphaME(2,ifreq,i,j),& !orbital, energy-term
                        alphaME(3,ifreq,i,j),& !spin
                        sum(alphaME(:,ifreq,i,j))
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

         do ifreq=1,nfreq
             freq=optics_energy_min+(ifreq-1)*d_freq
             do p=1,10
                !*******************************
                ! eventually change to 1,2 below
                !*******************************
                write(gamma_unit(p),'(5E18.8)') freq,&
                     imgamma(1,ifreq,p),& !orbital, matrix-element-term
                     imgamma(2,ifreq,p),& !orbital, energy-term
                     imgamma(3,ifreq,p),& !spin (***should vanish***)
                     sum(imgamma(:,ifreq,p))
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
                do ifreq=1,nfreq_cut
                   write(tot_unit,'(2(E18.8,1x))') (ifreq-1)*d_freq_cut,&
                        alphaspn_cut(ifreq,i,j)
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
       

       if (timing_level>1) call io_stopwatch('berry: kpts',2)

    end if !on_root

  end subroutine berry_main


  subroutine get_imf_k(kpt,imf_k)
  !======================================!
  !                                      !
  ! Calculates -2Im[f_{alpha beta}(k)],  !
  ! Eq.33 CTVR06, Eq.6 LVTS12, stores it !
  ! in axial-vector form                 !
  !                                      !
  !===================================== !

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_diagonalize,utility_re_tr,&
                                   utility_im_tr
    use w90_parameters, only     : num_wann
    use w90_postw90_common, only : fourier_R_to_k
    use w90_wan_ham, only        : get_eig_UU_JJ_HH,get_occ_mat
    use w90_get_oper, only       : HH_R,AA_R,FF_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: imf_k(3,3)

    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f(:,:)
    complex(kind=dp), allocatable :: g(:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: OOmega_i(:,:)
    complex(kind=dp), allocatable :: JJp(:,:,:)
    complex(kind=dp), allocatable :: JJm(:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i

    allocate(UU(num_wann,num_wann))
    allocate(f(num_wann,num_wann))
    allocate(g(num_wann,num_wann))
    allocate(JJp(num_wann,num_wann,3))
    allocate(JJm(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(OOmega_i(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    ! Gather W-gauge matrix objects
    !
    call get_eig_UU_JJ_HH(kpt,eig,UU,JJp,JJm,mdum)
    ! occupation f, and g=1-f
    call get_occ_mat(eig,UU,f,g)
    ! AA_a = i<u|del_a u> [Eq. (28) LVTS12]
    call fourier_R_to_k(kpt,AA_R(:,:,:,1),AA(:,:,1),0)
    call fourier_R_to_k(kpt,AA_R(:,:,:,2),AA(:,:,2),0)
    call fourier_R_to_k(kpt,AA_R(:,:,:,3),AA(:,:,3),0)

    ! Cartesian components of the Berry curvature pseudovector
    !
    do i=1,3
       !
       call fourier_R_to_k(kpt,AA_R(:,:,:,bet(i)),mdum,alp(i))
       Oomega_i=mdum
       call fourier_R_to_k(kpt,AA_R(:,:,:,alp(i)),mdum,bet(i))
       OOmega_i=Oomega_i-mdum
       !
       ! Trace formula, Eq.(51) LVTS12
       !
       ! J0 (Omega_bar) term
       mdum=matmul(f,OOmega_i)
       imf_k(1,i)=utility_re_tr(mdum)
       !
       ! J1 (DA) term
       mdum =matmul(AA(:,:,alp(i)),JJp(:,:,bet(i)))&
            +matmul(JJm(:,:,alp(i)),AA(:,:,bet(i)))
       imf_k(2,i)=-2.0_dp*utility_im_tr(mdum)
       !
       ! J2 (DD) term
       mdum=matmul(JJm(:,:,alp(i)),JJp(:,:,bet(i)))
       imf_k(3,i)=-2.0_dp*utility_im_tr(mdum)
       !
    end do

  end subroutine get_imf_k


  subroutine get_imh_k(kpt,imh_k)
  !======================================!
  !                                      !
  ! Calculates -2Im[h_{alpha beta}(k)],  !
  ! Eq.35 CTVR06, Eq.8 LVTS12, stores it !
  ! in axial-vector form                 !
  !                                      !
  !======================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_diagonalize,utility_re_tr,&
                                   utility_im_tr
    use w90_parameters, only     : num_wann
    use w90_postw90_common, only : fourier_R_to_k
    use w90_wan_ham, only        : get_eig_UU_JJ_HH,get_occ_mat
    use w90_get_oper, only       : HH_R,AA_R,FF_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: imh_k(3,3)

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f(:,:)
    complex(kind=dp), allocatable :: g(:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
!    complex(kind=dp), allocatable :: FF(:,:)
    complex(kind=dp), allocatable :: OOmega_i(:,:)
    complex(kind=dp), allocatable :: JJp(:,:,:)
    complex(kind=dp), allocatable :: JJm(:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i

    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
    allocate(f(num_wann,num_wann))
    allocate(g(num_wann,num_wann))
    allocate(JJp(num_wann,num_wann,3))
    allocate(JJm(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(OOmega_i(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    ! Gather W-gauge matrix objects
    !
    call get_eig_UU_JJ_HH(kpt,eig,UU,JJp,JJm,HH)
    ! occupation f, and g=1-f
    call get_occ_mat(eig,UU,f,g)
    ! AA_a = i<u|del_a u> [Eq. (28) LVTS12]
    call fourier_R_to_k(kpt,AA_R(:,:,:,1),AA(:,:,1),0)
    call fourier_R_to_k(kpt,AA_R(:,:,:,2),AA(:,:,2),0)
    call fourier_R_to_k(kpt,AA_R(:,:,:,3),AA(:,:,3),0)

    do i=1,3
       call fourier_R_to_k(kpt,AA_R(:,:,:,bet(i)),mdum,alp(i))
       Oomega_i=mdum
       call fourier_R_to_k(kpt,AA_R(:,:,:,alp(i)),mdum,bet(i))
       OOmega_i=Oomega_i-mdum
       !
       ! Trace formula, Eq.(56) LVTS12
       !
       ! J0 term
       mdum=matmul(f,matmul(HH,OOmega_i))
       imh_k(1,i)=utility_re_tr(mdum)
       mdum=matmul(f,matmul(HH,matmul(AA(:,:,alp(i)),&
            matmul(f,AA(:,:,bet(i))))))
       imh_k(1,i)=imh_k(1,i)+2.0_dp*utility_im_tr(mdum)
       !
       ! J1 term
       mdum=matmul(HH,matmul(AA(:,:,alp(i)),JJp(:,:,bet(i))))&
           +matmul(HH,matmul(JJm(:,:,alp(i)),AA(:,:,bet(i))))
       imh_k(2,i)=-2.0_dp*utility_im_tr(mdum)
       !
       ! J2 term
       mdum=matmul(HH,matmul(JJm(:,:,alp(i)),JJp(:,:,bet(i))))
       imh_k(3,i)=-2.0_dp*utility_im_tr(mdum)
       !
    enddo

  end subroutine get_imh_k


  subroutine get_img_k(kpt,img_k)
  !======================================!
  !                                      !
  ! Calculates -2Im[g_{alpha beta}(k)],  !
  ! Eq.34 CTVR06, Eq.7 LVTS12, stores it !
  ! in axial-vector form                 !
  !                                      !
  !======================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_diagonalize,utility_re_tr,&
                                   utility_im_tr
    use w90_parameters, only     : num_wann
    use w90_postw90_common, only : fourier_R_to_k
    use w90_wan_ham, only        : get_eig_UU_JJ_HH,get_occ_mat
    use w90_get_oper, only       : HH_R,AA_R,BB_R,CC_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: img_k(3,3)

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f(:,:)
    complex(kind=dp), allocatable :: g(:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: BB(:,:,:)
    complex(kind=dp), allocatable :: CC(:,:,:,:)
    complex(kind=dp), allocatable :: LLambda_i(:,:)
    complex(kind=dp), allocatable :: JJp(:,:,:)
    complex(kind=dp), allocatable :: JJm(:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i,j

    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
    allocate(f(num_wann,num_wann))
    allocate(g(num_wann,num_wann))
    allocate(JJp(num_wann,num_wann,3))
    allocate(JJm(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(BB(num_wann,num_wann,3))
    allocate(CC(num_wann,num_wann,3,3))
    allocate(LLambda_i(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    ! Gather W-gauge matrix objects
    !
    call get_eig_UU_JJ_HH(kpt,eig,UU,JJp,JJm,HH)
    ! occupation f, and g=1-f
    call get_occ_mat(eig,UU,f,g)
    do j=1,3
       !
       ! AA_j = i<u|del_j u> [Eq. (28) LVTS12]
       call fourier_R_to_k(kpt,AA_R(:,:,:,j),AA(:,:,j),0)
       ! BB_j = i<u|H|del_j u> [Eq. (34) LVTS12]
       call fourier_R_to_k(kpt,BB_R(:,:,:,j),BB(:,:,j),0)
       do i=1,3
          ! CC_ij = i<del_i u|H|del_j u> [Eq. (35) LVTS12]
          call fourier_R_to_k(kpt,CC_R(:,:,:,i,j),CC(:,:,i,j),0)
       enddo
    enddo

    ! Trace formula, Eq.(66) LVTS12
    !
    do i=1,3
       !
       ! J0 term
       ! LLambda_ij [Eq. (37) LVTS12] in pseudovector form
       LLambda_i=cmplx_i*(CC(:,:,alp(i),bet(i))&
              -conjg(transpose(CC(:,:,alp(i),bet(i)))))
       mdum=matmul(f,LLambda_i)
       img_k(1,i)=utility_re_tr(mdum)
       mdum=matmul(f,matmul(HH,matmul(AA(:,:,alp(i)),&
            matmul(f,AA(:,:,bet(i))))))
       img_k(1,i)=img_k(1,i)-2.0_dp*utility_im_tr(mdum)
       !
       ! J1 term
       mdum=matmul(JJm(:,:,alp(i)),BB(:,:,bet(i)))&
           -matmul(JJm(:,:,bet(i)),BB(:,:,alp(i)))
       img_k(2,i)=-2.0_dp*utility_im_tr(mdum)
       !
       ! J2 term
       mdum=matmul(JJm(:,:,alp(i)),matmul(HH,JJp(:,:,bet(i))))
       img_k(3,i)=-2.0_dp*utility_im_tr(mdum)
    enddo

  end subroutine get_img_k

  !===========================================================!
  !                   PRIVATE PROCEDURES                      ! 
  !===========================================================!

  subroutine get_sig_k(kpt,sig_k,jdos_k)
  !=======================================================!
  !                                                       !
  ! Interband optical conductivity (absorptive), in the   !
  ! independent-particle approximation (no local fields), !
  ! and joint density of states (JDOS)                    !
  !                                                       !
  !====================================================== !

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_diagonalize,utility_rotate,w0gauss
    use w90_parameters, only     : num_wann,optics_energy_min,fermi_energy,&
                                   optics_adpt_smr,optics_adpt_smr_factor,&
                                   optics_smr_index,optics_smr_fixed_en_width,&
                                   optics_smr_max,&
                                   berry_interp_mesh,spn_decomp
    use w90_postw90_common, only : get_occ,kmesh_spacing,fourier_R_to_k
    use w90_wan_ham, only        : get_D_h,get_eig_deleig
    use w90_get_oper, only       : HH_R,AA_R
    use w90_spin, only           : get_spn_nk

    ! Arguments
    !
    real(kind=dp),                   intent(in)    :: kpt(3)
    real(kind=dp), dimension(:,:,:,:), intent(out) :: sig_k
    real(kind=dp), dimension(:,:),   intent(out)   :: jdos_k

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: AA_bar(:,:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    
    ! Adaptive smearing
    !
    real(kind=dp) :: del_eig(num_wann,3),joint_level_spacing,&
                     smear,Delta_k,arg
 
    integer          :: i,n,m,ifreq,case,ncomp
    real(kind=dp)    :: rdum,rvdum(3),eig(num_wann),occ(num_wann),occ_prod,&
                        freq,spn_nk(num_wann)

    real(kind=dp), allocatable :: Aprod(:,:)

    if(T_odd) then
       ncomp=3
    else
       ncomp=6
    endif
    allocate(Aprod(3,ncomp))

    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(D_h(num_wann,num_wann,3))
    allocate(AA_bar(num_wann,num_wann,3))
    allocate(mdum(num_wann,num_wann))

    if(optics_adpt_smr) then
       call get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
       Delta_k=kmesh_spacing(berry_interp_mesh)
    else
       call fourier_R_to_k(kpt,HH_R,HH,0) 
       call utility_diagonalize(HH,num_wann,eig,UU)
    endif

    call get_occ(eig,occ,fermi_energy)
    call get_D_h(delHH,UU,eig,D_h)
    do i=1,3
       call fourier_R_to_k(kpt,AA_R(:,:,:,i),mdum(:,:),0)
       AA_bar(:,:,i)=utility_rotate(mdum,UU,num_wann)
    enddo

    if(spn_decomp) call get_spn_nk(kpt,spn_nk)

    sig_k=0.0_dp
    jdos_k=0.0_dp
    case=0
    do n=1,num_wann
       do m=1,num_wann
          occ_prod=occ(n)*(1.0_dp-occ(m))
          if(occ_prod < 1.0e-7_dp) cycle
          if(spn_decomp) then
             if(spn_nk(n)>=0 .and. spn_nk(m)>=0) then 
                case=1 ! up --> up transition 
             elseif(spn_nk(n)<0 .and. spn_nk(m)<0) then 
                case=2 ! down --> down
             else
                case=3 ! spin-flip
             end if
          end if
          !
          ! Optical matrix elements
          !
          if(T_odd) then ! Eqs.(39,42) YWVS07
             do i=1,ncomp !ncomp=3
                Aprod(1,i)=-aimag(D_h(n,m,alp(i))*D_h(m,n,bet(i)))
                Aprod(2,i)=aimag(cmplx_i*(AA_bar(n,m,alp(i))*D_h(m,n,bet(i))&
                     +D_h(n,m,alp(i))*AA_bar(m,n,bet(i))))
                Aprod(3,i)=aimag(AA_bar(n,m,alp(i))*AA_bar(m,n,bet(i)))
             enddo
          else
             do i=1,ncomp !ncomp=6
                Aprod(1,i)=-aimag(cmplx_i*D_h(n,m,alps(i))*D_h(m,n,bets(i)))
                Aprod(2,i)=-aimag(D_h(n,m,alps(i))*AA_bar(m,n,bets(i))&
                                 +AA_bar(n,m,alps(i))*D_h(m,n,bets(i)))
                Aprod(3,i)=aimag(cmplx_i*AA_bar(n,m,alps(i))*AA_bar(m,n,bets(i)))
             enddo
          end if
             
          if(optics_adpt_smr) then
             !
             ! Eq.(35) YWVS07, except for the factor 1/sqrt(2) (understand!)
             !
             rvdum(:)=del_eig(m,:)-del_eig(n,:)
             joint_level_spacing=sqrt(dot_product(rvdum(:),rvdum(:)))*Delta_k
             smear=min(joint_level_spacing*optics_adpt_smr_factor/sqrt(2.0_dp),&
                  optics_smr_max)
          else
             smear=optics_smr_fixed_en_width
          endif

          do ifreq=1,nfreq   
             
             freq=optics_energy_min+(ifreq-1)*d_freq
             arg=(eig(m)-eig(n)-freq)/smear
             if(abs(arg)>10.0_dp) then !optimisation
                cycle
             else
                !
                ! Adaptive broadening of the delta-function in Eq.(39) YWVS07
                ! 
                ! Ivo: previously hard-coded for M-P, w0gauss(arg,1)
                !
                rdum=w0gauss(arg,optics_smr_index)/smear
                !
                rdum=rdum*occ_prod ! Fermi occupancy factor
             end if
             !
             ! Joint density of states
             !
             jdos_k(ifreq,1)=jdos_k(ifreq,1)+rdum
             if(spn_decomp) jdos_k(ifreq,1+case)&
                  =jdos_k(ifreq,1+case)+rdum
             !
             ! Optical conductivity sigma(omega)/omega 
             ! (exclude factor of omega in Eq.(39) YWVS07)
             ! NOTE: Minus sign therein is a typo
             !
             do i=1,ncomp
                sig_k(1,i,ifreq,1)=sig_k(1,i,ifreq,1)+Aprod(1,i)*rdum
                sig_k(2,i,ifreq,1)=sig_k(2,i,ifreq,1)+Aprod(2,i)*rdum
                sig_k(3,i,ifreq,1)=sig_k(3,i,ifreq,1)+Aprod(3,i)*rdum
                if(spn_decomp) then
                   sig_k(1,i,ifreq,1+case)&
                        =sig_k(1,i,ifreq,1+case)+Aprod(1,i)*rdum
                   sig_k(2,i,ifreq,1+case)&
                        =sig_k(2,i,ifreq,1+case)+Aprod(2,i)*rdum
                   sig_k(3,i,ifreq,1+case)&
                        =sig_k(3,i,ifreq,1+case)+Aprod(3,i)*rdum
                end if
             enddo
          end do !ifreq
       end do !m
    end do !n

  end subroutine get_sig_k


  !################################################
  ! EVERYTHING BELOW THIS POINT IS WORK IN PROGRESS
  !################################################


  subroutine get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
  !==================================================================!
  !                                                                  !
  ! Optical conductivity of an insulator at first order in           !
  ! the wavevector of light (independent-particle, no local fields). !
  ! The reactive component is evaluated, for frequencies assumed     !
  ! to be in the gap.                                                !
  !                                                                  !
  ! * If T_odd==F, on output sig_abc_k contains Im[sigma_{S,abc}(k)] !
  !   ("S" stands for symmetric part under a <--> b)                 !
  !                                                                  !
  ! * If T_odd=T, on output sig_abc_k contains                       !
  !   Re[sigma_{A,abc}(k)]/(hbar.omega)                              !
  !   ("A" stands for antisymmetric under a <--> b)                  !
  !                                                                  !
  ! NOTE: Here a=alpha, b=beta, c=gamma                              !
  !                                                                  !
  !==================================================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i,bohr,eV_au
    use w90_utility, only       : utility_diagonalize,utility_rotate
    use w90_parameters, only    : num_wann,fermi_energy,&
                                  optics_energy_min,optics_energy_max,&
                                  alpha,beta,gamma,sigma_abc_onlyorb
    use w90_postw90_common, only : fourier_R_to_k
    use w90_wan_ham, only       : get_D_h_a,get_deleig_a
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

    integer                       :: i,j,n,m,ifreq
    real(kind=dp)                 :: eig(num_wann),del_eig(num_wann,3),&
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
    
    sig_abc_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann !conduction, inside frozen window
          if(eig(m)<fermi_energy.or.eig(m)-eig(n)>optics_energy_max) cycle
          omega_mn=eig(m)-eig(n)
          do ifreq=1,nfreq
             freq=optics_energy_min+(ifreq-1)*d_freq
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
             sig_abc_k(ifreq,1)=sig_abc_k(ifreq,1)+sig_matel(1) !orb (m)-term
             sig_abc_k(ifreq,2)=sig_abc_k(ifreq,2)+sig_en !orb (e)-term
             sig_abc_k(ifreq,3)=sig_abc_k(ifreq,3)+sig_matel(2) !spin
          enddo !ifreq
       enddo !m
    enddo !n

    ! Now recalculate the dc (static) value sig_abc_k(omega=0), 
    ! including only transitions *below* a given energy, scanned 
    ! between zero and optics_energy_max
    !
    sig_abc_cut_k=0.0_dp
    do n=1,num_wann
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann
          if(eig(m)<fermi_energy.or.eig(m)-eig(n)>optics_energy_max) cycle
          omega_mn=eig(m)-eig(n)
          do ifreq=nfreq_cut,1,-1
             !---------------------------------------!
             if(omega_mn>(ifreq-1)*d_freq_cut) exit !
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
             sig_abc_cut_k(ifreq,1)=sig_abc_cut_k(ifreq,1)+sig_matel(1)
             sig_abc_cut_k(ifreq,2)=sig_abc_cut_k(ifreq,2)+sig_en
             sig_abc_cut_k(ifreq,3)=sig_abc_cut_k(ifreq,3)+sig_matel(2)
          enddo !ifreq
       enddo !m
    enddo !n
    
  end subroutine get_sig_abc_k


  subroutine get_ME_EQ_k(kpt,alphaME_k,imgamma_k)
  !==================================================================!
  !                                                                  !
  ! Traceless optical magnetoelectric (ME) tensor alpha_ij           !
  ! and totally symmetric electric quadrupolar (EQ) tensor gamma_ijl !
  !                                                                  ! 
  !==================================================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i,bohr,eV_au
    use w90_utility, only       : utility_diagonalize,utility_rotate
    use w90_parameters, only    : num_wann,fermi_energy,optics_energy_min,&
                                  optics_energy_max,sigma_abc_onlyorb
    use w90_postw90_common, only : fourier_R_to_k
    use w90_wan_ham, only       : get_D_h_a,get_deleig_a
    use w90_get_oper, only      : HH_R,AA_R,BB_R,CC_R,FF_R,SS_R

    ! Arguments
    !
    real(kind=dp),                     intent(in)  :: kpt(3)
    real(kind=dp), dimension(:,:,:,:), intent(out) :: alphaME_k
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

    integer                       :: i,j,n,m,ifreq,a,b,c,d,p
    real(kind=dp)                 :: eig(num_wann),del_eig(num_wann,3),&
                                     !eig_cut,
                                     sig_matel(2),sig_en,freq,omega_mn,eps_abc,&
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
       enddo !i
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
    
    ! -----------------
    ! Im[sigma_{S,abc}] (symmetric under a <--> b)
    ! -----------------
    !
    imsig_s_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann !conduction, inside frozen window
          if(eig(m)<fermi_energy.or.eig(m)-eig(n)>optics_energy_max) cycle
          omega_mn=eig(m)-eig(n)
          do c=1,3
          do b=1,3
          do a=1,b ! Will symmetrize at the end
             do ifreq=1,nfreq
                freq=optics_energy_min+(ifreq-1)*d_freq
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
                imsig_s_k(1,ifreq,a,b,c)=imsig_s_k(1,ifreq,a,b,c)&
                     +sig_matel(1) !orb (m)-term
                imsig_s_k(2,ifreq,a,b,c)=imsig_s_k(2,ifreq,a,b,c)&
                     +sig_en !orb (e)-term
                imsig_s_k(3,ifreq,a,b,c)=imsig_s_k(3,ifreq,a,b,c)&
                     +sig_matel(2) !spin
             enddo !ifreq
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
    alphaME_k=0.0_dp
    do a=1,3
       do d=1,3
          do b=1,3
             do c=1,3
                eps_abc=real((b-a)*(c-a)*(c-b)/2,kind=dp) ! Levi-Civita symbol
                do ifreq=1,nfreq
                   do i=1,3 ! orb (m)-term, orb (e)-term, and spin
                      alphaME_k(i,ifreq,d,a)=alphaME_k(i,ifreq,d,a)&
                           +eps_abc*imsig_s_k(i,ifreq,d,b,c)/3.0_dp
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
       do ifreq=1,nfreq
          !***************************
          ! eventually change to i=1,2
          !***************************
          do i=1,3
             imgamma_k(i,ifreq,p)=imsig_s_k(i,ifreq,aa(p),bb(p),cc(p))&
                                  +imsig_s_k(i,ifreq,cc(p),aa(p),bb(p))&
                                  +imsig_s_k(i,ifreq,bb(p),cc(p),aa(p))
          enddo
       enddo
    enddo
    imgamma_k=imgamma_k/3.0_dp

  end subroutine get_ME_EQ_k


  subroutine get_MEspn_k(kpt,alphaspn_k,alphaspn_cut_k)
  !===========================================================!
  !                                                           !
  ! Spin-electronic static magnetoelectric tensor alphaspn_ij !
  !                                                           !
  !===========================================================!

    use w90_constants, only     : dp,cmplx_0,cmplx_i,bohr,eV_au
    use w90_utility, only       : utility_diagonalize,utility_rotate
    use w90_parameters, only    : num_wann,fermi_energy,optics_energy_max,&
                                  cell_volume
    use w90_postw90_common, only: fourier_R_to_k
    use w90_wan_ham, only       : get_D_h_a
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

    integer                       :: i,j,n,m,ifreq
    real(kind=dp)                 :: eig(num_wann),fac,rdum

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

    ! hbar^2/(m_e.V_cell) after converting hbar^2/m_e from a.u. to eV.Angstrom^2
    !
    fac=bohr**2/eV_au/cell_volume ! units of ev/Angstrom

    ! Compute alphaspn_k (dimensionless due to multiplication by fac). Also,
    ! recalculate it including only transitions below a given energy, scanned 
    ! between zero and optics_energy_max
    !
    alphaspn_k=0.0_dp
    alphaspn_cut_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy) cycle
       do m=1,num_wann !conduction, inside frozen window
!         if(eig(m)<fermi_energy.or.eig(m)>eig_max) cycle
          if(eig(m)<fermi_energy.or.eig(m)-eig(n)>optics_energy_max) cycle
          do j=1,3
             do i=1,3
                rdum=fac*aimag(cmplx_i*SS_h(n,m,j)*AA_h(m,n,i))/(eig(m)-eig(n))
                alphaspn_k(i,j)=alphaspn_k(i,j)+rdum
          freq: do ifreq=nfreq_cut,1,-1
                   if((eig(m)-eig(n))<(ifreq-1)*d_freq_cut) then
                      alphaspn_cut_k(ifreq,i,j)=alphaspn_cut_k(ifreq,i,j)+rdum
                   else
                      exit freq
                   endif
                enddo freq
             enddo !i
          enddo !j
       enddo !m
    enddo !n

  end subroutine get_MEspn_k

end module w90_berry
