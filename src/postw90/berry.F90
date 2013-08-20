!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2007-13 Jonathan Yates, Arash Mostofi,       !
!                Giovanni Pizzi, Young-Su Lee,               !
!                Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

! ---------------------------------------------------------------
! REFERENCES
!
!   WYSV06 = PRB 74, 195118 (2006)  (dc anomalous Hall conductivity - AHC)
!   YWVS07 = PRB 75, 195121 (2007)  (frequency-dependent conductivity - Kubo)
!   LVTS12 = PRB 85, 014435 (2012)  (orbital magnetization and AHC)
!   CTVR06 = PRB 74, 024408 (2006)  (  "          "       )
!
! ---------------------------------------------------------------
!
! * Implementation in progress (undocumented, largely untested):
!
!   MS10   = PRB 82, 245118 (2010)  (Spatial dispersion)
!     -- Add logical input keyword 'kubo_spatial_dispersion' (default F)
!        If 'T', calculate sig_abc, and the conversion into
!        optical rotatory power (if time-even), or into
!        magnetoelectric + electric-quadrupole response tensors (if time-odd)
!
! * Also undocumented, works for limited purposes only: 
!                                 reading k-points and weights from file

module w90_berry

  use w90_constants, only : dp

  implicit none

  private

  public :: berry_main,get_imf_k_list,get_imfgh_k_list,get_curv_k

  ! Pseudovector <--> Antisymmetric tensor
  !
  ! x <--> (y,z)
  ! y <--> (z,x)
  ! z <--> (x,y)
  !
  integer, dimension(3), parameter :: alpha_A=(/ 2,3,1 /)
  integer, dimension(3), parameter ::  beta_A=(/ 3,1,2 /)

  ! Independent components of a symmetric tensor
  !
  ! 1 <--> xx
  ! 2 <--> yy
  ! 3 <--> zz
  ! 4 <--> xy
  ! 5 <--> xz
  ! 6 <--> yz
  !
  integer, dimension(6), parameter :: alpha_S=(/ 1,2,3,1,1,2 /)
  integer, dimension(6), parameter ::  beta_S=(/ 1,2,3,2,3,3 /)

  integer       :: nfreq_cut,aa(10),bb(10),cc(10)
  real(kind=dp) :: d_freq,d_freq_cut
  logical       :: T_odd

  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine berry_main
  !============================================================!
  !                                                            !
  ! Computes the following quantities:                         !
  !                                                            !
  !   (i) Anomalous Hall conductivity (Berry curvature)        ! 
  !  (ii) Complex optical conductivity (Kubo-Greenwood) & JDOS ! 
  ! (iii) Orbital magnetization                                !
  !                                                            !
  !   -------------- WORK IN PROGRESS ---------------          !
  !  (iv) Optical conductivity at order q^1                    !
  !                                                            !
  !============================================================!

    use w90_constants, only     : dp,pi,cmplx_0,cmplx_i,elem_charge_SI,hbar_SI,&
                                  elec_mass_SI,eV_au,bohr,bohr_magn_SI,&
                                  eps0_SI,speedlight_SI
    use w90_comms, only         : on_root,num_nodes,my_node_id,comms_reduce
    use w90_io, only            : io_error,stdout,io_file_unit,seedname,&
                                  io_stopwatch
    use w90_postw90_common, only : nrpts,irvec,num_int_kpts_on_node,int_kpts,&
                                   weight
    use w90_parameters, only    : timing_level,iprint,num_wann,berry_kmesh,&
                                  !---------remove eventually---------
                                  alpha,beta,gamma,sigma_abc_onlyorb,&
                                  !-----------------------------------
                                  ahc_adpt_kmesh,ahc_adpt_kmesh_thresh,&
                                  wanint_kpoint_file,cell_volume,transl_inv,&
                                  berry_task,berry_curv_unit,spin_decomp,&
                                  kubo_nfreq,kubo_freq_list,nfermi,&
                                  fermi_energy_list
    use w90_get_oper, only      : get_HH_R,get_AA_R,get_BB_R,get_CC_R,&
                                  get_FF_R,get_SS_R

    real(kind=dp), allocatable    :: adkpt(:,:) ! Adaptive refinement mesh

    ! AHC and orbital magnetization, calculated for a list of Fermi levels
    !
    ! First index labels J0,J1,J2 terms, second labels the Cartesian component 
    !
    real(kind=dp) :: imf_k_list(3,3,nfermi),imf_list(3,3,nfermi)
    real(kind=dp) :: img_k_list(3,3,nfermi),img_list(3,3,nfermi)
    real(kind=dp) :: imh_k_list(3,3,nfermi),imh_list(3,3,nfermi)
    real(kind=dp) :: ahc_list(3,3,nfermi)
    real(kind=dp) :: LCtil_list(3,3,nfermi),ICtil_list(3,3,nfermi),&
                     Morb_list(3,3,nfermi) 
    real(kind=dp) :: imf_k_list_dummy(3,3,nfermi) ! adaptive refinement of AHC
    
    ! Complex optical conductivity, dividided into Hermitean
    ! (absorptive) and anti-Hermitean (reactive) parts
    !
    complex(kind=dp), allocatable :: kubo_H_k(:,:,:)
    complex(kind=dp), allocatable :: kubo_H(:,:,:)
    complex(kind=dp), allocatable :: kubo_AH_k(:,:,:)
    complex(kind=dp), allocatable :: kubo_AH(:,:,:)
    ! decomposition into contribs from up-up, down-down and spin-flip transitions
    complex(kind=dp), allocatable :: kubo_H_k_spn(:,:,:,:)
    complex(kind=dp), allocatable :: kubo_H_spn(:,:,:,:)
    complex(kind=dp), allocatable :: kubo_AH_k_spn(:,:,:,:)
    complex(kind=dp), allocatable :: kubo_AH_spn(:,:,:,:)

    ! Joint density of states
    !
    real(kind=dp), allocatable :: jdos_k(:)
    real(kind=dp), allocatable :: jdos(:)
    ! decomposition into contribs from up-up, down-down and spin-flip transitions
    real(kind=dp), allocatable :: jdos_k_spn(:,:)
    real(kind=dp), allocatable :: jdos_spn(:,:)
    
    ! Reactive spatial-dispersion optical conductivity
    !
    ! sig_abc_(:,1) is the "matrix element term" part of the orbital contrib
    ! sig_abc_(:,2) is the "energy term" part of the orbital contribution
    ! sig_abc_(:,3) is the spin contribution ("matrix element term"-only)
    !
    real(kind=dp), allocatable :: sig_abc_k(:,:)
    real(kind=dp), allocatable :: sig_abc(:,:)

    ! Same as above, but contrib to the static value from optical
    ! transitions below a given frequency, scanned between zero and
    ! the maximum value of the real part of kubo_freq_list
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
    ! between zero and the maximum value of the real part of kubo_freq_list
    !
    real(kind=dp) :: alphaspn_k(3,3)
    real(kind=dp) :: alphaspn(3,3)
    !
    real(kind=dp), allocatable :: alphaspn_cut_k(:,:,:)
    real(kind=dp), allocatable :: alphaspn_cut(:,:,:)

    real(kind=dp)     :: kweight,kweight_adpt,kpt(3),fac,vdum(3),&
                         freq,kpt_ad(3),rdum
    integer           :: n,i,j,k,ikpt,if,ispn,ierr,p,loop_x,loop_y,loop_z,&
                         loop_tot,loop_adpt,adpt_counter_list(nfermi),ifreq,&
                         ndim,jdos_unit,file_unit,tot_unit,alpha_unit(3,3),&
                         gamma_unit(10)
    character(len=24) :: file_name
    logical           :: eval_ahc,eval_morb,eval_sig_ab,eval_sig_abc,&
                         eval_ME_EQ,eval_MEspn,not_scannable,&
                         eval_kubo !***NEW***

    if(nfermi==0) call io_error(&
         'Must set either "fermi_energy," "num_valence_bands," or '&
         //'(fermi_energy_min,fermi_energy_max,nfermi) for optical properties')

    if (timing_level>1.and.on_root) call io_stopwatch('berry: prelims',1)

    ! NOTE: This assumes that real(kubo_freq_list) is a list of
    ! uniformly spaced increasing frequencies
    !
    d_freq=(maxval(real(kubo_freq_list,dp))&
           -minval(real(kubo_freq_list,dp)))/(kubo_nfreq-1)

    rdum=(maxval(real(kubo_freq_list,dp))&
                       -minval(real(kubo_freq_list,dp)))/(kubo_nfreq-1)
    nfreq_cut=max(nint(maxval(real(kubo_freq_list,dp))/rdum),2)
    d_freq_cut=maxval(real(kubo_freq_list,dp))/(nfreq_cut-1)
    
    ! Initialize to .false. all eval_ flags
    !
    eval_ahc=.false.
    eval_morb=.false.
    eval_kubo=.false.
    eval_sig_abc=.false.
    eval_ME_EQ=.false.
    eval_MEspn=.false.
    T_odd=.false.
    if(index(berry_task,'ahc')>0) eval_ahc=.true.
    if(index(berry_task,'morb')>0) eval_morb=.true.
    if(index(berry_task,'kubo')>0) eval_kubo=.true.
    if(index(berry_task,'gyro')>0) then
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

    ! Wannier matrix elements, allocations and initializations
    !
    if(eval_ahc) then
       call get_HH_R 
       call get_AA_R
       imf_list=0.0_dp
       adpt_counter_list=0
    endif

    if(eval_morb) then
       call get_HH_R 
       call get_AA_R
       call get_BB_R
       call get_CC_R
       imf_list=0.0_dp
       img_list=0.0_dp
       imh_list=0.0_dp
    endif

    ! List here berry_tasks that assume nfermi=1
    not_scannable=eval_kubo.or.eval_sig_ab.or.eval_sig_abc&
         .or.eval_ME_EQ.or.eval_MEspn
    if(not_scannable .and. nfermi.ne.1) call io_error(&
         'The berry_task(s) you chose require that you specify a single '&
         //'Fermi energy: scanning the Fermi energy is not implemented')

    if(eval_kubo) then
       call get_HH_R 
       call get_AA_R
       allocate(kubo_H_k(3,3,kubo_nfreq))
       allocate(kubo_H(3,3,kubo_nfreq)) 
       allocate(kubo_AH_k(3,3,kubo_nfreq))
       allocate(kubo_AH(3,3,kubo_nfreq)) 
       allocate(jdos_k(kubo_nfreq))
       allocate(jdos(kubo_nfreq)) 
       kubo_H=cmplx_0
       kubo_AH=cmplx_0
       jdos=0.0_dp
       if(spin_decomp) then
          call get_SS_R
          allocate(kubo_H_k_spn(3,3,3,kubo_nfreq))
          allocate(kubo_H_spn(3,3,3,kubo_nfreq)) 
          allocate(kubo_AH_k_spn(3,3,3,kubo_nfreq))
          allocate(kubo_AH_spn(3,3,3,kubo_nfreq)) 
          allocate(jdos_k_spn(3,kubo_nfreq))
          allocate(jdos_spn(3,kubo_nfreq)) 
          kubo_H_spn=cmplx_0
          kubo_AH_spn=cmplx_0
          jdos_spn=0.0_dp
       endif
    endif

    if(eval_sig_abc) then
       call get_HH_R 
       call get_AA_R
       call get_BB_R
       call get_CC_R
       call get_FF_R
       if(.not.sigma_abc_onlyorb) call get_SS_R
       allocate(sig_abc_k(kubo_nfreq,3))
       allocate(sig_abc(kubo_nfreq,3))
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
       allocate(alphaME_k(3,kubo_nfreq,3,3))
       allocate(alphaME(3,kubo_nfreq,3,3))
       !
       ! N.B.: There should be no spin contribution to gamma, so ultimately want
       !       to use the dimensions (2,kubo_nfreq,10). For the moment include 
       !       spin to check whether its contribution vanishes numerically
       !
       allocate(imgamma_k(3,kubo_nfreq,10))
       allocate(imgamma(3,kubo_nfreq,10))
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

       if(eval_morb) write(stdout,'(/,3x,a)') '* Orbital magnetization'

       if(eval_kubo) then
          if(spin_decomp) then
             write(stdout,'(/,3x,a)')&
                  '* Complex optical conductivity and its spin-decomposition'
             write(stdout,'(/,3x,a)')&
                  '* Joint density of states and its spin-decomposition'
          else
             write(stdout,'(/,3x,a)') '* Complex optical conductivity'
             write(stdout,'(/,3x,a)') '* Joint density of states'
          endif
       endif

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

       if (timing_level>1) then
          call io_stopwatch('berry: prelims',2)
          call io_stopwatch('berry: k-interpolation',1)
       endif

    end if !on_root

    ! Set up adaptive refinement mesh
    !
    allocate(adkpt(3,ahc_adpt_kmesh**3),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating adkpt in berry')
    ikpt=0
    do i=-(ahc_adpt_kmesh-1)/2,(ahc_adpt_kmesh-1)/2
       do j=-(ahc_adpt_kmesh-1)/2,(ahc_adpt_kmesh-1)/2
          do k=-(ahc_adpt_kmesh-1)/2,(ahc_adpt_kmesh-1)/2
             ikpt=ikpt+1 
             adkpt(1,ikpt)=real(i,dp)/(berry_kmesh(1)*ahc_adpt_kmesh)
             adkpt(2,ikpt)=real(j,dp)/(berry_kmesh(2)*ahc_adpt_kmesh)
             adkpt(3,ikpt)=real(k,dp)/(berry_kmesh(3)*ahc_adpt_kmesh)
          end do
       end do
    end do
           
    ! Loop over interpolation k-points
    !
    if(wanint_kpoint_file) then 

       ! NOTE: still need to specify berry_kmesh in the input file
       !  
       !        - Must use the correct nominal value in order to
       !          correctly set up adaptive smearing in kubo
       
       
       if(on_root) write(stdout,'(/,1x,a,i10,a)')&
            'Reading interpolation grid from file kpoint.dat: ',&
            sum(num_int_kpts_on_node),' points'

       ! Loop over k-points on the irreducible wedge of the Brillouin zone,
       ! read from file 'kpoint.dat'
       !
       do loop_tot=1,num_int_kpts_on_node(my_node_id)
          kpt(:)=int_kpts(:,loop_tot)
          kweight=weight(loop_tot)
          kweight_adpt=kweight/ahc_adpt_kmesh**3
          !               .
          ! ***BEGIN COPY OF CODE BLOCK 1***
          !
          if(eval_ahc) then
             call get_imf_k_list(kpt,imf_k_list)
             do if=1,nfermi
                do i=1,3
                   vdum(i)=sum(imf_k_list(1:3,i,if))
                enddo
                if(berry_curv_unit=='bohr2') vdum=vdum/bohr**2
                if(sqrt(dot_product(vdum,vdum))>ahc_adpt_kmesh_thresh) then
                   adpt_counter_list(if)=adpt_counter_list(if)+1
                   do loop_adpt=1,ahc_adpt_kmesh**3
                      ! Using imf_k_list here would corrupt values for other
                      ! frequencies, hence dummy. Only if-th element is used
                      call get_imf_k_list(kpt(:)+adkpt(:,loop_adpt),&
                           imf_k_list_dummy)
                      imf_list(:,:,if)=imf_list(:,:,if)&
                           +imf_k_list_dummy(:,:,if)*kweight_adpt
                   end do
                else
                   imf_list(:,:,if)=imf_list(:,:,if)&
                        +imf_k_list(:,:,if)*kweight
                endif
             enddo
          end if

          if(eval_morb) then
             call get_imfgh_k_list(kpt,imf_k_list,img_k_list,imh_k_list)
             imf_list=imf_list+imf_k_list*kweight
             img_list=img_list+img_k_list*kweight
             imh_list=imh_list+imh_k_List*kweight
          endif

          if(eval_kubo) then
             if(spin_decomp) then
                call get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k,&
                     kubo_H_k_spn,kubo_AH_k_spn,jdos_k_spn)
             else
                call get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k)
             endif
             kubo_H=kubo_H+kubo_H_k*kweight
             kubo_AH=kubo_AH+kubo_AH_k*kweight
             jdos=jdos+jdos_k*kweight
             if(spin_decomp) then
                kubo_H_spn=kubo_H_spn+kubo_H_k_spn*kweight
                kubo_AH_spn=kubo_AH_spn+kubo_AH_k_spn*kweight
                jdos_spn=jdos_spn+jdos_k_spn*kweight
             endif
          endif

          if(eval_sig_abc) then
             call get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
             sig_abc=sig_abc+sig_abc_k*kweight
             sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight
          endif

          if(eval_ME_EQ) then
             call get_ME_EQ_k(kpt,alphaME_k,imgamma_k)
             alphaME=alphaME+alphaME_k*kweight
             imgamma=imgamma+imgamma_k*kweight
          endif

          if(eval_MEspn) then
             call get_MEspn_k(kpt,alphaspn_k,alphaspn_cut_k)
             alphaspn=alphaspn+alphaspn_k*kweight
             alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight
          endif
          !
          ! ***END COPY OF CODE BLOCK 1***

       end do !loop_tot

    else ! Do not read 'kpoint.dat'. Loop over a regular grid in the full BZ

       kweight = 1.0_dp/real(PRODUCT(berry_kmesh),kind=dp)
       kweight_adpt=kweight/ahc_adpt_kmesh**3

       do loop_tot=my_node_id,PRODUCT(berry_kmesh)-1,num_nodes
          loop_x= loop_tot/(berry_kmesh(2)*berry_kmesh(3))
          loop_y=(loop_tot-loop_x*(berry_kmesh(2)&
               *berry_kmesh(3)))/berry_kmesh(3)
          loop_z=loop_tot-loop_x*(berry_kmesh(2)*berry_kmesh(3))&
                -loop_y*berry_kmesh(3)
          kpt(1)=real(loop_x,dp)/real(berry_kmesh(1),dp)
          kpt(2)=real(loop_y,dp)/real(berry_kmesh(2),dp)
          kpt(3)=real(loop_z,dp)/real(berry_kmesh(3),dp)

          ! ***BEGIN CODE BLOCK 1***
          !
          if(eval_ahc) then
             call get_imf_k_list(kpt,imf_k_list)
             do if=1,nfermi
                do i=1,3
                   vdum(i)=sum(imf_k_list(1:3,i,if))
                enddo
                if(berry_curv_unit=='bohr2') vdum=vdum/bohr**2
                if(sqrt(dot_product(vdum,vdum))>ahc_adpt_kmesh_thresh) then
                   adpt_counter_list(if)=adpt_counter_list(if)+1
                   do loop_adpt=1,ahc_adpt_kmesh**3
                      ! Using imf_k_list here would corrupt values for other
                      ! frequencies, hence dummy. Only if-th element is used
                      call get_imf_k_list(kpt(:)+adkpt(:,loop_adpt),&
                           imf_k_list_dummy)
                      imf_list(:,:,if)=imf_list(:,:,if)&
                           +imf_k_list_dummy(:,:,if)*kweight_adpt
                   end do
                else
                   imf_list(:,:,if)=imf_list(:,:,if)&
                        +imf_k_list(:,:,if)*kweight
                endif
             enddo
          end if

          if(eval_morb) then
             call get_imfgh_k_list(kpt,imf_k_list,img_k_list,imh_k_list)
             imf_list=imf_list+imf_k_list*kweight
             img_list=img_list+img_k_list*kweight
             imh_list=imh_list+imh_k_List*kweight
          endif

          if(eval_kubo) then
             if(spin_decomp) then
                call get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k,&
                     kubo_H_k_spn,kubo_AH_k_spn,jdos_k_spn)
             else
                call get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k)
             endif
             kubo_H=kubo_H+kubo_H_k*kweight
             kubo_AH=kubo_AH+kubo_AH_k*kweight
             jdos=jdos+jdos_k*kweight
             if(spin_decomp) then
                kubo_H_spn=kubo_H_spn+kubo_H_k_spn*kweight
                kubo_AH_spn=kubo_AH_spn+kubo_AH_k_spn*kweight
                jdos_spn=jdos_spn+jdos_k_spn*kweight
             endif
          endif

          if(eval_sig_abc) then
             call get_sig_abc_k(kpt,sig_abc_k,sig_abc_cut_k)
             sig_abc=sig_abc+sig_abc_k*kweight
             sig_abc_cut=sig_abc_cut+sig_abc_cut_k*kweight
          endif

          if(eval_ME_EQ) then
             call get_ME_EQ_k(kpt,alphaME_k,imgamma_k)
             alphaME=alphaME+alphaME_k*kweight
             imgamma=imgamma+imgamma_k*kweight
          endif

          if(eval_MEspn) then
             call get_MEspn_k(kpt,alphaspn_k,alphaspn_cut_k)
             alphaspn=alphaspn+alphaspn_k*kweight
             alphaspn_cut=alphaspn_cut+alphaspn_cut_k*kweight
          endif
          !
          ! ***END CODE BLOCK 1***
          
       end do !loop_tot
   
    end if !wanint_kpoint_file

    ! Collect contributions from all nodes    
    !
    if(eval_ahc) then
       call comms_reduce(imf_list(1,1,1),3*3*nfermi,'SUM')
       call comms_reduce(adpt_counter_list(1),nfermi,'SUM')
    endif

    if(eval_morb) then
       call comms_reduce(imf_list(1,1,1),3*3*nfermi,'SUM')
       call comms_reduce(img_list(1,1,1),3*3*nfermi,'SUM')
       call comms_reduce(imh_list(1,1,1),3*3*nfermi,'SUM')
    end if

    if(eval_kubo) then
       call comms_reduce(kubo_H(1,1,1),3*3*kubo_nfreq,'SUM')
       call comms_reduce(kubo_AH(1,1,1),3*3*kubo_nfreq,'SUM')
       call comms_reduce(jdos(1),kubo_nfreq,'SUM')
       if(spin_decomp) then
          call comms_reduce(kubo_H_spn(1,1,1,1),3*3*3*kubo_nfreq,'SUM')
          call comms_reduce(kubo_AH_spn(1,1,1,1),3*3*3*kubo_nfreq,'SUM')
          call comms_reduce(jdos_spn(1,1),3*kubo_nfreq,'SUM')
       endif
    endif

    if(eval_sig_abc) then
       call comms_reduce(sig_abc(1,1),kubo_nfreq*3,'SUM')
       call comms_reduce(sig_abc_cut(1,1),nfreq_cut*3,'SUM')
    endif

    if(eval_ME_EQ) then
       call comms_reduce(alphaME(1,1,1,1),3*kubo_nfreq*3*3,'SUM')
       !
       ! ***************************************************
       ! Eventually change 3*kubo_nfreq*3 --> 2*kubo_nfreq*3
       ! ***************************************************
       !
       call comms_reduce(imgamma(1,1,1),3*kubo_nfreq*10,'SUM')
    endif
    if(eval_MEspn) then
       call comms_reduce(alphaspn(1,1),3*3,'SUM')
       call comms_reduce(alphaspn_cut(1,1,1),nfreq_cut*3*3,'SUM')
    endif
    
    if(on_root) then

       if (timing_level>1) call io_stopwatch('berry: k-interpolation',2)

       write(stdout,'(1x,a)') ' '

       if(eval_ahc .and. ahc_adpt_kmesh.ne.1) then
          if(.not.wanint_kpoint_file) write(stdout, '(1x,a28,3(i0,1x))')&
               'Regular interpolation grid: ',berry_kmesh
          write(stdout, '(1x,a28,3(i0,1x))') 'Adaptive refinement grid: ',&
               ahc_adpt_kmesh,ahc_adpt_kmesh,ahc_adpt_kmesh
          if(berry_curv_unit=='ang2') then
             write(stdout, '(1x,a28,a17,f6.2,a)')&
                  'Refinement threshold: ','Berry curvature >',&
                  ahc_adpt_kmesh_thresh,' Ang^2'
          elseif(berry_curv_unit=='bohr2') then
             write(stdout, '(1x,a28,a17,f6.2,a)')&
                  'Refinement threshold: ','Berry curvature >',&
                  ahc_adpt_kmesh_thresh,' bohr^2'
          endif
          if(nfermi==1) then
             if(wanint_kpoint_file) then
                write(stdout,'(1x,a30,i5,a,f5.2,a)')&
                     ' Points triggering refinement: ',&
                     adpt_counter_list(1),'(',&
                     100*real(adpt_counter_list(1),dp)&
                     /sum(num_int_kpts_on_node),'%)'
             else
                write(stdout,'(1x,a30,i5,a,f5.2,a)')&
                     ' Points triggering refinement: ',&
                     adpt_counter_list(1),'(',&
                     100*real(adpt_counter_list(1),dp)/product(berry_kmesh),'%)'
             endif
          endif
       else
          if(.not.wanint_kpoint_file) write(stdout, '(1x,a20,3(i0,1x))')&
               'Interpolation grid: ',berry_kmesh(1:3)
       endif

       if(eval_ahc) then
          !
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
          ! (i)   Divide by V_c to obtain (1/N) sum_k omega(k)/V_c, with units 
          !       of [L]^{-1} (Berry curvature Omega(k) has units of [L]^2)
          ! (ii)  [L] = Angstrom. Multiply by 10^8 to convert to (cm)^{-1}
          ! (iii) Multiply by -e^2/hbar in SI, with has units ofconductance, 
          !       (Ohm)^{-1}, or Siemens (S), to get the final result in S/cm
          !
          ! ===========================
          ! fac = -e^2/(hbar.V_c*10^-8) 
          ! ===========================
          !
          ! with 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units
          ! --------------------------------------------------------------------
          !
          fac=-1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
          ahc_list(:,:,:)=imf_list(:,:,:)*fac
          if(nfermi>1) then
             write(stdout,'(/,1x,a)')&
                  '---------------------------------'
             write(stdout,'(1x,a)')&
                  'Output data files related to AHC:'
             write(stdout,'(1x,a)')&
                  '---------------------------------'
             file_name=trim(seedname)//'-ahc-fermiscan.dat'
             write(stdout,'(/,3x,a)') '* '//file_name
             file_unit=io_file_unit()
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
          endif
          do if=1,nfermi
             if(nfermi>1) write(file_unit,'(4(F12.6,1x))')&
                  fermi_energy_list(if),sum(ahc_list(:,1,if)),&
                  sum(ahc_list(:,2,if)),sum(ahc_list(:,3,if))
             write(stdout,'(/,1x,a18,F10.4)') 'Fermi energy (ev):',&
                  fermi_energy_list(if)
             if(nfermi>1) then
                if(wanint_kpoint_file) then
                   write(stdout,'(1x,a30,i5,a,f5.2,a)')&
                        ' Points triggering refinement: ',&
                        adpt_counter_list(if),'(',&
                        100*real(adpt_counter_list(if),dp)&
                        /sum(num_int_kpts_on_node),'%)'
                else
                   write(stdout,'(1x,a30,i5,a,f5.2,a)')&
                        ' Points triggering refinement: ',&
                        adpt_counter_list(if),'(',&
                        100*real(adpt_counter_list(if),dp)/product(berry_kmesh),'%)'
                endif
             endif
             write(stdout,'(/,1x,a)')&
                  'AHC (S/cm)       x          y          z'
             if(iprint>1) then
             write(stdout,'(1x,a)')&
                  '=========='
                write(stdout,'(1x,a9,2x,3(f10.4,1x))') 'J0 term :',&
                     ahc_list(1,1,if),ahc_list(1,2,if),ahc_list(1,3,if)
                write(stdout,'(1x,a9,2x,3(f10.4,1x))') 'J1 term :',&
                     ahc_list(2,1,if),ahc_list(2,2,if),ahc_list(2,3,if)
                write(stdout,'(1x,a9,2x,3(f10.4,1x))') 'J2 term :',&
                     ahc_list(3,1,if),ahc_list(3,2,if),ahc_list(3,3,if)
                write(stdout,'(1x,a)')&
                     '-------------------------------------------'
                write(stdout,'(1x,a9,2x,3(f10.4,1x),/)') 'Total   :',&
                     sum(ahc_list(:,1,if)),sum(ahc_list(:,2,if)),&
                     sum(ahc_list(:,3,if))
             else
                write(stdout,'(1x,a10,1x,3(f10.4,1x),/)') '==========',&
                     sum(ahc_list(:,1,if)),sum(ahc_list(:,2,if)),&
                     sum(ahc_list(:,3,if))
             endif
          enddo
          if(nfermi>1) close(file_unit)
       endif
       
       if(eval_morb) then
          !
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
          if(nfermi>1) then
             write(stdout,'(/,1x,a)')&
                  '---------------------------------'
             write(stdout,'(1x,a)')&
                  'Output data files related to the orbital magnetization:'
             write(stdout,'(1x,a)')&
                  '---------------------------------'
             file_name=trim(seedname)//'-morb-fermiscan.dat'
             write(stdout,'(/,3x,a)') '* '//file_name
             file_unit=io_file_unit()
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
          endif
          do if=1,nfermi
             LCtil_list(:,:,if)=(img_list(:,:,if)&
                  -fermi_energy_list(if)*imf_list(:,:,if))*fac
             ICtil_list(:,:,if)=(imh_list(:,:,if)&
                  -fermi_energy_list(if)*imf_list(:,:,if))*fac
             Morb_list(:,:,if)=LCtil_list(:,:,if)+ICtil_list(:,:,if)
             if(nfermi>1) write(file_unit,'(4(F12.6,1x))')&
                  fermi_energy_list(if),sum(Morb_list(1:3,1,if)),&
                  sum(Morb_list(1:3,2,if)),sum(Morb_list(1:3,3,if))
             write(stdout,'(/,/,1x,a,F12.6)') 'Fermi energy (ev) =',&
                  fermi_energy_list(if)
             write(stdout,'(/,/,1x,a)')&
                  'M_orb (bohr magn/cell)        x          y          z'
             if(iprint>1) then
                write(stdout,'(1x,a)')&
                     '======================'
                write(stdout,'(1x,a22,2x,3(f10.4,1x))') 'Local circulation :',&
                     sum(LCtil_list(1:3,1,if)),sum(LCtil_list(1:3,2,if)),&
                     sum(LCtil_list(1:3,3,if))
                write(stdout,'(1x,a22,2x,3(f10.4,1x))') 'Itinerant circulation:',&
                     sum(ICtil_list(1:3,1,if)),sum(ICtil_list(1:3,2,if)),&
                     sum(ICtil_list(1:3,3,if))
                write(stdout,'(1x,a)')&
                     '--------------------------------------------------------'
                write(stdout,'(1x,a22,2x,3(f10.4,1x),/)') 'Total   :',&
                     sum(Morb_list(1:3,1,if)),sum(Morb_list(1:3,2,if)),&
                     sum(Morb_list(1:3,3,if))
             else
                write(stdout,'(1x,a22,2x,3(f10.4,1x),/)')&
                     '======================',&
                     sum(Morb_list(1:3,1,if)),sum(Morb_list(1:3,2,if)),&
                     sum(Morb_list(1:3,3,if))
             endif
          enddo
          if(nfermi>1) close(file_unit)
       endif

       ! -----------------------------!
       ! Complex optical conductivity !
       ! -----------------------------!
       !
       if(eval_kubo) then

          ! Convert to S/cm
          fac=1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
          kubo_H=kubo_H*fac
          kubo_AH=kubo_AH*fac
          if(spin_decomp) then
             kubo_H_spn=kubo_H_spn*fac
             kubo_AH_spn=kubo_AH_spn*fac
          endif

          write(stdout,'(/,1x,a)')&
               '----------------------------------------------------------'
          write(stdout,'(1x,a)')&
               'Output data files related to complex optical conductivity:'
          write(stdout,'(1x,a)')&
               '----------------------------------------------------------'

          ! Symmetric: real (imaginary) part is Hermitean (anti-Hermitean)
          !
          do n=1,6
             i=alpha_S(n)
             j=beta_S(n)
             file_name= trim(seedname)//'-kubo_S_'//&
                  achar(119+i)//achar(119+j)//'.dat'
             file_name=trim(file_name)
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
             do ifreq=1,kubo_nfreq
                if(spin_decomp) then
                   write(file_unit,'(9E16.8)') real(kubo_freq_list(ifreq),dp),&
                        real(0.5_dp*(kubo_H(i,j,ifreq)+kubo_H(j,i,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_AH(i,j,ifreq)+kubo_AH(j,i,ifreq))),&
                        real(0.5_dp*(kubo_H_spn(i,j,1,ifreq)&
                                    +kubo_H_spn(j,i,1,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_AH_spn(i,j,1,ifreq)&
                                      +kubo_AH_spn(j,i,1,ifreq))),&
                        real(0.5_dp*(kubo_H_spn(i,j,2,ifreq)&
                                    +kubo_H_spn(j,i,2,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_AH_spn(i,j,2,ifreq)&
                                      +kubo_AH_spn(j,i,2,ifreq))),&
                        real(0.5_dp*(kubo_H_spn(i,j,3,ifreq)&
                                    +kubo_H_spn(j,i,3,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_AH_spn(i,j,3,ifreq)&
                                      +kubo_AH_spn(j,i,3,ifreq)))
                else
                   write(file_unit,'(3E16.8)') real(kubo_freq_list(ifreq),dp),&
                        real(0.5_dp*(kubo_H(i,j,ifreq)+kubo_H(j,i,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_AH(i,j,ifreq)+kubo_AH(j,i,ifreq)))
                endif
             enddo
             close(file_unit)
          enddo

          ! Antisymmetric: real (imaginary) part is anti-Hermitean (Hermitean)
          !
          do n=1,3
             i=alpha_A(n)
             j=beta_A(n)
             file_name= trim(seedname)//'-kubo_A_'//&
                  achar(119+i)//achar(119+j)//'.dat'
             file_name=trim(file_name)
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
             do ifreq=1,kubo_nfreq
                if(spin_decomp) then
                   write(file_unit,'(9E16.8)') real(kubo_freq_list(ifreq),dp),&
                        real(0.5_dp*(kubo_AH(i,j,ifreq)-kubo_AH(j,i,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_H(i,j,ifreq)-kubo_H(j,i,ifreq))),&
                        real(0.5_dp*(kubo_AH_spn(i,j,1,ifreq)&
                                    -kubo_AH_spn(j,i,1,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_H_spn(i,j,1,ifreq)&
                                     -kubo_H_spn(j,i,1,ifreq))),&
                        real(0.5_dp*(kubo_AH_spn(i,j,2,ifreq)&
                                     -kubo_AH_spn(j,i,2,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_H_spn(i,j,2,ifreq)&
                                     -kubo_H_spn(j,i,2,ifreq))),&
                        real(0.5_dp*(kubo_AH_spn(i,j,3,ifreq)&
                                    -kubo_AH_spn(j,i,3,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_H_spn(i,j,3,ifreq)&
                                     -kubo_H_spn(j,i,3,ifreq)))
                else
                   write(file_unit,'(3E16.8)') real(kubo_freq_list(ifreq),dp),&
                        real(0.5_dp*(kubo_AH(i,j,ifreq)-kubo_AH(j,i,ifreq)),dp),&
                        aimag(0.5_dp*(kubo_H(i,j,ifreq)-kubo_H(j,i,ifreq)))
                endif
             enddo
             close(file_unit)
          enddo

          ! Joint density of states
          !
          file_name=trim(seedname)//'-jdos.dat'
          write(stdout,'(/,3x,a)') '* '//file_name
          file_unit=io_file_unit()
          open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
          do ifreq=1,kubo_nfreq
             if(spin_decomp) then
                write(file_unit,'(5E16.8)') real(kubo_freq_list(ifreq),dp),&
                     jdos(ifreq),jdos_spn(:,ifreq)
             else
                write(file_unit,'(2E16.8)') real(kubo_freq_list(ifreq),dp),&
                     jdos(ifreq)
             endif
          enddo
          close(file_unit)

       endif

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

         do ifreq=1,kubo_nfreq
             write(tot_unit,'(5E18.8)') real(kubo_freq_list(ifreq),dp),&
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

          do ifreq=1,kubo_nfreq
             do j=1,3
                do i=1,3
                   write(alpha_unit(i,j),'(5E18.8)') kubo_freq_list(ifreq),&
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

         do ifreq=1,kubo_nfreq
             do p=1,10
                !*******************************
                ! eventually change to 1,2 below
                !*******************************
                write(gamma_unit(p),'(5E18.8)')&
                     real(kubo_freq_list(ifreq),dp),&
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
       
    end if !on_root

  end subroutine berry_main


  subroutine get_imf_k_list(kpt,imf_k_list)
  !======================================!
  !                                      !
  ! Calculates -2Im[f(k)]                !
  ! [Eq.33 CTVR06, Eq.6 LVTS12] for a    !
  ! list of Fermi energies, and stores   !
  ! it in axial-vector form              !
  !                                      !
  !===================================== !

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_diagonalize,utility_re_tr,&
                                   utility_im_tr
    use w90_parameters, only     : num_wann,nfermi
    use w90_postw90_common, only : fourier_R_to_k_vec
    use w90_wan_ham, only        : get_eig_UU_HH_JJlist,get_occ_mat_list
    use w90_get_oper, only       : AA_R

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), intent(out), dimension(:,:,:) :: imf_k_list

    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f_list(:,:,:)
    complex(kind=dp), allocatable :: g_list(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: OOmega(:,:,:)
    complex(kind=dp), allocatable :: JJp_list(:,:,:,:)
    complex(kind=dp), allocatable :: JJm_list(:,:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i,if

    allocate(UU(num_wann,num_wann))
    allocate(f_list(num_wann,num_wann,nfermi))
    allocate(g_list(num_wann,num_wann,nfermi))
    allocate(JJp_list(num_wann,num_wann,nfermi,3))
    allocate(JJm_list(num_wann,num_wann,nfermi,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(OOmega(num_wann,num_wann,3))
    allocate(mdum(num_wann,num_wann))

    ! Gather W-gauge matrix objects
    !
    call get_eig_UU_HH_JJlist(kpt,eig,UU,mdum,JJp_list,JJm_list)
    ! list of occupations f, and g=1-f
    call get_occ_mat_list(eig,UU,f_list,g_list)

    ! Eqs.(39-40) WYSV06
    !
    call fourier_R_to_k_vec(kpt,AA_R,OO_true=AA,OO_pseudo=OOmega)

    ! Trace formula, Eq.(51) LVTS12
    !
    do if=1,nfermi
       do i=1,3
          !
          ! J0 term (Omega_bar term of WYSV06)
          mdum=matmul(f_list(:,:,if),OOmega(:,:,i))
          imf_k_list(1,i,if)=utility_re_tr(mdum)
          !
          ! J1 term (DA term of WYSV06)
          mdum =matmul(AA(:,:,alpha_A(i)),JJp_list(:,:,if,beta_A(i)))&
               +matmul(JJm_list(:,:,if,alpha_A(i)),AA(:,:,beta_A(i)))
          imf_k_list(2,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
          ! J2 term (DD of WYSV06)
          mdum=matmul(JJm_list(:,:,if,alpha_A(i)),JJp_list(:,:,if,beta_A(i)))
          imf_k_list(3,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
       end do
    enddo

  end subroutine get_imf_k_list


  subroutine get_imfgh_k_list(kpt,imf_k_list,img_k_list,imh_k_list)
  !==============================================================!
  !                                                              !
  ! Calculates together (to reduce the number of Fourier calls), !
  ! the three quantities entering the orbital magnetization:     !
  !                                                              !
  ! * -2Im[f(k)] [Eq.33 CTVR06, Eq.6 LVTS12]                     !
  ! * -2Im[g(k)] [Eq.34 CTVR06, Eq.7 LVTS12]                     !
  ! * -2Im[h(k)] [Eq.35 CTVR06, Eq.8 LVTS12]                     !
  !                                                              !
  ! They are calculated for a list of Fermi energies, and stored !
  ! in axial-vector form.                                        !
  !                                                              !
  !==============================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_diagonalize,utility_re_tr,&
                                   utility_im_tr
    use w90_parameters, only     : num_wann,nfermi
    use w90_postw90_common, only : fourier_R_to_k_vec,fourier_R_to_k
    use w90_wan_ham, only        : get_eig_UU_HH_JJlist,get_occ_mat_list
    use w90_get_oper, only       : AA_R,BB_R,CC_R

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), intent(out), dimension(:,:,:) :: imf_k_list
    real(kind=dp), intent(out), dimension(:,:,:) :: img_k_list
    real(kind=dp), intent(out), dimension(:,:,:) :: imh_k_list

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f_list(:,:,:)
    complex(kind=dp), allocatable :: g_list(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: BB(:,:,:)
    complex(kind=dp), allocatable :: CC(:,:,:,:)
    complex(kind=dp), allocatable :: OOmega(:,:,:)
    complex(kind=dp), allocatable :: LLambda_i(:,:)
    complex(kind=dp), allocatable :: JJp_list(:,:,:,:)
    complex(kind=dp), allocatable :: JJm_list(:,:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i,j,if

    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
    allocate(f_list(num_wann,num_wann,nfermi))
    allocate(g_list(num_wann,num_wann,nfermi))
    allocate(JJp_list(num_wann,num_wann,nfermi,3))
    allocate(JJm_list(num_wann,num_wann,nfermi,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(BB(num_wann,num_wann,3))
    allocate(CC(num_wann,num_wann,3,3))
    allocate(OOmega(num_wann,num_wann,3))
    allocate(LLambda_i(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    ! Gather W-gauge matrix objects
    !
    call get_eig_UU_HH_JJlist(kpt,eig,UU,HH,JJp_list,JJm_list)
    call get_occ_mat_list(eig,UU,f_list,g_list)

    call fourier_R_to_k_vec(kpt,AA_R,OO_true=AA,OO_pseudo=OOmega) 
    call fourier_R_to_k_vec(kpt,BB_R,OO_true=BB)
    do j=1,3
       do i=1,j
          call fourier_R_to_k(kpt,CC_R(:,:,:,i,j),CC(:,:,i,j),0)
          CC(:,:,j,i)=-conjg(transpose(CC(:,:,i,j)))
       enddo
    enddo

    ! Trace formula for -2Im[f], Eq.(51) LVTS12
    !
    do if=1,nfermi
       do i=1,3
          !
          ! J0 term (Omega_bar term of WYSV06)
          mdum=matmul(f_list(:,:,if),OOmega(:,:,i))
          imf_k_list(1,i,if)=utility_re_tr(mdum)
          !
          ! J1 term (DA term of WYSV06)
          mdum =matmul(AA(:,:,alpha_A(i)),JJp_list(:,:,if,beta_A(i)))&
               +matmul(JJm_list(:,:,if,alpha_A(i)),AA(:,:,beta_A(i)))
          imf_k_list(2,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
          ! J2 term (DD of WYSV06)
          mdum=matmul(JJm_list(:,:,if,alpha_A(i)),JJp_list(:,:,if,beta_A(i)))
          imf_k_list(3,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
       end do
    enddo

    ! Trace formula for -2Im[g], Eq.(66) LVTS12
    !
    do i=1,3
       !
       ! J0 term
       ! LLambda_ij [Eq. (37) LVTS12] in pseudovector form
       LLambda_i=cmplx_i*(CC(:,:,alpha_A(i),beta_A(i))&
              -conjg(transpose(CC(:,:,alpha_A(i),beta_A(i)))))
       do if=1,nfermi
          mdum=matmul(f_list(:,:,if),LLambda_i)
          img_k_list(1,i,if)=utility_re_tr(mdum)
          mdum=matmul(f_list(:,:,if),matmul(HH,matmul(AA(:,:,alpha_A(i)),&
               matmul(f_list(:,:,if),AA(:,:,beta_A(i))))))
          img_k_list(1,i,if)=img_k_list(1,i,if)-2.0_dp*utility_im_tr(mdum)
          !
          ! J1 term
          mdum=matmul(JJm_list(:,:,if,alpha_A(i)),BB(:,:,beta_A(i)))&
               -matmul(JJm_list(:,:,if,beta_A(i)),BB(:,:,alpha_A(i)))
          img_k_list(2,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
          ! J2 term
          mdum=matmul(JJm_list(:,:,if,alpha_A(i)),&
               matmul(HH,JJp_list(:,:,if,beta_A(i))))
          img_k_list(3,i,if)=-2.0_dp*utility_im_tr(mdum)
       enddo
    enddo

    ! Trace formula for -2Im[h], Eq.(56) LVTS12
    !
    do if=1,nfermi
       do i=1,3
          !
          ! J0 term
          mdum=matmul(f_list(:,:,if),matmul(HH,OOmega(:,:,i)))
          imh_k_list(1,i,if)=utility_re_tr(mdum)
          mdum=matmul(f_list(:,:,if),matmul(HH,matmul(AA(:,:,alpha_A(i)),&
               matmul(f_list(:,:,if),AA(:,:,beta_A(i))))))
          imh_k_list(1,i,if)=imh_k_list(1,i,if)+2.0_dp*utility_im_tr(mdum)
          !
          ! J1 term
          mdum=matmul(HH,matmul(AA(:,:,alpha_A(i)),JJp_list(:,:,if,beta_A(i))))&
              +matmul(HH,matmul(JJm_list(:,:,if,alpha_A(i)),AA(:,:,beta_A(i))))
          imh_k_list(2,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
          ! J2 term
          mdum=matmul(HH,matmul(JJm_list(:,:,if,alpha_A(i)),JJp_list(:,:,if,beta_A(i))))
          imh_k_list(3,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
       enddo
    enddo

  end subroutine get_imfgh_k_list


  subroutine get_curv_k(kpt,curv_k,UU)
  !====================================================================!
  !                                                                    !
  ! Calculates the non-abelian Berry curvature of the individual Bloch !
  ! eigenstates at k. Optionally, returns the unitary matrix UU.       !
  !                                                                    !
  !====================================================================!


    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_diagonalize,utility_rotate
    use w90_parameters, only     : num_wann,fermi_energy_list
    use w90_postw90_common, only : get_occ,fourier_R_to_k_vec,fourier_R_to_k_new
    use w90_wan_ham, only        : get_D_h
    use w90_get_oper, only       : HH_R,AA_R

    ! Arguments
    !
    real(kind=dp),              intent(in)                  :: kpt(3)
    ! curv_k has dimensions (3,num_wann), for Cartesian component and band index
    real(kind=dp),              intent(out), dimension(:,:) :: curv_k
    complex(kind=dp), optional, intent(out), dimension(:,:) :: UU       

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU_loc(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: OOmega(:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)

    real(kind=dp) :: eig(num_wann),occ(num_wann),&
                     curv_DD(3),curv_DA(3),curv_AA(3)
    integer       :: n,m,i

    allocate(HH(num_wann,num_wann))
    allocate(UU_loc(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(D_h(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(OOmega(num_wann,num_wann,3))
    allocate(mdum(num_wann,num_wann))

    call fourier_R_to_k_new(kpt,HH_R,OO=HH,&
                                OO_dx=delHH(:,:,1),&
                                OO_dy=delHH(:,:,2),&
                                OO_dz=delHH(:,:,3))
    call utility_diagonalize(HH,num_wann,eig,UU_loc)
    call get_occ(eig,occ,fermi_energy_list(1))
    call get_D_h(delHH,UU_loc,eig,D_h)

    call fourier_R_to_k_vec(kpt,AA_R,OO_true=AA,OO_pseudo=OOmega)
    do i=1,3
       AA(:,:,i)=utility_rotate(AA(:,:,i),UU_loc,num_wann) ! AA_bar, Eq.(21) WYSV06
       OOmega(:,:,i)=utility_rotate(OOmega(:,:,i),UU_loc,num_wann) ! Omega_bar
    enddo

    if(present(UU)) UU=UU_loc

    curv_k=0.0_dp
    do n=1,num_wann
       curv_DD(:)=0.0_dp
       curv_DA(:)=0.0_dp
       curv_AA(:)=occ(n)*real(OOmega(n,n,:),dp)
       if(occ(n)<1.0e-7_dp) cycle
       do m=1,num_wann
          do i=1,3
             curv_DD(i)=curv_DD(i)+2.0_dp*occ(n)*(1.0_dp-occ(m))&
                    *aimag(D_h(n,m,alpha_A(i))*D_h(m,n,beta_A(i)))
             curv_DA(i)=curv_DA(i)-2.0_dp*occ(n)*(1.0_dp-occ(m))&
                    *real(D_h(n,m,alpha_A(i))*AA(m,n,beta_A(i))&
                         -D_h(n,m,beta_A(i))*AA(m,n,alpha_A(i)),dp)
             curv_AA(i)=curv_AA(i)&
                       +2.0_dp*occ(n)*occ(m)*aimag(AA(n,m,alpha_A(i))*AA(m,n,beta_A(i)))
          enddo
       enddo
       curv_k(:,n)=curv_DD+curv_DA+curv_AA
    enddo

  end subroutine get_curv_k


  !===========================================================!
  !                   PRIVATE PROCEDURES                      ! 
  !===========================================================!

  subroutine get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k,&
                        kubo_H_k_spn,kubo_AH_k_spn,jdos_k_spn)
  !====================================================================!
  !                                                                    !
  ! Contribution from point k to the complex interband optical         !
  ! conductivity, separated into Hermitian (H) and anti-Hermitian (AH) ! 
  ! parts. Also returns the joint density of states (JDOS).            !
  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_diagonalize,utility_rotate,w0gauss
    use w90_parameters, only     : num_wann,kubo_nfreq,kubo_freq_list,&
                                   fermi_energy_list,kubo_eigval_max,&
                                   kubo_adpt_smr,kubo_smr_fixed_en_width,&
                                   kubo_adpt_smr_max,kubo_adpt_smr_fac,&
                                   kubo_smr_index,berry_kmesh,spin_decomp
    use w90_postw90_common, only : get_occ,fourier_R_to_k_new,&
                                   fourier_R_to_k_vec,kmesh_spacing
    use w90_wan_ham, only        : get_D_h,get_eig_deleig
    use w90_get_oper, only       : HH_R,AA_R
    use w90_spin, only           : get_spin_nk

    ! Arguments
    !
    ! Last three arguments should be present iff spin_decomp=T (but
    ! this is not checked: do it?)
    !
    real(kind=dp),                                  intent(in)  :: kpt(3)
    complex(kind=dp),           dimension(:,:,:),   intent(out) :: kubo_H_k
    complex(kind=dp),           dimension(:,:,:),   intent(out) :: kubo_AH_k
    real(kind=dp),              dimension(:),       intent(out) :: jdos_k
    complex(kind=dp), optional, dimension(:,:,:,:), intent(out) :: kubo_H_k_spn
    complex(kind=dp), optional, dimension(:,:,:,:), intent(out) :: kubo_AH_k_spn
    real(kind=dp),    optional, dimension(:,:),     intent(out) :: jdos_k_spn

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
     
    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann,3),joint_level_spacing,&
                        eta_smr,Delta_k,arg,rvdum(3)
    
    integer          :: i,j,n,m,ifreq,ispn
    real(kind=dp)    :: eig(num_wann),occ(num_wann),delta,&
                        rfac1,rfac2,occ_prod,spn_nk(num_wann)
    complex(kind=dp) :: cfac,omega

    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(D_h(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))

    if(kubo_adpt_smr) then
       call get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
       Delta_k=kmesh_spacing(berry_kmesh)
    else
       call fourier_R_to_k_new(kpt,HH_R,OO=HH,&
                                        OO_dx=delHH(:,:,1),&
                                        OO_dy=delHH(:,:,2),&
                                        OO_dz=delHH(:,:,3))
       call utility_diagonalize(HH,num_wann,eig,UU)
    endif
    call get_occ(eig,occ,fermi_energy_list(1))
    call get_D_h(delHH,UU,eig,D_h)

    call fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
    do i=1,3
       AA(:,:,i)=utility_rotate(AA(:,:,i),UU,num_wann)
    enddo
    AA=AA+cmplx_i*D_h ! Eq.(25) WYSV06

    ! Replace imaginary part of frequency with a fixed value
    if(.not.kubo_adpt_smr .and. kubo_smr_fixed_en_width/=0.0_dp)&
         kubo_freq_list=real(kubo_freq_list,dp)&
                          +cmplx_i*kubo_smr_fixed_en_width

    kubo_H_k=cmplx_0  
    kubo_AH_k=cmplx_0
    jdos_k=0.0_dp
    if(spin_decomp) then
       call get_spin_nk(kpt,spn_nk)
       kubo_H_k_spn=cmplx_0
       kubo_AH_k_spn=cmplx_0
       jdos_k_spn=0.0_dp
    end if
    do m=1,num_wann
       do n=1,num_wann
          if(n==m) cycle
          if(eig(m)>kubo_eigval_max .or. eig(n)>kubo_eigval_max) cycle
          if(spin_decomp) then
             if(spn_nk(n)>=0 .and. spn_nk(m)>=0) then 
                ispn=1 ! up --> up transition 
             elseif(spn_nk(n)<0 .and. spn_nk(m)<0) then 
                ispn=2 ! down --> down
             else
                ispn=3 ! spin-flip
             end if
          end if
          if(kubo_adpt_smr) then
             ! Eq.(35) YWVS07 
             rvdum(:)=del_eig(m,:)-del_eig(n,:)
             joint_level_spacing=sqrt(dot_product(rvdum(:),rvdum(:)))*Delta_k
             eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,&
                  kubo_adpt_smr_max)
          endif
          rfac1=(occ(m)-occ(n))*(eig(m)-eig(n))
          occ_prod=occ(n)*(1.0_dp-occ(m))
          do ifreq=1,kubo_nfreq   
             !
             ! Complex frequency for the anti-Hermitian conductivity
             !
!             if(kubo_adpt_smr) kubo_freq_list(ifreq)=&
!                  real(kubo_freq_list(ifreq),dp)+cmplx_i*eta_smr
!             omega=kubo_freq_list(ifreq)
             if(kubo_adpt_smr) then
                omega=real(kubo_freq_list(ifreq),dp)+cmplx_i*eta_smr
             else
                omega=kubo_freq_list(ifreq)
             endif
             !
             ! Broadened delta function for the Hermitian conductivity and JDOS
             !
             arg=(eig(m)-eig(n)-real(omega,dp))/eta_smr
             ! If only Hermitean part were computed, could speed up
             ! by inserting here 'if(abs(arg)>10.0_dp) cycle'
             delta=w0gauss(arg,kubo_smr_index)/eta_smr
             !
             ! Lorentzian shape (for testing purposes)
!             delta=1.0_dp/(1.0_dp+arg*arg)/pi
!             delta=delta/eta_smr
             !
             jdos_k(ifreq)=jdos_k(ifreq)+occ_prod*delta
             if(spin_decomp)&
                  jdos_k_spn(ispn,ifreq)=jdos_k_spn(ispn,ifreq)+occ_prod*delta
             cfac=cmplx_i*rfac1/(eig(m)-eig(n)-omega)
             rfac2=-pi*rfac1*delta
             do j=1,3
                do i=1,3
                   kubo_H_k(i,j,ifreq)=kubo_H_k(i,j,ifreq)&
                                    +rfac2*AA(n,m,i)*AA(m,n,j)
                   kubo_AH_k(i,j,ifreq)=kubo_AH_k(i,j,ifreq)&
                                    +cfac*AA(n,m,i)*AA(m,n,j)
                   if(spin_decomp) then
                      kubo_H_k_spn(i,j,ispn,ifreq)=&
                           kubo_H_k_spn(i,j,ispn,ifreq)&
                           +rfac2*AA(n,m,i)*AA(m,n,j)
                      kubo_AH_k_spn(i,j,ispn,ifreq)=&
                           kubo_AH_k_spn(i,j,ispn,ifreq)&
                           +cfac*AA(n,m,i)*AA(m,n,j)
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine get_kubo_k


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
    use w90_parameters, only    : num_wann,fermi_energy_list,&
                                  kubo_nfreq,kubo_freq_list,&
                                  kubo_eigval_max,&
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
       call get_D_h_a(mdum,UU,eig,fermi_energy_list(1),D_h(:,:,i))
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

    ! Evaluate, in the H-gauge of the projected subspace, the matrices 
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
       if(eig(n)>fermi_energy_list(1)) cycle
       do m=1,num_wann !conduction, inside frozen window
!          if(eig(m)<fermi_energy_list(1).or.&
!             eig(m)-eig(n)>maxval(real(kubo_freq_list,dp))) cycle
          if(eig(m)<fermi_energy_list(1).or.&
             eig(m)>kubo_eigval_max) cycle
          omega_mn=eig(m)-eig(n)
          do ifreq=1,kubo_nfreq
             if(T_odd) then
                !
                ! ------------------------------
                ! Im[sigma_{S,alpha beta gamma}]
                ! ------------------------------
                !
                omega_fac_matel=2.0_dp*omega_mn&
                     /(omega_mn**2-kubo_freq_list(ifreq)**2)
                omega_fac_en=2.0_dp*omega_mn**3&
                     /(omega_mn**2-kubo_freq_list(ifreq)**2)**2
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
                omega_fac_matel=2.0_dp/(omega_mn**2-kubo_freq_list(ifreq)**2)
                omega_fac_en=(3.0_dp*omega_mn**2-kubo_freq_list(ifreq)**2)&
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
                !
             endif !T_odd
             sig_abc_k(ifreq,1)=sig_abc_k(ifreq,1)+sig_matel(1) !orb (m)-term
             sig_abc_k(ifreq,2)=sig_abc_k(ifreq,2)+sig_en !orb (e)-term
             sig_abc_k(ifreq,3)=sig_abc_k(ifreq,3)+sig_matel(2) !spin
          enddo !ifreq
       enddo !m
    enddo !n

    ! Now recalculate the dc (static) value sig_abc_k(omega=0),
    ! including only transitions *below* a given energy, scanned
    ! between zero and the maximum value of the real part of
    ! kubo_freq_list
    !
    sig_abc_cut_k=0.0_dp
    do n=1,num_wann
       if(eig(n)>fermi_energy_list(1)) cycle
       do m=1,num_wann
          if(eig(m)<fermi_energy_list(1) .or.&
             eig(m)-eig(n)>maxval(real(kubo_freq_list,dp))) cycle
          omega_mn=eig(m)-eig(n)
          do ifreq=nfreq_cut,1,-1
             !---------------------------------------!
             if(omega_mn>(ifreq-1)*d_freq_cut) exit  !
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
    use w90_parameters, only    : num_wann,fermi_energy_list,&
                                  kubo_nfreq,kubo_freq_list,sigma_abc_onlyorb
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

    allocate(imsig_s_k(3,kubo_nfreq,3,3,3)) ! MOVE TO ARGUMENT LIST?
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
       call get_D_h_a(mdum,UU,eig,fermi_energy_list(1),D_h(:,:,i))
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
       if(eig(n)>fermi_energy_list(1)) cycle
       do m=1,num_wann !conduction, inside frozen window
          if(eig(m)<fermi_energy_list(1).or.&
             eig(m)-eig(n)>maxval(real(kubo_freq_list,dp))) cycle
          omega_mn=eig(m)-eig(n)
          do c=1,3
          do b=1,3
          do a=1,b ! Will symmetrize at the end
             do ifreq=1,kubo_nfreq
                omega_fac_matel=2.0_dp*omega_mn/(omega_mn**2-kubo_freq_list(ifreq)**2)
                omega_fac_en=2.0_dp*omega_mn**3/(omega_mn**2-kubo_freq_list(ifreq)**2)**2
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
                do ifreq=1,kubo_nfreq
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
       do ifreq=1,kubo_nfreq
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
    use w90_parameters, only    : num_wann,fermi_energy_list,kubo_freq_list,&
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
       call get_D_h_a(mdum,UU,eig,fermi_energy_list(1),D_h_i(:,:))
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
    ! between zero and the maximum value of the real part of kubo_freq_list
    !
    alphaspn_k=0.0_dp
    alphaspn_cut_k=0.0_dp
    do n=1,num_wann !valence
       if(eig(n)>fermi_energy_list(1)) cycle
       do m=1,num_wann !conduction, inside frozen window
!         if(eig(m)<fermi_energy_list(1).or.eig(m)>eig_max) cycle
          if(eig(m)<fermi_energy_list(1).or.&
             eig(m)-eig(n)>maxval(real(kubo_freq_list,dp))) cycle
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
