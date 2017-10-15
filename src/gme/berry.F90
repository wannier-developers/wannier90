!-*- mode: F90 -*-!
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
!   WYSV06 = PRB 74, 195118 (2006)  (anomalous Hall conductivity - AHC)
!   YWVS07 = PRB 75, 195121 (2007)  (Kubo frequency-dependent conductivity)
!   LVTS12 = PRB 85, 014435 (2012)  (orbital magnetization and AHC)
!   CTVR06 = PRB 74, 024408 (2006)  (  "          "       )
!   ZMS16  = PRL 116, 077201 (2016) (gyrotropic magnetelectric effect- GME)
!
! ---------------------------------------------------------------
!
! * Undocumented, works for limited purposes only: 
!                                 reading k-points and weights from file

module w90_berry

  use w90_constants, only : dp, twopi, elec_mass_SI,pi,eV_au

  implicit none

  private

  public :: berry_main,berry_get_imf_klist,berry_get_imfgh_klist,&
	berry_get_gme_k_list,berry_get_morb_nk,berry_get_morb_k,&
	    berry_get_curv_k,  berry_get_curv_w_k

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

  integer, dimension(9), parameter ::  alpha_xyz=(/ 1,1,1,2,2,2,3,3,3 /)
  integer, dimension(9), parameter ::  beta_xyz=(/  1,2,3,1,2,3,1,2,3 /)


  integer, dimension(9), parameter :: alpha_Axyz=(/ 1,3,2 ,1,3,2 ,1,3,2 /)
  integer, dimension(9), parameter ::  beta_Axyz=(/ 2,1,3 ,2,1,3 ,2,1,3 /)
  integer, dimension(9), parameter :: gamma_Axyz=(/ 1,1,1 ,2,2,2 ,3,3,3 /)


  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine berry_main
  !============================================================!
  !                                                            !
  ! Computes the following quantities:                         !
  !                                                            !
  !   (i) Anomalous Hall conductivity (from Berry curvature)   ! 
  !  (ii) Complex optical conductivity (Kubo-Greenwood) & JDOS ! 
  ! (iii) Orbital magnetization                                !
  !                                                            !
  !============================================================!

    use w90_constants, only      : dp,cmplx_0,elem_charge_SI,hbar_SI,&
                                   eV_au,bohr
    use w90_comms, only          : on_root,num_nodes,my_node_id,comms_reduce,&
					comms_send,comms_recv,root_id
    use w90_io, only             : io_error,stdout,io_file_unit,seedname,&
                                   io_stopwatch
    use w90_postw90_common, only : nrpts,irvec,num_int_kpts_on_node,int_kpts,&
                                   weight
    use w90_wan_ham, only        : wham_get_eig
    use w90_utility, only        : utility_det3
    use w90_parameters, only     : timing_level,iprint,num_wann,berry_kmesh,&
                                   berry_curv_adpt_kmesh,&
                                   berry_curv_adpt_kmesh_thresh,&
                                   gme_berry_adpt_kmesh_thresh,& 
                                   gme_spin_adpt_kmesh_thresh,&  
                                   gme_orb_adpt_kmesh_thresh,&   
                                   wanint_kpoint_file,cell_volume,transl_inv,&
                                   berry_task,gme_task,berry_curv_unit,spin_decomp,&
                                   kubo_nfreq,kubo_freq_list,nfermi,&
                                   fermi_energy_list,num_berry_bands,berry_band_list,&
                                   berry_box,berry_box_corner,spinors


    use w90_get_oper, only       : get_HH_R,get_AA_R,get_BB_R,get_CC_R,&
                                   get_SS_R

    real(kind=dp), allocatable    :: adkpt(:,:)

    ! AHC and orbital magnetization and GME calculated for a list of Fermi levels
    !
    ! First index labels J0,J1,J2 terms, second labels the Cartesian component 
    !
    real(kind=dp) :: imf_k_list(3,3,nfermi),imf_list(3,3,nfermi)
    real(kind=dp) :: img_k_list(3,3,nfermi),img_list(3,3,nfermi)
    real(kind=dp) :: imh_k_list(3,3,nfermi),imh_list(3,3,nfermi)
    real(kind=dp) :: ahc_list(3,3,nfermi)
    real(kind=dp) :: total_spin_k_list(3,nfermi),total_spin_list(3,nfermi)
    
    !
    ! Added Feb 2016 (Ivo)
    !
    real(kind=dp) :: gme_spn_k_list(3,3,nfermi),gme_spn_list(3,3,nfermi),&
                     mu0alpha_gme_spn_list(3,3,nfermi)
    real(kind=dp) :: omega  !for kubo_current
    !
    ! Added April 2016 (Ivo): alternative Fermi-sea calculation of spin GME
    !
    real(kind=dp) :: gme_spn_k_list_alt(3,3,nfermi),&
                     gme_spn_list_alt(3,3,nfermi),&
                     mu0alpha_gme_spn_list_alt(3,3,nfermi)
    !
    real(kind=dp) :: gme_orb_k_list(3,3,nfermi),gme_orb_list(3,3,nfermi),&
                     mu0alpha_gme_orb_list(3,3,nfermi)
    real(kind=dp) :: gme_ber_k_list(3,3,nfermi),gme_berE_k_list(3,3,nfermi),&
		     gme_ber_list(3,3,nfermi),gme_berE_list(3,3,nfermi),&
                     pge_list(3,3,nfermi),&
                     gme_ohm_list(3,3,nfermi),gme_ohm_k_list(3,3,nfermi) 
                     
    logical       :: gme_dorefine  ! Tsirkin
    real(kind=dp) :: xx,yy,zz,xy,xz,yz ! symmetric part of GME tensor
    real(kind=dp) :: x,y,z ! antisymmetric part, as an axial vector

    
    real(kind=dp) :: LCtil_list(3,3,nfermi),ICtil_list(3,3,nfermi),&
                     Morb_list(3,3,nfermi) 
    real(kind=dp) :: imf_k_list_dummy(3,3,nfermi) ! adaptive refinement of AHC
    
    ! Complex optical conductivity, dividided into Hermitean and
    ! anti-Hermitean parts
    !
    complex(kind=dp), allocatable :: kubo_H_k(:,:,:)
    complex(kind=dp), allocatable :: kubo_H(:,:,:)
    real(kind=dp), allocatable :: NOA_sigma_e(:,:,:)
    real(kind=dp), allocatable :: NOA_sigma_e_k(:,:,:)
    real(kind=dp), allocatable :: NOA_sigma_m_orb(:,:,:)
    real(kind=dp), allocatable :: NOA_sigma_m_orb_k(:,:,:)
    real(kind=dp), allocatable :: NOA_sigma_m_spin(:,:,:)
    real(kind=dp), allocatable :: NOA_sigma_m_spin_k(:,:,:)
    real(kind=dp), allocatable :: kubo_current_re_k(:,:,:)
    real(kind=dp), allocatable :: kubo_current_re(:,:,:)
    real(kind=dp), allocatable :: kubo_current_im_k(:,:,:)
    real(kind=dp), allocatable :: kubo_current_im(:,:,:)
    real(kind=dp), allocatable :: beta_CPGE_im_k(:,:,:)
    real(kind=dp), allocatable :: beta_CPGE_im(:,:,:)
    complex(kind=dp), allocatable :: kubo_AH_k(:,:,:)
    complex(kind=dp), allocatable :: kubo_AH(:,:,:)
    ! decomposition into up-up, down-down and spin-flip transitions
    complex(kind=dp), allocatable :: kubo_H_k_spn(:,:,:,:)
    complex(kind=dp), allocatable :: kubo_H_spn(:,:,:,:)
    complex(kind=dp), allocatable :: kubo_AH_k_spn(:,:,:,:)
    complex(kind=dp), allocatable :: kubo_AH_spn(:,:,:,:)

    ! Joint density of states
    !
    real(kind=dp), allocatable :: jdos_k(:)
    real(kind=dp), allocatable :: jdos(:)
    ! decomposition into up-up, down-down and spin-flip transitions
    real(kind=dp), allocatable :: jdos_k_spn(:,:)
    real(kind=dp), allocatable :: jdos_spn(:,:)
    real(kind=dp), allocatable :: eigenval(:,:),eigenval_all(:)
   
    real(kind=dp)     :: kweight,kweight_adpt,kpt(3),dkpt(3),kpt_ad(3),&
                         db1,db2,db3,fac,freq,rdum,vdum(3)

    integer           :: n,n1,i,j,k,ikpt,ixyz,if,ispn,ierr,loop_x,loop_y,loop_z,&
                         loop_xyz,loop_adpt,adpt_counter_list(nfermi),ifreq,&
                         adpt_counter_list_gme_berry(nfermi),&  ! Tsirkin
                         adpt_counter_list_gme_spin(nfermi),&   ! Tsirkin
                         adpt_counter_list_gme_orb(nfermi),&    ! Tsirkin
                         adpt_counter_list_gme_spin_alt(nfermi),&   ! Tsirkin
                         adpt_counter_gme,&    ! Tsirkin
                         file_unit,ikpt_node,nkpt_tot,nkpt_node
    character(len=24) :: file_name
    logical           :: eval_ahc,eval_gme,eval_CPGE,eval_CPGE_2b,eval_berry_sphere,&
    			 eval_morb,eval_kubo,eval_spin,eval_morg,eval_kubo_current,&
			 not_scannable,eval_wpsearch,eval_eigval,eval_noa,eval_noa_spin


    if(nfermi==0) call io_error(&
         'Must specify one or more Fermi levels when berry=true')

    if (timing_level>1.and.on_root) call io_stopwatch('berry: prelims',1)

    ! Mesh spacing in reduced coordinates
    !
    db1=1.0_dp/real(berry_kmesh(1),dp)
    db2=1.0_dp/real(berry_kmesh(2),dp)
    db3=1.0_dp/real(berry_kmesh(3),dp)
    
    eval_ahc=.false.
    eval_noa=.false.
    eval_noa_spin=.false.
    eval_gme=.false.
    eval_CPGE=.false.
    eval_CPGE_2b=.false.
    eval_berry_sphere=.false.
    eval_eigval=.false.
    eval_morb=.false.
    eval_morg=.false.
    eval_spin=.false.
    eval_kubo=.false.
    eval_kubo_current=.false.
    ! Added DGM Jan 2015
    eval_wpsearch=.false.
    if(index(berry_task,'gme')>0) eval_gme=.true.
    if(index(berry_task,'noa')>0) eval_noa=.true.
    if( (index(berry_task,'noa+spin')>0).and.spinors) eval_noa_spin=.true.
    if(index(berry_task,'cpge')>0) eval_CPGE=.true.
    if(index(berry_task,'cpge2b')>0) eval_CPGE_2b=.true.
    if(index(berry_task,'ahc')>0) eval_ahc=.true.
    if(index(berry_task,'morb')>0) eval_morb=.true.
    if(index(berry_task,'morg')>0) eval_morg=.true.
    if(index(berry_task,'spin')>0) eval_spin=.true.
    if(index(berry_task,'sphere')>0) eval_berry_sphere=.true.
    if(index(berry_task,'kubo')>0) eval_kubo=.true.
    if(index(berry_task,'eigval')>0) eval_eigval=.true.
    if(index(berry_task,'kubcur')>0) eval_kubo_current=.true.
    ! Added DGM Jan 2015
    if(index(berry_task,'wpsearch')>0) eval_wpsearch=.true.

    if ( (nfermi==0).and. .not. eval_wpsearch) call io_error(&
         'Must specify one or more Fermi levels when berry=true')

    ! Wannier matrix elements, allocations and initializations
    !
    if(eval_ahc) then
       call get_HH_R 
       call get_AA_R
       imf_list=0.0_dp
       adpt_counter_list=0
    endif

    if(eval_morb.or.eval_morg) then
       call get_HH_R 
       call get_AA_R
       call get_BB_R
       call get_CC_R
       imf_list=0.0_dp
       img_list=0.0_dp
       imh_list=0.0_dp
    endif
    
    if(eval_spin) then
       call get_HH_R 
       call get_SS_R
    endif
    
    
    if(eval_gme) then
       call get_HH_R
       if(index(gme_task,'spin')>0) call get_SS_R
       if(index(gme_task,'berry')>0.or.index(gme_task,'orb')>0) then   ! Tsirkin
          call get_AA_R
       endif
       if(index(gme_task,'orb')>0) then
          call get_BB_R
          call get_CC_R
       endif
       gme_spn_list=0.0_dp
       gme_spn_list_alt=0.0_dp
       gme_orb_list=0.0_dp
       gme_ber_list=0.0_dp
       gme_berE_list=0.0_dp
       gme_ohm_list=0.0_dp
       adpt_counter_gme=0                  !
       adpt_counter_list_gme_berry=0       ! 
       adpt_counter_list_gme_orb=0         ! Tsirkin
       adpt_counter_list_gme_spin=0        !
       adpt_counter_list_gme_spin_alt=0    !
    endif


    ! List here berry_tasks that assume nfermi=1
    !
    not_scannable=eval_kubo.or.eval_noa
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

    if(eval_kubo_current) then
       call get_HH_R 
       call get_AA_R
       allocate(kubo_current_re_k(9,kubo_nfreq,nfermi))
       allocate(kubo_current_re(9,kubo_nfreq,nfermi)) 
       allocate(kubo_current_im_k(9,kubo_nfreq,nfermi))
       allocate(kubo_current_im(9,kubo_nfreq,nfermi)) 
       kubo_current_re=0.0_dp
       kubo_current_im=0.0_dp
    endif

    if(eval_CPGE) then
       call get_HH_R 
       call get_AA_R
       allocate(beta_CPGE_im_k(9,kubo_nfreq,nfermi))
       allocate(beta_CPGE_im(9,kubo_nfreq,nfermi))
       beta_CPGE_im=0.0_dp
    endif


    if(eval_noa) then
       if(on_root) write(stdout,*) "eval_noa - preparation - start"
       call get_HH_R 
       call get_AA_R
       allocate(NOA_sigma_e(3,3,kubo_nfreq))
       allocate(NOA_sigma_e_k(3,3,kubo_nfreq))
       allocate(NOA_sigma_m_orb(3,3,kubo_nfreq))
       allocate(NOA_sigma_m_orb_k(3,3,kubo_nfreq))
       if (eval_noa_spin) then
        if(on_root) write(stdout,*) "eval_noa - get_SS_R - start"
	call get_SS_R
        if(on_root) write(stdout,*) "eval_noa - get_SS_R - end"
	allocate(NOA_sigma_m_spin(3,3,kubo_nfreq))
	allocate(NOA_sigma_m_spin_k(3,3,kubo_nfreq))
        NOA_sigma_m_spin=0.0_dp
        if(on_root) write(stdout,*) "eval_noa -  - start"
       endif
       NOA_sigma_e=0.0_dp
       NOA_sigma_m_orb=0.0_dp
       if(on_root) write(stdout,*) "eval_noa - preparation - OK"
    endif



    if(eval_berry_sphere) then
       call get_HH_R 
       call get_AA_R
    endif
    
    if(eval_wpsearch) then
      
	call berry_eval_wp_search()
	
    endif !wp_search


    if (on_root) then

       write(stdout,'(/,/,1x,a)')&
            'Properties calculated in module  b e r r y'
       write(stdout,'(1x,a)')&
            '------------------------------------------'

       if(eval_ahc) write(stdout,'(/,3x,a)')&
            '* Anomalous Hall conductivity'

       if(eval_morb) write(stdout,'(/,3x,a)') '* Orbital magnetization'
       if(eval_morg) write(stdout,'(/,3x,a)') '* Orbital magnetization (in GME tensor)'
       if(eval_spin) write(stdout,'(/,3x,a)') '* Spin magnetization (in GME tensor)'

       if(eval_gme) write(stdout,'(/,3x,a)') '* Gyrotropic magnetic effect'

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

       if(eval_kubo_current) then
             write(stdout,'(/,3x,a)') '* Current-induced optical conductivity'
       endif


       if(eval_CPGE) then
             write(stdout,'(/,3x,a)') '* Circular Photogalvanic effect'
       endif


       if(eval_NOA) then
             write(stdout,'(/,3x,a)') '* natural optical activity'
       endif

       if(eval_CPGE_2b) then
             write(stdout,'(/,3x,a)') '* Circular Photogalvanic effect calculated from 2-band model with berry-curvature integral'
       endif

       if(eval_berry_sphere) then
             write(stdout,'(/,3x,a)') '* Flux through berry sphere'
    	    call berry_sphere_flux()
       endif

	
       if(transl_inv) then
          if(eval_morb)&
            call io_error('transl_inv=T disabled for morb')
          if(eval_morb)&
            call io_error('transl_inv=T disabled for morb')
          if(eval_morg)&
            call io_error('transl_inv=T disabled for morg')
          if(eval_spin)&
            call io_error('transl_inv=T disabled for spin')
          write(stdout,'(/,1x,a)')&
               'Using a translationally-invariant discretization for the'
          write(stdout,'(1x,a)')&
               'band-diagonal Wannier matrix elements of r, etc.'
       endif

       if (timing_level>1) then
          call io_stopwatch('berry: prelims',2)
          call io_stopwatch('berry: k-interpolation',1)
       endif

    end if !on_root

    ! Set up adaptive refinement mesh
    !
    allocate(adkpt(3,berry_curv_adpt_kmesh**3),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating adkpt in berry')
    ikpt=0
    !
    ! OLD VERSION (only works correctly for odd grids including original point)
    !
    ! do i=-(berry_curv_adpt_kmesh-1)/2,(berry_curv_adpt_kmesh-1)/2
    !    do j=-(berry_curv_adpt_kmesh-1)/2,(berry_curv_adpt_kmesh-1)/2
    !       do k=-(berry_curv_adpt_kmesh-1)/2,(berry_curv_adpt_kmesh-1)/2
    !          ikpt=ikpt+1 
    !          adkpt(1,ikpt)=i*db1/berry_curv_adpt_kmesh
    !          adkpt(2,ikpt)=j*db2/berry_curv_adpt_kmesh
    !          adkpt(3,ikpt)=k*db3/berry_curv_adpt_kmesh
    !       end do
    !    end do
    ! end do
    !
    ! NEW VERSION (both even and odd grids) 
    !
    do i=0,berry_curv_adpt_kmesh-1
       do j=0,berry_curv_adpt_kmesh-1
          do k=0,berry_curv_adpt_kmesh-1
             ikpt=ikpt+1 
             adkpt(1,ikpt)=db1*((i+0.5_dp)/berry_curv_adpt_kmesh-0.5_dp)
             adkpt(2,ikpt)=db2*((j+0.5_dp)/berry_curv_adpt_kmesh-0.5_dp)
             adkpt(3,ikpt)=db3*((k+0.5_dp)/berry_curv_adpt_kmesh-0.5_dp)
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

       ! Loop over k-points on the irreducible wedge of the Brillouin
       ! zone, read from file 'kpoint.dat'
       !
       do loop_xyz=1,num_int_kpts_on_node(my_node_id)
          kpt(:)=int_kpts(:,loop_xyz)
          kweight=weight(loop_xyz)
          kweight_adpt=kweight/berry_curv_adpt_kmesh**3
          !               .
          ! ***BEGIN COPY OF CODE BLOCK 1***
          !
          if(eval_ahc) then
             call berry_get_imf_klist(kpt,imf_k_list)
             do if=1,nfermi
                vdum(1)=sum(imf_k_list(:,1,if))
                vdum(2)=sum(imf_k_list(:,2,if))
                vdum(3)=sum(imf_k_list(:,3,if))
                if(berry_curv_unit=='bohr2') vdum=vdum/bohr**2
                rdum=sqrt(dot_product(vdum,vdum))
                if(rdum>berry_curv_adpt_kmesh_thresh) then
                   adpt_counter_list(if)=adpt_counter_list(if)+1
                   do loop_adpt=1,berry_curv_adpt_kmesh**3
                      ! Using imf_k_list here would corrupt values for other
                      ! frequencies, hence dummy. Only if-th element is used
                      call berry_get_imf_klist(kpt(:)+adkpt(:,loop_adpt),&
                           imf_k_list_dummy)
                      imf_list(:,:,if)=imf_list(:,:,if)&
                           +imf_k_list_dummy(:,:,if)*kweight_adpt
                   end do
                else
                   imf_list(:,:,if)=imf_list(:,:,if)+imf_k_list(:,:,if)*kweight
                endif
             enddo
          end if
	    
	  if (eval_eigval) then 
	    call io_error("eval_eigval is not implemented for wanint_kpoint_file=.true.")
	  endif
	  
          if(eval_morb.or.eval_morg) then
             call berry_get_imfgh_klist(kpt,imf_k_list,img_k_list,imh_k_list)
             imf_list=imf_list+imf_k_list*kweight
             img_list=img_list+img_k_list*kweight
             imh_list=imh_list+imh_k_List*kweight
          endif

          if(eval_spin) then
            call berry_get_total_spin(kpt,total_spin_k_list)
	    total_spin_list=total_spin_list+total_spin_k_list*kweight
          endif


	  if(eval_gme) then
             ! call get_gme_k_list(kpt,gme_spn_k_list,gme_orb_k_list)
             call berry_get_gme_k_list(kpt,gme_spn_k_list,gme_orb_k_list,&
                  gme_ber_k_list,gme_berE_k_list,gme_spn_k_list_alt,gme_ohm_k_list)

        !!  Tsirkin : adaptive refinement
            gme_dorefine=.false.
             do if=1,nfermi
        	
                if (maxval(abs(gme_spn_k_list(:,:,if)))>=gme_spin_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_spin(if)=adpt_counter_list_gme_spin(if)+1
                   gme_dorefine=.true.
                endif
                if (maxval(abs(gme_spn_k_list_alt(:,:,if)))>=gme_spin_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_spin_alt(if)=adpt_counter_list_gme_spin_alt(if)+1
                   gme_dorefine=.true.
                endif
                if (maxval(abs(gme_ber_k_list(:,:,if)))>=gme_berry_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_berry(if)=adpt_counter_list_gme_berry(if)+1
                   gme_dorefine=.true.
                endif
                if (maxval(abs(gme_orb_k_list(:,:,if)))>=gme_orb_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_orb(if)=adpt_counter_list_gme_orb(if)+1
                   gme_dorefine=.true.
                endif
             enddo
             
                if (gme_dorefine) then
            	   adpt_counter_gme=adpt_counter_gme+1
                   do loop_adpt=1,berry_curv_adpt_kmesh**3
                      call berry_get_gme_k_list(kpt(:)+adkpt(:,loop_adpt),gme_spn_k_list,gme_orb_k_list,&
                          gme_ber_k_list,gme_berE_k_list,gme_spn_k_list_alt,gme_ohm_k_list)
                      gme_spn_list=gme_spn_list+gme_spn_k_list*kweight_adpt
                      gme_orb_list=gme_orb_list+gme_orb_k_list*kweight_adpt
                      gme_ber_list=gme_ber_list+gme_ber_k_list*kweight_adpt
                      gme_berE_list=gme_berE_list+gme_berE_k_list*kweight_adpt
                      gme_ohm_list=gme_ohm_list+gme_ohm_k_list*kweight_adpt
                      gme_spn_list_alt=gme_spn_list_alt+gme_spn_k_list_alt*kweight_adpt
                   end do
                else
                   gme_spn_list=gme_spn_list+gme_spn_k_list*kweight
                   gme_orb_list=gme_orb_list+gme_orb_k_list*kweight
                   gme_ber_list=gme_ber_list+gme_ber_k_list*kweight
                   gme_berE_list=gme_berE_list+gme_berE_k_list*kweight
                   gme_ohm_list=gme_ohm_list+gme_ohm_k_list*kweight
                   gme_spn_list_alt=gme_spn_list_alt+gme_spn_k_list_alt*kweight
                endif
!!!!!!!!!!  Tsirkin : end of adaptive refinement

          endif


          if(eval_kubo) then
             if(spin_decomp) then
                call berry_get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k,&
                     kubo_H_k_spn,kubo_AH_k_spn,jdos_k_spn)
             else
                call berry_get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k)
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


          if(eval_kubo_current) then
             call berry_get_kubo_current_k(kpt,kubo_current_re_k,kubo_current_im_k)
             kubo_current_re=kubo_current_re+kubo_current_re_k*kweight
             kubo_current_im=kubo_current_im+kubo_current_im_k*kweight
          endif

          if(eval_CPGE) then
             call berry_get_CPGE_k(kpt,beta_CPGE_im_k,eval_CPGE_2b)
             beta_CPGE_im=beta_CPGE_im+beta_CPGE_im_k*kweight
          endif


          if(eval_noa) then
             if(eval_noa_spin) then
                call berry_get_NOA_sigma_k(kpt,NOA_sigma_e_k,NOA_sigma_m_orb_k,NOA_sigma_m_spin_k)
                NOA_sigma_m_spin=NOA_sigma_m_spin+NOA_sigma_m_spin_k*kweight
             else
                call berry_get_NOA_sigma_k(kpt,NOA_sigma_e_k,NOA_sigma_m_orb_k)
             endif
             NOA_sigma_e=NOA_sigma_e+NOA_sigma_e_k*kweight
             NOA_sigma_m_orb =NOA_sigma_m_orb +NOA_sigma_m_orb_k*kweight
          endif
         



          !
          ! ***END COPY OF CODE BLOCK 1***

       end do !loop_xyz

    else ! Do not read 'kpoint.dat'. Loop over a regular grid in the full BZ

       kweight=db1*db2*db3*utility_det3(berry_box)
       kweight_adpt=kweight/berry_curv_adpt_kmesh**3



    	nkpt_tot=PRODUCT(berry_kmesh)
    	nkpt_node=nkpt_tot/num_nodes
        if (mod(nkpt_tot,num_nodes).ge.1) nkpt_node=nkpt_node+1
    	ikpt_node=0

	if (eval_eigval) then	
    	    allocate(eigenval(nkpt_node,num_wann))
    	    if(on_root) allocate(eigenval_all(nkpt_node*num_nodes))
    	endif

       do loop_xyz=my_node_id,PRODUCT(berry_kmesh)-1,num_nodes
	  ikpt_node=ikpt_node+1
          loop_x= loop_xyz/(berry_kmesh(2)*berry_kmesh(3))
          loop_y=(loop_xyz-loop_x*(berry_kmesh(2)&
               *berry_kmesh(3)))/berry_kmesh(3)
          loop_z=loop_xyz-loop_x*(berry_kmesh(2)*berry_kmesh(3))&
                -loop_y*berry_kmesh(3)

          
          dkpt(1)=loop_x*db1
          dkpt(2)=loop_y*db2
          dkpt(3)=loop_z*db3
          
          kpt(:)=berry_box_corner(:)+matmul(dkpt,berry_box)

          ! ***BEGIN CODE BLOCK 1***
          !
          if(eval_ahc) then
             call berry_get_imf_klist(kpt,imf_k_list)
             do if=1,nfermi
                vdum(1)=sum(imf_k_list(:,1,if))
                vdum(2)=sum(imf_k_list(:,2,if))
                vdum(3)=sum(imf_k_list(:,3,if))
                if(berry_curv_unit=='bohr2') vdum=vdum/bohr**2
                rdum=sqrt(dot_product(vdum,vdum))
                if(rdum>berry_curv_adpt_kmesh_thresh) then
                   adpt_counter_list(if)=adpt_counter_list(if)+1
                   do loop_adpt=1,berry_curv_adpt_kmesh**3
                      ! Using imf_k_list here would corrupt values for other
                      ! frequencies, hence dummy. Only if-th element is used
                      call berry_get_imf_klist(kpt(:)+adkpt(:,loop_adpt),&
                           imf_k_list_dummy)
                      imf_list(:,:,if)=imf_list(:,:,if)&
                           +imf_k_list_dummy(:,:,if)*kweight_adpt
                   end do
                else
                   imf_list(:,:,if)=imf_list(:,:,if)+imf_k_list(:,:,if)*kweight
                endif
             enddo
          end if

          if(eval_morb.or.eval_morg) then
             call berry_get_imfgh_klist(kpt,imf_k_list,img_k_list,imh_k_list)
             imf_list=imf_list+imf_k_list*kweight
             img_list=img_list+img_k_list*kweight
             imh_list=imh_list+imh_k_List*kweight
          endif

	if (eval_eigval) then	
    	    call wham_get_eig(kpt,eigenval(ikpt_node,:))
    	endif

          if(eval_spin) then
            call berry_get_total_spin(kpt,total_spin_k_list)
	    total_spin_list=total_spin_list+total_spin_k_list*kweight
          endif

          if(eval_gme) then
             ! call get_gme_k_list(kpt,gme_spn_k_list,gme_orb_k_list)
             call berry_get_gme_k_list(kpt,gme_spn_k_list,gme_orb_k_list,&
                  gme_ber_k_list,gme_berE_k_list,gme_spn_k_list_alt,gme_ohm_k_list)

        !!  Tsirkin : adaptive refinement

            gme_dorefine=.false.
             do if=1,nfermi
        	
                if (maxval(abs(gme_spn_k_list(:,:,if)))>=gme_spin_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_spin(if)=adpt_counter_list_gme_spin(if)+1
                   gme_dorefine=.true.
                endif
                if (maxval(abs(gme_spn_k_list_alt(:,:,if)))>=gme_spin_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_spin_alt(if)=adpt_counter_list_gme_spin_alt(if)+1
                   gme_dorefine=.true.
                endif
                if (maxval(abs(gme_ber_k_list(:,:,if)))>=gme_berry_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_berry(if)=adpt_counter_list_gme_berry(if)+1
                   gme_dorefine=.true.
                endif
                if (maxval(abs(gme_orb_k_list(:,:,if)))>=gme_orb_adpt_kmesh_thresh) then
                   adpt_counter_list_gme_orb(if)=adpt_counter_list_gme_orb(if)+1
                   gme_dorefine=.true.
                endif
             enddo
             
                if (gme_dorefine) then
            	   adpt_counter_gme=adpt_counter_gme+1
                   do loop_adpt=1,berry_curv_adpt_kmesh**3
                      call berry_get_gme_k_list(kpt(:)+adkpt(:,loop_adpt),gme_spn_k_list,gme_orb_k_list,&
                          gme_ber_k_list,gme_berE_k_list,gme_spn_k_list_alt,gme_ohm_k_list)
                      gme_spn_list=gme_spn_list+gme_spn_k_list*kweight_adpt
                      gme_orb_list=gme_orb_list+gme_orb_k_list*kweight_adpt
                      gme_ber_list=gme_ber_list+gme_ber_k_list*kweight_adpt
                      gme_berE_list=gme_berE_list+gme_berE_k_list*kweight_adpt
                      gme_ohm_list=gme_ohm_list+gme_ohm_k_list*kweight_adpt
                      gme_spn_list_alt=gme_spn_list_alt+gme_spn_k_list_alt*kweight_adpt
                   end do
                else
                   gme_spn_list=gme_spn_list+gme_spn_k_list*kweight
                   gme_orb_list=gme_orb_list+gme_orb_k_list*kweight
                   gme_ber_list=gme_ber_list+gme_ber_k_list*kweight
                   gme_berE_list=gme_berE_list+gme_berE_k_list*kweight
                   gme_ohm_list=gme_ohm_list+gme_ohm_k_list*kweight
                   gme_spn_list_alt=gme_spn_list_alt+gme_spn_k_list_alt*kweight
                endif
!!!!!!!!        !!  Tsirkin : end of adaptive refinement

          endif


          if(eval_kubo) then
             if(spin_decomp) then
                call berry_get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k,&
                     kubo_H_k_spn,kubo_AH_k_spn,jdos_k_spn)
             else
                call berry_get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k)
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


          if(eval_kubo_current) then
             call berry_get_kubo_current_k(kpt,kubo_current_re_k,kubo_current_im_k)
             kubo_current_re=kubo_current_re+kubo_current_re_k*kweight
             kubo_current_im=kubo_current_im+kubo_current_im_k*kweight
          endif

          if(eval_CPGE) then
             call berry_get_CPGE_k(kpt,beta_CPGE_im_k,eval_CPGE_2b)
             beta_CPGE_im=beta_CPGE_im+beta_CPGE_im_k*kweight
          endif

          if(eval_noa) then
             if(eval_noa_spin) then
                call berry_get_NOA_sigma_k(kpt,NOA_sigma_e_k,NOA_sigma_m_orb_k,NOA_sigma_m_spin_k)
                NOA_sigma_m_spin=NOA_sigma_m_spin+NOA_sigma_m_spin_k*kweight
             else
                call berry_get_NOA_sigma_k(kpt,NOA_sigma_e_k,NOA_sigma_m_orb_k)
             endif
             NOA_sigma_e=NOA_sigma_e+NOA_sigma_e_k*kweight
             NOA_sigma_m_orb =NOA_sigma_m_orb +NOA_sigma_m_orb_k*kweight
          endif


          !
          ! ***END CODE BLOCK 1***
          
       end do !loop_xyz
   
    end if !wanint_kpoint_file

    ! Collect contributions from all nodes    
    !
    if(eval_ahc) then
       call comms_reduce(imf_list(1,1,1),3*3*nfermi,'SUM')
       call comms_reduce(adpt_counter_list(1),nfermi,'SUM')
    endif


    if (eval_eigval) then
      do n1=1,num_berry_bands
        n=berry_band_list(n1)
	if (on_root) then 
	    eigenval_all(::num_nodes)=eigenval(:,n)
	endif
	
	if (.not.on_root) then
	    call comms_send(eigenval(1,n),nkpt_node,root_id)
	else
	    do i=1,num_nodes-1
		call comms_recv(eigenval(1,n),nkpt_node,i)
		eigenval_all(1+i::num_nodes)=eigenval(:,n)
	    enddo
	endif
	
	if (on_root) then
	     write(file_name,201) trim(seedname)//'-eigenval',n
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=trim(file_name),STATUS='UNKNOWN',FORM='FORMATTED')   

             write(file_unit,'(a,3I6)') '# berry kmesh : ', berry_kmesh
                write(file_unit,*) eigenval_all(1:nkpt_tot)
             close(file_unit)
	endif
      enddo
    endif

201 format (a,'_',i3.3,'.dat')


    if(eval_morb.or.eval_morg) then
       call comms_reduce(imf_list(1,1,1),3*3*nfermi,'SUM')
       call comms_reduce(img_list(1,1,1),3*3*nfermi,'SUM')
       call comms_reduce(imh_list(1,1,1),3*3*nfermi,'SUM')
    end if

    if(eval_spin) then
       call comms_reduce(total_spin_list(1,1),3*nfermi,'SUM')
    end if


    if(eval_gme) then
       call comms_reduce(gme_spn_list(1,1,1),3*3*nfermi,'SUM')       
       call comms_reduce(gme_orb_list(1,1,1),3*3*nfermi,'SUM')     
       call comms_reduce(gme_ber_list(1,1,1),3*3*nfermi,'SUM')     
       call comms_reduce(gme_berE_list(1,1,1),3*3*nfermi,'SUM')     
       call comms_reduce(gme_ohm_list(1,1,1),3*3*nfermi,'SUM')     
       call comms_reduce(gme_spn_list_alt(1,1,1),3*3*nfermi,'SUM') 
       call comms_reduce(adpt_counter_gme,1,'SUM')                      !
       call comms_reduce(adpt_counter_list_gme_spin(1),nfermi,'SUM')	!
       call comms_reduce(adpt_counter_list_gme_spin_alt(1),nfermi,'SUM')!  Tsirkin
       call comms_reduce(adpt_counter_list_gme_berry(1),nfermi,'SUM')   !
       call comms_reduce(adpt_counter_list_gme_orb(1),nfermi,'SUM')     !
    endif


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

    if(eval_kubo_current) then
       call comms_reduce(kubo_current_re(1,1,1),9*kubo_nfreq*nfermi,'SUM')
       call comms_reduce(kubo_current_im(1,1,1),9*kubo_nfreq*nfermi,'SUM')
    endif

    if(eval_CPGE) then
       call comms_reduce(beta_CPGE_im(1,1,1),9*kubo_nfreq*nfermi,'SUM')
    endif

    if(eval_NOA) then
       if (eval_noa_spin) call comms_reduce(NOA_sigma_m_spin(1,1,1),9*kubo_nfreq,'SUM')
       call comms_reduce(NOA_sigma_m_orb(1,1,1),9*kubo_nfreq,'SUM')
       call comms_reduce(NOA_sigma_e(1,1,1),9*kubo_nfreq,'SUM')
    endif

    
    if(on_root) then

       if (timing_level>1) call io_stopwatch('berry: k-interpolation',2)
       write(stdout,'(1x,a)') ' '
       if(eval_ahc .and. berry_curv_adpt_kmesh.ne.1) then
          if(.not.wanint_kpoint_file) write(stdout, '(1x,a28,3(i0,1x))')&
               'Regular interpolation grid: ',berry_kmesh
          write(stdout, '(1x,a28,3(i0,1x))') 'Adaptive refinement grid: ',&
               berry_curv_adpt_kmesh,berry_curv_adpt_kmesh,berry_curv_adpt_kmesh
          if(berry_curv_unit=='ang2') then
             write(stdout, '(1x,a28,a17,f6.2,a)')&
                  'Refinement threshold: ','Berry curvature >',&
                  berry_curv_adpt_kmesh_thresh,' Ang^2'
          elseif(berry_curv_unit=='bohr2') then
             write(stdout, '(1x,a28,a17,f6.2,a)')&
                  'Refinement threshold: ','Berry curvature >',&
                  berry_curv_adpt_kmesh_thresh,' bohr^2'
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
                        100*real(adpt_counter_list(if),dp)&
                        /product(berry_kmesh),'%)'
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
                write(stdout,'(1x,a22,2x,3(f10.4,1x))')&
                     'Itinerant circulation:',&
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


       if(eval_morg) then
         fac=eV_au/bohr**2
          if(nfermi>1) then
             write(stdout,'(/,1x,a)')&
                  '---------------------------------'
             write(stdout,'(1x,a)')&
                  'Output data files related to the orbital magnetization:'
             write(stdout,'(1x,a)')&
                  '---------------------------------'
             file_name=trim(seedname)//'-morb-gme-fermiscan.dat'
             write(stdout,'(/,3x,a)') '* '//file_name
             file_unit=io_file_unit()
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
          endif
          do if=1,nfermi
             Morb_list(:,:,if)=fac*(imh_list(:,:,if)-img_list(:,:,if))
             if(nfermi>1) write(file_unit,'(4(F12.6,1x))')&
                  fermi_energy_list(if),sum(Morb_list(1:3,1,if)),&
                  sum(Morb_list(1:3,2,if)),sum(Morb_list(1:3,3,if))
             write(stdout,'(/,/,1x,a,F12.6)') 'Fermi energy (ev) =',&
                  fermi_energy_list(if)
             write(stdout,'(/,/,1x,a)')&
                  'M_orb (bohr magn/cell)        x          y          z'
             write(stdout,'(1x,a22,2x,3(f10.4,1x),/)')&
                     '======================',&
                     sum(Morb_list(1:3,1,if)),sum(Morb_list(1:3,2,if)),&
                     sum(Morb_list(1:3,3,if))
          enddo
          if(nfermi>1) close(file_unit)
       endif


       if(eval_spin) then
          if(nfermi>1) then
             write(stdout,'(/,1x,a)')&
                  '---------------------------------'
             write(stdout,'(1x,a)')&
                  'Output data files related to the spin magnetization:'
             write(stdout,'(1x,a)')&
                  '---------------------------------'
             file_name=trim(seedname)//'-mspin-fermiscan.dat'
             write(stdout,'(/,3x,a)') '* '//file_name
             file_unit=io_file_unit()
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
          endif
          do if=1,nfermi
             if(nfermi>1) write(file_unit,'(4(F12.6,1x))')&
                  fermi_energy_list(if),-total_spin_list(:,if)
             write(stdout,'(/,/,1x,a,F12.6)') 'Fermi energy (ev) =',&
                  fermi_energy_list(if)
             write(stdout,'(/,/,1x,a)')&
                  'Spin magnetization (bohr magn/cell)        x          y          z'
             write(stdout,'(1x,a22,2x,3(f10.4,1x),/)')&
                     '======================',-total_spin_list(:,if)
          enddo
          if(nfermi>1) close(file_unit)
       endif

       ! ---------------------------!
       ! Gyrotropic magnetic effect !
       ! ---------------------------!
       !
       if(eval_gme) then
       
	!!!   Tsirkin - refinement
    	if  (berry_curv_adpt_kmesh.ne.1) then
          if(.not.wanint_kpoint_file) write(stdout, '(1x,a28,3(i0,1x))')&
               'Regular interpolation grid: ',berry_kmesh
          write(stdout, '(1x,a28,3(i0,1x))') 'Adaptive refinement grid: ',&
               berry_curv_adpt_kmesh,berry_curv_adpt_kmesh,berry_curv_adpt_kmesh
          write(stdout, '(1x,a28,a22,f15.2,a)')&
            'Refinement threshold for:','Berry curvature part>',&
                      gme_berry_adpt_kmesh_thresh
          write(stdout, '(1x,a28,a22,f15.2,a)')&
            'Refinement threshold for:','spin part >',&
                      gme_spin_adpt_kmesh_thresh
          write(stdout, '(1x,a28,a22,f15.2,a)')&
            'Refinement threshold for:','orbital part >',&
                      gme_orb_adpt_kmesh_thresh

        write(stdout,'(1x,a20,i20,a,f10.2,a)')&
                     ' Points refined: ',&
                     adpt_counter_gme,'(',&
                     100*real(adpt_counter_gme,dp)/product(berry_kmesh),'%)'
	write(stdout,*) 'Points refined due to :'
	write(stdout,*) " EF   |  orbital  |   spin    |  berry   |  tot  |"
	do if=1,nfermi
	    write(stdout,'(f8.4,a,4(f20.2,a3))') fermi_energy_list(if),"|",&
		    100*real(adpt_counter_list_gme_orb(if),dp)/product(berry_kmesh),"% |",&
		    100*real(adpt_counter_list_gme_spin(if),dp)/product(berry_kmesh),"% |",&
		    100*real(adpt_counter_list_gme_berry(if),dp)/product(berry_kmesh),"% |",&
		    100*real(adpt_counter_gme,dp)/product(berry_kmesh),"% "
	enddo
        endif  !! Tsirkin - refinement

          !
          ! -------------------------------
          ! Conversion of spin contribution
          ! -------------------------------
          !
          ! At this point gme_spn_list contains
          !
          ! (1/N) sum_k delta(E_kn-E_f).(d E_{kn}/d k_i).sigma_{kn,j}
          !
          ! (units of length) in Angstroms. 
          !
          !        ====================================
          ! To get mu0.alpha^GME in units of meter^{-1} do the following:
          !        ====================================
          !
          !   * Divide by V_c in Ang^3 to get a quantity with units of [L]^{-2}
          !   * Multiply by 10^20 to convert to SI
          !   * Multiply by -g_s.e^2/(4m_e) \simeq e^2/(2.m_e) in SI units 
          !     to get alpha^GME in SI
          !   * Multiply by vacuum magnetic permeability in SI, 
          !     mu0 = 4.pi.10^{-7}
          !
          ! ==============================
          ! fac = -2.pi.10^13*e^2/(m_e.V_c)
          ! ==============================
          !
          ! with 'V_c' in Ang^3, and 'e' and 'm_e' in SI units
          !
          fac=-twopi*1.0e13_dp*elem_charge_SI**2/(elec_mass_SI*cell_volume)
          mu0alpha_gme_spn_list(:,:,:)=gme_spn_list(:,:,:)*fac
          mu0alpha_gme_spn_list_alt(:,:,:)=gme_spn_list_alt(:,:,:)*fac
          !
          !
          ! ----------------------------------
          ! Conversion of orbital contribution
          ! ----------------------------------
          !
          ! At this point gme_orb_list contains
          !
          ! (1/N)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i)
          !                           .Im[<del_k u_kn| x (H_k-E_kn)|del_k u_kn>]
          !
          ! (units of energy times length^3) in eV.Ang^3.
          !
          !        ====================================
          ! To get mu0.alpha^GME in units of meter^{-1} do the following:
          !        ====================================
          !
          !   * Divide by V_c in Ang^3 to get a quantity with units of eV
          !   * Multiply by 'e' in SI to convert to SI (Joules)
          !   * Multiply by e^2/(2.hbar^2) to get -alpha^GME in SI
          !   * Multiply by vacuum magnetic permeability in SI, 
          !     mu0 = 4.pi.10^{-7}
          !
          ! ====================================
          ! fac = 2.pi.10^{-7}*e^3/(hbar^2.V_c)
          ! ====================================
          !
          ! with 'V_c' in Ang^3, and 'e' and 'hbar' in SI units
          !
          fac=twopi*1.0e-7_dp*elem_charge_SI**3/(hbar_SI**2*cell_volume)
          mu0alpha_gme_orb_list(:,:,:)=gme_orb_list(:,:,:)*fac
          !
          do if=1,nfermi
             write(stdout,'(/,/,1x,a,F12.6)') 'Fermi energy (eV) =',&
                  fermi_energy_list(if)
             if(index(gme_task,'spin')>0) then
                !
                ! Spin part obtained by integrating over the Fermi surface
                !
                ! Symmetric part
                !
                xx=mu0alpha_gme_spn_list(1,1,if)
                yy=mu0alpha_gme_spn_list(2,2,if)
                zz=mu0alpha_gme_spn_list(3,3,if)
                xy=(mu0alpha_gme_spn_list(1,2,if)&
                   +mu0alpha_gme_spn_list(2,1,if))/2.0_dp
                xz=(mu0alpha_gme_spn_list(1,3,if)&
                   +mu0alpha_gme_spn_list(3,1,if))/2.0_dp
                yz=(mu0alpha_gme_spn_list(2,3,if)&
                   +mu0alpha_gme_spn_list(3,2,if))/2.0_dp
                !
                ! Antisymmetric part, in axial-vector form
                !
                x=(mu0alpha_gme_spn_list(2,3,if)&
                  -mu0alpha_gme_spn_list(3,2,if))/2.0_dp
                y=(mu0alpha_gme_spn_list(3,1,if)&
                  -mu0alpha_gme_spn_list(1,3,if))/2.0_dp
                z=(mu0alpha_gme_spn_list(1,2,if)&
                  -mu0alpha_gme_spn_list(2,1,if))/2.0_dp
                !
                write(stdout,'(/,/,1x,a)')&
                     '-mu0.alpha^GME (1/meter)  [spin part, Fermi-surface]'
                write(stdout,'(1x,a)')&
                     '========================'
                write(stdout,'(/,/,1x,a30)') 'Symmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xx',-xx
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yy',-yy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      zz',-zz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xy',-xy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xz',-xz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yz',-yz
                write(stdout,'(/,1x,a30,2x,e12.4)') ' (1/3)Tr',-(xx+yy+zz)/3._dp
                write(stdout,'(/,/,1x,a30)') 'Antisymmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '       x',-x
                write(stdout,'(/,1x,a30,2x,e12.4)') '       y',-y
                write(stdout,'(/,1x,a30,2x,e12.4)') '       z',-z
                !
                ! Spin part obtained by integrating over the Fermi sea
                !
                ! Symmetric part
                !
                xx=mu0alpha_gme_spn_list_alt(1,1,if)
                yy=mu0alpha_gme_spn_list_alt(2,2,if)
                zz=mu0alpha_gme_spn_list_alt(3,3,if)
                xy=(mu0alpha_gme_spn_list_alt(1,2,if)&
                   +mu0alpha_gme_spn_list_alt(2,1,if))/2.0_dp
                xz=(mu0alpha_gme_spn_list_alt(1,3,if)&
                   +mu0alpha_gme_spn_list_alt(3,1,if))/2.0_dp
                yz=(mu0alpha_gme_spn_list_alt(2,3,if)&
                   +mu0alpha_gme_spn_list_alt(3,2,if))/2.0_dp
                !
                ! Antisymmetric part, in polar-vector form
                !
                x=(mu0alpha_gme_spn_list_alt(2,3,if)&
                  -mu0alpha_gme_spn_list_alt(3,2,if))/2.0_dp
                y=(mu0alpha_gme_spn_list_alt(3,1,if)&
                  -mu0alpha_gme_spn_list_alt(1,3,if))/2.0_dp
                z=(mu0alpha_gme_spn_list_alt(1,2,if)&
                  -mu0alpha_gme_spn_list_alt(2,1,if))/2.0_dp
                !
                write(stdout,'(/,/,1x,a)')&
                     '-mu0.alpha^GME (1/meter)  [spin part, Fermi-sea]'
                write(stdout,'(1x,a)')&
                     '========================'
                write(stdout,'(/,/,1x,a30)') 'Symmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xx',-xx
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yy',-yy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      zz',-zz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xy',-xy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xz',-xz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yz',-yz
                write(stdout,'(/,1x,a30,2x,e12.4)') ' (1/3)Tr',-(xx+yy+zz)/3._dp
                write(stdout,'(/,/,1x,a30)') 'Antisymmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '       x',-x
                write(stdout,'(/,1x,a30,2x,e12.4)') '       y',-y
                write(stdout,'(/,1x,a30,2x,e12.4)') '       z',-z

             endif
             if(index(gme_task,'orb')>0) then
                !
                ! Symmetric part
                !
                !
                ! Symmetric part
                !
                xx=mu0alpha_gme_orb_list(1,1,if)
                yy=mu0alpha_gme_orb_list(2,2,if)
                zz=mu0alpha_gme_orb_list(3,3,if)
                xy=(mu0alpha_gme_orb_list(1,2,if)&
                   +mu0alpha_gme_orb_list(2,1,if))/2.0_dp
                xz=(mu0alpha_gme_orb_list(1,3,if)&
                   +mu0alpha_gme_orb_list(3,1,if))/2.0_dp
                yz=(mu0alpha_gme_orb_list(2,3,if)&
                   +mu0alpha_gme_orb_list(3,2,if))/2.0_dp
                !
                ! Antisymmetric part, in polar-vector form
                !
                x=(mu0alpha_gme_orb_list(2,3,if)&
                  -mu0alpha_gme_orb_list(3,2,if))/2.0_dp
                y=(mu0alpha_gme_orb_list(3,1,if)&
                  -mu0alpha_gme_orb_list(1,3,if))/2.0_dp
                z=(mu0alpha_gme_orb_list(1,2,if)&
                  -mu0alpha_gme_orb_list(2,1,if))/2.0_dp
                !
                write(stdout,'(/,/,1x,a)')&
                     '-mu0.alpha^GME (1/meter)  [orbital part, Fermi-surface]'
                write(stdout,'(1x,a)')&
                     '========================'
                write(stdout,'(/,/,1x,a30)') 'Symmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xx',-xx
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yy',-yy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      zz',-zz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xy',-xy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xz',-xz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yz',-yz
                write(stdout,'(/,1x,a30,2x,e12.4)') ' (1/3)Tr',-(xx+yy+zz)/3._dp
                write(stdout,'(/,/,1x,a30)') 'Antisymmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '       x',-x
                write(stdout,'(/,1x,a30,2x,e12.4)') '       y',-y
                write(stdout,'(/,1x,a30,2x,e12.4)') '       z',-z
             endif
          enddo !if
          !
          !
          ! ----------------------------------------
          ! Conversion of Berry-curvature (PGE) term
          ! ----------------------------------------
          !
          ! At this point gme_ber_list contains
          !
          ! (1/N)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i).Omega_{kn}
	  !
          ! (units of length^3) in Ang^3
          !
          ! and point gme_berE_list contains
          !
          ! (1/N)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i).Omega_{kn}*E_kn
          !
          ! (units of energy*length^3) in Ang^3*eV
          !
          ! To get the dimensionless PGE tensor D_ij in Eq.(10) SF15,
          !  (and D^E in eV)
          ! divide by V_c in Ang^3
          !
          if(index(gme_task,'berry')>0) then
             pge_list(:,:,:)=gme_ber_list(:,:,:)/cell_volume
             do if=1,nfermi
                write(stdout,'(/,/,1x,a,F12.6)') 'Fermi energy (eV) =',&
                     fermi_energy_list(if)
                !
                ! Symmetric part
                !
                xx=pge_list(1,1,if)
                yy=pge_list(2,2,if)
                zz=pge_list(3,3,if)
                xy=(pge_list(1,2,if)+pge_list(2,1,if))/2.0_dp
                xz=(pge_list(1,3,if)+pge_list(3,1,if))/2.0_dp
                yz=(pge_list(2,3,if)+pge_list(3,2,if))/2.0_dp
                !
                ! Antisymmetric part, in polar-vector form
                !
                x=(pge_list(2,3,if)-pge_list(3,2,if))/2.0_dp
                y=(pge_list(3,1,if)-pge_list(1,3,if))/2.0_dp
                z=(pge_list(1,2,if)-pge_list(2,1,if))/2.0_dp
                !
                write(stdout,'(/,/,1x,a)')&
                     'D_ij (dimensionless)  [Berry curvature]'
                write(stdout,'(1x,a)')&
                     '===================='
                write(stdout,'(/,/,1x,a30)') 'Symmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xx',xx
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yy',yy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      zz',zz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xy',xy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xz',xz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yz',yz
                write(stdout,'(/,1x,a30,2x,e12.4)') ' (1/3)Tr',(xx+yy+zz)/3._dp
                write(stdout,'(/,/,1x,a30)') 'Antisymmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '       x',x
                write(stdout,'(/,1x,a30,2x,e12.4)') '       y',y
                write(stdout,'(/,1x,a30,2x,e12.4)') '       z',z
             enddo !if

             pge_list(:,:,:)=gme_berE_list(:,:,:)/cell_volume
             do if=1,nfermi
                write(stdout,'(/,/,1x,a,F12.6)') 'Fermi energy (eV) =',&
                     fermi_energy_list(if)
                !
                ! Symmetric part
                !
                xx=pge_list(1,1,if)
                yy=pge_list(2,2,if)
                zz=pge_list(3,3,if)
                xy=(pge_list(1,2,if)+pge_list(2,1,if))/2.0_dp
                xz=(pge_list(1,3,if)+pge_list(3,1,if))/2.0_dp
                yz=(pge_list(2,3,if)+pge_list(3,2,if))/2.0_dp
                !
                ! Antisymmetric part, in polar-vector form
                !
                x=(pge_list(2,3,if)-pge_list(3,2,if))/2.0_dp
                y=(pge_list(3,1,if)-pge_list(1,3,if))/2.0_dp
                z=(pge_list(1,2,if)-pge_list(2,1,if))/2.0_dp
                !
                write(stdout,'(/,/,1x,a)')&
                     'DE_ij (eV)  [EnBerry curvature ]'
                write(stdout,'(1x,a)')&
                     '===================='
                write(stdout,'(/,/,1x,a30)') 'Symmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xx',xx
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yy',yy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      zz',zz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xy',xy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xz',xz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yz',yz
                write(stdout,'(/,1x,a30,2x,e12.4)') ' (1/3)Tr',(xx+yy+zz)/3._dp
                write(stdout,'(/,/,1x,a30)') 'Antisymmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '       x',x
                write(stdout,'(/,1x,a30,2x,e12.4)') '       y',y
                write(stdout,'(/,1x,a30,2x,e12.4)') '       z',z
             enddo !if
          endif


          ! ----------------------------------------
          ! Conversion of ohmic conductance term
          ! ----------------------------------------
          !
          ! At this point gme_ber_list contains
          !
          ! (1/N)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i).(d E_{kn}/d k_j)
          !
          ! (units of energy*length^2) in eV*Ang^2
          !
          ! To get it in Cab = e/h * (1/N*V_cell)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i).(d E_{kn}/d k_j)
          ! in units Ampere/cm
          ! divide by V_c 
          !
          if(index(gme_task,'ohmic')>0) then
             fac=1.0e+8_dp*elem_charge_SI**2/(twopi*hbar_SI*cell_volume)
             gme_ohm_list(:,:,:)=gme_ohm_list(:,:,:)*fac
             do if=1,nfermi
                write(stdout,'(/,/,1x,a,F12.6)') 'Fermi energy (eV) =',&
                     fermi_energy_list(if)
                !
                ! Symmetric part
                !
                xx=gme_ohm_list(1,1,if)
                yy=gme_ohm_list(2,2,if)
                zz=gme_ohm_list(3,3,if)
                xy=(gme_ohm_list(1,2,if)+gme_ohm_list(2,1,if))/2.0_dp
                xz=(gme_ohm_list(1,3,if)+gme_ohm_list(3,1,if))/2.0_dp
                yz=(gme_ohm_list(2,3,if)+gme_ohm_list(3,2,if))/2.0_dp
                !
                ! Antisymmetric part, in polar-vector form
                !
                x=(gme_ohm_list(2,3,if)-gme_ohm_list(3,2,if))/2.0_dp
                y=(gme_ohm_list(3,1,if)-gme_ohm_list(1,3,if))/2.0_dp
                z=(gme_ohm_list(1,2,if)-gme_ohm_list(2,1,if))/2.0_dp
                !
                write(stdout,'(/,/,1x,a)')&
                     'C_ij (dimensionless)  [Ohmic conductance]'
                write(stdout,'(1x,a)')&
                     '===================='
                write(stdout,'(/,/,1x,a30)') 'Symmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xx',xx
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yy',yy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      zz',zz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xy',xy
                write(stdout,'(/,1x,a30,2x,e12.4)') '      xz',xz
                write(stdout,'(/,1x,a30,2x,e12.4)') '      yz',yz
                write(stdout,'(/,1x,a30,2x,e12.4)') ' (1/3)Tr',(xx+yy+zz)/3._dp
                write(stdout,'(/,/,1x,a30)') 'Antisymmetric part'
                write(stdout,'(/,1x,a30,2x,e12.4)') '       x',x
                write(stdout,'(/,1x,a30,2x,e12.4)') '       y',y
                write(stdout,'(/,1x,a30,2x,e12.4)') '       z',z
             enddo !if
          endif




       endif





       ! -----------------------------!
       ! Complex optical conductivity !
       ! -----------------------------!
       !
       if(eval_kubo) then
          !
          ! Convert to S/cm
          fac=1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
          kubo_H=kubo_H*fac
          kubo_AH=kubo_AH*fac
          if(spin_decomp) then
             kubo_H_spn=kubo_H_spn*fac
             kubo_AH_spn=kubo_AH_spn*fac
          endif
          !
          write(stdout,'(/,1x,a)')&
               '----------------------------------------------------------'
          write(stdout,'(1x,a)')&
               'Output data files related to complex optical conductivity:'
          write(stdout,'(1x,a)')&
               '----------------------------------------------------------'
          !
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
          !
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
          !
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

       !-----------------------------!
       ! Current-induced optical conductivity !
       ! -----------------------------!
       !
       if(eval_kubo_current) then


          ! at this point kubo_current contains 
          ! sigma_k(iab,ifreq,ifermi) =                                        !
          !	SUM_{n,m} delta(E_kn-E_f) (d E_{kn}/d k_i)*                    !
          !      [-2*Im(A_nm_a A_mn_b)]*w_mn^2/(w_mn^2-w^2)                    !
          ! [units of length^3]                                                !
          !  where w_mn=E_{km}-E_{kn}                                          !
          !  We want the result in Dimensionless units   !
          ! 
          !   * Divide by V_c in Ang^3 to make it dimensionless
          fac=1./cell_volume
          kubo_current_re=kubo_current_re*fac
          kubo_current_im=kubo_current_im*fac
          !
          write(stdout,'(/,1x,a)')&
            '------------------------------------------------------------------'
          write(stdout,'(1x,a)')&
            'Output data files related to current-induced optical conductivity:'
          write(stdout,'(1x,a)')&
            '------------------------------------------------------------------'

          file_name= trim(seedname)//'-kubo_cur_freq.dat'
          file_name=trim(file_name)
          file_unit=io_file_unit()
          write(stdout,'(/,3x,a)') '* '//file_name
          open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
          write(file_unit,*) '# FREQUENCIES='
          do ifreq=1,kubo_nfreq
             write(file_unit,'(E16.8)') real(kubo_freq_list(ifreq),dp)
          enddo ! ifreq
          close(file_unit)

          file_name= trim(seedname)//'-kubo_cur_efermi.dat'
          file_name=trim(file_name)
          file_unit=io_file_unit()
          write(stdout,'(/,3x,a)') '* '//file_name
          open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
          write(file_unit,*) '# FERMI_LEVELS='
          do if=1,nfermi
             write(file_unit,'(E16.8)') fermi_energy_list(if)
          enddo ! if
          close(file_unit)

          do n=1,9
             i=alpha_Axyz(n)
             j=beta_Axyz(n)
             k=gamma_Axyz(n)
             file_name= trim(seedname)//'-kubo_cur_'//achar(119+k)//'_'&
                  //achar(119+i)//achar(119+j)//'.dat'
             file_name=trim(file_name)
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
             
             do if=1,nfermi
                write(file_unit,'(a,F10.5)') '# EFERMI=',fermi_energy_list(if)
                do ifreq=1,kubo_nfreq
                   write(file_unit,'(2E16.8)') kubo_current_re(n,ifreq,if),&
                        kubo_current_im(n,ifreq,if)
                enddo ! ifreq
                write(file_unit,*)
                write(file_unit,*)
             enddo !if
             close(file_unit)
          enddo ! n
       endif  ! eval_kubo_current



       !-----------------------------!
       ! Natural optical conductivity !
       ! -----------------------------!
       !
       if(eval_noa) then


          ! at this point NOA_sigma_m, NOA_sigma_e  are in ang^3   !
          !  We want the result in Dimensionless units   !
          ! 
          !   * Divide by V_c in Ang^3 to make it dimensionless
          fac=1./cell_volume
          NOA_sigma_m_orb=NOA_sigma_m_orb*fac
          if (eval_noa_spin) NOA_sigma_m_spin=NOA_sigma_m_spin*fac
          NOA_sigma_e=NOA_sigma_e*fac
          !
          write(stdout,'(/,1x,a)')&
            '------------------------------------------------------------------'
          write(stdout,'(1x,a)')&
            'Output data files related to natural optical conductivity:'
          write(stdout,'(1x,a)')&
            '------------------------------------------------------------------'
             file_name= trim(seedname)//'-NOA-m-orb.dat'
             file_name=trim(file_name)
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
             write(file_unit,'(10a12)') "omega", &
                	(((achar(119+alpha_A(i))//achar(119+beta_A(i))//achar(119+k) ), i=1,3), k=1,3)
             do ifreq=1,kubo_nfreq
                  write(file_unit,'(10E12.4)') real(kubo_freq_list(ifreq),dp),&
                	((NOA_sigma_m_orb(i,k,ifreq), i=1,3), k=1,3)
             enddo ! ifreq
             close(file_unit)

	if (eval_noa_spin) then
             file_name= trim(seedname)//'-NOA-m-spin.dat'
             file_name=trim(file_name)
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
             write(file_unit,'(10a12)') "omega", &
                	(((achar(119+alpha_A(i))//achar(119+beta_A(i))//achar(119+k) ), i=1,3), k=1,3)
             do ifreq=1,kubo_nfreq
                  write(file_unit,'(10E12.4)') real(kubo_freq_list(ifreq),dp),&
                	((NOA_sigma_m_spin(i,k,ifreq), i=1,3), k=1,3)
             enddo ! ifreq
             close(file_unit)
        endif


             file_name= trim(seedname)//'-NOA-e.dat'
             file_name=trim(file_name)
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')

             write(file_unit,'(10a12)') "omega", &
                	(((achar(119+alpha_A(i))//achar(119+beta_A(i))//achar(119+k) ), i=1,3), k=1,3)
             do ifreq=1,kubo_nfreq
                  write(file_unit,'(10E12.4)') real(kubo_freq_list(ifreq),dp),&
                	((NOA_sigma_e(i,k,ifreq), i=1,3), k=1,3)
             enddo ! ifreq
             close(file_unit)


       endif  ! eval_NOA




       !-------------------------------!
       ! Circular photogalvanic effect !
       ! ------------------------------!
       !
       if(eval_CPGE) then

      ! at this point kubo_current contains 
      ! Im[beta0_CPGE(ij,ifreq,ifermi)] =                !
      !	1/Nk*SUM_{k,n,m,a,b} eps_{jab} (f_nk-f_mk) (d E_{knm}/d k_i) Im[A_nm_a A_mn_b] * delta(W-Emn)
      ! OR  (CPGE_2b)
      ! Im[beta0_CPGE(ij,ifreq,ifermi)] =                !
      !	1/Nk*SUM_{k}  (f_1k-f_2k) (d E_{k21}/d k_i) Omega^j_1 * delta(W-E_{k21})
      ! [units of Ang^3]                                                !
      !  We want the quantity  
      !  beta_CPGE=pi*e^3/hbar^2 *beta0_CPGE/V_uc
      !  in units of the quantum pi*e^3/(2*pi*hbar)^2  !
      ! 
      !   * Divide by V_c in Ang^3 to get a quantity in dimensionless units
      !   * multiply  by (2*pi)^2
         fac=(4.*pi*pi)/cell_volume
         
          beta_CPGE_im=beta_CPGE_im*fac
          !
          write(stdout,'(/,1x,a)')&
               '----------------------------------------------------------'
          write(stdout,'(1x,a)')&
               'Output data files related to circular photogalvanic effect'
          write(stdout,'(1x,a)')&
               '----------------------------------------------------------'


             file_name= trim(seedname)//'-CPGE.dat'
             if(eval_CPGE_2b)  file_name= trim(seedname)//'-CPGE-2band.dat'
             file_name=trim(file_name)
             file_unit=io_file_unit()
             write(stdout,'(/,3x,a)') '* '//file_name
             open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')   
	     
             if(eval_CPGE_2b) then 
                write(file_unit,'(a,I3)')&
                     '#  Im(beta)/beta0=(2pi)^2/(Nk*Vuc)*SUM_{k,n} (f_1k-f_2k) &
                     &(d E_{k21}/d k_i) Omega^j_1 * delta(W-E21)'             
             else
                write(file_unit,'(a,I3)')&
                     '#  Im(beta)/beta0=(2pi)^2/(Nk*Vuc)*SUM_{k,n,m,a,b} eps_{jab} &
                     &(f_nk-f_mk) (d E_{knm}/d k_i) Im[A_nm_a A_mn_b] * delta(W-Emn)'
	     endif

             write(file_unit,'(a,I3)') '#          beta0=pi*e^3/(2*pi*hbar)^2  '
             write(file_unit,'(a,I3)') '# Number of Fermi levels : ',nfermi
             do  if=1,nfermi
               write(file_unit,'(a,F10.5,a)') '# EFERMI=',fermi_energy_list(if),' eV '
               write(file_unit,'(a10,9a16)') , '#  W (eV) ',(achar(119+alpha_xyz(n))//achar(119+beta_xyz(n)) , n=1,9 )
               do ifreq=1,kubo_nfreq
                   write(file_unit,'(F10.5,9E16.8)') real(kubo_freq_list(ifreq),dp),(beta_CPGE_im(n,ifreq,if), n=1, 9)
               enddo ! ifreq
               write(file_unit,*)
               write(file_unit,*)
             enddo !if
             close(file_unit)
       endif  ! eval_CPGE


        
    end if !on_root

  end subroutine berry_main


  subroutine berry_get_imf_klist(kpt,imf_k_list,occ)
  !============================================================!
  !                                                            !
  ! Calculates the Berry curvature traced over the occupied    !
  ! states, -2Im[f(k)] [Eq.33 CTVR06, Eq.6 LVTS12] for a list  !
  ! of Fermi energies, and stores it in axial-vector form      !
  !                                                            !
  !============================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_re_tr,utility_im_tr
    use w90_parameters, only     : num_wann,nfermi
    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec
    use w90_wan_ham, only        : wham_get_eig_UU_HH_JJlist,wham_get_occ_mat_list
    use w90_get_oper, only       : AA_R

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), intent(out), dimension(:,:,:) :: imf_k_list
    real(kind=dp), intent(in), optional, dimension(:) :: occ
    
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f_list(:,:,:)
    complex(kind=dp), allocatable :: g_list(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: OOmega(:,:,:)
    complex(kind=dp), allocatable :: JJp_list(:,:,:,:)
    complex(kind=dp), allocatable :: JJm_list(:,:,:,:)
    complex(kind=dp), allocatable :: mdum(:,:)
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i,if,nfermi_loc

    allocate(UU(num_wann,num_wann))
    allocate(AA(num_wann,num_wann,3))
    allocate(OOmega(num_wann,num_wann,3))
    allocate(mdum(num_wann,num_wann))

     if(present(occ)) then
        nfermi_loc=1
     else
        nfermi_loc=nfermi
     endif
     allocate(f_list(num_wann,num_wann,nfermi_loc))
     allocate(g_list(num_wann,num_wann,nfermi_loc))
     allocate(JJp_list(num_wann,num_wann,nfermi_loc,3))
     allocate(JJm_list(num_wann,num_wann,nfermi_loc,3))
    


    ! Gather W-gauge matrix objects
    !
    !call wham_get_eig_UU_HH_JJlist(kpt,eig,UU,mdum,JJp_list,JJm_list)
    ! occupation matrices f and g=1-f
    !call wham_get_occ_mat_list(eig,UU,f_list,g_list)
    
     if(present(occ)) then
        call wham_get_eig_UU_HH_JJlist(kpt,eig,UU,mdum,JJp_list,JJm_list,occ=occ)
        call wham_get_occ_mat_list(UU,f_list,g_list,occ=occ)
     else
        call wham_get_eig_UU_HH_JJlist(kpt,eig,UU,mdum,JJp_list,JJm_list)
        call wham_get_occ_mat_list(UU,f_list,g_list,eig=eig)
     endif
    

    ! Eqs.(39-40) WYSV06
    !
    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA,OO_pseudo=OOmega)

    ! Trace formula, Eq.(51) LVTS12
    !
    do if=1,nfermi_loc
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

  end subroutine berry_get_imf_klist


   subroutine berry_get_imfgh_klist(kpt,imf_k_list,img_k_list,imh_k_list,occ)
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
   ! in axial-vector form. Optionally, can specify by hand over   !
   ! which bands to trace (overrides list of fermi energies).     !
   !                                                              !
   !==============================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_re_tr,utility_im_tr,utility_rotate
    use w90_parameters, only     : num_wann,nfermi
    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec,pw90common_fourier_R_to_k
    use w90_wan_ham, only        : wham_get_eig_UU_HH_JJlist,wham_get_occ_mat_list
    use w90_get_oper, only       : AA_R,BB_R,CC_R

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), intent(out), dimension(:,:,:) :: imf_k_list
    real(kind=dp), intent(out), dimension(:,:,:) :: img_k_list
    real(kind=dp), intent(out), dimension(:,:,:) :: imh_k_list
    real(kind=dp), intent(in), optional, dimension(:) :: occ

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
    integer                       :: i,j,if,nfermi_loc,n

    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
     if(present(occ)) then
        nfermi_loc=1
     else
        nfermi_loc=nfermi
     endif
     allocate(f_list(num_wann,num_wann,nfermi_loc))
     allocate(g_list(num_wann,num_wann,nfermi_loc))
     allocate(JJp_list(num_wann,num_wann,nfermi_loc,3))
     allocate(JJm_list(num_wann,num_wann,nfermi_loc,3))    
    
    allocate(AA(num_wann,num_wann,3))
    allocate(BB(num_wann,num_wann,3))
    allocate(CC(num_wann,num_wann,3,3))
    allocate(OOmega(num_wann,num_wann,3))
    allocate(LLambda_i(num_wann,num_wann))
    allocate(mdum(num_wann,num_wann))

    ! Gather W-gauge matrix objects
    !
     if(present(occ)) then
        call wham_get_eig_UU_HH_JJlist(kpt,eig,UU,HH,JJp_list,JJm_list,occ=occ)
        call wham_get_occ_mat_list(UU,f_list,g_list,occ=occ)
     else
        call wham_get_eig_UU_HH_JJlist(kpt,eig,UU,HH,JJp_list,JJm_list)
        call wham_get_occ_mat_list(UU,f_list,g_list,eig=eig)
     endif

    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA,OO_pseudo=OOmega) 
    call pw90common_fourier_R_to_k_vec(kpt,BB_R,OO_true=BB)
    do j=1,3
       do i=1,j
          call pw90common_fourier_R_to_k(kpt,CC_R(:,:,:,i,j),CC(:,:,i,j),0)
          CC(:,:,j,i)=conjg(transpose(CC(:,:,i,j)))
       enddo
    enddo

    ! Trace formula for -2Im[f], Eq.(51) LVTS12
    !
    do if=1,nfermi_loc
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
       ! LLambda_ij [Eq. (37) LVTS12] expressed as a pseudovector
       LLambda_i=cmplx_i*(CC(:,:,alpha_A(i),beta_A(i))&
              -conjg(transpose(CC(:,:,alpha_A(i),beta_A(i)))))
       do if=1,nfermi_loc
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
    do if=1,nfermi_loc
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
          mdum=matmul(HH,matmul(JJm_list(:,:,if,alpha_A(i)),&
               JJp_list(:,:,if,beta_A(i))))
          imh_k_list(3,i,if)=-2.0_dp*utility_im_tr(mdum)
          !
       enddo
    enddo

  end subroutine berry_get_imfgh_klist


  !===========================================================!
  !                   PRIVATE PROCEDURES                      ! 
  !===========================================================!

  subroutine berry_get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k,&
                        kubo_H_k_spn,kubo_AH_k_spn,jdos_k_spn)
  !====================================================================!
  !                                                                    !
  ! Contribution from point k to the complex interband optical         !
  ! conductivity, separated into Hermitian (H) and anti-Hermitian (AH) ! 
  ! parts. Also returns the joint density of states                    !
  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_diagonalize,utility_rotate,utility_w0gauss
    use w90_parameters, only     : num_wann,kubo_nfreq,kubo_freq_list,&
                                   fermi_energy_list,kubo_eigval_max,&
                                   kubo_adpt_smr,kubo_smr_fixed_en_width,&
                                   kubo_adpt_smr_max,kubo_adpt_smr_fac,&
                                   kubo_smr_index,berry_kmesh,spin_decomp
    use w90_postw90_common, only : pw90common_get_occ,pw90common_fourier_R_to_k_new,&
                                   pw90common_fourier_R_to_k_vec,pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : HH_R,AA_R
    use w90_spin, only           : spin_get_nk

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
                        eta_smr,Delta_k,arg,vdum(3)
    
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
       call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
       Delta_k=pw90common_kmesh_spacing(berry_kmesh)
    else
       call pw90common_fourier_R_to_k_new(kpt,HH_R,OO=HH,&
                                        OO_dx=delHH(:,:,1),&
                                        OO_dy=delHH(:,:,2),&
                                        OO_dz=delHH(:,:,3))
       call utility_diagonalize(HH,num_wann,eig,UU)
    endif
    call pw90common_get_occ(eig,occ,fermi_energy_list(1))
    call wham_get_D_h(delHH,UU,eig,D_h)

    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
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
       call spin_get_nk(kpt,spn_nk)
       kubo_H_k_spn=cmplx_0
       kubo_AH_k_spn=cmplx_0
       jdos_k_spn=0.0_dp
    end if
    do m=1,num_wann
       do n=1,num_wann
          if( abs(eig(n)-eig(m))<0.0001 ) cycle
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
             vdum(:)=del_eig(m,:)-del_eig(n,:)
             joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
             eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,&
                  kubo_adpt_smr_max)
	  else
	    eta_smr=imag(omega)
          endif
          rfac1=(occ(m)-occ(n))*(eig(m)-eig(n))
          occ_prod=occ(n)*(1.0_dp-occ(m))
          do ifreq=1,kubo_nfreq   
             !
             ! Complex frequency for the anti-Hermitian conductivity
             !
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
             delta=utility_w0gauss(arg,kubo_smr_index)/eta_smr
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

  end subroutine berry_get_kubo_k


  subroutine berry_get_kubo_current_k(kpt,resigma_k,imsigma_k) 
  !====================================================================!
  !                                                                    !
  ! Contribution from point k to the current-induced complex interband !
  ! optical conductivity                                               !
  !                                                                    !
  ! sigma_k(iab,ifreq,ifermi) =                                        !
  !	SUM_{n,m} delta(E_kn-E_f) (d E_{kn}/d k_i)*                    !
  !      [-2*Im(A_nm_a A_mn_b)]*w_mn^2/(w_mn^2-w^2)                    !
  ! [units of length^3]                                                !
  !  where w_mn=E_{km}-E_{kn}                                          !
  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_rotate,utility_w0gauss
    use w90_parameters, only     : num_wann,kubo_nfreq,kubo_freq_list,&
                                   fermi_energy_list,nfermi,kubo_eigval_max,&
                                   kubo_adpt_smr,kubo_smr_fixed_en_width,&
                                   kubo_adpt_smr_max,kubo_adpt_smr_fac,&
                                   kubo_smr_index,kubo_fermi_adpt_smr,&
                                   kubo_fermi_smr_fixed_en_width,&
                                   kubo_fermi_adpt_smr_max,&
                                   kubo_fermi_adpt_smr_fac,&
                                   kubo_fermi_smr_index,&
                                   berry_kmesh,kubo_smr_max_arg,&
                                   berry_band_list,num_berry_bands
    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec,&
                                   pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : AA_R

    ! Arguments
    !
    real(kind=dp),                   intent(in)  :: kpt(3)
    real(kind=dp), dimension(:,:,:), intent(out) :: resigma_k
    real(kind=dp), dimension(:,:,:), intent(out) :: imsigma_k

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    real(kind=dp)                 :: iAAAA(9)
    complex(kind=dp), allocatable :: multW(:)
    real(kind=dp), allocatable    :: multWre(:)
    real(kind=dp), allocatable    :: multWim(:)
     
    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann,3),joint_level_spacing,level_spacing,&
                        eta_smr_fermi,eta_smr,Delta_k,arg,vdum(3)
    
    integer          :: i,j,n,m,n1,m1,ifreq,ifermi
    real(kind=dp)    :: eig(num_wann),deltaF,wmn
    logical          :: AAcalculated,AAAAcalculated

    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(multW(kubo_nfreq))
    allocate(multWre(kubo_nfreq))
    allocate(multWim(kubo_nfreq))
  
    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)

    if(kubo_adpt_smr.or.kubo_fermi_adpt_smr) &
         Delta_k=pw90common_kmesh_spacing(berry_kmesh)
    
    if (.not. kubo_adpt_smr) eta_smr=kubo_smr_fixed_en_width
    if (.not. kubo_fermi_adpt_smr) eta_smr_fermi=kubo_fermi_smr_fixed_en_width
    

    ! Replace imaginary part of frequency with a fixed value
    if(.not.kubo_adpt_smr .and. kubo_smr_fixed_en_width/=0.0_dp)&
         kubo_freq_list=real(kubo_freq_list,dp)+cmplx_i*kubo_smr_fixed_en_width

    resigma_k=0.0_dp
    imsigma_k=0.0_dp
    
    AAcalculated=.false.
    do n1=1,num_berry_bands
       n=berry_band_list(n1)
       if(eig(n)>kubo_eigval_max) cycle
       if(kubo_fermi_adpt_smr) then
          level_spacing=sqrt(dot_product(del_eig(n,:),del_eig(n,:)))*Delta_k
          eta_smr_fermi=min(level_spacing*kubo_fermi_adpt_smr_fac,&
               kubo_fermi_adpt_smr_max)
       endif
       do m1=1,num_berry_bands
    	  m=berry_band_list(m1)
          if(n==m) cycle
          if(eig(m)>kubo_eigval_max) cycle
          
          wmn=eig(m)-eig(n)
          
          
          if(kubo_adpt_smr) then
             ! Eq.(35) YWVS07 
             vdum(:)=del_eig(m,:)-del_eig(n,:)
             joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
             eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,&
                  kubo_adpt_smr_max)
          endif
          
          ! multW(:)=2.0_dp/(wnm**2-kubo_freq_list(:)**2)
          multW(:)=wmn**2/(wmn**2-kubo_freq_list(:)**2) !Ivo
          multWre(:)=real(multW(:))
          multWim(:)=imag(multW(:))
          
          AAAAcalculated=.false.
          do ifermi=1,nfermi
             arg=(eig(n)-fermi_energy_list(ifermi))/eta_smr_fermi
             if(abs(arg)>kubo_smr_max_arg) cycle
             if (.not.AAAAcalculated) then
                if(.not.AAcalculated) then
		    allocate(AA(num_wann,num_wann,3))
		    allocate(D_h(num_wann,num_wann,3))
		    call wham_get_D_h(delHH,UU,eig,D_h)
		    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
		    do i=1,3
    			AA(:,:,i)=utility_rotate(AA(:,:,i),UU,num_wann)
		    enddo
		    AA=AA+cmplx_i*D_h ! Eq.(25) WYSV06
		    AAcalculated=.true.
		endif                
        	do i=1,9
            	    ! iAAAA(i)=imag(AA(n,m,alpha_Axyz(i))*AA(m,n,beta_Axyz(i)))*&
            	    !      del_eig(n,gamma_Axyz(i))
            	    !
            	    ! ---------------------------
            	    ! Ivo: Note added minus sign!
            	    ! ---------------------------
            	    !
            	    iAAAA(i)=-2.0_dp*imag(AA(n,m,alpha_Axyz(i))*AA(m,n,beta_Axyz(i)))&
                	*del_eig(n,gamma_Axyz(i))
        	enddo
		! iAAAA(:)=iAAAA(:)*wnm**2
                AAAAcalculated=.true.
             endif  ! .not.AAAAcalculated
             deltaF=utility_w0gauss(arg,kubo_fermi_smr_index)/eta_smr_fermi
             call dger(9,kubo_nfreq,deltaF,iAAAA,1,multWim,1,&
                  imsigma_k(:,:,ifermi),9)
             call dger(9,kubo_nfreq,deltaF,iAAAA,1,multWre,1,&
                  resigma_k(:,:,ifermi),9)
          enddo  ! ifermi

       enddo !m
    enddo !n

  end subroutine berry_get_kubo_current_k


  subroutine berry_get_NOA_Bnl_orb(eig,del_eig,AA,&
	num_occ,occ_list,num_unocc,unocc_list,Bnl) 
  !====================================================================!
  !                                                                    !
  ! Calculating the matrix                                             !
  ! 2*hbar*B_{nl,ac}(num_occ,num_unocc,3,3)=                    !
  !      -sum_m(  (E_m-E_n)A_nma*Amlc +(E_l-E_m)A_nmc*A_mla -          !
  !      -i( del_a (E_n+E_l) A_nlc                                     !
  !   in units eV*Ang^2                                                !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_i,cmplx_0
    use w90_parameters, only     : num_berry_bands,berry_band_list
                                   
    use w90_comms, only          : on_root
    use w90_io, only          : stdout,io_time,io_error

    ! Arguments
    !
    integer,                              intent(in) ::num_occ,num_unocc
    integer,          dimension(:),       intent(in) ::occ_list , unocc_list
    real(kind=dp),    dimension(:),       intent(in) ::eig     !  n
    real(kind=dp),    dimension(:,:),     intent(in) ::del_eig !  n
    complex(kind=dp), dimension(:,:,:),   intent(in) ::AA      !  n,l,a
    complex(kind=dp), dimension(:,:,:,:), intent(out)::Bnl     !   n,l,a,c
    integer n,m,l,a,c,n1,m1,l1

    Bnl(:,:,:,:)=cmplx_0
    
    
    do a=1,3
    do c=1,3
	do n1=1,num_occ
	    n=occ_list(n1)
	    do l1=1,num_unocc
		l=unocc_list(l1)
		Bnl(n1,l1,a,c)=-cmplx_i*(del_eig(n,a)+del_eig(l,a))*AA(n,l,c)
		do m1=1,num_berry_bands
		    m=berry_band_list(m1)
		    Bnl(n1,l1,a,c)=Bnl(n1,l1,a,c)+ &
		        (eig(n)-eig(m))*AA(n,m,a)*AA(m,l,c)- &
			(eig(l)-eig(m))*AA(n,m,c)*AA(m,l,a)
		enddo ! m1
	    enddo !l1
	enddo !n1
    enddo !c
    enddo !a

  end subroutine berry_get_NOA_Bnl_orb


  subroutine berry_get_NOA_Bnl_spin( S_h,&
	num_occ,occ_list,num_unocc,unocc_list,Bnl) 
  !====================================================================!
  !                                                                    !
  ! Calculating the matrix                                             !
  ! 2*hbar*B_{nl,ac}^spin(num_occ,num_unocc,3,3)=                    !
  !      -i hbar^2/m_e eps_{abc}  < u_n | sigma_b | u_l >
  !   in units eV*Ang^2                                                !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_i,cmplx_0,elec_mass_SI,hbar_SI,elem_charge_SI
                                   
    use w90_comms, only          : on_root
    use w90_io, only          : stdout,io_time,io_error

    ! Arguments
    !
    integer,                              intent(in) ::num_occ,num_unocc
    integer,          dimension(:),       intent(in) ::occ_list , unocc_list
    complex(kind=dp), dimension(:,:,:),   intent(in) ::S_h     !  n,l,a
    complex(kind=dp), dimension(:,:,:,:), intent(out)::Bnl     !   n,l,a,c
    integer n,l,a,b,c,n1,l1

    Bnl(:,:,:,:)=cmplx_0
    
    
    do b=1,3
!	if (on_root) write(stdout,*) b
	c=alpha_A(b)
	a=beta_A(b)
!	if (on_root) write(stdout,*) a,b,c
	do n1=1,num_occ
	    n=occ_list(n1)
	    do l1=1,num_unocc
		l=unocc_list(l1)
!		if (on_root) write(stdout,*) a,b,c,n1,n,l1,l
		Bnl(n1,l1,a,c)=S_h(n,l,b)
	    enddo !l1
	enddo !n1
    enddo !b
    
    Bnl=Bnl*(-cmplx_i)*(1e+20*hbar_SI**2/(elec_mass_SI*elem_charge_SI))

  end subroutine berry_get_NOA_Bnl_spin


  subroutine berry_get_NOA_sigma_k(kpt,sigma_e_k,sigma_m_orb_k,sigma_m_spin_k) 
  !====================================================================!
  !                                                                    !
  ! Contribution from point k to the real (antisymmetric) part         !
  ! of the natural  complex interband optical conductivity             !
  !                                                                    !
  ! Re sigma_m  =  SUM_{n,l}^{oe}  w/(w_nl^2-w^2) *  
  !   Re (  A_lnb Bnlac -Alna Bnlbc)
  ! [units of Ang^3]                                                !
  !
  ! Re sigma_e  =  -SUM_{n,l}^{oe}  w*(w_ln^2-w^2)/(w_nl^2-w^2)^2 *  
  ! Im (  A_lnb Bnlac -Alna Bnlac)nm_a A_mn_b )
  ! [units of Ang^3]                                                !
  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_rotate,utility_w0gauss,utility_wgauss
    use w90_parameters, only     : num_wann,kubo_nfreq,kubo_freq_list,kubo_smr_fixed_en_width,&
                                   fermi_energy_list,nfermi,kubo_eigval_max,&
                                   num_berry_bands,berry_band_list
                                   
    use w90_comms, only          : on_root
    use w90_io, only          : stdout,io_time,io_error

    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec,&
				   pw90common_fourier_R_to_k_new,&
                                   pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : AA_R,SS_R
    use w90_spin, only           : get_S

    ! Arguments
    !
    real(kind=dp),                     intent(in)  :: kpt(3)
    real(kind=dp), dimension(:,:,:), intent(out) :: sigma_m_orb_k,sigma_e_k ! ifreq,ab,c
    real(kind=dp), dimension(:,:,:),optional, intent(out) :: sigma_m_spin_k ! ifreq,ab,c


    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: Bnl_orb(:,:,:,:)
    complex(kind=dp), allocatable :: Bnl_spin(:,:,:,:)
    complex(kind=dp), allocatable :: SS(:,:,:)
    complex(kind=dp), allocatable :: S_h(:,:,:)
    complex(kind=dp)              :: W(kubo_nfreq)
!    real(kind=dp)              :: W(kubo_nfreq)
    
    real(kind=dp) :: multWe(kubo_nfreq),multWm(kubo_nfreq)

    integer ::  num_occ,num_unocc,occ_list(num_wann),unocc_list(num_wann)
     
    real(kind=dp)    :: del_eig(num_wann,3),joint_level_spacing,level_spacing,&
                    	time0,time1,time2,wmn
    
    integer          :: i,j,n,l,n1,l1,a,b,c,ab
    real(kind=dp)    :: eig(num_wann),wln
    logical :: deltaWzero,AAcalculated
    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(SS(num_wann,num_wann,3))
    allocate(S_h(num_wann,num_wann,3))

    if (on_root) write(stdout,*) "kpt=",kpt

    W=real(kubo_freq_list(:))+cmplx_i*kubo_smr_fixed_en_width
!    W=real(kubo_freq_list(:)) ! +cmplx_i*kubo_smr_fixed_en_width

    if (on_root) write(stdout,*) "W[10]=",W(10)

    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
  
    num_occ=0
    num_unocc=0
    do n=1,num_berry_bands
	if (eig(berry_band_list(n))<fermi_energy_list(1) ) then
	    num_occ=num_occ+1
	    occ_list(num_occ)=berry_band_list(n)
	elseif (eig(berry_band_list(n))<kubo_eigval_max) then
	    num_unocc=num_unocc+1
	    unocc_list(num_unocc)=berry_band_list(n)
	endif
    enddo

    if (num_occ==0) then 
	call io_error("no occupied bands included in the calculation")
    endif

    if (num_unocc==0) then 
	call io_error("no unoccupied bands included in the calculation")
    endif

!    if (on_root) write(stdout,*) "occupied bands : ",(occ_list(n),n=1,num_occ)
!    if (on_root) write(stdout,*) "unoccupied bands : ",(unocc_list(n),n=1,num_unocc)

    allocate(D_h(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))
    call wham_get_D_h(delHH,UU,eig,D_h)
    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
    do i=1,3
        AA(:,:,i)=utility_rotate(AA(:,:,i),UU,num_wann)
    enddo
    AA=AA+cmplx_i*D_h ! Eq.(25) WYSV06

!    if (on_root) write(stdout,*) "AA calculated "

    allocate (Bnl_orb(num_occ,num_unocc,3,3))    
    call berry_get_NOA_Bnl_orb(eig,del_eig,AA,num_occ,occ_list,num_unocc,unocc_list,Bnl_orb)
    sigma_m_orb_k=cmplx_0
    sigma_e_k=cmplx_0


    if (present(sigma_m_spin_k)) then
       do j=1,3 ! spin direction
          call pw90common_fourier_R_to_k_new(kpt,SS_R(:,:,:,j),OO=SS(:,:,j))
          S_h(:,:,j)=utility_rotate(SS(:,:,j),UU,num_wann)
       enddo
	allocate (Bnl_spin(num_occ,num_unocc,3,3))
	call berry_get_NOA_Bnl_spin(S_h,num_occ,occ_list,num_unocc,unocc_list,Bnl_spin)
	sigma_m_spin_k=cmplx_0
    endif




    
    do n1=1,num_occ
	n=occ_list(n1)
	do l1=1,num_unocc
	    l=unocc_list(l1)

!	    if (on_root) write(stdout,*) n1,l1,n,l

            wln=eig(l)-eig(n)

	    multWm=real(W(:))*real(1./(wln*wln-W(:)*W(:)))
	    multWe=real(-multWm*(3*wln*wln-W(:)*W(:))/(wln*wln-W(:)*W(:)))
	    do ab=1,3
	      a=alpha_A(ab)
	      b=beta_A(ab)
	      do c=1,3
        	sigma_m_orb_k(ab,c,:)=sigma_m_orb_k(ab,c,:)+&
        	    multWm*real(AA(l,n,b)*Bnl_orb(n1,l1,a,c)-AA(l,n,a)*Bnl_orb(n1,l1,b,c))
        	if (present(sigma_m_spin_k)) then
        	    sigma_m_spin_k(ab,c,:)=sigma_m_spin_k(ab,c,:)+&
        		multWm*real(AA(l,n,b)*Bnl_spin(n1,l1,a,c)-AA(l,n,a)*Bnl_spin(n1,l1,b,c))
        	endif
        	sigma_e_k(ab,c,:)=sigma_e_k(ab,c,:)+&
        	    multWe*(del_eig(n,c)+del_eig(l,c))*imag(AA(n,l,a)*AA(l,n,b))
    	      enddo
    	    enddo
    	enddo  ! l1
    enddo ! n1

  end subroutine berry_get_NOA_sigma_k


  subroutine berry_get_CPGE_k(kpt,beta_CPGE_k,two_band_model) 
  !====================================================================!
  !                                                                    !
  ! Contribution from point k to the current-induced                   !
  ! complex interband optical conductivity                             !
  !                                                                    !
  ! beta_CPGE_k(ij,ifreq,ifermi) =                		       !
  !	SUM_{n,m,a,b} eps_{jab} (f_nk-f_mk) (d E_{knm}/d k_i)          !
  !        A_nm_a A_mn_b  *  delta(W-Emn)                              !
  ! [units of Ang^3]                                                !
  !  OR in a two-band model
  ! beta_CPGE_k(ij,ifreq,ifermi) =                			!
  !	(f_1k-f_2k) (d E_{k21}/d k_i)           !
  !        Omaga^j_1  *  delta(W-E21)                              !

  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_rotate,utility_w0gauss,utility_wgauss
    use w90_parameters, only     : num_wann,kubo_nfreq,kubo_freq_list,&
                                   fermi_energy_list,nfermi,kubo_eigval_max,&
                                   kubo_adpt_smr,kubo_smr_fixed_en_width,&
                                   kubo_adpt_smr_max,kubo_adpt_smr_fac,&
                                   kubo_smr_index,&
                                   kubo_fermi_adpt_smr,&
                                   kubo_fermi_smr_fixed_en_width,&
                                   kubo_fermi_adpt_smr_max,&
                                   kubo_fermi_adpt_smr_fac,&
                                   kubo_fermi_smr_index,&
                                   berry_kmesh,kubo_smr_max_arg,&
                                   num_berry_bands,berry_band_list
                                   
    use w90_comms, only          : on_root
    use w90_io, only          : stdout,io_time,io_error

    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec,&
                                   pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : AA_R

    ! Arguments
    !
    real(kind=dp),                                  intent(in)  :: kpt(3)
    real(kind=dp),           dimension(:,:,:),   intent(out) :: beta_CPGE_k   !   iab,ifreq,ifermi
    logical,    optional       ,   intent(in) :: two_band_model   !   iab,ifreq,ifermi

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    real(kind=dp) :: iAAAA(9)
    real(kind=dp) :: iAAA(3) 
    real(kind=dp), allocatable :: deltaW(:)
    integer :: num_bands
    integer,allocatable :: band_list_n(:),band_list_m(:)
     
    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann,3),joint_level_spacing,level_spacing,&
                        eta_smr_fermi_n,eta_smr_fermi_m,eta_smr,Delta_k,arg,vdum(3),& !enmax,enmin,&
                    	time0,time1,time2,wmn,fermifactor
    
    integer          :: i,j,ij,n,m,n1,m1,ifreq,ifermi
    real(kind=dp)    :: eig(num_wann),deltaF,omega
    logical :: deltaWzero,AAcalculated,twoband
    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(deltaW(kubo_nfreq))
  
!        write(*,*) "num_berry_bands= ",num_berry_bands
!       write(*,*) "berry bands are : ",(berry_band_list(i) , i=1,num_berry_bands)

    twoband=.false.
    if (present(two_band_model)) twoband=two_band_model
    if (two_band_model) then
       if (num_berry_bands.ne.2) call io_error("Exaxctly 2 bands should be provided for a two-band model of CPGE")
       allocate(band_list_n(1))
       allocate(band_list_m(1))
       band_list_n(1)=berry_band_list(1)
       band_list_m(1)=berry_band_list(2)
       num_bands=1
    else
       allocate(band_list_n(num_berry_bands))
       allocate(band_list_m(num_berry_bands))
       num_bands=num_berry_bands
       band_list_n(:)=berry_band_list(:)
       band_list_m(:)=berry_band_list(:)
    endif
    

    
    beta_CPGE_k=0.0_dp
    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)

    if(kubo_adpt_smr.or.kubo_fermi_adpt_smr) \
       Delta_k=pw90common_kmesh_spacing(berry_kmesh)
    
    if (.not. kubo_adpt_smr )  eta_smr=kubo_smr_fixed_en_width
    if (.not. kubo_fermi_adpt_smr ) then
	 eta_smr_fermi_n=kubo_fermi_smr_fixed_en_width
	 eta_smr_fermi_m=kubo_fermi_smr_fixed_en_width
    endif
    
    AAcalculated=.false.

    do n1=1,num_bands
       n=band_list_n(n1)
       if(eig(n)>kubo_eigval_max)  cycle
       if(kubo_fermi_adpt_smr) then
            level_spacing=sqrt(dot_product(del_eig(n,:),del_eig(n,:)))*Delta_k
            eta_smr_fermi_n=min(level_spacing*kubo_fermi_adpt_smr_fac,&
                 kubo_fermi_adpt_smr_max)
         endif
       do m1=1,num_bands
          m=band_list_m(m1) 
          if(n==m) cycle
          if(eig(m)>kubo_eigval_max)  cycle
!          write(*,*) "bands:",n,m
          if(kubo_adpt_smr) then
             ! Eq.(35) YWVS07 
             vdum(:)=del_eig(m,:)-del_eig(n,:)
             joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
             if(kubo_adpt_smr) eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,kubo_adpt_smr_max)
          endif
          if(kubo_fermi_adpt_smr) then
              level_spacing=sqrt(dot_product(del_eig(m,:),del_eig(m,:)))*Delta_k
              eta_smr_fermi_n=min(level_spacing*kubo_fermi_adpt_smr_fac,&
                    kubo_fermi_adpt_smr_max)
          endif

          wmn=eig(m)-eig(n)
!	  write(*,*) "calculating deltaW for wmn=",wmn,"eta_smr=",eta_smr
          deltaWzero=.true.
	  deltaW=0.0_dp
          do i=1,kubo_nfreq
             arg=(wmn-real(kubo_freq_list(i)))/eta_smr
             if (abs(arg)>kubo_smr_max_arg) cycle 
	     deltaWzero=.false.
	     deltaW(i)=utility_w0gauss(arg,kubo_smr_index)/eta_smr
          enddo
	  if (deltaWzero) cycle
!	  write(*,*) "calculating iAAA"

	if (twoband) then
		call berry_get_berry_curv_nk(kpt,n,iAAA)
		iAAA(:)=-iAAA(:)
	else
	    if (.not.AAcalculated) then 
		allocate(D_h(num_wann,num_wann,3))
		allocate(AA(num_wann,num_wann,3))
		call wham_get_D_h(delHH,UU,eig,D_h)
		call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
		do i=1,3
    		AA(:,:,i)=utility_rotate(AA(:,:,i),UU,num_wann)
		enddo
		AA=AA+cmplx_i*D_h ! Eq.(25) WYSV06
		AAcalculated=.true.
	    endif
    	    do j=1,3
                iAAA(j)=imag(AA(n,m,alpha_A(j))*AA(m,n,beta_A(j))-AA(n,m,beta_A(j))*AA(m,n,alpha_A(j)))
    	    enddo
	endif
	  
	  do ij=1,9
              iAAAA(ij)=(del_eig(n,alpha_xyz(ij))-del_eig(m,alpha_xyz(ij))) * iAAA(beta_xyz(ij))
	  enddo
          
          do ifermi=1,nfermi 
	     fermifactor= & 
	        utility_wgauss((fermi_energy_list(ifermi)-eig(n))/eta_smr_fermi_n,kubo_fermi_smr_index) &
	          -utility_wgauss((fermi_energy_list(ifermi)-eig(m))/eta_smr_fermi_m,kubo_fermi_smr_index)
!	     write(*,*) "fermi factor for Efermi",fermi_energy_list(ifermi),"E=",eig(n),eig(m)," smearing ",eta_smr_fermi_n,eta_smr_fermi_m," i ",fermifactor
	     if (abs(fermifactor)<1e-5) cycle
!	     write(*,*) "calling dger!"
	     call  dger(9,kubo_nfreq,fermifactor,iAAAA,1,deltaW, 1 , &
			    beta_CPGE_k(:,:,ifermi) , 9 )
           enddo ! ifermi 
       enddo
    enddo

  end subroutine berry_get_CPGE_k





  subroutine berry_sphere_flux() 
  !====================================================================!
  !                                                                    !
  ! 
  ! calculates the fluxes through a sphere around a k-point
  !	Omega_n^j=  i*eps^{jab}  A_nm_a A_mn_b                         !
  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_rotate,utility_w0gauss,utility_wgauss
    use w90_parameters, only     : num_wann,&
                                   num_berry_bands,berry_band_list,&
                            	   berry_sphere_rad,berry_sphere_nk, berry_sphere_center,&
                            	   berry_sphere_band_list,num_berry_sphere_bands,real_lattice
                                   
    use w90_comms, only          : on_root
    use w90_io, only             : stdout,io_time

    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec,&
                                   pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : AA_R

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    real(kind=dp), allocatable :: fluxes(:)
     
    real(kind=dp)    :: time0,time1,time2
    
    integer          :: i,j,ij,n,m,n1,m1,it,iphi,nt,nphi
    real(kind=dp)    :: eig(num_wann),dk(3),dk_frac(3),kpt(3),weight,weightsum,t,phi,dt,dphi,&
			del_eig(num_wann,3)
    
    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(D_h(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(fluxes(num_berry_sphere_bands))

    nt=berry_sphere_nk+1
    dt=pi/(nt-1)
    fluxes=0
    weightsum=0

	write(stdout,*) "berry_sphere_band_list(",num_berry_sphere_bands,")=",berry_sphere_band_list
	write(stdout,*) "berry_band_list(",num_berry_bands,")=",berry_band_list


    do it=1,nt
	write(stdout,*) "it=",it," of ",nt
	t=dt*(it-1)
	if (it==1 .or. it==nt) then
	    nphi=1
	    weight=2*pi*(1-cos(dt/2))
	else
	    nphi=NINT(berry_sphere_nk*2*sin(t))
	    weight=2*pi*(cos(dt*(it-1.5))-cos(dt*(it-0.5)))/nphi
	endif
	
	do iphi=1,nphi
!	    write(stdout,*) "iphi=",iphi," of ",nphi
	    phi=2*pi*iphi/nphi
	    dk(1)=cos(phi)*sin(t) 
	    dk(2)=sin(phi)*sin(t) 
	    dk(3)=cos(t)          
	    dk_frac=matmul(real_lattice(:,:),dk(:))/(2.*pi)
	    kpt=berry_sphere_center(:)+dk_frac(:)*berry_sphere_rad
!	    write (stdout,*) "kpt=",kpt
	    weightsum=weightsum+weight
!	    write(stdout,*) " 1 - OK"
	    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
!	    write(stdout,*) " 2 - OK"
	    call wham_get_D_h(delHH,UU,eig,D_h)
!	    write(stdout,*) " 3 - OK"
	    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
!	    write(stdout,*) " 4 - OK"
	    do i=1,3
    		AA(:,:,i)=utility_rotate(AA(:,:,i),UU,num_wann)
!		write(stdout,*) " 5 - OK",i,"/3"
	    enddo
	    AA=AA+cmplx_i*D_h ! Eq.(25) WYSV06
!	    write(stdout,*) " AA - OK"
	    do n1=1,num_berry_sphere_bands
		n=berry_sphere_band_list(n1)
    		do m1=1,num_berry_bands
        	    m=berry_band_list(m1) 
        	    if(n==m) cycle
        	    do j=1,3
        	        fluxes(n1)=fluxes(n1)-imag(AA(n,m,alpha_A(j))*AA(m,n,beta_A(j))-AA(n,m,beta_A(j))*AA(m,n,alpha_A(j)))*dk(j)*weight
        	    enddo
        	enddo !m1
    	    enddo ! n1
    	enddo   !iphi
    enddo !it

    fluxes=fluxes*berry_sphere_rad**2/(2*pi)
    
    write(stdout,*) "Fluxes through the sphere of center ",berry_sphere_center," and radius ",berry_sphere_rad
    do n1=1,num_berry_sphere_bands
	n=berry_sphere_band_list(n1)
	write(stdout,'(a5,I3,a12,F9.6)') "band ",n," flx/2pi = ",fluxes(n1)
    enddo
    write(stdout,*) "---------------------------------------------------"
    write(stdout,'(a,F6.4)') "total integrated 3D angle : 4pi*",weightsum/(4.*pi)
  end subroutine berry_sphere_flux



  subroutine berry_get_berry_curv_nk(kpt,n,Omega)
    
    use w90_parameters, only     : num_wann

    real(kind=dp), intent(in)                    :: kpt(3)
    integer      , intent(in) 			 :: n
    real(kind=dp), dimension(:), intent(out) :: Omega

    real(kind=dp)    ::  imf_k(3,3,1),occ(num_wann)
    integer i


    occ=0.0_dp
    occ(n)=1.0_dp
    call berry_get_imf_klist(kpt,imf_k,occ)
    do i=1,3
        Omega(i)=sum(imf_k(:,i,1))
    enddo

  end subroutine berry_get_berry_curv_nk
    

  subroutine berry_get_gme_k_list(kpt,gme_spn_k_list,gme_orb_k_list,&
       gme_ber_k_list,gme_berE_k_list,gme_spn_k_list_alt,gme_ohm_k_list)
  !======================================================================!
  !                                                                      !
  ! Contribution from point k to the GME tensor, Eq.(9) of ZMS16,        !
  ! evaluated in the clean limit of omega.tau >> 1 where  it is real.    !
  ! The following two quantities are calculated (sigma = Pauli matrix):  !
  !                                                                      !
  ! gme_spn_k = delta(E_kn-E_f).(d E_{kn}/d k_i).sigma_{kn,j}            !
  ! [units of length]                                                    !
  !                                                                      !
  ! gme_orb_k = delta(E_kn-E_f).(d E_{kn}/d k_i).(2.hbar/e).m^orb_{kn,j} !
  ! [units of (length^3)*energy]                                         !
  !                                                                      !
  ! gme_ber_k = delta(E_kn-E_f).(d E_{kn}/d k_i).Omega_{kn,j}            !
  ! [units of length^3]                                                  ! 
  !                                                                      !
  ! gme_berE_k = delta(E_kn-E_f).(d E_{kn}/d k_i).Omega_{kn,j}*E_{kn}    !
  ! [units of Energy*length^3]                                                  ! 
  !                                                                      !
  ! gme_ohm_k = delta(E_kn-E_f).(d E_{kn}/d k_i).(d E_{kn}/d k_j)        !
  ! [units of energy*length^2]                                           ! 
  !                                                                      !
  !======================================================================!

    use w90_constants, only      : dp,cmplx_0
    use w90_utility, only        : utility_rotate,utility_rotate_diag,utility_w0gauss
    use w90_parameters, only     : num_wann,fermi_energy_list,&
                                   kubo_smr_fixed_en_width,kubo_adpt_smr,&
                                   kubo_adpt_smr_max,kubo_adpt_smr_fac,&
                                   kubo_smr_index,berry_kmesh,nfermi,gme_task,&
                                   gme_degen_thresh ,kubo_smr_max_arg,&  !  Tsirkin
                                   berry_band_list,num_berry_bands
    use w90_postw90_common, only : pw90common_kmesh_spacing,pw90common_get_occ,&
				    pw90common_fourier_R_to_k_new
    use w90_wan_ham, only        : wham_get_eig_deleig,wham_get_D_h
    use w90_get_oper, only       : HH_R,SS_R
    use w90_spin, only           : get_S,get_del_S
    use w90_io, only             : stdout

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), dimension(:,:,:), intent(out) :: gme_spn_k_list,&
                                                    gme_orb_k_list,&
                                                    gme_ber_k_list,&
                                                    gme_berE_k_list,&
                                                    gme_ohm_k_list
    real(kind=dp), dimension(:,:,:), intent(out), optional :: gme_spn_k_list_alt

    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: SS(:,:,:)
    complex(kind=dp), allocatable :: delSS(:,:,:,:)
    complex(kind=dp), allocatable :: S_h(:,:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)

    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann,3),level_spacing,&
                        eta_smr,Delta_k,arg,rvdum(3)

    integer          :: i,j,n,n1,m1,m,ifermi
    real(kind=dp)    :: delta,eig(num_wann),occ(num_wann),&
                        S(num_wann,3),del_S(num_wann,3,3),&
                        orb_nk(3),ber_nk(3),&
                        imf_k(3,3,1),img_k(3,3,1),imh_k(3,3,1)
    logical          :: got_spin,got_orb_n

    allocate(UU(num_wann,num_wann))
    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(SS(num_wann,num_wann,3))
    allocate(delSS(num_wann,num_wann,3,3))
    allocate(S_h(num_wann,num_wann,3))
    allocate(D_h(num_wann,num_wann,3))

    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
    if(kubo_adpt_smr) Delta_k=pw90common_kmesh_spacing(berry_kmesh)

    gme_spn_k_list=0.0_dp
    gme_orb_k_list=0.0_dp
    gme_ber_k_list=0.0_dp
    gme_berE_k_list=0.0_dp
    gme_ohm_k_list=0.0_dp
    got_spin=.false.
    do n1=1,num_berry_bands
	n=berry_band_list(n1)
       !
       ! ***ADJUSTABLE PARAMETER***
       !
       ! avoid degeneracies
       !---------------------------------------------------
        if (n>1) then
    	    if (eig(n)-eig(n-1) <= gme_degen_thresh) cycle
    	endif
        if (n<num_wann) then
    	    if (eig(n+1)-eig(n) <= gme_degen_thresh) cycle
    	endif
!       if(min(eig(n+1)-eig(n),eig(n)-eig(n-1))<1.e-4) cycle 
       !---------------------------------------------------
       got_orb_n=.false.
       if(kubo_adpt_smr) then
          !
          ! Eq.(35) YWVS07 
          !
          level_spacing=sqrt(dot_product(del_eig(n,:),del_eig(n,:)))*Delta_k
          eta_smr=min(level_spacing*kubo_adpt_smr_fac,kubo_adpt_smr_max)
       else
          eta_smr=kubo_smr_fixed_en_width
       endif
       do ifermi=1,nfermi
          arg=(eig(n)-fermi_energy_list(ifermi))/eta_smr
          !
          ! To save time: far from the Fermi surface, negligible contribution
          !
          !-------------------------
          if(abs(arg)>kubo_smr_max_arg) cycle 
          !-------------------------
          !
          ! Spin is computed for all bands simultaneously
          !
          if(index(gme_task,'spin')>0 .and. .not.got_spin) then
             call get_S(kpt,S)
             got_spin=.true. ! Do it for only one value of ifermi and n
          endif
          !
          ! Orbital quantities are computed for each band separately
          !
          if (.not.got_orb_n) then
             ! Tsirkin  -  do not calculate orbital quantities, when not needed (spin-only or Berry-curvature-only mode)
             if(index(gme_task,'orb')>0)  then   
             !
             ! Fake occupations: band n occupied, others empty
             !
                occ=0.0_dp
                occ(n)=1.0_dp
                call berry_get_imfgh_klist(kpt,imf_k,img_k,imh_k,occ)
                do i=1,3
                   orb_nk(i)=sum(imh_k(:,i,1))-sum(img_k(:,i,1))
                   ber_nk(i)=sum(imf_k(:,i,1))
                enddo
                got_orb_n=.true. ! Do it for only one value of ifermi
             else if(index(gme_task,'berry')>0)  then
                occ=0.0_dp
                occ(n)=1.0_dp
                call berry_get_imf_klist(kpt,imf_k,occ)
                do i=1,3
                    ber_nk(i)=sum(imf_k(:,i,1))
                enddo
                got_orb_n=.true. ! Do it for only one value of ifermi
             endif
          endif
          !
          delta=utility_w0gauss(arg,kubo_smr_index)/eta_smr ! Broadened delta(E_nk-E_f)
          !
          ! Loop over Cartesian tensor components 
          !
          do j=1,3
             do i=1,3
                !
                ! Spin GME tensor, integrand of Eq.(9) ZMS16
                !
                if(index(gme_task,'spin')>0) gme_spn_k_list(i,j,ifermi)=&
                     gme_spn_k_list(i,j,ifermi)+del_eig(n,i)*S(n,j)*delta
                !
                ! Orbital GME tensor, integrand of Eq.(9) ZMS16
                !
                if(index(gme_task,'orb')>0) gme_orb_k_list(i,j,ifermi)=&
                     gme_orb_k_list(i,j,ifermi)+del_eig(n,i)*orb_nk(j)*delta
                !
                ! PGE tensor, integrand of Eq.(8) SF15 
                !
                if(index(gme_task,'berry')>0) then
            	    gme_ber_k_list(i,j,ifermi)=&
                     gme_ber_k_list(i,j,ifermi)+del_eig(n,i)*ber_nk(j)*delta
            	    gme_berE_k_list(i,j,ifermi)=&
                     gme_berE_k_list(i,j,ifermi)+del_eig(n,i)*ber_nk(j)*delta*eig(n)                     
                endif
                if(index(gme_task,'ohmic')>0) gme_ohm_k_list(i,j,ifermi)=&
                     gme_ohm_k_list(i,j,ifermi)+del_eig(n,i)*del_eig(n,j)*delta
             enddo
          enddo
       enddo !ifermi
    enddo !n

    ! Recalculate spin GME using a Fermi-sea integration
    !
    if(index(gme_task,'spin')>0 .and. present(gme_spn_k_list_alt)) then
       gme_spn_k_list_alt(:,:,:)=0.0_dp
       !
       ! OLD VERSION: Eqs.(1)+(2) in notes of 2016-05-12 (Ivo)
       !
       ! call get_del_S(kpt,del_S)
       ! do ifermi=1,nfermi
       !    call pw90common_get_occ(eig,occ,fermi_energy_list(ifermi))
       !    do n=1,num_wann
       !       gme_spn_k_list_alt(:,:,ifermi)=gme_spn_k_list_alt(:,:,ifermi)&
       !            +occ(n)*del_S(n,:,:)
       !    enddo
       ! enddo
       !
       ! NEW VERSION: Eq.(5) in notes of 2016-05-12 (Ivo)
       !
       call wham_get_D_h(delHH,UU,eig,D_h)
       do j=1,3 ! spin direction
          call pw90common_fourier_R_to_k_new(kpt,SS_R(:,:,:,j),OO=SS(:,:,j),&
                                                    OO_dx=delSS(:,:,1,j),&
                                                    OO_dy=delSS(:,:,2,j),&
                                                    OO_dz=delSS(:,:,3,j))
          ! spin matrix in H-gauge
          S_h(:,:,j)=utility_rotate(SS(:,:,j),UU,num_wann)
       enddo
       do j=1,3 ! spin direction
          do i=1,3 ! k-derivative direction
             ! U^dagger.[del_i.sigma^(W)_j].U
             del_S(:,i,j)=&
                  real(utility_rotate_diag(delSS(:,:,i,j),UU,num_wann),dp)
          enddo
       enddo

       do ifermi=1,nfermi
          call pw90common_get_occ(eig,occ,fermi_energy_list(ifermi))
          do n1=1,num_berry_bands
             n=berry_band_list(n1)
             ! 1st term in Eq.(5) of 2016-05-12 notes
             gme_spn_k_list_alt(:,:,ifermi)=gme_spn_k_list_alt(:,:,ifermi)&
                  +occ(n)*del_S(n,:,:)
             do m1=1,num_berry_bands
                m=berry_band_list(m1)
                !------------------------------------
                if(abs(occ(n)-occ(m))<1.e-7_dp) cycle ! not needed?
                !------------------------------------
                do j=1,3
                   do i=1,3
                      ! 2nd term in Eq.(5) of 2016-05-12 notes
                      gme_spn_k_list_alt(i,j,ifermi)=&
                           gme_spn_k_list_alt(i,j,ifermi)&
                           +(occ(n)-occ(m))*S_h(n,m,j)*D_h(m,n,i)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
       
  end subroutine berry_get_gme_k_list




  subroutine berry_eval_wp_search

	!--------------------------------------------------------------
	!          
	! this subroutine searches for Weyl points  --------------------
	!                   
	!-----------------------------------------------------------

    use w90_comms, only          : on_root
    use w90_io, only             : stdout,io_file_unit

    ! Added Jan 2016 (DGM)      
    use w90_parameters, only    : wp_gap_thresh,kslice_coor,recip_lattice,&
                                  wp_det_thresh,wp_dis_thresh,wp_corner,&
                                  kslice_num_slices,kslice_k3min,kslice_k3max

    ! new
    use w90_parameters, only    : wp_box_corner,wp_box_b1,wp_box_b2,&
                                  wp_box_b3,wp_box_k3min,wp_box_k3max,&
                                  wp_box_kmesh,wp_num_filled_bands
    !! new way to get the double and triple wp
    use w90_parameters, only    : num_c4_axis,c4_axis,num_c6_axis,c6_axis
    use w90_wan_ham, only       : get_kdotp_2d,get_kvec_wp
    use w90_get_oper,only       : get_hh_r
    use w90_utility,only        : utility_frac_to_cart


    real(kind=dp), allocatable :: wp_pt(:,:),E0_wp(:),gap_wp(:),deter_wp(:),modulus(:)
    integer, allocatable :: charge_wp(:)

    real(kind=dp)     :: kpt(3),kpt_min(3),kpt_gap,deter,kpt_cart(3),vec_dum(3),& ! wpsearch parameters ,
                         kslice_corner_cart(3),b1(3),b2(3),b3(3),curv_kp(3),&
                         kpoint(3),k_wp(3),energy_wp,&
                         wp_corner_cart(3)
    real(kind=dp)     :: zvec(3),yvec(3),b1mod,b2mod,cosb1b2,cosyb2,kpt_x,&
                         kpt_y,ymod,kaux(3),kslice3_min,kslice3_max,k3_aux

    real(kind=dp)     :: wp_box_db1,wp_box_db2,wp_box_db3
    real(kind=dp)     :: b3_box(3),wp_box_corner_cart(3)
    real(kind=dp)     :: wp_k3_list(0:wp_box_kmesh(3)) ! the third component of the wp_box_kmesh define the wp_list

    real(kind=dp)     :: J_vec(0:3),J_mat(3,0:3) ! k.p parameters

    integer           :: nslices3  ! wpsearch parameters 
    integer           :: counter_wp,chivato  ! wpsearch parameters
    logical           :: inbox !wpsearch
    logical           :: multiple_wp

    integer           :: num_walkers
    integer           :: file_unit,file_unit_wp
    character(len=24) :: file_name,file_name_wp

    integer           ::   i,n,i1,i2,i3
    real(kind=dp)     :: k1,k2,k3
     if(on_root) then
          write(stdout,'(/,/,1x,a,I4,a,I4,/)') '* Steepest descent search of Weyl points &
                           between bands:',wp_num_filled_bands,',',wp_num_filled_bands+1
          write(stdout,'(3x,a,/)') '- Outputs files: degeneracies.dat and wp_list.dat'
       end if !on_root

       ! Definition of the step along wp_box vectors and wp_box_k3 list
       wp_box_db1=1.0_dp/real(wp_box_kmesh(1),dp)
       wp_box_db2=1.0_dp/real(wp_box_kmesh(2),dp)
       wp_box_db3=(wp_box_k3max-wp_box_k3min)/wp_box_kmesh(3)

       do i3=0,wp_box_kmesh(3)
          wp_k3_list(i3)=wp_box_k3min+i3*wp_box_db3
       enddo
   
       ! Allocate all the arrays to the maximun number of walkers
       num_walkers=wp_box_kmesh(1)*wp_box_kmesh(2)*wp_box_kmesh(3)

	   allocate(wp_pt(3,num_walkers))
       allocate(E0_wp(num_walkers))
       allocate(gap_wp(num_walkers))
       allocate(deter_wp(num_walkers))
       allocate(charge_wp(num_walkers))
       allocate(modulus(num_walkers))

       modulus=0.0_dp
       charge_wp=0

       call get_HH_R

	   ! Vectors spanning the slice and corner (Cartesian)
       b1(:)=matmul(wp_box_b1(:),recip_lattice(:,:))
!       b1=b1/recip_lattice(1,1)
       b2(:)=matmul(wp_box_b2(:),recip_lattice(:,:))
!       b2=b2/recip_lattice(1,1)
       b3(:)=matmul(wp_box_b3(:),recip_lattice(:,:))
!	   b3=b3/recip_lattice(1,1)

       ! the box is restricted by kslice_k3max and kslice_k3min
       b3_box(:)=b3(:)*(kslice_k3max-kslice_k3min)

       ! z_vec (orthogonal to b1 and b2)
       zvec(1)=b1(2)*b2(3)-b1(3)*b2(2)
       zvec(2)=b1(3)*b2(1)-b1(1)*b2(3)
       zvec(3)=b1(1)*b2(2)-b1(2)*b2(1)
       ! y_vec (orthogonal to b1=x_vec)
       yvec(1)=zvec(2)*b1(3)-zvec(3)*b1(2)
       yvec(2)=zvec(3)*b1(1)-zvec(1)*b1(3)
       yvec(3)=zvec(1)*b1(2)-zvec(2)*b1(1)
       ! Moduli b_1,b_2,y_vec
       b1mod=sqrt(b1(1)**2+b1(2)**2+b1(3)**2)
       b2mod=sqrt(b2(1)**2+b2(2)**2+b2(3)**2)
       ymod=sqrt(yvec(1)**2+yvec(2)**2+yvec(3)**2)

       ! Cosine of the angle between b_1=x_vec and b_2
       cosb1b2=b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3)
       cosb1b2=cosb1b2/(b1mod*b2mod)      
       ! Cosine of the angle between y_vec and b_2
       cosyb2=yvec(1)*b2(1)+yvec(2)*b2(2)+yvec(3)*b2(3)
       cosyb2=cosyb2/(ymod*b2mod)
 
       wp_box_corner_cart=matmul(wp_box_corner,recip_lattice(:,:))
!       wp_box_corner_cart=wp_box_corner_cart/recip_lattice(1,1)

       wp_corner_cart=matmul(wp_corner,recip_lattice(:,:))
!	   wp_corner_cart=wp_corner_cart/recip_lattice(1,1)

       if (on_root) then
          file_name='degeneracies.dat'
	      file_unit=io_file_unit()
          open(file_unit,FILE=file_name,STATUS='UNKNOWN',&
                                                FORM='FORMATTED') 

          write(file_unit,'(4(A,x))') '  *         degeneracy point         *',&
                             '      GAP      ','WP energy ','   det(nu) '
          write(file_unit,'(4(A,x))') '  ------------------------------------',&
                             '   ---------  ','---------- ',' ----------'
       end if! on_root  

       counter_wp=0
!       write(stdout,*) "wp_box_b1=",wp_box_b1
!       write(stdout,*) "wp_box_b2=",wp_box_b2
!       write(stdout,*) "wp_box_b3=",wp_box_b3
!       write(stdout,*) "wp_box_b1=",wp_box_b1
       do i3=0,wp_box_kmesh(3)

          k3=wp_k3_list(i3)
            
          do i1=0,wp_box_kmesh(1)
             do i2=0,wp_box_kmesh(2)

                k1=i1*wp_box_db1
                k2=i2*wp_box_db2
                kpt=wp_box_corner+k1*wp_box_b1+k2*wp_box_b2+k3*wp_box_b3  
               write(stdout,*) i1,i2,i3," -  ",k1,k2,k3," -  ",kpt

                ! Obtain the kpt in Cartesian coordiantes
                ! In units of twopi/a [***NOTE*** Hard-wired for bcc structure!]
                !
                call utility_frac_to_cart(kpt,kpt_cart,recip_lattice)
!                kpt_cart=kpt_cart/recip_lattice(1,1)

                
                call get_kvec_wp(wp_num_filled_bands,kpt,kpt_min,gap=kpt_gap)         !minimiza

                ! Obtain the minimized point in Cartesian coordinates
                ! In units of twopi/a [***NOTE*** Hard-wired for bcc structure!]
                !
                call utility_frac_to_cart(kpt_min,kpt_cart,recip_lattice)
!                kpt_cart=kpt_cart/recip_lattice(1,1)

	            if (abs(kpt_gap) < wp_gap_thresh) then  

 		   ! Get the energy where the degeneracy is and the determiant of the
                   ! coefficients of the generalized Weyl Hamiltonian
                   call get_kdotp_2d(kpt_min,wp_num_filled_bands,J_vec,J_mat,det_out=deter)
                   energy_wp=J_vec(0)

                   ! Evaluate if the degeneracy point is inside the box
	               call in_or_out(kpt_cart,b1,b2,b3_box,wp_box_corner_cart,inbox)
 	               if(inbox) then 
                      
                      ! Plot in file all the degeneracies inside the box
                      if (on_root) then
                         if(kslice_coor) then     
	 		                !Fix this to write well points
                            ! Convert to (kpt_x,kpt_y), the 2D Cartesian coordinates
                            ! with x along x_vec=b1 and y along y_vec
                            !***************************************
                            kpt_cart=kpt_cart-wp_corner_cart

                            kpt_x=(kpt_cart(1)*b1(1)+kpt_cart(2)*b1(2)+kpt_cart(3)*b1(3))/b1mod
                            kpt_y=(kpt_cart(1)*yvec(1)+kpt_cart(2)*yvec(2)+kpt_cart(3)*yvec(3))/ymod
              
                            write(file_unit,'(6(f10.7,1x))')&
                            !in vector axis
                            kpt_x*recip_lattice(1,1),kpt_y*recip_lattice(1,1),kpt_cart(3),kpt_gap,energy_wp,deter
                         else
                            ! in cartesianas
                            write(file_unit,'(6(F12.7,x))')&
                            kpt_cart(1),kpt_cart(2),kpt_cart(3),kpt_gap,energy_wp,deter
                         end if !kslice_coor
                      end if !on_root

                      ! check chirality
	                  if(abs(deter) >= wp_det_thresh) then

                         if (counter_wp==0) then               ! The first wp is always included in the Weyl point list
                            counter_wp = 1
	                        wp_pt( : ,1)=kpt_cart 
                            charge_wp(1)=int(dsign(1.0_dp,deter))
                            E0_wp(1)=energy_wp
                            gap_wp(1)=kpt_gap
                            deter_wp(1)=deter
 	                     else          
                            ! Check if the new wp is already in the Weyl points list.
                            ! Check if the distance to each degeneracy in the list is bigger than wp_dis_thresh
                            ! if it is true then we add a new WP to the list.
	                        do n=1, counter_wp	           
                               vec_dum=wp_pt( : ,n)-kpt_cart
	                           modulus(n)=sqrt(DOT_PRODUCT(vec_dum,vec_dum)) !
                            end do !n

                            chivato=0
                            do n=1, counter_wp
                               if (modulus(n) > wp_dis_thresh) then
	                              chivato=chivato+1
 		                       end if
	                        end do
  
                            ! If the distace to all the WP in the list is bigger than wp_dis_thresh
                            ! i.e., the variable chivato has the same counts than counter_wp, then
                            ! we add the kpt_cart to the WP list as a new wp and increment the counter_wp 
                            ! in one unit. 

	                        if (chivato == counter_wp) then 
                               counter_wp=counter_wp+1    
	                           wp_pt( : ,counter_wp)=kpt_cart
                               charge_wp(counter_wp)=int(dsign(1.0_dp, deter))
	                           E0_wp(counter_wp)=energy_wp
                               gap_wp(counter_wp)=kpt_gap
                               deter_wp(counter_wp)=deter
                            end if ! chivato	            
	                     end if ! counter_wp

                      else 

                         ! In the case that the determinand of the Hamiltonian coefficients is smaller than
                         ! wp_det_thresh we can have two scenarios: i) the point is a multiple Weyl node
                         ! or ii) it is a degeneracy without chirality. We want to avoid this last group
                         ! since we want to provide a list of possible candidates to WP's. 

                         ! If they are multiple Weyl nodes we have to belong to a high symmetry rotation axis
                         ! with c4 or c6 symmetry, that should be provided in the input file. 
                         ! We check if they belong to shuch axis.

                         ! We check if kpt_cart belongs to c4 and c6 we call to check_c4_rot
                         ! To do: check_c4_rot
                                                 
                         call check_c4_rot(kpt_cart,multiple_wp)
                         if (multiple_wp) then
                            if(counter_wp==0) then
                               counter_wp = 1
	                           wp_pt( : ,1)=kpt_cart               
                               charge_wp(1)=int(dsign(2.0_dp, deter))
                               E0_wp(1)=energy_wp
                               gap_wp(1)=kpt_gap
				 			   deter_wp(1)=deter
 	                        else
	                           do n=1, counter_wp	           
                                  vec_dum=wp_pt( : ,n)-kpt_cart
	                              modulus(n)=sqrt(DOT_PRODUCT(vec_dum,vec_dum)) !
                               end do !n
                      
                               chivato=0
                      
                               do n=1, counter_wp
                                  if (modulus(n) > wp_dis_thresh) then
	                                 chivato=chivato+1
 		                          end if
	                           end do
  
	                           if (chivato == counter_wp) then

                                  counter_wp=counter_wp+1
	                              wp_pt( : ,counter_wp)=kpt_cart
                                  charge_wp(counter_wp)=int(dsign(2.0_dp, deter))
	                              E0_wp(counter_wp)=energy_wp
                                  gap_wp(counter_wp)=kpt_gap
                                  deter_wp(counter_wp)=deter

                                  if(on_root) then
                                     if(kslice_coor) then 
                                        write(file_unit,'(6(f10.7,1x))')&
                                        kpt_x,kpt_y,kpt_cart(3),kpt_gap,deter,energy_wp
                                     else
                                        write(file_unit,'(6(F12.7,x))')&
                                        kpt_cart(1),kpt_cart(2),kpt_cart(3),kpt_gap,energy_wp,deter
                                     end if
                                  end if !on_root

	                           end if ! chivato	            
	                        end if ! counter_wp
                         end if! multiple_wp
		              end if ! deter
	               end if ! inside
                end if ! kpt_gap
	         end do !i2
	      end do !i1
       end do !i3  
       ! to indicate the end of the file we add a line with all zeros
       if(on_root) then
          write(file_unit,'(6(f12.7,1x))')&
                     0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0
          close(file_unit)
       end if

       if(on_root) then

          file_name_wp='wp_list.dat'
          file_unit_wp=io_file_unit()
          open(file_unit_wp,FILE=file_name_wp,STATUS='UNKNOWN',&
                                                FORM='FORMATTED') 

          write(file_unit_wp,'(A,I4)') "Number of Weyl point candidates:",counter_wp
          write(file_unit_wp,*)        
          write(file_unit_wp,'(5(A,x))') 'Chirality', '*         degeneracy point         *',&
                             '      GAP      ','WP energy ','   det(nu) '
          write(file_unit_wp,'(5(A,x))') '---------','------------------------------------',&
                             '   ---------  ','---------- ',' ----------'

	      do n=1,counter_wp
             if(kslice_coor) then 
	            !Fix this to write well points
                ! Convert to (kpt_x,kpt_y), the 2D Cartesian coordinates
                ! with x along x_vec=b1 and y along y_vec
                !Since wp_pt contains the cartesian we have to write them in the appropiate coordinate system
                kpt_x=(wp_pt(1,n)*b1(1)+wp_pt(2,n)*b1(2)+wp_pt(3,n)*b1(3))/b1mod
                kpt_y=(wp_pt(1,n)*yvec(1)+wp_pt(2,n)*yvec(2)+wp_pt(3,n)*yvec(3))/ymod

                write(file_unit_wp,'(I3,1x,3(f12.7,1x),f4.1,x,f10.7)')n,&
                kpt_x*recip_lattice(1,1),kpt_y*recip_lattice(1,1),&
                wp_pt(3,n),charge_wp(n),E0_wp(n)

             else
                !write(file_unit,'(6(F12.7,x))')&
                !                        kpt_cart(1),kpt_cart(2),kpt_cart(3),wp_gap,energy_wp,deter
                write(file_unit_wp,'(I5,3x,6(f12.7,x))') charge_wp(n),&
                    (wp_pt(i,n),i=1,3),gap_wp(n),E0_wp(n),deter_wp(n)         
             end if !kslice_coor
       
          end do
          close(file_unit_wp)

       end if !on_root     
       deallocate(wp_pt)

  
  end subroutine berry_eval_wp_search

  subroutine cross_product(vec1,vec2,vec_out) 
    use w90_constants, only     : dp
    real(kind=dp), intent(in)   :: vec1(3),vec2(3)
    real(kind=dp), intent(out)  :: vec_out(3)
    vec_out=0.0_dp
    vec_out(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
    vec_out(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
    vec_out(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
  end subroutine cross_product


  subroutine in_or_out(pnt,t1,t2,t3,corner,eval)
  !===========================================================!
  !                                                           !
  ! eval=.true. if pnt is inside the Parallelepiped         !
  !                                    of t1,t2 and t3        !
  !                                                           !
  !===========================================================!

    use w90_constants, only     : dp
    
    use w90_parameters, only    : wp_box_delta


    real(kind=dp), intent(in) :: pnt(3),t1(3),t2(3),t3(3),corner(3)
    logical ,intent(out)     :: eval

    real(kind=dp)      :: n1(3),n2(3),n3(3),vdum1(3),vdum2(3)
    real(kind=dp)      :: t1_prime(3),t2_prime(3),t3_prime(3),corner_prime(3)

    ! debug variables
    real(kind=dp)      :: n1_prime(3),n2_prime(3),n3_prime(3)
    real(kind=dp)      :: rdum1,rdum2,rdum3
    real(kind=dp)      :: debug1(3),debug2(3),debug3(3)
    real(kind=dp)      :: delta

    ! definition of the normal vectors plane 1-2, 1-3 and 2-3

    delta=0.001_dp

    corner_prime(:)=corner(:)-wp_box_delta*0.5_dp*(t1(:)+t2(:)+t3(:))


    t1_prime(:)=(1.0d0+wp_box_delta)*t1(:)
    t2_prime(:)=(1.0d0+wp_box_delta)*t2(:)
    t3_prime(:)=(1.0d0+wp_box_delta)*t3(:)

    call cross_product(t1,t2,n1)
    call cross_product(t3,t1,n2) 
    call cross_product(t2,t3,n3)

    call cross_product(t1_prime,t2_prime,n1_prime)
    call cross_product(t3_prime,t1_prime,n2_prime) 
    call cross_product(t2_prime,t3_prime,n3_prime)

   ! print*, "+++++"
   ! print*, t1_prime
   ! print*, t2_prime
   ! print*, n1_prime
   ! print*, "+++++"

    call cross_product(n1,n1_prime,debug1)
    call cross_product(n2,n2_prime,debug2)
    call cross_product(n3,n3_prime,debug3)

   ! print*,
   ! print*, n1
   ! print*, n1_prime
   ! print*, debug1
   ! print*,
   ! print*, n2
   ! print*, n2_prime
   ! print*, debug2
   ! print*,
   ! print*, n3
   ! print*, n3_prime
   ! print*, debug3

    vdum1=pnt-corner_prime
    vdum2=pnt-(corner_prime+t3_prime)

    eval=.false.

    if (dot_product(n1_prime,vdum1) >= 0.0_dp .and. dot_product(n1_prime,vdum2)&
                                                         <= 0.0_dp) then
       vdum2=pnt-(corner_prime+t2_prime)
       if (dot_product(n2_prime,vdum1) >= 0.0_dp .and. dot_product(n2_prime,vdum2)&
                                                         <= 0.0_dp) then
          vdum2=pnt-(corner_prime+t1_prime)
          if (dot_product(n3_prime,vdum1) >= 0.0_dp .and. dot_product(n3_prime,vdum2)&
                                                         <= 0.0_dp) then
            eval=.true.
          end if	
       end if
    end if

  end subroutine in_or_out

  subroutine check_c4_rot(pnt,eval)
  !===========================================================!
  !                                                           
  !===========================================================!

    use w90_constants, only     : dp,pi,twopi,eps3
    
    use w90_parameters, only    : num_c4_axis,c4_axis,num_c6_axis,c6_axis


    real(kind=dp), intent(in) :: pnt(3)
    logical ,intent(out)     :: eval

    integer            :: i,j,counter

    real(kind=dp)      :: n(3),theta,R(3,3),rot_pnt(3)
    real(kind=dp)      :: rdum,vdum(3)

    eval=.false.
    do i=1,num_c4_axis
       vdum(:)=c4_axis(:,i)
       rdum=dsqrt(dot_product(vdum,vdum))
       n(:)=vdum(:)/rdum
       counter=0
       do j=1,3
          theta=real(j,dp)*(pi/2.0_dp)
          !
          ! Definition of the rotation matrix
          ! 
          R(1,1)=dcos(theta)+n(1)**2*(1.0_dp-dcos(theta))
          R(1,2)=n(1)*n(2)*(1.0_dp-dcos(theta))-n(3)*dsin(theta)
          R(1,3)=n(1)*n(3)*(1.0_dp-dcos(theta))+n(2)*dsin(theta)
          !
          R(2,1)=n(1)*n(2)*(1.0_dp-dcos(theta))+n(3)*dsin(theta)
          R(2,2)=dcos(theta)+n(2)**2*(1.0_dp-dcos(theta))
          R(2,3)=n(2)*n(3)*(1.0_dp-dcos(theta))-n(1)*dsin(theta)
          !
          R(3,1)=n(1)*n(3)*(1.0_dp-dcos(theta))-n(2)*dsin(theta)
          R(3,2)=n(2)*n(3)*(1.0_dp-dcos(theta))+n(1)*dsin(theta)
          R(3,3)=dcos(theta)+n(3)**2*(1.0_dp-dcos(theta))
          !
          rot_pnt=matmul(R,pnt)
          ! debug
          !print*, 
          !print*, "axis:",i," theta=",j*90
          !print*, "pnt=    ",pnt
          !print*, "rot_pnt=",rot_pnt
          vdum(:)=rot_pnt(:)-pnt(:)
          rdum=dsqrt(dot_product(vdum,vdum))
          !print*, "Distance between pnt and rot_pnt:",rdum
          if(rdum < eps3) then
             counter=counter+1
          end if
       end do !j
       if (counter == 3) then
          eval=.true.
       end if
       ! if a point belongs to an axis then we do not need to continue
       ! looking for another axis of rotation
       if (eval) then
          exit
       end if

    end do


  end subroutine check_c4_rot

  subroutine berry_get_total_spin(kpt,total_spin_k_list)
  !======================================================================!
  !                                                                      !
  ! Contribution from point k to the total spin of the system            !
  !                                                                      !
  !======================================================================!

    use w90_constants, only      : dp,cmplx_0
    use w90_utility, only        : utility_wgauss
    use w90_parameters, only     : num_wann,fermi_energy_list,&
                                   kubo_smr_fixed_en_width,kubo_adpt_smr,&
                                   kubo_adpt_smr_max,kubo_adpt_smr_fac,&
                                   kubo_smr_index,berry_kmesh,nfermi
    use w90_postw90_common, only : pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_eig_deleig,wham_get_D_h
    use w90_get_oper, only       : HH_R
    use w90_spin, only           : get_S
    use w90_io, only             : stdout

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), dimension(:,:), intent(out) :: total_spin_k_list

    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)

    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann,3),level_spacing,&
                        eta_smr,Delta_k,arg

    integer          :: n,ifermi
    real(kind=dp)    :: theta,eig(num_wann),S(num_wann,3)
    logical          :: got_spin,got_orb_n

    allocate(UU(num_wann,num_wann))
    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))

    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
    if(kubo_adpt_smr) Delta_k=pw90common_kmesh_spacing(berry_kmesh)

    total_spin_k_list=0.0_dp

    call get_S(kpt,S)
    
    do n=1,num_wann
       write(stdout,*) "spin of state ",n," : ",S(n,:)
       if(kubo_adpt_smr) then
          ! Eq.(35) YWVS07 
          level_spacing=sqrt(dot_product(del_eig(n,:),del_eig(n,:)))*Delta_k
          eta_smr=min(level_spacing*kubo_adpt_smr_fac,kubo_adpt_smr_max)
       else
          eta_smr=kubo_smr_fixed_en_width
       endif
       do ifermi=1,nfermi
          arg=(fermi_energy_list(ifermi)-eig(n))/eta_smr
          theta=utility_wgauss(arg,kubo_smr_index) ! Broadened theta( E_f-E_nk)
          total_spin_k_list(:,ifermi)=total_spin_k_list(:,ifermi)+S(n,:)*theta          
       enddo !ifermi
    enddo !n

       
  end subroutine berry_get_total_spin


  subroutine berry_get_morb_nk(kpt,morb_nk)
  !======================================================================!
  !                                                                      !
  ! calculation of the projection of band-resolved orbital moment        !
  !                                                                      !
  ! (2.hbar/e).m^orb_{kn,j}                                              !
  !                                                                      !
  !======================================================================!

    use w90_constants, only      : dp,cmplx_0

    use w90_parameters, only     : num_wann,fermi_energy_list,&
				spin_axis_polar, spin_axis_azimuth

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), dimension(:), intent(out), optional :: morb_nk


    ! Adaptive smearing
    !

    integer          :: i,n
    real(kind=dp)    :: morb_k(num_wann,3),conv,&
                        orb_nk(3),alpha(3)

    ! Unit vector along the magnetization direction
    !
    conv=180.0_dp/pi
    alpha(1)=sin(spin_axis_polar/conv)*cos(spin_axis_azimuth/conv)
    alpha(2)=sin(spin_axis_polar/conv)*sin(spin_axis_azimuth/conv)
    alpha(3)=cos(spin_axis_polar/conv)


    call berry_get_morb_k(kpt,morb_k)
    do n=1,num_wann
    	morb_nk(n)=alpha(1)*morb_k(n,1)+alpha(2)*morb_k(n,2)+alpha(3)*morb_k(n,3)
    enddo !n
       
  end subroutine berry_get_morb_nk


  subroutine berry_get_morb_k(kpt,morb_k,curv_k)
  !======================================================================!
  !                                                                      !
  ! calculation of the band-resolved orbital moment                      !
  !                                                                      !
  ! (2.hbar/e).m^orb_{kn,j}                                              !
  !  and optionall berry-curvature                                       !
  !======================================================================!

    use w90_constants, only      : dp,cmplx_0,eV_au,bohr

    use w90_parameters, only     : num_wann,fermi_energy_list,&
				spin_axis_polar, spin_axis_azimuth

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), dimension(:,:), intent(out)   :: morb_k
    real(kind=dp), dimension(:,:), intent(out), optional :: curv_k


    ! Adaptive smearing
    !

    integer          :: i,n
    real(kind=dp)    :: occ(num_wann),&
                        orb_nk(3),alpha(3),fac,conv,&
                        imf_k(3,3,1),img_k(3,3,1),imh_k(3,3,1)

     fac=eV_au/bohr**2
    do n=1,num_wann
        occ=0.0_dp
        occ(n)=1.0_dp
        call berry_get_imfgh_klist(kpt,imf_k,img_k,imh_k,occ)
        do i=1,3
           morb_k(n,i)=sum(imh_k(:,i,1))-sum(img_k(:,i,1))
	   if (present(curv_k)) then
		curv_k(n,i)=sum(imf_k(:,i,1))
	   endif
        enddo
    enddo !n
    morb_k(:,:)=morb_k(:,:)*fac
       
  end subroutine berry_get_morb_k





  subroutine berry_get_curv_k(kpt,curv_k)
  !======================================================================!
  !                                                                      !
  ! calculation of the band-resolved berry curvature                     !
  !                                                                      !
  !======================================================================!

    use w90_constants, only      : dp,cmplx_0,eV_au,bohr

    use w90_parameters, only     : num_wann,fermi_energy_list,&
				spin_axis_polar, spin_axis_azimuth
    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), dimension(:,:), intent(out):: curv_k

    integer          :: i,n
    real(kind=dp)    :: occ(num_wann),imf_k(3,3,1)

    do n=1,num_wann
        occ=0.0_dp
        occ(n)=1.0_dp
        call berry_get_imf_klist(kpt,imf_k,occ)
        do i=1,3
	    curv_k(n,i)=sum(imf_k(:,i,1))
        enddo
    enddo !n
       
  end subroutine berry_get_curv_k


  subroutine berry_get_curv_w_k(kpt,curv_w_k)
  !======================================================================!
  !                                                                      !
  ! calculation of the band-resolved                                     !
  ! frequency-dependent   berry curvature                                !
  !                                                                      !
  ! tildeOmeg(w)=                                                        !
  !     -eps_{bcd}sum_m ( w_mn^2/(wmn^2-w^2)) *Im[A_{nm,c}A_{mn,d}       !
  !  and optionall berry-curvature                                       !
  !======================================================================!


    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_rotate,utility_w0gauss
    use w90_parameters, only     : num_wann,kubo_nfreq,kubo_freq_list,&
				   kubo_smr_fixed_en_width,&
				  kubo_adpt_smr,&
                                   berry_band_list,num_berry_bands
    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec,&
                                   pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : AA_R
    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), dimension(:,:,:), intent(out):: curv_w_k  ! (num_wann,kubo_n_freq,3)


    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    real(kind=dp), allocatable    :: multWre(:)
     
    real(kind=dp)    :: del_eig(num_wann,3)
    
    integer          :: i,n,m,n1,m1
    real(kind=dp)    :: eig(num_wann),wmn

    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(multWre(kubo_nfreq))
    allocate(AA(num_wann,num_wann,3))
    allocate(D_h(num_wann,num_wann,3))
  

    if(.not.kubo_adpt_smr .and. kubo_smr_fixed_en_width/=0.0_dp)&
         kubo_freq_list=real(kubo_freq_list,dp)+cmplx_i*kubo_smr_fixed_en_width


    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
    call wham_get_D_h(delHH,UU,eig,D_h)
    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
    do i=1,3
	AA(:,:,i)=utility_rotate(AA(:,:,i),UU,num_wann)
    enddo
    AA=AA+cmplx_i*D_h ! Eq.(25) WYSV06

    curv_w_k(:,:,:)=0_dp
    
    do n1=1,num_berry_bands
       n=berry_band_list(n1)
       do m1=1,num_berry_bands
          m=berry_band_list(m1)
          if(n==m) cycle
          wmn=eig(m)-eig(n)
          multWre(:)=real(wmn**2/(wmn**2-kubo_freq_list(:)**2))

          do i=1,3
            	curv_w_k(n,:,i)=curv_w_k(n,:,i)-&
            	    2_dp*imag(AA(n,m,alpha_A(i))*AA(m,n,beta_A(i)))*multWre
          enddo
          
       enddo !m
    enddo !n

       
  end subroutine berry_get_curv_w_k



end module w90_berry
