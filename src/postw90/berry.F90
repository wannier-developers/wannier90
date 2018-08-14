!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

! ---------------------------------------------------------------

module w90_berry
  !! This module computes various "Berry phase" related properties
  !!
  !! Key REFERENCES
  !!
  !! *  WYSV06 = PRB 74, 195118 (2006)  (anomalous Hall conductivity - AHC)
  !! *  YWVS07 = PRB 75, 195121 (2007)  (Kubo frequency-dependent conductivity)
  !! *  LVTS12 = PRB 85, 014435 (2012)  (orbital magnetization and AHC)
  !! *  CTVR06 = PRB 74, 024408 (2006)  (  "          "       )
  ! ---------------------------------------------------------------
  !
  ! * Undocumented, works for limited purposes only: 
  !                                 reading k-points and weights from file

  use w90_constants, only : dp

  implicit none

  private

  public :: berry_main,berry_get_imf_klist,berry_get_imfgh_klist

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


  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine berry_main
  !============================================================!
  !                                                            !
  !! Computes the following quantities:                        
  !!   (i) Anomalous Hall conductivity (from Berry curvature)   
  !!  (ii) Complex optical conductivity (Kubo-Greenwood) & JDOS 
  !! (iii) Orbital magnetization                               
  !                                                            !
  !============================================================!

    use w90_constants, only      : dp,cmplx_0,elem_charge_SI,hbar_SI,&
                                   eV_au,bohr
    use w90_comms, only          : on_root,num_nodes,my_node_id,comms_reduce
    use w90_io, only             : io_error,stdout,io_file_unit,seedname,&
                                   io_stopwatch
    use w90_postw90_common, only : nrpts,irvec,num_int_kpts_on_node,int_kpts,&
                                   weight
    use w90_parameters, only     : timing_level,iprint,num_wann,berry_kmesh,&
                                   berry_curv_adpt_kmesh,&
                                   berry_curv_adpt_kmesh_thresh,&
                                   wanint_kpoint_file,cell_volume,transl_inv,&
                                   berry_task,berry_curv_unit,spin_decomp,&
                                   kubo_nfreq,kubo_freq_list,nfermi,&
                                   fermi_energy_list,shc_freq_scan
    use w90_get_oper, only       : get_HH_R,get_AA_R,get_BB_R,get_CC_R,&
                                   get_SS_R,get_SHC_R

    real(kind=dp), allocatable    :: adkpt(:,:)

    ! AHC and orbital magnetization, calculated for a list of Fermi levels
    !
    ! First index labels J0,J1,J2 terms, second labels the Cartesian component 
    !
    real(kind=dp) :: imf_k_list(3,3,nfermi),imf_list(3,3,nfermi),imf_list2(3,3,nfermi)
    real(kind=dp) :: img_k_list(3,3,nfermi),img_list(3,3,nfermi)
    real(kind=dp) :: imh_k_list(3,3,nfermi),imh_list(3,3,nfermi)
    real(kind=dp) :: ahc_list(3,3,nfermi)
    real(kind=dp) :: LCtil_list(3,3,nfermi),ICtil_list(3,3,nfermi),&
                     Morb_list(3,3,nfermi) 
    real(kind=dp) :: imf_k_list_dummy(3,3,nfermi) ! adaptive refinement of AHC
    
    ! Complex optical conductivity, dividided into Hermitean and
    ! anti-Hermitean parts
    !
    complex(kind=dp), allocatable :: kubo_H_k(:,:,:)
    complex(kind=dp), allocatable :: kubo_H(:,:,:)
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

    ! Spin Hall conductivity
    real(kind=dp), allocatable   :: shc_fermi(:),shc_k_fermi(:)
    complex(kind=dp),allocatable :: shc_freq(:),shc_k_freq(:)
    ! for fermi energy scan, adaptive kmesh
    real(kind=dp), allocatable   :: shc_k_fermi_dummy(:)

    real(kind=dp)     :: kweight,kweight_adpt,kpt(3),kpt_ad(3),&
                         db1,db2,db3,fac,freq,rdum,vdum(3)
    integer           :: n,i,j,k,ikpt,if,ispn,ierr,loop_x,loop_y,loop_z,&
                         loop_xyz,loop_adpt,adpt_counter_list(nfermi),ifreq,&
                         file_unit
    character(len=24) :: file_name
    logical           :: eval_ahc,eval_morb,eval_kubo,not_scannable,eval_shc

    if(nfermi==0) call io_error(&
         'Must specify one or more Fermi levels when berry=true')

    if (timing_level>1.and.on_root) call io_stopwatch('berry: prelims',1)

    ! Mesh spacing in reduced coordinates
    !
    db1=1.0_dp/real(berry_kmesh(1),dp)
    db2=1.0_dp/real(berry_kmesh(2),dp)
    db3=1.0_dp/real(berry_kmesh(3),dp)
    
    eval_ahc=.false.
    eval_morb=.false.
    eval_kubo=.false.
    eval_shc=.false.
    if(index(berry_task,'ahc')>0) eval_ahc=.true.
    if(index(berry_task,'morb')>0) eval_morb=.true.
    if(index(berry_task,'kubo')>0) eval_kubo=.true.
    if(index(berry_task,'shc')>0) eval_shc=.true.

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
       imf_list2=0.0_dp
       img_list=0.0_dp
       imh_list=0.0_dp
    endif

    ! List here berry_tasks that assume nfermi=1
    !
    not_scannable=eval_kubo
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

    if(eval_shc) then
        call get_HH_R
        call get_AA_R
        call get_SS_R
        call get_SHC_R

        if (shc_freq_scan) then
            allocate(shc_freq(kubo_nfreq))
            allocate(shc_k_freq(kubo_nfreq))
            shc_freq=0.0_dp
            shc_k_freq=0.0_dp
        else
            allocate(shc_fermi(nfermi))
            allocate(shc_k_fermi(nfermi))
            allocate(shc_k_fermi_dummy(nfermi))
            shc_fermi=0.0_dp
            shc_k_fermi=0.0_dp
            !only used for fermiscan
            shc_k_fermi_dummy=0.0_dp
            adpt_counter_list=0
        end if

    endif


    if(on_root) then

       write(stdout,'(/,/,1x,a)')&
            'Properties calculated in module  b e r r y'
       write(stdout,'(1x,a)')&
            '------------------------------------------'

       if(eval_ahc) write(stdout,'(/,3x,a)')&
            '* Anomalous Hall conductivity'

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

       if(eval_shc) then
           write(stdout,'(/,3x,a)') '* Spin Hall Conductivity'
           if (shc_freq_scan) then
               write(stdout,'(/,3x,a)') '  frequency scan'
           else
               write(stdout,'(/,3x,a)') '  fermi energy scan'
           end if
       endif

       if(transl_inv) then
          if(eval_morb)&
            call io_error('transl_inv=T disabled for morb')
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

          if(eval_morb) then
             call berry_get_imfgh_klist(kpt,imf_k_list,img_k_list,imh_k_list)
             imf_list2=imf_list2+imf_k_list*kweight
             img_list=img_list+img_k_list*kweight
             imh_list=imh_list+imh_k_List*kweight
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

          if(eval_shc) then
              ! be aware that index starts from 1

              if (.not. shc_freq_scan) then
                  call berry_get_shc_k(kpt,shc_k_fermi=shc_k_fermi)
                  do if=1,nfermi
                      rdum=abs(shc_k_fermi(if))
                      if(rdum>berry_curv_adpt_kmesh_thresh) then
                          adpt_counter_list(if)=adpt_counter_list(if)+1
                          do loop_adpt=1,berry_curv_adpt_kmesh**3
                              ! Using shc_k here would corrupt values for other
                              ! kpt, hence dummy. Only if-th element is used
                              call berry_get_shc_k(kpt(:)+adkpt(:,loop_adpt),&
                                      shc_k_fermi_dummy)
                              shc_fermi(if)=shc_fermi(if)+shc_k_fermi_dummy(if)*kweight_adpt
                          end do
                      else
                          shc_fermi(if)=shc_fermi(if)+shc_k_fermi(if)*kweight
                      endif
                  enddo
              else ! freq_scan, no adaptive kmesh
                  call berry_get_shc_k(kpt,shc_k_freq=shc_k_freq)
                  shc_freq = shc_freq + kweight*shc_k_freq
              end if

              !write(stdout,'(a6,i4,a6,3(f9.6,1x))') 'node ',my_node_id,&
              !        ' kpt ',kpt(1:3)
          end if
          !
          ! ***END COPY OF CODE BLOCK 1***

       end do !loop_xyz

    else ! Do not read 'kpoint.dat'. Loop over a regular grid in the full BZ

       kweight=db1*db2*db3
       kweight_adpt=kweight/berry_curv_adpt_kmesh**3

       do loop_xyz=my_node_id,PRODUCT(berry_kmesh)-1,num_nodes
          loop_x= loop_xyz/(berry_kmesh(2)*berry_kmesh(3))
          loop_y=(loop_xyz-loop_x*(berry_kmesh(2)&
               *berry_kmesh(3)))/berry_kmesh(3)
          loop_z=loop_xyz-loop_x*(berry_kmesh(2)*berry_kmesh(3))&
                -loop_y*berry_kmesh(3)
          kpt(1)=loop_x*db1
          kpt(2)=loop_y*db2
          kpt(3)=loop_z*db3

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

          if(eval_morb) then
             call berry_get_imfgh_klist(kpt,imf_k_list,img_k_list,imh_k_list)
             imf_list2=imf_list2+imf_k_list*kweight
             img_list=img_list+img_k_list*kweight
             imh_list=imh_list+imh_k_List*kweight
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

          if(eval_shc) then
            ! be aware that index starts from 1 

             if (.not. shc_freq_scan) then
                 call berry_get_shc_k(kpt,shc_k_fermi=shc_k_fermi)
                 do if=1,nfermi
                     rdum=abs(shc_k_fermi(if))
                     !TODO if adpt_kmesh=1, will this do the same kpt calculation again? add check adpt_kmesh!=1?
                     if(rdum>berry_curv_adpt_kmesh_thresh) then
                         adpt_counter_list(if)=adpt_counter_list(if)+1
                         do loop_adpt=1,berry_curv_adpt_kmesh**3
                             ! Using shc_k here would corrupt values for other
                             ! kpt, hence dummy. Only if-th element is used
                             call berry_get_shc_k(kpt(:)+adkpt(:,loop_adpt),&
                                     shc_k_fermi_dummy)
                             shc_fermi(if)=shc_fermi(if)+shc_k_fermi_dummy(if)*kweight_adpt
                         end do
                     else
                        shc_fermi(if)=shc_fermi(if)+shc_k_fermi(if)*kweight
                     endif
                 enddo
             else ! freq_scan, no adaptive kmesh
                 call berry_get_shc_k(kpt,shc_k_freq=shc_k_freq)
                 shc_freq = shc_freq + kweight*shc_k_freq
             end if

             !write(stdout,'(a6,i4,a6,3(f9.6,1x))') 'node ',my_node_id,&
             !        ' kpt ',kpt(1:3)
          end if
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

    if(eval_morb) then
       call comms_reduce(imf_list2(1,1,1),3*3*nfermi,'SUM')
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

    if(eval_shc) then
       if (shc_freq_scan) then
           call comms_reduce(shc_freq(1),kubo_nfreq,'SUM')
       else
           call comms_reduce(shc_fermi(1),nfermi,'SUM')
       end if
    end if
    
    if(on_root) then

       if (timing_level>1) call io_stopwatch('berry: k-interpolation',2)
       write(stdout,'(1x,a)') ' '
       if((eval_ahc) .and. berry_curv_adpt_kmesh.ne.1) then
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
       else if (eval_shc) then
           if(berry_curv_adpt_kmesh.ne.1) then
               if(.not.wanint_kpoint_file ) write(stdout, '(1x,a28,3(i0,1x))')&
                       'Regular interpolation grid: ',berry_kmesh
               if(.not. shc_freq_scan) then
                   write(stdout, '(1x,a28,3(i0,1x))') &
                           'Adaptive refinement grid: ',&
                            berry_curv_adpt_kmesh,berry_curv_adpt_kmesh,berry_curv_adpt_kmesh
                   write(stdout, '(1x,a28,a17,f6.2,a)')&
                       'Refinement threshold: ','SHC >',&
                       berry_curv_adpt_kmesh_thresh,' S/cm'
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
               endif
           else
               if(.not.wanint_kpoint_file) write(stdout, '(1x,a20,3(i0,1x))')&
                       'Interpolation grid: ',berry_kmesh(1:3)
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
                  -fermi_energy_list(if)*imf_list2(:,:,if))*fac
             ICtil_list(:,:,if)=(imh_list(:,:,if)&
                  -fermi_energy_list(if)*imf_list2(:,:,if))*fac
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

       ! -----------------------!
       ! Spin Hall conductivity !
       ! -----------------------!
       !
       if(eval_shc) then
           !
           ! Convert to S/cm
           !fac=1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
           !fac=elem_charge_SI*hbar_SI/cell_volume
           fac=1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
           if (shc_freq_scan) then
               shc_freq=shc_freq*fac
           else
               shc_fermi=shc_fermi*fac
           end if
           !
           write(stdout,'(/,1x,a)')&
                   '----------------------------------------------------------'
           write(stdout,'(1x,a)')&
                   'Output data files related to Spin Hall conductivity:'
           write(stdout,'(1x,a)')&
                   '----------------------------------------------------------'
           !
           if (.not. shc_freq_scan) then
               file_name= trim(seedname)//'-shc-fermiscan'//'.dat'
           else
               file_name= trim(seedname)//'-shc-freqscan'//'.dat'
           end if
           file_name=trim(file_name)
           file_unit=io_file_unit()
           write(stdout,'(/,3x,a)') '* '//file_name
           open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
           if (.not. shc_freq_scan) then
               write(file_unit,'(a,3x,a,3x,a)') 'No.','Fermi energy(eV)','SHC(S/cm)'
               do n=1,nfermi
                   write(file_unit,'(I4,1x,F12.6,1x,E17.8)') &
                       n,fermi_energy_list(n),shc_fermi(n)
               enddo
           else
               write(file_unit,'(a,3x,a,3x,a,3x,a)') 'No.','Frequency(eV)','Re(sigma)(S/cm)','Im(sigma)(S/cm)'
               do n=1,kubo_nfreq
                   write(file_unit,'(I4,1x,F12.6,1x,1x,2(E17.8,1x))') &
                           n,real(kubo_freq_list(n),dp),real(shc_freq(n),dp),aimag(shc_freq(n))
               end do
           end if
           close(file_unit)

       endif
       
    end if !on_root

    !write(stdout,'(a,i4)') 'berry_main end, my_node_id',my_node_id
    !if (eval_shc .and. wanint_kpoint_file) then
    !    write(file_name,'(a,a,I5.5,a)') trim(seedname),'-omega_xy_sz',my_node_id,'.dat'
    !    file_name=trim(file_name)
    !    file_unit=io_file_unit()
    !    write(stdout,'(/,3x,a)') '* '//file_name
    !    open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
    !    write(file_unit,'(a,3x,a,3x,a)') 'kpt(x,y,z)','nbnd','Omega'
    !    do loop_xyz=1,num_int_kpts_on_node(my_node_id)
    !        kpt(:)=int_kpts(:,loop_xyz)
    !        do n=1,num_wann
    !            write(file_unit,'(3E8.5,1I4,1E16.8)') &
    !                kpt(1),kpt(2),kpt(3),shc_k_list(n)
    !    enddo
    !    close(file_unit)
    !end if

  end subroutine berry_main


  subroutine berry_get_imf_klist(kpt,imf_k_list,occ)
  !============================================================!
  !                                                            !
  !! Calculates the Berry curvature traced over the occupied
  !! states, -2Im[f(k)] [Eq.33 CTVR06, Eq.6 LVTS12] for a list 
  !! of Fermi energies, and stores it in axial-vector form 
  !                                                            !
  !============================================================!
    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), intent(out), dimension(:,:,:) :: imf_k_list
    real(kind=dp), intent(in), optional, dimension(:) :: occ

    if (present(occ)) then 
        call berry_get_imfgh_klist(kpt,imf_k_list,occ=occ)
    else
        call berry_get_imfgh_klist(kpt,imf_k_list)
    endif

  end subroutine berry_get_imf_klist


  subroutine berry_get_imfgh_klist(kpt,imf_k_list,img_k_list,imh_k_list,occ)
  !=========================================================!
  !
  !! Calculates the three quantities needed for the orbital
  !! magnetization: 
  !!  
  !! * -2Im[f(k)] [Eq.33 CTVR06, Eq.6 LVTS12] 
  !! * -2Im[g(k)] [Eq.34 CTVR06, Eq.7 LVTS12] 
  !! * -2Im[h(k)] [Eq.35 CTVR06, Eq.8 LVTS12] 
  !! They are calculated together (to reduce the number of  
  !! Fourier calls) for a list of Fermi energies, and stored
  !! in axial-vector form. 
  !
  ! The two optional output parameters 'imh_k_list' and
  ! 'img_k_list' are only calculated if both of them are
  ! present.
  !
  !=========================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i
    use w90_utility, only        : utility_re_tr_prod,utility_im_tr_prod
    use w90_parameters, only     : num_wann,nfermi
    use w90_postw90_common, only : pw90common_fourier_R_to_k_vec,pw90common_fourier_R_to_k
    use w90_wan_ham, only        : wham_get_eig_UU_HH_JJlist,wham_get_occ_mat_list
    use w90_get_oper, only       : AA_R,BB_R,CC_R
    use w90_utility, only        : utility_zgemm_new

    ! Arguments
    !
    real(kind=dp), intent(in)     :: kpt(3)
    real(kind=dp), intent(out), dimension(:,:,:), optional & 
                                  :: imf_k_list, img_k_list, imh_k_list
    real(kind=dp), intent(in), optional, dimension(:) :: occ 

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: f_list(:,:,:)
    complex(kind=dp), allocatable :: g_list(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)
    complex(kind=dp), allocatable :: BB(:,:,:)
    complex(kind=dp), allocatable :: CC(:,:,:,:)
    complex(kind=dp), allocatable :: OOmega(:,:,:)
    complex(kind=dp), allocatable :: JJp_list(:,:,:,:)
    complex(kind=dp), allocatable :: JJm_list(:,:,:,:)
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i,j,ife,nfermi_loc
    real(kind=dp)                 :: s

    ! Temporary space for matrix products
    complex(kind=dp), allocatable, dimension(:,:,:) :: tmp

     if(present(occ)) then
        nfermi_loc=1
     else
        nfermi_loc=nfermi
     endif


    allocate(HH(num_wann,num_wann))
    allocate(UU(num_wann,num_wann))
    allocate(f_list(num_wann,num_wann,nfermi_loc))
    allocate(g_list(num_wann,num_wann,nfermi_loc))
    allocate(JJp_list(num_wann,num_wann,nfermi_loc,3))
    allocate(JJm_list(num_wann,num_wann,nfermi_loc,3))
    allocate(AA(num_wann,num_wann,3))
    allocate(OOmega(num_wann,num_wann,3))

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

    if(present(imf_k_list)) then
      ! Trace formula for -2Im[f], Eq.(51) LVTS12
      !
      do ife=1,nfermi_loc
        do i=1,3
          !
          ! J0 term (Omega_bar term of WYSV06)
          imf_k_list(1,i,ife) = &
             utility_re_tr_prod(f_list(:,:,ife), OOmega(:,:,i))
          !
          ! J1 term (DA term of WYSV06)
          imf_k_list(2,i,ife) = -2.0_dp * &
          ( &
              utility_im_tr_prod(AA(:,:,alpha_A(i)), JJp_list(:,:,ife,beta_A(i))) &
            + utility_im_tr_prod(JJm_list(:,:,ife,alpha_A(i)), AA(:,:,beta_A(i))) &
          )
          !
          ! J2 term (DD of WYSV06)
          imf_k_list(3,i,ife) = -2.0_dp * &
            utility_im_tr_prod(JJm_list(:,:,ife,alpha_A(i)),JJp_list(:,:,ife,beta_A(i)))
        end do
      end do
    end if

    if(present(img_k_list)) img_k_list = 0.0_dp
    if(present(imh_k_list)) imh_k_list = 0.0_dp

    if(present(img_k_list) .and. present(imh_k_list)) then
      allocate(BB(num_wann,num_wann,3))
      allocate(CC(num_wann,num_wann,3,3))

      allocate(tmp(num_wann,num_wann,5))
      ! tmp(:,:,1:3) ... not dependent on inner loop variables
      ! tmp(:,:,1) ..... HH . AA(:,:,alpha_A(i))
      ! tmp(:,:,2) ..... LLambda_ij [Eq. (37) LVTS12] expressed as a pseudovector
      ! tmp(:,:,3) ..... HH . OOmega(:,:,i)
      ! tmp(:,:,4:5) ... working matrices for matrix products of inner loop

      call pw90common_fourier_R_to_k_vec(kpt,BB_R,OO_true=BB)
      do j=1,3
        do i=1,j
          call pw90common_fourier_R_to_k(kpt,CC_R(:,:,:,i,j),CC(:,:,i,j),0)
          CC(:,:,j,i)=conjg(transpose(CC(:,:,i,j)))
        end do
      end do

      ! Trace formula for -2Im[g], Eq.(66) LVTS12
      ! Trace formula for -2Im[h], Eq.(56) LVTS12
      !
      do i=1,3
        call utility_zgemm_new(HH, AA(:,:,alpha_A(i)), tmp(:,:,1))
        call utility_zgemm_new(HH, OOmega(:,:,i), tmp(:,:,3))
        !
        ! LLambda_ij [Eq. (37) LVTS12] expressed as a pseudovector
        tmp(:,:,2) = cmplx_i*(CC(:,:,alpha_A(i),beta_A(i))&
                     -conjg(transpose(CC(:,:,alpha_A(i),beta_A(i)))))

        do ife=1,nfermi_loc
          !
          ! J0 terms for -2Im[g] and -2Im[h]
          !
          ! tmp(:,:,5) = HH . AA(:,:,alpha_A(i)) . f_list(:,:,ife) . AA(:,:,beta_A(i))
          call utility_zgemm_new(tmp(:,:,1), f_list(:,:,ife),   tmp(:,:,4))
          call utility_zgemm_new(tmp(:,:,4), AA(:,:,beta_A(i)), tmp(:,:,5))

          s = 2.0_dp * utility_im_tr_prod(f_list(:,:,ife), tmp(:,:,5));
          img_k_list(1,i,ife) = utility_re_tr_prod(f_list(:,:,ife), tmp(:,:,2)) - s
          imh_k_list(1,i,ife) = utility_re_tr_prod(f_list(:,:,ife), tmp(:,:,3)) + s

          !
          ! J1 terms for -2Im[g] and -2Im[h]
          !
          ! tmp(:,:,1) = HH . AA(:,:,alpha_A(i))
          ! tmp(:,:,4) = HH . JJm_list(:,:,ife,alpha_A(i))
          call utility_zgemm_new(HH, JJm_list(:,:,ife,alpha_A(i)), tmp(:,:,4))

          img_k_list(2,i,ife) = -2.0_dp * &
          ( &
              utility_im_tr_prod(JJm_list(:,:,ife,alpha_A(i)),BB(:,:,beta_A(i))) &
            - utility_im_tr_prod(JJm_list(:,:,ife,beta_A(i)),BB(:,:,alpha_A(i))) &
          )
          imh_k_list(2,i,ife) = -2.0_dp * &
          ( &
              utility_im_tr_prod(tmp(:,:,1), JJp_list(:,:,ife,beta_A(i))) &
            + utility_im_tr_prod(tmp(:,:,4), AA(:,:,beta_A(i))) &
          )

          !
          ! J2 terms for -2Im[g] and -2Im[h]
          !
          ! tmp(:,:,4) = JJm_list(:,:,ife,alpha_A(i)) . HH
          ! tmp(:,:,5) = HH . JJm_list(:,:,ife,alpha_A(i))
          call utility_zgemm_new(JJm_list(:,:,ife,alpha_A(i)), HH, tmp(:,:,4))
          call utility_zgemm_new(HH, JJm_list(:,:,ife,alpha_A(i)), tmp(:,:,5))

          img_k_list(3,i,ife) = -2.0_dp * &
             utility_im_tr_prod(tmp(:,:,4), JJp_list(:,:,ife,beta_A(i)))
          imh_k_list(3,i,ife) = -2.0_dp * &
             utility_im_tr_prod(tmp(:,:,5), JJp_list(:,:,ife,beta_A(i)))
        end do
      end do
      deallocate(tmp)
    end if

  end subroutine berry_get_imfgh_klist




  !===========================================================!
  !                   PRIVATE PROCEDURES                      ! 
  !===========================================================!

  subroutine berry_get_kubo_k(kpt,kubo_H_k,kubo_AH_k,jdos_k,&
                        kubo_H_k_spn,kubo_AH_k_spn,jdos_k_spn)
  !====================================================================!
  !                                                                    !
  !! Contribution from point k to the complex interband optical 
  !! conductivity, separated into Hermitian (H) and anti-Hermitian (AH)
  !! parts. Also returns the joint density of states 
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
             vdum(:)=del_eig(m,:)-del_eig(n,:)
             joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
             eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,&
                  kubo_adpt_smr_max)
          else
             eta_smr=kubo_smr_fixed_en_width
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




  !===========================================================!
  !                   PRIVATE PROCEDURES                      !
  !===========================================================!

  subroutine berry_get_shc_k(kpt,shc_k_fermi,shc_k_freq)
  !====================================================================!
  !                                                                    !
  !! Contribution from point k to the Spin Hall conductivity
  !! \sigma_{ipol,jpol}^{spol}(k), ipol,jpol,spol=1,2,3 (x,y,z,respectively)
  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi,hbar_SI
    use w90_utility, only        : utility_diagonalize,utility_rotate,utility_w0gauss
    use w90_parameters, only     : num_wann,kubo_eigval_max,&
                                   kubo_nfreq,kubo_freq_list,&
                                   kubo_adpt_smr,kubo_smr_fixed_en_width,&
                                   kubo_adpt_smr_max,kubo_adpt_smr_fac,&
                                   kubo_smr_index,berry_kmesh,&
                                   fermi_energy_list,nfermi,shc_freq_scan,&
                                   ipol=>shc_ipol,jpol=>shc_jpol,spol=>shc_spol
    use w90_postw90_common, only : pw90common_get_occ,pw90common_fourier_R_to_k_new,&
                                   pw90common_fourier_R_to_k_vec,pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : HH_R,AA_R
    use w90_comms,only           : my_node_id
    !!!
    use w90_io, only             : io_error,io_file_unit

    ! args
    real(kind=dp),              intent(in)  :: kpt(3)
    real(kind=dp), allocatable, optional, intent(out)    :: shc_k_fermi(:)
    complex(kind=dp), allocatable, optional, intent(out) :: shc_k_freq(:)

    ! internal vars
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    complex(kind=dp), allocatable :: D_h(:,:,:)
    complex(kind=dp), allocatable :: AA(:,:,:)

    complex(kind=dp)              :: Jsi_k(num_wann,num_wann)

    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann,3),joint_level_spacing,&
                        eta_smr,Delta_k,arg,vdum(3)

    integer          :: n,m,ispn,i,ifreq
    real(kind=dp)    :: eig(num_wann)
    real(kind=dp),allocatable    :: occ_list(:,:)
    complex(kind=dp),allocatable :: omega_list(:)
    real(kind=dp)    :: omega,rfac
    complex(kind=dp) :: prod,cdum,cfac

    write(*,'((a),1x,1I4,1x,3(f9.6,1x))') &
            'node',my_node_id,kpt(1:3)

    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))
    allocate(D_h(num_wann,num_wann,3))
    allocate(AA(num_wann,num_wann,3))
    if (shc_freq_scan) then
        allocate(shc_k_freq(kubo_nfreq))
        allocate(occ_list(num_wann,1))
        allocate(omega_list(kubo_nfreq))
        shc_k_freq=0.0_dp
    else
        allocate(shc_k_fermi(nfermi))
        allocate(occ_list(num_wann,nfermi))
        shc_k_fermi=0.0_dp
    end if


    call wham_get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
    call wham_get_D_h(delHH,UU,eig,D_h)

    call pw90common_fourier_R_to_k_vec(kpt,AA_R,OO_true=AA)
    do i=1,3
       AA(:,:,i)=utility_rotate(AA(:,:,i),UU,num_wann)
    enddo
    AA=AA+cmplx_i*D_h ! Eq.(25) WYSV06

    call berry_get_js_k(kpt,&
            eig,del_eig(:,ipol),D_h(:,:,ipol),UU,&
            Jsi_k)

    if (shc_freq_scan) then
        if(kubo_adpt_smr) then
           Delta_k=pw90common_kmesh_spacing(berry_kmesh)
        end if
        call pw90common_get_occ(eig,occ_list(:,1),fermi_energy_list(1))
    else
      ! get occ for different fermi_energy
      do i=1,nfermi
         call pw90common_get_occ(eig,occ_list(:,i),fermi_energy_list(i))
      end do
    end if

    !write(*,'((a),3(f9.6,1x))') 'kpt', kpt(1),kpt(2),kpt(3)
    do n=1,num_wann
       ! get \Omega_{n,ipol jpol}^{spol}
       omega=0.0_dp
       do m=1,num_wann
          if(m==n) cycle
          if(eig(m)>kubo_eigval_max .or. eig(n)>kubo_eigval_max) cycle

          rfac = eig(m)-eig(n)
          !this will calculate AHC
          !prod = -rfac*cmplx_i*AA(n,m,ipol) * rfac*cmplx_i*AA(m,n,jpol)
          prod = Jsi_k(n,m) * cmplx_i*rfac*AA(m,n,jpol)
          if (kubo_adpt_smr) then
              ! Eq.(35) YWVS07
              vdum(:)=del_eig(m,:)-del_eig(n,:)
              joint_level_spacing=sqrt(dot_product(vdum(:),vdum(:)))*Delta_k
              eta_smr=min(joint_level_spacing*kubo_adpt_smr_fac,&
                      kubo_adpt_smr_max)
          else
              eta_smr=kubo_smr_fixed_en_width
          endif
          if (shc_freq_scan) then
              do ifreq=1,kubo_nfreq
                  cdum = real(kubo_freq_list(ifreq),dp)+cmplx_i*eta_smr
                  cfac = (occ_list(m,1)-occ_list(n,1))/(rfac**2-cdum**2)
                  omega_list(ifreq) = omega_list(ifreq) + cfac*aimag(prod)
              end do
          else
            rfac = -2.0_dp/(rfac**2+eta_smr**2)
            omega = omega + rfac*aimag(prod)
          end if
       enddo

       !write(*,'(3(f9.6,1x),1I4,1E16.8)') &
       !        kpt(1),kpt(2),kpt(3),n,omega

       if (.not. shc_freq_scan) then
         do i=1,nfermi
            shc_k_fermi(i) = shc_k_fermi(i) + occ_list(n,i)*omega
         end do
       else
         shc_k_freq = shc_k_freq + omega_list
       end if
    enddo

    !do i = 1,nfermi
    !  write(*,'(3(f9.6,1x),1I4,1E16.8)') &
    !        kpt(1),kpt(2),kpt(3),i,shc_k(i)
    !end do

    return

  contains

  !===========================================================!
  !                   PRIVATE PROCEDURES                      !
  !===========================================================!

  subroutine berry_get_js_k(kpt,eig,del_i_eig,D_i_h,UU,js_k)
  !====================================================================!
  !                                                                    !
  !! Contribution from point k to the
  !!    < \psi_k | 1/2 * (s_z v_x + v_x s_z) | \psi_k >
  !
  !  Junfeng Qiao (8/7/2018)
  !                                                                    !
  !====================================================================!

    use w90_constants, only      : dp,cmplx_0,cmplx_i,pi
    use w90_utility, only        : utility_diagonalize,utility_rotate,utility_w0gauss
    use w90_parameters, only     : num_wann,kubo_eigval_max,&
                                   ipol=>shc_ipol,spol=>shc_spol
    use w90_postw90_common, only : pw90common_get_occ,pw90common_fourier_R_to_k_new,&
                                   pw90common_fourier_R_to_k_vec,pw90common_kmesh_spacing
    use w90_wan_ham, only        : wham_get_D_h,wham_get_eig_deleig
    use w90_get_oper, only       : HH_R,AA_R,SS_R,SR_R,SHR_R,SH_R

    ! args
    real(kind=dp),                                  intent(in)  :: kpt(3)
    real(kind=dp),              dimension(:),       intent(in)  :: eig
    real(kind=dp),              dimension(:),       intent(in)  :: del_i_eig
    complex(kind=dp),           dimension(:,:),     intent(in)  :: D_i_h
    complex(kind=dp),           dimension(:,:),     intent(in)  :: UU
    complex(kind=dp),           dimension(:,:),     intent(out) :: js_k

    ! internal vars
    complex(kind=dp)    :: K_k(num_wann,num_wann)
    complex(kind=dp)    :: L_k(num_wann,num_wann)
    complex(kind=dp)    :: B_k(num_wann,num_wann)
    complex(kind=dp)    :: S_k(num_wann,num_wann)
    complex(kind=dp)    :: SR(num_wann,num_wann,3)
    complex(kind=dp)    :: SR_i(num_wann,num_wann)
    complex(kind=dp)    :: SHR(num_wann,num_wann,3)
    complex(kind=dp)    :: SHR_i(num_wann,num_wann)
    complex(kind=dp)    :: SH(num_wann,num_wann,3)
    complex(kind=dp)    :: eig_mat(num_wann,num_wann)
    complex(kind=dp)    :: del_eig_mat(num_wann,num_wann)
    complex(kind=dp)    :: S_w(num_wann,num_wann)

    !===========
    js_k = cmplx_0

    !=========== S_k ===========
    ! < u_k | s_spol | u_k >
    !call spin_get_SS_k(kpt,spol,S_k)
    call pw90common_fourier_R_to_k_new(kpt,SS_R(:,:,:,spol),OO=S_w)
    S_k(:,:)=utility_rotate(S_w,UU,num_wann)

    !=========== K_k ===========
    ! < u_k | s_spol | \partial_ipol u_k >
    call pw90common_fourier_R_to_k_vec(kpt,SR_R(:,:,:,spol,:),OO_true=SR)
    SR_i = -cmplx_i * utility_rotate(SR(:,:,ipol),UU,num_wann)
    K_k = SR_i + matmul(S_k,D_i_h) !


    !=========== L_k ===========
    ! < u_k | s_spol.H | \partial_ipol u_k >
    call pw90common_fourier_R_to_k_vec(kpt,SHR_R(:,:,:,spol,:),OO_true=SHR)
    SHR_i = -cmplx_i * utility_rotate(SHR(:,:,ipol),UU,num_wann)
    call pw90common_fourier_R_to_k_vec(kpt,SH_R,OO_true=SH)
    SH(:,:,spol) = utility_rotate(SH(:,:,spol),UU,num_wann)
    L_k = SHR_i + matmul(SH(:,:,spol),D_i_h) !


    !=========== B_k ===========
    ! < \psi_nk | s_z v_x | \psi_mk >
    B_k = cmplx_0
    do i=1,num_wann
        eig_mat(i,:) = eig(:)
        del_eig_mat(i,:) = del_i_eig(:)
    end do
    ! note * is not matmul
    B_k =  del_eig_mat*S_k + eig_mat*K_k - L_k


    !=========== Js_k ===========
    ! note the S in SR_R,SHR_R,SH_R of get_SHC_R is sigma,
    ! we need to multiply it by hbar/2, and this is
    ! cancelled by 1/hbar so only 1/2 is left
    js_k = 1.0_dp/4.0_dp * (B_k + conjg(transpose(B_k)))

  end subroutine berry_get_js_k

  end subroutine berry_get_shc_k

end module w90_berry
