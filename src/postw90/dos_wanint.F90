!-*- mode: F90; mode: font-lock -*-!

module w90_dos_wanint

  use w90_constants, only : dp

  implicit none

  private

  public :: dos,find_fermi_level

  integer       :: num_freq
  real(kind=dp) :: d_omega

contains

  !=========================================================!
  !                   PUBLIC PROCEDURES                     ! 
  !=========================================================!

  subroutine dos
    !=========================================================!
    !                                                         !
    ! Computes the electronic density of states using         !
    ! adaptive broadening: PRB 75, 195121 (2007) [YWVS07].    !
    ! Can resolve the DOS into up-spin and down-spin parts    !
    !                                                         !
    !=========================================================!

    use w90_io, only            : io_error,io_file_unit,io_date,io_stopwatch,&
         seedname,stdout
    use w90_comms, only         : on_root,num_nodes,my_node_id,comms_reduce
    use w90_wanint_common, only : num_int_kpts_on_node,int_kpts,weight
    use w90_parameters, only    : num_wann,dos_num_points,dos_min_energy,&
         dos_max_energy,dos_energy_step,timing_level,&
         wanint_kpoint_file,adpt_smr_steps,&
         adpt_smr_width,spn_decomp
    use w90_get_oper, only      : get_HH_R,get_SS_R


    ! 'dos_k' contains contrib. from one k-point, 'dos_node' from k-points
    ! in one node, 'dos_all' from all nodes/k-points
    !
    real(kind=dp), allocatable :: dos_k(:,:,:)
    real(kind=dp), allocatable :: dos_all(:,:,:)

    real(kind=dp)    :: kweight,kpt(3),omega
    integer          :: i,loop_x,loop_y,loop_z,loop_kpt,loop_f
    integer          :: dos_unit,ndim,ierr

    num_freq=nint((dos_max_energy-dos_min_energy)/dos_energy_step)
    if(num_freq==0) num_freq=1
    d_omega=(dos_max_energy-dos_min_energy)/(num_freq-1)

    call get_HH_R
    if(spn_decomp) then
       ndim=3
       call get_SS_R
    else
       ndim=1
    end if

    allocate(dos_k(num_freq,adpt_smr_steps,ndim))
    allocate(dos_all(num_freq,adpt_smr_steps,ndim))

    if(on_root) then

       if (timing_level>1) call io_stopwatch('dos_wanint: dos',1)

       write(stdout,'(/,1x,a)') '============'
       write(stdout,'(1x,a)')   'Calculating:'
       write(stdout,'(1x,a)')   '============'

       write(stdout,'(/,3x,a)') '* Density of states (_dos)'

       write(stdout,'(/,5x,a,f9.4,a,f9.4,a)')&
            'Energy range: [',dos_min_energy,',',dos_max_energy,'] eV'

       write(stdout,'(/,5x,a,6(f6.3,1x))')&
            'Adaptive smearing width prefactors: ',&
            (adpt_smr_width(i),i=1,adpt_smr_steps)

       if(dos_num_points < 10) then
          write(stdout,'(5x,a,i1,a,i1,a,i1,a,i12,a)')&
               'Interpolation mesh in full BZ: ',&
               dos_num_points,'x',dos_num_points,'x',&
               dos_num_points,'=',dos_num_points**3,' points'
       elseif(dos_num_points < 100) then
          write(stdout,'(5x,a,i2,a,i2,a,i2,a,i10,a)')&
               'Interpolation mesh in full BZ: ',&
               dos_num_points,'x',dos_num_points,'x',&
               dos_num_points,'=',dos_num_points**3,' points'
       elseif(dos_num_points < 1000) then
          write(stdout,'(5x,a,i3,a,i3,a,i3,a,i12,a)')&
               'Interpolation mesh in full BZ: ',&
               dos_num_points,'x',dos_num_points,'x',&
               dos_num_points,'=',dos_num_points**3,' points'
       else
          write(stdout,'(5x,a,i4,a,i4,a,i4,a,i12,a)')&
               'Interpolation mesh in full BZ: ',&
               dos_num_points,'x',dos_num_points,'x',&
               dos_num_points,'=',dos_num_points**3,' points'
       end if

    end if

    dos_all=0.0_dp

    if(wanint_kpoint_file) then
       !
       ! Unlike for optical properties, this should always work for DOS
       !
       if(on_root)  write(stdout,'(/,1x,a)') 'Sampling the irreducible BZ only'

       ! Loop over k-points on the irreducible wedge of the Brillouin zone,
       ! read from file 'kpoint.dat'
       !
       ! ---------------------------------------------------------------------
       ! NOTE: Still need to set the variable 'dos_num_points' in the .wanint
       !       file to the linear dimensions of the corresponding nominal 
       !       interpolation mesh in the full BZ. Ideally that information 
       !       should be contained in the file 'kpoint.dat', and if the 
       !       variable is also set in the .wanint file, the code should check 
       !       that they are equal
       ! ---------------------------------------------------------------------
       !
       do loop_kpt=1,num_int_kpts_on_node(my_node_id)
          kpt(:)=int_kpts(:,loop_kpt)
          call dos_kpt(kpt,dos_k)
          dos_all=dos_all+dos_k*weight(loop_kpt)
       end do

    else

       if (on_root) write(stdout,'(/,1x,a)') 'Sampling the full BZ'

       kweight=1.0_dp/dos_num_points**3
       do loop_kpt=my_node_id,dos_num_points**3-1,num_nodes
          loop_x=loop_kpt/dos_num_points**2
          loop_y=(loop_kpt-loop_x*dos_num_points**2)/dos_num_points
          loop_z=loop_kpt-loop_x*dos_num_points**2-loop_y*dos_num_points
          kpt(1)=real(loop_x,dp)/dos_num_points
          kpt(2)=real(loop_y,dp)/dos_num_points
          kpt(3)=real(loop_z,dp)/dos_num_points
          call dos_kpt(kpt,dos_k)
          dos_all=dos_all+dos_k*kweight
       end do

    end if

    ! Collect contributions from all nodes
    !
    call comms_reduce(dos_all(1,1,1),num_freq*adpt_smr_steps*ndim,'SUM')

    if(on_root) then
       write(stdout,'(/,/,1x,a)') '------------------'
       write(stdout,'(1x,a)')     'Output data files:'
       write(stdout,'(1x,a)')     '------------------'
       write(stdout,'(/,3x,a)') trim(seedname)//'_dos.dat'
       dos_unit=io_file_unit()
       open(dos_unit,FILE=trim(seedname)//'_dos.dat',STATUS='UNKNOWN',&
            FORM='FORMATTED')
       do loop_f=1,num_freq
          omega=dos_min_energy+(loop_f-1)*d_omega
          write(dos_unit,'(20E16.8)') omega,dos_all(loop_f,:,:)
       enddo
       close(dos_unit)
       if (timing_level>1) call io_stopwatch('dos_wanint: dos',2)
    end if

  end subroutine dos

  ! =========================================================================

  subroutine find_fermi_level
    !==============================================!
    !                                              !
    ! Finds the Fermi level by integrating the DOS !
    !                                              !
    !==============================================!

    use w90_io, only            : stdout,io_error
    use w90_comms
    use w90_wanint_common, only : max_int_kpts_on_node,num_int_kpts_on_node,&
         int_kpts,weight
    use w90_parameters, only    : fermi_energy,found_fermi_energy,&
         num_elec_cell,&
         num_wann,dos_num_points,dos_min_energy,&
         dos_max_energy,dos_energy_step,&
         wanint_kpoint_file,&
         adpt_smr_steps,adpt_smr_width

#ifdef MPI 
    include 'mpif.h'
#endif

    real(kind=dp) :: kpt(3),sum_max_node,sum_max_all,&
         sum_mid_node,sum_mid_all,emin,emax,emid,&
         emin_node(0:num_nodes-1),emax_node(0:num_nodes-1),&
         ef(adpt_smr_steps)
    integer       :: loop_x,loop_y,loop_z,loop_kpt,loop_s,loop_nodes,&
         loop_iter,ierr,num_int_kpts,ikp

    real(kind=dp), allocatable :: eig_node(:,:)
    real(kind=dp), allocatable :: levelspacing_node(:,:)

    if(on_root) write(stdout,'(/,a)') 'Finding the value of the Fermi level'

    if(.not.wanint_kpoint_file) then
       !
       ! Already done in wanint_get_kpoint_file if 
       ! wanint_kpoint_file=.true.
       !
       allocate(num_int_kpts_on_node(0:num_nodes-1))
       num_int_kpts=dos_num_points**3
       !
       ! Local k-point counter on each node (lazy way of doing it, there is
       ! probably a smarter way)
       !
       ikp=0
       do loop_kpt=my_node_id,num_int_kpts-1,num_nodes
          ikp=ikp+1
       end do
       num_int_kpts_on_node(my_node_id)=ikp
#ifdef MPI
       call MPI_reduce(ikp,max_int_kpts_on_node,1,MPI_integer,&
            MPI_MAX,0,MPI_COMM_WORLD,ierr)
#else
       max_int_kpts_on_node=ikp
#endif
       call comms_bcast(max_int_kpts_on_node,1)
    end if

    allocate(eig_node(num_wann,max_int_kpts_on_node),stat=ierr)
    if (ierr/=0)&
         call io_error('Error in allocating eig_node in find_fermi_level')
    eig_node=0.0_dp
    allocate(levelspacing_node(num_wann,max_int_kpts_on_node),stat=ierr)
    if (ierr/=0)&
         call io_error('Error in allocating levelspacing in find_fermi_level')
    levelspacing_node=0.0_dp

    if(wanint_kpoint_file) then
       if(on_root) write(stdout,'(/,1x,a)') 'Sampling the irreducible BZ only'
       do loop_kpt=1,num_int_kpts_on_node(my_node_id)
          kpt(:)=int_kpts(:,loop_kpt)
          call get_eig_levelspacing_k(kpt,eig_node(:,loop_kpt),&
               levelspacing_node(:,loop_kpt))
       end do
    else
       if (on_root)&
            write(stdout,'(/,1x,a)') 'Sampling the full BZ (not using symmetry)'
       allocate(weight(max_int_kpts_on_node),stat=ierr)
       if (ierr/=0)&
            call io_error('Error in allocating weight in find_fermi_level')
       weight=0.0_dp
       ikp=0
       do loop_kpt=my_node_id,num_int_kpts-1,num_nodes
          ikp=ikp+1
          loop_x=loop_kpt/dos_num_points**2
          loop_y=(loop_kpt-loop_x*dos_num_points**2)/dos_num_points
          loop_z=loop_kpt-loop_x*dos_num_points**2-loop_y*dos_num_points
          kpt(1)=real(loop_x,dp)/dos_num_points
          kpt(2)=real(loop_y,dp)/dos_num_points
          kpt(3)=real(loop_z,dp)/dos_num_points
          weight(ikp)=1.0_dp/dos_num_points**3
          call get_eig_levelspacing_k(kpt,eig_node(:,ikp),&
               levelspacing_node(:,ikp))
       end do
    end if

    do loop_s=1,adpt_smr_steps

       ! Find minimum and maximum band energies within projected subspace
       !
       emin_node(my_node_id)=&
            minval(eig_node(1,1:num_int_kpts_on_node(my_node_id)))
       emax_node(my_node_id)=&
            maxval(eig_node(num_wann,1:num_int_kpts_on_node(my_node_id)))
       if(.not.on_root) then
          call comms_send(emin_node(my_node_id),1,root_id)
          call comms_send(emax_node(my_node_id),1,root_id)
       else
          do loop_nodes=1,num_nodes-1
             call comms_recv(emin_node(loop_nodes),1,loop_nodes)
             call comms_recv(emax_node(loop_nodes),1,loop_nodes)
          end do
          emin=minval(emin_node)
          emax=maxval(emax_node)
       end if
       call comms_bcast(emin,1)
       call comms_bcast(emax,1)

       ! Check that the Fermi level lies within the projected subspace
       !
       sum_max_node=count_states(emax,eig_node,levelspacing_node,loop_s,&
            num_int_kpts_on_node(my_node_id))
#ifdef MPI
       call MPI_reduce(sum_max_node,sum_max_all,1,MPI_DOUBLE_PRECISION,&
            MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
       sum_max_all=sum_max_node
#endif
       if(on_root) then
          if(num_elec_cell>sum_max_all) then
             write(stdout,*) 'Something wrong in find_fermi_level:'
             write(stdout,*) '   Fermi level does not lie within projected subspace'
             write(stdout,*) 'num_elec_cell= ',num_elec_cell
             write(stdout,*) 'sum_max_all= ',sum_max_all
             stop 'Stopped: see output file'
          end if
       end if

       ! Now interval search for the Fermi level
       !
       do loop_iter=1,1000
          emid=(emin+emax)/2.0_dp
          sum_mid_node=count_states(emid,eig_node,levelspacing_node,loop_s,&
               num_int_kpts_on_node(my_node_id))
#ifdef MPI
          call MPI_reduce(sum_mid_node,sum_mid_all,1,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
          sum_mid_all=sum_mid_node
#endif
          ! This is needed because MPI_reduce only returns sum_mid_all to the 
          ! root (To understand: could we use MPI_Allreduce instead?)
          !
          call comms_bcast(sum_mid_all,1)
          if(abs(sum_mid_all-num_elec_cell) < 1.e-10_dp) then
             !
             ! NOTE: Here should assign a value to an entry in a fermi-level 
             !       vector. Then at the end average over adaptive smearing 
             !       widths and broadcast the result
             !
             ef(loop_s)=emid
             exit
          elseif((sum_mid_all-num_elec_cell) < -1.e-10_dp) then
             emin=emid
          else
             emax=emid
          end if
       end do

    end do

    fermi_energy=sum(ef)/real(adpt_smr_steps,dp)
    found_fermi_energy=.true.
    if(on_root) then
       if(adpt_smr_steps>1) then
          write(stdout,*) ' '
          do loop_s=1,adpt_smr_steps
             write(stdout,'(1x,a,f5.2,2x,a,f10.6)')&
                  'Adaptive smearing prefact= ',&
                  adpt_smr_width(loop_s),'Fermi energy= ',ef(loop_s)
          end do
          write(stdout,'(1x,a)')&
               '---------------------------------------------------------'
          write(stdout,'(39x,a,F10.6,1x)') 'Average= ',fermi_energy
       else
          write(stdout,'(/,1x,a,F10.6,1x,a)') 'Fermi energy= ',fermi_energy,'eV'
       end if
       write(stdout,'(/,1x,a,F10.6,1x,a)')&
            'In subsequent calculations Fermi energy=',fermi_energy,'eV'
    end if

  end subroutine find_fermi_level

  !=========================================================!
  !                   PRIVATE PROCEDURES                    ! 
  !=========================================================!

  subroutine dos_kpt(kpt,dos_k)
    !=========================================================!
    !                                                         !
    ! Calculates the contribution from one k-point to the DOS !
    !                                                         !
    !=========================================================!

    use w90_constants, only     : dp
    use w90_utility, only       : w0gauss
    use w90_parameters, only    : num_wann,dos_min_energy,dos_num_points,&
         adpt_smr_steps,adpt_smr_width,spn_decomp
    use w90_spin_wanint, only   : get_spn_nk

    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), dimension(:,:,:), intent(out) :: dos_k

    ! Adaptive smearing
    !
    real(kind=dp) :: eig_k(num_wann),levelspacing_k(num_wann),smear,arg

    ! Misc/Dummy
    !
    integer          :: i,loop_f,loop_s
    real(kind=dp)    :: rdum,omega,spn_nk(num_wann),alpha_sq,beta_sq 

    call get_eig_levelspacing_k(kpt,eig_k,levelspacing_k)

    ! Get spin projections for every band
    !
    if(spn_decomp) call get_spn_nk(kpt,spn_nk)

    dos_k=0.0_dp
    do i=1,num_wann
       do loop_s=1,adpt_smr_steps
          !
          ! Except for the factor 1/sqrt(2), this is Eq.(34) YWVS07
          ! !!!UNDERSTAND THAT FACTOR!!!
          !
          smear=levelspacing_k(i)*adpt_smr_width(loop_s)/sqrt(2.0_dp)
          do loop_f=1,num_freq
             omega=dos_min_energy+(loop_f-1)*d_omega
             arg=(omega-eig_k(i))/smear
             if(abs(arg) > 10.0_dp) then ! optimization
                cycle
             else
                !
                ! Adaptive broadening of the delta-function in Eq.(39) YWVS07
                !
                ! hard code for M-P (1)
                rdum=w0gauss(arg,1)/smear
             end if
             !
             ! Contribution to total DOS
             !
             dos_k(loop_f,loop_s,1)=dos_k(loop_f,loop_s,1)+rdum
             if(spn_decomp) then
                !
                ! Contribution to spin-up DOS of Bloch spinor with component 
                ! (alpha,beta) with respect to the chosen quantization axis
                !
                alpha_sq=(1.0_dp+spn_nk(i))/2.0_dp ! |alpha|^2
                dos_k(loop_f,loop_s,2)=dos_k(loop_f,loop_s,2)+rdum*alpha_sq
                !
                ! Contribution to spin-down DOS 
                !
                beta_sq=1.0_dp-alpha_sq ! |beta|^2 = 1 - |alpha|^2
                dos_k(loop_f,loop_s,3)=dos_k(loop_f,loop_s,3)+rdum*beta_sq
             end if
          end do
       end do
    end do !loop over bands

  end subroutine dos_kpt

  ! =========================================================================

  subroutine get_eig_levelspacing_k(kpt,eig,levelspacing)

    use w90_constants, only     : dp,cmplx_0,cmplx_i,twopi
    use w90_io, only            : io_error
    use w90_utility, only   : utility_diagonalize
    use w90_parameters, only    : num_wann,dos_num_points
    use w90_wanint_common, only : fourier_R_to_k,kmesh_spacing
    use w90_get_oper, only      : HH_R
    use w90_wan_ham, only   : get_deleig_a

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: eig(num_wann)
    real(kind=dp), intent(out) :: levelspacing(num_wann)

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)

    ! Adaptive smearing
    !
    real(kind=dp) :: del_eig(num_wann,3),Delta_k

    integer          :: i

    allocate(HH(num_wann,num_wann))
    allocate(delHH(num_wann,num_wann,3))
    allocate(UU(num_wann,num_wann))

    call fourier_R_to_k(kpt,HH_R,HH,0) 
    call utility_diagonalize(HH,num_wann,eig,UU) 
    call fourier_R_to_k(kpt,HH_R,delHH(:,:,1),1) 
    call fourier_R_to_k(kpt,HH_R,delHH(:,:,2),2) 
    call fourier_R_to_k(kpt,HH_R,delHH(:,:,3),3) 
    call get_deleig_a(del_eig(:,1),eig,delHH(:,:,1),UU)
    call get_deleig_a(del_eig(:,2),eig,delHH(:,:,2),UU)
    call get_deleig_a(del_eig(:,3),eig,delHH(:,:,3),UU)

    Delta_k=kmesh_spacing(dos_num_points)
    do i=1,num_wann
       levelspacing(i)=&
            sqrt(dot_product(del_eig(i,:),del_eig(i,:)))*Delta_k
    end do

  end subroutine get_eig_levelspacing_k

  ! =========================================================================

  function count_states(energy,eig,levelspacing,loop_s,npts)

    use w90_constants, only     : dp,cmplx_0,cmplx_i,twopi
    use w90_utility, only       : wgauss
    use w90_wanint_common, only : weight
    use w90_parameters, only    : num_wann,adpt_smr_steps,adpt_smr_width

    real(kind=dp) :: count_states

    ! Arguments
    !
    real(kind=dp)                  :: energy
    real(kind=dp), dimension (:,:) :: eig
    real(kind=dp), dimension (:,:) :: levelspacing
    integer                        :: loop_s
    integer                        :: npts

    ! Misc/Dummy
    !
    integer       :: loop_k,i
    real(kind=dp) :: sum,smear,arg

    count_states=0.0_dp
    do loop_k=1,npts
       sum=0.0_dp
       do i=1,num_wann
          smear=levelspacing(i,loop_k)*adpt_smr_width(loop_s)/sqrt(2.0_dp)
          arg=(energy-eig(i,loop_k))/smear
          !
          ! For Fe and a 125x125x125 interpolation mesh, E_f=12.6306 with M-P
          ! smearing, and E_f=12.6512 with F-D smearing
          !
          !          sum=sum+wgauss(arg,-99) ! Fermi-Dirac
          sum=sum+wgauss(arg,1)    ! Methfessel-Paxton case
       end do
       count_states=count_states+weight(loop_k)*sum
    end do

  end function count_states

end module w90_dos_wanint

