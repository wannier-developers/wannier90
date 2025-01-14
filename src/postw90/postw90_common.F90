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
!                                                            !
!  w90_postw90_common: routines used throughout postw90      !
!                                                            !
!------------------------------------------------------------!

module w90_postw90_common

  !! This contains the common variables and procedures needed to set up a Wannier
  !! interpolatation calculation for any physical property

  use w90_constants, only: dp
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
    set_error_input, set_error_fatal, set_error_file

  implicit none

  private

  public :: pw90common_fourier_R_to_k
  public :: pw90common_fourier_R_to_k_new
  public :: pw90common_fourier_R_to_k_new_second_d
  public :: pw90common_fourier_R_to_k_new_second_d_TB_conv
  public :: pw90common_fourier_R_to_k_vec
  public :: pw90common_fourier_R_to_k_vec_dadb
  public :: pw90common_fourier_R_to_k_vec_dadb_TB_conv
  public :: pw90common_get_occ
  public :: pw90common_kmesh_spacing
  public :: pw90common_wanint_data_dist
  public :: pw90common_wanint_get_kpoint_file
  public :: pw90common_wanint_setup
  public :: pw90common_wanint_w90_wannier90_readwrite_dist

  interface pw90common_kmesh_spacing
    module procedure kmesh_spacing_singleinteger
    module procedure kmesh_spacing_mesh
  end interface pw90common_kmesh_spacing

contains

  !================================================!
  !                   PUBLIC PROCEDURES
  ! Public procedures have names starting with wanint_
  !================================================!

  subroutine pw90common_wanint_setup(num_wann, print_output, real_lattice, mp_grid, &
                                     effective_model, ws_region, wigner_seitz, stdout, &
                                     seedname, timer, error, comm)
    !================================================!
    !
    !! Setup data ready for interpolation
    !
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: print_output_type, timer_list_type, ws_region_type
    use w90_comms, only: mpirank, w90_comm_type, comms_bcast
    use w90_postw90_types, only: wigner_seitz_type

    type(print_output_type), intent(in) :: print_output
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_region_type), intent(inout) :: ws_region

    real(kind=dp), intent(in) :: real_lattice(3, 3)
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: mp_grid(3)
    logical, intent(in) :: effective_model
    character(len=50), intent(in)  :: seedname

    integer :: ierr, ir, file_unit, num_wann_loc
    logical :: on_root = .false.

    if (mpirank(comm) == 0) on_root = .true.

    ! Find nrpts, the number of points in the Wigner-Seitz cell
    if (effective_model) then
      if (on_root) then
        ! nrpts is read from file, together with num_wann
        open (newunit=file_unit, file=trim(seedname)//'_HH_R.dat', form='formatted', &
              status='old', err=101)
        read (file_unit, *) !header
        read (file_unit, *) num_wann_loc
        if (num_wann_loc /= num_wann) then
          call set_error_fatal(error, 'Inconsistent values of num_wann in '//trim(seedname) &
                               //'_HH_R.dat and '//trim(seedname)//'.win', comm)
          return
        endif
        read (file_unit, *) wigner_seitz%nrpts
        close (file_unit)
      endif
      call comms_bcast(wigner_seitz%nrpts, 1, error, comm)
      if (allocated(error)) return
    else
      call wignerseitz(print_output, real_lattice, mp_grid, ws_region, wigner_seitz, stdout, &
                       .true., timer, error, comm)
      if (allocated(error)) return
    endif

    ! Now can allocate several arrays
    allocate (wigner_seitz%irvec(3, wigner_seitz%nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating irvec in pw90common_wanint_setup', comm)
      return
    endif
    wigner_seitz%irvec = 0
    allocate (wigner_seitz%crvec(3, wigner_seitz%nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating crvec in pw90common_wanint_setup', comm)
      return
    endif
    wigner_seitz%crvec = 0.0_dp
    allocate (wigner_seitz%ndegen(wigner_seitz%nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ndegen in pw90common_wanint_setup', comm)
      return
    endif
    wigner_seitz%ndegen = 0

    ! Also rpt_origin, so that when effective_model=.true it is not
    ! passed to get_HH_R without being initialized.
    wigner_seitz%rpt_origin = 0

    ! If effective_model, this is done in get_HH_R
    if (.not. effective_model) then
      ! Set up the lattice vectors on the Wigner-Seitz supercell
      ! where the Wannier functions live

      call wignerseitz(print_output, real_lattice, mp_grid, ws_region, wigner_seitz, stdout, &
                       .false., timer, error, comm)
      if (allocated(error)) return

      ! Convert from reduced to Cartesian coordinates

      do ir = 1, wigner_seitz%nrpts
        ! Note that 'real_lattice' stores the lattice vectors as *rows*
        wigner_seitz%crvec(:, ir) = matmul(transpose(real_lattice), wigner_seitz%irvec(:, ir))
      end do
    endif

    return

101 call set_error_file(error, 'Error in pw90common_wanint_setup: problem opening file '// &
                        trim(seedname)//'_HH_R.dat', comm)
    return !jj fixme restructure

  end subroutine pw90common_wanint_setup

  !================================================!
  subroutine pw90common_wanint_get_kpoint_file(kpoint_dist, error, comm)
    !================================================!
    !
    !! read kpoints from kpoint.dat and distribute
    !
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_date, io_time
    use w90_comms, only: mpirank, mpisize, w90_comm_type, comms_bcast
    use w90_postw90_types, only: kpoint_dist_type

    ! arguments
    type(kpoint_dist_type), intent(inout) :: kpoint_dist
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    ! local variables
    integer :: i, ierr, my_node_id, num_nodes, k_unit
    !real(kind=dp) :: sum
    logical :: on_root = .false.
    real(kind=dp), allocatable :: kt(:, :), wt(:) !temp for read/dist
    integer :: ik, nk, ir, nkloc, off

    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)

    if (my_node_id == 0) on_root = .true.

    if (on_root) then
      open (newunit=k_unit, file='kpoint.dat', status='old', form='formatted', err=106)
      read (k_unit, *) kpoint_dist%num_int_kpts
    end if
    call comms_bcast(kpoint_dist%num_int_kpts, 1, error, comm)
    if (allocated(error)) return

    allocate (kpoint_dist%num_int_kpts_on_node(0:num_nodes - 1))
    kpoint_dist%num_int_kpts_on_node(:) = kpoint_dist%num_int_kpts/num_nodes
    kpoint_dist%max_int_kpts_on_node = kpoint_dist%num_int_kpts &
                                       - (num_nodes - 1)*(kpoint_dist%num_int_kpts/num_nodes)
    kpoint_dist%num_int_kpts_on_node(0) = kpoint_dist%max_int_kpts_on_node
!   if(my_node_id < num_int_kpts- num_int_kpts_on_node*num_nodes)  num_int_kpts_on_node= num_int_kpts_on_node+1

    allocate (kpoint_dist%int_kpts(3, kpoint_dist%max_int_kpts_on_node), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating max_int_kpts_on_node in w90_wannier90_readwrite_read_um', comm)
      return
    endif
    kpoint_dist%int_kpts = 0.0_dp
    allocate (kpoint_dist%weight(kpoint_dist%max_int_kpts_on_node), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating weight in w90_wannier90_readwrite_read_um', comm)
      return
    endif
    kpoint_dist%weight = 0.0_dp

    !sum = 0.0_dp

    nk = kpoint_dist%num_int_kpts
    allocate (kt(3, nk))
    allocate (wt(nk))
    if (on_root) then
      do ik = 1, nk
        read (k_unit, *) (kt(i, ik), i=1, 3), wt(ik)
      enddo
    endif
    call comms_bcast(kt(1, 1), 3*nk, error, comm)
    if (allocated(error)) return
    call comms_bcast(wt(1), nk, error, comm)
    if (allocated(error)) return

    nkloc = kpoint_dist%num_int_kpts_on_node(my_node_id)
    off = 1
    do ir = 0, my_node_id - 1
      off = off + kpoint_dist%num_int_kpts_on_node(ir)
    enddo
    call dcopy(3*nkloc, kt(1, off), 1, kpoint_dist%int_kpts, 1)
    call dcopy(nkloc, wt(off), 1, kpoint_dist%weight, 1)
    deallocate (kt)
    deallocate (wt)

!    if (on_root) then
!      do loop_nodes = 1, num_nodes - 1
!        do loop_kpt = 1, kpoint_dist%num_int_kpts_on_node(loop_nodes)
!          read (k_unit, *) (kpoint_dist%int_kpts(i, loop_kpt), i=1, 3), kpoint_dist%weight(loop_kpt)
!          sum = sum + kpoint_dist%weight(loop_kpt)
!        end do
!
!        call comms_no_sync_send(kpoint_dist%int_kpts(1, 1), 3*kpoint_dist%num_int_kpts_on_node(loop_nodes), loop_nodes, error, comm)
!        if (allocated(error)) return
!        call comms_no_sync_send(kpoint_dist%weight(1), kpoint_dist%num_int_kpts_on_node(loop_nodes), loop_nodes, error, comm)
!        if (allocated(error)) return
!      end do
!      do loop_kpt = 1, kpoint_dist%num_int_kpts_on_node(0)
!        read (k_unit, *) (kpoint_dist%int_kpts(i, loop_kpt), i=1, 3), kpoint_dist%weight(loop_kpt)
!        sum = sum + kpoint_dist%weight(loop_kpt)
!      end do
!!       print*,'rsum',sum
!    end if
!
!    if (.not. on_root) then
!      call comms_no_sync_recv(kpoint_dist%int_kpts(1, 1), 3*kpoint_dist%num_int_kpts_on_node(my_node_id), 0, error, comm)
!      if (allocated(error)) return
!      call comms_no_sync_recv(kpoint_dist%weight(1), kpoint_dist%num_int_kpts_on_node(my_node_id), 0, error, comm)
!      if (allocated(error)) return
!    end if

    return

106 call set_error_file(error, 'Error: Problem opening file kpoint.dat in pw90common_wanint_get_kpoint_file', comm)
    return

  end subroutine pw90common_wanint_get_kpoint_file

  !================================================!
  subroutine pw90common_wanint_w90_wannier90_readwrite_dist(print_output, ws_region, kmesh_info, kpt_latt, num_kpts, &
                                                            dis_manifold, w90_system, fermi_energy_list, num_bands, &
                                                            num_wann, eigval, mp_grid, real_lattice, &
                                                            pw90_calculation, scissors_shift, effective_model, &
                                                            pw90_spin, pw90_band_deriv_degen, pw90_kpath, &
                                                            pw90_kslice, pw90_dos, pw90_berry, pw90_spin_hall, &
                                                            pw90_gyrotropic, pw90_geninterp, pw90_boltzwann, &
                                                            eig_found, error, comm)
    !================================================!
    !
    !! distribute the parameters across processors
    !! NOTE: we only send the ones postw90 uses, not all in w90
    !
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_date, io_time
    use w90_comms, only: mpirank, w90_comm_type, comms_bcast
    use w90_types
    use w90_postw90_types, only: pw90_calculation_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, pw90_kpath_mod_type, pw90_kslice_mod_type, pw90_dos_mod_type, &
      pw90_berry_mod_type, pw90_spin_hall_type, pw90_gyrotropic_type, pw90_geninterp_mod_type, &
      pw90_boltzwann_type

    type(print_output_type), intent(inout) :: print_output
    type(ws_region_type), intent(inout) :: ws_region
    type(w90_system_type), intent(inout) :: w90_system
    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(pw90_calculation_type), intent(inout) :: pw90_calculation
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    type(pw90_band_deriv_degen_type), intent(inout) :: pw90_band_deriv_degen
    type(pw90_kpath_mod_type), intent(inout) :: pw90_kpath
    type(pw90_kslice_mod_type), intent(inout) :: pw90_kslice
    type(pw90_dos_mod_type), intent(inout) :: pw90_dos
    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(pw90_spin_hall_type), intent(inout) :: pw90_spin_hall
    type(pw90_gyrotropic_type), intent(inout) :: pw90_gyrotropic
    type(pw90_geninterp_mod_type), intent(inout) :: pw90_geninterp
    type(pw90_boltzwann_type), intent(inout) :: pw90_boltzwann
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    real(kind=dp), allocatable, intent(inout) :: fermi_energy_list(:)
    real(kind=dp), pointer, intent(inout) :: eigval(:, :)
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: scissors_shift
    integer, intent(inout) :: num_kpts
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_wann
    integer, intent(inout) :: mp_grid(3)
    logical, intent(inout) :: eig_found
    logical, intent(inout) :: effective_model

    integer :: ierr
    integer :: iprintroot
    integer :: fermi_n, kdotp_nbands
    logical :: on_root = .false.

    if (mpirank(comm) == 0) on_root = .true.

    call comms_bcast(effective_model, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(eig_found, 1, error, comm)
    if (allocated(error)) return

    if (.not. effective_model) then
      call comms_bcast(mp_grid(1), 3, error, comm)
      if (allocated(error)) return
      call comms_bcast(num_kpts, 1, error, comm)
      if (allocated(error)) return
      call comms_bcast(num_bands, 1, error, comm)
      if (allocated(error)) return
    endif
    call comms_bcast(num_wann, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(print_output%timing_level, 1, error, comm)
    if (allocated(error)) return

    !______________________________________
    !JJ fixme maybe? not so pretty solution to setting iprint to zero on non-root processes
    iprintroot = print_output%iprint
    print_output%iprint = 0
    call comms_bcast(print_output%iprint, 1, error, comm)
    if (allocated(error)) return
    if (on_root) print_output%iprint = iprintroot
    !______________________________________

    call comms_bcast(ws_region%ws_distance_tol, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(ws_region%ws_search_size(1), 3, error, comm)
    if (allocated(error)) return
!    call comms_bcast(num_atoms,1)   ! Ivo: not used in postw90, right?
!    call comms_bcast(num_species,1) ! Ivo: not used in postw90, right?
    call comms_bcast(real_lattice(1, 1), 9, error, comm)
    if (allocated(error)) return
    !call comms_bcast(recip_lattice(1, 1), 9, error, comm)
    !call comms_bcast(real_metric(1, 1), 9)
    !call comms_bcast(recip_metric(1, 1), 9)
    !call comms_bcast(cell_volume, 1, error, comm)
    call comms_bcast(pw90_dos%energy_step, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%smearing%use_adaptive, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%smearing%type_index, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%kmesh%spacing, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%kmesh%mesh(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%smearing%adaptive_max_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%smearing%fixed_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%smearing%adaptive_prefactor, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%num_project, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(pw90_calculation%berry, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%task, len(pw90_berry%task), error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kmesh%spacing, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kmesh%mesh(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%curv_adpt_kmesh, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%curv_adpt_kmesh_thresh, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%curv_unit, len(pw90_berry%curv_unit), error, comm)
    if (allocated(error)) return

! Tsirkin
    call comms_bcast(pw90_calculation%gyrotropic, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%task, len(pw90_gyrotropic%task), error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%kmesh%spacing, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%kmesh%mesh(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%eigval_max, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%nfreq, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%degen_thresh, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%num_bands, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%box(1, 1), 9, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%box_corner(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%smearing%use_adaptive, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%smearing%fixed_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%smearing%type_index, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%smearing%max_arg, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(w90_system%spinors, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(pw90_spin_hall%freq_scan, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin_hall%alpha, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin_hall%beta, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin_hall%gamma, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin_hall%bandshift, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin_hall%bandshift_firstband, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin_hall%bandshift_energyshift, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin_hall%method, len(pw90_spin_hall%method), error, comm)
    if (allocated(error)) return

    call comms_bcast(pw90_berry%kubo_smearing%use_adaptive, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kubo_smearing%adaptive_prefactor, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kubo_smearing%adaptive_max_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kubo_smearing%fixed_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kubo_smearing%type_index, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kubo_eigval_max, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kubo_nfreq, 1, error, comm)
    if (allocated(error)) return
    fermi_n = 0
    if (on_root) then
      if (allocated(fermi_energy_list)) fermi_n = size(fermi_energy_list)
    endif
    call comms_bcast(fermi_n, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%energy_min, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%energy_max, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin%kmesh%spacing, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin%kmesh%mesh(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%wanint_kpoint_file, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(dis_manifold%win_min, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(dis_manifold%win_max, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%sc_eta, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%sc_w_thr, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%sc_phase_conv, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%sc_use_eta_corr, 1, error, comm)
    if (allocated(error)) return
! ----------------------------------------------
!
! New input variables in development
!
    !call comms_bcast(print_output%devel_flag, len(print_output%devel_flag), error, comm)
    call comms_bcast(pw90_calculation%spin_moment, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin%axis_polar, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_spin%axis_azimuth, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_calculation%spin_decomp, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_band_deriv_degen%use_degen_pert, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_band_deriv_degen%degen_thr, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(w90_system%num_valence_bands, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_calculation%dos, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%task, len(pw90_dos%task), error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_calculation%kpath, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kpath%task, len(pw90_kpath%task), error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kpath%bands_colour, len(pw90_kpath%bands_colour), error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_calculation%kslice, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kslice%task, len(pw90_kslice%task), error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kslice%corner(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kslice%b1(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kslice%b2(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kslice%kmesh2d(1), 2, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_kslice%fermi_lines_colour, len(pw90_kslice%fermi_lines_colour), error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%transl_inv, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(w90_system%num_elec_per_state, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(scissors_shift, 1, error, comm)
    if (allocated(error)) return

    ! Do these have to be broadcasted? (Plots done on root node only)
    !
!    call comms_bcast(bands_num_points,1)
!    call comms_bcast(bands_num_spec_points,1)
!    if(allocated(bands_spec_points)) &
!         call comms_bcast(bands_spec_points(1,1),3*bands_num_spec_points)
!    if(allocated(bands_label)) &
!         call comms_bcast(bands_label(:),len(bands_label(1))*bands_num_spec_points)
! ----------------------------------------------
    call comms_bcast(pw90_calculation%geninterp, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_geninterp%alsofirstder, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_geninterp%single_file, 1, error, comm)
    if (allocated(error)) return
    ! BoltzWann variables
    call comms_bcast(pw90_calculation%boltzwann, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%calc_also_dos, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dir_num_2d, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_energy_step, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_energy_min, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_energy_max, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_smearing%use_adaptive, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_smearing%fixed_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_smearing%adaptive_prefactor, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_smearing%adaptive_max_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%mu_min, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%mu_max, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%mu_step, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%temp_min, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%temp_max, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%temp_step, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%kmesh%spacing, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%kmesh%mesh(1), 3, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%tdf_energy_step, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%relax_time, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%tdf_smearing%use_adaptive, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%tdf_smearing%fixed_width, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%tdf_smearing%type_index, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%dos_smearing%type_index, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%bandshift, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%bandshift_firstband, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_boltzwann%bandshift_energyshift, 1, error, comm)
    if (allocated(error)) return
    call comms_bcast(ws_region%use_ws_distance, 1, error, comm)
    if (allocated(error)) return

    ! These variables are different from the ones above in that they are
    ! allocatable, and in w90_wannier90_readwrite_read they were allocated on the root node only

    if (.not. on_root) then
      allocate (fermi_energy_list(fermi_n), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating fermi_energy_read in postw90_w90_wannier90_readwrite_dist', comm)
        return
      endif
      allocate (pw90_berry%kubo_freq_list(pw90_berry%kubo_nfreq), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating kubo_freq_list in postw90_w90_wannier90_readwrite_dist', comm)
        return
      endif
      allocate (pw90_gyrotropic%band_list(pw90_gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating gyrotropic_band_list in postw90_w90_wannier90_readwrite_dist', comm)
        return
      endif
      allocate (pw90_gyrotropic%freq_list(pw90_gyrotropic%nfreq), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating gyrotropic_freq_list in postw90_w90_wannier90_readwrite_dist', comm)
        return
      endif
      allocate (pw90_dos%project(pw90_dos%num_project), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating dos_project in postw90_w90_wannier90_readwrite_dist', comm)
        return
      endif
      if (.not. effective_model) then
        if (eig_found) then
          allocate (eigval(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating eigval in postw90_w90_wannier90_readwrite_dist', comm)
            return
          endif
        end if
        allocate (kpt_latt(3, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating kpt_latt in postw90_w90_wannier90_readwrite_dist', comm)
          return
        endif
      endif
    end if

    kdotp_nbands = 0
    if (on_root) then
      if (allocated(pw90_berry%kdotp_bands)) kdotp_nbands = size(pw90_berry%kdotp_bands)
    endif
    call comms_bcast(kdotp_nbands, 1, error, comm)
    if (kdotp_nbands > 0) then
      if (.not. on_root) then
        allocate (pw90_berry%kdotp_bands(kdotp_nbands), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating kdotp_bands in postw90_param_dist', comm)
          return
        endif
      endif
      call comms_bcast(pw90_berry%kdotp_bands(1), kdotp_nbands, error, comm)
    endif

    if (fermi_n > 0) then
      call comms_bcast(fermi_energy_list(1), fermi_n, error, comm)
      if (allocated(error)) return
    endif
    call comms_bcast(pw90_gyrotropic%freq_list(1), pw90_gyrotropic%nfreq, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_gyrotropic%band_list(1), pw90_gyrotropic%num_bands, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_berry%kubo_freq_list(1), pw90_berry%kubo_nfreq, error, comm)
    if (allocated(error)) return
    call comms_bcast(pw90_dos%project(1), pw90_dos%num_project, error, comm)
    if (allocated(error)) return
    if (.not. effective_model) then
      if (eig_found) then
        call comms_bcast(eigval(1, 1), num_bands*num_kpts, error, comm)
        if (allocated(error)) return
      end if
      call comms_bcast(kpt_latt(1, 1), 3*num_kpts, error, comm)
      if (allocated(error)) return
    endif

    ! kmesh: only nntot,wb, and bk are needed to evaluate the WF matrix
    ! elements of the position operator in reciprocal space. For the
    ! extra matrix elements entering the orbital magnetization, also
    ! need nnlist. In principle could only broadcast those four variables

    if (.not. effective_model) then

      call comms_bcast(kmesh_info%nnh, 1, error, comm)
      if (allocated(error)) return
      call comms_bcast(kmesh_info%nntot, 1, error, comm)
      if (allocated(error)) return

      if (.not. on_root) then
        allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating nnlist in pw90common_wanint_w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%neigh(num_kpts, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating neigh in pw90common_wanint_w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating nncell in pw90common_wanint_w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%wb(kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating wb in pw90common_wanint_w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%bka(3, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating bka in pw90common_wanint_w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%bk(3, kmesh_info%nntot, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating bk in pw90common_wanint_w90_wannier90_readwrite_dist', comm)
          return
        endif
      end if

      call comms_bcast(kmesh_info%nnlist(1, 1), num_kpts*kmesh_info%nntot, error, comm)
      if (allocated(error)) return
      call comms_bcast(kmesh_info%neigh(1, 1), num_kpts*kmesh_info%nntot/2, error, comm)
      if (allocated(error)) return
      call comms_bcast(kmesh_info%nncell(1, 1, 1), 3*num_kpts*kmesh_info%nntot, error, comm)
      if (allocated(error)) return
      call comms_bcast(kmesh_info%wb(1), kmesh_info%nntot, error, comm)
      if (allocated(error)) return
      call comms_bcast(kmesh_info%bka(1, 1), 3*kmesh_info%nntot/2, error, comm)
      if (allocated(error)) return
      call comms_bcast(kmesh_info%bk(1, 1, 1), 3*kmesh_info%nntot*num_kpts, error, comm)
      if (allocated(error)) return

    endif

  end subroutine pw90common_wanint_w90_wannier90_readwrite_dist

  !================================================!
  subroutine pw90common_wanint_data_dist(num_wann, num_kpts, num_bands, u_matrix_opt, u_matrix, &
                                         dis_manifold, wannier_data, scissors_shift, v_matrix, &
                                         num_valence_bands, have_disentangled, error, comm)
    !================================================!
    !
    !! Distribute the um and chk files
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0
    use w90_io, only: io_date, io_time
    use w90_types, only: dis_manifold_type, wannier_data_type
    use w90_comms, only: w90_comm_type, mpirank, comms_bcast

    implicit none

    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(wannier_data_type), intent(inout) :: wannier_data
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_valence_bands
    integer, intent(in) :: num_wann, num_kpts, num_bands
    real(kind=dp), intent(in) :: scissors_shift
    complex(kind=dp), allocatable :: v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :), u_matrix(:, :, :)
    logical, intent(inout) :: have_disentangled

    integer :: ierr, loop_kpt, m, i, j
    logical :: on_root = .false.

    if (mpirank(comm) == 0) on_root = .true.

!    if (.not. on_root) then
!      ! wannier_centres is allocated in w90_wannier90_readwrite_read, so only on root node
!      ! It is then read in w90_wannier90_readwrite_read_chpkt
!      ! Therefore, now we need to allocate it on all nodes, and then broadcast it
!      allocate (wannier_data%centres(3, num_wann), stat=ierr)
!      if (ierr /= 0) then
!        call set_error_alloc(error, 'Error allocating wannier_centres in pw90common_wanint_data_dist', comm)
!        return
!      endif
!    end if
!    call comms_bcast(wannier_data%centres(1, 1), 3*num_wann, error, comm)
!    if (allocated(error)) return

    ! -------------------
    ! Ivo: added 8april11
    ! -------------------
    !
    ! Calculate the matrix that describes the combined effect of
    ! disentanglement and maximal localization. This is the combination
    ! that is most often needed for interpolation purposes
    !
    ! Allocate on all nodes
    allocate (v_matrix(num_bands, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating v_matrix in pw90common_wanint_data_dist', comm)
      return
    endif
    ! u_matrix and u_matrix_opt are stored on root only
    if (on_root) then
      if (.not. have_disentangled) then
        v_matrix(1:num_wann, :, :) = u_matrix(1:num_wann, :, :)
      else
        v_matrix = cmplx_0
        do loop_kpt = 1, num_kpts
          do j = 1, num_wann
            do m = 1, dis_manifold%ndimwin(loop_kpt)
              do i = 1, num_wann
                v_matrix(m, j, loop_kpt) = v_matrix(m, j, loop_kpt) &
                                           + u_matrix_opt(m, i, loop_kpt)*u_matrix(i, j, loop_kpt)
              enddo
            enddo
          enddo
        enddo
      endif
      if (allocated(u_matrix_opt)) deallocate (u_matrix_opt)
      if (.not. (num_valence_bands > 0 .and. abs(scissors_shift) > 1.0e-7_dp)) then
        if (allocated(u_matrix)) deallocate (u_matrix)
      endif
    endif
    call comms_bcast(v_matrix(1, 1, 1), num_bands*num_wann*num_kpts, error, comm)
    if (allocated(error)) return

    if (num_valence_bands > 0 .and. abs(scissors_shift) > 1.0e-7_dp) then
    if (.not. on_root .and. .not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating u_matrix in pw90common_wanint_data_dist', comm)
        return
      endif
    endif
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, error, comm)
    endif

    call comms_bcast(have_disentangled, 1, error, comm)
    if (allocated(error)) return

    if (have_disentangled) then
      if (.not. on_root) then

        if (.not. allocated(dis_manifold%lwindow)) then
          allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating lwindow in pw90common_wanint_data_dist', comm)
            return
          endif
        endif

        if (.not. allocated(dis_manifold%ndimwin)) then
          allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating ndimwin in pw90common_wanint_data_dist', comm)
            return
          endif
        endif

      end if

!       call comms_bcast(u_matrix_opt(1,1,1),num_bands*num_wann*num_kpts)
      call comms_bcast(dis_manifold%lwindow(1, 1), num_bands*num_kpts, error, comm)
      if (allocated(error)) return
      call comms_bcast(dis_manifold%ndimwin(1), num_kpts, error, comm)
      if (allocated(error)) return
    end if

  end subroutine pw90common_wanint_data_dist

!================================================

  subroutine pw90common_get_occ(ef, eig, occ, num_wann)
    !================================================!
    !
    !! Compute the electronic occupancy
    !
    !================================================!

    use w90_constants, only: dp

    ! arguments
    integer, intent(in) :: num_wann

    real(kind=dp), intent(in)  :: eig(num_wann)
    !! Eigenvalues
    real(kind=dp), intent(in)  :: ef
    !! Fermi level
    real(kind=dp), intent(out) :: occ(num_wann)
    !! Occupancy of states

    ! local variables
    integer :: i

    ! State occupancies
!    if(smear_temp < eps7) then
    !
    ! Use a step function occupancy (T=0)
    !
    occ(:) = 0.0_dp
    do i = 1, num_wann
      if (eig(i) < ef) occ(i) = 1.0_dp
    end do
!    else
    !
    ! Use a Fermi-Dirac occupancy (T=smear_temp, in Kelvin)
    !
    ! k_B.T in electron-volts
    !
!       kt=k_B_SI*smear_temp/elem_charge_SI
!       do i=1,num_wann
!          occ(i)=1.0_dp/(exp((eig(i)-ef)/kt)+1.0_dp)
!       end do
!    end if

  end subroutine pw90common_get_occ

  !================================================
  function kmesh_spacing_singleinteger(num_points, recip_lattice)
    !================================================
    !! Set up the value of the interpolation mesh spacing, needed for
    !! adaptive smearing [see Eqs. (34-35) YWVS07]. Choose it as the largest of
    !! the three Delta_k's for each of the primitive translations b1, b2, and b3
    !================================================

    integer, intent(in) :: num_points
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp)       :: kmesh_spacing_singleinteger

    integer        :: i
    real(kind=dp) :: Delta_k_i(3)

    ! NOTE: The vectors b_i are stored as *rows* in recip_lattice (argh!).
    ! Hence I believe Jonathan's original code confused rows with columns
    ! when computing Delta_k, which he called 'rspace'
    ! (See my e-mail of 20Sept07)
    !
    do i = 1, 3
      Delta_k_i(i) = sqrt(dot_product(recip_lattice(i, :), recip_lattice(i, :)))/num_points
    end do
    kmesh_spacing_singleinteger = maxval(Delta_k_i)

  end function kmesh_spacing_singleinteger
  function kmesh_spacing_mesh(mesh, recip_lattice)
    !! Same as kmesh_spacing_singleinteger, but for a kmesh with three
    !! different mesh samplings along the three directions

    integer, dimension(3), intent(in) :: mesh
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp) :: kmesh_spacing_mesh

    integer        :: i
    real(kind=dp) :: Delta_k_i(3)

    do i = 1, 3
      Delta_k_i(i) = sqrt(dot_product(recip_lattice(i, :), recip_lattice(i, :)))/mesh(i)
    end do
    kmesh_spacing_mesh = maxval(Delta_k_i)

  end function kmesh_spacing_mesh

  !================================================!
  subroutine pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, OO, &
                                       OO_R, kpt, real_lattice, mp_grid, alpha, num_wann, &
                                       error, comm)
    !================================================!
    !
    !! For alpha=0:
    !! O_ij(R) --> O_ij(k) = sum_R e^{+ik.R}*O_ij(R)
    !!
    !! For alpha=1,2,3:
    !! sum_R [cmplx_i*R_alpha*e^{+ik.R}*O_ij(R)]
    !! where R_alpha is a Cartesian component of R
    !! ***REMOVE EVENTUALLY*** (replace with pw90common_fourier_R_to_k_new)
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_types, only: wannier_data_type, ws_region_type, ws_distance_type
    use w90_ws_distance, only: ws_translate_dist
    use w90_postw90_types, only: wigner_seitz_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: alpha

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), intent(out) :: OO(:, :)

    ! local variables
    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (ws_region%use_ws_distance) then
      call ws_translate_dist(ws_distance, ws_region, num_wann, wannier_data%centres, real_lattice, &
                             mp_grid, wigner_seitz%nrpts, wigner_seitz%irvec, error, comm)
      if (allocated(error)) return
    endif

    OO(:, :) = cmplx_0
    do ir = 1, wigner_seitz%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (ws_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(wigner_seitz%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (alpha == 0) then
              OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
              OO(i, j) = OO(i, j) + cmplx_i*ws_distance%crdist(alpha, ideg, i, j, ir) &
                         *phase_fac*OO_R(i, j, ir)
            else
              stop 'wrong value of alpha in pw90common_fourier_R_to_k'
            endif
          enddo
        enddo
        enddo
      else
        ! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), wigner_seitz%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir), dp)
        if (alpha == 0) then
          OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
          OO(:, :) = OO(:, :) + &
                     cmplx_i*wigner_seitz%crvec(alpha, ir)*phase_fac*OO_R(:, :, ir)
        else
          stop 'wrong value of alpha in pw90common_fourier_R_to_k'
        endif
      endif

    enddo

  end subroutine pw90common_fourier_R_to_k

  !================================================!
  subroutine pw90common_fourier_R_to_k_new(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           OO_R, kpt, real_lattice, mp_grid, num_wann, error, &
                                           comm, OO, OO_dx, OO_dy, OO_dz)
    !================================================!
    !
    !! For OO:
    !! $$O_{ij}(k) = \sum_R e^{+ik.R}.O_{ij}(R)$$
    !! For $$OO_{dx,dy,dz}$$:
    !! $$\sum_R [i.R_{dx,dy,dz}.e^{+ik.R}.O_{ij}(R)]$$
    !! where R_{x,y,z} are the Cartesian components of R
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_types, only: ws_region_type, wannier_data_type, ws_distance_type
    use w90_ws_distance, only: ws_translate_dist
    use w90_postw90_types, only: wigner_seitz_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO(:, :)
    complex(kind=dp), optional, intent(out) :: OO_dx(:, :)
    complex(kind=dp), optional, intent(out) :: OO_dy(:, :)
    complex(kind=dp), optional, intent(out) :: OO_dz(:, :)

    ! local variables
    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (ws_region%use_ws_distance) then
      call ws_translate_dist(ws_distance, ws_region, num_wann, wannier_data%centres, real_lattice, &
                             mp_grid, wigner_seitz%nrpts, wigner_seitz%irvec, error, comm)
      if (allocated(error)) return
    endif

    if (present(OO)) OO = cmplx_0
    if (present(OO_dx)) OO_dx = cmplx_0
    if (present(OO_dy)) OO_dy = cmplx_0
    if (present(OO_dz)) OO_dz = cmplx_0
    do ir = 1, wigner_seitz%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (ws_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(wigner_seitz%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_dx)) OO_dx(i, j) = OO_dx(i, j) + &
                                              cmplx_i*ws_distance%crdist(1, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
            if (present(OO_dy)) OO_dy(i, j) = OO_dy(i, j) + &
                                              cmplx_i*ws_distance%crdist(2, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
            if (present(OO_dz)) OO_dz(i, j) = OO_dz(i, j) + &
                                              cmplx_i*ws_distance%crdist(3, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), wigner_seitz%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir), dp)
        if (present(OO)) OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        if (present(OO_dx)) OO_dx(:, :) = OO_dx(:, :) + &
                                          cmplx_i*wigner_seitz%crvec(1, ir)*phase_fac*OO_R(:, :, ir)
        if (present(OO_dy)) OO_dy(:, :) = OO_dy(:, :) + &
                                          cmplx_i*wigner_seitz%crvec(2, ir)*phase_fac*OO_R(:, :, ir)
        if (present(OO_dz)) OO_dz(:, :) = OO_dz(:, :) + &
                                          cmplx_i*wigner_seitz%crvec(3, ir)*phase_fac*OO_R(:, :, ir)
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new

  !================================================!
  subroutine pw90common_fourier_R_to_k_new_second_d(kpt, OO_R, num_wann, ws_region, wannier_data, &
                                                    real_lattice, mp_grid, ws_distance, &
                                                    wigner_seitz, error, comm, OO, OO_da, OO_dadb)
    !================================================!
    !
    !! For OO:
    !! $$O_{ij}(k) = \sum_R e^{+ik.R}.O_{ij}(R)$$
    !! For $$OO_{dx,dy,dz}$$:
    !! $$\sum_R [i.R_{dx,dy,dz}.e^{+ik.R}.O_{ij}(R)]$$
    !! where R_{x,y,z} are the Cartesian components of R
    !! For $$OO_{dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$-\sum_R [R_{dx1,dy1,dz1}.R_{dx2,dy2,dz2}.e^{+ik.R}.O_{ij}(R)]$$
    !! where R_{xi,yi,zi} are the Cartesian components of R
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_types, only: ws_region_type, wannier_data_type, ws_distance_type
    use w90_ws_distance, only: ws_translate_dist
    use w90_postw90_types, only: wigner_seitz_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO(:, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_dadb(:, :, :, :)

    ! local variables
    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (ws_region%use_ws_distance) then
      call ws_translate_dist(ws_distance, ws_region, num_wann, wannier_data%centres, real_lattice, &
                             mp_grid, wigner_seitz%nrpts, wigner_seitz%irvec, error, comm)
      if (allocated(error)) return
    endif

    if (present(OO)) OO = cmplx_0
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, wigner_seitz%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (ws_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(wigner_seitz%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i*ws_distance%crdist(a, ideg, i, j, ir)* &
                                 phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) &
                                        - ws_distance%crdist(a, ideg, i, j, ir) &
                                        *ws_distance%crdist(b, ideg, i, j, ir) &
                                        *phase_fac*OO_R(i, j, ir)
                enddo
              enddo
            end if

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), wigner_seitz%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir), dp)
        if (present(OO)) OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        if (present(OO_da)) then
          do a = 1, 3
            OO_da(:, :, a) = OO_da(:, :, a) + cmplx_i*wigner_seitz%crvec(a, ir)*phase_fac &
                             *OO_R(:, :, ir)
          enddo
        endif
        if (present(OO_dadb)) then
          do a = 1, 3
            do b = 1, 3
              OO_dadb(:, :, a, b) = OO_dadb(:, :, a, b) - &
                                    wigner_seitz%crvec(a, ir)*wigner_seitz%crvec(b, ir)*phase_fac &
                                    *OO_R(:, :, ir)
            enddo
          enddo
        end if
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new_second_d

  !================================================!
  subroutine pw90common_fourier_R_to_k_new_second_d_TB_conv(kpt, OO_R, oo_a_R, num_wann, &
                                                            ws_region, wannier_data, real_lattice, &
                                                            mp_grid, ws_distance, wigner_seitz, &
                                                            error, comm, OO, OO_da, OO_dadb)
    !================================================!
    ! modified version of pw90common_fourier_R_to_k_new_second_d, includes wannier centres in
    ! the exponential inside the sum (so called TB convention)
    !
    !! For OO:
    !! $$O_{ij}(k) = \sum_R e^{+ik.(R+tau_ij)}.O_{ij}(R)$$
    !! For $$OO_{dx,dy,dz}$$:
    !! $$\sum_R [i.(R+tau_ij)_{dx,dy,dz}.e^{+ik.(R+tau_ij)}.O_{ij}(R)]$$
    !! where R_{x,y,z} are the Cartesian components of R
    !! For $$OO_{dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$-\sum_R [(R+tau_ij)_{dx1,dy1,dz1}.(R+tau_ij)_{dx2,dy2,dz2}.e^{+ik.(R+tau_ij)}.O_{ij}(R)]$$
    !! where {xi,yi,zi} denote the Cartesian components and
    !  tau_ij = tau_j - tau_i, being tau_i=<0i|r|0i> the individual wannier centres
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_types, only: ws_region_type, wannier_data_type, ws_distance_type
    use w90_ws_distance, only: ws_translate_dist
    use w90_utility, only: utility_cart_to_frac, utility_inverse_mat
    use w90_postw90_types, only: wigner_seitz_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: oo_a_R(:, :, :, :)
    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO(:, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_dadb(:, :, :, :)

    ! local variables
    real(kind=dp) :: inv_lattice(3, 3)
    integer :: ir, i, j, ideg, a, b
    real(kind=dp) :: rdotk
    real(kind=dp) :: local_wannier_centres(3, num_wann), wannier_centres_frac(3, num_wann)
    real(kind=dp) :: r_sum(3)
    complex(kind=dp) :: phase_fac

    r_sum = 0.d0

    if (ws_region%use_ws_distance) then
      call ws_translate_dist(ws_distance, ws_region, num_wann, wannier_data%centres, real_lattice, &
                             mp_grid, wigner_seitz%nrpts, wigner_seitz%irvec, error, comm)
      if (allocated(error)) return
    endif

    ! calculate wannier centres in cartesian
    local_wannier_centres(:, :) = 0.d0
    do j = 1, num_wann
      do ir = 1, wigner_seitz%nrpts
        if ((wigner_seitz%irvec(1, ir) .eq. 0) .and. (wigner_seitz%irvec(2, ir) .eq. 0) .and. &
            (wigner_seitz%irvec(3, ir) .eq. 0)) then
          local_wannier_centres(1, j) = real(oo_a_R(j, j, ir, 1))
          local_wannier_centres(2, j) = real(oo_a_R(j, j, ir, 2))
          local_wannier_centres(3, j) = real(oo_a_R(j, j, ir, 3))
        endif
      enddo
    enddo
    ! rotate wannier centres from cartesian to fractional coordinates
    wannier_centres_frac(:, :) = 0.d0
    call utility_inverse_mat(real_lattice, inv_lattice)
    do ir = 1, num_wann
      call utility_cart_to_frac(local_wannier_centres(:, ir), wannier_centres_frac(:, ir), &
                                inv_lattice)
    enddo

    if (present(OO)) OO = cmplx_0
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, wigner_seitz%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (ws_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir) + &
                                                   wannier_centres_frac(:, j) - wannier_centres_frac(:, i), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(wigner_seitz%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i* &
                                 (ws_distance%crdist(a, ideg, i, j, ir) &
                                  + local_wannier_centres(a, j) &
                                  - local_wannier_centres(a, i))*phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) - &
                                        (ws_distance%crdist(a, ideg, i, j, ir) &
                                         + local_wannier_centres(a, j) &
                                         - local_wannier_centres(a, i)) &
                                        *(ws_distance%crdist(b, ideg, i, j, ir) &
                                          + local_wannier_centres(b, j) - &
                                          local_wannier_centres(b, i))*phase_fac*OO_R(i, j, ir)
                enddo
              enddo
            end if

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        do j = 1, num_wann
          do i = 1, num_wann
            r_sum(:) = real(wigner_seitz%irvec(:, ir)) &
                       + wannier_centres_frac(:, j) - wannier_centres_frac(:, i)
            rdotk = twopi*dot_product(kpt(:), r_sum(:))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i* &
                                 (wigner_seitz%crvec(a, ir) + local_wannier_centres(a, j) - &
                                  local_wannier_centres(a, i))*phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = &
                    OO_dadb(i, j, a, b) - &
                    (wigner_seitz%crvec(a, ir) + local_wannier_centres(a, j) - &
                     local_wannier_centres(a, i))*(wigner_seitz%crvec(b, ir) + &
                                                   local_wannier_centres(b, j) - local_wannier_centres(b, i))* &
                    phase_fac*OO_R(i, j, ir)
                enddo
              enddo
            end if
          enddo
        enddo
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new_second_d_TB_conv

  !================================================!
  subroutine pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           OO_R, kpt, real_lattice, mp_grid, num_wann, error, &
                                           comm, OO_true, OO_pseudo)
    !================================================!
    !
    !! For OO_true (true vector):
    !! $${\vec O}_{ij}(k) = \sum_R e^{+ik.R} {\vec O}_{ij}(R)$$
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_types, only: ws_region_type, wannier_data_type, ws_distance_type
    use w90_ws_distance, only: ws_translate_dist
    use w90_postw90_types, only: wigner_seitz_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in)  :: OO_R(:, :, :, :)
    complex(kind=dp), optional, intent(out) :: OO_true(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_pseudo(:, :, :)

    ! local variables
    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (ws_region%use_ws_distance) then
      call ws_translate_dist(ws_distance, ws_region, num_wann, wannier_data%centres, real_lattice, &
                             mp_grid, wigner_seitz%nrpts, wigner_seitz%irvec, error, comm)
      if (allocated(error)) return
    endif

    if (present(OO_true)) OO_true = cmplx_0
    if (present(OO_pseudo)) OO_pseudo = cmplx_0
    do ir = 1, wigner_seitz%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (ws_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir) &
                                                               *ws_distance%ndeg(i, j, ir), dp)
            if (present(OO_true)) then
              OO_true(i, j, 1) = OO_true(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_true(i, j, 2) = OO_true(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_true(i, j, 3) = OO_true(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
            if (present(OO_pseudo)) then
              OO_pseudo(i, j, 1) = OO_pseudo(i, j, 1) &
                                   + cmplx_i*ws_distance%crdist(2, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 3) &
                                   - cmplx_i*ws_distance%crdist(3, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 2)
              OO_pseudo(i, j, 2) = OO_pseudo(i, j, 2) &
                                   + cmplx_i*ws_distance%crdist(3, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 1) &
                                   - cmplx_i*ws_distance%crdist(1, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 3)
              OO_pseudo(i, j, 3) = OO_pseudo(i, j, 3) &
                                   + cmplx_i*ws_distance%crdist(1, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 2) &
                                   - cmplx_i*ws_distance%crdist(2, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 1)
            endif
          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), wigner_seitz%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir), dp)
        if (present(OO_true)) then
          OO_true(:, :, 1) = OO_true(:, :, 1) + phase_fac*OO_R(:, :, ir, 1)
          OO_true(:, :, 2) = OO_true(:, :, 2) + phase_fac*OO_R(:, :, ir, 2)
          OO_true(:, :, 3) = OO_true(:, :, 3) + phase_fac*OO_R(:, :, ir, 3)
        endif
        if (present(OO_pseudo)) then
          OO_pseudo(:, :, 1) = OO_pseudo(:, :, 1) &
                               + cmplx_i*wigner_seitz%crvec(2, ir)*phase_fac*OO_R(:, :, ir, 3) &
                               - cmplx_i*wigner_seitz%crvec(3, ir)*phase_fac*OO_R(:, :, ir, 2)
          OO_pseudo(:, :, 2) = OO_pseudo(:, :, 2) &
                               + cmplx_i*wigner_seitz%crvec(3, ir)*phase_fac*OO_R(:, :, ir, 1) &
                               - cmplx_i*wigner_seitz%crvec(1, ir)*phase_fac*OO_R(:, :, ir, 3)
          OO_pseudo(:, :, 3) = OO_pseudo(:, :, 3) &
                               + cmplx_i*wigner_seitz%crvec(1, ir)*phase_fac*OO_R(:, :, ir, 2) &
                               - cmplx_i*wigner_seitz%crvec(2, ir)*phase_fac*OO_R(:, :, ir, 1)
        endif
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec

  !================================================!
  subroutine pw90common_fourier_R_to_k_vec_dadb(ws_region, wannier_data, ws_distance, &
                                                wigner_seitz, OO_R, kpt, real_lattice, mp_grid, &
                                                num_wann, error, comm, OO_da, OO_dadb)
    !================================================!
    !
    !! For $$OO_{ij;dx,dy,dz}$$:
    !! $$O_{ij;dx,dy,dz}(k) = \sum_R e^{+ik.R} O_{ij;dx,dy,dz}(R)$$
    !! For $$OO_{ij;dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$O_{ij;dx1,dy1,dz1;dx2,dy2,dz2}(k) = \sum_R e^{+ik.R} i.R_{dx2,dy2,dz2}
    !!                                       .O_{ij;dx1,dy1,dz1}(R)$$
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_types, only: ws_region_type, wannier_data_type, ws_distance_type
    use w90_ws_distance, only: ws_translate_dist
    use w90_postw90_types, only: wigner_seitz_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in)  :: OO_R(:, :, :, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_dadb(:, :, :, :)

    ! local variables
    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (ws_region%use_ws_distance) then
      call ws_translate_dist(ws_distance, ws_region, num_wann, wannier_data%centres, real_lattice, &
                             mp_grid, wigner_seitz%nrpts, wigner_seitz%irvec, error, comm)
      if (allocated(error)) return
    endif

    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, wigner_seitz%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (ws_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(wigner_seitz%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO_da)) then
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) &
                                        + cmplx_i*ws_distance%crdist(b, ideg, i, j, ir) &
                                        *phase_fac*OO_R(i, j, ir, a)
                enddo
              enddo
            endif

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), wigner_seitz%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir), dp)
        if (present(OO_da)) then
          OO_da(:, :, 1) = OO_da(:, :, 1) + phase_fac*OO_R(:, :, ir, 1)
          OO_da(:, :, 2) = OO_da(:, :, 2) + phase_fac*OO_R(:, :, ir, 2)
          OO_da(:, :, 3) = OO_da(:, :, 3) + phase_fac*OO_R(:, :, ir, 3)
        endif
        if (present(OO_dadb)) then
          do a = 1, 3
            do b = 1, 3
              OO_dadb(:, :, a, b) = OO_dadb(:, :, a, b) &
                                    + cmplx_i*wigner_seitz%crvec(b, ir)*phase_fac*OO_R(:, :, ir, a)
            enddo
          enddo
        endif
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec_dadb

  !================================================!
  subroutine pw90common_fourier_R_to_k_vec_dadb_TB_conv(ws_region, wannier_data, ws_distance, &
                                                        wigner_seitz, OO_R, kpt, real_lattice, &
                                                        mp_grid, num_wann, error, comm, OO_da, &
                                                        OO_dadb)
    !================================================!
    !
    ! modified version of pw90common_fourier_R_to_k_vec_dadb, includes wannier centres in
    ! the exponential inside the sum (so called TB convention)
    !
    !! For $$OO_{ij;dx,dy,dz}$$:
    !! $$O_{ij;dx,dy,dz}(k) = \sum_R e^{+ik.(R+tau_ij)} O_{ij;dx,dy,dz}(R)$$
    !! For $$OO_{ij;dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$O_{ij;dx1,dy1,dz1;dx2,dy2,dz2}(k) = \sum_R e^{+ik.(R+tau_ij)} i.(R+tau_ij)_{dx2,dy2,dz2}
    !!                                       .O_{ij;dx1,dy1,dz1}(R)$$
    ! with tau_ij = tau_j - tau_i, being tau_i=<0i|r|0i> the individual wannier centres
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_types, only: ws_region_type, wannier_data_type, ws_distance_type
    use w90_ws_distance, only: ws_translate_dist
    use w90_utility, only: utility_cart_to_frac, utility_inverse_mat
    use w90_postw90_types, only: wigner_seitz_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: OO_R(:, :, :, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out)   :: OO_dadb(:, :, :, :)

    ! local variables
    real(kind=dp)    :: inv_lattice(3, 3)
    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac
    real(kind=dp)    :: local_wannier_centres(3, num_wann), wannier_centres_frac(3, num_wann)
    real(kind=dp)                                                 :: r_sum(3)

    r_sum = 0.d0

    if (ws_region%use_ws_distance) then
      call ws_translate_dist(ws_distance, ws_region, num_wann, wannier_data%centres, real_lattice, &
                             mp_grid, wigner_seitz%nrpts, wigner_seitz%irvec, error, comm)
      if (allocated(error)) return
    endif

    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0

    ! calculate wannier centres in cartesian
    local_wannier_centres(:, :) = 0.d0
    do j = 1, num_wann
      do ir = 1, wigner_seitz%nrpts
        if ((wigner_seitz%irvec(1, ir) .eq. 0) .and. (wigner_seitz%irvec(2, ir) .eq. 0) &
            .and. (wigner_seitz%irvec(3, ir) .eq. 0)) then
          local_wannier_centres(1, j) = real(OO_R(j, j, ir, 1))
          local_wannier_centres(2, j) = real(OO_R(j, j, ir, 2))
          local_wannier_centres(3, j) = real(OO_R(j, j, ir, 3))
        endif
      enddo
    enddo
    ! rotate wannier centres from cartesian to fractional coordinates
    wannier_centres_frac(:, :) = 0.d0
    call utility_inverse_mat(real_lattice, inv_lattice)
    do ir = 1, num_wann
      call utility_cart_to_frac(local_wannier_centres(:, ir), wannier_centres_frac(:, ir), &
                                inv_lattice)
    enddo

!    print *, 'wannier_centres_frac'
!    do ir = 1,num_wann
!     print *, wannier_centres_frac(:,ir)
!    enddo
!    stop
!
!    print *, 'crvec'
!    do ir = 1,nrpts
!     print *, crvec(:,ir)
!    enddo
!    stop
!    print *, 'wannier_centres'
!    do ir = 1,num_wann
!     print *, wannier_centres(:,ir)
!    enddo
!    stop

    do ir = 1, wigner_seitz%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (ws_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir) + &
                                                   wannier_centres_frac(:, j) &
                                                   - wannier_centres_frac(:, i), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(wigner_seitz%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO_da)) then
              ! if we are at the origin and at the same band, then the
              ! matrix element is zero in this convention
              if ((wigner_seitz%irvec(1, ir) .eq. 0) .and. &
                  (wigner_seitz%irvec(2, ir) .eq. 0) .and. &
                  (wigner_seitz%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
                cycle
              else
                OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
                OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
                OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
              endif
            endif
            if (present(OO_dadb)) then
              ! same skip as before
              if ((wigner_seitz%irvec(1, ir) .eq. 0) .and. &
                  (wigner_seitz%irvec(2, ir) .eq. 0) .and. &
                  (wigner_seitz%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
                cycle
              else
                do a = 1, 3
                  do b = 1, 3
                    OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + cmplx_i* &
                                          (ws_distance%crdist(b, ideg, i, j, ir) &
                                           + local_wannier_centres(b, j) &
                                           - local_wannier_centres(b, i)) &
                                          *phase_fac*OO_R(i, j, ir, a)
                  enddo
                enddo
              endif
            endif

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        do j = 1, num_wann
        do i = 1, num_wann
          r_sum(:) = real(wigner_seitz%irvec(:, ir)) &
                     + wannier_centres_frac(:, j) - wannier_centres_frac(:, i)
          rdotk = twopi*dot_product(kpt(:), r_sum(:))
          phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wigner_seitz%ndegen(ir), dp)
          if (present(OO_da)) then
            ! if we are at the origin and at the same band, then the
            ! matrix element is zero in this convention
            if ((wigner_seitz%irvec(1, ir) .eq. 0) .and. (wigner_seitz%irvec(2, ir) .eq. 0) .and. &
                (wigner_seitz%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*(OO_R(i, j, ir, 1) &
                                                           - local_wannier_centres(1, j))
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*(OO_R(i, j, ir, 2) &
                                                           - local_wannier_centres(2, j))
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*(OO_R(i, j, ir, 3) &
                                                           - local_wannier_centres(3, j))
!            print *, 'OO_R(i,j,ir,1)', OO_R(i,j,ir,1)
!            print *, 'local_wannier_centres(1,j)', local_wannier_centres(1,j)
!            print *, 'OO_R(i,j,ir,2)', OO_R(i,j,ir,2)
!            print *, 'local_wannier_centres(2,j)', local_wannier_centres(2,j)
              cycle
            else
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
          endif
          if (present(OO_dadb)) then
            ! same skip as before
            if ((wigner_seitz%irvec(1, ir) .eq. 0) .and. (wigner_seitz%irvec(2, ir) .eq. 0) .and. &
                (wigner_seitz%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + &
                                        cmplx_i*(wigner_seitz%crvec(b, ir) + &
                                                 local_wannier_centres(b, j) &
                                                 - local_wannier_centres(b, i))*phase_fac* &
                                        (OO_R(i, j, ir, a) - local_wannier_centres(a, j))
                enddo
              enddo
!           cycle
            else
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + &
                                        cmplx_i*(wigner_seitz%crvec(b, ir) + &
                                                 local_wannier_centres(b, j) &
                                                 - local_wannier_centres(b, i))*phase_fac*OO_R(i, j, ir, a)
                enddo
              enddo
            endif
          endif
        enddo
        enddo
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec_dadb_TB_conv

  !================================================!
  !                   PRIVATE PROCEDURES
  !================================================!

  !================================================!
  subroutine wignerseitz(print_output, real_lattice, mp_grid, ws_region, wigner_seitz, stdout, &
                         count_pts, timer, error, comm)
    !================================================!
    !! Calculates a grid of lattice vectors r that fall inside (and eventually
    !! on the surface of) the Wigner-Seitz supercell centered on the
    !! origin of the Bravais superlattice with primitive translations
    !! mp_grid(1)*a_1, mp_grid(2)*a_2, and mp_grid(3)*a_3
    !================================================!

    use w90_constants, only: eps8, dp
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: print_output_type, timer_list_type, ws_region_type
    use w90_utility, only: utility_metric
    use w90_comms, only: w90_comm_type, mpirank
    use w90_postw90_types, only: wigner_seitz_type

    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points

    implicit none
    ! arguments
    type(ws_region_type), intent(in)    :: ws_region
    type(print_output_type), intent(in) :: print_output
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout
    logical, intent(in) :: count_pts
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    ! local variables
    integer       :: ierr, dist_dim
    integer       :: n1, n2, n3, i1, i2, i3, icnt, i, j, ir
    integer       :: ndiff(3)
    logical :: on_root = .false.
    real(kind=dp), allocatable :: dist(:)
    real(kind=dp) :: real_metric(3, 3)
    real(kind=dp) :: tot, dist_min

    if (mpirank(comm) == 0) on_root = .true.

    if (print_output%timing_level > 1 .and. on_root) &
      call io_stopwatch_start('postw90_common: wigner_seitz', timer)

    call utility_metric(real_lattice, real_metric)

    dist_dim = 1
    do i = 1, 3
      dist_dim = dist_dim*((ws_region%ws_search_size(i) + 1)*2 + 1)
    end do
    allocate (dist(dist_dim), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating dist in wigner_seitz', comm)
      return
    endif

    ! The Wannier functions live in a periodic supercell of the real space unit
    ! cell. This supercell is mp_grid(i) unit cells long along each primitive
    ! translation vector a_i of the unit cell
    !
    ! We loop over grid points r on a unit cell that is (2*ws_search_size+1)**3 times
    ! larger than this primitive supercell.
    !
    ! One of these points is in the W-S supercell if it is closer to R=0 than any
    ! of the other points R (where R are the translation vectors of the
    ! supercell).

    ! In the end, nrpts contains the total number of grid points that have been
    ! found in the Wigner-Seitz cell

    wigner_seitz%nrpts = 0
    ! Loop over the lattice vectors of the primitive cell
    ! that live in a supercell which is (2*ws_search_size+1)**2
    ! larger than the Born-von Karman supercell.
    ! We need to find which among these live in the Wigner-Seitz cell
    do n1 = -ws_region%ws_search_size(1)*mp_grid(1), ws_region%ws_search_size(1)*mp_grid(1)
      do n2 = -ws_region%ws_search_size(2)*mp_grid(2), ws_region%ws_search_size(2)*mp_grid(2)
        do n3 = -ws_region%ws_search_size(3)*mp_grid(3), ws_region%ws_search_size(3)*mp_grid(3)
          ! Loop over the lattice vectors R of the Born-von Karman supercell
          ! that contains all the points of the previous loop.
          ! There are (2*(ws_search_size+1)+1)**3 points R. R=0 corresponds to
          ! i1=i2=i3=0, or icnt=((2*(ws_search_size+1)+1)**3 + 1)/2
          icnt = 0
          do i1 = -ws_region%ws_search_size(1) - 1, ws_region%ws_search_size(1) + 1
            do i2 = -ws_region%ws_search_size(2) - 1, ws_region%ws_search_size(2) + 1
              do i3 = -ws_region%ws_search_size(3) - 1, ws_region%ws_search_size(3) + 1
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1*mp_grid(1)
                ndiff(2) = n2 - i2*mp_grid(2)
                ndiff(3) = n3 - i3*mp_grid(3)
                dist(icnt) = 0.0_dp
                do i = 1, 3
                  do j = 1, 3
                    dist(icnt) = dist(icnt) + &
                                 real(ndiff(i), dp)*real_metric(i, j)*real(ndiff(j), dp)
                  enddo
                enddo
              enddo
            enddo
          enddo
          dist_min = minval(dist)
          if (abs(dist((dist_dim + 1)/2) - dist_min) .lt. ws_region%ws_distance_tol**2) then
            wigner_seitz%nrpts = wigner_seitz%nrpts + 1
            if (.not. count_pts) then
              wigner_seitz%ndegen(wigner_seitz%nrpts) = 0
              do i = 1, dist_dim
                if (abs(dist(i) - dist_min) .lt. ws_region%ws_distance_tol**2) &
                  wigner_seitz%ndegen(wigner_seitz%nrpts) = &
                  wigner_seitz%ndegen(wigner_seitz%nrpts) + 1
              end do
              wigner_seitz%irvec(1, wigner_seitz%nrpts) = n1
              wigner_seitz%irvec(2, wigner_seitz%nrpts) = n2
              wigner_seitz%irvec(3, wigner_seitz%nrpts) = n3

              ! Remember which grid point r is at the origin

              if (n1 == 0 .and. n2 == 0 .and. n3 == 0) wigner_seitz%rpt_origin = wigner_seitz%nrpts
            endif
          end if

          !n3
        enddo
        !n2
      enddo
      !n1
    enddo
    !
    deallocate (dist, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating dist wigner_seitz', comm)
      return
    end if

    if (count_pts) then
      if (print_output%timing_level > 1 .and. on_root) &
        call io_stopwatch_stop('postw90_common: wigner_seitz', timer)
      return
    end if

    if (print_output%iprint >= 3 .and. on_root) then
      write (stdout, '(1x,i4,a,/)') wigner_seitz%nrpts, &
        ' lattice points in Wigner-Seitz supercell:'
      do ir = 1, wigner_seitz%nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', wigner_seitz%irvec(1, ir), &
          wigner_seitz%irvec(2, ir), wigner_seitz%irvec(3, ir), '  degeneracy: ', &
          wigner_seitz%ndegen(ir)
      enddo
      write (stdout, '(1x,a,f12.3)') ' tot = ', tot
      write (stdout, '(1x,a,i12)') ' mp_grid product = ', mp_grid(1)*mp_grid(2)*mp_grid(3)
    endif
    ! Check the "sum rule"
    tot = 0.0_dp
    do ir = 1, wigner_seitz%nrpts
      !
      ! Corrects weights in Fourier sums for R-vectors on the boundary of the
      ! W-S supercell
      !
      tot = tot + 1.0_dp/real(wigner_seitz%ndegen(ir), dp)
    enddo
    if (abs(tot - real(mp_grid(1)*mp_grid(2)*mp_grid(3), dp)) > eps8) then
      call set_error_fatal(error, 'ERROR in wigner_seitz: error in finding Wigner-Seitz points', comm)
      return
    endif
    if (print_output%timing_level > 1 .and. on_root) &
      call io_stopwatch_stop('postw90_common: wigner_seitz', timer)
    return
  end subroutine wignerseitz

end module w90_postw90_common
