
module w90_helper_types

  ! as with fortran routines like allocate, the status variable indicates an error if non-zero
  ! positive is an error, negative is a warning (such as non-convergence) which is recoverable.

  ! note on the naming of output and error-output streams
  !   the specific names stdout and stderr cannot be used because
  !   the python wrapper uses an intermediate c representation with named arguments
  !   in this wrapper stdout/err are already defined as structs inconsistent
  !   which are inconsistent with the integer unit numbers used here.
  !   (The slightly 77 style naming here has no particula significance.)

  use w90_constants
  use w90_types
  use w90_wannier90_types

  implicit none

  ! datatype encapsulating types of use in both wannier90 and postw90
  type lib_global_type
    type(atom_data_type) :: atom_data
    type(dis_manifold_type) :: dis_manifold
    type(kmesh_info_type) :: kmesh_info
    type(kmesh_input_type) :: kmesh_input
    type(kpoint_path_type) :: kpoint_path
    type(print_output_type) :: print_output
    type(settings_type) :: settings
    !! container for input file (.win) data and options set via library interface
    type(timer_list_type) :: timer
    type(w90_physical_constants_type) :: physics
    type(w90_system_type) :: w90_system
    type(wannier_data_type) :: wannier_data
    type(ws_region_type) :: ws_region

    complex(kind=dp), pointer :: u_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: u_opt(:, :, :) => null()
    real(kind=dp), pointer :: eigval(:, :) => null()
    !! matrices

    real(kind=dp), allocatable :: fermi_energy_list(:)
    real(kind=dp), allocatable :: kpt_latt(:, :)
    real(kind=dp) :: real_lattice(3, 3)

    integer, pointer :: dist_kpoints(:) => null()
    !! dist_kpoints(i) = rank operating on k-point i

    integer, allocatable :: exclude_bands(:)
    integer :: mp_grid(3)
    integer :: num_bands
    integer :: num_kpts
    integer :: num_wann

    logical :: gamma_only
    logical :: have_disentangled = .false.

    character(len=128) :: seedname
  end type lib_global_type

  ! datatype encapsulating types of use in wannier90 exclusively
  type lib_w90_type
    type(dis_control_type) :: dis_control
    type(dis_spheres_type) :: dis_spheres
    type(ham_logical_type) :: ham_logical
    type(output_file_type) :: output_file
    type(proj_input_type) :: proj
    type(real_space_ham_type) :: real_space_ham
    type(select_projection_type) :: select_proj
    type(sitesym_type) :: sitesym
    type(w90_calculation_type) :: w90_calculation ! separate this? ... maybe yes (JJ)
    type(wann_control_type) :: wann_control
    type(wann_omega_type) :: omega
    type(wann_omega_type) :: wann_omega

    complex(kind=dp), allocatable :: ham_k(:, :, :)
    complex(kind=dp), allocatable :: ham_r(:, :, :)
    complex(kind=dp), pointer :: a_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: m_matrix_local(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_matrix(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_orig(:, :, :, :) => null()  !m_matrix_orig_local

    real(kind=dp), allocatable :: wannier_centres_translated(:, :)

    integer, allocatable :: irvec(:, :)
    integer, allocatable :: ndegen(:)
    integer, allocatable :: shift_vec(:, :)

    integer :: nrpts
    integer :: num_proj = 0
    integer :: optimisation = 3
    integer :: rpt_origin
    logical :: lhasproj = .false.

    ! plot
    type(band_plot_type) :: band_plot
    type(wannier_plot_type) :: wann_plot
    type(fermi_surface_plot_type) :: fermi_surface_data

    ! transport
    type(transport_type) :: tran
    type(wvfn_read_type) :: wvfn_read

    ! symmetry-adapted Wannier functions
    logical :: calc_only_A = .false.
    logical :: lsitesymmetry = .false.
    logical :: use_bloch_phases = .false.
  end type lib_w90_type

  private :: prterr
  public :: create_kmesh
  public :: get_fortran_stderr
  public :: get_fortran_stdout
  public :: input_reader
  public :: input_setopt
  public :: overlaps
  public :: plot_files
  public :: print_times
  public :: projovlp
  public :: read_chkpt
  public :: transport
  public :: wannierise
  public :: write_chkpt
  public :: write_kmesh

  public :: set_option ! interface for setting vars using library
  interface set_option
    module procedure set_option_logical
    !module procedure set_option_b1d
    module procedure set_option_text
    module procedure set_option_int
    module procedure set_option_i1d
    module procedure set_option_i2d
    module procedure set_option_r1d
    module procedure set_option_r2d
    module procedure set_option_real
  end interface set_option

contains

  subroutine get_fortran_stdout(istdout)
    use iso_fortran_env, only: output_unit
    implicit none
    integer, intent(out) :: istdout

    istdout = output_unit
  end subroutine get_fortran_stdout

  subroutine get_fortran_stderr(istdout)
    use iso_fortran_env, only: error_unit
    implicit none
    integer, intent(out) :: istdout

    istdout = error_unit
  end subroutine get_fortran_stderr

  subroutine get_fortran_file(output, name)
    implicit none
    integer, intent(out) :: output
    character(len=*), intent(in) :: name

    open (newunit=output, file=name, form='formatted', status='unknown')
  end subroutine get_fortran_file

  subroutine write_chkpt(helper, wan90, label, seedname, istdout, istderr, ierr, comm)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write_chkpt
    use w90_comms, only: comms_reduce, mpirank, w90_comm_type
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal

    implicit none

    ! arguments
    character(len=*), intent(in) :: label ! e.g. 'postdis' or 'postwann' after disentanglement, wannierisation
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(inout) :: ierr
    type(lib_global_type), target, intent(in) :: helper
    type(lib_w90_type), intent(in) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    complex(kind=dp), allocatable :: u(:, :, :), uopt(:, :, :), m(:, :, :, :)
    integer, allocatable :: global_k(:)
    integer, pointer :: nw, nb, nk, nn
    integer :: rank, nkl, ikg, ikl
    type(w90_error_type), allocatable :: error

    ierr = 0

    rank = mpirank(comm)

    nb => helper%num_bands
    nk => helper%num_kpts
    nn => helper%kmesh_info%nntot
    nw => helper%num_wann

    if (.not. associated(helper%u_opt)) then
      call set_error_fatal(error, 'u_opt not set for write_chkpt call', comm)
    else if (.not. associated(helper%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for write_chkpt call', comm)
    else if (.not. associated(wan90%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for write_chkpt call', comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    nkl = count(helper%dist_kpoints == rank)
    allocate (global_k(nkl))
    global_k = huge(1); ikl = 1
    do ikg = 1, nk
      if (helper%dist_kpoints(ikg) == rank) then
        global_k(ikl) = ikg
        ikl = ikl + 1
      endif
    enddo

    ! allocating and partially assigning the full matrix on all ranks and reducing is a terrible idea at scale
    ! alternatively, allocate on root and use point-to-point
    ! or, if required only for checkpoint file writing, then use mpi-io (but needs to be ordered io, alas)
    ! or, even better, use parallel hdf5
    ! fixme.  JJ Nov 22
    allocate (u(nw, nw, nk)) ! all kpts
    allocate (uopt(nb, nw, nk)) ! all kpts
    allocate (m(nw, nw, nn, nk)) ! all kpts
    u(:, :, :) = 0.d0
    uopt(:, :, :) = 0.d0
    m(:, :, :, :) = 0.d0
    do ikl = 1, nkl
      ikg = global_k(ikl)
      u(:, :, ikg) = helper%u_matrix(:, :, ikl)
      uopt(:, :, ikg) = helper%u_opt(:, :, ikl)
      m(:, :, :, ikg) = wan90%m_matrix_local(:, :, :, ikl)
    enddo
    call comms_reduce(u(1, 1, 1), nw*nw*nk, 'SUM', error, comm)
    call comms_reduce(uopt(1, 1, 1), nb*nw*nk, 'SUM', error, comm)
    call comms_reduce(m(1, 1, 1, 1), nw*nw*nn*nk, 'SUM', error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    if (rank == 0) then
      call w90_wannier90_readwrite_write_chkpt(label, helper%exclude_bands, helper%wannier_data, &
                                               helper%kmesh_info, helper%kpt_latt, nk, &
                                               helper%dis_manifold, nb, nw, u, uopt, m, &
                                               helper%mp_grid, helper%real_lattice, &
                                               wan90%omega%invariant, helper%have_disentangled, &
                                               istdout, seedname)
    endif
    deallocate (u)
    deallocate (uopt)
    deallocate (m)
  end subroutine write_chkpt

  subroutine read_chkpt(helper, wan90, checkpoint, seedname, istdout, istderr, ierr, comm)
    use w90_comms, only: w90_comm_type, mpirank
    use w90_error_base, only: w90_error_type
    use w90_readwrite, only: w90_readwrite_read_chkpt, w90_readwrite_chkpt_dist

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    character(len=*), intent(out) :: checkpoint
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), target, intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    complex(kind=dp), allocatable :: m(:, :, :, :)
    integer, pointer :: nw, nb, nk, nn
    integer :: rank, nexclude = 0
    logical :: ispostw90 = .false. ! ispostw90 is used to print a different error message in case the chk file is missing (did you run w90 first?)
    type(w90_error_type), allocatable :: error

    ierr = 0

    rank = mpirank(comm)

    nb => helper%num_bands
    nk => helper%num_kpts
    nn => helper%kmesh_info%nntot
    nw => helper%num_wann

    ! allocating and partially assigning the full matrix on all ranks and reducing is a terrible idea at scale
    ! alternatively, allocate on root and use point-to-point
    ! or, if required only for checkpoint file writing, then use mpi-io (but needs to be ordered io, alas)
    ! or, even better, use parallel hdf5
    ! fixme.  JJ Nov 22
    allocate (m(nw, nw, nn, nk)) ! all kpts

    if (rank == 0) then
      if (allocated(helper%exclude_bands)) nexclude = size(helper%exclude_bands)

      call w90_readwrite_read_chkpt(helper%dis_manifold, helper%exclude_bands, helper%kmesh_info, helper%kpt_latt, &
                                    helper%wannier_data, m, helper%u_matrix, helper%u_opt, &
                                    helper%real_lattice, wan90%omega%invariant, helper%mp_grid, nb, &
                                    nexclude, nk, nw, checkpoint, &
                                    helper%have_disentangled, ispostw90, seedname, istdout, error, comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, comm)
        return
      endif
    endif

    ! scatter from m_matrix to m_matrix_local
    ! normally achieved in overlap_read
    call w90_readwrite_chkpt_dist(helper%dis_manifold, helper%wannier_data, helper%u_matrix, &
                                  helper%u_opt, m, wan90%m_matrix_local, wan90%omega%invariant, &
                                  nb, nk, nw, nn, checkpoint, helper%have_disentangled, &
                                  helper%dist_kpoints, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
    deallocate (m)
  end subroutine read_chkpt

  subroutine input_setopt(helper, wan90, seedname, istdout, istderr, ierr, comm)
    use w90_readwrite, only: w90_readwrite_uppercase, w90_readwrite_read_final_alloc
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    !use w90_error, only: set_error_input, set_error_fatal
    use w90_comms, only: w90_comm_type, mpirank, comms_sync_err

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    logical :: cp_pp, disentanglement

    ierr = 0

    call w90_wannier90_readwrite_read(helper%settings, helper%atom_data, wan90%band_plot, wan90%dis_control, &
                                      wan90%dis_spheres, helper%dis_manifold, &
                                      helper%exclude_bands, helper%fermi_energy_list, &
                                      wan90%fermi_surface_data, helper%kmesh_input, &
                                      helper%kmesh_info, helper%kpt_latt, wan90%output_file, &
                                      wan90%wvfn_read, wan90%wann_control, wan90%proj, &
                                      wan90%real_space_ham, wan90%select_proj, &
                                      helper%kpoint_path, helper%w90_system, wan90%tran, &
                                      helper%print_output, wan90%wann_plot, io_params, &
                                      helper%ws_region, wan90%w90_calculation, &
                                      helper%real_lattice, helper%physics%bohr, &
                                      wan90%sitesym%symmetrize_eps, helper%mp_grid, &
                                      helper%num_bands, helper%num_kpts, wan90%num_proj, &
                                      helper%num_wann, wan90%optimisation, wan90%calc_only_A, &
                                      cp_pp, helper%gamma_only, wan90%lhasproj, &
                                      wan90%lsitesymmetry, wan90%use_bloch_phases, seedname, &
                                      istdout, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    else
      ! For aesthetic purposes, convert some things to uppercase
      call w90_readwrite_uppercase(helper%settings, helper%atom_data, helper%kpoint_path, &
                                   helper%print_output%length_unit)

      disentanglement = (helper%num_bands > helper%num_wann)

      call w90_readwrite_read_final_alloc(disentanglement, helper%dis_manifold, &
                                          helper%wannier_data, helper%num_wann, &
                                          helper%num_bands, helper%num_kpts, error, comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, comm)
        return
      endif
    endif
    helper%seedname = seedname ! maybe not keep this separate from "blob"? JJ

    if (mpirank(comm) /= 0) helper%print_output%iprint = 0 ! supress printing non-rank-0
  end subroutine input_setopt

  subroutine input_reader(helper, wan90, seedname, istdout, istderr, ierr, comm)
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_uppercase, &
      w90_readwrite_clean_infile, w90_readwrite_read_final_alloc
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_input, set_error_fatal
    use w90_comms, only: w90_comm_type, mpirank, comms_sync_err

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error
    !logical :: cp_pp ! ? when used?

    ierr = 0

    call w90_readwrite_in_file(helper%settings, seedname, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    call input_setopt(helper, wan90, seedname, istdout, istderr, ierr, comm)
    if (ierr /= 0) then
      return
    endif

    ! test mpi error handling using "unlucky" input token
    ! this machinery used to sit in w90_wannier90_readwrite_dist
    ! but that routine is obsolete if input file is read on all ranks
    ! fixme, this should be moved to wannier90 main routine (definately doesn't belong here)
    if (helper%print_output%timing_level < 0 .and. mpirank(comm) == abs(helper%print_output%timing_level)) then
      call set_error_input(error, 'received unlucky_rank', comm)
    else
      call comms_sync_err(comm, error, 0) ! this is necessary since non-root may never enter an mpi collective if root has exited here
    endif
    if (allocated(error)) then ! applies (is t) for all ranks now
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
    !!!!! end unlucky code

    call w90_readwrite_clean_infile(helper%settings, istdout, seedname, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
  end subroutine input_reader

  subroutine create_kmesh(helper, istdout, istderr, ierr, comm)
    use w90_kmesh, only: kmesh_get
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    ierr = 0
    call kmesh_get(helper%kmesh_input, helper%kmesh_info, helper%print_output, helper%kpt_latt, &
                   helper%real_lattice, helper%num_kpts, helper%gamma_only, istdout, helper%timer, &
                   error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
  end subroutine create_kmesh

  subroutine write_kmesh(helper, wan90, seedname, istdout, istderr, ierr, comm)
    use w90_kmesh, only: kmesh_get, kmesh_write
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank, comms_sync_err

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname ! needed for nnkp filename
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    type(w90_error_type), allocatable :: error
    logical :: calc_only_A = .false.! what does this do?

    ierr = 0

    if (mpirank(comm) == 0) then
      call kmesh_write(helper%exclude_bands, helper%kmesh_info, wan90%proj, helper%print_output, &
                       helper%kpt_latt, helper%real_lattice, helper%num_kpts, wan90%num_proj, &
                       calc_only_A, helper%w90_system%spinors, seedname, helper%timer)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, comm)
        return
      endif
    endif
    call comms_sync_err(comm, error, 0) ! this is necessary since non-root may never enter an mpi collective if root has exited here
  end subroutine write_kmesh

  subroutine overlaps(helper, wan90, istdout, istderr, ierr, comm)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_comms, only: w90_comm_type
    use w90_overlap, only: overlap_read

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    logical :: cp_pp = .false.
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. associated(helper%dist_kpoints)) then
      call set_error_fatal(error, 'dist_kpoints not set for overlap call', comm)
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    ! fixme jj checkit
    if (helper%num_bands > helper%num_wann) then ! disentanglement case
      if ((.not. associated(wan90%a_matrix)) .or. (.not. associated(wan90%m_orig))) then
        write (istderr, *) 'Matrices not set for overlap call (disentanglement case)'
        ierr = 1
        return
      endif
      call overlap_read(helper%kmesh_info, wan90%select_proj, wan90%sitesym, wan90%a_matrix, &
                        wan90%m_orig, helper%num_bands, helper%num_kpts, wan90%num_proj, &
                        helper%num_wann, helper%print_output%timing_level, cp_pp, &
                        helper%gamma_only, wan90%lsitesymmetry, wan90%use_bloch_phases, &
                        helper%seedname, istdout, helper%timer, helper%dist_kpoints, error, comm)
    else
      if ((.not. associated(helper%u_matrix)) .or. (.not. associated(wan90%m_matrix_local))) then
        write (istderr, *) 'Matrices not set for overlap call'
        ierr = 1
        return
      endif
      call overlap_read(helper%kmesh_info, wan90%select_proj, wan90%sitesym, helper%u_matrix, &
                        wan90%m_matrix_local, helper%num_bands, helper%num_kpts, wan90%num_proj, &
                        helper%num_wann, helper%print_output%timing_level, cp_pp, &
                        helper%gamma_only, wan90%lsitesymmetry, wan90%use_bloch_phases, &
                        helper%seedname, istdout, helper%timer, helper%dist_kpoints, error, comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
  end subroutine overlaps

  subroutine disentangle(helper, wan90, istdout, istderr, ierr, comm)
    use w90_disentangle, only: dis_main, setup_m_loc
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_comms, only: w90_comm_type, mpirank

    implicit none

    ! arguments
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error
    integer :: optimisation = 3

    ierr = 0

    if (.not. associated(wan90%m_orig)) then  ! m_matrix_orig_local (nband*nwann for disentangle)
      call set_error_fatal(error, 'm_orig not set for disentangle call', comm)
    else if (.not. associated(wan90%m_matrix_local)) then ! (nband*nwann*nknode for wannierise)
      call set_error_fatal(error, 'm_matrix_local not set for disentangle call', comm)
    else if (.not. associated(wan90%a_matrix)) then
      call set_error_fatal(error, 'a_matrix not set for disentangle call', comm)
    else if (.not. associated(helper%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for disentangle call', comm)
    else if (.not. associated(helper%u_opt)) then
      call set_error_fatal(error, 'u_opt not set for disentangle call', comm)
    else if (.not. associated(helper%dist_kpoints)) then
      call set_error_fatal(error, 'kpt decomp not set for disentangle call', comm)
    else if (.not. associated(helper%eigval)) then
      call set_error_fatal(error, 'eigval not set for disentangle call', comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    call dis_main(wan90%dis_control, wan90%dis_spheres, helper%dis_manifold, helper%kmesh_info, &
                  helper%kpt_latt, wan90%sitesym, helper%print_output, wan90%a_matrix, &
                  wan90%m_orig, helper%u_matrix, helper%u_opt, helper%eigval, &
                  helper%real_lattice, wan90%omega%invariant, helper%num_bands, helper%num_kpts, &
                  helper%num_wann, helper%gamma_only, wan90%lsitesymmetry, istdout, helper%timer, &
                  helper%dist_kpoints, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    ! copy to m_matrix_local from m_matrix_orig_local (aka m_orig)
    call setup_m_loc(helper%kmesh_info, helper%print_output, wan90%m_matrix_local, wan90%m_orig, &
                     helper%u_matrix, helper%num_bands, helper%num_kpts, helper%num_wann, &
                     optimisation, helper%timer, helper%dist_kpoints, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    helper%have_disentangled = .true.
  end subroutine disentangle

  subroutine projovlp(helper, wan90, istdout, istderr, ierr, comm) !fixme(jj) rename more sensibly
    use w90_comms, only: w90_comm_type
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_overlap, only: overlap_project

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. associated(wan90%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for disentangle call', comm)
    else if (.not. associated(helper%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for disentangle call', comm)
    else if (.not. associated(helper%dist_kpoints)) then
      call set_error_fatal(error, 'dist_kpoints not set for disentangle call', comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    if (.not. helper%have_disentangled) then
      ! fixme (jj) there is also a gamma only specialisation of this
      call overlap_project(wan90%sitesym, wan90%m_matrix_local, helper%u_matrix, &
                           helper%kmesh_info%nnlist, helper%kmesh_info%nntot, &
                           helper%num_wann, helper%num_kpts, helper%num_wann, &
                           helper%print_output%timing_level, wan90%lsitesymmetry, &
                           istdout, helper%timer, helper%dist_kpoints, error, comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, comm)
        return
      endif
    endif
  end subroutine projovlp

  subroutine wannierise(helper, wan90, istdout, istderr, ierr, comm)
    use w90_comms, only: w90_comm_type
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_wannierise, only: wann_main, wann_main_gamma

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. associated(wan90%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for disentangle call', comm)
    else if (.not. associated(helper%u_opt)) then
      call set_error_fatal(error, 'u_opt not set for disentangle call', comm)
    else if (.not. associated(helper%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for disentangle call', comm)
    else if (.not. associated(helper%dist_kpoints)) then
      call set_error_fatal(error, 'dist_kpoints not set for disentangle call', comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    if (helper%gamma_only) then
      call wann_main_gamma(helper%kmesh_info, wan90%wann_control, wan90%omega, &
                           helper%print_output, helper%wannier_data, wan90%m_matrix_local, &
                           helper%u_matrix, helper%real_lattice, helper%num_kpts, helper%num_wann, &
                           istdout, helper%timer, error, comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, comm)
        return
      endif
    else
      call wann_main(wan90%ham_logical, helper%kmesh_info, helper%kpt_latt, wan90%wann_control, &
                     wan90%omega, wan90%sitesym, helper%print_output, helper%wannier_data, &
                     helper%ws_region, wan90%w90_calculation, wan90%ham_k, wan90%ham_r, &
                     wan90%m_matrix_local, helper%u_matrix, helper%real_lattice, &
                     wan90%wannier_centres_translated, wan90%irvec, helper%mp_grid, wan90%ndegen, &
                     wan90%nrpts, helper%num_kpts, wan90%num_proj, helper%num_wann, &
                     wan90%optimisation, wan90%rpt_origin, wan90%band_plot%mode, wan90%tran%mode, &
                     wan90%lsitesymmetry, istdout, helper%timer, helper%dist_kpoints, error, comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
  end subroutine wannierise

  subroutine plot_files(helper, wan90, istdout, istderr, ierr, comm)
    use w90_comms, only: w90_comm_type
    use w90_error_base, only: w90_error_type
    use w90_plot, only: plot_main

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper ! inout due to ham_logical -- JJ: eugh, that's nasty, can we change it?
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    ierr = 0

    ! fixme(jj) what are our preconditions?

    call plot_main(helper%atom_data, wan90%band_plot, helper%dis_manifold, helper%fermi_energy_list, &
                   wan90%fermi_surface_data, wan90%ham_logical, helper%kmesh_info, helper%kpt_latt, &
                   wan90%output_file, wan90%wvfn_read, wan90%real_space_ham, helper%kpoint_path, &
                   helper%print_output, helper%wannier_data, wan90%wann_plot, helper%ws_region, &
                   wan90%w90_calculation, wan90%ham_k, wan90%ham_r, wan90%m_matrix_local, helper%u_matrix, &
                   helper%u_opt, helper%eigval, helper%real_lattice, wan90%wannier_centres_translated, &
                   helper%physics%bohr, wan90%irvec, helper%mp_grid, wan90%ndegen, wan90%shift_vec, wan90%nrpts, &
                   helper%num_bands, helper%num_kpts, helper%num_wann, wan90%rpt_origin, &
                   wan90%tran%mode, helper%have_disentangled, wan90%lsitesymmetry, helper%w90_system, &
                   helper%seedname, istdout, helper%timer, helper%dist_kpoints, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
  end subroutine plot_files

  subroutine transport(helper, wan90, istdout, istderr, ierr, comm)
    use w90_comms, only: w90_comm_type, mpirank
    use w90_error_base, only: w90_error_type
    use w90_transport, only: tran_main

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper ! because of ham_logical
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    ierr = 0

    ! fixme(jj) what are our preconditions?

    ! fixme(JJ) is tran_main pllel at all?
    if (mpirank(comm) == 0) then
      call tran_main(helper%atom_data, helper%dis_manifold, helper%fermi_energy_list, &
                     wan90%ham_logical, helper%kpt_latt, wan90%output_file, wan90%real_space_ham, &
                     wan90%tran, helper%print_output, helper%wannier_data, helper%ws_region, &
                     wan90%w90_calculation, wan90%ham_k, wan90%ham_r, helper%u_matrix, helper%u_opt, &
                     helper%eigval, helper%real_lattice, wan90%wannier_centres_translated, wan90%irvec, &
                     helper%mp_grid, wan90%ndegen, wan90%shift_vec, wan90%nrpts, helper%num_bands, &
                     helper%num_kpts, helper%num_wann, wan90%rpt_origin, wan90%band_plot%mode, &
                     helper%have_disentangled, wan90%lsitesymmetry, helper%seedname, istdout, &
                     helper%timer, error, comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, comm)
        return
      endif
    endif
  end subroutine transport

  subroutine print_times(helper, istdout)
    use w90_io, only: io_print_timings
    implicit none
    type(lib_global_type), intent(in) :: helper
    integer, intent(in) :: istdout

    if (helper%print_output%iprint > 0) call io_print_timings(helper%timer, istdout)
  end subroutine print_times

  subroutine set_a_matrix(helper, a_matrix)
    implicit none
    type(lib_w90_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: a_matrix(:, :, :)

    helper%a_matrix => a_matrix
  end subroutine set_a_matrix

  subroutine set_m_matrix(helper, m_matrix)
    implicit none
    type(lib_w90_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: m_matrix(:, :, :, :)

    helper%m_matrix => m_matrix
  end subroutine set_m_matrix

  subroutine set_m_matrix_local(helper, m_matrix_local) ! scattered m-matrix
    implicit none
    type(lib_w90_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: m_matrix_local(:, :, :, :)

    helper%m_matrix_local => m_matrix_local
  end subroutine set_m_matrix_local

  subroutine set_m_orig(helper, m_orig) ! m_matrix_local_orig
    implicit none
    type(lib_w90_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: m_orig(:, :, :, :)

    helper%m_orig => m_orig
  end subroutine set_m_orig

  subroutine set_u_matrix(helper, u_matrix)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: u_matrix(:, :, :)

    helper%u_matrix => u_matrix
  end subroutine set_u_matrix

  subroutine set_u_opt(helper, u_opt)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: u_opt(:, :, :)

    helper%u_opt => u_opt
  end subroutine set_u_opt

  subroutine set_eigval(helper, eigval)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    real(kind=dp), intent(in), target :: eigval(:, :)

    helper%eigval => eigval
    ! if not already initialised, set disentanglement window to limits of spectrum
    if (helper%dis_manifold%win_min == -huge(0.0_dp)) helper%dis_manifold%win_min = minval(helper%eigval)
    if (helper%dis_manifold%win_max == huge(0.0_dp)) helper%dis_manifold%win_max = maxval(helper%eigval)
  end subroutine set_eigval

  subroutine set_kpoint_distribution(helper, dist, istdout, istderr, ierr, comm)
    use w90_comms, only: w90_comm_type
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal

    implicit none

    ! arguments
    integer, intent(in) :: istderr, istdout
    integer, intent(inout), target :: dist(:)
    integer, intent(out) :: ierr
    type(lib_global_type), intent(inout) :: helper
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (size(dist) < 1) call set_error_fatal(error, 'Error in k-point distribution, mpisize < 1', comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    helper%dist_kpoints => dist
  end subroutine set_kpoint_distribution

  subroutine set_option_text(helper, keyword, text)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    character(*), intent(in) :: text
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%txtdata = text
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_text

  subroutine set_option_logical(helper, keyword, bool)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    logical, intent(in) :: bool
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%ldata = bool
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_logical

  subroutine set_option_i1d(helper, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    integer, intent(in) :: arr(:)
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%i1d = arr ! this causes an automatic allocation
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_i1d

  subroutine set_option_i2d(helper, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    integer, intent(in) :: arr(:, :)
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%i2d = arr
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_i2d

  subroutine set_option_int(helper, keyword, ival)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    integer, intent(in) :: ival
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%idata = ival
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_int

  subroutine set_option_r1d(helper, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    real(kind=dp), intent(in) :: arr(:)
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%r1d = arr
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_r1d

  subroutine set_option_r2d(helper, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    real(kind=dp), intent(in) :: arr(:, :)
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%r2d = arr
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_r2d

  subroutine set_option_real(helper, keyword, rval)
    use w90_readwrite, only: init_settings, expand_settings
    implicit none
    character(*), intent(in) :: keyword
    real(kind=dp), intent(in) :: rval
    type(lib_global_type), intent(inout) :: helper
    integer :: i

    if (.not. allocated(helper%settings%entries)) call init_settings(helper%settings)
    i = helper%settings%num_entries + 1
    helper%settings%entries(i)%keyword = keyword
    helper%settings%entries(i)%rdata = rval
    helper%settings%num_entries = i + 1
    if (helper%settings%num_entries == helper%settings%num_entries_max) call expand_settings(helper%settings)
  endsubroutine set_option_real

  subroutine prterr(error, ie, istdout, istderr, comm)
    use w90_comms, only: comms_no_sync_send, comms_no_sync_recv, w90_comm_type, mpirank, mpisize
    use w90_error_base, only: code_remote, w90_error_type

    ! arguments
    integer, intent(inout) :: ie! global error value to be returned
    integer, intent(in) :: istderr, istdout
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(in) :: error

    ! local
    type(w90_error_type), allocatable :: le ! unchecked error state for calls made in this routine
    integer :: je ! error value on remote ranks
    integer :: j ! rank index
    integer :: failrank ! lowest rank reporting an error
    character(len=128) :: mesg ! only print 128 chars of error

    ie = 0
    mesg = 'not set'

    if (mpirank(comm) == 0) then
      ! fixme, report all failing ranks
      do j = mpisize(comm) - 1, 1, -1
        call comms_no_sync_recv(je, 1, j, le, comm)

        if (je /= code_remote .and. je /= 0) then
          failrank = j
          ie = je
          call comms_no_sync_recv(mesg, 128, j, le, comm)
        endif
      enddo
      ! if the error is on rank0
      if (error%code /= code_remote .and. error%code /= 0) then
        failrank = 0
        ie = error%code
        mesg = error%message
      endif

      !if (ie == 0) write (istderr, *) "logic error" ! to arrive here requires this

      write (istderr, *) 'Exiting.......'
      write (istderr, '(1x,a)') trim(mesg)
      write (istderr, '(1x,a,i0,a)') '(rank: ', failrank, ')'
      write (istdout, '(1x,a)') ' error encountered; check error .werr error log'

    else ! non 0 ranks
      je = error%code
      call comms_no_sync_send(je, 1, 0, le, comm)
      if (je /= code_remote .and. je /= 0) then
        mesg = error%message
        call comms_no_sync_send(mesg, 128, 0, le, comm)
      endif
    endif
  end subroutine prterr

  ! this routine needs revising/moving fixme(jj)
  ! this routine assigns to the array passed as an argument, not to the internal module pointer
  ! this array must subsequently be the subject of associating the internal module pointer
  subroutine read_eigvals(w90main, w90dat, eigval, seedname, istdout, istderr, ierr, comm)
    use w90_comms, only: w90_comm_type
    use w90_error, only: w90_error_type, set_error_fatal
    use w90_readwrite, only: w90_readwrite_read_eigvals

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    real(kind=dp), intent(inout) :: eigval(:, :)
    type(lib_global_type), intent(inout) :: w90main
    type(lib_w90_type), intent(in) :: w90dat
    type(w90_comm_type), intent(in) :: comm

    ! local vars
    type(w90_error_type), allocatable :: error
    logical :: eig_found
    integer :: ierr

    ierr = 0

    if (size(eigval, 1) /= w90main%num_bands) then
      call set_error_fatal(error, 'eigval not dimensioned correctly (num_bands,num_kpts) in read_eigvals', comm)
      call prterr(error, ierr, istdout, istderr, comm)
      return
    elseif (size(eigval, 2) /= w90main%num_kpts) then
      call set_error_fatal(error, 'eigval not dimensioned correctly (num_bands,num_kpts) in read_eigvals', comm)
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif

    call w90_readwrite_read_eigvals(w90main%settings, eig_found, eigval, w90main%num_bands, w90main%num_kpts, &
                                    istdout, seedname, error, comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, comm)
      return
    else if (.not. eig_found) then
      call set_error_fatal(error, 'failed to read eigenvalues file in read_eigvals call', comm)
      call prterr(error, ierr, istdout, istderr, comm)
      return
    endif
  end subroutine read_eigvals
end module w90_helper_types
