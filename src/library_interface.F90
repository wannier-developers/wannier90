
module w90_library

  ! as with fortran routines like allocate, the status variable indicates an error if non-zero
  ! positive is an error, negative is a warning (such as non-convergence) which is recoverable.

  ! note on the naming of output and error-output streams
  !   the specific names stdout and stderr cannot be used because
  !   the python wrapper uses an intermediate c representation with named arguments
  !   in this wrapper stdout/err are already defined as structs inconsistent
  !   which are inconsistent with the integer unit numbers used here.
  !   (The slightly 77 style naming here has no particular significance.)

  use w90_constants
  use w90_types
  use w90_wannier90_types
  use w90_comms, only: w90_comm_type
  use iso_c_binding

  implicit none

  ! datatype encapsulating types used by wannier90
  type lib_common_type
    character(len=128) :: seedname

    !! matrices
    complex(kind=dp), allocatable :: ham_k(:, :, :)
    complex(kind=dp), allocatable :: ham_r(:, :, :)
    complex(kind=dp), pointer :: a_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: m_matrix_local(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_matrix(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_orig(:, :, :, :) => null()  !m_matrix_orig_local
    complex(kind=dp), pointer :: u_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: u_opt(:, :, :) => null()
    real(kind=dp), pointer :: eigval(:, :) => null()

    integer, allocatable :: dist_kpoints(:) ! dist_kpoints(i) = rank operating on k-point i
    integer, allocatable :: exclude_bands(:)

    integer, allocatable :: irvec(:, :)
    integer, allocatable :: ndegen(:)
    integer, allocatable :: shift_vec(:, :)

    integer :: mp_grid(3)
    integer :: nrpts
    integer :: num_bands
    integer :: num_kpts
    integer :: num_proj = 0
    integer :: num_wann
    integer :: optimisation = 3
    integer :: rpt_origin

    logical :: calc_only_A = .false. ! used only in kmesh
    logical :: gamma_only = .false.
    logical :: have_disentangled = .false.
    logical :: lhasproj = .false.
    logical :: lsitesymmetry = .false.
    logical :: use_bloch_phases = .false.
    logical :: setup_complete = .false.

    real(kind=dp), allocatable :: fermi_energy_list(:)
    real(kind=dp), allocatable :: kpt_latt(:, :)
    real(kind=dp), allocatable :: wannier_centres_translated(:, :)
    real(kind=dp) :: real_lattice(3, 3)

    type(atom_data_type) :: atom_data
    type(band_plot_type) :: band_plot
    type(dis_control_type) :: dis_control
    type(dis_manifold_type) :: dis_manifold
    type(dis_spheres_type) :: dis_spheres
    type(fermi_surface_plot_type) :: fermi_surface_data
    type(ham_logical_type) :: ham_logical
    type(kmesh_info_type) :: kmesh_info
    type(kmesh_input_type) :: kmesh_input
    type(kpoint_path_type) :: kpoint_path
    type(output_file_type) :: output_file
    type(print_output_type) :: print_output
    type(proj_type), allocatable :: proj(:), proj_input(:)
    type(real_space_ham_type) :: real_space_ham
    type(select_projection_type) :: select_proj
    type(settings_type) :: settings !! container for input file (.win) data and options set via library interface
    type(sitesym_type) :: sitesym ! symmetry-adapted Wannier functions
    type(timer_list_type) :: timer
    type(transport_type) :: tran
    type(w90_calculation_type) :: w90_calculation
    type(w90_comm_type) :: comm
    type(w90_physical_constants_type) :: physics
    type(w90_system_type) :: w90_system
    type(wann_control_type) :: wann_control
    type(wannier_data_type) :: wannier_data
    type(wannier_plot_type) :: wann_plot
    type(wann_omega_type) :: omega
    type(wann_omega_type) :: wann_omega
    type(ws_region_type) :: ws_region
    type(wvfn_read_type) :: wvfn_read
  end type lib_common_type

  public :: prterr
  public :: create_kmesh
  public :: disentangle
  public :: get_centres
  public :: get_fortran_file
  public :: get_fortran_stderr
  public :: get_fortran_stdout
  public :: get_proj
  public :: get_spreads
  public :: input_reader
  public :: input_reader_special
  public :: input_setopt
  public :: overlaps !to standalone driver (just reads a file)
  public :: plot_files
  public :: print_times
  public :: projovlp
  public :: read_chkpt
  public :: read_eigvals
  !public :: set_a_matrix
  public :: set_constant_bohr_to_ang
  public :: set_eigval
  public :: set_kpoint_distribution
  public :: set_m_orig
  public :: set_option ! interface for setting vars using library
  public :: set_parallel_comms
  public :: set_u_matrix
  public :: set_u_opt
  public :: transport
  public :: wannierise
  public :: write_chkpt
  public :: write_kmesh

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

  subroutine write_chkpt(common_data, label, istdout, istderr, ierr)
    use w90_comms, only: comms_reduce, mpirank
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_alloc, set_error_dealloc, set_error_fatal
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write_chkpt

    implicit none

    ! arguments
    character(len=*), intent(in) :: label ! e.g. 'postdis' or 'postwann' after disentanglement, wannierisation
    integer, intent(in) :: istdout, istderr
    integer, intent(inout) :: ierr
    type(lib_common_type), target, intent(in) :: common_data

    ! local variables
    complex(kind=dp), allocatable :: u(:, :, :), uopt(:, :, :), m(:, :, :, :)
    integer, allocatable :: global_k(:)
    integer, pointer :: nw, nb, nk, nn
    integer :: rank, nkrank, ikg, ikl, istat
    type(w90_error_type), allocatable :: error

    ierr = 0
    rank = mpirank(common_data%comm)

    nb => common_data%num_bands
    nk => common_data%num_kpts
    nn => common_data%kmesh_info%nntot
    nw => common_data%num_wann

    if (.not. associated(common_data%u_opt)) then
      call set_error_fatal(error, 'u_opt not set for write_chkpt call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for write_chkpt call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    else if (.not. associated(common_data%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for write_chkpt call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif

    nkrank = count(common_data%dist_kpoints == rank)
    allocate (global_k(nkrank), stat=istat)
    if (istat /= 0) then
      call set_error_alloc(error, 'Error allocating global_k in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
    global_k = huge(1); ikl = 1
    do ikg = 1, nk
      if (rank == common_data%dist_kpoints(ikg)) then
        global_k(ikl) = ikg
        ikl = ikl + 1
      endif
    enddo

    ! allocating and partially assigning the full matrix on all ranks and reducing is a terrible idea
    ! alternatively, allocate on root and use point-to-point
    ! or, if required only for checkpoint file writing, then use mpi-io (but needs to be ordered io, alas)
    ! or, even better, use parallel hdf5. JJ Nov 22
    allocate (u(nw, nw, nk), stat=istat) ! all kpts
    if (istat /= 0) then
      call set_error_alloc(error, 'Error allocating u in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
    allocate (uopt(nb, nw, nk), stat=istat) ! all kpts
    if (istat /= 0) then
      call set_error_alloc(error, 'Error allocating uopt in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
    allocate (m(nw, nw, nn, nk), stat=istat) ! all kpts
    if (istat /= 0) then
      call set_error_alloc(error, 'Error allocating m in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif

    u(:, :, :) = 0.d0
    uopt(:, :, :) = 0.d0
    m(:, :, :, :) = 0.d0

    do ikl = 1, nkrank
      ikg = global_k(ikl)
      u(:, :, ikg) = common_data%u_matrix(:, :, ikl)
      uopt(:, :, ikg) = common_data%u_opt(:, :, ikl)
      m(:, :, :, ikg) = common_data%m_matrix_local(1:nw, 1:nw, :, ikl)
    enddo

    call comms_reduce(u(1, 1, 1), nw*nw*nk, 'SUM', error, common_data%comm)
    call comms_reduce(uopt(1, 1, 1), nb*nw*nk, 'SUM', error, common_data%comm)
    call comms_reduce(m(1, 1, 1, 1), nw*nw*nn*nk, 'SUM', error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (rank == 0) then
      call w90_wannier90_readwrite_write_chkpt(label, common_data%exclude_bands, &
                                               common_data%wannier_data, common_data%kmesh_info, &
                                               common_data%kpt_latt, nk, common_data%dis_manifold, &
                                               nb, nw, u, uopt, m, common_data%mp_grid, &
                                               common_data%real_lattice, &
                                               common_data%omega%invariant, &
                                               common_data%have_disentangled, &
                                               common_data%print_output%iprint, istdout, &
                                               common_data%seedname)
    endif

    deallocate (u, stat=istat)
    if (istat /= 0) then
      call set_error_dealloc(error, 'Error deallocating u in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
    deallocate (uopt, stat=istat)
    if (istat /= 0) then
      call set_error_dealloc(error, 'Error deallocating uopt in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
    deallocate (m, stat=istat)
    if (istat /= 0) then
      call set_error_dealloc(error, 'Error deallocating m in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
  end subroutine write_chkpt

  subroutine read_chkpt(common_data, checkpoint, istdout, istderr, ierr)
    use w90_comms, only: mpirank
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_alloc, set_error_dealloc
    use w90_readwrite, only: w90_readwrite_read_chkpt, w90_readwrite_chkpt_dist

    implicit none

    ! arguments
    character(len=*), intent(out) :: checkpoint
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), target, intent(inout) :: common_data

    ! local variables
    complex(kind=dp), allocatable :: m(:, :, :, :)
    integer, pointer :: nw, nb, nk, nn
    integer :: rank, nexclude = 0, istat
    logical :: ispostw90 = .false. ! ispostw90 is used to print a different error message in case the chk file is missing (did you run w90 first?)
    type(w90_error_type), allocatable :: error

    ierr = 0
    rank = mpirank(common_data%comm)

    nb => common_data%num_bands
    nk => common_data%num_kpts
    nn => common_data%kmesh_info%nntot
    nw => common_data%num_wann

    ! allocating and partially assigning the full matrix on all ranks and reducing is a terrible idea at scale
    ! alternatively, allocate on root and use point-to-point
    ! or, if required only for checkpoint file writing, then use mpi-io (but needs to be ordered io, alas)
    ! or, even better, use parallel hdf5
    ! fixme, check allocation status of u, uopt?
    allocate (m(nw, nw, nn, nk), stat=istat) ! all kpts
    if (istat /= 0) then
      call set_error_alloc(error, 'Error allocating m in read_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif

    if (rank == 0) then
      if (allocated(common_data%exclude_bands)) nexclude = size(common_data%exclude_bands)

      call w90_readwrite_read_chkpt(common_data%dis_manifold, common_data%exclude_bands, &
                                    common_data%kmesh_info, common_data%kpt_latt, &
                                    common_data%wannier_data, m, common_data%u_matrix, &
                                    common_data%u_opt, common_data%real_lattice, &
                                    common_data%omega%invariant, common_data%mp_grid, nb, &
                                    nexclude, nk, nw, checkpoint, common_data%have_disentangled, &
                                    ispostw90, common_data%seedname, istdout, error, &
                                    common_data%comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif

    ! scatter from m_matrix to m_matrix_local (cf overlap_read)
    call w90_readwrite_chkpt_dist(common_data%dis_manifold, common_data%wannier_data, &
                                  common_data%u_matrix, common_data%u_opt, m, &
                                  common_data%m_matrix_local, common_data%omega%invariant, &
                                  nb, nk, nw, nn, checkpoint, common_data%have_disentangled, &
                                  common_data%dist_kpoints, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    deallocate (m, stat=istat)
    if (istat /= 0) then
      call set_error_alloc(error, 'Error deallocating m in read_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
  end subroutine read_chkpt

  subroutine input_setopt(common_data, seedname, comm, istdout, istderr, ierr)
#ifdef MPI08
    use mpi_f08
#endif
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_alloc, set_error_dealloc, set_error_fatal
    use w90_kmesh, only: kmesh_get
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, &
      w90_wannier90_readwrite_read_special, w90_extra_io_type

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data
#ifdef MPI08
    type(mpi_comm), intent(in) :: comm
#else
    integer, intent(in) :: comm
#endif

    ! local variables
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    logical :: cp_pp, disentanglement

    ierr = 0

    if (allocated(common_data%settings%in_data)) then
      call set_error_fatal(error, ' readinput and setopt clash at input_setopt call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    else if (.not. allocated(common_data%settings%entries)) then
      call set_error_fatal(error, ' input_setopt called with no input', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif

    common_data%comm%comm = comm ! set communicator
    common_data%seedname = seedname ! set seedname for input/output files

    call w90_wannier90_readwrite_read_special(common_data%settings, common_data%atom_data, &
                                              common_data%kmesh_input, common_data%kmesh_info, &
                                              common_data%kpt_latt, common_data%wann_control, &
                                              common_data%proj, common_data%proj_input, &
                                              common_data%select_proj, common_data%w90_system, &
                                              common_data%w90_calculation, &
                                              common_data%real_lattice, common_data%physics%bohr, &
                                              common_data%mp_grid, common_data%num_bands, &
                                              common_data%num_kpts, common_data%num_proj, &
                                              common_data%num_wann, common_data%gamma_only, &
                                              common_data%lhasproj, &
                                              common_data%use_bloch_phases, &
                                              common_data%dist_kpoints, istdout, error, &
                                              common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    disentanglement = (common_data%num_bands > common_data%num_wann)

    if (disentanglement) then
      allocate (common_data%dis_manifold%ndimwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating ndimwin in input_setopt() library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%nfirstwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nfirstwin in input_setopt() library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%lwindow(common_data%num_bands, common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating lwindow in input_setopt() library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif

    allocate (common_data%wannier_data%centres(3, common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wannier_centres in input_setopt() library call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%centres = 0.0_dp

    allocate (common_data%wannier_data%spreads(common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wannier_spreads in input_setopt() library call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%spreads = 0.0_dp

    ! read all other variables
    call w90_wannier90_readwrite_read(common_data%settings, common_data%band_plot, &
                                      common_data%dis_control, common_data%dis_spheres, &
                                      common_data%dis_manifold, common_data%exclude_bands, &
                                      common_data%fermi_energy_list, &
                                      common_data%fermi_surface_data, common_data%output_file, &
                                      common_data%wvfn_read, common_data%wann_control, &
                                      common_data%real_space_ham, common_data%kpoint_path, &
                                      common_data%w90_system, common_data%tran, &
                                      common_data%print_output, common_data%wann_plot, io_params, &
                                      common_data%ws_region, common_data%real_lattice, &
                                      common_data%w90_calculation, common_data%physics%bohr, &
                                      common_data%sitesym%symmetrize_eps, common_data%num_bands, &
                                      common_data%num_kpts, common_data%num_wann, &
                                      common_data%optimisation, common_data%calc_only_A, cp_pp, &
                                      common_data%gamma_only, common_data%lsitesymmetry, &
                                      common_data%use_bloch_phases, common_data%seedname, istdout, &
                                      error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    ! clear any settings (from settings interface not .win file)
    deallocate (common_data%settings%entries, stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in deallocating entries data in input_setopt() library call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
    endif
  end subroutine input_setopt

  subroutine input_reader(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_input, set_error_fatal, set_error_alloc
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_clean_infile
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    logical :: cp_pp

    ierr = 0

    if (allocated(common_data%settings%entries)) then
      call set_error_fatal(error, &
                           'Error: input reader called when unspent options present (setopt must be called first)', &
                           common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    ! read data from .win file to internal string array
    call w90_readwrite_in_file(common_data%settings, common_data%seedname, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    ! set options corresponding to string array from .win file
    call w90_wannier90_readwrite_read(common_data%settings, common_data%band_plot, &
                                      common_data%dis_control, common_data%dis_spheres, &
                                      common_data%dis_manifold, common_data%exclude_bands, &
                                      common_data%fermi_energy_list, &
                                      common_data%fermi_surface_data, common_data%output_file, &
                                      common_data%wvfn_read, common_data%wann_control, &
                                      common_data%real_space_ham, common_data%kpoint_path, &
                                      common_data%w90_system, common_data%tran, &
                                      common_data%print_output, common_data%wann_plot, io_params, &
                                      common_data%ws_region, common_data%real_lattice, &
                                      common_data%w90_calculation, common_data%physics%bohr, &
                                      common_data%sitesym%symmetrize_eps, common_data%num_bands, &
                                      common_data%num_kpts, common_data%num_wann, &
                                      common_data%optimisation, common_data%calc_only_A, cp_pp, &
                                      common_data%gamma_only, common_data%lsitesymmetry, &
                                      common_data%use_bloch_phases, common_data%seedname, &
                                      istdout, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    ! remove any remaining acceptable keywords; anything that remains is an input error
    call w90_readwrite_clean_infile(common_data%settings, istdout, common_data%seedname, error, &
                                    common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (allocated(common_data%settings%in_data)) deallocate (common_data%settings%in_data)
  end subroutine input_reader

  subroutine input_reader_special(common_data, seedname, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_input, set_error_fatal, set_error_alloc
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_clean_infile
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read_special, w90_extra_io_type

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    logical :: disentanglement

    ierr = 0

    ! read data from .win file to internal string array
    call w90_readwrite_in_file(common_data%settings, seedname, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    call w90_wannier90_readwrite_read_special(common_data%settings, common_data%atom_data, &
                                              common_data%kmesh_input, common_data%kmesh_info, &
                                              common_data%kpt_latt, common_data%wann_control, &
                                              common_data%proj, common_data%proj_input, &
                                              common_data%select_proj, common_data%w90_system, &
                                              common_data%w90_calculation, &
                                              common_data%real_lattice, common_data%physics%bohr, &
                                              common_data%mp_grid, common_data%num_bands, &
                                              common_data%num_kpts, common_data%num_proj, &
                                              common_data%num_wann, common_data%gamma_only, &
                                              common_data%lhasproj, &
                                              common_data%use_bloch_phases, &
                                              common_data%dist_kpoints, istdout, error, &
                                              common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    common_data%seedname = seedname

    disentanglement = (common_data%num_bands > common_data%num_wann)

    if (disentanglement) then
      allocate (common_data%dis_manifold%ndimwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating ndimwin in input_reader_special() call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%nfirstwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nfirstwin in input_reader_special() call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%lwindow(common_data%num_bands, common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating lwindow in input_reader_special() call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif
    allocate (common_data%wannier_data%centres(3, common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wannier_centres in input_reader_special() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%centres = 0.0_dp
    allocate (common_data%wannier_data%spreads(common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wannier_spreads in input_reader_special() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%spreads = 0.0_dp

    ! remove any remaining acceptable keywords; anything that remains is an input error
    call w90_readwrite_clean_infile(common_data%settings, istdout, seedname, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (allocated(common_data%settings%in_data)) deallocate (common_data%settings%in_data)
  end subroutine input_reader_special

  subroutine input_print_details(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_readwrite, only: w90_readwrite_write_header
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write, w90_extra_io_type

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params ! what is this? fixme

    ierr = 0

    ! write jazzy header info
    call w90_readwrite_write_header(common_data%physics%bohr_version_str, &
                                    common_data%physics%constants_version_str1, &
                                    common_data%physics%constants_version_str2, istdout)

    ! write simulation details
    call w90_wannier90_readwrite_write(common_data%atom_data, common_data%band_plot, &
                                       common_data%dis_control, common_data%dis_spheres, &
                                       common_data%fermi_energy_list, &
                                       common_data%fermi_surface_data, common_data%kpt_latt, &
                                       common_data%output_file, common_data%wvfn_read, &
                                       common_data%wann_control, common_data%proj, &
                                       common_data%proj_input, common_data%real_space_ham, &
                                       common_data%select_proj, common_data%kpoint_path, &
                                       common_data%tran, common_data%print_output, &
                                       common_data%wannier_data, &
                                       common_data%wann_plot, io_params, &
                                       common_data%w90_calculation, common_data%real_lattice, &
                                       common_data%sitesym%symmetrize_eps, common_data%mp_grid, &
                                       common_data%num_bands, common_data%num_kpts, &
                                       common_data%num_proj, common_data%num_wann, &
                                       common_data%optimisation, .false., &
                                       common_data%gamma_only, common_data%lsitesymmetry, &
                                       common_data%w90_system%spinors, &
                                       common_data%use_bloch_phases, istdout)

    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine input_print_details

  subroutine create_kmesh(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_kmesh, only: kmesh_get

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    call kmesh_get(common_data%kmesh_input, common_data%kmesh_info, common_data%print_output, &
                   common_data%kpt_latt, common_data%real_lattice, common_data%num_kpts, &
                   common_data%gamma_only, istdout, common_data%timer, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%setup_complete = .true.
  end subroutine create_kmesh

  subroutine write_kmesh(common_data, istdout, istderr, ierr)
    use w90_comms, only: mpirank, comms_sync_err
    use w90_error_base, only: w90_error_type
    use w90_kmesh, only: kmesh_get, kmesh_write

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. common_data%setup_complete) then
      call create_kmesh(common_data, istdout, istderr, ierr)
      if (ierr > 0) return
    endif

    if (mpirank(common_data%comm) == 0) then
      call kmesh_write(common_data%exclude_bands, common_data%kmesh_info, &
                       common_data%select_proj%auto_projections, common_data%proj_input, &
                       common_data%print_output, common_data%kpt_latt, common_data%real_lattice, &
                       common_data%num_kpts, common_data%num_proj, common_data%calc_only_A, &
                       common_data%w90_system%spinors, common_data%seedname, common_data%timer)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif
    call comms_sync_err(common_data%comm, error, 0) ! this is necessary since non-root may never enter an mpi collective if root has exited here
  end subroutine write_kmesh

  subroutine overlaps(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_overlap, only: overlap_read

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    logical :: cp_pp = .false.
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. common_data%setup_complete) then
      call create_kmesh(common_data, istdout, istderr, ierr)
      if (ierr > 0) return
    endif

    if (common_data%num_bands > common_data%num_wann) then ! disentanglement case
      call overlap_read(common_data%kmesh_info, common_data%select_proj, common_data%sitesym, common_data%u_opt, &
                        common_data%m_orig, common_data%num_bands, common_data%num_kpts, common_data%num_proj, &
                        common_data%num_wann, common_data%print_output, common_data%print_output%timing_level, cp_pp, &
                        common_data%gamma_only, common_data%lsitesymmetry, common_data%use_bloch_phases, &
                        common_data%seedname, istdout, common_data%timer, common_data%dist_kpoints, error, &
                        common_data%comm)
    else
      if (.not. associated(common_data%u_matrix)) then
        call set_error_fatal(error, 'u_matrix not associated at overlaps library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      if (.not. associated(common_data%m_matrix_local)) then
        call set_error_fatal(error, 'm_matrix_local not associated at overlaps library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif

      call overlap_read(common_data%kmesh_info, common_data%select_proj, common_data%sitesym, &
                        common_data%u_matrix, common_data%m_matrix_local, common_data%num_bands, &
                        common_data%num_kpts, common_data%num_proj, common_data%num_wann, &
                        common_data%print_output, common_data%print_output%timing_level, cp_pp, &
                        common_data%gamma_only, common_data%lsitesymmetry, &
                        common_data%use_bloch_phases, common_data%seedname, istdout, &
                        common_data%timer, common_data%dist_kpoints, error, common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine overlaps

  subroutine disentangle(common_data, istdout, istderr, ierr)
    use w90_disentangle, only: dis_main, setup_m_loc
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal

    implicit none

    ! arguments
    type(lib_common_type), intent(inout) :: common_data
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. associated(common_data%m_orig)) then  ! m_matrix_orig_local (nband*nwann for disentangle)
      call set_error_fatal(error, 'm_orig not associated for disentangle call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%m_matrix_local)) then ! (nband*nwann*nknode for wannierise)
      call set_error_fatal(error, 'm_matrix_local not associated for disentangle call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not associated for disentangle call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%u_opt)) then
      call set_error_fatal(error, 'u_opt not associated for disentangle call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%eigval)) then
      call set_error_fatal(error, 'eigval not associated for disentangle call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    call dis_main(common_data%dis_control, common_data%dis_spheres, common_data%dis_manifold, &
                  common_data%kmesh_info, common_data%kpt_latt, common_data%sitesym, &
                  common_data%print_output, common_data%m_orig, common_data%u_matrix, &
                  common_data%u_opt, common_data%eigval, common_data%real_lattice, &
                  common_data%omega%invariant, common_data%num_bands, common_data%num_kpts, &
                  common_data%num_wann, common_data%gamma_only, common_data%lsitesymmetry, &
                  istdout, common_data%timer, common_data%dist_kpoints, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    ! copy to m_matrix_local from m_matrix_orig_local (aka m_orig)
    call setup_m_loc(common_data%kmesh_info, common_data%print_output, common_data%m_matrix_local, &
                     common_data%m_orig, common_data%u_matrix, common_data%num_bands, &
                     common_data%num_kpts, common_data%num_wann, common_data%timer, &
                     common_data%dist_kpoints, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    common_data%have_disentangled = .true.
  end subroutine disentangle

  subroutine projovlp(common_data, istdout, istderr, ierr) !fixme(jj) rename more sensibly
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_overlap, only: overlap_project

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. associated(common_data%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for projovlp call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for projovlp call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (.not. common_data%have_disentangled) then
      ! fixme (jj) there is also a gamma only specialisation of this
      call overlap_project(common_data%sitesym, common_data%m_matrix_local, common_data%u_matrix, &
                           common_data%kmesh_info%nnlist, common_data%kmesh_info%nntot, &
                           common_data%num_wann, common_data%num_kpts, common_data%num_wann, &
                           common_data%print_output%timing_level, common_data%lsitesymmetry, &
                           istdout, common_data%timer, common_data%dist_kpoints, error, &
                           common_data%comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif
  end subroutine projovlp

  subroutine wannierise(common_data, istdout, istderr, ierr)
    use w90_comms, only: mpirank, comms_barrier
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_wannierise, only: wann_main, wann_main_gamma

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. associated(common_data%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for wannierise', common_data%comm)
    else if (.not. associated(common_data%u_opt)) then
      call set_error_fatal(error, 'u_opt not set for wannierise', common_data%comm)
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for wannierise', common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (common_data%gamma_only) then
      if (mpirank(common_data%comm) == 0) then
        call wann_main_gamma(common_data%kmesh_info, common_data%wann_control, &
                             common_data%omega, common_data%print_output, &
                             common_data%wannier_data, common_data%m_matrix_local, &
                             common_data%u_matrix, common_data%real_lattice, common_data%num_kpts, &
                             common_data%num_wann, istdout, common_data%timer, error, &
                             common_data%comm)
        if (allocated(error)) then
          call prterr(error, ierr, istdout, istderr, common_data%comm)
          return
        endif
      endif
      call comms_barrier(error, common_data%comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    else
      call wann_main(common_data%ham_logical, common_data%kmesh_info, common_data%kpt_latt, &
                     common_data%wann_control, common_data%omega, common_data%sitesym, &
                     common_data%print_output, common_data%wannier_data, common_data%ws_region, &
                     common_data%w90_calculation, common_data%ham_k, common_data%ham_r, &
                     common_data%m_matrix_local, common_data%u_matrix, common_data%real_lattice, &
                     common_data%wannier_centres_translated, common_data%irvec, &
                     common_data%mp_grid, common_data%ndegen, common_data%nrpts, &
                     common_data%num_kpts, common_data%num_proj, common_data%num_wann, &
                     common_data%optimisation, common_data%rpt_origin, &
                     common_data%band_plot%mode, common_data%tran%mode, &
                     common_data%lsitesymmetry, istdout, common_data%timer, &
                     common_data%dist_kpoints, error, common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine wannierise

  subroutine plot_files(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_plot, only: plot_main

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data ! inout due to ham_logical only

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    ! fixme(jj) what are our preconditions?

    call plot_main(common_data%atom_data, common_data%band_plot, common_data%dis_manifold, &
                   common_data%fermi_energy_list, common_data%fermi_surface_data, &
                   common_data%ham_logical, common_data%kmesh_info, common_data%kpt_latt, &
                   common_data%output_file, common_data%wvfn_read, common_data%real_space_ham, &
                   common_data%kpoint_path, common_data%print_output, common_data%wannier_data, &
                   common_data%wann_plot, common_data%ws_region, common_data%w90_calculation, &
                   common_data%ham_k, common_data%ham_r, common_data%m_matrix_local, &
                   common_data%u_matrix, common_data%u_opt, common_data%eigval, &
                   common_data%real_lattice, common_data%wannier_centres_translated, &
                   common_data%physics%bohr, common_data%irvec, common_data%mp_grid, &
                   common_data%ndegen, common_data%shift_vec, common_data%nrpts, &
                   common_data%num_bands, common_data%num_kpts, common_data%num_wann, &
                   common_data%rpt_origin, common_data%tran%mode, common_data%have_disentangled, &
                   common_data%lsitesymmetry, common_data%w90_system, common_data%seedname, &
                   istdout, common_data%timer, common_data%dist_kpoints, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine plot_files

  subroutine transport(common_data, istdout, istderr, ierr)
    use w90_comms, only: mpirank
    use w90_error_base, only: w90_error_type
    use w90_transport, only: tran_main

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data ! because of ham_logical

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    ! fixme(jj) what are our preconditions?

    ! currently tran_main is entirely serial
    ! test throwing an error!
    if (mpirank(common_data%comm) == 0) then
      call tran_main(common_data%atom_data, common_data%dis_manifold, common_data%fermi_energy_list, &
                     common_data%ham_logical, common_data%kpt_latt, common_data%output_file, &
                     common_data%real_space_ham, common_data%tran, common_data%print_output, &
                     common_data%wannier_data, common_data%ws_region, common_data%w90_calculation, &
                     common_data%ham_k, common_data%ham_r, common_data%u_matrix, common_data%u_opt, &
                     common_data%eigval, common_data%real_lattice, common_data%wannier_centres_translated, &
                     common_data%irvec, common_data%mp_grid, common_data%ndegen, common_data%shift_vec, &
                     common_data%nrpts, common_data%num_bands, common_data%num_kpts, common_data%num_wann, &
                     common_data%rpt_origin, common_data%band_plot%mode, common_data%have_disentangled, &
                     common_data%lsitesymmetry, common_data%seedname, istdout, common_data%timer, error, &
                     common_data%comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif
  end subroutine transport

  subroutine print_times(common_data, istdout)
    use w90_io, only: io_print_timings
    implicit none

    type(lib_common_type), intent(in) :: common_data
    integer, intent(in) :: istdout

    if (common_data%print_output%iprint > 0) call io_print_timings(common_data%timer, istdout)
  end subroutine print_times

  subroutine set_m_orig(common_data, m_orig) ! m_matrix_local_orig
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: m_orig(:, :, :, :)

    common_data%m_matrix_local => m_orig
    common_data%m_orig => m_orig
  end subroutine set_m_orig

  subroutine set_u_matrix(common_data, u_matrix)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_matrix(:, :, :)

    common_data%u_matrix => u_matrix
  end subroutine set_u_matrix

  subroutine set_u_opt(common_data, u_opt)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_opt(:, :, :)

    common_data%u_opt => u_opt
  end subroutine set_u_opt

  subroutine set_eigval(common_data, eigval)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    real(kind=dp), intent(in), target :: eigval(:, :)

    common_data%eigval => eigval
    ! if not already initialised, set disentanglement window to limits of spectrum
    if (common_data%dis_manifold%win_min == -huge(0.0_dp)) common_data%dis_manifold%win_min = minval(common_data%eigval)
    if (common_data%dis_manifold%win_max == huge(0.0_dp)) common_data%dis_manifold%win_max = maxval(common_data%eigval)
    if (common_data%dis_manifold%frozen_states) then
      if (common_data%dis_manifold%froz_min == -huge(0.0_dp)) &
        common_data%dis_manifold%froz_min = common_data%dis_manifold%win_min
    endif
  end subroutine set_eigval

  subroutine set_parallel_comms(common_data, comm)
#ifdef MPI08
    use mpi_f08
#endif
    implicit none
    type(lib_common_type), intent(inout) :: common_data
#ifdef MPI08
    type(mpi_comm), intent(in) :: comm
#else
    integer, intent(in) :: comm
#endif
    common_data%comm%comm = comm
  end subroutine set_parallel_comms

  subroutine set_kpoint_distribution(common_data, dist, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal

    implicit none

    ! arguments
    integer, intent(in) :: istderr, istdout, dist(:)
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (size(dist) < 1) then
      call set_error_fatal(error, 'Error in k-point distribution, mpisize < 1', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%dist_kpoints = dist
  end subroutine set_kpoint_distribution

  subroutine set_constant_bohr_to_ang(common_data, bohr_to_angstrom)
    ! used to set the bohr_to_angstrom value as used in the SCF code

    implicit none

    type(lib_common_type), intent(inout) :: common_data
    real(kind=dp), intent(in) :: bohr_to_angstrom

    common_data%physics%bohr = bohr_to_angstrom
    common_data%physics%bohr_version_str = "-> Using Bohr value from linked main code"
  end subroutine set_constant_bohr_to_ang

  subroutine set_option_text(common_data, keyword, text)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword, text
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%txtdata = text
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_text

  subroutine set_option_logical(common_data, keyword, bool)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    logical, intent(in) :: bool
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%ldata = bool
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_logical

  subroutine set_option_i1d(common_data, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    integer, intent(in) :: arr(:)
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%i1d = arr ! this causes an automatic allocation
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_i1d

  subroutine set_option_i2d(common_data, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    integer, intent(in) :: arr(:, :)
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%i2d = arr
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_i2d

  subroutine set_option_int(common_data, keyword, ival)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    integer, intent(in) :: ival
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%idata = ival
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_int

  subroutine set_option_r1d(common_data, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    real(kind=dp), intent(in) :: arr(:)
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%r1d = arr
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_r1d

  subroutine set_option_r2d(common_data, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    real(kind=dp), intent(in) :: arr(:, :)
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%r2d = arr
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_r2d

  subroutine set_option_real(common_data, keyword, rval)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    real(kind=dp), intent(in) :: rval
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%rdata = rval
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) call expand_settings(common_data%settings)
  endsubroutine set_option_real

  subroutine prterr(error, ie, istdout, istderr, comm)
    use w90_comms, only: comms_no_sync_send, comms_no_sync_recv, w90_comm_type, mpirank, mpisize
    use w90_error_base, only: code_remote, w90_error_type

    ! arguments
    integer, intent(inout) :: ie! global error value to be returned
    integer, intent(in) :: istderr, istdout
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(in) :: error

    ! local variables
    type(w90_error_type), allocatable :: le ! unchecked error state for calls made in this routine
    integer :: je ! error value on remote ranks
    integer :: j ! rank index
    integer :: failrank ! lowest rank reporting an error
    character(len=128) :: mesg ! only print 128 chars of error

    ie = 0
    mesg = 'not set'

    if (mpirank(comm) == 0) then
      ! fixme, report all failing ranks instead of lowest failing rank (current stand)
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

  subroutine read_eigvals(common_data, eigval, istdout, istderr, ierr)
    use w90_error, only: w90_error_type, set_error_fatal
    use w90_readwrite, only: w90_readwrite_read_eigvals

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    real(kind=dp), intent(inout) :: eigval(:, :)
    type(lib_common_type), intent(inout) :: common_data
    integer, intent(out) :: ierr

    ! local vars
    type(w90_error_type), allocatable :: error
    logical :: eig_found

    ierr = 0

    if (size(eigval, 1) /= common_data%num_bands) then
      call set_error_fatal(error, 'eigval not dimensioned correctly (num_bands,num_kpts) in read_eigvals', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    elseif (size(eigval, 2) /= common_data%num_kpts) then
      call set_error_fatal(error, 'eigval not dimensioned correctly (num_bands,num_kpts) in read_eigvals', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    call w90_readwrite_read_eigvals(eig_found, eigval, common_data%num_bands, common_data%num_kpts, &
                                    istdout, common_data%seedname, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. eig_found) then
      call set_error_fatal(error, 'failed to read eigenvalues file in read_eigvals call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine read_eigvals

  subroutine get_nn(common_data, istdout, istderr, nn, ierr)
    implicit none

    ! arguments
    integer, intent(out) :: nn, ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data

    if (.not. common_data%setup_complete) then
      call create_kmesh(common_data, istdout, istderr, ierr)
      ! fixme check ierr
    endif

    nn = common_data%kmesh_info%nntot
  end subroutine get_nn

  subroutine get_nnkp(common_data, istdout, istderr, nnkp, ierr)
    implicit none

    integer, intent(out) :: nnkp(:, :), ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data

    if (.not. common_data%setup_complete) then
      call create_kmesh(common_data, istdout, istderr, ierr)
      ! fixme check ierr
    endif

    nnkp = common_data%kmesh_info%nnlist
  end subroutine get_nnkp

  subroutine get_centres(common_data, centres)
    implicit none

    real(kind=dp), intent(out) :: centres(:, :)
    type(lib_common_type), intent(in) :: common_data

    if (.not. common_data%setup_complete) then
      ! fixme error condition
    endif

    centres = common_data%wannier_data%centres
  endsubroutine get_centres

  subroutine get_spreads(common_data, spreads)
    implicit none

    real(kind=dp), intent(out) :: spreads(:)
    type(lib_common_type), intent(in) :: common_data

    if (.not. common_data%setup_complete) then
      ! fixme error condition
    endif

    spreads = common_data%wannier_data%spreads
  endsubroutine get_spreads

  subroutine get_proj(common_data, n, site, l, m, s, istdout, istderr, ierr, sqa, z, x, rad, zona)
    use w90_error, only: w90_error_type, set_error_fatal
    implicit none

    ! returns arrays describing a list of projections derived from the .win input file
    ! n specifies the number of projections returned
    ! all other arguments should be of length >= num_band (the maximum possible size?)
    integer, intent(inout) :: n, l(:), m(:), s(:)
    real(kind=dp), intent(inout) :: site(:, :)

    type(lib_common_type), intent(in), target :: common_data
    integer, intent(out) :: ierr
    integer, intent(in) :: istdout, istderr

    ! probably the remaining variables find limited use, allow them to be absent
    integer, intent(inout), optional :: rad(:)
    real(kind=dp), intent(inout), optional :: sqa(:, :), z(:, :), x(:, :), zona(:)

    ! local variables
    integer :: ip
    type(proj_type), pointer :: proj
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. allocated(common_data%proj_input)) then
      call set_error_fatal(error, 'Projectors are not setup in Wannier90 library when requested via get_proj()', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    n = size(common_data%proj_input)

    ! check allocation of main output arrays
    if (size(l) < n) then
      call set_error_fatal(error, 'Array argument l in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (size(m) < n) then
      call set_error_fatal(error, 'Array argument m in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (size(s) < n) then
      call set_error_fatal(error, 'Array argument s in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (size(site, 2) < n) then
      call set_error_fatal(error, 'Array argument site in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (present(sqa) .and. size(sqa, 2) < n) then
      call set_error_fatal(error, 'Optional array argument sqa in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (present(sqa) .and. size(sqa, 2) < n) then
      call set_error_fatal(error, 'Optional array argument sqa in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (present(z) .and. size(z, 2) < n) then
      call set_error_fatal(error, 'Optional array argument z in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (present(x) .and. size(x, 2) < n) then
      call set_error_fatal(error, 'Optional array argument x in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (present(rad) .and. size(rad) < n) then
      call set_error_fatal(error, 'Optional array argument rad in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (present(zona) .and. size(zona) < n) then
      call set_error_fatal(error, 'Optional array argument zona in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    do ip = 1, n
      proj => common_data%proj_input(ip)
      l(ip) = proj%l
      m(ip) = proj%m
      s(ip) = proj%s
      site(:, ip) = proj%site(:)
      if (present(sqa)) sqa(:, ip) = proj%s_qaxis(:)
      if (present(z)) z(:, ip) = proj%z(:)
      if (present(x)) x(:, ip) = proj%x(:)
      if (present(rad)) rad(ip) = proj%radial
      if (present(zona)) zona(ip) = proj%zona
    enddo
  end subroutine get_proj
end module w90_library
