
module w90_library

  ! as with fortran routines like allocate, the status variable indicates an error if non-zero
  ! positive is an error, negative is a warning (such as non-convergence) which is recoverable.

  ! note on the naming of output and error-output streams
  !   the specific names stdout and stderr cannot be used because
  !   the python wrapper uses an intermediate c representation with named arguments
  !   in this wrapper stdout/err are already defined as structs inconsistent
  !   which are inconsistent with the integer unit numbers used here.
  !   (The slightly 77 style naming here has no particular significance.)

  use w90_constants, only: dp, w90_physical_constants_type
  use w90_types
  use w90_wannier90_types
  use w90_comms, only: w90_comm_type
  use w90_io, only: prterr
  use iso_c_binding

  implicit none

  private :: dp ! avoid polluting calling program's namespace (dp is defined in w90_constants)
  private :: prterr

  ! datatype encapsulating types used by wannier90
  type lib_common_type
    character(len=128) :: seedname

    !! matrices
    complex(kind=dp), allocatable :: ham_k(:, :, :)
    complex(kind=dp), allocatable :: ham_r(:, :, :)
    complex(kind=dp), pointer :: a_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: m_matrix_local(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_matrix(:, :, :, :) => null()
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

  !private :: w90_create_kmesh ! trigers the generation of k-mesh info (as do get_nn*)
  public :: w90_create_kmesh ! trigers the generation of k-mesh info (as do get_nn*)
  public :: w90_disentangle ! perform disentanglement
  public :: w90_get_centres ! get wannier centers
  public :: w90_get_fortran_file
  public :: w90_get_fortran_stderr
  public :: w90_get_fortran_stdout
  public :: w90_get_proj ! get projection info
  public :: w90_get_spreads ! get spreads
  public :: w90_get_nnkp ! get k' indexes
  public :: w90_get_nn ! get number of b-vectors
  public :: w90_get_gkpb ! get g offsets of k'
  public :: w90_input_reader ! extra variables from .win
  public :: w90_input_setopt ! set options specified by set_option
  public :: w90_plot ! plot functions
  public :: w90_project_overlap ! transform overlaps and initial projections
  public :: w90_set_constant_bohr_to_ang ! set value of angstrom
  public :: w90_set_eigval ! set (pointer to) eigenvalues
  public :: w90_set_m_local ! set (pointer to) m (decomposed)
  public :: w90_set_option ! set options
  public :: w90_set_u_matrix ! u matrix (nw,nw)
  public :: w90_set_u_opt ! u (nb,nb)
  public :: w90_transport ! transport functions
  public :: w90_wannierise ! perform wannierisation

  interface w90_set_option
    module procedure w90_set_option_logical
    !module procedure w90_set_option_b1d
    module procedure w90_set_option_text
    module procedure w90_set_option_int
    module procedure w90_set_option_i1d
    module procedure w90_set_option_i2d
    module procedure w90_set_option_r1d
    module procedure w90_set_option_r2d
    module procedure w90_set_option_real
  end interface w90_set_option

contains

  subroutine w90_get_fortran_stdout(istdout)
    use iso_fortran_env, only: output_unit
    implicit none
    integer, intent(out) :: istdout
    istdout = output_unit
  end subroutine w90_get_fortran_stdout

  subroutine w90_get_fortran_stderr(istdout)
    use iso_fortran_env, only: error_unit
    implicit none
    integer, intent(out) :: istdout
    istdout = error_unit
  end subroutine w90_get_fortran_stderr

  subroutine w90_get_fortran_file(output, name)
    implicit none
    integer, intent(out) :: output
    character(len=*), intent(in) :: name
    open (newunit=output, file=name, form='formatted', status='unknown')
  end subroutine w90_get_fortran_file

  subroutine w90_input_setopt(common_data, seedname, comm, istdout, istderr, ierr)
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
      return
    else if (.not. allocated(common_data%settings%entries)) then
      call set_error_fatal(error, ' input_setopt called with no input', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
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
      return
    endif
  end subroutine w90_input_setopt

  subroutine w90_input_reader(common_data, istdout, istderr, ierr)
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
  end subroutine w90_input_reader

  subroutine w90_disentangle(common_data, istdout, istderr, ierr)
    use w90_disentangle_mod, only: dis_main, setup_m_loc
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

    ! m_matrix_orig_local (nband*nwann for disentangle)
    if (.not. associated(common_data%m_matrix_local)) then ! (nband*nwann*nknode for wannierise)
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

    if (common_data%num_bands > common_data%num_wann) then
      call dis_main(common_data%dis_control, common_data%dis_spheres, common_data%dis_manifold, &
                    common_data%kmesh_info, common_data%kpt_latt, common_data%sitesym, &
                    common_data%print_output, common_data%m_matrix_local, common_data%u_matrix, &
                    common_data%u_opt, common_data%eigval, common_data%real_lattice, &
                    common_data%omega%invariant, common_data%num_bands, common_data%num_kpts, &
                    common_data%num_wann, common_data%gamma_only, common_data%lsitesymmetry, &
                    istdout, common_data%timer, common_data%dist_kpoints, error, common_data%comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif

      call setup_m_loc(common_data%kmesh_info, common_data%print_output, common_data%m_matrix_local, &
                       common_data%m_matrix_local, common_data%u_matrix, common_data%num_bands, &
                       common_data%num_kpts, common_data%num_wann, common_data%timer, &
                       common_data%dist_kpoints, error, common_data%comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif

      common_data%have_disentangled = .true.
    endif
  end subroutine w90_disentangle

  subroutine w90_project_overlap(common_data, istdout, istderr, ierr) !fixme(jj) rename more sensibly
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
  end subroutine w90_project_overlap

  subroutine w90_wannierise(common_data, istdout, istderr, ierr)
    use w90_comms, only: mpirank, comms_barrier
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_wannierise_mod, only: wann_main, wann_main_gamma

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
  end subroutine w90_wannierise

  subroutine w90_plot(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_plot_mod, only: plot_main

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
  end subroutine w90_plot

  subroutine w90_transport(common_data, istdout, istderr, ierr)
    use w90_comms, only: mpirank
    use w90_error_base, only: w90_error_type
    use w90_transport_mod, only: tran_main

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
  end subroutine w90_transport

  subroutine w90_set_m_local(common_data, m_orig) ! m_matrix_local_orig
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: m_orig(:, :, :, :)

    common_data%m_matrix_local => m_orig
    !common_data%m_orig => m_orig
  end subroutine w90_set_m_local

  subroutine w90_set_u_matrix(common_data, u_matrix)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_matrix(:, :, :)

    common_data%u_matrix => u_matrix
  end subroutine w90_set_u_matrix

  subroutine w90_set_u_opt(common_data, u_opt)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_opt(:, :, :)

    common_data%u_opt => u_opt
  end subroutine w90_set_u_opt

  subroutine w90_set_eigval(common_data, eigval)
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
  end subroutine w90_set_eigval

  subroutine w90_set_constant_bohr_to_ang(common_data, bohr_to_angstrom)
    ! used to set the bohr_to_angstrom value as used in the SCF code

    implicit none

    type(lib_common_type), intent(inout) :: common_data
    real(kind=dp), intent(in) :: bohr_to_angstrom

    common_data%physics%bohr = bohr_to_angstrom
    common_data%physics%bohr_version_str = "-> Using Bohr value from linked main code"
  end subroutine w90_set_constant_bohr_to_ang

  subroutine w90_set_option_text(common_data, keyword, text)
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
  endsubroutine w90_set_option_text

  subroutine w90_set_option_logical(common_data, keyword, bool)
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
  endsubroutine w90_set_option_logical

  subroutine w90_set_option_i1d(common_data, keyword, arr)
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
  endsubroutine w90_set_option_i1d

  subroutine w90_set_option_i2d(common_data, keyword, arr)
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
  endsubroutine w90_set_option_i2d

  subroutine w90_set_option_int(common_data, keyword, ival)
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
  endsubroutine w90_set_option_int

  subroutine w90_set_option_r1d(common_data, keyword, arr)
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
  endsubroutine w90_set_option_r1d

  subroutine w90_set_option_r2d(common_data, keyword, arr)
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
  endsubroutine w90_set_option_r2d

  subroutine w90_set_option_real(common_data, keyword, rval)
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
  endsubroutine w90_set_option_real

  subroutine w90_create_kmesh(common_data, istdout, istderr, ierr)
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
  end subroutine w90_create_kmesh

  subroutine w90_get_nn(common_data, nn, istdout, istderr, ierr)
    implicit none

    ! arguments
    integer, intent(out) :: nn, ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data

    if (.not. common_data%setup_complete) then
      call w90_create_kmesh(common_data, istdout, istderr, ierr)
      if (ierr > 0) return
    endif

    nn = common_data%kmesh_info%nntot
  end subroutine w90_get_nn

  subroutine w90_get_nnkp(common_data, nnkp, istdout, istderr, ierr)
    implicit none

    integer, intent(out) :: nnkp(:, :), ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data

    if (.not. common_data%setup_complete) then
      call w90_create_kmesh(common_data, istdout, istderr, ierr)
      if (ierr > 0) return
    endif

    nnkp = common_data%kmesh_info%nnlist
  end subroutine w90_get_nnkp

  subroutine w90_get_gkpb(common_data, gkpb, istdout, istderr, ierr)
    ! gkpb is the triple of reciprocal lattice translations to correct phase in k'= k+b
    ! gkpb should be dimensioned (3,nk,nnb)
    implicit none

    integer, intent(out) :: gkpb(:, :, :), ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data

    if (.not. common_data%setup_complete) then
      call w90_create_kmesh(common_data, istdout, istderr, ierr)
      if (ierr > 0) return
    endif

    gkpb = common_data%kmesh_info%nncell
  end subroutine w90_get_gkpb

  subroutine w90_get_centres(common_data, centres)
    implicit none

    real(kind=dp), intent(out) :: centres(:, :)
    type(lib_common_type), intent(in) :: common_data

    centres = common_data%wannier_data%centres
  endsubroutine w90_get_centres

  subroutine w90_get_spreads(common_data, spreads)
    implicit none

    real(kind=dp), intent(out) :: spreads(:)
    type(lib_common_type), intent(in) :: common_data

    spreads = common_data%wannier_data%spreads
  endsubroutine w90_get_spreads

  subroutine w90_get_proj(common_data, n, site, l, m, s, istdout, istderr, ierr, sqa, z, x, rad, zona)
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
  end subroutine w90_get_proj

end module w90_library
