!-*- mode: F90 -*-!
!------------------------------------------------------------!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Please cite                                                !
!                                                            !
!  [ref] "Wannier90 as a community code:                     !
!        new features and applications",                     !
!        G. Pizzi et al.,  J. Phys. Cond. Matt. 32,          !
!        165902 (2020).                                      !
!        http://doi.org/10.1088/1361-648X/ab51ff             !
!                                                            !
! in any publications arising from the use of this code.     !
!                                                            !
! Wannier90 is based on Wannier77, written by N. Marzari,    !
! I. Souza and D. Vanderbilt. For the method please cite     !
!                                                            !
! [ref] N. Marzari and D. Vanderbilt,                        !
!       Phys. Rev. B 56 12847 (1997)                         !
!       http://dx.doi.org/10.1103/PhysRevB.56.12847          !
!                                                            !
! [ref] I. Souza, N. Marzari and D. Vanderbilt,              !
!       Phys. Rev. B 65 035109 (2001)                        !
!       http://dx.doi.org/10.1103/PhysRevB.65.035109         !
!                                                            !
! [ref] N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza,    !
!       D. Vanderbilt, "Maximally localized Wannier          !
!       functions: theory and applications",                 !
!       Rev. Mod. Phys. 84, 1419 (2012)                      !
!       http://dx.doi.org/10.1103/RevModPhys.84.1419         !
!                                                            !
! For a full list of authors and contributors, please        !
! see the README file in the root directory of the           !
! distribution.                                              !
!                                                            !
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
    !! base name for reading/writing of files

    ! matrices
    complex(kind=dp), allocatable :: ham_k(:, :, :)
    !! KS Hamiltonian, takes size (number of WF, number of WF, number of k-points)
    complex(kind=dp), allocatable :: ham_r(:, :, :)
    !! RS Hamiltonian, takes size (number of WF, number of WF, number of RS-points)
    complex(kind=dp), pointer :: a_matrix(:, :, :) => null()
    !! pointer to matrix "a" (projections), dimensioned (number of bands, number WF, number of k-points)
    complex(kind=dp), pointer :: m_matrix_local(:, :, :, :) => null()
    !! pointer to rank local part of "m", dimensioned (number of bands, number of bands, FD neighbours, number of k-points on this rank)
    complex(kind=dp), pointer :: u_matrix(:, :, :) => null()
    !! pointer to matrix "u", dimensionsed (number of WF, number of WF, number of k-points)
    complex(kind=dp), pointer :: u_matrix_opt(:, :, :) => null()
    !! pointer to matrix "u_matrix_opt", dimensioned (number of bands, number of WF, number of k-points)
    real(kind=dp), pointer :: eigval(:, :) => null()
    !! pointer to eigenvalue array dimensioned (number of bands, number of k-points)

    integer, allocatable :: dist_kpoints(:)
    !! distribution of k-points; dist_kpoints(i) = rank operating on k-point i
    integer, allocatable :: exclude_bands(:)
    !! index of bands to exclude from the calculation: passed to DFT code to define space for overlaps, projections, etc
    integer, allocatable :: irvec(:, :)
    integer, allocatable :: ndegen(:)
    integer, allocatable :: shift_vec(:, :)

    integer :: mp_grid(3)
    !! Dimensions of the Monkhorst-Pack grid
    integer :: nrpts
    !! Dimensions of the RS grid
    integer :: num_bands
    !! Number of bands (greater or equal to number of WFs)
    integer :: num_kpts
    !! Total number of k-points
    integer :: num_proj = 0
    !! Number of projectors defined
    integer :: num_wann
    !! Number of WFs
    integer :: optimisation = 3
    !! Algorithm control (affects memory use/speed)
    integer :: rpt_origin

    logical :: calc_only_A = .false.
    logical :: gamma_only = .false.
    !! Select gamma_only branch of some algorithms
    logical :: have_disentangled = .false.
    !! Flag that disentanglment has been performed
    logical :: lhasproj = .false.
    !! Flag that projectors are defined
    logical :: lsitesymmetry = .false.
    !! Flag that symmetry-adapted WFs are to be calculated
    logical :: use_bloch_phases = .false.
    !! Flag to bypass disentanglement
    logical :: setup_complete = .false.
    !! Internal flag to indicate completion of setup

    real(kind=dp), allocatable :: fermi_energy_list(:)
    !! Array of energies around the Fermi energy
    real(kind=dp), allocatable :: kpt_latt(:, :)
    !! Tabulation of k-points in crystal coordinates
    real(kind=dp), allocatable :: wannier_centres_translated(:, :)
    real(kind=dp) :: real_lattice(3, 3)

    ! See types.F90 and wannier90_types.F90 for the use and composition of the different types
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
    type(settings_type) :: settings
    !! container for input file (.win) data and options set via library interface
    type(sitesym_type) :: sitesym
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
  !! container type for all variables used by the Wannier90 library

  public :: w90_print_info
  !! prints a wide variety of simulation parameters to stdout
  public :: w90_create_kmesh
  ! trigers the generation of k-mesh info (as do get_nn*)
  ! this is called by get_nnkp and get_gkpb
  public :: w90_disentangle
  !! perform disentanglement
  public :: w90_distribute_kpts
  !! provides an MPI k-point distribution for codes that don't have one
  public :: w90_get_centres
  !! get wannier centers
  public :: w90_get_fortran_file
  !! open a file and get the corresponding unit number
  public :: w90_get_fortran_stderr
  !! get a fortran unit number corresponding to standard error
  public :: w90_get_fortran_stdout
  !! get a fortran unit number corresponding to standard output
  public :: w90_get_proj
  !! get projection info (after interpreting projector specification strings)
  public :: w90_get_spreads
  !! get spreads
  public :: w90_get_nnkp
  !! get k' indexes for finite-difference scheme
  public :: w90_get_nn
  !! get number of b-vectors (finite-difference points)
  public :: w90_get_gkpb
  !! get g offsets of k'
  public :: w90_input_reader
  !! optionally read additional input variables from .win file
  public :: w90_input_setopt
  !! act upon (interpret & setup) options specified by set_option interface
  public :: w90_plot
  !! performs plot functions
  public :: w90_project_overlap
  !! transform overlaps and initial projections
  public :: w90_set_comm
  !! setup MPI communicator in parallel case
  public :: w90_set_constant_bohr_to_ang
  !! set value of Bohr/Angstrom conversion
  public :: w90_set_eigval
  !! set (pointer to) eigenvalues
  public :: w90_set_m_local
  !! set (pointer to) m (potentially MPI decomposed by k-points)
  public :: w90_set_option
  !! specify options to the library; set_option is overloaded to accept arguments of various types
  public :: w90_set_u_matrix
  !! set (pointer to)  u matrix (nw,nw)
  public :: w90_set_u_opt
  !! set (pointer to) optimised (disentangled) u matrix (nb,nb)
  public :: w90_transport
  !! perform transport functions
  public :: w90_wannierise
  !! perform wannierisation

  interface w90_set_option
    module procedure w90_set_option_logical
    !module procedure w90_set_option_b1d
    module procedure w90_set_option_text
    module procedure w90_set_option_int
    module procedure w90_set_option_i1d
    module procedure w90_set_option_i2d
    module procedure w90_set_option_r1d
    module procedure w90_set_option_r2d
    module procedure w90_set_option_c2d
    module procedure w90_set_option_real
  end interface w90_set_option

contains

  subroutine w90_get_fortran_stdout(istdout)
    !! fortran unit number for stdout
    use iso_fortran_env, only: output_unit
    implicit none
    integer, intent(out) :: istdout
    istdout = output_unit
  end subroutine w90_get_fortran_stdout

  subroutine w90_get_fortran_stderr(istdout)
    !! fortran unit number for stderr
    use iso_fortran_env, only: error_unit
    implicit none
    integer, intent(out) :: istdout
    istdout = error_unit
  end subroutine w90_get_fortran_stderr

  subroutine w90_get_fortran_file(output, name)
    !! open a (formatted) file and return the corresponding fortran unit number
    implicit none
    integer, intent(out) :: output
    character(len=*), intent(in) :: name
    open (newunit=output, file=name, form='formatted', status='unknown')
  end subroutine w90_get_fortran_file

  subroutine w90_input_setopt(common_data, seedname, istdout, istderr, ierr)
    !! mechanism to act upon options supplied to the library
    !! input is parsed and interpreted (any errors are identified) and
    !! the library data structure (variable common_data) is populated ready for use

    ! w90_input_setopt() processes options stored in common_data%settings
    ! w90_input_reader() processes options stored in common_data%in_data (from .win file, should be empty here)

#ifdef MPI08
    use mpi_f08
#endif
#ifdef MPI90
    use mpi
#endif

    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_alloc, set_error_fatal, code_mpi
    use w90_comms, only: w90_comm_type, valid_communicator
    use w90_kmesh, only: kmesh_get
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, &
      w90_wannier90_readwrite_read_special

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data
    !! instance of library data object, modified here

    ! local variables
    type(w90_error_type), allocatable :: error
    logical :: cp_pp

    ierr = 0

    if (.not. valid_communicator(common_data%comm)) then ! when MPI is not defined, valid_communicator returns true (see comms.F90)
      ! this is a problem: how do we exit using the parallel error handler when the communicator is unknown?
      write (istderr, *) ' Error: parallel Wannier90 library invoked with invalid communicator, exiting.  Use w90_set_comm()!'
      ierr = code_mpi
      return
    endif

    if (allocated(common_data%settings%in_data)) then
      call set_error_fatal(error, &
                           ' Error: w90_read_input() and w90_set_option() clash at w90_input_setopt() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. allocated(common_data%settings%entries)) then
      call set_error_fatal(error, ' Error: w90_input_setopt() called with no input set', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
      !else if (.not. allocated(common_data%dist_kpoints)) then
      !  call set_error_fatal(error, ' input_setopt called but distk unallocated', common_data%comm)
      !  call prterr(error, ierr, istdout, istderr, common_data%comm)
      !  return
    endif

    common_data%seedname = seedname ! set seedname for input/output files

    ! read_special can only be executed once per library object
    ! it sets key variables (eg, number of k-points sizing allocations, etc)
    call w90_wannier90_readwrite_read_special(common_data%settings, common_data%atom_data, &
                                              common_data%kmesh_input, common_data%kmesh_info, &
                                              common_data%kpt_latt, common_data%wann_control, &
                                              common_data%proj, common_data%proj_input, &
                                              common_data%select_proj, common_data%w90_system, &
                                              common_data%w90_calculation, &
                                              common_data%real_lattice, common_data%physics%bohr, &
                                              common_data%mp_grid, common_data%num_bands, &
                                              common_data%exclude_bands, &
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

    ! condition for disentanglement is number of bands > number of WF
    if (common_data%num_bands > common_data%num_wann) then
      allocate (common_data%dis_manifold%ndimwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating ndimwin in w90_input_setopt() library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%nfirstwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nfirstwin in w90_input_setopt() library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%lwindow(common_data%num_bands, common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating lwindow in w90_input_setopt() library call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif

    allocate (common_data%wannier_data%centres(3, common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wannier_centres in w90_input_setopt() library call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%centres = 0.0_dp

    allocate (common_data%wannier_data%spreads(common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wannier_spreads in w90_input_setopt() library call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%spreads = 0.0_dp

    ! read all other variables; mostly simple variables can be set and/or reset
    call w90_wannier90_readwrite_read(common_data%settings, common_data%band_plot, &
                                      common_data%dis_control, common_data%dis_spheres, &
                                      common_data%dis_manifold, common_data%fermi_energy_list, &
                                      common_data%fermi_surface_data, common_data%output_file, &
                                      common_data%wvfn_read, common_data%wann_control, &
                                      common_data%real_space_ham, common_data%kpoint_path, &
                                      common_data%w90_system, common_data%tran, &
                                      common_data%print_output, common_data%wann_plot, &
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

    ! clear settings container (from settings interface not .win file)
    deallocate (common_data%settings%entries, stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in deallocating entries data in w90_input_setopt() library call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine w90_input_setopt

  subroutine w90_input_reader(common_data, istdout, istderr, ierr)
    !! mechanism to act upon tokens read from the input ".win" file
    !! input is parsed and interpreted (any errors are identified) and
    !! the library data structure (variable common_data), already setup by w90_input_setopt, is modified

    ! w90_input_setopt() processes options stored in common_data%settings (must be empty here; emptied by w90_input_setopt call)
    ! w90_input_reader() processes options stored in common_data%in_data (from .win file)

    use w90_comms, only: valid_communicator
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_input, set_error_fatal, set_error_alloc, code_mpi
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_clean_infile
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    logical :: cp_pp

    ierr = 0

    if (.not. valid_communicator(common_data%comm)) then
      ! this is a problem: how do we exit using the parallel error handler when the communicator is unknown?
      write (istderr, *) ' Error: parallel Wannier90 library invoked with invalid communicator, exiting.  Use w90_set_comm()!'
      ierr = code_mpi
      return
    endif

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
                                      common_data%dis_manifold, common_data%fermi_energy_list, &
                                      common_data%fermi_surface_data, common_data%output_file, &
                                      common_data%wvfn_read, common_data%wann_control, &
                                      common_data%real_space_ham, common_data%kpoint_path, &
                                      common_data%w90_system, common_data%tran, &
                                      common_data%print_output, common_data%wann_plot, &
                                      common_data%ws_region, common_data%real_lattice, &
                                      common_data%w90_calculation, common_data%physics%bohr, &
                                      common_data%sitesym%symmetrize_eps, common_data%num_bands, &
                                      common_data%num_kpts, common_data%num_wann, &
                                      common_data%optimisation, common_data%calc_only_A, cp_pp, &
                                      common_data%gamma_only, common_data%lsitesymmetry, &
                                      common_data%use_bloch_phases, common_data%seedname, &
                                      istdout, error, common_data%comm)

    ! check
    if ((.not. common_data%kmesh_input%order_b_vectors) .and. &
            common_data%wann_control%use_ss_functional) then
            call set_error_input(error, 'Error: If use_ss_functional is true, &
                   order_b_vector must be true. ', common_data%comm)
      return 
    endif

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
    !! perform disentanglement; assumes library data object is already setup after w90_input_setopt() call
    !! no effect if number of bands == number of WF

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
      call set_error_fatal(error, 'Error: m_matrix_local not associated for w90_disentangle() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'Error: u_matrix not associated for w90_disentangle() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%u_matrix_opt)) then
      call set_error_fatal(error, 'Error: u_matrix_opt not associated for w90_disentangle() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%eigval)) then
      call set_error_fatal(error, 'Error: eigval not associated for w90_disentangle() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    ! if not already initialised, set disentanglement window to limits of spectrum
    if (common_data%dis_manifold%win_min == -huge(0.0_dp)) common_data%dis_manifold%win_min = minval(common_data%eigval)
    if (common_data%dis_manifold%win_max == huge(0.0_dp)) common_data%dis_manifold%win_max = maxval(common_data%eigval)
    if (common_data%dis_manifold%frozen_states) then
      if (common_data%dis_manifold%froz_min == -huge(0.0_dp)) then
        common_data%dis_manifold%froz_min = minval(common_data%eigval(:, :))
      endif
    endif

    ! condition for disentanglement is number of bands > number of WF
    if (common_data%num_bands > common_data%num_wann) then
      call dis_main(common_data%dis_control, common_data%dis_spheres, common_data%dis_manifold, &
                    common_data%kmesh_info, common_data%kpt_latt, common_data%sitesym, &
                    common_data%print_output, common_data%m_matrix_local, common_data%u_matrix, &
                    common_data%u_matrix_opt, common_data%eigval, common_data%real_lattice, &
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

  subroutine w90_project_overlap(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_overlap, only: overlap_project, overlap_project_gamma

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    integer :: ik, iw

    ierr = 0

    if (.not. associated(common_data%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for w90_project_overlap call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%u_matrix_opt)) then
      call set_error_fatal(error, 'u_matrix_opt not set for w90_project_overlap call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'u_matrixt not set for w90_project_overlap call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (.not. common_data%have_disentangled) then
      if (common_data%num_wann /= common_data%num_bands) then
        call set_error_fatal(error, 'Error: w90_project_overlap(): num_bands /= num_wann but disentanglement() was not called', &
                             common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif

      ! fixme, document!
      common_data%u_matrix(:, :, :) = common_data%u_matrix_opt(:, :, :) ! u_matrix_opt contains initial projections
      common_data%u_matrix_opt(:, :, :) = 0.d0
      do ik = 1, common_data%num_kpts
        do iw = 1, common_data%num_wann
          common_data%u_matrix_opt(iw, iw, ik) = 1.d0
        enddo
      enddo

      if (common_data%gamma_only) then
        call overlap_project_gamma(common_data%m_matrix_local, common_data%u_matrix, &
                                   common_data%kmesh_info%nntot, common_data%num_wann, &
                                   common_data%print_output%timing_level, istdout, &
                                   common_data%timer, error, common_data%comm)
      else
        call overlap_project(common_data%sitesym, common_data%m_matrix_local, common_data%u_matrix, &
                             common_data%kmesh_info%nnlist, common_data%kmesh_info%nntot, &
                             common_data%num_wann, common_data%num_kpts, common_data%num_wann, &
                             common_data%print_output%timing_level, common_data%lsitesymmetry, &
                             istdout, common_data%timer, common_data%dist_kpoints, error, &
                             common_data%comm)
      endif
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif
  end subroutine w90_project_overlap

  subroutine w90_wannierise(common_data, istdout, istderr, ierr)
    ! perform MLWF algorithm

    use w90_comms, only: mpirank, comms_sync_error
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
      call set_error_fatal(error, 'Error: m_matrix_local not set for call to w90_wannierise()', common_data%comm)
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'Error: u_matrix not set for w90_wannierise()', common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (common_data%gamma_only) then
      if (mpirank(common_data%comm) == 0) then
        call wann_main_gamma(common_data%kmesh_info, common_data%wann_control, common_data%omega, &
                             common_data%print_output, common_data%wannier_data, &
                             common_data%m_matrix_local, common_data%u_matrix, &
                             common_data%real_lattice, common_data%num_kpts, common_data%num_wann, &
                             istdout, common_data%timer, error, common_data%comm)
        if (allocated(error)) then
          call prterr(error, ierr, istdout, istderr, common_data%comm)
          return
        endif
      endif
      call comms_sync_error(common_data%comm, error, 0) ! this is necessary after root's excursion alone
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
                     common_data%optimisation, common_data%rpt_origin, common_data%band_plot%mode, &
                     common_data%tran%mode, common_data%lsitesymmetry, istdout, common_data%timer, &
                     common_data%dist_kpoints, error, common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine w90_wannierise

  subroutine w90_plot(common_data, istdout, istderr, ierr)
    !! performs a variety of plotting functions

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
                   common_data%u_matrix, common_data%u_matrix_opt, common_data%eigval, &
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
    !! performs a variety of transport calculations

    use w90_comms, only: mpirank, comms_sync_error
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
    if (mpirank(common_data%comm) == 0) then
      call tran_main(common_data%atom_data, common_data%dis_manifold, &
                     common_data%fermi_energy_list, common_data%ham_logical, common_data%kpt_latt, &
                     common_data%output_file, common_data%real_space_ham, common_data%tran, &
                     common_data%print_output, common_data%wannier_data, common_data%ws_region, &
                     common_data%w90_calculation, common_data%ham_k, common_data%ham_r, &
                     common_data%u_matrix, common_data%u_matrix_opt, common_data%eigval, &
                     common_data%real_lattice, common_data%wannier_centres_translated, &
                     common_data%irvec, common_data%mp_grid, common_data%ndegen, &
                     common_data%shift_vec, common_data%nrpts, common_data%num_bands, &
                     common_data%num_kpts, common_data%num_wann, common_data%rpt_origin, &
                     common_data%band_plot%mode, common_data%have_disentangled, &
                     common_data%lsitesymmetry, common_data%seedname, istdout, common_data%timer, &
                     error, common_data%comm)
    endif
    call comms_sync_error(common_data%comm, error, 0) ! this is necessary after root's excursion alone
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine w90_transport

  subroutine w90_set_m_local(common_data, m_matrix_local) ! m_matrix_local_orig
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: m_matrix_local(:, :, :, :)

    common_data%m_matrix_local => m_matrix_local
  end subroutine w90_set_m_local

  subroutine w90_set_u_matrix(common_data, u_matrix)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_matrix(:, :, :)

    common_data%u_matrix => u_matrix
  end subroutine w90_set_u_matrix

  subroutine w90_set_u_opt(common_data, u_matrix_opt)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_matrix_opt(:, :, :)

    common_data%u_matrix_opt => u_matrix_opt
  end subroutine w90_set_u_opt

  subroutine w90_set_eigval(common_data, eigval)
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    real(kind=dp), intent(in), target :: eigval(:, :)

    common_data%eigval => eigval
  end subroutine w90_set_eigval

  subroutine w90_set_constant_bohr_to_ang(common_data, bohr_to_angstrom)
    !! used to set the bohr_to_angstrom value as used in the SCF code
    implicit none

    type(lib_common_type), intent(inout) :: common_data
    real(kind=dp), intent(in) :: bohr_to_angstrom

    common_data%physics%bohr = bohr_to_angstrom
    common_data%physics%bohr_version_str = "-> Using Bohr value from linked main code"
  end subroutine w90_set_constant_bohr_to_ang

  subroutine w90_create_kmesh(common_data, istdout, istderr, ierr)
    !! causes w90 to calculate finite difference neighbour lists
    use w90_error_base, only: w90_error_type
    use w90_kmesh, only: kmesh_get, kmesh_sort

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0
    if (common_data%setup_complete) return

    call kmesh_get(common_data%kmesh_input, common_data%kmesh_info, common_data%print_output, &
                   common_data%kpt_latt, common_data%real_lattice, common_data%num_kpts, &
                   common_data%gamma_only, istdout, common_data%timer, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (.not. common_data%gamma_only) then
      if (common_data%kmesh_input%order_b_vectors) then
        call kmesh_sort(common_data%kmesh_info, common_data%num_kpts, error, common_data%comm)
      endif
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif

    common_data%setup_complete = .true.
  end subroutine w90_create_kmesh

  subroutine w90_get_nn(common_data, nn, istdout, istderr, ierr)
    !! probe w90 library for number of finite difference k-point neighbours
    implicit none

    ! arguments
    integer, intent(out) :: nn, ierr
    !! nn is the number of neighbours in F.D. scheme
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data
    !! library data object

    if (.not. common_data%setup_complete) then
      call w90_create_kmesh(common_data, istdout, istderr, ierr)
      !! setup k-mesh (b vectors) if not already done (sets setup_complete)
      if (ierr > 0) return
    endif

    nn = common_data%kmesh_info%nntot
  end subroutine w90_get_nn

  subroutine w90_get_nnkp(common_data, nnkp, istdout, istderr, ierr)
    !! probe w90 library for finite difference k-point neighbour indices
    implicit none

    integer, intent(out) :: nnkp(:, :), ierr
    !! nnkp must be dimensioned (n_neighbours, n_fbz)
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data
    !! library data object

    if (.not. common_data%setup_complete) then
      call w90_create_kmesh(common_data, istdout, istderr, ierr)
      !! setup k-mesh (b vectors) if not already done (sets setup_complete)
      if (ierr > 0) return
    endif

    nnkp = common_data%kmesh_info%nnlist
  end subroutine w90_get_nnkp

  subroutine w90_get_gkpb(common_data, gkpb, istdout, istderr, ierr)
    !! probe w90 library for the triple of reciprocal lattice translations determining phase in k'= k+b
    implicit none

    integer, intent(out) :: gkpb(:, :, :), ierr
    !! gkpb must be dimensioned (3,nk,nnb)
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data
    !! library data object

    if (.not. common_data%setup_complete) then
      call w90_create_kmesh(common_data, istdout, istderr, ierr)
      !! setup k-mesh (b vectors) if not already done (sets setup_complete)
      if (ierr > 0) return
    endif

    gkpb = common_data%kmesh_info%nncell
  end subroutine w90_get_gkpb

  subroutine w90_get_centres(common_data, centres)
    !! probes w90 library for (current) wannier centres
    implicit none

    real(kind=dp), intent(out) :: centres(:, :)
    !! must be allocated with size >= n_wannier
    type(lib_common_type), intent(in) :: common_data
    !! library data object

    centres = common_data%wannier_data%centres
  endsubroutine w90_get_centres

  subroutine w90_get_spreads(common_data, spreads)
    !! probes w90 library for (current) wannier spreads
    implicit none

    real(kind=dp), intent(out) :: spreads(:)
    !! must be allocated with size >= n_wannier
    type(lib_common_type), intent(in) :: common_data
    !! library data object

    spreads = common_data%wannier_data%spreads
  endsubroutine w90_get_spreads

  subroutine w90_get_proj(common_data, n, site, l, m, s, rad, x, z, sqa, zona, istdout, istderr, ierr)
    !! probes library data object and returns arrays describing a list of projections
    !! projectors defined either in .win file or passed (using same syntax) through w90_setopt
    !! array arguments assumed allocated at call
    !! array arguments must have length >= number of projectors (checked here)
    use w90_error, only: w90_error_type, set_error_fatal
    implicit none
    integer, intent(in) :: istdout, istderr
    integer, intent(inout) :: n, l(:), m(:), s(:)
    !! number of projectors, angular, orbital and spin numbers
    integer, intent(inout) :: rad(:)
    !! radial function defining projector
    integer, intent(out) :: ierr
    !! ierr returned > 0 in case of error
    real(kind=dp), intent(inout) :: site(:, :)
    !! projector origin (kind, site of kind)
    real(kind=dp), intent(inout) :: sqa(:, :), z(:, :), x(:, :), zona(:)
    !! spin quantisation axis, z- an x-axes
    type(lib_common_type), intent(in), target :: common_data
    !! library data object

    ! local variables
    integer :: ip
    type(proj_type), pointer :: proj
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. allocated(common_data%proj_input)) then
      call set_error_fatal(error, &
                           'Error: projectors are not setup in Wannier90 library when requested via get_proj()', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    n = size(common_data%proj_input)

    ! check allocation of main output arrays
    if (size(l) < n) then
      call set_error_fatal(error, 'Error: array argument l in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (size(m) < n) then
      call set_error_fatal(error, 'Error: array argument m in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (size(s) < n) then
      call set_error_fatal(error, 'Error: array argument s in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (size(site, 2) < n) then
      call set_error_fatal(error, 'Error: array argument site in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (size(sqa, 2) < n) then
      call set_error_fatal(error, 'Error: array argument sqa in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (size(sqa, 2) < n) then
      call set_error_fatal(error, 'Error: array argument sqa in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (size(z, 2) < n) then
      call set_error_fatal(error, 'Error: array argument z in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    elseif (size(x, 2) < n) then
      call set_error_fatal(error, 'Error: array argument x in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (size(rad) < n) then
      call set_error_fatal(error, 'Error: array argument rad in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (size(zona) < n) then
      call set_error_fatal(error, 'Error: array argument zona in get_proj() call is insufficiently sized', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    do ip = 1, n
      proj => common_data%proj_input(ip)
      l(ip) = proj%l
      m(ip) = proj%m
      s(ip) = proj%s
      site(1:3, ip) = proj%site(1:3)
      sqa(1:3, ip) = proj%s_qaxis(1:3)
      z(1:3, ip) = proj%z(1:3)
      x(1:3, ip) = proj%x(1:3)
      rad(ip) = proj%radial
      zona(ip) = proj%zona
    enddo
  end subroutine w90_get_proj

  subroutine w90_set_comm(common_data, comm)
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
  end subroutine w90_set_comm

  subroutine w90_print_info(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_readwrite, only: w90_readwrite_write_header
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write
    use w90_comms, only: mpisize, mpirank

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    integer :: mpi_size

    ierr = 0

    if (mpirank(common_data%comm) == 0) then
      mpi_size = mpisize(common_data%comm)
      ! write jazzy header info
      call w90_readwrite_write_header(common_data%physics%bohr_version_str, &
                                      common_data%physics%constants_version_str1, &
                                      common_data%physics%constants_version_str2, &
                                      mpi_size, istdout)
    endif

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
                                       common_data%wannier_data, common_data%wann_plot, &
                                       common_data%w90_calculation, common_data%real_lattice, &
                                       common_data%sitesym%symmetrize_eps, common_data%mp_grid, &
                                       common_data%num_bands, common_data%num_kpts, &
                                       common_data%num_proj, common_data%num_wann, &
                                       common_data%optimisation, .false., common_data%gamma_only, &
                                       common_data%lsitesymmetry, common_data%w90_system%spinors, &
                                       common_data%use_bloch_phases, istdout)

    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine w90_print_info

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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
  endsubroutine w90_set_option_r2d

  subroutine w90_set_option_c2d(common_data, keyword, arr)
    use w90_readwrite, only: init_settings, expand_settings

    implicit none

    character(*), intent(in) :: keyword
    character(len=*), intent(in) :: arr(:)
    type(lib_common_type), intent(inout) :: common_data
    integer :: i

    if (.not. allocated(common_data%settings%entries)) call init_settings(common_data%settings)
    i = common_data%settings%num_entries + 1
    common_data%settings%entries(i)%keyword = keyword
    common_data%settings%entries(i)%c2d = arr
    common_data%settings%num_entries = i + 1
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
  endsubroutine w90_set_option_c2d

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
    if (common_data%settings%num_entries == common_data%settings%num_entries_max) then
      call expand_settings(common_data%settings)
    endif
  endsubroutine w90_set_option_real

  subroutine w90_distribute_kpts(common_data, num_kpts, mpi_size, dist_k, istdout, istderr, ierr)
    !! provide a distribution of num_kpts k-points across mpi_size MPI ranks
    ! should be called from all ranks in a parallel environment for error propagation
    use w90_comms, only: comms_sync_error
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal

    implicit none

    ! arguments
    integer, intent(in) :: num_kpts
    !! number of k-points
    integer, intent(in) :: mpi_size
    !! number of ranks in MPI communicator
    integer, intent(in) :: istdout, istderr
    !! destination for error messages
    integer, intent(inout), allocatable :: dist_k(:)
    !! already allocated array
    !! assigned here such that dist_k(i) = rank handling kpt i
    !! size and allocation status are tested
    integer, intent(out) :: ierr
    !! return code, nonzero in case of error
    type(lib_common_type), intent(in) :: common_data
    !! library object: only the communicator type is referenced

    ! local variables
    type(w90_error_type), allocatable :: error
    integer :: ctr, i, nkl

    ierr = 0

    if (mpi_size < 1) then
      call set_error_fatal(error, 'Error: mpi_size < 1 in w90_distribute_kpts call.', common_data%comm)
    elseif (num_kpts < 1) then
      call set_error_fatal(error, 'Error: num_kpts < 1 in w90_distribute_kpts call.', common_data%comm)
    elseif (.not. allocated(dist_k)) then
      call set_error_fatal(error, 'Error: dist_k not allocated in w90_distribute_kpts call.', common_data%comm)
    elseif (size(dist_k) < num_kpts) then
      call set_error_fatal(error, 'Error: size(dist_k) < num_kpts in w90_distribute_kpts call.', common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    ctr = 0
    do i = 0, mpi_size - 1
      nkl = num_kpts/mpi_size ! number of kpoints per rank
      if (mod(num_kpts, mpi_size) > i) nkl = nkl + 1
      if (nkl > 0) then
        dist_k(ctr + 1:ctr + nkl) = i
        ctr = ctr + nkl
      endif
    enddo
  end subroutine w90_distribute_kpts

end module w90_library
