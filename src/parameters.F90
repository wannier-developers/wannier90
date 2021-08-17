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

module w90_param_types
  !! This module contains parameters to control the actions of wannier90.
  !! Also routines to read the parameters and write them out again.

  use w90_constants, only: dp
  use w90_io, only: maxlen

  implicit none

  public

  type print_output_type
    !! ==============================
    !! Contains variables to control output file formatting and verbosity.
    !! ==============================
    ! verbosity flags - param_read_verbosity
    integer :: iprint
    ! Controls the verbosity of the output
    integer :: timing_level
    ! REVIEW_2021-07-22: optimisation doesn't fit in this type (TO FINISH)
    integer :: optimisation !wannierise and disentangle
    ! REVIEW_2021-07-22: we agree that we don't need both length_unit and lenconfac;
    ! REVIEW_2021-07-22: instead could have a utility function.
    character(len=20) :: length_unit ! MAYBE, just have a separate variable?
    ! Units for length
    real(kind=dp) :: lenconfac !lots of write statements in wannier90
  end type print_output_type

  type w90_system_type
    !! ==============================
    !! Contains physical information about the material being calculated.
    !! ==============================
    integer :: num_valence_bands !wannierise, postw90/postw90_common, get_oper and berry
    integer :: num_elec_per_state !wannierise and postw90 dos and boltzwann
    logical :: spinors   !are our WF spinors? !kmesh, plot, wannier_lib, postw90/gyrotropic
  end type w90_system_type

  type ws_region_type
    logical :: use_ws_distance !ws_distance, plot and postw90_common
    real(kind=dp) :: ws_distance_tol !ws_distance, hamiltonian and postw90_common
    !! absolute tolerance for the distance to equivalent positions
    integer :: ws_search_size(3) ! ws_distance, hamiltonian
    !! maximum extension in each direction of the supercell of the BvK cell
    !! to search for points inside the Wigner-Seitz cell
  end type ws_region_type

  ! setup in wannierise, but used by plot, ws_distance etc
  type wannier_data_type
    !! =========================================
    !! Contains the centres and spreads of the MLWFs
    !! =========================================
    ! Wannier centres and spreads
    real(kind=dp), allocatable :: centres(:, :)
    real(kind=dp), allocatable :: spreads(:)
    ! REVIEW_2021-07-22: Do we want to expose other related variables such as the decomposition
    ! REVIEW_2021-07-22: of the spread, matrix elements of r and r^2, etc. (TO FINISH)
  end type wannier_data_type

  ! used in kmesh, and to allocate in parameters
  ! The maximum number of shells we need to satisfy B1 condition in kmesh
  integer, parameter :: max_shells = 6
  integer, parameter :: num_nnmax = 12

  type kmesh_input_type
    !! ================================
    !! Contains information that can be provided by the user about determining the kmesh
    !! ================================
    integer :: num_shells
    !! no longer an input keyword
    logical :: skip_B1_tests
    !! do not check the B1 condition
    integer, allocatable :: shell_list(:)
    integer :: search_shells
    real(kind=dp) :: tol
  end type kmesh_input_type

  !AAM: There are a number of ways one can handle the initial guess. (i) specify explicit
  !AAM: projections; (ii) use random (s-orbital) projections; (iii) a combination of (i) and
  !AAM: (ii); (iv) use the phases from the Bloch functions directly; (v) SCDM method. (i), (ii)
  !AAM: and (iii) require the arrays defined in the "projection_type" below. (iv) and (v) do not.
  !AAM: (vi) An external code may also simply supply to w90 an Amn(k) matrix that is has independently
  !AAM: generated, in which case projection_type is not needed.
  !AAM: It makes sense to keep the projection sites separate from the projection_type data below.
  type proj_input_type
    !! ========================
    !! Contains information that can be provided by the user about the projections
    !! ========================
    ! Projections, mainly used in parameters, maybe written by kmesh
    ! REVIEW_2021-07-22: site(:,:) has dual usage: for projections and for guiding centres.
    ! REVIEW_2021-07-22: Make a new type for guiding centres that only contains sites.
    ! REVIEW_2021-07-22: In the future this can be logically distinct from the projection sites.
    ! REVIEW_2021-07-22: For now, when defining proj_input_type, also define sites inside the
    ! REVIEW_2021-07-22: new guiding centres type.
    real(kind=dp), allocatable :: site(:, :)
    integer, allocatable :: l(:)
    integer, allocatable :: m(:)
    integer, allocatable :: s(:)
    real(kind=dp), allocatable :: s_qaxis(:, :)
    real(kind=dp), allocatable :: z(:, :)
    real(kind=dp), allocatable :: x(:, :)
    integer, allocatable :: radial(:)
    real(kind=dp), allocatable :: zona(:)
    ! a u t o m a t i c   p r o j e c t i o n s
    ! vv: Writes a new block in .nnkp
    logical :: auto_projections
  end type proj_input_type

  ! kmesh parameters (set in kmesh)
  type kmesh_info_type
    !! =======================
    !! Contains derived information about the kmesh
    !! =======================
    integer              :: nnh           ! the number of b-directions (bka)
    integer              :: nntot         ! total number of neighbours for each k-point
    integer, allocatable :: nnlist(:, :)   ! list of neighbours for each k-point
    integer, allocatable :: neigh(:, :)
    integer, allocatable :: nncell(:, :, :) ! gives BZ of each neighbour of each k-point
    real(kind=dp)              :: wbtot
    real(kind=dp), allocatable :: wb(:)         ! weights associated with neighbours of each k-point
    real(kind=dp), allocatable :: bk(:, :, :)     ! the b-vectors that go from each k-point to its neighbours
    real(kind=dp), allocatable :: bka(:, :)      ! the b-directions from 1st k-point to its neighbours
    logical :: explicit_nnkpts
    !! nnkpts block is in the input file (allowed only for post-proc setup)
  end type kmesh_info_type

  ! used in wannierise, hamiltonian, plot and others (postw90 also)
  type k_points_type
    !! =====================
    !! Contains information about the kpoints used in the calculation.
    !! =====================
    real(kind=dp), allocatable :: kpt_latt(:, :) !! kpoints in lattice vecs
    ! REVIEW_2021-07-22: we can generate kpt_cart from kpt_latt as and when
    ! REVIEW_2021-07-22: we need it (usage is very localised in the code).
    ! REVIEW_2021-07-22: We have a utility that does the conversion already.
    ! REVIEW_2021-07-22: Then it doesn't make sense to have a type for just kpt_latt.
    real(kind=dp), allocatable :: kpt_cart(:, :) !kpoints in cartesians - kmesh and transport
  end type k_points_type

  ! this contains data which described the disentangled manifold, also used in postw90
  type dis_manifold_type
    !! ===========================
    !! Contains information about the manifold of states from which the MLWFs are to be disentangled.
    !! ===========================
    real(kind=dp) :: win_min
    !! lower bound of the disentanglement outer window
    real(kind=dp) :: win_max
    !! upper bound of the disentanglement outer window
    real(kind=dp) :: froz_min
    !! lower bound of the disentanglement inner (frozen) window
    real(kind=dp) :: froz_max
    !! upper bound of the disentanglement inner (frozen) window
    logical :: frozen_states
    ! disentangle parameters
    ! Used by plot, hamiltonian, wannierise, postw90_common, get_oper - not read
    integer, allocatable :: ndimwin(:)
    logical, allocatable :: lwindow(:, :)
  end type dis_manifold_type

  ! Atom sites - often used in the write_* routines
  ! hamiltonian, wannierise, plot, transport, wannier_lib
  type atom_data_type
    !! =======================
    !! Contains information about the atoms (and maybe the cell...) of the system being calculated.
    !! =======================
    ! REVIEW_2021-07-22: Keep only one of pos_frac and pos_cart to avoid possible inconsistency
    ! REVIEW_2021-07-22: on passing these variables in from the external code. We can generate
    ! REVIEW_2021-07-22: the other internally using a utility function.
    ! REVIEW_2021-07-22: Shall we include unit cell data here too? Good arguments for and
    ! REVIEW_2021-07-22: against -- if most subroutines only need cell and not atomic data,
    ! REVIEW_2021-07-22: then keep separate. In any case, also only pass one of
    ! REVIEW_2021-07-22: real_cell / recip_cell and generate the other internally.
    real(kind=dp), allocatable :: pos_frac(:, :, :)
    real(kind=dp), allocatable :: pos_cart(:, :, :)
    integer, allocatable :: species_num(:)
    character(len=maxlen), allocatable :: label(:)
    character(len=2), allocatable :: symbol(:)
    integer :: num_atoms
    integer :: num_species
  end type atom_data_type

  ! plot.F90 and postw90/kpath
  type kpoint_path_type
    !! ============================
    !! Contains information that specifies the k-point path for plotting and other purposes.
    !! Note: The length of bands_label and the second index of bands_spec_points is twice the
    !! number of segments specified by the user. Each pair of special points defines a segment.
    !! ============================
    integer num_points_first_segment
    character(len=20), allocatable :: labels(:)
    real(kind=dp), allocatable :: points(:, :)
  end type kpoint_path_type

end module w90_param_types

module w90_param_methods
  ! very few of these use save, so may actually be local to subroutines

  use w90_constants, only: dp
  use w90_io, only: maxlen
  use w90_param_types

  implicit none

  private

  ! Private data for processing input file
  integer                            :: num_lines
  character(len=maxlen), allocatable :: in_data(:)

  public :: param_dealloc
  public :: param_write_header
  public :: param_read_chkpt
  public :: param_lib_set_atoms
  public :: param_get_smearing_type
  public :: param_get_convention_type
  public :: param_chkpt_dist
  ! for postw90 parameters
  public :: param_in_file
  public :: param_get_keyword
  public :: param_get_keyword_block
  public :: param_get_block_length
  public :: param_get_vector_length
  public :: param_get_keyword_vector
  public :: param_get_range_vector
  public :: param_get_centre_constraints
  public :: param_get_projections
  public :: get_smearing_index
  public :: internal_set_kmesh
  public :: param_uppercase
  ! common read routines
  public :: param_read_gamma_only
  public :: param_read_verbosity
  public :: param_read_num_wann
  public :: param_read_exclude_bands
  public :: param_read_lattice
  public :: param_read_atoms
  public :: param_read_units
  public :: param_read_devel
  public :: param_read_mp_grid
  public :: param_read_kpath
  public :: param_read_dis_manifold
  public :: param_read_num_bands
  public :: param_read_system
  public :: param_read_fermi_energy
  public :: param_read_ws_data
  public :: param_read_eigvals
  public :: param_read_kmesh_data
  public :: param_read_kpoints
  public :: param_read_global_kmesh
  public :: param_clean_infile
  public :: param_read_final_alloc
  public :: get_all_keywords
  private :: param_clear_block

contains

  subroutine param_read_verbosity(verbose, stdout, seedname)
    !%%%%%%%%%%%%%%%%
    !System variables
    !%%%%%%%%%%%%%%%%
    implicit none
    type(print_output_type), intent(inout) :: verbose
    logical :: found
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    verbose%timing_level = 1             ! Verbosity of timing output info
    call param_get_keyword(stdout, seedname, 'timing_level', found, i_value=verbose%timing_level)

    verbose%iprint = 1             ! Verbosity
    call param_get_keyword(stdout, seedname, 'iprint', found, i_value=verbose%iprint)

    verbose%optimisation = 3             ! Verbosity
    call param_get_keyword(stdout, seedname, 'optimisation', found, i_value=verbose%optimisation)

  end subroutine param_read_verbosity

  subroutine param_read_units(lenconfac, length_unit, energy_unit, bohr, stdout, seedname)
    !use w90_constants, only: bohr
    use w90_io, only: io_error
    implicit none
    real(kind=dp), intent(out) :: lenconfac
    integer, intent(in) :: stdout
    character(len=*), intent(out) :: length_unit
    character(len=*), intent(out) :: energy_unit
    character(len=50), intent(in)  :: seedname
    real(kind=dp), intent(in) :: bohr
    logical :: found

    energy_unit = 'ev'          !
    call param_get_keyword(stdout, seedname, 'energy_unit', found, c_value=energy_unit)

    length_unit = 'ang'         !
    lenconfac = 1.0_dp
    call param_get_keyword(stdout, seedname, 'length_unit', found, c_value=length_unit)
    if (length_unit .ne. 'ang' .and. length_unit .ne. 'bohr') &
      call io_error('Error: value of length_unit not recognised in param_read', stdout, seedname)
    if (length_unit .eq. 'bohr') lenconfac = 1.0_dp/bohr
  end subroutine param_read_units

  subroutine param_read_num_wann(num_wann, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    integer, intent(out) :: num_wann
    character(len=50), intent(in)  :: seedname

    logical :: found

    num_wann = -99
    call param_get_keyword(stdout, seedname, 'num_wann', found, i_value=num_wann)
    if (.not. found) call io_error('Error: You must specify num_wann', stdout, seedname)
    if (num_wann <= 0) call io_error('Error: num_wann must be greater than zero', stdout, seedname)
  end subroutine param_read_num_wann

  subroutine param_read_exclude_bands(exclude_bands, num_exclude_bands, stdout, seedname)
    use w90_io, only: io_error
    implicit none

    integer, allocatable, intent(inout) :: exclude_bands(:)
    integer, intent(out) :: num_exclude_bands
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: ierr
    logical :: found

    num_exclude_bands = 0
    call param_get_range_vector(stdout, seedname, 'exclude_bands', found, &
                                num_exclude_bands, lcount=.true.)
    if (found) then
      if (num_exclude_bands < 1) call io_error('Error: problem reading exclude_bands', stdout, seedname)
      if (allocated(exclude_bands)) deallocate (exclude_bands)
      allocate (exclude_bands(num_exclude_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating exclude_bands in param_read', stdout, seedname)
      call param_get_range_vector(stdout, seedname, 'exclude_bands', found, &
                                  num_exclude_bands, .false., exclude_bands)
      if (any(exclude_bands < 1)) &
        call io_error('Error: exclude_bands must contain positive numbers', stdout, seedname)
    end if
  end subroutine param_read_exclude_bands

  subroutine param_read_num_bands(pw90_effective_model, library, num_exclude_bands, num_bands, &
                                  num_wann, library_param_read_first_pass, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: pw90_effective_model, library
    integer, intent(in) :: num_exclude_bands
    integer, intent(inout) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    logical, intent(in) :: library_param_read_first_pass
    character(len=50), intent(in)  :: seedname

    integer :: i_temp
    logical :: found

    ! AAM_2016-09-16: some changes to logic to patch a problem with uninitialised num_bands in library mode
!    num_bands       =   -1
    call param_get_keyword(stdout, seedname, 'num_bands', found, i_value=i_temp)
    if (found .and. library) write (stdout, '(/a)') ' Ignoring <num_bands> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      if (found) num_bands = i_temp
      if (.not. found) num_bands = num_wann
    end if
    ! GP: I subtract it here, but only the first time when I pass the total number of bands
    ! In later calls, I need to pass instead num_bands already subtracted.
    if (library .and. library_param_read_first_pass) num_bands = num_bands - num_exclude_bands
    if (.not. pw90_effective_model) then
      if (found .and. num_bands < num_wann) then
        write (stdout, *) 'num_bands', num_bands
        write (stdout, *) 'num_wann', num_wann
        call io_error('Error: num_bands must be greater than or equal to num_wann', stdout, seedname)
      endif
    endif
  end subroutine param_read_num_bands

  subroutine param_read_devel(devel_flag, stdout, seedname)
!   use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    character(len=*), intent(out) :: devel_flag
    character(len=50), intent(in)  :: seedname

    logical :: found

    devel_flag = ' '          !
    call param_get_keyword(stdout, seedname, 'devel_flag', found, c_value=devel_flag)
  end subroutine param_read_devel

  subroutine param_read_gamma_only(gamma_only, num_kpts, library, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(inout) :: gamma_only
    integer, intent(in) :: num_kpts
    logical, intent(in) :: library
    character(len=50), intent(in)  :: seedname

    logical :: found, ltmp

    ltmp = .false.
    call param_get_keyword(stdout, seedname, 'gamma_only', found, l_value=ltmp)
    if (.not. library) then
      gamma_only = ltmp
      if (gamma_only .and. (num_kpts .ne. 1)) &
        call io_error('Error: gamma_only is true, but num_kpts > 1', stdout, seedname)
    else
      if (found) write (stdout, '(a)') ' Ignoring <gamma_only> in input file'
    endif
  end subroutine param_read_gamma_only

  subroutine param_read_mp_grid(pw90_effective_model, library, mp_grid, num_kpts, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: pw90_effective_model, library
    integer, intent(inout) :: mp_grid(3), num_kpts
    character(len=50), intent(in)  :: seedname

    integer :: iv_temp(3)
    logical :: found

!    mp_grid=-99
    call param_get_keyword_vector(stdout, seedname, 'mp_grid', found, 3, i_value=iv_temp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <mp_grid> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      if (found) mp_grid = iv_temp
      if (.not. found) then
        call io_error('Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid', stdout, seedname)
      elseif (any(mp_grid < 1)) then
        call io_error('Error: mp_grid must be greater than zero', stdout, seedname)
      end if
      num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
    end if
  end subroutine param_read_mp_grid

  subroutine param_read_system(library, system, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: library
    type(w90_system_type), intent(inout) :: system
    character(len=50), intent(in)  :: seedname

    logical :: found, ltmp

    ltmp = .false.  ! by default our WF are not spinors
    call param_get_keyword(stdout, seedname, 'spinors', found, l_value=ltmp)
    if (.not. library) then
      system%spinors = ltmp
    else
      if (found) write (stdout, '(a)') ' Ignoring <spinors> in input file'
    endif
!    if(spinors .and. (2*(num_wann/2))/=num_wann) &
!       call io_error('Error: For spinor WF num_wann must be even')

    ! We need to know if the bands are double degenerate due to spin, e.g. when
    ! calculating the DOS
    if (system%spinors) then
      system%num_elec_per_state = 1
    else
      system%num_elec_per_state = 2
    endif
    call param_get_keyword(stdout, seedname, 'num_elec_per_state', found, &
                           i_value=system%num_elec_per_state)
    if ((system%num_elec_per_state /= 1) .and. (system%num_elec_per_state /= 2)) &
      call io_error('Error: num_elec_per_state can be only 1 or 2', stdout, seedname)
    if (system%spinors .and. system%num_elec_per_state /= 1) &
      call io_error('Error: when spinors = T num_elec_per_state must be 1', stdout, seedname)

    ! set to a negative default value
    system%num_valence_bands = -99
    call param_get_keyword(stdout, seedname, 'num_valence_bands', found, i_value=system%num_valence_bands)
    if (found .and. (system%num_valence_bands .le. 0)) &
      call io_error('Error: num_valence_bands should be greater than zero', stdout, seedname)
    ! there is a check on this parameter later

  end subroutine param_read_system

  subroutine param_read_kpath(library, spec_points, ok, bands_plot, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: library, bands_plot
    type(kpoint_path_type), intent(out) :: spec_points
    integer, intent(in) :: stdout
    logical, intent(out) :: ok
    character(len=50), intent(in)  :: seedname

    integer :: i_temp, ierr, bands_num_spec_points
    logical :: found

    bands_num_spec_points = 0
    call param_get_block_length(stdout, seedname, 'kpoint_path', found, i_temp, library)
    if (found) then
      ok = .true.
      bands_num_spec_points = i_temp*2
      if (allocated(spec_points%labels)) deallocate (spec_points%labels)
      allocate (spec_points%labels(bands_num_spec_points), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating labels in param_read', stdout, seedname)
      if (allocated(spec_points%points)) deallocate (spec_points%points)
      allocate (spec_points%points(3, bands_num_spec_points), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating points in param_read', stdout, seedname)
      call param_get_keyword_kpath(spec_points, stdout, seedname)
    else
      ok = .false.
    end if
    spec_points%num_points_first_segment = 100
    call param_get_keyword(stdout, seedname, 'bands_num_points', found, &
                           i_value=spec_points%num_points_first_segment)
    ! checks
    if (bands_plot) then
      if (spec_points%num_points_first_segment < 0) &
        call io_error('Error: bands_num_points must be positive', stdout, seedname)
    endif
  end subroutine param_read_kpath

  subroutine param_read_fermi_energy(found_fermi_energy, fermi_energy_list, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    logical, intent(out) :: found_fermi_energy
    real(kind=dp), allocatable, intent(out) :: fermi_energy_list(:)
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    real(kind=dp) :: fermi_energy
    logical :: fermi_energy_scan
    real(kind=dp) :: fermi_energy_min
    real(kind=dp) :: fermi_energy_max
    real(kind=dp) :: fermi_energy_step
    integer :: i, ierr, n
    logical :: found

    n = 0
    found_fermi_energy = .false.
    call param_get_keyword(stdout, seedname, 'fermi_energy', found, r_value=fermi_energy)
    if (found) then
      found_fermi_energy = .true.
      n = 1
    endif
    !
    fermi_energy_scan = .false.
    call param_get_keyword(stdout, seedname, 'fermi_energy_min', found, r_value=fermi_energy_min)
    if (found) then
      if (found_fermi_energy) call io_error( &
        'Error: Cannot specify both fermi_energy and fermi_energy_min', stdout, seedname)
      fermi_energy_scan = .true.
      fermi_energy_max = fermi_energy_min + 1.0_dp
      call param_get_keyword(stdout, seedname, 'fermi_energy_max', found, &
                             r_value=fermi_energy_max)
      if (found .and. fermi_energy_max <= fermi_energy_min) call io_error( &
        'Error: fermi_energy_max must be larger than fermi_energy_min', stdout, seedname)
      fermi_energy_step = 0.01_dp
      call param_get_keyword(stdout, seedname, 'fermi_energy_step', found, &
                             r_value=fermi_energy_step)
      if (found .and. fermi_energy_step <= 0.0_dp) call io_error( &
        'Error: fermi_energy_step must be positive', stdout, seedname)
      n = nint(abs((fermi_energy_max - fermi_energy_min)/fermi_energy_step)) + 1
    endif
    !
    if (found_fermi_energy) then
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(1), stat=ierr)
      fermi_energy_list(1) = fermi_energy
    elseif (fermi_energy_scan) then
      if (n .eq. 1) then
        fermi_energy_step = 0.0_dp
      else
        fermi_energy_step = (fermi_energy_max - fermi_energy_min)/real(n - 1, dp)
      endif
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(n), stat=ierr)
      do i = 1, n
        fermi_energy_list(i) = fermi_energy_min + (i - 1)*fermi_energy_step
      enddo
!!    elseif(nfermi==0) then
!!        ! This happens when both found_fermi_energy=.false. and
!!        ! fermi_energy_scan=.false. Functionalities that require
!!        ! specifying a Fermi level should give an error message
!!        allocate(fermi_energy_list(1),stat=ierr) ! helps streamline things
!!
!! AAM_2017-03-27: if nfermi is zero (ie, fermi_energy* parameters are not set in input file)
!! then allocate fermi_energy_list with length 1 and set to zero as default.
    else
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(1), stat=ierr)
      fermi_energy_list(1) = 0.0_dp
    endif
    if (ierr /= 0) call io_error( &
      'Error allocating fermi_energy_list in param_read', stdout, seedname)
  end subroutine param_read_fermi_energy

  subroutine param_read_ws_data(param_input, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    type(ws_region_type), intent(inout) :: param_input
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: i
    logical :: found

    param_input%use_ws_distance = .true.
    call param_get_keyword(stdout, seedname, 'use_ws_distance', found, l_value=param_input%use_ws_distance)

    param_input%ws_distance_tol = 1.e-5_dp
    call param_get_keyword(stdout, seedname, 'ws_distance_tol', found, r_value=param_input%ws_distance_tol)

    param_input%ws_search_size = 2

    call param_get_vector_length(stdout, seedname, 'ws_search_size', found, length=i)
    if (found) then
      if (i .eq. 1) then
        call param_get_keyword_vector(stdout, seedname, 'ws_search_size', found, 1, &
                                      i_value=param_input%ws_search_size)
        param_input%ws_search_size(2) = param_input%ws_search_size(1)
        param_input%ws_search_size(3) = param_input%ws_search_size(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector(stdout, seedname, 'ws_search_size', found, 3, &
                                      i_value=param_input%ws_search_size)
      else
        call io_error('Error: ws_search_size must be provided as either one integer or a vector of three integers', &
                      stdout, seedname)
      end if
      if (any(param_input%ws_search_size <= 0)) &
        call io_error('Error: ws_search_size elements must be greater than zero', stdout, seedname)
    end if
  end subroutine param_read_ws_data

  subroutine param_read_eigvals(pw90_effective_model, pw90_boltzwann, pw90_geninterp, w90_plot, &
                                disentanglement, eig_found, eigval, library, postproc_setup, &
                                num_bands, num_kpts, stdout, seedname)

    use w90_io, only: io_file_unit, io_error

    implicit none
    logical, intent(in) :: pw90_effective_model, pw90_boltzwann, &
                           pw90_geninterp, w90_plot, disentanglement, &
                           library, postproc_setup
    logical, intent(out) :: eig_found
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    integer, intent(in) :: num_bands, num_kpts
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    !w90_plot = bands_plot or dos_plot or fermi_surface_plot or write_hr
    integer :: i, j, k, n, eig_unit, ierr
    !integer, allocatable, dimension(:, :) :: nnkpts_block
    !integer, allocatable, dimension(:) :: nnkpts_idx

    ! Read the eigenvalues from wannier.eig
    eig_found = .false.
    if (.not. library .and. .not. pw90_effective_model) then

      if (.not. postproc_setup) then
        inquire (file=trim(seedname)//'.eig', exist=eig_found)
        if (.not. eig_found) then
          if (disentanglement) then
            call io_error('No '//trim(seedname)//'.eig file found. Needed for disentanglement', stdout, seedname)
          else if ((w90_plot .or. pw90_boltzwann .or. pw90_geninterp)) then
            call io_error('No '//trim(seedname)//'.eig file found. Needed for interpolation', stdout, seedname)
          end if
        else
          ! Allocate only here
          allocate (eigval(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating eigval in param_read', stdout, seedname)

          eig_unit = io_file_unit()
          open (unit=eig_unit, file=trim(seedname)//'.eig', form='formatted', status='old', err=105)
          do k = 1, num_kpts
            do n = 1, num_bands
              read (eig_unit, *, err=106, end=106) i, j, eigval(n, k)
              if ((i .ne. n) .or. (j .ne. k)) then
                write (stdout, '(a)') 'Found a mismatch in '//trim(seedname)//'.eig'
                write (stdout, '(a,i0,a,i0)') 'Wanted band  : ', n, ' found band  : ', i
                write (stdout, '(a,i0,a,i0)') 'Wanted kpoint: ', k, ' found kpoint: ', j
                write (stdout, '(a)') ' '
                write (stdout, '(a)') 'A common cause of this error is using the wrong'
                write (stdout, '(a)') 'number of bands. Check your input files.'
                write (stdout, '(a)') 'If your pseudopotentials have shallow core states remember'
                write (stdout, '(a)') 'to account for these electrons.'
                write (stdout, '(a)') ' '
                call io_error('param_read: mismatch in '//trim(seedname)//'.eig', stdout, seedname)
              end if
            enddo
          end do
          close (eig_unit)
        end if
      end if
    end if

    if (library .and. allocated(eigval)) eig_found = .true.

    return

105 call io_error('Error: Problem opening eigenvalue file '//trim(seedname)//'.eig', stdout, seedname)
106 call io_error('Error: Problem reading eigenvalue file '//trim(seedname)//'.eig', stdout, seedname)

  end subroutine param_read_eigvals

  subroutine param_read_dis_manifold(eig_found, dis_window, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: eig_found
    !real(kind=dp), intent(in) :: eigval(:, :)
    type(dis_manifold_type), intent(inout) :: dis_window
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    !integer, intent(in) :: num_bands, num_wann
    !integer :: nkp, ierr
    logical :: found, found2

    !dis_window%win_min = -1.0_dp; dis_window%win_max = 0.0_dp
    !if (eig_found) dis_window%win_min = minval(eigval)
    call param_get_keyword(stdout, seedname, 'dis_win_min', found, r_value=dis_window%win_min)

    !if (eig_found) dis_window%win_max = maxval(eigval)
    call param_get_keyword(stdout, seedname, 'dis_win_max', found, r_value=dis_window%win_max)
    if (eig_found .and. (dis_window%win_max .lt. dis_window%win_min)) &
      call io_error('Error: param_read: check disentanglement windows', stdout, seedname)

    dis_window%froz_min = -1.0_dp; dis_window%froz_max = 0.0_dp
    ! no default for dis_froz_max
    dis_window%frozen_states = .false.
    call param_get_keyword(stdout, seedname, 'dis_froz_max', found, r_value=dis_window%froz_max)
    if (found) then
      dis_window%frozen_states = .true.
      dis_window%froz_min = dis_window%win_min ! default value for the bottom of frozen window
    end if
    call param_get_keyword(stdout, seedname, 'dis_froz_min', found2, r_value=dis_window%froz_min)
    if (eig_found) then
      if (dis_window%froz_max .lt. dis_window%froz_min) &
        call io_error('Error: param_read: check disentanglement frozen windows', stdout, seedname)
      if (found2 .and. .not. found) &
        call io_error('Error: param_read: found dis_froz_min but not dis_froz_max', stdout, seedname)
    endif
    ! ndimwin/lwindow are not read
  end subroutine param_read_dis_manifold

  subroutine param_read_kmesh_data(kmesh_data, stdout, seedname)
    use w90_io, only: io_error
!   use w90_utility, only: utility_recip_lattice
    implicit none
    type(kmesh_input_type), intent(out) :: kmesh_data
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    !real(kind=dp) :: real_lattice_tmp(3, 3), cell_volume
    integer :: itmp, ierr
    logical :: found

    kmesh_data%search_shells = 36
    call param_get_keyword(stdout, seedname, 'search_shells', found, i_value=kmesh_data%search_shells)
    if (kmesh_data%search_shells < 0) call io_error('Error: search_shells must be positive', stdout, seedname)

    kmesh_data%tol = 0.000001_dp
    call param_get_keyword(stdout, seedname, 'kmesh_tol', found, r_value=kmesh_data%tol)
    if (kmesh_data%tol < 0.0_dp) call io_error('Error: kmesh_tol must be positive', stdout, seedname)

    kmesh_data%num_shells = 0
    call param_get_range_vector(stdout, seedname, 'shell_list', found, kmesh_data%num_shells, lcount=.true.)
    if (found) then
      if (kmesh_data%num_shells < 0 .or. kmesh_data%num_shells > max_shells) &
        call io_error('Error: number of shell in shell_list must be between zero and six', stdout, seedname)
      if (allocated(kmesh_data%shell_list)) deallocate (kmesh_data%shell_list)
      allocate (kmesh_data%shell_list(kmesh_data%num_shells), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating shell_list in param_read', stdout, seedname)
      call param_get_range_vector(stdout, seedname, 'shell_list', found, kmesh_data%num_shells, .false., kmesh_data%shell_list)
      if (any(kmesh_data%shell_list < 1)) &
        call io_error('Error: shell_list must contain positive numbers', stdout, seedname)
    else
      if (allocated(kmesh_data%shell_list)) deallocate (kmesh_data%shell_list)
      allocate (kmesh_data%shell_list(max_shells), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating shell_list in param_read', stdout, seedname)
    end if

    call param_get_keyword(stdout, seedname, 'num_shells', found, i_value=itmp)
    if (found .and. (itmp /= kmesh_data%num_shells)) &
      call io_error('Error: Found obsolete keyword num_shells. Its value does not agree with shell_list', stdout, seedname)

    ! If .true., does not perform the check of B1 of
    ! Marzari, Vanderbild, PRB 56, 12847 (1997)
    ! in kmesh.F90
    ! mainly needed for the interaction with Z2PACK
    ! By default: .false. (perform the tests)
    kmesh_data%skip_B1_tests = .false.
    call param_get_keyword(stdout, seedname, 'skip_b1_tests', found, l_value=kmesh_data%skip_B1_tests)

  end subroutine param_read_kmesh_data

  subroutine param_read_kpoints(pw90_effective_model, library, k_points, num_kpts, &
                                recip_lattice, bohr, stdout, seedname)
    use w90_io, only: io_error
!   use w90_utility, only: utility_recip_lattice
    implicit none
    logical, intent(in) :: pw90_effective_model, library
    type(k_points_type), intent(inout) :: k_points
    integer, intent(in) :: num_kpts
    integer, intent(in) :: stdout
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname
    !real(kind=dp) :: real_lattice_tmp(3, 3), cell_volume
    integer :: nkp, ierr
    logical :: found

    if (.not. pw90_effective_model) allocate (k_points%kpt_cart(3, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating kpt_cart in param_read', stdout, seedname)
    if (.not. library) then
      allocate (k_points%kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpt_latt in param_read', stdout, seedname)
    end if

    call param_get_keyword_block(stdout, seedname, 'kpoints', found, num_kpts, 3, bohr, r_value=k_points%kpt_cart)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <kpoints> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      k_points%kpt_latt = k_points%kpt_cart
      if (.not. found) call io_error('Error: Did not find the kpoint information in the input file', stdout, seedname)
    end if

    ! Calculate the kpoints in cartesian coordinates
    if (.not. pw90_effective_model) then
      do nkp = 1, num_kpts
        k_points%kpt_cart(:, nkp) = matmul(k_points%kpt_latt(:, nkp), recip_lattice(:, :))
      end do
    endif

  end subroutine param_read_kpoints

  subroutine param_read_lattice(library, real_lattice, recip_lattice, bohr, stdout, seedname)
    use w90_io, only: io_error
    use w90_utility, only: utility_recip_lattice
    implicit none
    logical, intent(in) :: library
    integer, intent(in) :: stdout
    real(kind=dp), intent(out) :: real_lattice(3, 3), recip_lattice(3, 3)
    real(kind=dp) :: real_lattice_tmp(3, 3), cell_volume
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname

    logical :: found

    call param_get_keyword_block(stdout, seedname, 'unit_cell_cart', found, 3, 3, bohr, r_value=real_lattice_tmp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <unit_cell_cart> in input file'
    if (.not. library) then
      real_lattice = transpose(real_lattice_tmp)
      if (.not. found) call io_error('Error: Did not find the cell information in the input file', stdout, seedname)
    end if

    if (.not. library) &
      call utility_recip_lattice(real_lattice, recip_lattice, cell_volume, stdout, seedname)
    !call utility_metric(real_lattice, recip_lattice, real_metric, recip_metric)
  end subroutine param_read_lattice

  subroutine param_read_global_kmesh(global_kmesh_set, kmesh_spacing, kmesh, recip_lattice, stdout, seedname)
    use w90_io, only: io_error
    !use w90_utility, only: utility_recip_lattice
    implicit none
    integer, intent(in) :: stdout
    logical, intent(out) :: global_kmesh_set
    real(kind=dp), intent(out) :: kmesh_spacing
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    integer, intent(out) :: kmesh(3)
    character(len=50), intent(in)  :: seedname

    integer :: i
    logical :: found

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! k meshes                                                                                 !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! [GP-begin, Apr13, 2012]
    ! Global interpolation k-mesh; this is overridden by "local" meshes of a given submodule
    ! This bit of code must appear *before* all other codes for the local interpolation meshes,
    ! BUT *after* having calculated the reciprocal-space vectors.
    global_kmesh_set = .false.
    kmesh_spacing = -1._dp
    kmesh = 0
    call param_get_keyword(stdout, seedname, 'kmesh_spacing', found, r_value=kmesh_spacing)
    if (found) then
      if (kmesh_spacing .le. 0._dp) &
        call io_error('Error: kmesh_spacing must be greater than zero', stdout, seedname)
      global_kmesh_set = .true.

      call internal_set_kmesh(kmesh_spacing, recip_lattice, kmesh)
    end if
    call param_get_vector_length(stdout, seedname, 'kmesh', found, length=i)
    if (found) then
      if (global_kmesh_set) &
        call io_error('Error: cannot set both kmesh and kmesh_spacing', stdout, seedname)
      if (i .eq. 1) then
        global_kmesh_set = .true.
        call param_get_keyword_vector(stdout, seedname, 'kmesh', found, 1, i_value=kmesh)
        kmesh(2) = kmesh(1)
        kmesh(3) = kmesh(1)
      elseif (i .eq. 3) then
        global_kmesh_set = .true.
        call param_get_keyword_vector(stdout, seedname, 'kmesh', found, 3, i_value=kmesh)
      else
        call io_error('Error: kmesh must be provided as either one integer or a vector of three integers', stdout, seedname)
      end if
      if (any(kmesh <= 0)) &
        call io_error('Error: kmesh elements must be greater than zero', stdout, seedname)
    end if
    ! [GP-end]
  end subroutine param_read_global_kmesh

  subroutine param_read_atoms(library, atoms, real_lattice, recip_lattice, bohr, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: library
    integer, intent(in) :: stdout
    type(atom_data_type), intent(inout) :: atoms
    real(kind=dp), intent(in) :: real_lattice(3, 3), recip_lattice(3, 3)
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname

    integer :: i_temp, i_temp2
    logical :: found, found2, lunits

    ! Atoms
    if (.not. library) atoms%num_atoms = 0
    call param_get_block_length(stdout, seedname, 'atoms_frac', found, i_temp, library)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <atoms_frac> in input file'
    call param_get_block_length(stdout, seedname, 'atoms_cart', found2, i_temp2, library, lunits)
    if (found2 .and. library) write (stdout, '(a)') ' Ignoring <atoms_cart> in input file'
    if (.not. library) then
      if (found .and. found2) call io_error('Error: Cannot specify both atoms_frac and atoms_cart', stdout, seedname)
      if (found .and. i_temp > 0) then
        lunits = .false.
        atoms%num_atoms = i_temp
      elseif (found2 .and. i_temp2 > 0) then
        atoms%num_atoms = i_temp2
        if (lunits) atoms%num_atoms = atoms%num_atoms - 1
      end if
      if (atoms%num_atoms > 0) then
        call param_get_atoms(atoms, library, lunits, real_lattice, recip_lattice, bohr, stdout, seedname)
      end if
    endif
  end subroutine param_read_atoms

  subroutine get_all_keywords(stdout, seedname)
    ! wannier90.x and postw90.x now only read their own subset of the valid tokens in the ctrl file
    ! checking of the ctrl file is by testing for the presence of any remaining strings in the file
    ! after removing all valid keys.
    !
    ! this routine hoovers up any remaining keys by scanning the ctrl file for (the union of) all
    ! wannier90.x and postw90.x keywords.  The param_get_keyword* functions only assign to optional
    ! arguments: here we call without any, which has the side effect of clearing the input stream.
    !
    ! these lists have been populated using a grep command on the source; it needs to be updated by
    ! hand when the code changes.  There are a lot of keywords; it's not an ideal solution.
    !
    ! (for _vector: just specify zero length)
    ! (for _block: small modification to skip checking/failure when rows=0 )
    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in) :: seedname

    !type(kpoint_path_type) :: spec_points ! for the special case of param_get_keyword_kpath
    logical :: found

    ! keywords for wannier.x
    call param_get_keyword_block(stdout, seedname, 'dis_spheres', found, 0, 0, 0.0_dp)
    call param_get_keyword_block(stdout, seedname, 'kpoints', found, 0, 0, 0.0_dp)
    call param_get_keyword_block(stdout, seedname, 'nnkpts', found, 0, 0, 0.0_dp)
    call param_get_keyword_block(stdout, seedname, 'unit_cell_cart', found, 0, 0, 0.0_dp)
    call param_clear_block(stdout, seedname, 'projections')
    call param_clear_block(stdout, seedname, 'kpoint_path')
    call param_get_keyword(stdout, seedname, 'auto_projections', found)
    call param_get_keyword(stdout, seedname, 'bands_num_points', found)
    call param_get_keyword(stdout, seedname, 'bands_plot_dim', found)
    call param_get_keyword(stdout, seedname, 'bands_plot_format', found)
    call param_get_keyword(stdout, seedname, 'bands_plot', found)
    call param_get_keyword(stdout, seedname, 'bands_plot_mode', found)
    call param_get_keyword(stdout, seedname, 'calc_only_A', found)
    call param_get_keyword(stdout, seedname, 'conv_noise_amp', found)
    call param_get_keyword(stdout, seedname, 'conv_noise_num', found)
    call param_get_keyword(stdout, seedname, 'conv_tol', found)
    call param_get_keyword(stdout, seedname, 'conv_window', found)
    call param_get_keyword(stdout, seedname, 'cp_pp', found)
    call param_get_keyword(stdout, seedname, 'devel_flag', found)
    call param_get_keyword(stdout, seedname, 'dis_conv_tol', found)
    call param_get_keyword(stdout, seedname, 'dis_conv_window', found)
    call param_get_keyword(stdout, seedname, 'dis_froz_max', found)
    call param_get_keyword(stdout, seedname, 'dis_froz_min', found)
    call param_get_keyword(stdout, seedname, 'dis_mix_ratio', found)
    call param_get_keyword(stdout, seedname, 'dis_num_iter', found)
    call param_get_keyword(stdout, seedname, 'dis_spheres_first_wann', found)
    call param_get_keyword(stdout, seedname, 'dis_spheres_num', found)
    call param_get_keyword(stdout, seedname, 'dist_cutoff', found)
    call param_get_keyword(stdout, seedname, 'dist_cutoff_hc', found)
    call param_get_keyword(stdout, seedname, 'dist_cutoff_mode', found)
    call param_get_keyword(stdout, seedname, 'dis_win_max', found)
    call param_get_keyword(stdout, seedname, 'dis_win_min', found)
    call param_get_keyword(stdout, seedname, 'energy_unit', found)
    call param_get_keyword(stdout, seedname, 'fermi_energy', found)
    call param_get_keyword(stdout, seedname, 'fermi_energy_max', found)
    call param_get_keyword(stdout, seedname, 'fermi_energy_min', found)
    call param_get_keyword(stdout, seedname, 'fermi_energy_step', found)
    call param_get_keyword(stdout, seedname, 'fermi_surface_num_points', found)
    call param_get_keyword(stdout, seedname, 'fermi_surface_plot_format', found)
    call param_get_keyword(stdout, seedname, 'fermi_surface_plot', found)
    call param_get_keyword(stdout, seedname, 'fixed_step', found)
    call param_get_keyword(stdout, seedname, 'gamma_only', found)
    call param_get_keyword(stdout, seedname, 'guiding_centres', found)
    call param_get_keyword(stdout, seedname, 'hr_cutoff', found)
    call param_get_keyword(stdout, seedname, 'hr_plot', found)
    call param_get_keyword(stdout, seedname, 'iprint', found)
    call param_get_keyword(stdout, seedname, 'kmesh_spacing', found)
    call param_get_keyword(stdout, seedname, 'kmesh_tol', found)
    call param_get_keyword(stdout, seedname, 'length_unit', found)
    call param_get_keyword(stdout, seedname, 'num_bands', found)
    call param_get_keyword(stdout, seedname, 'num_cg_steps', found)
    call param_get_keyword(stdout, seedname, 'num_dump_cycles', found)
    call param_get_keyword(stdout, seedname, 'num_elec_per_state', found)
    call param_get_keyword(stdout, seedname, 'num_guide_cycles', found)
    call param_get_keyword(stdout, seedname, 'num_iter', found)
    call param_get_keyword(stdout, seedname, 'num_no_guide_iter', found)
    call param_get_keyword(stdout, seedname, 'num_print_cycles', found)
    call param_get_keyword(stdout, seedname, 'num_shells', found)
    call param_get_keyword(stdout, seedname, 'num_valence_bands', found)
    call param_get_keyword(stdout, seedname, 'num_wann', found)
    call param_get_keyword(stdout, seedname, 'one_dim_axis', found)
    call param_get_keyword(stdout, seedname, 'optimisation', found)
    call param_get_keyword(stdout, seedname, 'postproc_setup', found)
    call param_get_keyword(stdout, seedname, 'precond', found)
    call param_get_keyword(stdout, seedname, 'restart', found)
    call param_get_keyword(stdout, seedname, 'search_shells', found)
    call param_get_keyword(stdout, seedname, 'site_symmetry', found)
    call param_get_keyword(stdout, seedname, 'skip_b1_tests', found)
    call param_get_keyword(stdout, seedname, 'slwf_constrain', found)
    call param_get_keyword(stdout, seedname, 'slwf_lambda', found)
    call param_get_keyword(stdout, seedname, 'slwf_num', found)
    call param_get_keyword(stdout, seedname, 'spin', found)
    call param_get_keyword(stdout, seedname, 'spinors', found)
    call param_get_keyword(stdout, seedname, 'symmetrize_eps', found)
    call param_get_keyword(stdout, seedname, 'timing_level', found)
    call param_get_keyword(stdout, seedname, 'tran_easy_fix', found)
    call param_get_keyword(stdout, seedname, 'tran_energy_step', found)
    call param_get_keyword(stdout, seedname, 'tran_group_threshold', found)
    call param_get_keyword(stdout, seedname, 'tran_num_bandc', found)
    call param_get_keyword(stdout, seedname, 'tran_num_bb', found)
    call param_get_keyword(stdout, seedname, 'tran_num_cc', found)
    call param_get_keyword(stdout, seedname, 'tran_num_cell_ll', found)
    call param_get_keyword(stdout, seedname, 'tran_num_cell_rr', found)
    call param_get_keyword(stdout, seedname, 'tran_num_cr', found)
    call param_get_keyword(stdout, seedname, 'tran_num_lc', found)
    call param_get_keyword(stdout, seedname, 'tran_num_ll', found)
    call param_get_keyword(stdout, seedname, 'tran_num_rr', found)
    call param_get_keyword(stdout, seedname, 'tran_read_ht', found)
    call param_get_keyword(stdout, seedname, 'translate_home_cell', found)
    call param_get_keyword(stdout, seedname, 'transport', found)
    call param_get_keyword(stdout, seedname, 'transport_mode', found)
    call param_get_keyword(stdout, seedname, 'tran_use_same_lead', found)
    call param_get_keyword(stdout, seedname, 'tran_win_max', found)
    call param_get_keyword(stdout, seedname, 'tran_win_min', found)
    call param_get_keyword(stdout, seedname, 'tran_write_ht', found)
    call param_get_keyword(stdout, seedname, 'trial_step', found)
    call param_get_keyword(stdout, seedname, 'use_bloch_phases', found)
    call param_get_keyword(stdout, seedname, 'use_ws_distance', found)
    call param_get_keyword(stdout, seedname, 'wannier_plot_format', found)
    call param_get_keyword(stdout, seedname, 'wannier_plot', found)
    call param_get_keyword(stdout, seedname, 'wannier_plot_mode', found)
    call param_get_keyword(stdout, seedname, 'wannier_plot_radius', found)
    call param_get_keyword(stdout, seedname, 'wannier_plot_scale', found)
    call param_get_keyword(stdout, seedname, 'wannier_plot_spinor_mode', found)
    call param_get_keyword(stdout, seedname, 'wannier_plot_spinor_phase', found)
    call param_get_keyword(stdout, seedname, 'write_bvec', found)
    call param_get_keyword(stdout, seedname, 'write_hr_diag', found)
    call param_get_keyword(stdout, seedname, 'write_hr', found)
    call param_get_keyword(stdout, seedname, 'write_proj', found)
    call param_get_keyword(stdout, seedname, 'write_r2mn', found)
    call param_get_keyword(stdout, seedname, 'write_rmn', found)
    call param_get_keyword(stdout, seedname, 'write_tb', found)
    call param_get_keyword(stdout, seedname, 'write_u_matrices', found)
    call param_get_keyword(stdout, seedname, 'write_vdw_data', found)
    call param_get_keyword(stdout, seedname, 'write_xyz', found)
    call param_get_keyword(stdout, seedname, 'ws_distance_tol', found)
    call param_get_keyword(stdout, seedname, 'wvfn_formatted', found)
    call param_get_keyword_vector(stdout, seedname, 'kmesh', found, 0) ! the absent arrays have zero length ;-)
    call param_get_keyword_vector(stdout, seedname, 'mp_grid', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'translation_centre_frac', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'wannier_plot_supercell', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'ws_search_size', found, 0)
    ! ends list of wannier.x keywords

    ! keywords for postw90.x
    call param_get_keyword(stdout, seedname, 'adpt_smr_fac', found)
    call param_get_keyword(stdout, seedname, 'adpt_smr', found)
    call param_get_keyword(stdout, seedname, 'adpt_smr_max', found)
    call param_get_keyword(stdout, seedname, 'berry_curv_adpt_kmesh', found)
    call param_get_keyword(stdout, seedname, 'berry_curv_adpt_kmesh_thresh', found)
    call param_get_keyword(stdout, seedname, 'berry_curv_unit', found)
    call param_get_keyword(stdout, seedname, 'berry', found)
    call param_get_keyword(stdout, seedname, 'berry_kmesh_spacing', found)
    call param_get_keyword(stdout, seedname, 'berry_task', found)
    call param_get_keyword(stdout, seedname, 'boltz_2d_dir', found)
    call param_get_keyword(stdout, seedname, 'boltz_bandshift_energyshift', found)
    call param_get_keyword(stdout, seedname, 'boltz_bandshift_firstband', found)
    call param_get_keyword(stdout, seedname, 'boltz_bandshift', found)
    call param_get_keyword(stdout, seedname, 'boltz_calc_also_dos', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr_fac', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr_max', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_energy_max', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_energy_min', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_energy_step', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_smr_fixed_en_width', found)
    call param_get_keyword(stdout, seedname, 'boltz_dos_smr_type', found)
    call param_get_keyword(stdout, seedname, 'boltz_kmesh_spacing', found)
    call param_get_keyword(stdout, seedname, 'boltz_mu_max', found)
    call param_get_keyword(stdout, seedname, 'boltz_mu_min', found)
    call param_get_keyword(stdout, seedname, 'boltz_mu_step', found)
    call param_get_keyword(stdout, seedname, 'boltz_relax_time', found)
    call param_get_keyword(stdout, seedname, 'boltz_tdf_energy_step', found)
    call param_get_keyword(stdout, seedname, 'boltz_tdf_smr_fixed_en_width', found)
    call param_get_keyword(stdout, seedname, 'boltz_tdf_smr_type', found)
    call param_get_keyword(stdout, seedname, 'boltz_temp_max', found)
    call param_get_keyword(stdout, seedname, 'boltz_temp_min', found)
    call param_get_keyword(stdout, seedname, 'boltz_temp_step', found)
    call param_get_keyword(stdout, seedname, 'boltzwann', found)
    call param_get_keyword(stdout, seedname, 'degen_thr', found)
    call param_get_keyword(stdout, seedname, 'dos_adpt_smr_fac', found)
    call param_get_keyword(stdout, seedname, 'dos_adpt_smr', found)
    call param_get_keyword(stdout, seedname, 'dos_adpt_smr_max', found)
    call param_get_keyword(stdout, seedname, 'dos_energy_max', found)
    call param_get_keyword(stdout, seedname, 'dos_energy_min', found)
    call param_get_keyword(stdout, seedname, 'dos_energy_step', found)
    call param_get_keyword(stdout, seedname, 'dos', found)
    call param_get_keyword(stdout, seedname, 'dos_kmesh_spacing', found)
    call param_get_keyword(stdout, seedname, 'dos_smr_fixed_en_width', found)
    call param_get_keyword(stdout, seedname, 'dos_smr_type', found)
    call param_get_keyword(stdout, seedname, 'dos_task', found)
    call param_get_keyword(stdout, seedname, 'effective_model', found)
    call param_get_keyword(stdout, seedname, 'geninterp_alsofirstder', found)
    call param_get_keyword(stdout, seedname, 'geninterp', found)
    call param_get_keyword(stdout, seedname, 'geninterp_single_file', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_degen_thresh', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_eigval_max', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_freq_max', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_freq_min', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_freq_step', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_kmesh_spacing', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_smr_fixed_en_width', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_smr_max_arg', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_smr_type', found)
    call param_get_keyword(stdout, seedname, 'gyrotropic_task', found)
    call param_get_keyword(stdout, seedname, 'kpath_bands_colour', found)
    call param_get_keyword(stdout, seedname, 'kpath', found)
    call param_get_keyword(stdout, seedname, 'kpath_num_points', found)
    call param_get_keyword(stdout, seedname, 'kpath_task', found)
    call param_get_keyword(stdout, seedname, 'kslice_fermi_lines_colour', found)
    call param_get_keyword(stdout, seedname, 'kslice', found)
    call param_get_keyword(stdout, seedname, 'kslice_task', found)
    call param_get_keyword(stdout, seedname, 'kubo_adpt_smr_fac', found)
    call param_get_keyword(stdout, seedname, 'kubo_adpt_smr', found)
    call param_get_keyword(stdout, seedname, 'kubo_adpt_smr_max', found)
    call param_get_keyword(stdout, seedname, 'kubo_eigval_max', found)
    call param_get_keyword(stdout, seedname, 'kubo_freq_max', found)
    call param_get_keyword(stdout, seedname, 'kubo_freq_min', found)
    call param_get_keyword(stdout, seedname, 'kubo_freq_step', found)
    call param_get_keyword(stdout, seedname, 'kubo_smr_fixed_en_width', found)
    call param_get_keyword(stdout, seedname, 'kubo_smr_type', found)
    call param_get_keyword(stdout, seedname, 'sc_eta', found)
    call param_get_keyword(stdout, seedname, 'scissors_shift', found)
    call param_get_keyword(stdout, seedname, 'sc_phase_conv', found)
    call param_get_keyword(stdout, seedname, 'sc_w_thr', found)
    call param_get_keyword(stdout, seedname, 'shc_alpha', found)
    call param_get_keyword(stdout, seedname, 'shc_bandshift_energyshift', found)
    call param_get_keyword(stdout, seedname, 'shc_bandshift_firstband', found)
    call param_get_keyword(stdout, seedname, 'shc_bandshift', found)
    call param_get_keyword(stdout, seedname, 'shc_beta', found)
    call param_get_keyword(stdout, seedname, 'shc_freq_scan', found)
    call param_get_keyword(stdout, seedname, 'shc_gamma', found)
    call param_get_keyword(stdout, seedname, 'smr_fixed_en_width', found)
    call param_get_keyword(stdout, seedname, 'smr_max_arg', found)
    call param_get_keyword(stdout, seedname, 'smr_type', found)
    call param_get_keyword(stdout, seedname, 'spin_axis_azimuth', found)
    call param_get_keyword(stdout, seedname, 'spin_axis_polar', found)
    call param_get_keyword(stdout, seedname, 'spin_decomp', found)
    call param_get_keyword(stdout, seedname, 'spin_kmesh_spacing', found)
    call param_get_keyword(stdout, seedname, 'spin_moment', found)
    call param_get_keyword(stdout, seedname, 'spn_formatted', found)
    call param_get_keyword(stdout, seedname, 'transl_inv', found)
    call param_get_keyword(stdout, seedname, 'uhu_formatted', found)
    call param_get_keyword(stdout, seedname, 'use_degen_pert', found)
    call param_get_keyword(stdout, seedname, 'wanint_kpoint_file', found)
    call param_get_keyword_vector(stdout, seedname, 'berry_kmesh', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'boltz_kmesh', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'dos_kmesh', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'gyrotropic_box_b1', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'gyrotropic_box_b2', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'gyrotropic_box_b3', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'gyrotropic_box_center', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'gyrotropic_kmesh', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'kslice_2dkmesh', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'kslice_b1', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'kslice_b2', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'kslice_corner', found, 0)
    call param_get_keyword_vector(stdout, seedname, 'spin_kmesh', found, 0)
    ! ends list of postw90 keywords

  end subroutine get_all_keywords

  subroutine param_clean_infile(stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: loop, ierr

    ! JJ, filter out any remaining accepted keywords
    call get_all_keywords(stdout, seedname)

    if (any(len_trim(in_data(:)) > 0)) then
      write (stdout, '(1x,a)') 'The following section of file '//trim(seedname)//'.win contained unrecognised keywords'
      write (stdout, *)
      do loop = 1, num_lines
        if (len_trim(in_data(loop)) > 0) then
          write (stdout, '(1x,a)') trim(in_data(loop))
        end if
      end do
      write (stdout, *)
      call io_error('Unrecognised keyword(s) in input file, see also output file', stdout, seedname)
    end if

    deallocate (in_data, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating in_data in param_read', stdout, seedname)

  end subroutine param_clean_infile

  subroutine param_read_final_alloc(disentanglement, dis_window, wann_data, &
                                    num_wann, num_bands, num_kpts, stdout, seedname)
    ! =============================== !
    ! Some checks and initialisations !
    ! =============================== !
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: disentanglement
    !type(parameter_input_type), intent(inout) :: param_input
    type(dis_manifold_type), intent(inout) :: dis_window
    !type(param_wannierise_type), intent(inout) :: param_wannierise
    type(wannier_data_type), intent(inout) :: wann_data
    integer, intent(in) :: num_wann, num_bands, num_kpts
    character(len=50), intent(in)  :: seedname

    integer :: ierr

!    if (restart.ne.' ') disentanglement=.false.

    if (disentanglement) then
      if (allocated(dis_window%ndimwin)) deallocate (dis_window%ndimwin)
      allocate (dis_window%ndimwin(num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating ndimwin in param_read', stdout, seedname)
      if (allocated(dis_window%lwindow)) deallocate (dis_window%lwindow)
      allocate (dis_window%lwindow(num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating lwindow in param_read', stdout, seedname)
    endif

!    if ( wannier_plot .and. (index(wannier_plot_format,'cub').ne.0) ) then
!       cosa(1)=dot_product(real_lattice(1,:),real_lattice(2,:))
!       cosa(2)=dot_product(real_lattice(1,:),real_lattice(3,:))
!       cosa(3)=dot_product(real_lattice(2,:),real_lattice(3,:))
!       cosa = abs(cosa)
!       if (any(cosa.gt.eps6)) &
!            call io_error('Error: plotting in cube format requires orthogonal lattice vectors')
!    endif

    if (allocated(wann_data%centres)) deallocate (wann_data%centres)
    allocate (wann_data%centres(3, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wannier_centres in param_read', stdout, seedname)
    wann_data%centres = 0.0_dp
    if (allocated(wann_data%spreads)) deallocate (wann_data%spreads)
    allocate (wann_data%spreads(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating wannier_spreads in param_read', stdout, seedname)
    wann_data%spreads = 0.0_dp
  end subroutine param_read_final_alloc

  subroutine internal_set_kmesh(spacing, reclat, mesh)
    !! This routines returns the three integers that define the interpolation k-mesh, satisfying
    !! the condition that the spacing between two neighboring points along each of the three
    !! k_x, k_y and k_z directions is at smaller than a given spacing.
    !!
    !! The reclat is defined as:
    !!   * 'b_1' = (recip_lattice(1,I), i=1,3)
    !!   * 'b_2' = (recip_lattice(2,I), i=1,3)
    !!   * 'b_3' = (recip_lattice(3,I), i=1,3)
    !!
    !!  spacing must be > 0 (and in particular different from zero). We don't check this here.
    !!
    implicit none
    real(kind=dp), intent(in) :: spacing
    !! Minimum spacing between neighboring points, in angstrom^(-1)
    real(kind=dp), dimension(3, 3), intent(in) :: reclat
    !! Matrix of the reciprocal lattice vectors in cartesian coordinates, in angstrom^(-1)
    integer, dimension(3), intent(out) :: mesh
    !! Will contain the three integers defining the interpolation k-mesh

    real(kind=dp), dimension(3) :: blen
    integer :: i

    do i = 1, 3
      blen(i) = sqrt(sum(reclat(i, :)**2))
    end do

    do i = 1, 3
      mesh(i) = int(floor(blen(i)/spacing)) + 1
    end do

  end subroutine internal_set_kmesh

  function param_get_smearing_type(smearing_index)
    !! This function returns a string describing the type of smearing
    !! associated to a given smr_index integer value.
    integer, intent(in) :: smearing_index
    !! The integer index for which we want to get the string
    character(len=80)   :: param_get_smearing_type

    character(len=4)   :: orderstr

    if (smearing_index > 0) then
      write (orderstr, '(I0)') smearing_index
      param_get_smearing_type = "Methfessel-Paxton of order "//trim(orderstr)
    else if (smearing_index .eq. 0) then
      param_get_smearing_type = "Gaussian"
    else if (smearing_index .eq. -1) then
      param_get_smearing_type = "Marzari-Vanderbilt cold smearing"
    else if (smearing_index .eq. -99) then
      param_get_smearing_type = "Fermi-Dirac smearing"
    else
      param_get_smearing_type = "Unknown type of smearing"
    end if

  end function param_get_smearing_type

  function param_get_convention_type(sc_phase_conv)
    !! This function returns a string describing the convention
    !! associated to a sc_phase_conv integer value.
    integer, intent(in) :: sc_phase_conv
    !! The integer index for which we want to get the string
    character(len=80)   :: param_get_convention_type

    !character(len=4)   :: orderstr

    if (sc_phase_conv .eq. 1) then
      param_get_convention_type = "Tight-binding convention"
    else if (sc_phase_conv .eq. 2) then
      param_get_convention_type = "Wannier90 convention"
    else
      param_get_convention_type = "Unknown type of convention"
    end if

  end function param_get_convention_type

  function get_smearing_index(string, keyword, stdout, seedname)
    !! This function parses a string containing the type of
    !! smearing and returns the correct index for the smearing_index variable
    !
    !! If the string is not valid, an io_error is issued
    use w90_io, only: io_error
    integer, intent(in) :: stdout
    character(len=*), intent(in) :: string
    !! The string read from input
    character(len=*), intent(in) :: keyword
    !! The keyword that was read (e.g., smr_type), so that we can print a more useful error message
    character(len=50), intent(in)  :: seedname
    integer :: get_smearing_index

    integer :: pos

    get_smearing_index = 0 ! To avoid warnings of unset variables

    if (index(string, 'm-v') > 0) then
      get_smearing_index = -1
    elseif (index(string, 'm-p') > 0) then
      pos = index(string, 'm-p')
      if (len(trim(string(pos + 3:))) .eq. 0) then
        ! If the string is only 'm-p', we assume that 'm-p1' was intended
        get_smearing_index = 1
      else
        read (string(pos + 3:), *, err=337) get_smearing_index
        if (get_smearing_index < 0) &
          call io_error('Wrong m-p smearing order in keyword '//trim(keyword), stdout, seedname)
      end if
    elseif (index(string, 'f-d') > 0) then
      get_smearing_index = -99
      ! Some aliases
    elseif (index(string, 'cold') > 0) then
      get_smearing_index = -1
    elseif (index(string, 'gauss') > 0) then
      get_smearing_index = 0
      ! Unrecognised keyword
    else
      call io_error('Unrecognised value for keyword '//trim(keyword), stdout, seedname)
    end if

    return

337 call io_error('Wrong m-p smearing order in keyword '//trim(keyword), stdout, seedname)

  end function get_smearing_index

!===================================================================
  subroutine param_uppercase(atoms, spec_points, length_unit)
    !===================================================================
    !                                                                  !
    !! Convert a few things to uppercase to look nice in the output
    !                                                                  !
    !===================================================================

    implicit none

    type(atom_data_type), intent(inout) :: atoms
    type(kpoint_path_type), intent(inout) :: spec_points
    character(len=*), intent(inout) :: length_unit
    integer :: nsp, ic, loop, inner_loop

    ! Atom labels (eg, si --> Si)
    do nsp = 1, atoms%num_species
      ic = ichar(atoms%label(nsp) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atoms%label(nsp) (1:1) = char(ic + ichar('Z') - ichar('z'))
    enddo

    do nsp = 1, atoms%num_species
      ic = ichar(atoms%symbol(nsp) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atoms%symbol(nsp) (1:1) = char(ic + ichar('Z') - ichar('z'))
    enddo

    ! Bands labels (eg, x --> X)
    if (allocated(spec_points%labels)) then
      do loop = 1, size(spec_points%labels)
        do inner_loop = 1, len(spec_points%labels(loop))
          ic = ichar(spec_points%labels(loop) (inner_loop:inner_loop))
          if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
            spec_points%labels(loop) (inner_loop:inner_loop) = char(ic + ichar('Z') - ichar('z'))
        enddo
      enddo
    endif

    ! Length unit (ang --> Ang, bohr --> Bohr)
    ic = ichar(length_unit(1:1))
    if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
      length_unit(1:1) = char(ic + ichar('Z') - ichar('z'))

    return

  end subroutine param_uppercase

  subroutine param_write_header(bohr_version_str, constants_version_str1, constants_version_str2, stdout)
    !! Write a suitable header for the calculation - version authors etc
    use w90_io, only: io_date, w90_version
    !use w90_constants, only: bohr_version_str, constants_version_str1, constants_version_str2
    implicit none

    integer, intent(in) :: stdout
    character(len=*), intent(in) :: bohr_version_str, constants_version_str1, constants_version_str2
    character(len=9) :: cdate, ctime

    call io_date(cdate, ctime)

    write (stdout, *)
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                   WANNIER90                       |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |        Welcome to the Maximally-Localized         |'
    write (stdout, *) '            |        Generalized Wannier Functions code         |'
    write (stdout, *) '            |            http://www.wannier.org                 |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  Wannier90 Developer Group:                       |'
    write (stdout, *) '            |    Giovanni Pizzi    (EPFL)                       |'
    write (stdout, *) '            |    Valerio Vitale    (Cambridge)                  |'
    write (stdout, *) '            |    David Vanderbilt  (Rutgers University)         |'
    write (stdout, *) '            |    Nicola Marzari    (EPFL)                       |'
    write (stdout, *) '            |    Ivo Souza         (Universidad del Pais Vasco) |'
    write (stdout, *) '            |    Arash A. Mostofi  (Imperial College London)    |'
    write (stdout, *) '            |    Jonathan R. Yates (University of Oxford)       |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  For the full list of Wannier90 3.x authors,      |'
    write (stdout, *) '            |  please check the code documentation and the      |'
    write (stdout, *) '            |  README on the GitHub page of the code            |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  Please cite                                      |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  [ref] "Wannier90 as a community code:            |'
    write (stdout, *) '            |        new features and applications",            |'
    write (stdout, *) '            |        G. Pizzi et al., J. Phys. Cond. Matt. 32,  |'
    write (stdout, *) '            |        165902 (2020).                             |'
    write (stdout, *) '            |        http://doi.org/10.1088/1361-648X/ab51ff    |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  in any publications arising from the use of      |'
    write (stdout, *) '            |  this code. For the method please cite            |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  [ref] "Maximally Localized Generalised Wannier   |'
    write (stdout, *) '            |         Functions for Composite Energy Bands"     |'
    write (stdout, *) '            |         N. Marzari and D. Vanderbilt              |'
    write (stdout, *) '            |         Phys. Rev. B 56 12847 (1997)              |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  [ref] "Maximally Localized Wannier Functions     |'
    write (stdout, *) '            |         for Entangled Energy Bands"               |'
    write (stdout, *) '            |         I. Souza, N. Marzari and D. Vanderbilt    |'
    write (stdout, *) '            |         Phys. Rev. B 65 035109 (2001)             |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | Copyright (c) 1996-2020                           |'
    write (stdout, *) '            |        The Wannier90 Developer Group and          |'
    write (stdout, *) '            |        individual contributors                    |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |      Release: ', adjustl(w90_version), '   5th March    2020      |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | This program is free software; you can            |'
    write (stdout, *) '            | redistribute it and/or modify it under the terms  |'
    write (stdout, *) '            | of the GNU General Public License as published by |'
    write (stdout, *) '            | the Free Software Foundation; either version 2 of |'
    write (stdout, *) '            | the License, or (at your option) any later version|'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | This program is distributed in the hope that it   |'
    write (stdout, *) '            | will be useful, but WITHOUT ANY WARRANTY; without |'
    write (stdout, *) '            | even the implied warranty of MERCHANTABILITY or   |'
    write (stdout, *) '            | FITNESS FOR A PARTICULAR PURPOSE. See the GNU     |'
    write (stdout, *) '            | General Public License for more details.          |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | You should have received a copy of the GNU General|'
    write (stdout, *) '            | Public License along with this program; if not,   |'
    write (stdout, *) '            | write to the Free Software Foundation, Inc.,      |'
    write (stdout, *) '            | 675 Mass Ave, Cambridge, MA 02139, USA.           |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) '            |    Execution started on ', cdate, ' at ', ctime, '    |'
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) ''
    write (stdout, '(1X,A)') '******************************************************************************'
    write (stdout, '(1X,A)') '* '//constants_version_str1//'*'
    write (stdout, '(1X,A)') '* '//constants_version_str2//'*'
    write (stdout, '(1X,A)') '* '//bohr_version_str//'*'
    write (stdout, '(1X,A)') '******************************************************************************'
    write (stdout, *) ''

  end subroutine param_write_header

!==================================================================!
  subroutine param_dealloc(exclude_bands, wann_data, input_proj, kmesh_data, k_points, &
                           dis_window, atoms, eigval, spec_points, stdout, seedname)
    !==================================================================!
    !                                                                  !
    !! release memory from allocated parameters
    !                                                                  !
    !===================================================================
    use w90_io, only: io_error

    implicit none
    !data from parameters module
    !type(param_driver_type), intent(inout) :: driver
    integer, allocatable, intent(inout) :: exclude_bands(:)
    type(wannier_data_type), intent(inout) :: wann_data
    type(proj_input_type), intent(inout) :: input_proj
    type(kmesh_input_type), intent(inout) :: kmesh_data
    type(k_points_type), intent(inout) :: k_points
    type(dis_manifold_type), intent(inout) :: dis_window
    type(atom_data_type), intent(inout) :: atoms
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    type(kpoint_path_type), intent(inout) :: spec_points
    character(len=50), intent(in)  :: seedname
    integer, intent(in) :: stdout

    integer :: ierr

    if (allocated(dis_window%ndimwin)) then
      deallocate (dis_window%ndimwin, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ndimwin in param_dealloc', stdout, seedname)
    end if
    if (allocated(dis_window%lwindow)) then
      deallocate (dis_window%lwindow, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating lwindow in param_dealloc', stdout, seedname)
    end if
    if (allocated(eigval)) then
      deallocate (eigval, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating eigval in param_dealloc', stdout, seedname)
    endif
    if (allocated(kmesh_data%shell_list)) then
      deallocate (kmesh_data%shell_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating shell_list in param_dealloc', stdout, seedname)
    endif
    if (allocated(k_points%kpt_latt)) then
      deallocate (k_points%kpt_latt, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kpt_latt in param_dealloc', stdout, seedname)
    endif
    if (allocated(k_points%kpt_cart)) then
      deallocate (k_points%kpt_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kpt_cart in param_dealloc', stdout, seedname)
    endif
    if (allocated(spec_points%labels)) then
      deallocate (spec_points%labels, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating labels in param_dealloc', stdout, seedname)
    end if
    if (allocated(spec_points%points)) then
      deallocate (spec_points%points, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating points in param_dealloc', stdout, seedname)
    end if
    if (allocated(atoms%label)) then
      deallocate (atoms%label, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_label in param_dealloc', stdout, seedname)
    end if
    if (allocated(atoms%symbol)) then
      deallocate (atoms%symbol, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_symbol in param_dealloc', stdout, seedname)
    end if
    if (allocated(atoms%pos_frac)) then
      deallocate (atoms%pos_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atom_pos_frac in param_dealloc', stdout, seedname)
    end if
    if (allocated(atoms%pos_cart)) then
      deallocate (atoms%pos_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_pos_cart in param_dealloc', stdout, seedname)
    end if
    if (allocated(atoms%species_num)) then
      deallocate (atoms%species_num, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_species_num in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%site)) then
      deallocate (input_proj%site, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_site in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%l)) then
      deallocate (input_proj%l, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_l in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%m)) then
      deallocate (input_proj%m, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_m in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%s)) then
      deallocate (input_proj%s, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_s in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%s_qaxis)) then
      deallocate (input_proj%s_qaxis, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_s_qaxis in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%z)) then
      deallocate (input_proj%z, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_z in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%x)) then
      deallocate (input_proj%x, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_x in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%radial)) then
      deallocate (input_proj%radial, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_radial in param_dealloc', stdout, seedname)
    end if
    if (allocated(input_proj%zona)) then
      deallocate (input_proj%zona, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_zona in param_dealloc', stdout, seedname)
    end if
    if (allocated(exclude_bands)) then
      deallocate (exclude_bands, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating exclude_bands in param_dealloc', stdout, seedname)
    end if
    if (allocated(wann_data%centres)) then
      deallocate (wann_data%centres, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_centres in param_dealloc', stdout, seedname)
    end if
    if (allocated(wann_data%spreads)) then
      deallocate (wann_data%spreads, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_spreads in param_dealloc', stdout, seedname)
    endif
    return

  end subroutine param_dealloc

!~  !================================!
!~  subroutine param_write_um
!~    !================================!
!~    !                                !
!~    ! Dump the U and M to *_um.dat   !
!~    !                                !
!~    !================================!
!~
!~
!~    use w90_io,        only : io_file_unit,io_error,seedname,io_date
!~    implicit none
!~
!~    integer :: i,j,k,l,um_unit
!~    character (len=9) :: cdate, ctime
!~    character(len=33) :: header
!~
!~    call io_date(cdate, ctime)
!~    header='written on '//cdate//' at '//ctime
!~
!~    um_unit=io_file_unit()
!~    open(unit=um_unit,file=trim(seedname)//'_um.dat',form='unformatted')
!~    write(um_unit) header
!~    write(um_unit) omega_invariant
!~    write(um_unit) num_wann,num_kpts,num_nnmax
!~    write(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
!~    write(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
!~    close(um_unit)
!~
!~    return
!~
!~  end subroutine param_write_um

!~  !================================!
!~  subroutine param_read_um
!~    !================================!
!~    !                                !
!~    ! Restore U and M from file      !
!~    !                                !
!~    !================================!
!~
!~    use w90_io,        only : io_file_unit,io_error,seedname
!~    implicit none
!~
!~    integer       :: tmp_num_wann,tmp_num_kpts,tmp_num_nnmax
!~    integer       :: i,j,k,l,um_unit,ierr
!~    character(len=33) :: header
!~    real(kind=dp) :: tmp_omi
!~
!~    um_unit=io_file_unit()
!~    open(unit=um_unit,file=trim(seedname)//'_um.dat',status="old",form='unformatted',err=105)
!~    read(um_unit) header
!~    write(stdout,'(1x,4(a))') 'Reading U and M from file ',trim(seedname),'_um.dat ', header
!~    read(um_unit) tmp_omi
!~    if ( have_disentangled ) then
!~       if ( abs(tmp_omi-omega_invariant).gt.1.0e-10_dp )  &
!~            call io_error('Error in restart: omega_invariant in .chk and um.dat files do not match')
!~    endif
!~    read(um_unit) tmp_num_wann,tmp_num_kpts,tmp_num_nnmax
!~    if(tmp_num_wann/=num_wann) call io_error('Error in param_read_um: num_wann mismatch')
!~    if(tmp_num_kpts/=num_kpts) call io_error('Error in param_read_um: num_kpts mismatch')
!~    if(tmp_num_nnmax/=num_nnmax) call io_error('Error in param_read_um: num_nnmax mismatch')
!~    if (.not.allocated(u_matrix)) then
!~       allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
!~       if (ierr/=0) call io_error('Error allocating u_matrix in param_read_um')
!~    endif
!~    read(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
!~    if (.not.allocated(m_matrix)) then
!~       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!~       if (ierr/=0) call io_error('Error allocating m_matrix in param_read_um')
!~    endif
!~    read(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
!~    close(um_unit)
!~
!~    return
!~
!~105 call io_error('Error: Problem opening file '//trim(seedname)//'_um.dat in param_read_um')
!~
! $  end subroutine param_read_um

!=================================================!
  subroutine param_read_chkpt(dis_data, exclude_bands, kmesh_info, k_points, wann_data, m_matrix, &
                              u_matrix, u_matrix_opt, real_lattice, recip_lattice, &
                              omega_invariant, mp_grid, num_bands, num_exclude_bands, num_kpts, &
                              num_wann, checkpoint, have_disentangled, ispostw90, seedname, stdout)
    !=================================================!
    !! Read checkpoint file
    !! IMPORTANT! If you change the chkpt format, adapt
    !! accordingly also the w90chk2chk.x utility!
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !! This function should be called
    !=================================================!

    use w90_constants, only: eps6
    use w90_io, only: io_file_unit, io_error

    implicit none

    !data from parameters module
    integer, allocatable, intent(inout) :: exclude_bands(:)
    type(wannier_data_type), intent(inout) :: wann_data
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(k_points_type), intent(in) :: k_points
    type(dis_manifold_type), intent(inout) :: dis_data

    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_exclude_bands

    complex(kind=dp), allocatable, intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix(:, :, :, :)

    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(inout) :: omega_invariant

    character(len=50), intent(in)  :: seedname
    character(len=*), intent(inout) :: checkpoint

    logical, intent(in) :: ispostw90 ! Are we running postw90?
    logical, intent(out) :: have_disentangled

!   local variables
    integer :: chk_unit, nkp, i, j, k, l, ntmp, ierr
    character(len=33) :: header
    real(kind=dp) :: tmp_latt(3, 3), tmp_kpt_latt(3, num_kpts)
    integer :: tmp_excl_bands(1:num_exclude_bands), tmp_mp_grid(1:3)

    write (stdout, '(1x,3a)') 'Reading restart information from file ', trim(seedname), '.chk :'

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk', status='old', form='unformatted', err=121)

    ! Read comment line
    read (chk_unit) header
    write (stdout, '(1x,a)', advance='no') trim(header)

    ! Consistency checks
    read (chk_unit) ntmp                           ! Number of bands
    if (ntmp .ne. num_bands) call io_error('param_read_chk: Mismatch in num_bands', stdout, seedname)
    read (chk_unit) ntmp                           ! Number of excluded bands
    if (ntmp .ne. num_exclude_bands) &
      call io_error('param_read_chk: Mismatch in num_exclude_bands', stdout, seedname)
    read (chk_unit) (tmp_excl_bands(i), i=1, num_exclude_bands) ! Excluded bands
    do i = 1, num_exclude_bands
      if (tmp_excl_bands(i) .ne. exclude_bands(i)) &
        call io_error('param_read_chk: Mismatch in exclude_bands', stdout, seedname)
    enddo
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Real lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - real_lattice(i, j)) .gt. eps6) &
          call io_error('param_read_chk: Mismatch in real_lattice', stdout, seedname)
      enddo
    enddo
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Reciprocal lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - recip_lattice(i, j)) .gt. eps6) &
          call io_error('param_read_chk: Mismatch in recip_lattice', stdout, seedname)
      enddo
    enddo
    read (chk_unit) ntmp                ! K-points
    if (ntmp .ne. num_kpts) &
      call io_error('param_read_chk: Mismatch in num_kpts', stdout, seedname)
    read (chk_unit) (tmp_mp_grid(i), i=1, 3)         ! M-P grid
    do i = 1, 3
      if (tmp_mp_grid(i) .ne. mp_grid(i)) &
        call io_error('param_read_chk: Mismatch in mp_grid', stdout, seedname)
    enddo
    read (chk_unit) ((tmp_kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts)
    do nkp = 1, num_kpts
      do i = 1, 3
        if (abs(tmp_kpt_latt(i, nkp) - k_points%kpt_latt(i, nkp)) .gt. eps6) &
          call io_error('param_read_chk: Mismatch in kpt_latt', stdout, seedname)
      enddo
    enddo
    read (chk_unit) ntmp                ! nntot
    if (ntmp .ne. kmesh_info%nntot) &
      call io_error('param_read_chk: Mismatch in nntot', stdout, seedname)
    read (chk_unit) ntmp                ! num_wann
    if (ntmp .ne. num_wann) &
      call io_error('param_read_chk: Mismatch in num_wann', stdout, seedname)
    ! End of consistency checks

    read (chk_unit) checkpoint             ! checkpoint
    checkpoint = adjustl(trim(checkpoint))

    read (chk_unit) have_disentangled      ! whether a disentanglement has been performed

    if (have_disentangled) then

      read (chk_unit) omega_invariant     ! omega invariant

      ! lwindow
      if (.not. allocated(dis_data%lwindow)) then
        allocate (dis_data%lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating lwindow in param_read_chkpt', stdout, seedname)
      endif
      read (chk_unit, err=122) ((dis_data%lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)

      ! ndimwin
      if (.not. allocated(dis_data%ndimwin)) then
        allocate (dis_data%ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ndimwin in param_read_chkpt', stdout, seedname)
      endif
      read (chk_unit, err=123) (dis_data%ndimwin(nkp), nkp=1, num_kpts)

      ! U_matrix_opt
      if (.not. allocated(u_matrix_opt)) then
        allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating u_matrix_opt in param_read_chkpt', stdout, seedname)
      endif
      read (chk_unit, err=124) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)

    endif

    ! U_matrix
    if (.not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating u_matrix in param_read_chkpt', stdout, seedname)
    endif
    read (chk_unit, err=125) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)

    ! M_matrix
    if (.not. allocated(m_matrix)) then
      allocate (m_matrix(num_wann, num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating m_matrix in param_read_chkpt', stdout, seedname)
    endif
    read (chk_unit, err=126) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, kmesh_info%nntot), l=1, num_kpts)

    ! wannier_centres
    read (chk_unit, err=127) ((wann_data%centres(i, j), i=1, 3), j=1, num_wann)

    ! wannier spreads
    read (chk_unit, err=128) (wann_data%spreads(i), i=1, num_wann)

    close (chk_unit)

    write (stdout, '(a/)') ' ... done'

    return

121 if (ispostw90) then
      call io_error('Error opening '//trim(seedname)//'.chk in param_read_chkpt: did you run wannier90.x first?', stdout, seedname)
    else
      call io_error('Error opening '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)
    end if
122 call io_error('Error reading lwindow from '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)
123 call io_error('Error reading ndimwin from '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)
124 call io_error('Error reading u_matrix_opt from '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)
125 call io_error('Error reading u_matrix from '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)
126 call io_error('Error reading m_matrix from '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)
127 call io_error('Error reading wannier_centres from '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)
128 call io_error('Error reading wannier_spreads from '//trim(seedname)//'.chk in param_read_chkpt', stdout, seedname)

  end subroutine param_read_chkpt

!===========================================================!
  subroutine param_chkpt_dist(dis_data, wann_data, u_matrix, u_matrix_opt, omega_invariant, &
                              num_bands, num_kpts, num_wann, checkpoint, have_disentangled, &
                              seedname, stdout, comm)
    !===========================================================!
    !                                                           !
    !! Distribute the chk files
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp !, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, io_date, io_time, io_stopwatch
    use w90_comms, only: comms_bcast, w90commtype, mpirank

    implicit none

    !data from parameters module
    type(wannier_data_type), intent(inout) :: wann_data
    type(dis_manifold_type), intent(inout) :: dis_data
    type(w90commtype), intent(in) :: comm

    integer, intent(in) :: stdout
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_wann
    integer, intent(inout) :: num_kpts

    complex(kind=dp), allocatable, intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :)
    real(kind=dp), intent(inout) :: omega_invariant

    character(len=50), intent(in)  :: seedname
    character(len=*), intent(inout) :: checkpoint
    logical, intent(inout) :: have_disentangled

!   local variables
    integer :: ierr

    logical :: on_root = .false.

    if (mpirank(comm) == 0) on_root = .true.

    call comms_bcast(checkpoint, len(checkpoint), stdout, seedname, comm)

    if (.not. on_root .and. .not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating u_matrix in param_chkpt_dist', stdout, seedname)
    endif
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, stdout, seedname, comm)

!    if (.not.on_root .and. .not.allocated(m_matrix)) then
!       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!       if (ierr/=0)&
!            call io_error('Error allocating m_matrix in param_chkpt_dist')
!    endif
!    call comms_bcast(m_matrix(1,1,1,1),num_wann*num_wann*nntot*num_kpts)

    call comms_bcast(have_disentangled, 1, stdout, seedname, comm)

    if (have_disentangled) then
      if (.not. on_root) then

        if (.not. allocated(u_matrix_opt)) then
          allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating u_matrix_opt in param_chkpt_dist', stdout, seedname)
        endif

        if (.not. allocated(dis_data%lwindow)) then
          allocate (dis_data%lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating lwindow in param_chkpt_dist', stdout, seedname)
        endif

        if (.not. allocated(dis_data%ndimwin)) then
          allocate (dis_data%ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating ndimwin in param_chkpt_dist', stdout, seedname)
        endif

      end if

      call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts, stdout, seedname, comm)
      call comms_bcast(dis_data%lwindow(1, 1), num_bands*num_kpts, stdout, seedname, comm)
      call comms_bcast(dis_data%ndimwin(1), num_kpts, stdout, seedname, comm)
      call comms_bcast(omega_invariant, 1, stdout, seedname, comm)
    end if
    call comms_bcast(wann_data%centres(1, 1), 3*num_wann, stdout, seedname, comm)
    call comms_bcast(wann_data%spreads(1), num_wann, stdout, seedname, comm)

  end subroutine param_chkpt_dist

!=======================================!
  subroutine param_in_file(seedname, stdout)
    !=======================================!
    !! Load the *.win file into a character
    !! array in_file, ignoring comments and
    !! blank lines and converting everything
    !! to lowercase characters
    !=======================================!

    use w90_utility, only: utility_lowercase
    use w90_io, only: io_file_unit, io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer           :: in_unit, tot_num_lines, ierr, line_counter, loop, in1, in2
    character(len=maxlen) :: dummy
    integer           :: pos
    character, parameter :: TABCHAR = char(9)

    in_unit = io_file_unit()
    open (in_unit, file=trim(seedname)//'.win', form='formatted', status='old', err=101)

    num_lines = 0; tot_num_lines = 0
    do
      read (in_unit, '(a)', iostat=ierr, err=200, end=210) dummy
      ! [GP-begin, Apr13, 2012]: I convert all tabulation characters to spaces
      pos = index(dummy, TABCHAR)
      do while (pos .ne. 0)
        dummy(pos:pos) = ' '
        pos = index(dummy, TABCHAR)
      end do
      ! [GP-end]
      dummy = adjustl(dummy)
      tot_num_lines = tot_num_lines + 1
      if (.not. dummy(1:1) == '!' .and. .not. dummy(1:1) == '#') then
        if (len(trim(dummy)) > 0) num_lines = num_lines + 1
      endif

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.win', stdout, seedname)
200 call io_error('Error: Problem reading input file '//trim(seedname)//'.win', stdout, seedname)
210 continue
    rewind (in_unit)

    allocate (in_data(num_lines), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating in_data in param_in_file', stdout, seedname)

    line_counter = 0
    do loop = 1, tot_num_lines
      read (in_unit, '(a)', iostat=ierr, err=200) dummy
      ! [GP-begin, Apr13, 2012]: I convert all tabulation characters to spaces
      pos = index(dummy, TABCHAR)
      do while (pos .ne. 0)
        dummy(pos:pos) = ' '
        pos = index(dummy, TABCHAR)
      end do
      ! [GP-end]
      dummy = utility_lowercase(dummy)
      dummy = adjustl(dummy)
      if (dummy(1:1) == '!' .or. dummy(1:1) == '#') cycle
      if (len(trim(dummy)) == 0) cycle
      line_counter = line_counter + 1
      in1 = index(dummy, '!')
      in2 = index(dummy, '#')
      if (in1 == 0 .and. in2 == 0) in_data(line_counter) = dummy
      if (in1 == 0 .and. in2 > 0) in_data(line_counter) = dummy(:in2 - 1)
      if (in2 == 0 .and. in1 > 0) in_data(line_counter) = dummy(:in1 - 1)
      if (in2 > 0 .and. in1 > 0) in_data(line_counter) = dummy(:min(in1, in2) - 1)
    end do

    close (in_unit)

  end subroutine param_in_file

!===========================================================================!
  subroutine param_get_keyword(stdout, seedname, keyword, found, c_value, l_value, i_value, r_value)
    !===========================================================================!
    !                                                                           !
    !! Finds the value of the required keyword.
    !                                                                           !
    !===========================================================================!

    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    character(*), optional, intent(inout) :: c_value
    !! Keyword value
    logical, optional, intent(inout) :: l_value
    !! Keyword value
    integer, optional, intent(inout) :: i_value
    !! Keyword value
    real(kind=dp), optional, intent(inout) :: r_value
    !! Keyword value

    integer           :: kl, in, loop, itmp
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      itmp = in + len(trim(keyword))
      if (in_data(loop) (itmp:itmp) /= '=' &
          .and. in_data(loop) (itmp:itmp) /= ':' &
          .and. in_data(loop) (itmp:itmp) /= ' ') cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      in_data(loop) (1:maxlen) = ' '
      dummy = adjustl(dummy)
      if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
        dummy = dummy(2:)
        dummy = adjustl(dummy)
      end if
    end do

    if (found) then
      if (present(c_value)) c_value = dummy
      if (present(l_value)) then
        if (index(dummy, 't') > 0) then
          l_value = .true.
        elseif (index(dummy, 'f') > 0) then
          l_value = .false.
        else
          call io_error('Error: Problem reading logical keyword '//trim(keyword), stdout, seedname)
        endif
      endif
      if (present(i_value)) read (dummy, *, err=220, end=220) i_value
      if (present(r_value)) read (dummy, *, err=220, end=220) r_value
    end if

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword), stdout, seedname)

  end subroutine param_get_keyword

!=========================================================================================!
  subroutine param_get_keyword_vector(stdout, seedname, keyword, found, length, c_value, l_value, &
                                      i_value, r_value)
    !=========================================================================================!
    !                                                                                         !
    !! Finds the values of the required keyword vector
    !                                                                                         !
    !=========================================================================================!

    use w90_io, only: io_error

    implicit none

    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(in)  :: length
    !! Length of vecotr to read
    character(*), optional, intent(inout) :: c_value(length)
    !! Keyword data
    logical, optional, intent(inout) :: l_value(length)
    !! Keyword data
    integer, optional, intent(inout) :: i_value(length)
    !! Keyword data
    real(kind=dp), optional, intent(inout) :: r_value(length)
    !! Keyword data
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer           :: kl, in, loop, i
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      in_data(loop) (1:maxlen) = ' '
      dummy = adjustl(dummy)
      if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
        dummy = dummy(2:)
        dummy = adjustl(dummy)
      end if
    end do

    if (found) then
      if (present(c_value)) read (dummy, *, err=230, end=230) (c_value(i), i=1, length)
      if (present(l_value)) then
        ! I don't think we need this. Maybe read into a dummy charater
        ! array and convert each element to logical
        call io_error('param_get_keyword_vector unimplemented for logicals', stdout, seedname)
      endif
      if (present(i_value)) read (dummy, *, err=230, end=230) (i_value(i), i=1, length)
      if (present(r_value)) read (dummy, *, err=230, end=230) (r_value(i), i=1, length)
    end if

    return

230 call io_error('Error: Problem reading keyword '//trim(keyword)//' in param_get_keyword_vector', stdout, seedname)

  end subroutine param_get_keyword_vector

!========================================================!
  subroutine param_get_vector_length(stdout, seedname, keyword, found, length)
    !======================================================!
    !                                                      !
    !! Returns the length of a keyword vector
    !                                                      !
    !======================================================!

    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(out)  :: length
    !! length of vector

    integer           :: kl, in, loop, pos
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      dummy = adjustl(dummy)
      if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
        dummy = dummy(2:)
        dummy = adjustl(dummy)
      end if
    end do

    length = 0
    if (found) then
      if (len_trim(dummy) == 0) call io_error('Error: keyword '//trim(keyword)//' is blank', stdout, seedname)
      length = 1
      dummy = adjustl(dummy)
      do
        pos = index(dummy, ' ')
        dummy = dummy(pos + 1:)
        dummy = adjustl(dummy)
        if (len_trim(dummy) > 0) then
          length = length + 1
        else
          exit
        endif

      end do

    end if

    return

  end subroutine param_get_vector_length

!==============================================================================================!
  subroutine param_get_keyword_block(stdout, seedname, keyword, found, rows, columns, bohr, c_value, &
                                     l_value, i_value, r_value)
    !==============================================================================================!
    !                                                                                              !
    !!   Finds the values of the required data block
    !                                                                                              !
    !==============================================================================================!

    !use w90_constants, only: bohr
    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(in)  :: rows
    !! Number of rows
    integer, intent(in)  :: columns
    !! Number of columns
    character(*), optional, intent(inout) :: c_value(columns, rows)
    !! keyword block data
    logical, optional, intent(inout) :: l_value(columns, rows)
    !! keyword block data
    integer, optional, intent(inout) :: i_value(columns, rows)
    !! keyword block data
    real(kind=dp), optional, intent(inout) :: r_value(columns, rows)
    !! keyword block data
    real(kind=dp), intent(in) :: bohr

    integer           :: in, ins, ine, loop, i, line_e, line_s, counter, blen
    logical           :: found_e, found_s, lconvert
    character(len=maxlen) :: dummy, end_st, start_st

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    if (.not. found_s) then
      found = .false.
      return
    end if

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    ! number of lines of data in block
    blen = line_e - line_s - 1

    !    if( blen /= rows) then
    !       if ( index(trim(keyword),'unit_cell_cart').ne.0 ) then
    !          if ( blen /= rows+1 ) call io_error('Error: Wrong number of lines in block '//trim(keyword))
    !       else
    !          call io_error('Error: Wrong number of lines in block '//trim(keyword))
    !       endif
    !    endif

    if ((blen .ne. rows) .and. (blen .ne. rows + 1) .and. (rows .gt. 0)) &
      call io_error('Error: Wrong number of lines in block '//trim(keyword), stdout, seedname)

    if ((blen .eq. rows + 1) .and. (rows .gt. 0) .and. (index(trim(keyword), 'unit_cell_cart') .eq. 0)) &
      call io_error('Error: Wrong number of lines in block '//trim(keyword), stdout, seedname)

    found = .true.

    lconvert = .false.
    if (blen == rows + 1) then
      dummy = in_data(line_s + 1)
      if (index(dummy, 'ang') .ne. 0) then
        lconvert = .false.
      elseif (index(dummy, 'bohr') .ne. 0) then
        lconvert = .true.
      else
        call io_error('Error: Units in block '//trim(keyword)//' not recognised', stdout, seedname)
      endif
      in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    endif

!    r_value=1.0_dp
    counter = 0
    do loop = line_s + 1, line_e - 1
      dummy = in_data(loop)
      counter = counter + 1
      if (present(c_value)) read (dummy, *, err=240, end=240) (c_value(i, counter), i=1, columns)
      if (present(l_value)) then
        ! I don't think we need this. Maybe read into a dummy charater
        ! array and convert each element to logical
        call io_error('param_get_keyword_block unimplemented for logicals', stdout, seedname)
      endif
      if (present(i_value)) read (dummy, *, err=240, end=240) (i_value(i, counter), i=1, columns)
      if (present(r_value)) read (dummy, *, err=240, end=240) (r_value(i, counter), i=1, columns)
    end do

    if (lconvert) then
      if (present(r_value)) then
        r_value = r_value*bohr
      endif
    endif

    in_data(line_s:line_e) (1:maxlen) = ' '

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword), stdout, seedname)

  end subroutine param_get_keyword_block

!=====================================================!
  subroutine param_get_block_length(stdout, seedname, keyword, found, rows, library, lunits)
    !=====================================================!
    !                                                     !
    !! Finds the length of the data block
    !                                                     !
    !=====================================================!

    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(out) :: rows
    !! Number of rows
    logical, intent(in) :: library
    logical, optional, intent(out) :: lunits
    !! Have we found a unit specification

    integer           :: i, in, ins, ine, loop, line_e, line_s
    logical           :: found_e, found_s
    character(len=maxlen) :: end_st, start_st, dummy
    character(len=2)  :: atsym
    real(kind=dp)     :: atpos(3)

    rows = 0
    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    if (.not. found_s) then
      found = .false.
      return
    end if

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    rows = line_e - line_s - 1

    found = .true.

    ! Ignore atoms_cart and atoms_frac blocks if running in library mode
    if (library) then
      if (trim(keyword) .eq. 'atoms_cart' .or. trim(keyword) .eq. 'atoms_frac') then
        in_data(line_s:line_e) (1:maxlen) = ' '
      endif
    endif

    if (present(lunits)) then
      dummy = in_data(line_s + 1)
      !       write(stdout,*) dummy
      !       write(stdout,*) trim(dummy)
      read (dummy, *, end=555) atsym, (atpos(i), i=1, 3)
      lunits = .false.
    endif

    if (rows <= 0) then !cope with empty blocks
      found = .false.
      in_data(line_s:line_e) (1:maxlen) = ' '
    end if

    return

555 lunits = .true.

    if (rows <= 1) then !cope with empty blocks
      found = .false.
      in_data(line_s:line_e) (1:maxlen) = ' '
    end if

    return

  end subroutine param_get_block_length

!===================================!
  subroutine param_get_atoms(atoms, library, lunits, real_lattice, recip_lattice, bohr, stdout, seedname)
    !===================================!
    !                                   !
    !!   Fills the atom data block
    !                                   !
    !===================================!

    !use w90_constants, only: bohr
    use w90_utility, only: utility_frac_to_cart, utility_cart_to_frac
    use w90_io, only: io_error
    implicit none

    type(atom_data_type), intent(inout) :: atoms
    integer, intent(in) :: stdout
    logical, intent(in) :: library
    logical, intent(in) :: lunits
    !! Do we expect a first line with the units
    real(kind=dp), intent(in) :: real_lattice(3, 3), recip_lattice(3, 3)
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname

    real(kind=dp)     :: atoms_pos_frac_tmp(3, atoms%num_atoms)
    real(kind=dp)     :: atoms_pos_cart_tmp(3, atoms%num_atoms)
    character(len=20) :: keyword
    integer           :: in, ins, ine, loop, i, line_e, line_s, counter
    integer           :: i_temp, loop2, max_sites, ierr, ic
    logical           :: found_e, found_s, found, frac
    character(len=maxlen) :: dummy, end_st, start_st
    character(len=maxlen) :: ctemp(atoms%num_atoms)
    character(len=maxlen) :: atoms_label_tmp(atoms%num_atoms)
    logical           :: lconvert

    keyword = "atoms_cart"
    frac = .false.
    call param_get_block_length(stdout, seedname, "atoms_frac", found, i_temp, library)
    if (found) then
      keyword = "atoms_frac"
      frac = .true.
    end if

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    lconvert = .false.
    if (lunits) then
      dummy = in_data(line_s + 1)
      if (index(dummy, 'ang') .ne. 0) then
        lconvert = .false.
      elseif (index(dummy, 'bohr') .ne. 0) then
        lconvert = .true.
      else
        call io_error('Error: Units in block atoms_cart not recognised in param_get_atoms', stdout, seedname)
      endif
      in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    endif

    counter = 0
    do loop = line_s + 1, line_e - 1
      dummy = in_data(loop)
      counter = counter + 1
      if (frac) then
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_frac_tmp(i, counter), i=1, 3)
      else
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_cart_tmp(i, counter), i=1, 3)
      end if
    end do

    if (lconvert) atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr

    in_data(line_s:line_e) (1:maxlen) = ' '

    if (frac) then
      do loop = 1, atoms%num_atoms
        call utility_frac_to_cart(atoms_pos_frac_tmp(:, loop), atoms_pos_cart_tmp(:, loop), real_lattice)
      end do
    else
      do loop = 1, atoms%num_atoms
        call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), atoms_pos_frac_tmp(:, loop), recip_lattice)
      end do
    end if

    ! Now we sort the data into the proper structures
    atoms%num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, atoms%num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          atoms%num_species = atoms%num_species + 1
          ctemp(atoms%num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atoms%species_num(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in param_get_atoms', stdout, seedname)
    allocate (atoms%label(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in param_get_atoms', stdout, seedname)
    allocate (atoms%symbol(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in param_get_atoms', stdout, seedname)
    atoms%species_num(:) = 0

    do loop = 1, atoms%num_species
      atoms%label(loop) = ctemp(loop)
      do loop2 = 1, atoms%num_atoms
        if (trim(atoms%label(loop)) == trim(atoms_label_tmp(loop2))) then
          atoms%species_num(loop) = atoms%species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atoms%species_num)
    allocate (atoms%pos_frac(3, max_sites, atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in param_get_atoms', stdout, seedname)
    allocate (atoms%pos_cart(3, max_sites, atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in param_get_atoms', stdout, seedname)

    do loop = 1, atoms%num_species
      counter = 0
      do loop2 = 1, atoms%num_atoms
        if (trim(atoms%label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          atoms%pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atoms%pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, atoms%num_species
      atoms%symbol(loop) (1:2) = atoms%label(loop) (1:2)
      ic = ichar(atoms%symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atoms%symbol(loop) (2:2) = ' '
    end do

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword), stdout, seedname)

  end subroutine param_get_atoms

!=====================================================!
  subroutine param_lib_set_atoms(atoms, atoms_label_tmp, atoms_pos_cart_tmp, recip_lattice, &
                                 stdout, seedname)
    !=====================================================!
    !                                                     !
    !!   Fills the atom data block during a library call
    !                                                     !
    !=====================================================!

    use w90_utility, only: utility_cart_to_frac, utility_lowercase
    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    type(atom_data_type), intent(inout) :: atoms
    character(len=*), intent(in) :: atoms_label_tmp(atoms%num_atoms)
    !! Atom labels
    real(kind=dp), intent(in)      :: atoms_pos_cart_tmp(3, atoms%num_atoms)
    !! Atom positions
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    character(len=50), intent(in)  :: seedname

    real(kind=dp)     :: atoms_pos_frac_tmp(3, atoms%num_atoms)
    integer           :: loop2, max_sites, ierr, ic, loop, counter
    character(len=maxlen) :: ctemp(atoms%num_atoms)
    character(len=maxlen) :: tmp_string

    do loop = 1, atoms%num_atoms
      call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), &
                                atoms_pos_frac_tmp(:, loop), recip_lattice)
    enddo

    ! Now we sort the data into the proper structures
    atoms%num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, atoms%num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          atoms%num_species = atoms%num_species + 1
          ctemp(atoms%num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atoms%species_num(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in param_lib_set_atoms', stdout, seedname)
    allocate (atoms%label(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in param_lib_set_atoms', stdout, seedname)
    allocate (atoms%symbol(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in param_lib_set_atoms', stdout, seedname)
    atoms%species_num(:) = 0

    do loop = 1, atoms%num_species
      atoms%label(loop) = ctemp(loop)
      do loop2 = 1, atoms%num_atoms
        if (trim(atoms%label(loop)) == trim(atoms_label_tmp(loop2))) then
          atoms%species_num(loop) = atoms%species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atoms%species_num)
    allocate (atoms%pos_frac(3, max_sites, atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in param_lib_set_atoms', stdout, seedname)
    allocate (atoms%pos_cart(3, max_sites, atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in param_lib_set_atoms', stdout, seedname)

    do loop = 1, atoms%num_species
      counter = 0
      do loop2 = 1, atoms%num_atoms
        if (trim(atoms%label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          atoms%pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atoms%pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, atoms%num_species
      atoms%symbol(loop) (1:2) = atoms%label(loop) (1:2)
      ic = ichar(atoms%symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atoms%symbol(loop) (2:2) = ' '
      tmp_string = trim(adjustl(utility_lowercase(atoms%symbol(loop))))
      atoms%symbol(loop) (1:2) = tmp_string(1:2)
      tmp_string = trim(adjustl(utility_lowercase(atoms%label(loop))))
      atoms%label(loop) (1:2) = tmp_string(1:2)
    end do

    return

  end subroutine param_lib_set_atoms

!====================================================================!
  subroutine param_get_range_vector(stdout, seedname, keyword, found, length, lcount, i_value)
    !====================================================================!
    !!   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100
    !!   if(lcount) we return the number of states in length
    !====================================================================!
    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)    :: keyword
    !! Keyword to examine
    logical, intent(out)   :: found
    !! Is keyword found
    integer, intent(inout) :: length
    !! Number of states
    logical, intent(in)    :: lcount
    !! If T only count states
    integer, optional, intent(out)   :: i_value(length)
    !! States specified in range vector

    integer   :: kl, in, loop, num1, num2, i_punc
    integer   :: counter, i_digit, loop_r, range_size
    character(len=maxlen) :: dummy
    character(len=10), parameter :: c_digit = "0123456789"
    character(len=2), parameter :: c_range = "-:"
    character(len=3), parameter :: c_sep = " ,;"
    character(len=5), parameter :: c_punc = " ,;-:"
    character(len=5)  :: c_num1, c_num2

    if (lcount .and. present(i_value)) call io_error('param_get_range_vector: incorrect call', stdout, seedname)

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      dummy = adjustl(dummy)
      if (.not. lcount) in_data(loop) (1:maxlen) = ' '
      if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
        dummy = dummy(2:)
        dummy = adjustl(dummy)
      end if
    end do

    if (.not. found) return

    counter = 0
    if (len_trim(dummy) == 0) call io_error('Error: keyword '//trim(keyword)//' is blank', stdout, seedname)
    dummy = adjustl(dummy)
    do
      i_punc = scan(dummy, c_punc)
      if (i_punc == 0) call io_error('Error parsing keyword '//trim(keyword), stdout, seedname)
      c_num1 = dummy(1:i_punc - 1)
      read (c_num1, *, err=101, end=101) num1
      dummy = adjustl(dummy(i_punc:))
      !look for range
      if (scan(dummy, c_range) == 1) then
        i_digit = scan(dummy, c_digit)
        dummy = adjustl(dummy(i_digit:))
        i_punc = scan(dummy, c_punc)
        c_num2 = dummy(1:i_punc - 1)
        read (c_num2, *, err=101, end=101) num2
        dummy = adjustl(dummy(i_punc:))
        range_size = abs(num2 - num1) + 1
        do loop_r = 1, range_size
          counter = counter + 1
          if (.not. lcount) i_value(counter) = min(num1, num2) + loop_r - 1
        end do
      else
        counter = counter + 1
        if (.not. lcount) i_value(counter) = num1
      end if

      if (scan(dummy, c_sep) == 1) dummy = adjustl(dummy(2:))
      if (scan(dummy, c_range) == 1) call io_error('Error parsing keyword '//trim(keyword)//' incorrect range', stdout, seedname)
      if (index(dummy, ' ') == 1) exit
    end do

    if (lcount) length = counter
    if (.not. lcount) then
      do loop = 1, counter - 1
        do loop_r = loop + 1, counter
          if (i_value(loop) == i_value(loop_r)) &
            call io_error('Error parsing keyword '//trim(keyword)//' duplicate values', stdout, seedname)
        end do
      end do
    end if

    return

101 call io_error('Error parsing keyword '//trim(keyword), stdout, seedname)

  end subroutine param_get_range_vector

  subroutine param_get_centre_constraints(ccentres_frac, ccentres_cart, &
                                          proj_site, num_wann, real_lattice, stdout, seedname)
    !=============================================================================!
    !                                                                             !
    !!  assigns projection centres as default centre constraints and global
    !!  Lagrange multiplier as individual Lagrange multipliers then reads
    !!  the centre_constraints block for individual centre constraint parameters
    !                                                                             !
    !=============================================================================!
    use w90_io, only: io_error
    use w90_utility, only: utility_frac_to_cart
    integer, intent(in) :: stdout
    real(kind=dp), intent(inout) :: ccentres_frac(:, :), ccentres_cart(:, :)
    real(kind=dp), intent(in) :: proj_site(:, :)
    integer, intent(in) :: num_wann
    character(len=50), intent(in)  :: seedname
    !type(param_wannierise_type), intent(inout) :: param_wannierise
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    integer           :: loop1, index1, constraint_num, loop2
    integer           :: column, start, finish, wann
    !logical           :: found
    character(len=maxlen) :: dummy

    do loop1 = 1, num_wann
      do loop2 = 1, 3
        ccentres_frac(loop1, loop2) = proj_site(loop2, loop1)
      end do
    end do

    constraint_num = 0
    do loop1 = 1, num_lines
      dummy = in_data(loop1)
      if (constraint_num > 0) then
        if (trim(dummy) == '') cycle
        index1 = index(dummy, 'begin')
        if (index1 > 0) call io_error("slwf_centres block hasn't ended yet", stdout, seedname)
        index1 = index(dummy, 'end')
        if (index1 > 0) then
          index1 = index(dummy, 'slwf_centres')
          if (index1 == 0) call io_error('Wrong ending of block (need to end slwf_centres)', stdout, seedname)
          in_data(loop1) (1:maxlen) = ' '
          exit
        end if
        column = 0
        start = 1
        finish = 1
        do loop2 = 1, len_trim(dummy)
          if (start == loop2 .and. dummy(loop2:loop2) == ' ') then
            start = loop2 + 1
          end if
          if (start < loop2) then
            if (dummy(loop2:loop2) == ' ') then
              finish = loop2 - 1
              call param_get_centre_constraint_from_column(column, start, finish, &
                                                           wann, dummy, ccentres_frac, stdout, seedname)
              start = loop2 + 1
              finish = start
            end if
          end if
          if (loop2 == len_trim(dummy) .and. dummy(loop2:loop2) /= ' ') then
            finish = loop2
            call param_get_centre_constraint_from_column(column, start, finish, &
                                                         wann, dummy, ccentres_frac, stdout, seedname)
            start = loop2 + 1
            finish = start
          end if
        end do
        in_data(loop1) (1:maxlen) = ' '
        constraint_num = constraint_num + 1
      end if
      index1 = index(dummy, 'slwf_centres')
      if (index1 > 0) then
        index1 = index(dummy, 'begin')
        if (index1 > 0) then
          constraint_num = 1
          in_data(loop1) (1:maxlen) = ' '
        end if
      end if
    end do
    do loop1 = 1, num_wann
      call utility_frac_to_cart(ccentres_frac(loop1, :), &
                                ccentres_cart(loop1, :), real_lattice)
    end do
  end subroutine param_get_centre_constraints

  subroutine param_get_centre_constraint_from_column(column, start, finish, &
                                                     wann, dummy, ccentres_frac, stdout, seedname)
    !===================================!
    !                                   !
    !!  assigns value read to constraint
    !!  parameters based on column
    !                                   !
    !===================================!
    use w90_io, only: io_error
    integer, intent(in) :: stdout
    integer, intent(inout):: column, start, finish, wann
    character(len=maxlen), intent(inout):: dummy
    character(len=50), intent(in)  :: seedname
    real(kind=dp), intent(inout) :: ccentres_frac(:, :)

    if (column == 0) then
      read (dummy(start:finish), '(i3)') wann
    end if
    if (column > 0) then
      if (column > 4) call io_error("Didn't expect anything else after Lagrange multiplier", stdout, seedname)
      if (column < 4) read (dummy(start:finish), '(f10.10)') ccentres_frac(wann, column)
    end if
    column = column + 1
  end subroutine param_get_centre_constraint_from_column

!===================================!
  subroutine param_get_projections(num_proj, atoms, num_wann, input_proj, proj, &
                                   recip_lattice, lcount, spinors, bohr, stdout, seedname)
    !===================================!
    !                                   !
    !!  Fills the projection data block
    !                                   !
    !===================================!

    use w90_constants, only: eps6, eps2
    use w90_utility, only: utility_cart_to_frac, &
      utility_string_to_coord, utility_strip
    use w90_io, only: io_error

    implicit none

    integer, intent(inout) :: num_proj
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(atom_data_type), intent(in) :: atoms
    ! projection data
    !type(param_driver_type), intent(inout) :: driver
    !type(param_wannierise_type), intent(inout) :: param_wannierise
    integer, intent(in) :: num_wann
    !type(select_projection_type), intent(inout) :: select_proj
    type(proj_input_type), intent(inout) :: proj ! intent(out)?
    type(proj_input_type), intent(inout) :: input_proj
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    logical, intent(in)    :: lcount
    logical, intent(in) :: spinors
    real(kind=dp), intent(in) :: bohr

    real(kind=dp)     :: pos_frac(3)
    real(kind=dp)     :: pos_cart(3)
    character(len=20) :: keyword
    integer           :: in, ins, ine, loop, line_e, line_s, counter
    integer           :: sites, species, line, pos1, pos2, pos3, m_tmp, l_tmp, mstate
    integer           :: loop_l, loop_m, loop_sites, ierr, loop_s, spn_counter
    logical           :: found_e, found_s
    character(len=maxlen) :: dummy, end_st, start_st
    character(len=maxlen) :: ctemp, ctemp2, ctemp3, ctemp4, ctemp5, m_string
    !
    integer, parameter :: min_l = -5
    integer, parameter :: max_l = 3
    integer, parameter :: min_m = 1
    integer, parameter :: max_m = 7
    integer            :: ang_states(min_m:max_m, min_l:max_l)
    ! default values for the optional part of the projection definitions
    real(kind=dp), parameter :: proj_z_def(3) = (/0.0_dp, 0.0_dp, 1.0_dp/)
    real(kind=dp), parameter :: proj_x_def(3) = (/1.0_dp, 0.0_dp, 0.0_dp/)
    real(kind=dp), parameter :: proj_s_qaxis_def(3) = (/0.0_dp, 0.0_dp, 1.0_dp/)
    real(kind=dp), parameter :: proj_zona_def = 1.0_dp
    integer, parameter       :: proj_radial_def = 1
    !
    real(kind=dp) :: proj_z_tmp(3)
    real(kind=dp) :: proj_x_tmp(3)
    real(kind=dp) :: proj_s_qaxis_tmp(3)
    real(kind=dp) :: proj_zona_tmp
    integer       :: proj_radial_tmp
    logical       :: lconvert, lrandom, proj_u_tmp, proj_d_tmp
    logical       :: lpartrandom
    !
    real(kind=dp) :: xnorm, znorm, cosphi, sinphi, xnorm_new, cosphi_new

    keyword = "projections"

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

!     if(spinors) num_proj=num_wann/2

    if (.not. lcount) then
      allocate (input_proj%site(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_site in param_get_projections', stdout, seedname)
      allocate (input_proj%l(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_l in param_get_projections', stdout, seedname)
      allocate (input_proj%m(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_m in param_get_projections', stdout, seedname)
      allocate (input_proj%z(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_z in param_get_projections', stdout, seedname)
      allocate (input_proj%x(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_x in param_get_projections', stdout, seedname)
      allocate (input_proj%radial(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_radial in param_get_projections', stdout, seedname)
      allocate (input_proj%zona(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_zona in param_get_projections', stdout, seedname)
      if (spinors) then
        allocate (input_proj%s(num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s in param_get_projections', stdout, seedname)
        allocate (input_proj%s_qaxis(3, num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s_qaxis in param_get_projections', stdout, seedname)
      endif

      allocate (proj%site(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_site in param_get_projections', stdout, seedname)
      allocate (proj%l(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_l in param_get_projections', stdout, seedname)
      allocate (proj%m(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_m in param_get_projections', stdout, seedname)
      allocate (proj%z(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_z in param_get_projections', stdout, seedname)
      allocate (proj%x(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_x in param_get_projections', stdout, seedname)
      allocate (proj%radial(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_radial in param_get_projections', stdout, seedname)
      allocate (proj%zona(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_zona in param_get_projections', stdout, seedname)
      if (spinors) then
        allocate (proj%s(num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_s in param_get_projections', stdout, seedname)
        allocate (proj%s_qaxis(3, num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_s_qaxis in param_get_projections', stdout, seedname)
      endif
    endif

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('param_get_projections: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('param_get_projections: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('param_get_projections: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    dummy = in_data(line_s + 1)
    lconvert = .false.
    lrandom = .false.
    lpartrandom = .false.
    if (index(dummy, 'ang') .ne. 0) then
      if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    elseif (index(dummy, 'bohr') .ne. 0) then
      if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
      lconvert = .true.
    elseif (index(dummy, 'random') .ne. 0) then
      if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
      if (index(in_data(line_s + 1), end_st) .ne. 0) then
        lrandom = .true.     ! all projections random
      else
        lpartrandom = .true. ! only some projections random
        if (index(in_data(line_s + 1), 'ang') .ne. 0) then
          if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
          line_s = line_s + 1
        elseif (index(in_data(line_s + 1), 'bohr') .ne. 0) then
          if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
          line_s = line_s + 1
          lconvert = .true.
        endif
      endif
    endif

    counter = 0
    if (.not. lrandom) then
      do line = line_s + 1, line_e - 1
        ang_states = 0
        !Assume the default values
        proj_z_tmp = proj_z_def
        proj_x_tmp = proj_x_def
        proj_zona_tmp = proj_zona_def
        proj_radial_tmp = proj_radial_def
        if (spinors) then
          proj_s_qaxis_tmp = proj_s_qaxis_def
          spn_counter = 2
          proj_u_tmp = .true.
          proj_d_tmp = .true.
        else
          spn_counter = 1
        endif
        ! Strip input line of all spaces
        dummy = utility_strip(in_data(line))
        dummy = adjustl(dummy)
        pos1 = index(dummy, ':')
        if (pos1 == 0) &
          call io_error('param_read_projection: malformed projection definition: '//trim(dummy), stdout, seedname)
        sites = 0
        ctemp = dummy(:pos1 - 1)
        ! Read the atomic site
        if (index(ctemp, 'c=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_cart, stdout, seedname)
          if (lconvert) pos_cart = pos_cart*bohr
          call utility_cart_to_frac(pos_cart(:), pos_frac(:), recip_lattice)
        elseif (index(ctemp, 'f=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_frac, stdout, seedname)
        else
          if (atoms%num_species == 0) &
            call io_error('param_get_projection: Atom centred projection requested but no atoms defined', stdout, seedname)
          do loop = 1, atoms%num_species
            if (trim(ctemp) == atoms%label(loop)) then
              species = loop
              sites = atoms%species_num(loop)
              exit
            end if
            if (loop == atoms%num_species) then
              call io_error('param_get_projection: Atom site not recognised '//trim(ctemp), &
                            stdout, seedname)
            endif
          end do
        end if

        dummy = dummy(pos1 + 1:)

        ! scan for quantisation direction
        pos1 = index(dummy, '[')
        if (spinors) then
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 1:))
            pos2 = index(ctemp, ']')
            if (pos2 == 0) call io_error &
              ('param_get_projections: no closing square bracket for spin quantisation dir', stdout, seedname)
            ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_s_qaxis_tmp, stdout, seedname)
            dummy = dummy(:pos1 - 1) ! remove [ ] section
          endif
        else
          if (pos1 > 0) call io_error('param_get_projections: spin qdir is defined but spinors=.false.', stdout, seedname)
        endif

        ! scan for up or down
        pos1 = index(dummy, '(')
        if (spinors) then
          if (pos1 > 0) then
            proj_u_tmp = .false.; proj_d_tmp = .false.
            ctemp = (dummy(pos1 + 1:))
            pos2 = index(ctemp, ')')
            if (pos2 == 0) call io_error('param_get_projections: no closing bracket for spin', stdout, seedname)
            ctemp = ctemp(:pos2 - 1)
            if (index(ctemp, 'u') > 0) proj_u_tmp = .true.
            if (index(ctemp, 'd') > 0) proj_d_tmp = .true.
            if (proj_u_tmp .and. proj_d_tmp) then
              spn_counter = 2
            elseif (.not. proj_u_tmp .and. .not. proj_d_tmp) then
              call io_error('param_get_projections: found brackets but neither u or d', stdout, seedname)
            else
              spn_counter = 1
            endif
            dummy = dummy(:pos1 - 1) ! remove ( ) section
          endif
        else
          if (pos1 > 0) call io_error('param_get_projections: spin is defined but spinors=.false.', stdout, seedname)
        endif

        !Now we know the sites for this line. Get the angular momentum states
        pos1 = index(dummy, ':')
        if (pos1 > 0) then
          ctemp = dummy(:pos1 - 1)
        else
          ctemp = dummy
        end if
        ctemp2 = ctemp
        do
          pos2 = index(ctemp2, ';')
          if (pos2 == 0) then
            ctemp3 = ctemp2
          else
            ctemp3 = ctemp2(:pos2 - 1)
          endif
          if (index(ctemp3, 'l=') == 1) then
            mstate = index(ctemp3, ',')
            if (mstate > 0) then
              read (ctemp3(3:mstate - 1), *, err=101, end=101) l_tmp
            else
              read (ctemp3(3:), *, err=101, end=101) l_tmp
            end if
            if (l_tmp < -5 .or. l_tmp > 3) call io_error('param_get_projection: Incorrect l state requested', stdout, seedname)
            if (mstate == 0) then
              if (l_tmp >= 0) then
                do loop_m = 1, 2*l_tmp + 1
                  ang_states(loop_m, l_tmp) = 1
                end do
              elseif (l_tmp == -1) then !sp
                ang_states(1:2, l_tmp) = 1
              elseif (l_tmp == -2) then !sp2
                ang_states(1:3, l_tmp) = 1
              elseif (l_tmp == -3) then !sp3
                ang_states(1:4, l_tmp) = 1
              elseif (l_tmp == -4) then !sp3d
                ang_states(1:5, l_tmp) = 1
              elseif (l_tmp == -5) then !sp3d2
                ang_states(1:6, l_tmp) = 1
              endif
            else
              if (index(ctemp3, 'mr=') /= mstate + 1) &
                call io_error('param_get_projection: Problem reading m state', stdout, seedname)
              ctemp4 = ctemp3(mstate + 4:)
              do
                pos3 = index(ctemp4, ',')
                if (pos3 == 0) then
                  ctemp5 = ctemp4
                else
                  ctemp5 = ctemp4(:pos3 - 1)
                endif
                read (ctemp5(1:), *, err=102, end=102) m_tmp
                if (l_tmp >= 0) then
                  if ((m_tmp > 2*l_tmp + 1) .or. (m_tmp <= 0)) call io_error('param_get_projection: m is > l !', stdout, seedname)
                elseif (l_tmp == -1 .and. (m_tmp > 2 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (1)', stdout, seedname)
                elseif (l_tmp == -2 .and. (m_tmp > 3 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (2)', stdout, seedname)
                elseif (l_tmp == -3 .and. (m_tmp > 4 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (3)', stdout, seedname)
                elseif (l_tmp == -4 .and. (m_tmp > 5 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (4)', stdout, seedname)
                elseif (l_tmp == -5 .and. (m_tmp > 6 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (5)', stdout, seedname)
                endif
                ang_states(m_tmp, l_tmp) = 1
                if (pos3 == 0) exit
                ctemp4 = ctemp4(pos3 + 1:)
              enddo
            end if
          else
            do
              pos3 = index(ctemp3, ',')
              if (pos3 == 0) then
                ctemp4 = ctemp3
              else
                ctemp4 = ctemp3(:pos3 - 1)
              endif
              read (ctemp4(1:), *, err=106, end=106) m_string
              select case (trim(adjustl(m_string)))
              case ('s')
                ang_states(1, 0) = 1
              case ('p')
                ang_states(1:3, 1) = 1
              case ('pz')
                ang_states(1, 1) = 1
              case ('px')
                ang_states(2, 1) = 1
              case ('py')
                ang_states(3, 1) = 1
              case ('d')
                ang_states(1:5, 2) = 1
              case ('dz2')
                ang_states(1, 2) = 1
              case ('dxz')
                ang_states(2, 2) = 1
              case ('dyz')
                ang_states(3, 2) = 1
              case ('dx2-y2')
                ang_states(4, 2) = 1
              case ('dxy')
                ang_states(5, 2) = 1
              case ('f')
                ang_states(1:7, 3) = 1
              case ('fz3')
                ang_states(1, 3) = 1
              case ('fxz2')
                ang_states(2, 3) = 1
              case ('fyz2')
                ang_states(3, 3) = 1
              case ('fxyz')
                ang_states(4, 3) = 1
              case ('fz(x2-y2)')
                ang_states(5, 3) = 1
              case ('fx(x2-3y2)')
                ang_states(6, 3) = 1
              case ('fy(3x2-y2)')
                ang_states(7, 3) = 1
              case ('sp')
                ang_states(1:2, -1) = 1
              case ('sp-1')
                ang_states(1, -1) = 1
              case ('sp-2')
                ang_states(2, -1) = 1
              case ('sp2')
                ang_states(1:3, -2) = 1
              case ('sp2-1')
                ang_states(1, -2) = 1
              case ('sp2-2')
                ang_states(2, -2) = 1
              case ('sp2-3')
                ang_states(3, -2) = 1
              case ('sp3')
                ang_states(1:4, -3) = 1
              case ('sp3-1')
                ang_states(1, -3) = 1
              case ('sp3-2')
                ang_states(2, -3) = 1
              case ('sp3-3')
                ang_states(3, -3) = 1
              case ('sp3-4')
                ang_states(4, -3) = 1
              case ('sp3d')
                ang_states(1:5, -4) = 1
              case ('sp3d-1')
                ang_states(1, -4) = 1
              case ('sp3d-2')
                ang_states(2, -4) = 1
              case ('sp3d-3')
                ang_states(3, -4) = 1
              case ('sp3d-4')
                ang_states(4, -4) = 1
              case ('sp3d-5')
                ang_states(5, -4) = 1
              case ('sp3d2')
                ang_states(1:6, -5) = 1
              case ('sp3d2-1')
                ang_states(1, -5) = 1
              case ('sp3d2-2')
                ang_states(2, -5) = 1
              case ('sp3d2-3')
                ang_states(3, -5) = 1
              case ('sp3d2-4')
                ang_states(4, -5) = 1
              case ('sp3d2-5')
                ang_states(5, -5) = 1
              case ('sp3d2-6')
                ang_states(6, -5) = 1
              case default
                call io_error('param_get_projection: Problem reading l state '//trim(ctemp3), stdout, seedname)
              end select
              if (pos3 == 0) exit
              ctemp3 = ctemp3(pos3 + 1:)
            enddo
          endif
          if (pos2 == 0) exit
          ctemp2 = ctemp2(pos2 + 1:)
        enddo
        ! check for non-default values
        if (pos1 > 0) then
          dummy = dummy(pos1 + 1:)
          ! z axis
          pos1 = index(dummy, 'z=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_z_tmp, stdout, seedname)
          endif
          ! x axis
          pos1 = index(dummy, 'x=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_x_tmp, stdout, seedname)
          endif
          ! diffusivity of orbital
          pos1 = index(dummy, 'zona=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 5:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            read (ctemp, *, err=104, end=104) proj_zona_tmp
          endif
          ! nodes for the radial part
          pos1 = index(dummy, 'r=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            read (ctemp, *, err=105, end=105) proj_radial_tmp
          endif
        end if
        ! if (sites == -1) then
        !   if (counter + spn_counter*sum(ang_states) > num_proj) &
        !     call io_error('param_get_projection: too many projections defined')
        ! else
        !   if (counter + spn_counter*sites*sum(ang_states) > num_proj) &
        !     call io_error('param_get_projection: too many projections defined')
        ! end if
        !
        if (sites == -1) then
          do loop_l = min_l, max_l
            do loop_m = min_m, max_m
              if (ang_states(loop_m, loop_l) == 1) then
                do loop_s = 1, spn_counter
                  counter = counter + 1
                  if (lcount) cycle
                  input_proj%site(:, counter) = pos_frac
                  input_proj%l(counter) = loop_l
                  input_proj%m(counter) = loop_m
                  input_proj%z(:, counter) = proj_z_tmp
                  input_proj%x(:, counter) = proj_x_tmp
                  input_proj%radial(counter) = proj_radial_tmp
                  input_proj%zona(counter) = proj_zona_tmp
                  if (spinors) then
                    if (spn_counter == 1) then
                      if (proj_u_tmp) input_proj%s(counter) = 1
                      if (proj_d_tmp) input_proj%s(counter) = -1
                    else
                      if (loop_s == 1) input_proj%s(counter) = 1
                      if (loop_s == 2) input_proj%s(counter) = -1
                    endif
                    input_proj%s_qaxis(:, counter) = proj_s_qaxis_tmp
                  endif
                end do
              endif
            end do
          end do
        else
          do loop_sites = 1, sites
            do loop_l = min_l, max_l
              do loop_m = min_m, max_m
                if (ang_states(loop_m, loop_l) == 1) then
                  do loop_s = 1, spn_counter
                    counter = counter + 1
                    if (lcount) cycle
                    input_proj%site(:, counter) = atoms%pos_frac(:, loop_sites, species)
                    input_proj%l(counter) = loop_l
                    input_proj%m(counter) = loop_m
                    input_proj%z(:, counter) = proj_z_tmp
                    input_proj%x(:, counter) = proj_x_tmp
                    input_proj%radial(counter) = proj_radial_tmp
                    input_proj%zona(counter) = proj_zona_tmp
                    if (spinors) then
                      if (spn_counter == 1) then
                        if (proj_u_tmp) input_proj%s(counter) = 1
                        if (proj_d_tmp) input_proj%s(counter) = -1
                      else
                        if (loop_s == 1) input_proj%s(counter) = 1
                        if (loop_s == 2) input_proj%s(counter) = -1
                      endif
                      input_proj%s_qaxis(:, counter) = proj_s_qaxis_tmp
                    endif
                  end do
                end if
              end do
            end do
          end do
        end if

      end do !end loop over projection block

      ! check there are enough projections and add random projections if required
      if (.not. lpartrandom) then
        if (counter .lt. num_wann) call io_error( &
          'param_get_projections: too few projection functions defined', stdout, seedname)
      end if
    end if ! .not. lrandom

    if (lcount) then
      if (counter .lt. num_wann) then
        num_proj = num_wann
      else
        num_proj = counter
      endif
      return
    endif

    if (lpartrandom .or. lrandom) then
      call random_seed()  ! comment out this line for reproducible random positions!
      do loop = counter + 1, num_wann
        call random_number(input_proj%site(:, loop))
        input_proj%l(loop) = 0
        input_proj%m(loop) = 1
        input_proj%z(:, loop) = proj_z_def
        input_proj%x(:, loop) = proj_x_def
        input_proj%zona(loop) = proj_zona_def
        input_proj%radial(loop) = proj_radial_def
        if (spinors) then
          if (modulo(loop, 2) == 1) then
            input_proj%s(loop) = 1
          else
            input_proj%s(loop) = -1
          end if
          input_proj%s_qaxis(1, loop) = 0.
          input_proj%s_qaxis(2, loop) = 0.
          input_proj%s_qaxis(3, loop) = 1.
        end if
      enddo
    endif

    ! I shouldn't get here, but just in case
    if (.not. lcount) in_data(line_s:line_e) (1:maxlen) = ' '

!~     ! Check
!~     do loop=1,num_wann
!~        if ( abs(sum(proj_z(:,loop)*proj_x(:,loop))).gt.1.0e-6_dp ) then
!~           write(stdout,*) ' Projection:',loop
!~           call io_error(' Error in projections: z and x axes are not orthogonal')
!~        endif
!~     enddo

    ! Normalise z-axis and x-axis and check/fix orthogonality
    do loop = 1, num_proj

      znorm = sqrt(sum(input_proj%z(:, loop)*input_proj%z(:, loop)))
      xnorm = sqrt(sum(input_proj%x(:, loop)*input_proj%x(:, loop)))
      input_proj%z(:, loop) = input_proj%z(:, loop)/znorm             ! normalise z
      input_proj%x(:, loop) = input_proj%x(:, loop)/xnorm             ! normalise x
      cosphi = sum(input_proj%z(:, loop)*input_proj%x(:, loop))

      ! Check whether z-axis and z-axis are orthogonal
      if (abs(cosphi) .gt. eps6) then

        ! Special case of circularly symmetric projections (pz, dz2, fz3)
        ! just choose an x-axis that is perpendicular to the given z-axis
        if ((input_proj%l(loop) .ge. 0) .and. (input_proj%m(loop) .eq. 1)) then
          proj_x_tmp(:) = input_proj%x(:, loop)            ! copy of original x-axis
          call random_seed()
          call random_number(proj_z_tmp(:))         ! random vector
          ! calculate new x-axis as the cross (vector) product of random vector with z-axis
          input_proj%x(1, loop) = proj_z_tmp(2)*input_proj%z(3, loop) &
                                  - proj_z_tmp(3)*input_proj%z(2, loop)
          input_proj%x(2, loop) = proj_z_tmp(3)*input_proj%z(1, loop) &
                                  - proj_z_tmp(1)*input_proj%z(3, loop)
          input_proj%x(3, loop) = proj_z_tmp(1)*input_proj%z(2, loop) &
                                  - proj_z_tmp(2)*input_proj%z(1, loop)
          xnorm_new = sqrt(sum(input_proj%x(:, loop)*input_proj%x(:, loop)))
          input_proj%x(:, loop) = input_proj%x(:, loop)/xnorm_new   ! normalise
          goto 555
        endif

        ! If projection axes non-orthogonal enough, then
        ! user may have made a mistake and should check
        if (abs(cosphi) .gt. eps2) then
          write (stdout, *) ' Projection:', loop
          call io_error(' Error in projections: z and x axes are not orthogonal', stdout, seedname)
        endif

        ! If projection axes are "reasonably orthogonal", project x-axis
        ! onto plane perpendicular to z-axis to make them more so
        sinphi = sqrt(1 - cosphi*cosphi)
        proj_x_tmp(:) = input_proj%x(:, loop)               ! copy of original x-axis
        ! calculate new x-axis:
        ! x = z \cross (x_tmp \cross z) / sinphi = ( x_tmp - z(z.x_tmp) ) / sinphi
        input_proj%x(:, loop) = (proj_x_tmp(:) - cosphi*input_proj%z(:, loop))/sinphi

        ! Final check
555     cosphi_new = sum(input_proj%z(:, loop)*input_proj%x(:, loop))
        if (abs(cosphi_new) .gt. eps6) then
          write (stdout, *) ' Projection:', loop
          call io_error(' Error: z and x axes are still not orthogonal after projection', stdout, seedname)
        endif

      endif

    enddo

    return

101 call io_error('param_get_projection: Problem reading l state into integer '//trim(ctemp3), stdout, seedname)
102 call io_error('param_get_projection: Problem reading m state into integer '//trim(ctemp3), stdout, seedname)
104 call io_error('param_get_projection: Problem reading zona into real '//trim(ctemp), stdout, seedname)
105 call io_error('param_get_projection: Problem reading radial state into integer '//trim(ctemp), stdout, seedname)
106 call io_error('param_get_projection: Problem reading m state into string '//trim(ctemp3), stdout, seedname)

  end subroutine param_get_projections

!===================================!
  subroutine param_get_keyword_kpath(spec_points, stdout, seedname)
    !===================================!
    !                                   !
    !!  Fills the kpath data block
    !                                   !
    !===================================!
    use w90_io, only: io_error

    implicit none

    type(kpoint_path_type), intent(inout) :: spec_points
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    character(len=20) :: keyword
    integer           :: in, ins, ine, loop, i, line_e, line_s, counter
    logical           :: found_e, found_s
    character(len=maxlen) :: dummy, end_st, start_st

    keyword = "kpoint_path"

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    counter = 0
    do loop = line_s + 1, line_e - 1

      counter = counter + 2
      dummy = in_data(loop)
      read (dummy, *, err=240, end=240) spec_points%labels(counter - 1), &
        (spec_points%points(i, counter - 1), i=1, 3), &
        spec_points%labels(counter), (spec_points%points(i, counter), i=1, 3)
    end do

    in_data(line_s:line_e) (1:maxlen) = ' '

    return

240 call io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy), stdout, seedname)

  end subroutine param_get_keyword_kpath

  subroutine param_clear_block(stdout, seedname, keyword)
    ! a dummy read routine to remove unused but legitimate input block from input stream
    ! needed to preserve input file error checking (i.e. input stream should be empty after all
    ! legitimate keywords/blocks are read)
    use w90_io, only: io_error

    implicit none

    ! passed variables
    integer, intent(in) :: stdout
    character(len=50), intent(in) :: seedname
    character(len=*), intent(in) :: keyword

    ! local variables
    integer :: in, ins, ine, loop, line_e, line_s
    logical :: found_e, found_s
    character(len=maxlen) :: end_st, start_st

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (found_s .and. (.not. found_e)) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (found_e .and. (.not. found_s)) then
      call io_error('Error: Found '//trim(end_st)//' but no '//trim(start_st)//' in input file', stdout, seedname)
    end if

    if (found_s .and. found_e) then
      if (line_e <= line_s) then
        call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
      end if

      in_data(line_s:line_e) (1:maxlen) = ' '  ! clear the block from the input stream
    end if ! found tags
  end subroutine param_clear_block

end module w90_param_methods
