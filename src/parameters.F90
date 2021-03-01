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

module w90_parameters
  !! This module contains parameters to control the actions of wannier90.
  !! Also routines to read the parameters and write them out again.

  use w90_constants, only: dp
  use w90_io, only: maxlen
  use w90_types, only: projection_type

  implicit none

  public

  ! GP: added a flag to check if this is the first run of param_read in library mode or not
  logical, save :: library_param_read_first_pass !BGS wannier_lib only
  !BGS flag this flag should eventually be removed, or put in driver_type?

  !Input
  type parameter_input_type
    integer :: iprint
    !! Controls the verbosity of the output
    integer :: timing_level
    character(len=20) :: length_unit
    !! Units for length
    character(len=50) :: devel_flag !kmesh, disentangle, postw90/postw90_common
    integer :: num_valence_bands !wannierise, postw90/postw90_common, get_oper and berry
    integer, allocatable :: exclude_bands(:) !kmesh, wannier_lib, w90chk2chk
    integer :: num_exclude_bands
    logical :: gamma_only !overlap, kmesh, disentangle, wannierise, wannier_prog
    !! Use the special Gamma-point routines
    character(len=20) :: bands_plot_mode !hamiltonian (setup only), plot
    !BGS - maybe a band_plot_type with band_num_points etc from plot_type?
    integer :: one_dim_dir ! transport and plot
    real(kind=dp) :: hr_cutoff !plot and transport
    real(kind=dp) :: dist_cutoff !plot and transport
    character(len=20) :: dist_cutoff_mode !plot and transport
    real(kind=dp) :: dist_cutoff_hc !plot and transport
    logical :: use_ws_distance !ws_distance, plot and postw90_common
    real(kind=dp) :: ws_distance_tol !ws_distance, hamiltonian and postw90_common
    !! absolute tolerance for the distance to equivalent positions
    integer :: ws_search_size(3) ! ws_distance, hamiltonian
    !! maximum extension in each direction of the supercell of the BvK cell
    !! to search for points inside the Wigner-Seitz cell
    logical :: spinors   !are our WF spinors? !kmesh, plot, wannier_lib, postw90/gyrotropic
    integer :: num_elec_per_state !wannierise and postw90 dos and boltzwann
    logical :: write_xyz !wannierise and transport
    integer :: optimisation !wannierise and disentangle
    real(kind=dp) :: omega_invariant !wannierise, disentangle and chk2chk
    logical :: have_disentangled !disentangle, plot, wannierise, postw90...
    real(kind=dp) :: lenconfac !lots of write statements in wannier90
  end type parameter_input_type
  type(parameter_input_type), save :: param_input

  ! setup in wannierise, but used by plot, ws_distance etc
  type wannier_data_type
    ! Wannier centres and spreads
    real(kind=dp), allocatable :: centres(:, :)
    real(kind=dp), allocatable :: spreads(:)
  end type wannier_data_type
  type(wannier_data_type), save :: wann_data

  ! used in kmesh, and to allocate in parameters
  ! The maximum number of shells we need to satisfy B1 condition in kmesh
  integer, parameter :: max_shells = 6
  integer, parameter :: num_nnmax = 12
  type param_kmesh_type
    integer :: num_shells
    !! no longer an input keyword
    logical :: skip_B1_tests
    !! do not check the B1 condition
    integer, allocatable :: shell_list(:)
    integer :: search_shells
    real(kind=dp) :: tol
    ! Projections, mainly used in parameters, maybe written by kmesh !BGS move?
    real(kind=dp), allocatable :: input_proj_site(:, :)
    type(projection_type) :: input_proj
    ! a u t o m a t i c   p r o j e c t i o n s
    ! vv: Writes a new block in .nnkp
    logical :: auto_projections
  end type param_kmesh_type
  type(param_kmesh_type), save :: kmesh_data

  ! kmesh parameters (set in kmesh)
  type kmesh_info_type
    integer              :: nnh           ! the number of b-directions (bka)
    integer              :: nntot         ! total number of neighbours for each k-point
    integer, allocatable :: nnlist(:, :)   ! list of neighbours for each k-point
    integer, allocatable :: neigh(:, :)
    integer, allocatable :: nncell(:, :, :) ! gives BZ of each neighbour of each k-point
    real(kind=dp)              :: wbtot
    real(kind=dp), allocatable :: wb(:)         ! weights associated with neighbours of each k-point
    real(kind=dp), allocatable :: bk(:, :, :)     ! the b-vectors that go from each k-point to its neighbours
    real(kind=dp), allocatable :: bka(:, :)      ! the b-directions from 1st k-point to its neighbours
  end type kmesh_info_type
  type(kmesh_info_type), save :: kmesh_info

  ! used in wannierise, hamiltonian, plot and others (postw90 also)
  type k_point_type
    real(kind=dp), allocatable :: kpt_latt(:, :) !! kpoints in lattice vecs
    real(kind=dp), allocatable :: kpt_cart(:, :) !kpoints in cartesians - kmesh and transport
  end type k_point_type
  type(k_point_type), save :: k_points
  integer, save :: num_kpts !BGS put in k_point_type?

  ! postw90/boltzwann use win_min/max, so maybe should move these?
  type disentangle_type
    real(kind=dp) :: win_min
    !! lower bound of the disentanglement outer window
    real(kind=dp) :: win_max
    !! upper bound of the disentanglement outer window
    real(kind=dp) :: froz_min
    !! lower bound of the disentanglement inner (frozen) window
    real(kind=dp) :: froz_max
    !! upper bound of the disentanglement inner (frozen) window
    integer :: num_iter
    !! number of disentanglement iteration steps
    real(kind=dp) :: mix_ratio
    !! Mixing ratio for the disentanglement routine
    real(kind=dp) :: conv_tol
    !! Convergence tolerance for the disentanglement
    integer :: conv_window
    !! Size of the convergence window for disentanglement
    ! GS-start
    integer :: spheres_first_wann
    integer :: spheres_num
    real(kind=dp), allocatable :: spheres(:, :)
    ! GS-end
    logical :: frozen_states
    ! disentangle parameters
    ! BGS also used by plot, hamiltonian, wannierise, postw90_common, get_oper
    integer, allocatable :: ndimwin(:)
    logical, allocatable :: lwindow(:, :)
  end type disentangle_type
  type(disentangle_type), save :: dis_data

  ! plot, transport, postw90: common, wan_ham, spin, berry, gyrotropic, kpath, kslice
  type fermi_data_type
    integer :: n
    real(kind=dp), allocatable :: energy_list(:)
  end type fermi_data_type
  type(fermi_data_type), save :: fermi

  ! Atom sites - often used in the write_* routines
  ! hamiltonian, wannierise, plot, transport, wannier_lib
  type atom_data_type
    real(kind=dp), allocatable :: pos_frac(:, :, :)
    real(kind=dp), allocatable :: pos_cart(:, :, :)
    integer, allocatable :: species_num(:)
    character(len=maxlen), allocatable :: label(:)
    character(len=2), allocatable :: symbol(:)
    integer :: num_atoms
    integer :: num_species
  end type atom_data_type
  type(atom_data_type), save :: atoms

  integer, save :: num_bands
  !! Number of bands

  integer, save :: num_wann
  !! number of wannier functions

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states
  !BGS disentangle, hamiltonian, a wannierise print, and postw90/get_oper
  real(kind=dp), allocatable, save :: eigval(:, :)

  !BGS u_matrix_opt in postw90 only for generation of v_matrix
  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states
  complex(kind=dp), allocatable, save :: u_matrix_opt(:, :, :)

  ! u_matrix gives the unitary rotations from the optimal subspace to the
  ! optimally smooth states.
  ! m_matrix we store here, becuase it is needed for restart of wannierise
  complex(kind=dp), allocatable, save :: u_matrix(:, :, :)
  !BGS is disentangle and overlap
  complex(kind=dp), allocatable, save :: m_matrix_local(:, :, :, :)

  integer, save :: mp_grid(3)
  !! Dimensions of the Monkhorst-Pack grid

  integer, save :: num_proj
  !BGS used by stuff in driver/kmesh/wannier - keep separate or duplicate?

  real(kind=dp), save :: real_lattice(3, 3)

  !parameters derived from input
  real(kind=dp), save :: recip_lattice(3, 3)

  ! plot.F90 and postw90/kpath
  type special_kpoints_type
    integer :: bands_num_spec_points
    character(len=20), allocatable ::bands_label(:)
    real(kind=dp), allocatable ::bands_spec_points(:, :)
  end type special_kpoints_type
  type(special_kpoints_type), save :: spec_points

end module w90_parameters

module w90_param_methods
  ! very few of these use save, so may actually be local to subroutines

  use w90_constants, only: dp
  use w90_io, only: stdout, maxlen
  use w90_parameters
  use wannier_parameters

  implicit none

  private

  ! Adaptive vs. fixed smearing stuff [GP, Jul 12, 2012]
  ! Only internal, always use the local variables defined by each module
  ! that take this value as default
  logical                         :: adpt_smr
  real(kind=dp)                   :: adpt_smr_fac
  real(kind=dp)                   :: adpt_smr_max
  real(kind=dp)                   :: smr_fixed_en_width

  real(kind=dp), save :: fermi_energy

  logical                                  :: fermi_energy_scan
  real(kind=dp)                            :: fermi_energy_min
  real(kind=dp)                            :: fermi_energy_max
  real(kind=dp)                            :: fermi_energy_step

  ! extra vars identified as private
  !logical, save :: berry_uHu_formatted
  !! Read the uHu from fortran formatted file

  ! Private data
  integer                            :: num_lines
  character(len=maxlen), allocatable :: in_data(:)
  character(len=maxlen)              :: ctmp
  logical                            :: ltmp

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
  ! common read routines
  public :: param_read_03
  public :: param_read_05
  public :: param_read_09
  public :: param_w90_read_10
  public :: param_read_11
  public :: param_read_13
  public :: param_read_16
  public :: param_w90_read_18a
  public :: param_read_21
  public :: param_read_23
  public :: param_read_25
  public :: param_read_28
  public :: param_read_32
  public :: param_w90_read_33 ! both
  public :: param_read_40a
  public :: param_read_40c
  public :: param_w90_read_41
  public :: param_read_44
  public :: param_w90_read_45

contains

  subroutine param_read_03
    !%%%%%%%%%%%%%%%%
    !System variables
    !%%%%%%%%%%%%%%%%
    implicit none
    logical :: found

    param_input%timing_level = 1             ! Verbosity of timing output info
    call param_get_keyword('timing_level', found, i_value=param_input%timing_level)

    param_input%iprint = 1             ! Verbosity
    call param_get_keyword('iprint', found, i_value=param_input%iprint)

    param_input%optimisation = 3             ! Verbosity
    call param_get_keyword('optimisation', found, i_value=param_input%optimisation)

  end subroutine param_read_03

  subroutine param_read_05(energy_unit)
    use w90_constants, only: bohr
    use w90_io, only: io_error
    implicit none
    character(len=*), intent(out) :: energy_unit
    logical :: found

    energy_unit = 'ev'          !
    call param_get_keyword('energy_unit', found, c_value=energy_unit)

    param_input%length_unit = 'ang'         !
    param_input%lenconfac = 1.0_dp
    call param_get_keyword('length_unit', found, c_value=param_input%length_unit)
    if (param_input%length_unit .ne. 'ang' .and. param_input%length_unit .ne. 'bohr') &
      call io_error('Error: value of length_unit not recognised in param_read')
    if (param_input%length_unit .eq. 'bohr') param_input%lenconfac = 1.0_dp/bohr
  end subroutine param_read_05

  subroutine param_read_09
    use w90_io, only: io_error
    implicit none
    logical :: found

    num_wann = -99
    call param_get_keyword('num_wann', found, i_value=num_wann)
    if (.not. found) call io_error('Error: You must specify num_wann')
    if (num_wann <= 0) call io_error('Error: num_wann must be greater than zero')
  end subroutine param_read_09

  subroutine param_w90_read_10
    use w90_io, only: io_error
    implicit none
    integer :: ierr
    logical :: found

    param_input%num_exclude_bands = 0
    call param_get_range_vector('exclude_bands', found, param_input%num_exclude_bands, lcount=.true.)
    if (found) then
      if (param_input%num_exclude_bands < 1) call io_error('Error: problem reading exclude_bands')
      if (allocated(param_input%exclude_bands)) deallocate (param_input%exclude_bands)
      allocate (param_input%exclude_bands(param_input%num_exclude_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating exclude_bands in param_read')
      call param_get_range_vector('exclude_bands', found, param_input%num_exclude_bands, .false., param_input%exclude_bands)
      if (any(param_input%exclude_bands < 1)) &
        call io_error('Error: exclude_bands must contain positive numbers')
    end if
  end subroutine param_w90_read_10

  subroutine param_read_11(pw90_effective_model, library)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: pw90_effective_model, library
    integer :: i_temp
    logical :: found

    ! AAM_2016-09-16: some changes to logic to patch a problem with uninitialised num_bands in library mode
!    num_bands       =   -1
    call param_get_keyword('num_bands', found, i_value=i_temp)
    if (found .and. library) write (stdout, '(/a)') ' Ignoring <num_bands> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      if (found) num_bands = i_temp
      if (.not. found) num_bands = num_wann
    end if
    ! GP: I subtract it here, but only the first time when I pass the total number of bands
    ! In later calls, I need to pass instead num_bands already subtracted.
    if (library .and. library_param_read_first_pass) num_bands = num_bands - param_input%num_exclude_bands
    if (.not. pw90_effective_model) then
      if (found .and. num_bands < num_wann) then
        write (stdout, *) 'num_bands', num_bands
        write (stdout, *) 'num_wann', num_wann
        call io_error('Error: num_bands must be greater than or equal to num_wann')
      endif
    endif
  end subroutine param_read_11

  subroutine param_read_13(pw90_effective_model, library)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: pw90_effective_model, library
    integer :: iv_temp(3)
    logical :: found

    param_input%devel_flag = ' '          !
    call param_get_keyword('devel_flag', found, c_value=param_input%devel_flag)

!    mp_grid=-99
    call param_get_keyword_vector('mp_grid', found, 3, i_value=iv_temp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <mp_grid> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      if (found) mp_grid = iv_temp
      if (.not. found) then
        call io_error('Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid')
      elseif (any(mp_grid < 1)) then
        call io_error('Error: mp_grid must be greater than zero')
      end if
      num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
    end if
  end subroutine param_read_13

![ysl-b]
  subroutine param_read_16(library)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: library
    logical :: found, ltmp

    ltmp = .false.  ! by default our WF are not spinors
    call param_get_keyword('spinors', found, l_value=ltmp)
    if (.not. library) then
      param_input%spinors = ltmp
    else
      if (found) write (stdout, '(a)') ' Ignoring <spinors> in input file'
    endif
!    if(spinors .and. (2*(num_wann/2))/=num_wann) &
!       call io_error('Error: For spinor WF num_wann must be even')

    ! We need to know if the bands are double degenerate due to spin, e.g. when
    ! calculating the DOS
    if (param_input%spinors) then
      param_input%num_elec_per_state = 1
    else
      param_input%num_elec_per_state = 2
    endif
    call param_get_keyword('num_elec_per_state', found, i_value=param_input%num_elec_per_state)
    if ((param_input%num_elec_per_state /= 1) .and. (param_input%num_elec_per_state /= 2)) &
      call io_error('Error: num_elec_per_state can be only 1 or 2')
    if (param_input%spinors .and. param_input%num_elec_per_state /= 1) &
      call io_error('Error: when spinors = T num_elec_per_state must be 1')
  end subroutine param_read_16

  subroutine param_w90_read_18a
    implicit none
    logical :: found

    param_input%write_xyz = .false.
    call param_get_keyword('write_xyz', found, l_value=param_input%write_xyz)
  end subroutine param_w90_read_18a

  subroutine param_read_21(bands_plot, library)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: bands_plot, library
    integer :: i_temp, ierr
    logical :: found

    spec_points%bands_num_spec_points = 0
    call param_get_block_length('kpoint_path', found, i_temp, library)
    if (found) then
      spec_points%bands_num_spec_points = i_temp*2
      if (allocated(spec_points%bands_label)) deallocate (spec_points%bands_label)
      allocate (spec_points%bands_label(spec_points%bands_num_spec_points), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating bands_label in param_read')
      if (allocated(spec_points%bands_spec_points)) deallocate (spec_points%bands_spec_points)
      allocate (spec_points%bands_spec_points(3, spec_points%bands_num_spec_points), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating bands_spec_points in param_read')
      call param_get_keyword_kpath
    end if
    if (.not. found .and. bands_plot) &
      call io_error('A bandstructure plot has been requested but there is no kpoint_path block')
  end subroutine param_read_21

  subroutine param_read_23(found_fermi_energy)
    use w90_io, only: io_error
    implicit none
    logical, intent(out) :: found_fermi_energy
    integer :: i, ierr
    logical :: found

    fermi%n = 0
    found_fermi_energy = .false.
    call param_get_keyword('fermi_energy', found, r_value=fermi_energy)
    if (found) then
      found_fermi_energy = .true.
      fermi%n = 1
    endif
    !
    fermi_energy_scan = .false.
    call param_get_keyword('fermi_energy_min', found, r_value=fermi_energy_min)
    if (found) then
      if (found_fermi_energy) call io_error( &
        'Error: Cannot specify both fermi_energy and fermi_energy_min')
      fermi_energy_scan = .true.
      fermi_energy_max = fermi_energy_min + 1.0_dp
      call param_get_keyword('fermi_energy_max', found, &
                             r_value=fermi_energy_max)
      if (found .and. fermi_energy_max <= fermi_energy_min) call io_error( &
        'Error: fermi_energy_max must be larger than fermi_energy_min')
      fermi_energy_step = 0.01_dp
      call param_get_keyword('fermi_energy_step', found, &
                             r_value=fermi_energy_step)
      if (found .and. fermi_energy_step <= 0.0_dp) call io_error( &
        'Error: fermi_energy_step must be positive')
      fermi%n = nint(abs((fermi_energy_max - fermi_energy_min)/fermi_energy_step)) + 1
    endif
    !
    if (found_fermi_energy) then
      if (allocated(fermi%energy_list)) deallocate (fermi%energy_list)
      allocate (fermi%energy_list(1), stat=ierr)
      fermi%energy_list(1) = fermi_energy
    elseif (fermi_energy_scan) then
      if (fermi%n .eq. 1) then
        fermi_energy_step = 0.0_dp
      else
        fermi_energy_step = (fermi_energy_max - fermi_energy_min)/real(fermi%n - 1, dp)
      endif
      if (allocated(fermi%energy_list)) deallocate (fermi%energy_list)
      allocate (fermi%energy_list(fermi%n), stat=ierr)
      do i = 1, fermi%n
        fermi%energy_list(i) = fermi_energy_min + (i - 1)*fermi_energy_step
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
      if (allocated(fermi%energy_list)) deallocate (fermi%energy_list)
      allocate (fermi%energy_list(1), stat=ierr)
      fermi%energy_list(1) = 0.0_dp
    endif
    if (ierr /= 0) call io_error( &
      'Error allocating fermi_energy_list in param_read')

  end subroutine param_read_23

  subroutine param_read_25(smr_index)
    use w90_io, only: io_error
    implicit none
    integer, intent(out) :: smr_index
    logical :: found
    ! [gp-begin, Apr 20, 2012]

    ! By default: Gaussian
    smr_index = 0
    call param_get_keyword('smr_type', found, c_value=ctmp)
    if (found) smr_index = get_smearing_index(ctmp, 'smr_type')

    ! By default: adaptive smearing
    adpt_smr = .true.
    call param_get_keyword('adpt_smr', found, l_value=adpt_smr)

    ! By default: a=sqrt(2)
    adpt_smr_fac = sqrt(2.0_dp)
    call param_get_keyword('adpt_smr_fac', found, r_value=adpt_smr_fac)
    if (found .and. (adpt_smr_fac <= 0._dp)) &
      call io_error('Error: adpt_smr_fac must be greater than zero')

    ! By default: 1 eV
    adpt_smr_max = 1.0_dp
    call param_get_keyword('adpt_smr_max', found, r_value=adpt_smr_max)
    if (adpt_smr_max <= 0._dp) &
      call io_error('Error: adpt_smr_max must be greater than zero')

    ! By default: if adpt_smr is manually set to false by the user, but he/she doesn't
    ! define smr_fixed_en_width: NO smearing, i.e. just the histogram
    smr_fixed_en_width = 0.0_dp
    call param_get_keyword('smr_fixed_en_width', found, r_value=smr_fixed_en_width)
    if (found .and. (smr_fixed_en_width < 0._dp)) &
      call io_error('Error: smr_fixed_en_width must be greater than or equal to zero')
    ! [gp-end]
  end subroutine param_read_25

  subroutine param_read_28
    use w90_io, only: io_error
    implicit none
    integer :: i
    logical :: found

    param_input%use_ws_distance = .true.
    call param_get_keyword('use_ws_distance', found, l_value=param_input%use_ws_distance)

    param_input%ws_distance_tol = 1.e-5_dp
    call param_get_keyword('ws_distance_tol', found, r_value=param_input%ws_distance_tol)

    param_input%ws_search_size = 2

    call param_get_vector_length('ws_search_size', found, length=i)
    if (found) then
      if (i .eq. 1) then
        call param_get_keyword_vector('ws_search_size', found, 1, &
                                      i_value=param_input%ws_search_size)
        param_input%ws_search_size(2) = param_input%ws_search_size(1)
        param_input%ws_search_size(3) = param_input%ws_search_size(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector('ws_search_size', found, 3, &
                                      i_value=param_input%ws_search_size)
      else
        call io_error('Error: ws_search_size must be provided as either one integer or a vector of three integers')
      end if
      if (any(param_input%ws_search_size <= 0)) &
        call io_error('Error: ws_search_size elements must be greater than zero')
    end if
  end subroutine param_read_28

  subroutine param_read_32(pw90_effective_model, pw90_boltzwann, &
                           pw90_geninterp, w90_plot, disentanglement, &
                           eig_found, library, postproc_setup)
    use w90_io, only: seedname, io_file_unit, io_error
    implicit none
    logical, intent(in) :: pw90_effective_model, pw90_boltzwann, &
                           pw90_geninterp, w90_plot, disentanglement, &
                           library, postproc_setup
    logical, intent(out) :: eig_found
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
            call io_error('No '//trim(seedname)//'.eig file found. Needed for disentanglement')
          else if ((w90_plot .or. pw90_boltzwann .or. pw90_geninterp)) then
            call io_error('No '//trim(seedname)//'.eig file found. Needed for interpolation')
          end if
        else
          ! Allocate only here
          allocate (eigval(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating eigval in param_read')

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
                call io_error('param_read: mismatch in '//trim(seedname)//'.eig')
              end if
            enddo
          end do
          close (eig_unit)
        end if
      end if
    end if

    if (library .and. allocated(eigval)) eig_found = .true.

105 call io_error('Error: Problem opening eigenvalue file '//trim(seedname)//'.eig')
106 call io_error('Error: Problem reading eigenvalue file '//trim(seedname)//'.eig')

  end subroutine param_read_32

  subroutine param_w90_read_33(eig_found)
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: eig_found
    integer :: nkp, ierr
    logical :: found, found2

    dis_data%win_min = -1.0_dp; dis_data%win_max = 0.0_dp
    if (eig_found) dis_data%win_min = minval(eigval)
    call param_get_keyword('dis_win_min', found, r_value=dis_data%win_min)

    if (eig_found) dis_data%win_max = maxval(eigval)
    call param_get_keyword('dis_win_max', found, r_value=dis_data%win_max)
    if (eig_found .and. (dis_data%win_max .lt. dis_data%win_min)) &
      call io_error('Error: param_read: check disentanglement windows')

    dis_data%froz_min = -1.0_dp; dis_data%froz_max = 0.0_dp
    ! no default for dis_froz_max
    dis_data%frozen_states = .false.
    call param_get_keyword('dis_froz_max', found, r_value=dis_data%froz_max)
    if (found) then
      dis_data%frozen_states = .true.
      dis_data%froz_min = dis_data%win_min ! default value for the bottom of frozen window
    end if
    call param_get_keyword('dis_froz_min', found2, r_value=dis_data%froz_min)
    if (eig_found) then
      if (dis_data%froz_max .lt. dis_data%froz_min) &
        call io_error('Error: param_read: check disentanglement frozen windows')
      if (found2 .and. .not. found) &
        call io_error('Error: param_read: found dis_froz_min but not dis_froz_max')
    endif

    dis_data%num_iter = 200
    call param_get_keyword('dis_num_iter', found, i_value=dis_data%num_iter)
    if (dis_data%num_iter < 0) call io_error('Error: dis_num_iter must be positive')

    dis_data%mix_ratio = 0.5_dp
    call param_get_keyword('dis_mix_ratio', found, r_value=dis_data%mix_ratio)
    if (dis_data%mix_ratio <= 0.0_dp .or. dis_data%mix_ratio > 1.0_dp) &
      call io_error('Error: dis_mix_ratio must be greater than 0.0 but not greater than 1.0')

    dis_data%conv_tol = 1.0e-10_dp
    call param_get_keyword('dis_conv_tol', found, r_value=dis_data%conv_tol)
    if (dis_data%conv_tol < 0.0_dp) call io_error('Error: dis_conv_tol must be positive')

    dis_data%conv_window = 3
    call param_get_keyword('dis_conv_window', found, i_value=dis_data%conv_window)
    if (dis_data%conv_window < 0) call io_error('Error: dis_conv_window must be positive')

    ! GS-start
    dis_data%spheres_first_wann = 1
    call param_get_keyword('dis_spheres_first_wann', found, i_value=dis_data%spheres_first_wann)
    if (dis_data%spheres_first_wann < 1) call io_error('Error: dis_spheres_first_wann must be greater than 0')
    if (dis_data%spheres_first_wann > num_bands - num_wann + 1) &
      call io_error('Error: dis_spheres_first_wann is larger than num_bands-num_wann+1')
    dis_data%spheres_num = 0
    call param_get_keyword('dis_spheres_num', found, i_value=dis_data%spheres_num)
    if (dis_data%spheres_num < 0) call io_error('Error: dis_spheres_num cannot be negative')
    if (dis_data%spheres_num > 0) then
      allocate (dis_data%spheres(4, dis_data%spheres_num), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dis_spheres in param_read')
      call param_get_keyword_block('dis_spheres', found, dis_data%spheres_num, 4, r_value=dis_data%spheres)
      if (.not. found) call io_error('Error: Did not find dis_spheres in the input file')
      do nkp = 1, dis_data%spheres_num
        if (dis_data%spheres(4, nkp) < 1.0e-15_dp) &
          call io_error('Error: radius for dis_spheres must be > 0')
      enddo
    endif
    ! GS-end
  end subroutine param_w90_read_33

  subroutine param_read_40a(pw90_effective_model, library)
    use w90_io, only: io_error
    use w90_utility, only: utility_recip_lattice
    implicit none
    logical, intent(in) :: pw90_effective_model, library
    real(kind=dp) :: real_lattice_tmp(3, 3), cell_volume
    integer :: itmp, nkp, ierr
    logical :: found

    kmesh_data%search_shells = 36
    call param_get_keyword('search_shells', found, i_value=kmesh_data%search_shells)
    if (kmesh_data%search_shells < 0) call io_error('Error: search_shells must be positive')

    kmesh_data%tol = 0.000001_dp
    call param_get_keyword('kmesh_tol', found, r_value=kmesh_data%tol)
    if (kmesh_data%tol < 0.0_dp) call io_error('Error: kmesh_tol must be positive')

    kmesh_data%num_shells = 0
    call param_get_range_vector('shell_list', found, kmesh_data%num_shells, lcount=.true.)
    if (found) then
      if (kmesh_data%num_shells < 0 .or. kmesh_data%num_shells > max_shells) &
        call io_error('Error: number of shell in shell_list must be between zero and six')
      if (allocated(kmesh_data%shell_list)) deallocate (kmesh_data%shell_list)
      allocate (kmesh_data%shell_list(kmesh_data%num_shells), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating shell_list in param_read')
      call param_get_range_vector('shell_list', found, kmesh_data%num_shells, .false., kmesh_data%shell_list)
      if (any(kmesh_data%shell_list < 1)) &
        call io_error('Error: shell_list must contain positive numbers')
    else
      if (allocated(kmesh_data%shell_list)) deallocate (kmesh_data%shell_list)
      allocate (kmesh_data%shell_list(max_shells), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating shell_list in param_read')
    end if

    call param_get_keyword('num_shells', found, i_value=itmp)
    if (found .and. (itmp /= kmesh_data%num_shells)) &
      call io_error('Error: Found obsolete keyword num_shells. Its value does not agree with shell_list')

    ! If .true., does not perform the check of B1 of
    ! Marzari, Vanderbild, PRB 56, 12847 (1997)
    ! in kmesh.F90
    ! mainly needed for the interaction with Z2PACK
    ! By default: .false. (perform the tests)
    kmesh_data%skip_B1_tests = .false.
    call param_get_keyword('skip_b1_tests', found, l_value=kmesh_data%skip_B1_tests)

    call param_get_keyword_block('unit_cell_cart', found, 3, 3, r_value=real_lattice_tmp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <unit_cell_cart> in input file'
    if (.not. library) then
      real_lattice = transpose(real_lattice_tmp)
      if (.not. found) call io_error('Error: Did not find the cell information in the input file')
    end if

    if (.not. library) &
      call utility_recip_lattice(real_lattice, recip_lattice, cell_volume)
    !call utility_metric(real_lattice, recip_lattice, real_metric, recip_metric)

    if (.not. pw90_effective_model) allocate (k_points%kpt_cart(3, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating kpt_cart in param_read')
    if (.not. library .and. .not. pw90_effective_model) then
      allocate (k_points%kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpt_latt in param_read')
    end if

    call param_get_keyword_block('kpoints', found, num_kpts, 3, r_value=k_points%kpt_cart)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <kpoints> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      k_points%kpt_latt = k_points%kpt_cart
      if (.not. found) call io_error('Error: Did not find the kpoint information in the input file')
    end if

    ! Calculate the kpoints in cartesian coordinates
    if (.not. pw90_effective_model) then
      do nkp = 1, num_kpts
        k_points%kpt_cart(:, nkp) = matmul(k_points%kpt_latt(:, nkp), recip_lattice(:, :))
      end do
    endif
  end subroutine param_read_40a

  subroutine param_read_40c(global_kmesh_set, kmesh_spacing, kmesh)
    use w90_io, only: io_error
    !use w90_utility, only: utility_recip_lattice
    implicit none
    logical, intent(out) :: global_kmesh_set
    real(kind=dp), intent(out) :: kmesh_spacing
    integer, intent(out) :: kmesh(3)
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
    call param_get_keyword('kmesh_spacing', found, r_value=kmesh_spacing)
    if (found) then
      if (kmesh_spacing .le. 0._dp) &
        call io_error('Error: kmesh_spacing must be greater than zero')
      global_kmesh_set = .true.

      call internal_set_kmesh(kmesh_spacing, recip_lattice, kmesh)
    end if
    call param_get_vector_length('kmesh', found, length=i)
    if (found) then
      if (global_kmesh_set) &
        call io_error('Error: cannot set both kmesh and kmesh_spacing')
      if (i .eq. 1) then
        global_kmesh_set = .true.
        call param_get_keyword_vector('kmesh', found, 1, i_value=kmesh)
        kmesh(2) = kmesh(1)
        kmesh(3) = kmesh(1)
      elseif (i .eq. 3) then
        global_kmesh_set = .true.
        call param_get_keyword_vector('kmesh', found, 3, i_value=kmesh)
      else
        call io_error('Error: kmesh must be provided as either one integer or a vector of three integers')
      end if
      if (any(kmesh <= 0)) &
        call io_error('Error: kmesh elements must be greater than zero')
    end if
    ! [GP-end]
  end subroutine param_read_40c

  subroutine param_w90_read_41(library)
    use w90_io, only: io_error, stdout
    implicit none
    logical, intent(in) :: library
    integer :: i_temp, i_temp2
    logical :: found, found2, lunits

    ! Atoms
    if (.not. library) atoms%num_atoms = 0
    call param_get_block_length('atoms_frac', found, i_temp, library)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <atoms_frac> in input file'
    call param_get_block_length('atoms_cart', found2, i_temp2, library, lunits)
    if (found2 .and. library) write (stdout, '(a)') ' Ignoring <atoms_cart> in input file'
    if (.not. library) then
      if (found .and. found2) call io_error('Error: Cannot specify both atoms_frac and atoms_cart')
      if (found .and. i_temp > 0) then
        lunits = .false.
        atoms%num_atoms = i_temp
      elseif (found2 .and. i_temp2 > 0) then
        atoms%num_atoms = i_temp2
        if (lunits) atoms%num_atoms = atoms%num_atoms - 1
      end if
      if (atoms%num_atoms > 0) then
        call param_get_atoms(lunits, library)
      end if
    endif
  end subroutine param_w90_read_41

  subroutine param_read_44(read_transport)
    use w90_io, only: seedname, stdout, io_error
    implicit none
    logical, intent(in) :: read_transport
    integer :: loop, ierr

    if (any(len_trim(in_data(:)) > 0)) then
      write (stdout, '(1x,a)') 'The following section of file '//trim(seedname)//'.win contained unrecognised keywords'
      write (stdout, *)
      do loop = 1, num_lines
        if (len_trim(in_data(loop)) > 0) then
          write (stdout, '(1x,a)') trim(in_data(loop))
        end if
      end do
      write (stdout, *)
      call io_error('Unrecognised keyword(s) in input file, see also output file')
    end if

    if (.not. read_transport) then

      ! For aesthetic purposes, convert some things to uppercase
      call param_uppercase()

    endif

    deallocate (in_data, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating in_data in param_read')

  end subroutine param_read_44

  subroutine param_w90_read_45(disentanglement)
    ! =============================== !
    ! Some checks and initialisations !
    ! =============================== !
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: disentanglement
    integer :: ierr

!    if (restart.ne.' ') disentanglement=.false.

    if (disentanglement) then
      if (allocated(dis_data%ndimwin)) deallocate (dis_data%ndimwin)
      allocate (dis_data%ndimwin(num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating ndimwin in param_read')
      if (allocated(dis_data%lwindow)) deallocate (dis_data%lwindow)
      allocate (dis_data%lwindow(num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating lwindow in param_read')
    endif

!    if ( wannier_plot .and. (index(wannier_plot_format,'cub').ne.0) ) then
!       cosa(1)=dot_product(real_lattice(1,:),real_lattice(2,:))
!       cosa(2)=dot_product(real_lattice(1,:),real_lattice(3,:))
!       cosa(3)=dot_product(real_lattice(2,:),real_lattice(3,:))
!       cosa = abs(cosa)
!       if (any(cosa.gt.eps6)) &
!            call io_error('Error: plotting in cube format requires orthogonal lattice vectors')
!    endif

    ! Initialise
    param_input%omega_invariant = -999.0_dp
    param_input%have_disentangled = .false.

    if (allocated(wann_data%centres)) deallocate (wann_data%centres)
    allocate (wann_data%centres(3, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wannier_centres in param_read')
    wann_data%centres = 0.0_dp
    if (allocated(wann_data%spreads)) deallocate (wann_data%spreads)
    allocate (wann_data%spreads(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating wannier_spreads in param_read')
    wann_data%spreads = 0.0_dp
  end subroutine param_w90_read_45

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

  function get_smearing_index(string, keyword)
    !! This function parses a string containing the type of
    !! smearing and returns the correct index for the smearing_index variable
    !
    !! If the string is not valid, an io_error is issued
    use w90_io, only: io_error
    character(len=*), intent(in) :: string
    !! The string read from input
    character(len=*), intent(in) :: keyword
    !! The keyword that was read (e.g., smr_type), so that we can print a more useful error message
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
          call io_error('Wrong m-p smearing order in keyword '//trim(keyword))
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
      call io_error('Unrecognised value for keyword '//trim(keyword))
    end if

    return

337 call io_error('Wrong m-p smearing order in keyword '//trim(keyword))

  end function get_smearing_index

!===================================================================
  subroutine param_uppercase
    !===================================================================
    !                                                                  !
    !! Convert a few things to uppercase to look nice in the output
    !                                                                  !
    !===================================================================

    implicit none

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
    do loop = 1, spec_points%bands_num_spec_points
      do inner_loop = 1, len(spec_points%bands_label(loop))
        ic = ichar(spec_points%bands_label(loop) (inner_loop:inner_loop))
        if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
          spec_points%bands_label(loop) (inner_loop:inner_loop) = char(ic + ichar('Z') - ichar('z'))
      enddo
    enddo

    ! Length unit (ang --> Ang, bohr --> Bohr)
    ic = ichar(param_input%length_unit(1:1))
    if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
      param_input%length_unit(1:1) = char(ic + ichar('Z') - ichar('z'))

    return

  end subroutine param_uppercase

  subroutine param_write_header
    !! Write a suitable header for the calculation - version authors etc
    use w90_io, only: io_date, w90_version
    use w90_constants, only: bohr_version_str, constants_version_str1, constants_version_str2
    implicit none

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
  subroutine param_dealloc
    !==================================================================!
    !                                                                  !
    !! release memory from allocated parameters
    !                                                                  !
    !===================================================================
    use w90_io, only: io_error

    implicit none
    integer :: ierr

    if (allocated(dis_data%ndimwin)) then
      deallocate (dis_data%ndimwin, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ndimwin in param_dealloc')
    end if
    if (allocated(dis_data%lwindow)) then
      deallocate (dis_data%lwindow, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating lwindow in param_dealloc')
    end if
    if (allocated(eigval)) then
      deallocate (eigval, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating eigval in param_dealloc')
    endif
    if (allocated(kmesh_data%shell_list)) then
      deallocate (kmesh_data%shell_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating shell_list in param_dealloc')
    endif
    if (allocated(k_points%kpt_latt)) then
      deallocate (k_points%kpt_latt, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kpt_latt in param_dealloc')
    endif
    if (allocated(k_points%kpt_cart)) then
      deallocate (k_points%kpt_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kpt_cart in param_dealloc')
    endif
    if (allocated(spec_points%bands_label)) then
      deallocate (spec_points%bands_label, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating bands_label in param_dealloc')
    end if
    if (allocated(spec_points%bands_spec_points)) then
      deallocate (spec_points%bands_spec_points, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating bands_spec_points in param_dealloc')
    end if
    if (allocated(atoms%label)) then
      deallocate (atoms%label, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_label in param_dealloc')
    end if
    if (allocated(atoms%symbol)) then
      deallocate (atoms%symbol, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_symbol in param_dealloc')
    end if
    if (allocated(atoms%pos_frac)) then
      deallocate (atoms%pos_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atom_pos_frac in param_dealloc')
    end if
    if (allocated(atoms%pos_cart)) then
      deallocate (atoms%pos_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_pos_cart in param_dealloc')
    end if
    if (allocated(atoms%species_num)) then
      deallocate (atoms%species_num, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_species_num in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj_site)) then
      deallocate (kmesh_data%input_proj_site, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_site in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%l)) then
      deallocate (kmesh_data%input_proj%l, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_l in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%m)) then
      deallocate (kmesh_data%input_proj%m, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_m in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%s)) then
      deallocate (kmesh_data%input_proj%s, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_s in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%s_qaxis)) then
      deallocate (kmesh_data%input_proj%s_qaxis, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_s_qaxis in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%z)) then
      deallocate (kmesh_data%input_proj%z, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_z in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%x)) then
      deallocate (kmesh_data%input_proj%x, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_x in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%radial)) then
      deallocate (kmesh_data%input_proj%radial, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_radial in param_dealloc')
    end if
    if (allocated(kmesh_data%input_proj%zona)) then
      deallocate (kmesh_data%input_proj%zona, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_zona in param_dealloc')
    end if
    if (allocated(param_input%exclude_bands)) then
      deallocate (param_input%exclude_bands, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating exclude_bands in param_dealloc')
    end if
    if (allocated(wann_data%centres)) then
      deallocate (wann_data%centres, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_centres in param_dealloc')
    end if
    if (allocated(wann_data%spreads)) then
      deallocate (wann_data%spreads, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_spreads in param_dealloc')
    endif
    if (allocated(dis_data%spheres)) then
      deallocate (dis_data%spheres, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating dis_spheres in param_dealloc')
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
  subroutine param_read_chkpt(ispostw90, checkpoint, m_matrix)
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
    use w90_io, only: io_error, io_file_unit, stdout, seedname

    implicit none
    logical, intent(in) :: ispostw90
    character(len=*), intent(out) :: checkpoint
    complex(kind=dp), allocatable, intent(inout) :: m_matrix(:, :, :, :)
    integer :: chk_unit, nkp, i, j, k, l, ntmp, ierr
    character(len=33) :: header
    real(kind=dp) :: tmp_latt(3, 3), tmp_kpt_latt(3, num_kpts)
    integer :: tmp_excl_bands(1:param_input%num_exclude_bands), tmp_mp_grid(1:3)
    complex(kind=dp) :: m_tmp ! postw90 dummy for m_matrix read

    write (stdout, '(1x,3a)') 'Reading restart information from file ', trim(seedname), '.chk :'

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk', status='old', form='unformatted', err=121)

    ! Read comment line
    read (chk_unit) header
    write (stdout, '(1x,a)', advance='no') trim(header)

    ! Consistency checks
    read (chk_unit) ntmp                           ! Number of bands
    if (ntmp .ne. num_bands) call io_error('param_read_chk: Mismatch in num_bands')
    read (chk_unit) ntmp                           ! Number of excluded bands
    if (ntmp .ne. param_input%num_exclude_bands) &
      call io_error('param_read_chk: Mismatch in num_exclude_bands')
    read (chk_unit) (tmp_excl_bands(i), i=1, param_input%num_exclude_bands) ! Excluded bands
    do i = 1, param_input%num_exclude_bands
      if (tmp_excl_bands(i) .ne. param_input%exclude_bands(i)) &
        call io_error('param_read_chk: Mismatch in exclude_bands')
    enddo
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Real lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - real_lattice(i, j)) .gt. eps6) &
          call io_error('param_read_chk: Mismatch in real_lattice')
      enddo
    enddo
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Reciprocal lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - recip_lattice(i, j)) .gt. eps6) &
          call io_error('param_read_chk: Mismatch in recip_lattice')
      enddo
    enddo
    read (chk_unit) ntmp                ! K-points
    if (ntmp .ne. num_kpts) &
      call io_error('param_read_chk: Mismatch in num_kpts')
    read (chk_unit) (tmp_mp_grid(i), i=1, 3)         ! M-P grid
    do i = 1, 3
      if (tmp_mp_grid(i) .ne. mp_grid(i)) &
        call io_error('param_read_chk: Mismatch in mp_grid')
    enddo
    read (chk_unit) ((tmp_kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts)
    do nkp = 1, num_kpts
      do i = 1, 3
        if (abs(tmp_kpt_latt(i, nkp) - k_points%kpt_latt(i, nkp)) .gt. eps6) &
          call io_error('param_read_chk: Mismatch in kpt_latt')
      enddo
    enddo
    read (chk_unit) ntmp                ! nntot
    if (ntmp .ne. kmesh_info%nntot) &
      call io_error('param_read_chk: Mismatch in nntot')
    read (chk_unit) ntmp                ! num_wann
    if (ntmp .ne. num_wann) &
      call io_error('param_read_chk: Mismatch in num_wann')
    ! End of consistency checks

    read (chk_unit) checkpoint             ! checkpoint
    checkpoint = adjustl(trim(checkpoint))

    read (chk_unit) param_input%have_disentangled      ! whether a disentanglement has been performed

    if (param_input%have_disentangled) then

      read (chk_unit) param_input%omega_invariant     ! omega invariant

      ! lwindow
      if (.not. allocated(dis_data%lwindow)) then
        allocate (dis_data%lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating lwindow in param_read_chkpt')
      endif
      read (chk_unit, err=122) ((dis_data%lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)

      ! ndimwin
      if (.not. allocated(dis_data%ndimwin)) then
        allocate (dis_data%ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ndimwin in param_read_chkpt')
      endif
      read (chk_unit, err=123) (dis_data%ndimwin(nkp), nkp=1, num_kpts)

      ! U_matrix_opt
      if (.not. allocated(u_matrix_opt)) then
        allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating u_matrix_opt in param_read_chkpt')
      endif
      read (chk_unit, err=124) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)

    endif

    ! U_matrix
    if (.not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating u_matrix in param_read_chkpt')
    endif
    read (chk_unit, err=125) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)

    ! M_matrix
    if (ispostw90) then
      read (chk_unit, err=126) ((((m_tmp, i=1, num_wann), j=1, num_wann), k=1, kmesh_info%nntot), l=1, num_kpts)
    else
      if (.not. allocated(m_matrix)) then
        allocate (m_matrix(num_wann, num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating m_matrix in param_read_chkpt')
      endif
      read (chk_unit, err=126) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, kmesh_info%nntot), l=1, num_kpts)
    endif

    ! wannier_centres
    read (chk_unit, err=127) ((wann_data%centres(i, j), i=1, 3), j=1, num_wann)

    ! wannier spreads
    read (chk_unit, err=128) (wann_data%spreads(i), i=1, num_wann)

    close (chk_unit)

    write (stdout, '(a/)') ' ... done'

    return

121 if (ispostw90) then
      call io_error('Error opening '//trim(seedname)//'.chk in param_read_chkpt: did you run wannier90.x first?')
    else
      call io_error('Error opening '//trim(seedname)//'.chk in param_read_chkpt')
    end if
122 call io_error('Error reading lwindow from '//trim(seedname)//'.chk in param_read_chkpt')
123 call io_error('Error reading ndimwin from '//trim(seedname)//'.chk in param_read_chkpt')
124 call io_error('Error reading u_matrix_opt from '//trim(seedname)//'.chk in param_read_chkpt')
125 call io_error('Error reading u_matrix from '//trim(seedname)//'.chk in param_read_chkpt')
126 call io_error('Error reading m_matrix from '//trim(seedname)//'.chk in param_read_chkpt')
127 call io_error('Error reading wannier_centres from '//trim(seedname)//'.chk in param_read_chkpt')
128 call io_error('Error reading wannier_spreads from '//trim(seedname)//'.chk in param_read_chkpt')

  end subroutine param_read_chkpt

!===========================================================!
  subroutine param_chkpt_dist(checkpoint)
    !===========================================================!
    !                                                           !
    !! Distribute the chk files
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, &
      io_date, io_time, io_stopwatch
    use w90_comms, only: on_root, comms_bcast

    implicit none
    character(len=*), intent(inout) :: checkpoint
    integer :: ierr

    call comms_bcast(checkpoint, len(checkpoint))

    if (.not. on_root .and. .not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating u_matrix in param_chkpt_dist')
    endif
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts)

!    if (.not.on_root .and. .not.allocated(m_matrix)) then
!       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!       if (ierr/=0)&
!            call io_error('Error allocating m_matrix in param_chkpt_dist')
!    endif
!    call comms_bcast(m_matrix(1,1,1,1),num_wann*num_wann*nntot*num_kpts)

    call comms_bcast(param_input%have_disentangled, 1)

    if (param_input%have_disentangled) then
      if (.not. on_root) then

        if (.not. allocated(u_matrix_opt)) then
          allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating u_matrix_opt in param_chkpt_dist')
        endif

        if (.not. allocated(dis_data%lwindow)) then
          allocate (dis_data%lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating lwindow in param_chkpt_dist')
        endif

        if (.not. allocated(dis_data%ndimwin)) then
          allocate (dis_data%ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating ndimwin in param_chkpt_dist')
        endif

      end if

      call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts)
      call comms_bcast(dis_data%lwindow(1, 1), num_bands*num_kpts)
      call comms_bcast(dis_data%ndimwin(1), num_kpts)
      call comms_bcast(param_input%omega_invariant, 1)
    end if
    call comms_bcast(wann_data%centres(1, 1), 3*num_wann)
    call comms_bcast(wann_data%spreads(1), num_wann)

  end subroutine param_chkpt_dist

!=======================================!
  subroutine param_in_file
    !=======================================!
    !! Load the *.win file into a character
    !! array in_file, ignoring comments and
    !! blank lines and converting everything
    !! to lowercase characters
    !=======================================!

    use w90_io, only: io_file_unit, io_error, seedname
    use w90_utility, only: utility_lowercase

    implicit none

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

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.win')
200 call io_error('Error: Problem reading input file '//trim(seedname)//'.win')
210 continue
    rewind (in_unit)

    allocate (in_data(num_lines), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating in_data in param_in_file')

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
  subroutine param_get_keyword(keyword, found, c_value, l_value, i_value, r_value)
    !===========================================================================!
    !                                                                           !
    !! Finds the value of the required keyword.
    !                                                                           !
    !===========================================================================!

    use w90_io, only: io_error

    implicit none

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
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
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
          call io_error('Error: Problem reading logical keyword '//trim(keyword))
        endif
      endif
      if (present(i_value)) read (dummy, *, err=220, end=220) i_value
      if (present(r_value)) read (dummy, *, err=220, end=220) r_value
    end if

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword))

  end subroutine param_get_keyword

!=========================================================================================!
  subroutine param_get_keyword_vector(keyword, found, length, c_value, l_value, i_value, r_value)
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

    integer           :: kl, in, loop, i
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
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
        call io_error('param_get_keyword_vector unimplemented for logicals')
      endif
      if (present(i_value)) read (dummy, *, err=230, end=230) (i_value(i), i=1, length)
      if (present(r_value)) read (dummy, *, err=230, end=230) (r_value(i), i=1, length)
    end if

    return

230 call io_error('Error: Problem reading keyword '//trim(keyword)//' in param_get_keyword_vector')

  end subroutine param_get_keyword_vector

!========================================================!
  subroutine param_get_vector_length(keyword, found, length)
    !======================================================!
    !                                                      !
    !! Returns the length of a keyword vector
    !                                                      !
    !======================================================!

    use w90_io, only: io_error

    implicit none

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
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
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
      if (len_trim(dummy) == 0) call io_error('Error: keyword '//trim(keyword)//' is blank')
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
  subroutine param_get_keyword_block(keyword, found, rows, columns, c_value, l_value, i_value, r_value)
    !==============================================================================================!
    !                                                                                              !
    !!   Finds the values of the required data block
    !                                                                                              !
    !==============================================================================================!

    use w90_constants, only: bohr
    use w90_io, only: io_error

    implicit none

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
        call io_error('Error: Found '//trim(start_st)//' more than once in input file')
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
        call io_error('Error: Found '//trim(end_st)//' more than once in input file')
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
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

    if ((blen .ne. rows) .and. (blen .ne. rows + 1)) &
      call io_error('Error: Wrong number of lines in block '//trim(keyword))

    if ((blen .eq. rows + 1) .and. (index(trim(keyword), 'unit_cell_cart') .eq. 0)) &
      call io_error('Error: Wrong number of lines in block '//trim(keyword))

    found = .true.

    lconvert = .false.
    if (blen == rows + 1) then
      dummy = in_data(line_s + 1)
      if (index(dummy, 'ang') .ne. 0) then
        lconvert = .false.
      elseif (index(dummy, 'bohr') .ne. 0) then
        lconvert = .true.
      else
        call io_error('Error: Units in block '//trim(keyword)//' not recognised')
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
        call io_error('param_get_keyword_block unimplemented for logicals')
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

240 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine param_get_keyword_block

!=====================================================!
  subroutine param_get_block_length(keyword, found, rows, library, lunits)
    !=====================================================!
    !                                                     !
    !! Finds the length of the data block
    !                                                     !
    !=====================================================!

    use w90_io, only: io_error

    implicit none

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
        call io_error('Error: Found '//trim(start_st)//' more than once in input file')
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
        call io_error('Error: Found '//trim(end_st)//' more than once in input file')
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
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
  subroutine param_get_atoms(lunits, library)
    !===================================!
    !                                   !
    !!   Fills the atom data block
    !                                   !
    !===================================!

    use w90_constants, only: bohr
    use w90_utility, only: utility_frac_to_cart, utility_cart_to_frac
    use w90_io, only: io_error
    implicit none

    logical, intent(in) :: lunits, library
    !! Do we expect a first line with the units

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
    call param_get_block_length("atoms_frac", found, i_temp, library)
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
        call io_error('Error: Found '//trim(start_st)//' more than once in input file')
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
        call io_error('Error: Found '//trim(end_st)//' more than once in input file')
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    lconvert = .false.
    if (lunits) then
      dummy = in_data(line_s + 1)
      if (index(dummy, 'ang') .ne. 0) then
        lconvert = .false.
      elseif (index(dummy, 'bohr') .ne. 0) then
        lconvert = .true.
      else
        call io_error('Error: Units in block atoms_cart not recognised in param_get_atoms')
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
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in param_get_atoms')
    allocate (atoms%label(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in param_get_atoms')
    allocate (atoms%symbol(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in param_get_atoms')
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
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in param_get_atoms')
    allocate (atoms%pos_cart(3, max_sites, atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in param_get_atoms')

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

240 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine param_get_atoms

!=====================================================!
  subroutine param_lib_set_atoms(atoms_label_tmp, atoms_pos_cart_tmp)
    !=====================================================!
    !                                                     !
    !!   Fills the atom data block during a library call
    !                                                     !
    !=====================================================!

    use w90_utility, only: utility_cart_to_frac, utility_lowercase
    use w90_io, only: io_error

    implicit none

    character(len=*), intent(in) :: atoms_label_tmp(atoms%num_atoms)
    !! Atom labels
    real(kind=dp), intent(in)      :: atoms_pos_cart_tmp(3, atoms%num_atoms)
    !! Atom positions

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
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in param_lib_set_atoms')
    allocate (atoms%label(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in param_lib_set_atoms')
    allocate (atoms%symbol(atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in param_lib_set_atoms')
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
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in param_lib_set_atoms')
    allocate (atoms%pos_cart(3, max_sites, atoms%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in param_lib_set_atoms')

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
  subroutine param_get_range_vector(keyword, found, length, lcount, i_value)
    !====================================================================!
    !!   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100
    !!   if(lcount) we return the number of states in length
    !====================================================================!
    use w90_io, only: io_error

    implicit none

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

    if (lcount .and. present(i_value)) call io_error('param_get_range_vector: incorrect call')

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
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
    if (len_trim(dummy) == 0) call io_error('Error: keyword '//trim(keyword)//' is blank')
    dummy = adjustl(dummy)
    do
      i_punc = scan(dummy, c_punc)
      if (i_punc == 0) call io_error('Error parsing keyword '//trim(keyword))
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
      if (scan(dummy, c_range) == 1) call io_error('Error parsing keyword '//trim(keyword)//' incorrect range')
      if (index(dummy, ' ') == 1) exit
    end do

    if (lcount) length = counter
    if (.not. lcount) then
      do loop = 1, counter - 1
        do loop_r = loop + 1, counter
          if (i_value(loop) == i_value(loop_r)) &
            call io_error('Error parsing keyword '//trim(keyword)//' duplicate values')
        end do
      end do
    end if

    return

101 call io_error('Error parsing keyword '//trim(keyword))

  end subroutine param_get_range_vector

  subroutine param_get_centre_constraints(ccentres_frac, ccentres_cart, &
                                          proj_site)
    !=============================================================================!
    !                                                                             !
    !!  assigns projection centres as default centre constraints and global
    !!  Lagrange multiplier as individual Lagrange multipliers then reads
    !!  the centre_constraints block for individual centre constraint parameters
    !                                                                             !
    !=============================================================================!
    use w90_io, only: io_error
    use w90_utility, only: utility_frac_to_cart
    real(kind=dp), intent(inout) :: ccentres_frac(:, :), ccentres_cart(:, :)
    real(kind=dp), intent(in) :: proj_site(:, :)
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
        if (index1 > 0) call io_error("slwf_centres block hasn't ended yet")
        index1 = index(dummy, 'end')
        if (index1 > 0) then
          index1 = index(dummy, 'slwf_centres')
          if (index1 == 0) call io_error('Wrong ending of block (need to end slwf_centres)')
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
              call param_get_centre_constraint_from_column(column, start, finish, wann, dummy, ccentres_frac)
              start = loop2 + 1
              finish = start
            end if
          end if
          if (loop2 == len_trim(dummy) .and. dummy(loop2:loop2) /= ' ') then
            finish = loop2
            call param_get_centre_constraint_from_column(column, start, finish, wann, dummy, ccentres_frac)
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

  subroutine param_get_centre_constraint_from_column(column, start, finish, wann, dummy, ccentres_frac)
    !===================================!
    !                                   !
    !!  assigns value read to constraint
    !!  parameters based on column
    !                                   !
    !===================================!
    use w90_io, only: io_error
    integer, intent(inout):: column, start, finish, wann
    character(len=maxlen), intent(inout):: dummy
    real(kind=dp), intent(inout) :: ccentres_frac(:, :)
    if (column == 0) then
      read (dummy(start:finish), '(i3)') wann
    end if
    if (column > 0) then
      if (column > 4) call io_error("Didn't expect anything else after Lagrange multiplier")
      if (column < 4) read (dummy(start:finish), '(f10.10)') ccentres_frac(wann, column)
    end if
    column = column + 1
  end subroutine param_get_centre_constraint_from_column

!===================================!
  subroutine param_get_projections(num_proj, lcount, proj_site, proj)
    !===================================!
    !                                   !
    !!  Fills the projection data block
    !                                   !
    !===================================!

    use w90_constants, only: bohr, eps6, eps2
    use w90_utility, only: utility_cart_to_frac, &
      utility_string_to_coord, utility_strip
    use w90_io, only: io_error

    implicit none

    integer, intent(inout) :: num_proj
    logical, intent(in)    :: lcount
    real(kind=dp), allocatable, dimension(:, :), intent(out) :: proj_site
    type(projection_type), intent(inout) :: proj ! intent(out)?

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
      allocate (kmesh_data%input_proj_site(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_site in param_get_projections')
      allocate (kmesh_data%input_proj%l(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_l in param_get_projections')
      allocate (kmesh_data%input_proj%m(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_m in param_get_projections')
      allocate (kmesh_data%input_proj%z(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_z in param_get_projections')
      allocate (kmesh_data%input_proj%x(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_x in param_get_projections')
      allocate (kmesh_data%input_proj%radial(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_radial in param_get_projections')
      allocate (kmesh_data%input_proj%zona(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_zona in param_get_projections')
      if (param_input%spinors) then
        allocate (kmesh_data%input_proj%s(num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s in param_get_projections')
        allocate (kmesh_data%input_proj%s_qaxis(3, num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s_qaxis in param_get_projections')
      endif

      allocate (proj_site(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_site in param_get_projections')
      allocate (proj%l(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_l in param_get_projections')
      allocate (proj%m(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_m in param_get_projections')
      allocate (proj%z(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_z in param_get_projections')
      allocate (proj%x(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_x in param_get_projections')
      allocate (proj%radial(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_radial in param_get_projections')
      allocate (proj%zona(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_zona in param_get_projections')
      if (param_input%spinors) then
        allocate (proj%s(num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_s in param_get_projections')
        allocate (proj%s_qaxis(3, num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_s_qaxis in param_get_projections')
      endif
    endif

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file')
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
        call io_error('param_get_projections: Found '//trim(end_st)//' more than once in input file')
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('param_get_projections: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e <= line_s) then
      call io_error('param_get_projections: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
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
        if (param_input%spinors) then
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
          call io_error('param_read_projection: malformed projection definition: '//trim(dummy))
        sites = 0
        ctemp = dummy(:pos1 - 1)
        ! Read the atomic site
        if (index(ctemp, 'c=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_cart)
          if (lconvert) pos_cart = pos_cart*bohr
          call utility_cart_to_frac(pos_cart(:), pos_frac(:), recip_lattice)
        elseif (index(ctemp, 'f=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_frac)
        else
          if (atoms%num_species == 0) &
            call io_error('param_get_projection: Atom centred projection requested but no atoms defined')
          do loop = 1, atoms%num_species
            if (trim(ctemp) == atoms%label(loop)) then
              species = loop
              sites = atoms%species_num(loop)
              exit
            end if
            if (loop == atoms%num_species) call io_error('param_get_projection: Atom site not recognised '//trim(ctemp))
          end do
        end if

        dummy = dummy(pos1 + 1:)

        ! scan for quantisation direction
        pos1 = index(dummy, '[')
        if (param_input%spinors) then
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 1:))
            pos2 = index(ctemp, ']')
            if (pos2 == 0) call io_error &
              ('param_get_projections: no closing square bracket for spin quantisation dir')
            ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_s_qaxis_tmp)
            dummy = dummy(:pos1 - 1) ! remove [ ] section
          endif
        else
          if (pos1 > 0) call io_error('param_get_projections: spin qdir is defined but spinors=.false.')
        endif

        ! scan for up or down
        pos1 = index(dummy, '(')
        if (param_input%spinors) then
          if (pos1 > 0) then
            proj_u_tmp = .false.; proj_d_tmp = .false.
            ctemp = (dummy(pos1 + 1:))
            pos2 = index(ctemp, ')')
            if (pos2 == 0) call io_error('param_get_projections: no closing bracket for spin')
            ctemp = ctemp(:pos2 - 1)
            if (index(ctemp, 'u') > 0) proj_u_tmp = .true.
            if (index(ctemp, 'd') > 0) proj_d_tmp = .true.
            if (proj_u_tmp .and. proj_d_tmp) then
              spn_counter = 2
            elseif (.not. proj_u_tmp .and. .not. proj_d_tmp) then
              call io_error('param_get_projections: found brackets but neither u or d')
            else
              spn_counter = 1
            endif
            dummy = dummy(:pos1 - 1) ! remove ( ) section
          endif
        else
          if (pos1 > 0) call io_error('param_get_projections: spin is defined but spinors=.false.')
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
            if (l_tmp < -5 .or. l_tmp > 3) call io_error('param_get_projection: Incorrect l state requested')
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
                call io_error('param_get_projection: Problem reading m state')
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
                  if ((m_tmp > 2*l_tmp + 1) .or. (m_tmp <= 0)) call io_error('param_get_projection: m is > l !')
                elseif (l_tmp == -1 .and. (m_tmp > 2 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (1)')
                elseif (l_tmp == -2 .and. (m_tmp > 3 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (2)')
                elseif (l_tmp == -3 .and. (m_tmp > 4 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (3)')
                elseif (l_tmp == -4 .and. (m_tmp > 5 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (4)')
                elseif (l_tmp == -5 .and. (m_tmp > 6 .or. m_tmp <= 0)) then
                  call io_error('param_get_projection: m has incorrect value (5)')
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
                call io_error('param_get_projection: Problem reading l state '//trim(ctemp3))
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
            call utility_string_to_coord(ctemp, proj_z_tmp)
          endif
          ! x axis
          pos1 = index(dummy, 'x=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_x_tmp)
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
                  kmesh_data%input_proj_site(:, counter) = pos_frac
                  kmesh_data%input_proj%l(counter) = loop_l
                  kmesh_data%input_proj%m(counter) = loop_m
                  kmesh_data%input_proj%z(:, counter) = proj_z_tmp
                  kmesh_data%input_proj%x(:, counter) = proj_x_tmp
                  kmesh_data%input_proj%radial(counter) = proj_radial_tmp
                  kmesh_data%input_proj%zona(counter) = proj_zona_tmp
                  if (param_input%spinors) then
                    if (spn_counter == 1) then
                      if (proj_u_tmp) kmesh_data%input_proj%s(counter) = 1
                      if (proj_d_tmp) kmesh_data%input_proj%s(counter) = -1
                    else
                      if (loop_s == 1) kmesh_data%input_proj%s(counter) = 1
                      if (loop_s == 2) kmesh_data%input_proj%s(counter) = -1
                    endif
                    kmesh_data%input_proj%s_qaxis(:, counter) = proj_s_qaxis_tmp
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
                    kmesh_data%input_proj_site(:, counter) = atoms%pos_frac(:, loop_sites, species)
                    kmesh_data%input_proj%l(counter) = loop_l
                    kmesh_data%input_proj%m(counter) = loop_m
                    kmesh_data%input_proj%z(:, counter) = proj_z_tmp
                    kmesh_data%input_proj%x(:, counter) = proj_x_tmp
                    kmesh_data%input_proj%radial(counter) = proj_radial_tmp
                    kmesh_data%input_proj%zona(counter) = proj_zona_tmp
                    if (param_input%spinors) then
                      if (spn_counter == 1) then
                        if (proj_u_tmp) kmesh_data%input_proj%s(counter) = 1
                        if (proj_d_tmp) kmesh_data%input_proj%s(counter) = -1
                      else
                        if (loop_s == 1) kmesh_data%input_proj%s(counter) = 1
                        if (loop_s == 2) kmesh_data%input_proj%s(counter) = -1
                      endif
                      kmesh_data%input_proj%s_qaxis(:, counter) = proj_s_qaxis_tmp
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
          'param_get_projections: too few projection functions defined')
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
        call random_number(kmesh_data%input_proj_site(:, loop))
        kmesh_data%input_proj%l(loop) = 0
        kmesh_data%input_proj%m(loop) = 1
        kmesh_data%input_proj%z(:, loop) = proj_z_def
        kmesh_data%input_proj%x(:, loop) = proj_x_def
        kmesh_data%input_proj%zona(loop) = proj_zona_def
        kmesh_data%input_proj%radial(loop) = proj_radial_def
        if (param_input%spinors) then
          if (modulo(loop, 2) == 1) then
            kmesh_data%input_proj%s(loop) = 1
          else
            kmesh_data%input_proj%s(loop) = -1
          end if
          kmesh_data%input_proj%s_qaxis(1, loop) = 0.
          kmesh_data%input_proj%s_qaxis(2, loop) = 0.
          kmesh_data%input_proj%s_qaxis(3, loop) = 1.
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

      znorm = sqrt(sum(kmesh_data%input_proj%z(:, loop)*kmesh_data%input_proj%z(:, loop)))
      xnorm = sqrt(sum(kmesh_data%input_proj%x(:, loop)*kmesh_data%input_proj%x(:, loop)))
      kmesh_data%input_proj%z(:, loop) = kmesh_data%input_proj%z(:, loop)/znorm             ! normalise z
      kmesh_data%input_proj%x(:, loop) = kmesh_data%input_proj%x(:, loop)/xnorm             ! normalise x
      cosphi = sum(kmesh_data%input_proj%z(:, loop)*kmesh_data%input_proj%x(:, loop))

      ! Check whether z-axis and z-axis are orthogonal
      if (abs(cosphi) .gt. eps6) then

        ! Special case of circularly symmetric projections (pz, dz2, fz3)
        ! just choose an x-axis that is perpendicular to the given z-axis
        if ((kmesh_data%input_proj%l(loop) .ge. 0) .and. (kmesh_data%input_proj%m(loop) .eq. 1)) then
          proj_x_tmp(:) = kmesh_data%input_proj%x(:, loop)            ! copy of original x-axis
          call random_seed()
          call random_number(proj_z_tmp(:))         ! random vector
          ! calculate new x-axis as the cross (vector) product of random vector with z-axis
          kmesh_data%input_proj%x(1, loop) = proj_z_tmp(2)*kmesh_data%input_proj%z(3, loop) &
                                             - proj_z_tmp(3)*kmesh_data%input_proj%z(2, loop)
          kmesh_data%input_proj%x(2, loop) = proj_z_tmp(3)*kmesh_data%input_proj%z(1, loop) &
                                             - proj_z_tmp(1)*kmesh_data%input_proj%z(3, loop)
          kmesh_data%input_proj%x(3, loop) = proj_z_tmp(1)*kmesh_data%input_proj%z(2, loop) &
                                             - proj_z_tmp(2)*kmesh_data%input_proj%z(1, loop)
          xnorm_new = sqrt(sum(kmesh_data%input_proj%x(:, loop)*kmesh_data%input_proj%x(:, loop)))
          kmesh_data%input_proj%x(:, loop) = kmesh_data%input_proj%x(:, loop)/xnorm_new   ! normalise
          goto 555
        endif

        ! If projection axes non-orthogonal enough, then
        ! user may have made a mistake and should check
        if (abs(cosphi) .gt. eps2) then
          write (stdout, *) ' Projection:', loop
          call io_error(' Error in projections: z and x axes are not orthogonal')
        endif

        ! If projection axes are "reasonably orthogonal", project x-axis
        ! onto plane perpendicular to z-axis to make them more so
        sinphi = sqrt(1 - cosphi*cosphi)
        proj_x_tmp(:) = kmesh_data%input_proj%x(:, loop)               ! copy of original x-axis
        ! calculate new x-axis:
        ! x = z \cross (x_tmp \cross z) / sinphi = ( x_tmp - z(z.x_tmp) ) / sinphi
        kmesh_data%input_proj%x(:, loop) = (proj_x_tmp(:) - cosphi*kmesh_data%input_proj%z(:, loop))/sinphi

        ! Final check
555     cosphi_new = sum(kmesh_data%input_proj%z(:, loop)*kmesh_data%input_proj%x(:, loop))
        if (abs(cosphi_new) .gt. eps6) then
          write (stdout, *) ' Projection:', loop
          call io_error(' Error: z and x axes are still not orthogonal after projection')
        endif

      endif

    enddo

    return

101 call io_error('param_get_projection: Problem reading l state into integer '//trim(ctemp3))
102 call io_error('param_get_projection: Problem reading m state into integer '//trim(ctemp3))
104 call io_error('param_get_projection: Problem reading zona into real '//trim(ctemp))
105 call io_error('param_get_projection: Problem reading radial state into integer '//trim(ctemp))
106 call io_error('param_get_projection: Problem reading m state into string '//trim(ctemp3))

  end subroutine param_get_projections

!===================================!
  subroutine param_get_keyword_kpath
    !===================================!
    !                                   !
    !!  Fills the kpath data block
    !                                   !
    !===================================!
    use w90_io, only: io_error

    implicit none

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
        call io_error('Error: Found '//trim(start_st)//' more than once in input file')
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
        call io_error('Error: Found '//trim(end_st)//' more than once in input file')
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    counter = 0
    do loop = line_s + 1, line_e - 1

      counter = counter + 2
      dummy = in_data(loop)
      read (dummy, *, err=240, end=240) spec_points%bands_label(counter - 1), &
        (spec_points%bands_spec_points(i, counter - 1), i=1, 3), &
        spec_points%bands_label(counter), (spec_points%bands_spec_points(i, counter), i=1, 3)
    end do

    in_data(line_s:line_e) (1:maxlen) = ' '

    return

240 call io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy))

  end subroutine param_get_keyword_kpath

end module w90_param_methods
