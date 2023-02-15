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
!  w90_types: derived types encapsulating data required by   !
!      both wannier90.x and postw90.x                        !
!                                                            !
!------------------------------------------------------------!

module w90_types

  !! Definition of types encapsulating various quantities, data and parameters.
  !! Variables are grouped according to physical meaning and their use in the Wannier90 project.
  !!
  !! Here are defined types used by both wannier90.x and postw90.x.
  !! Types specific to wannier90.x (not used by postw90.x) are defined in wannier90_types.F90.
  !! Types specific to postw90.x (not used by wannier90.x) are defined in postw90/postw90_types.F90.

  use w90_constants, only: dp, maxlen

  implicit none

  public

  type print_output_type
    !!==================================================
    !! Contains variables to control output file formatting and verbosity.
    !!==================================================
    ! verbosity flags - w90_readwrite_read_verbosity
    integer :: iprint = 1
    ! Controls the verbosity of the output
    integer :: timing_level = 1
    ! REVIEW_2021-07-22: we agree that we don't need both length_unit and lenconfac;
    ! REVIEW_2021-07-22: instead could have a utility function.
    character(len=20) :: length_unit = 'ang' ! MAYBE, just have a separate variable?
    ! Units for length
    real(kind=dp) :: lenconfac = 1.0_dp !lots of write statements in wannier90
  end type print_output_type

  type w90_system_type
    !!==================================================
    !! Contains physical information about the material being calculated.
    !!==================================================
    integer :: num_valence_bands !**no sensibe default**
    integer :: num_elec_per_state = 2 ! used in: wannierise and postw90 dos and boltzwann
    logical :: spinors = .false.  !are our WF spinors? !kmesh, plot, wannier_lib, postw90/gyrotropic
  end type w90_system_type

  ! timer from io.F90
  type timing_data_type                                          !! Data about each stopwatch - for timing routines
    integer :: ncalls                                            !! Number of times stopwatch has been called
    real(kind=DP) :: ctime                                       !! Total time on stopwatch
    real(kind=DP) :: ptime                                       !! Temporary record of time when watch is started
    character(len=60) :: label                                   !! What is this stopwatch timing
  end type timing_data_type

  integer, parameter :: nmax = 100                               !! Maximum number of stopwatches
  type timer_list_type
    type(timing_data_type) :: clocks(nmax)                         !! Data for the stopwatches
    integer :: nnames = 0                                          !! Number of active stopwatches
    logical :: overflow = .false.
  end type timer_list_type

  type ws_region_type
    logical :: use_ws_distance = .true. !ws_distance, plot and postw90_common
    real(kind=dp) :: ws_distance_tol = 1.e-5_dp !ws_distance, hamiltonian and postw90_common
    !! absolute tolerance for the distance to equivalent positions
    integer :: ws_search_size(3) = 2 ! ws_distance, hamiltonian
    !! maximum extension in each direction of the supercell of the BvK cell
    !! to search for points inside the Wigner-Seitz cell
  end type ws_region_type

  ! used by ws_distance
  type ws_distance_type
    integer, allocatable :: irdist(:, :, :, :, :)!(3,ndegenx,num_wann,num_wann,nrpts)
    !! The integer number of unit cells to shift Wannier function j to put its centre
    !! inside the Wigner-Seitz of wannier function i. If several shifts are
    !! equivalent (i.e. they take the function on the edge of the WS) they are
    !! all listed. First index: xyz, second index: number of degenerate shifts,
    !! third and fourth indices: i,j; fifth index: index on the R vector.
    real(DP), allocatable :: crdist(:, :, :, :, :)!(3,ndegenx,num_wann,num_wann,nrpts)
    !! Cartesian version of irdist_ws, in angstrom
    integer, allocatable :: ndeg(:, :, :)!(num_wann,num_wann,nrpts)
    !! The number of equivalent vectors for each set of (i,j,R) (that is, loops on
    !! the second index of irdist_ws(:,:,i,j,R) go from 1 to wdist_ndeg(i,j,R))
    !
    logical :: done = .false.
    !! Global variable to know if the properties were already calculated, and avoid
    !! recalculating them when the [[ws_translate_dist]] function is called multiple times
  end type ws_distance_type

  ! setup in wannierise, but used by plot, ws_distance etc
  type wannier_data_type
    !!==================================================
    !! Contains the centres and spreads of the MLWFs
    !!==================================================
    ! Wannier centres and spreads
    real(kind=dp), allocatable :: centres(:, :)
    real(kind=dp), allocatable :: spreads(:)
    ! REVIEW_2021-07-22: Do we want to expose other related variables such as the decomposition
    ! REVIEW_2021-07-22: of the spread, matrix elements of r and r^2, etc. (TO FINISH)
  end type wannier_data_type

  ! The maximum number of shells we need to satisfy B1 condition in kmesh
  integer, parameter :: max_shells = 6
  integer, parameter :: num_nnmax = 12

  type kmesh_input_type
    !!==================================================
    !! Contains information that can be provided by the user about determining the kmesh
    !!==================================================
    integer :: num_shells = 0
    !! no longer an input keyword
    logical :: skip_B1_tests = .false.
    !! do not check the B1 condition
    integer, allocatable :: shell_list(:)
    integer :: search_shells = 36
    real(kind=dp) :: tol = 0.000001_dp
  end type kmesh_input_type

  !AAM: There are a number of ways one can handle the initial guess. (i) specify explicit
  !AAM: projections; (ii) use random (s-orbital) projections; (iii) a combination of (i) and
  !AAM: (ii); (iv) use the phases from the Bloch functions directly; (v) SCDM method. (i), (ii)
  !AAM: and (iii) require the arrays defined in the "projection_type" below. (iv) and (v) do not.
  !AAM: (vi) An external code may also simply supply to w90 an Amn(k) matrix that is has independently
  !AAM: generated, in which case projection_type is not needed.
  !AAM: It makes sense to keep the projection sites separate from the projection_type data below.
  type proj_type
    !!==================================================
    !! Contains information that can be provided by the user about the projections
    !!==================================================
    ! REVIEW_2021-07-22: site(:,:) has dual usage: for projections and for guiding centres.
    ! REVIEW_2021-07-22: Make a new type for guiding centres that only contains sites.
    ! REVIEW_2021-07-22: In the future this can be logically distinct from the projection sites.
    ! REVIEW_2021-07-22: For now, when defining proj_input_type, also define sites inside the
    ! REVIEW_2021-07-22: new guiding centres type.

    ! regarding defaults: specific flow in readwrite means all values are assigned anyway
    ! using defaults here requires restructuring get_projections() (readwrite.F90)
    ! site, l, m don't have reasonable defaults
    integer :: l, m, s
    integer :: radial           != 1
    real(kind=dp) :: site(3)    != (/0.0_dp, 0.0_dp, 0.0_dp/)
    real(kind=dp) :: s_qaxis(3) != (/0.0_dp, 0.0_dp, 1.0_dp/)
    real(kind=dp) :: z(3)       != (/0.0_dp, 0.0_dp, 1.0_dp/)
    real(kind=dp) :: x(3)       != (/1.0_dp, 0.0_dp, 0.0_dp/)
    real(kind=dp) :: zona       != 1.0_dp
  end type proj_type

  ! kmesh information (set in kmesh)
  type kmesh_info_type
    !!==================================================
    !! Contains derived information about the kmesh
    !!==================================================
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

  ! this contains data which described the disentangled manifold, also used in postw90
  type dis_manifold_type
    !!==================================================
    !! Contains information about the manifold of states from which the MLWFs are to be disentangled.
    !!==================================================
    real(kind=dp) :: win_min = -huge(0.0_dp)
    !! lower bound of the disentanglement outer window
    real(kind=dp) :: win_max = huge(0.0_dp)
    !! upper bound of the disentanglement outer window
    real(kind=dp) :: froz_min = -1.0_dp
    !! lower bound of the disentanglement inner (frozen) window
    real(kind=dp) :: froz_max = 0.0_dp
    !! upper bound of the disentanglement inner (frozen) window
    logical :: frozen_states = .false.
    ! disentangle parameters
    ! Used by plot, hamiltonian, wannierise, postw90_common, get_oper - not read
    integer, allocatable :: ndimwin(:)
    logical, allocatable :: lwindow(:, :)
  end type dis_manifold_type

  ! Atom sites - often used in the write_* routines
  ! hamiltonian, wannierise, plot, transport, wannier_lib
  type atom_data_type
    !!==================================================
    !! Contains information about the atoms (and maybe the cell...) of the system being calculated.
    !!==================================================
    real(kind=dp), allocatable :: pos_cart(:, :, :)
    integer, allocatable :: species_num(:)
    character(len=maxlen), allocatable :: label(:)
    character(len=2), allocatable :: symbol(:)
    integer :: num_atoms = 0
    integer :: num_species = 0
  end type atom_data_type

  ! plot.F90 and postw90/kpath
  type kpoint_path_type
    !!==================================================
    !! Contains information that specifies the k-point path for plotting and other purposes.
    !! Note: The length of bands_label and the second index of bands_spec_points is twice the
    !! number of segments specified by the user. Each pair of special points defines a segment.
    !!==================================================
    integer :: num_points_first_segment = 100
    character(len=20), allocatable :: labels(:)
    real(kind=dp), allocatable :: points(:, :)
  end type kpoint_path_type

  type settings_data
    !!==================================================
    !! structure to hold a scalar and array settings
    !!==================================================
    ! for simplicity, consider arrays of different rank
    ! as different types; otherwise reshape, etc.
    character(len=:), allocatable :: keyword ! token
    character(len=:), allocatable :: txtdata ! text data item
    ! integer data
    integer, allocatable :: i1d(:)
    integer, allocatable :: i2d(:, :)
    integer :: idata
    ! logical data
    logical, allocatable :: l1d(:)
    logical :: ldata
    ! fp data
    real(kind=dp), allocatable :: r1d(:)
    real(kind=dp), allocatable :: r2d(:, :)
    real(kind=dp) :: rdata
  end type settings_data

  type settings_type
    !!==================================================
    !! structure to hold input var/values
    !! contents of .win file and settings set by library interface
    !!==================================================
    integer :: num_entries = 0, num_entries_max = 0 ! number of keywords stored and max
    type(settings_data), allocatable :: entries(:)
    ! data for processing input file
    integer :: num_lines
    character(len=maxlen), allocatable :: in_data(:) ! contents of .win file
  end type settings_type

end module w90_types
