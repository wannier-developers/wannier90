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
!  w90_postw90_types: data types specific to postw90.x       !
!                                                            !
!------------------------------------------------------------!

module w90_postw90_types

  ! Definition of types encapsulating various quantities, data and parameters.
  ! Variables are grouped according to physical meaning and their use in the Wannier90 project.
  !
  !! Here are defined types specific to postw90.x (not used by wannier90.x).
  !! Types used by both wannier90.x and postw90.x are defined in types.F90.
  !! Types specific to wannier90.x (not used by postw90.x) are defined in wannier90_types.F90.

  use w90_constants, only: dp
  use w90_comms, only: w90comm_type

  implicit none

  public

  type pw90_calculation_type
    !!==================================================
    !! Contains information about the high-level task that pw90 is being asked to do,
    !! including any global variables that affect all calculation branches.
    !!==================================================
    logical :: kpath = .false.
    logical :: kslice = .false.
    logical :: dos = .false.
    logical :: berry = .false.
    logical :: gyrotropic = .false.
    logical :: geninterp = .false.
    logical :: boltzwann = .false.
    logical :: spin_moment = .false.
    logical :: spin_decomp = .false.
  end type pw90_calculation_type

  type pw90_oper_read_type
    !!==================================================
    !! Contains variables for determining whether formatted or unformatted input is read by get_oper.F90
    !!==================================================
    logical :: spn_formatted = .false.
    !! Read the spin from fortran formatted file
    logical :: uHu_formatted = .false.
    !! Read the uHu from fortran formatted file
  end type pw90_oper_read_type

  type kmesh_spacing_type
    real(kind=dp) :: spacing = -1._dp
    integer :: mesh(3) = 0
  end type kmesh_spacing_type

  ! Module  s p i n
  type pw90_spin_mod_type
    !!==================================================
    !! Contains variables used in the spin module of postw90
    !!==================================================
    real(kind=dp) :: axis_polar = 0.0_dp
    real(kind=dp) :: axis_azimuth = 0.0_dp
    type(kmesh_spacing_type) :: kmesh
  end type pw90_spin_mod_type

  ! REVIEW_2021-08-09: A general comment, which we record here so we don't forget.
  ! REVIEW_2021-08-09: It would be good to try to keep type and variable names
  ! REVIEW_2021-08-09: as consistent as possible. As an example, at the top of berry.F90, we have
  ! REVIEW_2021-08-09: type(kmesh_info_type), intent(in) :: kmesh_info   ! consistent
  ! REVIEW_2021-08-09: type(kpoint_dist_type), intent(in) :: kdist       ! not consistent
  ! REVIEW_2021-08-09: type(k_point_type), intent(in) :: k_points        ! fixed
  ! REVIEW_2021-08-09: This would make reading and debugging the code easier in future.

  type pw90_band_deriv_degen_type
    !!==================================================
    !! Contains variables for doing degenerate perturbation theory when bands are degenerate
    !! and band derivatives are needed.
    !!==================================================
    logical :: use_degen_pert = .false.
    real(kind=dp) :: degen_thr = 1.0d-4
  end type pw90_band_deriv_degen_type

  ! module  k p a t h (used by postw90/kpath)
  type pw90_kpath_mod_type
    !!==================================================
    !! Contains control variables for the kpath module of postw90
    !!==================================================
    character(len=20) :: task = 'bands'
    integer :: num_points = 100
    character(len=20) :: bands_colour = 'none'
  end type pw90_kpath_mod_type

  ! module  k s l i c e (postw90/kslice)
  type pw90_kslice_mod_type
    !!==================================================
    !! Contains control variables for the kslice module of postw90
    !!==================================================
    character(len=20) :: task = 'fermi_lines'
    real(kind=dp) :: corner(3) != 0.0_dp
    real(kind=dp) :: b1(3) != [1.0_dp, 0.0_dp, 0.0_dp]
    real(kind=dp) :: b2(3) ! = [0.0_dp, 1.0_dp, 0.0_dp]
    integer :: kmesh2d(2) ! = 50
    character(len=20) :: fermi_lines_colour = 'none'
  end type pw90_kslice_mod_type

  type pw90_smearing_type
    !!==================================================
    !! Contains variables for controlling the smearing.
    !!==================================================
    logical    :: use_adaptive = .true.
    real(kind=dp)    :: adaptive_prefactor = sqrt(2.0_dp)
    integer    :: type_index = 0 ! Gaussian default
    real(kind=dp)    :: fixed_width = 0.0_dp !none
    real(kind=dp)    :: adaptive_max_width = 1.0_dp ! 1 eV
    ! REVIEW_2021-08-09: Is this a speed-up that could be applied more generally?
    ! BGS currently only implemented in gyrotropic
    real(kind=dp) :: max_arg
  end type pw90_smearing_type

  type pw90_dos_mod_type
    !!==================================================
    !! Contains variables for the dos module of postw90
    !!==================================================
    character(len=20)    :: task = 'dos_plot'
    type(pw90_smearing_type) :: smearing ! default from global version
    real(kind=dp)    :: energy_max
    real(kind=dp)    :: energy_min
    real(kind=dp)    :: energy_step = 0.01_dp
    integer    :: num_project
    integer, allocatable :: project(:)
    !  character(len=20)    :: plot_format
    type(kmesh_spacing_type) :: kmesh
    !  real(kind=dp) :: gaussian_width
  end type pw90_dos_mod_type

  ! Module  b e r r y (mainly postw90/berry)
  type pw90_berry_mod_type
    !!==================================================
    !! Contains variables for the berry module of postw90
    !!==================================================
    character(len=120) :: task = ' '
    type(kmesh_spacing_type) :: kmesh
    integer :: curv_adpt_kmesh = 1
    real(kind=dp) :: curv_adpt_kmesh_thresh = 100.0_dp
    character(len=20) :: curv_unit = 'ang2' ! postw90/kpath, kslice as well
    type(pw90_smearing_type) :: kubo_smearing ! defaults to global values
    integer :: sc_phase_conv = 1
    real(kind=dp) :: sc_eta = 0.04 ! also postw90/wan_ham
    real(kind=dp) :: sc_w_thr = 5.0d0
    logical :: sc_use_eta_corr = .true.
    logical :: wanint_kpoint_file = .false. ! also postw90/spin, postw90/dos, postw90.F90
    logical :: transl_inv = .false. !also used in postw90/get_oper, postw90/gyrotropic
    real(kind=dp) :: kdotp_kpoint(3) = 0.0_dp
    integer, allocatable :: kdotp_bands(:)
    integer :: kubo_nfreq
    complex(kind=dp), allocatable :: kubo_freq_list(:)
    real(kind=dp) :: kubo_eigval_max
  end type pw90_berry_mod_type

  ! spin Hall conductivity (postw90 - common, get_oper, berry, kpath)
  type pw90_spin_hall_type
    !!==================================================
    !! Contains variables controlling the calculation of spin hall conductivity in postw90
    !!==================================================
    logical :: freq_scan = .false.
    integer :: alpha = 1
    integer :: beta = 2
    integer :: gamma = 3
    logical :: bandshift = .false.
    integer :: bandshift_firstband = 0
    real(kind=dp) :: bandshift_energyshift = 0._dp
    character(len=120) :: method = ' '
  end type pw90_spin_hall_type

  type pw90_gyrotropic_type ! postw90 - common, gyrotropic
    !!==================================================
    !! Contains variables for the gyrotropic module of postw90
    !!==================================================
    character(len=120) :: task = 'all'
    type(kmesh_spacing_type) :: kmesh
    integer :: nfreq
    complex(kind=dp), allocatable :: freq_list(:)
    real(kind=dp) :: box_corner(3) = 0.0_dp, box(3, 3) = 0.0
    real(kind=dp) :: degen_thresh = 0.0_dp
    integer, allocatable :: band_list(:)
    integer :: num_bands
    type(pw90_smearing_type) :: smearing ! adaptive is .false.
    real(kind=dp) :: eigval_max
  end type pw90_gyrotropic_type

  ! [gp-begin, Jun 1, 2012]
  ! GeneralInterpolator variables - postw90/geninterp
  type pw90_geninterp_mod_type
    !!==================================================
    !! Contains variables for the geninterp module of postw90
    !!==================================================
    logical :: alsofirstder = .false.
    logical :: single_file = .true.
  end type pw90_geninterp_mod_type
  ! [gp-end, Jun 1, 2012]

  ! [gp-begin, Apr 12, 2012]
  ! BoltzWann variables (postw90/boltzwann.F90)
  type pw90_boltzwann_type
    !!==================================================
    !! Contains variables for the boltzwann module of postw90
    !!==================================================
    logical :: calc_also_dos = .false.
    integer :: dir_num_2d = 0
    real(kind=dp) :: dos_energy_step = 0.001_dp
    real(kind=dp) :: dos_energy_min = -1.0_dp
    real(kind=dp) :: dos_energy_max = 0.0_dp
    type(pw90_smearing_type) :: dos_smearing ! default to global smearing values
    real(kind=dp) :: mu_min = -999._dp
    real(kind=dp) :: mu_max = -999._dp
    real(kind=dp) :: mu_step = 0._dp
    real(kind=dp) :: temp_min = -999._dp
    real(kind=dp) :: temp_max = -999._dp
    real(kind=dp) :: temp_step = 0._dp
    type(kmesh_spacing_type) :: kmesh
    real(kind=dp) :: tdf_energy_step = 0.001_dp
    ! adaptive smearing .false.
    type(pw90_smearing_type) :: tdf_smearing ! TDF_smr_index and TDF_smr_fixed_en_width
    real(kind=dp) :: relax_time = 10._dp
    logical :: bandshift = .false.
    integer :: bandshift_firstband = 0
    real(kind=dp) :: bandshift_energyshift = 0._dp
  end type pw90_boltzwann_type
  ! [gp-end, Apr 12, 2012]

  ! Parameters describing the direct lattice points R on a
  ! Wigner-Seitz supercell
  ! these were in postw90_common
  !
  type wigner_seitz_type
    integer, allocatable       :: irvec(:, :)
    real(kind=dp), allocatable :: crvec(:, :)
    integer, allocatable       :: ndegen(:)
    integer                    :: nrpts
    integer                    :: rpt_origin
  end type wigner_seitz_type

  type kpoint_dist_type ! kpoints from file
    integer                       :: max_int_kpts_on_node
    integer                       :: num_int_kpts
    integer, allocatable          :: num_int_kpts_on_node(:)
    real(kind=dp), allocatable    :: int_kpts(:, :), weight(:)
  end type kpoint_dist_type

end module w90_postw90_types
