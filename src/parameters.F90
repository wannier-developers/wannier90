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
  use w90_io, only: stdout, maxlen

  implicit none

  private

  !Input
  integer, public, save :: iprint
  !! Controls the verbosity of the output
  character(len=20), public, save :: energy_unit
  !! Units for energy
  character(len=20), public, save :: length_unit
  !! Units for length
  logical, public, save :: wvfn_formatted
  !! Read the wvfn from fortran formatted file
  logical, public, save :: spn_formatted
  !! Read the spin from fortran formatted file
  logical, public, save :: uHu_formatted
  logical, public, save :: berry_uHu_formatted
  !! Read the uHu from fortran formatted file
  integer, public, save :: spin
  !! Spin up=1 down=2
  integer, public, save :: num_bands
  !! Number of bands
  integer, public, save :: num_dump_cycles
  !! Number of steps before writing checkpoint
  integer, public, save :: num_print_cycles
  !! Number of steps between writing output
  integer, public, save :: slwf_num
  !! Number of objective Wannier functions (others excluded from spread functional)
  logical, public, save :: selective_loc
  !! Selective localization
  logical, public, save :: slwf_constrain
  !! Constrained centres
  real(kind=dp), allocatable, public, save :: ccentres_frac(:, :)
  real(kind=dp), allocatable, public, save :: ccentres_cart(:, :)
  real(kind=dp), public, save :: slwf_lambda
  !! Centre constraints for each Wannier function. Co-ordinates of centre constraint defaults
  !! to centre of trial orbital. Individual Lagrange multipliers, lambdas, default to global Lagrange multiplier.
  character(len=50), public, save :: devel_flag
  ! Adaptive vs. fixed smearing stuff [GP, Jul 12, 2012]
  ! Only internal, always use the local variables defined by each module
  ! that take this value as default
  logical                         :: adpt_smr
  real(kind=dp)                   :: adpt_smr_fac
  real(kind=dp)                   :: adpt_smr_max
  real(kind=dp)                   :: smr_fixed_en_width
  ! GP: added a flag to check if this is the first run of param_read in library mode or not
  logical, public, save :: library_param_read_first_pass
  !IVO
  logical, public, save :: spin_moment
  real(kind=dp), public, save :: spin_axis_polar
  real(kind=dp), public, save :: spin_axis_azimuth
  logical, public, save :: use_degen_pert
  real(kind=dp), public, save :: degen_thr
  logical, public, save :: spin_decomp
  integer, public, save :: num_valence_bands
  logical                                  :: found_fermi_energy
  real(kind=dp), public, save :: scissors_shift
  !IVO_END
  ! [gp-begin, Apr 20, 2012] Smearing type
  ! The prefactor is given with the above parameters smr_...
  ! This is an internal variable, obtained from the input string smr_type
  ! Only internal, always use the local variables defined by each module
  ! that take this value as default
  integer                          :: smr_index
  ! [gp-end]
  integer, allocatable, public, save :: exclude_bands(:)
  integer, public, save :: num_wann
  !! number of wannier functions
  integer, public, save :: mp_grid(3)
  !! Dimensions of the Monkhorst-Pack grid
!  logical,           public, save :: automatic_mp_grid
  logical, public, save :: gamma_only
  !! Use the special Gamma-point routines
  real(kind=dp), public, save :: dis_win_min
  !! lower bound of the disentanglement outer window
  real(kind=dp), public, save :: dis_win_max
  !! upper bound of the disentanglement outer window
  real(kind=dp), public, save :: dis_froz_min
  !! lower bound of the disentanglement inner (frozen) window
  real(kind=dp), public, save :: dis_froz_max
  !! upper bound of the disentanglement inner (frozen) window
  integer, public, save :: dis_num_iter
  !! number of disentanglement iteration steps
  real(kind=dp), public, save :: dis_mix_ratio
  !! Mixing ratio for the disentanglement routine
  real(kind=dp), public, save :: dis_conv_tol
  !! Convergence tolerance for the disentanglement
  integer, public, save :: dis_conv_window
  !! Size of the convergence window for disentanglement
  ! GS-start
  integer, public, save :: dis_spheres_first_wann
  integer, public, save :: dis_spheres_num
  real(kind=dp), allocatable, public, save :: dis_spheres(:, :)
  ! GS-end
  integer, public, save :: num_iter
  !! Number of wannierisation iterations
  integer, public, save :: num_cg_steps
  !! Number of Conjugate Gradient steps
  real(kind=dp), public, save :: conv_tol
  integer, public, save :: conv_window
  logical, public, save :: wannier_plot
  integer, allocatable, public, save :: wannier_plot_list(:)
  integer, public, save :: wannier_plot_supercell(3)
  character(len=20), public, save :: wannier_plot_format
  character(len=20), public, save :: wannier_plot_mode
  character(len=20), public, save :: wannier_plot_spinor_mode
  logical, public, save :: wannier_plot_spinor_phase
  logical, public, save :: write_u_matrices
  logical, public, save :: bands_plot
  logical, public, save :: write_bvec
  integer, public, save :: bands_num_points
  character(len=20), public, save :: bands_plot_format
  character(len=20), public, save :: bands_plot_mode
  integer, allocatable, public, save :: bands_plot_project(:)
  integer, public, save :: bands_plot_dim
  logical, public, save :: write_hr
  logical, public, save :: write_rmn
  logical, public, save :: write_tb
  real(kind=dp), public, save :: hr_cutoff
  real(kind=dp), public, save :: dist_cutoff
  character(len=20), public, save :: dist_cutoff_mode
  real(kind=dp), public, save :: dist_cutoff_hc
  character(len=20), public, save :: one_dim_axis
  logical, public, save :: use_ws_distance
  real(kind=dp), public, save :: ws_distance_tol
  !! absolute tolerance for the distance to equivalent positions
  integer, public, save :: ws_search_size(3)
  !! maximum extension in each direction of the supercell of the BvK cell
  !! to search for points inside the Wigner-Seitz cell
  logical, public, save :: fermi_surface_plot
  integer, public, save :: fermi_surface_num_points
  character(len=20), public, save :: fermi_surface_plot_format
  real(kind=dp), save :: fermi_energy

  ! module  k p a t h
  logical, public, save :: kpath
  character(len=20), public, save :: kpath_task
  integer, public, save :: kpath_num_points
  character(len=20), public, save :: kpath_bands_colour

  ! module  k s l i c e
  logical, public, save :: kslice
  character(len=20), public, save :: kslice_task
  real(kind=dp), public, save :: kslice_corner(3)
  real(kind=dp), public, save :: kslice_b1(3)
  real(kind=dp), public, save :: kslice_b2(3)
  integer, public, save :: kslice_2dkmesh(2)
  character(len=20), public, save :: kslice_fermi_lines_colour

  ! module  d o s
  logical, public, save    :: dos
! No need to save 'dos_plot', only used here (introduced 'dos_task')
  logical, public          :: dos_plot
  character(len=20), public, save    :: dos_task
  logical, public, save    :: dos_adpt_smr
  real(kind=dp), public, save    :: dos_adpt_smr_fac
  integer, public, save    :: dos_smr_index
  real(kind=dp), public, save    :: dos_smr_fixed_en_width
  real(kind=dp), public, save    :: dos_adpt_smr_max
  real(kind=dp), public, save    :: dos_energy_max
  real(kind=dp), public, save    :: dos_energy_min
  real(kind=dp), public, save    :: dos_energy_step
  integer, public, save    :: num_dos_project
  integer, allocatable, public, save :: dos_project(:)
!  character(len=20), public, save    :: dos_plot_format
  real(kind=dp), public, save    :: dos_kmesh_spacing
  integer, public, save    :: dos_kmesh(3)
!  real(kind=dp),     public, save :: dos_gaussian_width

! Module  b e r r y
  logical, public, save :: berry
  character(len=120), public, save :: berry_task
  real(kind=dp), public, save :: berry_kmesh_spacing
  integer, public, save :: berry_kmesh(3)
  ! --------------remove eventually----------------
!  integer,           public, save :: alpha
!  integer,           public, save :: beta
!  integer,           public, save :: gamma
  ! --------------remove eventually----------------
  integer, public, save :: berry_curv_adpt_kmesh
  real(kind=dp), public, save :: berry_curv_adpt_kmesh_thresh
  character(len=20), public, save :: berry_curv_unit
  logical, public, save :: kubo_adpt_smr
  real(kind=dp), public, save :: kubo_adpt_smr_fac
  integer, public, save :: kubo_smr_index
  real(kind=dp), public, save :: kubo_smr_fixed_en_width
  real(kind=dp), public, save :: kubo_adpt_smr_max
  integer, public, save :: sc_phase_conv
  real(kind=dp), public, save :: sc_eta
  real(kind=dp), public, save :: sc_w_thr
  logical, public, save :: wanint_kpoint_file
!  logical,           public, save :: sigma_abc_onlyorb
  logical, public, save :: transl_inv

  ! spin Hall conductivity
  logical, public, save :: shc_freq_scan
  integer, public, save :: shc_alpha
  integer, public, save :: shc_beta
  integer, public, save :: shc_gamma
  logical, public, save :: shc_bandshift
  integer, public, save :: shc_bandshift_firstband
  real(kind=dp), public, save :: shc_bandshift_energyshift

  logical, public, save :: gyrotropic
  character(len=120), public, save :: gyrotropic_task
  integer, public, save :: gyrotropic_kmesh(3)
  real(kind=dp), public, save :: gyrotropic_kmesh_spacing
  integer, public, save :: gyrotropic_smr_index
  real(kind=dp), public, save :: gyrotropic_smr_fixed_en_width
  real(kind=dp)                               :: gyrotropic_freq_min
  real(kind=dp)                               :: gyrotropic_freq_max
  real(kind=dp)                               :: gyrotropic_freq_step
  integer, public, save :: gyrotropic_nfreq
  complex(kind=dp), allocatable, public, save :: gyrotropic_freq_list(:)
  real(kind=dp), public, save :: gyrotropic_box_corner(3), gyrotropic_box(3, 3)
  real(kind=dp)                   :: gyrotropic_box_tmp(3)
  real(kind=dp), public, save :: gyrotropic_degen_thresh
  integer, allocatable, public, save :: gyrotropic_band_list(:)
  integer, public, save :: gyrotropic_num_bands
  real(kind=dp)                   :: smr_max_arg
  real(kind=dp), public, save :: gyrotropic_smr_max_arg
  real(kind=dp), public, save :: gyrotropic_eigval_max

  logical                                  :: fermi_energy_scan
  real(kind=dp)                            :: fermi_energy_min
  real(kind=dp)                            :: fermi_energy_max
  real(kind=dp)                            :: fermi_energy_step
  integer, public, save :: nfermi
  real(kind=dp), allocatable, public, save :: fermi_energy_list(:)

  real(kind=dp)                               :: kubo_freq_min
  real(kind=dp)                               :: kubo_freq_max
  real(kind=dp)                               :: kubo_freq_step
  integer, public, save :: kubo_nfreq
  complex(kind=dp), allocatable, public, save :: kubo_freq_list(:)
  real(kind=dp), public, save :: kubo_eigval_max

! Module  s p i n
  real(kind=dp), public, save :: spin_kmesh_spacing
  integer, public, save :: spin_kmesh(3)

  ! [gp-begin, Apr 13, 2012]
  ! Global interpolation k mesh variables
  ! These don't need to be public, since their values are copied in the variables of the
  ! local interpolation meshes. JRY: added save attribute
  real(kind=dp), save             :: kmesh_spacing
  integer, save                   :: kmesh(3)
  logical, save                   :: global_kmesh_set
  ! [gp-end]

  ! [gp-begin, Jun 1, 2012]
  ! GeneralInterpolator variables
  logical, public, save :: geninterp
  logical, public, save :: geninterp_alsofirstder
  logical, public, save :: geninterp_single_file
  ! [gp-end, Jun 1, 2012]

  ! [gp-begin, Apr 12, 2012]
  ! BoltzWann variables
  logical, public, save :: boltzwann
  logical, public, save :: boltz_calc_also_dos
  integer, public, save :: boltz_2d_dir_num
  character(len=4), save :: boltz_2d_dir
  real(kind=dp), public, save :: boltz_dos_energy_step
  real(kind=dp), public, save :: boltz_dos_energy_min
  real(kind=dp), public, save :: boltz_dos_energy_max
  logical, public, save :: boltz_dos_adpt_smr
  real(kind=dp), public, save :: boltz_dos_smr_fixed_en_width
  real(kind=dp), public, save :: boltz_dos_adpt_smr_fac
  real(kind=dp), public, save :: boltz_dos_adpt_smr_max
  real(kind=dp), public, save :: boltz_mu_min
  real(kind=dp), public, save :: boltz_mu_max
  real(kind=dp), public, save :: boltz_mu_step
  real(kind=dp), public, save :: boltz_temp_min
  real(kind=dp), public, save :: boltz_temp_max
  real(kind=dp), public, save :: boltz_temp_step
  real(kind=dp), public, save :: boltz_kmesh_spacing
  integer, public, save :: boltz_kmesh(3)
  real(kind=dp), public, save :: boltz_tdf_energy_step
  integer, public, save :: boltz_TDF_smr_index
  integer, public, save :: boltz_dos_smr_index
  real(kind=dp), public, save :: boltz_relax_time
  real(kind=dp), public, save :: boltz_TDF_smr_fixed_en_width
  logical, public, save :: boltz_bandshift
  integer, public, save :: boltz_bandshift_firstband
  real(kind=dp), public, save :: boltz_bandshift_energyshift
  ! [gp-end, Apr 12, 2012]

  logical, public, save :: transport
  logical, public, save :: tran_easy_fix
  character(len=20), public, save :: transport_mode
  real(kind=dp), public, save :: tran_win_min
  real(kind=dp), public, save :: tran_win_max
  real(kind=dp), public, save :: tran_energy_step
  integer, public, save :: tran_num_bb
  integer, public, save :: tran_num_ll
  integer, public, save :: tran_num_rr
  integer, public, save :: tran_num_cc
  integer, public, save :: tran_num_lc
  integer, public, save :: tran_num_cr
  integer, public, save :: tran_num_bandc
  logical, public, save :: tran_write_ht
  logical, public, save :: tran_read_ht
  logical, public, save :: tran_use_same_lead
  integer, public, save :: tran_num_cell_ll
  integer, public, save :: tran_num_cell_rr
  real(kind=dp), public, save :: tran_group_threshold
  real(kind=dp), public, save :: translation_centre_frac(3)
  integer, public, save :: num_shells
  !! no longer an input keyword
  logical, public, save :: skip_B1_tests
  !! do not check the B1 condition
  logical, public, save :: explicit_nnkpts
  !! nnkpts block is in the input file (allowed only for post-proc setup)
  integer, allocatable, public, save :: shell_list(:)
  real(kind=dp), allocatable, public, save :: kpt_latt(:, :)
  !! kpoints in lattice vecs
  real(kind=dp), public, save :: real_lattice(3, 3)
  logical, public, save :: postproc_setup
  logical, public, save :: cp_pp
  !! Car-Parinello post-proc flag/transport

  logical, public, save :: calc_only_A
  logical, public, save :: use_bloch_phases
  character(len=20), public, save :: restart
  logical, public, save :: write_r2mn
  logical, public, save :: guiding_centres
  integer, public, save :: num_guide_cycles
  integer, public, save :: num_no_guide_iter
  real(kind=dp), public, save :: fixed_step
  real(kind=dp), public, save :: trial_step
  logical, public, save :: precond
  logical, public, save :: write_proj
  integer, public, save :: timing_level
  logical, public, save :: spinors   !are our WF spinors?
  integer, public, save :: num_elec_per_state
  logical, public, save :: translate_home_cell
  logical, public, save :: write_xyz
  logical, public, save :: write_hr_diag
  real(kind=dp), public, save :: conv_noise_amp
  integer, public, save :: conv_noise_num
  real(kind=dp), public, save :: wannier_plot_radius
  real(kind=dp), public, save :: wannier_plot_scale
  integer, public, save :: search_shells   !for kmesh
  real(kind=dp), public, save :: kmesh_tol
  integer, public, save :: optimisation
  ! aam: for WF-based calculation of vdW C6 coefficients
  logical, public, save :: write_vdw_data

  ! Restarts
  real(kind=dp), public, save :: omega_invariant
  character(len=20), public, save :: checkpoint
  logical, public, save :: have_disentangled

  ! Atom sites
  real(kind=dp), allocatable, public, save :: atoms_pos_frac(:, :, :)
  real(kind=dp), allocatable, public, save :: atoms_pos_cart(:, :, :)
  integer, allocatable, public, save :: atoms_species_num(:)
  character(len=maxlen), allocatable, public, save :: atoms_label(:)
  character(len=2), allocatable, public, save :: atoms_symbol(:)
  integer, public, save :: num_atoms
  integer, public, save :: num_species

  ! Projections
  logical, public, save :: lhasproj
  real(kind=dp), allocatable, public, save :: input_proj_site(:, :)
  integer, allocatable, public, save :: input_proj_l(:)
  integer, allocatable, public, save :: input_proj_m(:)
  integer, allocatable, public, save :: input_proj_s(:)
  real(kind=dp), allocatable, public, save :: input_proj_s_qaxis(:, :)
  real(kind=dp), allocatable, public, save :: input_proj_z(:, :)
  real(kind=dp), allocatable, public, save :: input_proj_x(:, :)
  integer, allocatable, public, save :: input_proj_radial(:)
  real(kind=dp), allocatable, public, save :: input_proj_zona(:)
  real(kind=dp), allocatable, public, save :: proj_site(:, :)
  integer, allocatable, public, save :: proj_l(:)
  integer, allocatable, public, save :: proj_m(:)
  integer, allocatable, public, save :: proj_s(:)
  real(kind=dp), allocatable, public, save :: proj_s_qaxis(:, :)
  real(kind=dp), allocatable, public, save :: proj_z(:, :)
  real(kind=dp), allocatable, public, save :: proj_x(:, :)
  integer, allocatable, public, save :: proj_radial(:)
  real(kind=dp), allocatable, public, save :: proj_zona(:)
  integer, public, save :: num_proj
  ! projections selection
  logical, public, save :: lselproj
  integer, public, save :: num_select_projections
  integer, allocatable, public, save :: select_projections(:)
  integer, allocatable, public, save :: proj2wann_map(:)
  ! a u t o m a t i c   p r o j e c t i o n s
  ! vv: Writes a new block in .nnkp
  logical, public, save :: auto_projections

  !parameters dervied from input
  integer, public, save :: num_kpts
  real(kind=dp), public, save :: recip_lattice(3, 3)
  real(kind=dp), public, save :: cell_volume
  real(kind=dp), public, save :: real_metric(3, 3)
  real(kind=dp), public, save :: recip_metric(3, 3)
  integer, public, save :: bands_num_spec_points
  character(len=20), allocatable, public, save ::bands_label(:)
  real(kind=dp), allocatable, public, save ::bands_spec_points(:, :)
  real(kind=dp), allocatable, public, save ::kpt_cart(:, :) !kpoints in cartesians
  logical, public, save :: disentanglement
  real(kind=dp), public, save :: lenconfac
  integer, public, save :: num_wannier_plot
  integer, public, save :: num_bands_project
  integer, public, save :: num_exclude_bands
  logical, public, save :: lfixstep

  ! kmesh parameters (set in kmesh)

  integer, public, save              :: nnh           ! the number of b-directions (bka)
  integer, public, save              :: nntot         ! total number of neighbours for each k-point
  integer, public, save, allocatable :: nnlist(:, :)   ! list of neighbours for each k-point
  integer, public, save, allocatable :: neigh(:, :)
  integer, public, save, allocatable :: nncell(:, :, :) ! gives BZ of each neighbour of each k-point
  real(kind=dp), public, save              :: wbtot
  real(kind=dp), public, save, allocatable :: wb(:)         ! weights associated with neighbours of each k-point
  real(kind=dp), public, save, allocatable :: bk(:, :, :)     ! the b-vectors that go from each k-point to its neighbours
  real(kind=dp), public, save, allocatable :: bka(:, :)      ! the b-directions from 1st k-point to its neighbours

  ! disentangle parameters
  integer, public, save, allocatable :: ndimwin(:)
  logical, public, save, allocatable :: lwindow(:, :)
  logical, public, save :: frozen_states

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states

  complex(kind=dp), allocatable, save, public :: a_matrix(:, :, :)
  complex(kind=dp), allocatable, save, public :: m_matrix_orig(:, :, :, :)
  complex(kind=dp), allocatable, save, public :: m_matrix_orig_local(:, :, :, :)
  real(kind=dp), allocatable, save, public :: eigval(:, :)
  logical, save, public :: eig_found

! $![ysl-b]
! $  ! ph_g = phase factor of Bloch functions at Gamma
! $  !  assuming that Bloch functions at Gamma are real except this phase factor
! $  complex(kind=dp), allocatable, save, public :: ph_g(:)
! $![ysl-e]

  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states

  complex(kind=dp), allocatable, save, public :: u_matrix_opt(:, :, :)

  ! u_matrix gives the unitary rotations from the optimal subspace to the
  ! optimally smooth states.
  ! m_matrix we store here, becuase it is needed for restart of wannierise

  complex(kind=dp), allocatable, save, public :: u_matrix(:, :, :)
  complex(kind=dp), allocatable, save, public :: m_matrix(:, :, :, :)
  complex(kind=dp), allocatable, save, public :: m_matrix_local(:, :, :, :)

  ! RS: symmetry-adapted Wannier functions
  logical, public, save :: lsitesymmetry = .false.
  real(kind=dp), public, save :: symmetrize_eps = 1.d-3

  ! The maximum number of shells we need to satisfy B1 condition in kmesh
  integer, parameter, public :: max_shells = 6
  integer, parameter, public :: num_nnmax = 12

  ! Are we running as a library
  logical, save, public :: library = .false.

  ! Are we running postw90?
  logical, save, public :: ispostw90 = .false.

  ! IVO
  ! Are we running postw90 starting from an effective model?
  logical, save, public :: effective_model = .false.

  ! Wannier centres and spreads
  real(kind=dp), public, save, allocatable :: wannier_centres(:, :)
  real(kind=dp), public, save, allocatable :: wannier_spreads(:)
  real(kind=dp), public, save :: omega_total
  real(kind=dp), public, save :: omega_tilde
  ! [ omega_invariant is declared above ]

  ! For Hamiltonian matrix in WF representation
  logical, public, save              :: automatic_translation
  integer, public, save              :: one_dim_dir

  ! Private data
  integer                            :: num_lines
  character(len=maxlen), allocatable :: in_data(:)
  character(len=maxlen)              :: ctmp
  logical                            :: ltmp
  ! AAM_2016-09-15: hr_plot is a deprecated input parameter. Replaced by write_hr.
  logical                            :: hr_plot

  public :: param_read
  public :: param_write
  public :: param_postw90_write
  public :: param_dealloc
  public :: param_write_header
  public :: param_write_chkpt
  public :: param_read_chkpt
  public :: param_lib_set_atoms
  public :: param_memory_estimate
  public :: param_get_smearing_type
  public :: param_get_convention_type
  public :: param_dist
  public :: param_chkpt_dist

contains

  !==================================================================!
  subroutine param_read()
    !==================================================================!
    !                                                                  !
    !! Read parameters and calculate derived values
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !                                                                  !
    !===================================================================
    use w90_constants, only: bohr, eps6, cmplx_i
    use w90_utility, only: utility_recip_lattice, utility_metric
    use w90_io, only: io_error, io_file_unit, seedname, post_proc_flag
    implicit none

    !local variables
    real(kind=dp)  :: real_lattice_tmp(3, 3)
    integer :: nkp, i, j, n, k, itmp, i_temp, i_temp2, eig_unit, loop, ierr, iv_temp(3), rows
    logical :: found, found2, lunits, chk_found
    character(len=6) :: spin_str
    real(kind=dp) :: cosa(3), rv_temp(3)
    integer, allocatable, dimension(:, :) :: nnkpts_block
    integer, allocatable, dimension(:) :: nnkpts_idx

    call param_in_file

    !%%%%%%%%%%%%%%%%
    ! Site symmetry
    !%%%%%%%%%%%%%%%%

    ! default value is lsitesymmetry=.false.
    call param_get_keyword('site_symmetry', found, l_value=lsitesymmetry)!YN:

    ! default value is symmetrize_eps=0.001
    call param_get_keyword('symmetrize_eps', found, r_value=symmetrize_eps)!YN:

    !%%%%%%%%%%%%%%%%
    ! Transport
    !%%%%%%%%%%%%%%%%

    transport = .false.
    call param_get_keyword('transport', found, l_value=transport)

    tran_read_ht = .false.
    call param_get_keyword('tran_read_ht', found, l_value=tran_read_ht)

    tran_easy_fix = .false.
    call param_get_keyword('tran_easy_fix', found, l_value=tran_easy_fix)

    if (transport .and. tran_read_ht) restart = ' '

    !%%%%%%%%%%%%%%%%
    !System variables
    !%%%%%%%%%%%%%%%%

    timing_level = 1             ! Verbosity of timing output info
    call param_get_keyword('timing_level', found, i_value=timing_level)

    iprint = 1             ! Verbosity
    call param_get_keyword('iprint', found, i_value=iprint)

    optimisation = 3             ! Verbosity
    call param_get_keyword('optimisation', found, i_value=optimisation)

    if (transport .and. tran_read_ht) goto 301

    !ivo
    call param_get_keyword('effective_model', found, l_value=effective_model)

    energy_unit = 'ev'          !
    call param_get_keyword('energy_unit', found, c_value=energy_unit)

    length_unit = 'ang'         !
    lenconfac = 1.0_dp
    call param_get_keyword('length_unit', found, c_value=length_unit)
    if (length_unit .ne. 'ang' .and. length_unit .ne. 'bohr') &
      call io_error('Error: value of length_unit not recognised in param_read')
    if (length_unit .eq. 'bohr') lenconfac = 1.0_dp/bohr

    wvfn_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword('wvfn_formatted', found, l_value=wvfn_formatted)

    spn_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword('spn_formatted', found, l_value=spn_formatted)

    uHu_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword('uhu_formatted', found, l_value=uHu_formatted)

    spin = 1
    call param_get_keyword('spin', found, c_value=spin_str)
    if (found) then
      if (index(spin_str, 'up') > 0) then
        spin = 1
      elseif (index(spin_str, 'down') > 0) then
        spin = 2
      else
        call io_error('Error: unrecognised value of spin found: '//trim(spin_str))
      end if
    end if

    num_wann = -99
    call param_get_keyword('num_wann', found, i_value=num_wann)
    if (.not. found) call io_error('Error: You must specify num_wann')
    if (num_wann <= 0) call io_error('Error: num_wann must be greater than zero')

    num_exclude_bands = 0
    call param_get_range_vector('exclude_bands', found, num_exclude_bands, lcount=.true.)
    if (found) then
      if (num_exclude_bands < 1) call io_error('Error: problem reading exclude_bands')
      if (allocated(exclude_bands)) deallocate (exclude_bands)
      allocate (exclude_bands(num_exclude_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating exclude_bands in param_read')
      call param_get_range_vector('exclude_bands', found, num_exclude_bands, .false., exclude_bands)
      if (any(exclude_bands < 1)) &
        call io_error('Error: exclude_bands must contain positive numbers')
    end if

    ! AAM_2016-09-16: some changes to logic to patch a problem with uninitialised num_bands in library mode
!    num_bands       =   -1
    call param_get_keyword('num_bands', found, i_value=i_temp)
    if (found .and. library) write (stdout, '(/a)') ' Ignoring <num_bands> in input file'
    if (.not. library .and. .not. effective_model) then
      if (found) num_bands = i_temp
      if (.not. found) num_bands = num_wann
    end if
    ! GP: I subtract it here, but only the first time when I pass the total number of bands
    ! In later calls, I need to pass instead num_bands already subtracted.
    if (library .and. library_param_read_first_pass) num_bands = num_bands - num_exclude_bands
    if (.not. effective_model) then
      if (found .and. num_bands < num_wann) then
        write (stdout, *) 'num_bands', num_bands
        write (stdout, *) 'num_wann', num_wann
        call io_error('Error: num_bands must be greater than or equal to num_wann')
      endif
    endif

    num_dump_cycles = 100          ! frequency to write backups at
    call param_get_keyword('num_dump_cycles', found, i_value=num_dump_cycles)
    if (num_dump_cycles < 0) call io_error('Error: num_dump_cycles must be positive')

    num_print_cycles = 1          ! frequency to write at
    call param_get_keyword('num_print_cycles', found, i_value=num_print_cycles)
    if (num_print_cycles < 0) call io_error('Error: num_print_cycles must be positive')

    devel_flag = ' '          !
    call param_get_keyword('devel_flag', found, c_value=devel_flag)

!    mp_grid=-99
    call param_get_keyword_vector('mp_grid', found, 3, i_value=iv_temp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <mp_grid> in input file'
    if (.not. library .and. .not. effective_model) then
      if (found) mp_grid = iv_temp
      if (.not. found) then
        call io_error('Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid')
      elseif (any(mp_grid < 1)) then
        call io_error('Error: mp_grid must be greater than zero')
      end if
      num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
    end if

![ysl-b]
    ltmp = .false.
    call param_get_keyword('gamma_only', found, l_value=ltmp)
    if (.not. library) then
      gamma_only = ltmp
      if (gamma_only .and. (num_kpts .ne. 1)) &
        call io_error('Error: gamma_only is true, but num_kpts > 1')
    else
      if (found) write (stdout, '(a)') ' Ignoring <gamma_only> in input file'
    endif
![ysl-e]

!    aam: automatic_mp_grid no longer used
!    automatic_mp_grid = .false.
!    call param_get_keyword('automatic_mp_grid',found,l_value=automatic_mp_grid)

    postproc_setup = .false.            ! set to true to write .nnkp file and exit
    call param_get_keyword('postproc_setup', found, l_value=postproc_setup)
    ! We allow this keyword to be overriden by a command line arg -pp
    if (post_proc_flag) postproc_setup = .true.

    cp_pp = .false.                  ! set to true if doing CP post-processing
    call param_get_keyword('cp_pp', found, l_value=cp_pp)

    calc_only_A = .false.
    call param_get_keyword('calc_only_A', found, l_value=calc_only_A)

    restart = ' '
    call param_get_keyword('restart', found, c_value=restart)
    if (found) then
      if ((restart .ne. 'default') .and. (restart .ne. 'wannierise') &
          .and. (restart .ne. 'plot') .and. (restart .ne. 'transport')) then
        call io_error('Error in input file: value of restart not recognised')
      else
        inquire (file=trim(seedname)//'.chk', exist=chk_found)
        if (.not. chk_found) &
          call io_error('Error: restart requested but '//trim(seedname)//'.chk file not found')
      endif
    endif
    !post processing takes priority (user is not warned of this)
    if (postproc_setup) restart = ' '

    write_r2mn = .false.
    call param_get_keyword('write_r2mn', found, l_value=write_r2mn)

    write_proj = .false.
    call param_get_keyword('write_proj', found, l_value=write_proj)

    ltmp = .false.  ! by default our WF are not spinors
    call param_get_keyword('spinors', found, l_value=ltmp)
    if (.not. library) then
      spinors = ltmp
    else
      if (found) write (stdout, '(a)') ' Ignoring <spinors> in input file'
    endif
!    if(spinors .and. (2*(num_wann/2))/=num_wann) &
!       call io_error('Error: For spinor WF num_wann must be even')

    ! We need to know if the bands are double degenerate due to spin, e.g. when
    ! calculating the DOS
    if (spinors) then
      num_elec_per_state = 1
    else
      num_elec_per_state = 2
    endif
    call param_get_keyword('num_elec_per_state', found, i_value=num_elec_per_state)
    if ((num_elec_per_state /= 1) .and. (num_elec_per_state /= 2)) &
      call io_error('Error: num_elec_per_state can be only 1 or 2')
    if (spinors .and. num_elec_per_state /= 1) &
      call io_error('Error: when spinors = T num_elec_per_state must be 1')

    translate_home_cell = .false.
    call param_get_keyword('translate_home_cell', found, l_value=translate_home_cell)

    write_xyz = .false.
    call param_get_keyword('write_xyz', found, l_value=write_xyz)

    write_hr_diag = .false.
    call param_get_keyword('write_hr_diag', found, l_value=write_hr_diag)

    !%%%%%%%%%%%
    ! Wannierise
    !%%%%%%%%%%%

    num_iter = 100
    call param_get_keyword('num_iter', found, i_value=num_iter)
    if (num_iter < 0) call io_error('Error: num_iter must be positive')

    num_cg_steps = 5
    call param_get_keyword('num_cg_steps', found, i_value=num_cg_steps)
    if (num_cg_steps < 0) call io_error('Error: num_cg_steps must be positive')

    conv_tol = 1.0e-10_dp
    call param_get_keyword('conv_tol', found, r_value=conv_tol)
    if (conv_tol < 0.0_dp) call io_error('Error: conv_tol must be positive')

    conv_noise_amp = -1.0_dp
    call param_get_keyword('conv_noise_amp', found, r_value=conv_noise_amp)

    conv_window = -1
    if (conv_noise_amp > 0.0_dp) conv_window = 5
    call param_get_keyword('conv_window', found, i_value=conv_window)

    conv_noise_num = 3
    call param_get_keyword('conv_noise_num', found, i_value=conv_noise_num)
    if (conv_noise_num < 0) call io_error('Error: conv_noise_num must be positive')

    guiding_centres = .false.
    call param_get_keyword('guiding_centres', found, l_value=guiding_centres)

    num_guide_cycles = 1
    call param_get_keyword('num_guide_cycles', found, i_value=num_guide_cycles)
    if (num_guide_cycles < 0) call io_error('Error: num_guide_cycles must be >= 0')

    num_no_guide_iter = 0
    call param_get_keyword('num_no_guide_iter', found, i_value=num_no_guide_iter)
    if (num_no_guide_iter < 0) call io_error('Error: num_no_guide_iter must be >= 0')

    fixed_step = -999.0_dp; lfixstep = .false.
    call param_get_keyword('fixed_step', found, r_value=fixed_step)
    if (found .and. (fixed_step < 0.0_dp)) call io_error('Error: fixed_step must be > 0')
    if (fixed_step > 0.0_dp) lfixstep = .true.

    trial_step = 2.0_dp
    call param_get_keyword('trial_step', found, r_value=trial_step)
    if (found .and. lfixstep) call io_error('Error: cannot specify both fixed_step and trial_step')

    precond = .false.
    call param_get_keyword('precond', found, l_value=precond)

    slwf_num = num_wann
    selective_loc = .false.
    call param_get_keyword('slwf_num', found, i_value=slwf_num)
    if (found) then
      if (slwf_num .gt. num_wann .or. slwf_num .lt. 1) then
        call io_error('Error: slwf_num must be an integer between 1 and num_wann')
      end if
      if (slwf_num .lt. num_wann) selective_loc = .true.
    end if

    slwf_constrain = .false.
    call param_get_keyword('slwf_constrain', found, l_value=slwf_constrain)
    if (found .and. slwf_constrain) then
      if (selective_loc) then
        allocate (ccentres_frac(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_frac in param_get_centre_constraints')
        allocate (ccentres_cart(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_cart in param_get_centre_constraints')
      else
        write (stdout, *) ' No selective localisation requested. Ignoring constraints on centres'
        slwf_constrain = .false.
      end if
    end if

    slwf_lambda = 1.0_dp
    call param_get_keyword('slwf_lambda', found, r_value=slwf_lambda)
    if (found) then
      if (slwf_lambda < 0.0_dp) call io_error('Error: slwf_lambda  must be positive.')
    endif

    !%%%%%%%%%
    ! Plotting
    !%%%%%%%%%

    wannier_plot = .false.
    call param_get_keyword('wannier_plot', found, l_value=wannier_plot)

    wannier_plot_supercell = 2

    call param_get_vector_length('wannier_plot_supercell', found, length=i)
    if (found) then
      if (i .eq. 1) then
        call param_get_keyword_vector('wannier_plot_supercell', found, 1, &
                                      i_value=wannier_plot_supercell)
        wannier_plot_supercell(2) = wannier_plot_supercell(1)
        wannier_plot_supercell(3) = wannier_plot_supercell(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector('wannier_plot_supercell', found, 3, &
                                      i_value=wannier_plot_supercell)
      else
        call io_error('Error: wannier_plot_supercell must be provided as either one integer or a vector of three integers')
      end if
      if (any(wannier_plot_supercell <= 0)) &
        call io_error('Error: wannier_plot_supercell elements must be greater than zero')
    end if

    wannier_plot_format = 'xcrysden'
    call param_get_keyword('wannier_plot_format', found, c_value=wannier_plot_format)

    wannier_plot_mode = 'crystal'
    call param_get_keyword('wannier_plot_mode', found, c_value=wannier_plot_mode)

    wannier_plot_spinor_mode = 'total'
    call param_get_keyword('wannier_plot_spinor_mode', found, c_value=wannier_plot_spinor_mode)
    wannier_plot_spinor_phase = .true.
    call param_get_keyword('wannier_plot_spinor_phase', found, l_value=wannier_plot_spinor_phase)

    call param_get_range_vector('wannier_plot_list', found, num_wannier_plot, lcount=.true.)
    if (found) then
      if (num_wannier_plot < 1) call io_error('Error: problem reading wannier_plot_list')
      if (allocated(wannier_plot_list)) deallocate (wannier_plot_list)
      allocate (wannier_plot_list(num_wannier_plot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_plot_list in param_read')
      call param_get_range_vector('wannier_plot_list', found, num_wannier_plot, .false., wannier_plot_list)
      if (any(wannier_plot_list < 1) .or. any(wannier_plot_list > num_wann)) &
        call io_error('Error: wannier_plot_list asks for a non-valid wannier function to be plotted')
    else
      ! we plot all wannier functions
      num_wannier_plot = num_wann
      if (allocated(wannier_plot_list)) deallocate (wannier_plot_list)
      allocate (wannier_plot_list(num_wannier_plot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_plot_list in param_read')
      do loop = 1, num_wann
        wannier_plot_list(loop) = loop
      end do
    end if

    wannier_plot_radius = 3.5_dp
    call param_get_keyword('wannier_plot_radius', found, r_value=wannier_plot_radius)

    wannier_plot_scale = 1.0_dp
    call param_get_keyword('wannier_plot_scale', found, r_value=wannier_plot_scale)

    ! checks
    if (wannier_plot) then
      if ((index(wannier_plot_format, 'xcrys') .eq. 0) .and. (index(wannier_plot_format, 'cub') .eq. 0)) &
        call io_error('Error: wannier_plot_format not recognised')
      if ((index(wannier_plot_mode, 'crys') .eq. 0) .and. (index(wannier_plot_mode, 'mol') .eq. 0)) &
        call io_error('Error: wannier_plot_mode not recognised')
      if ((index(wannier_plot_spinor_mode, 'total') .eq. 0) .and. (index(wannier_plot_spinor_mode, 'up') .eq. 0) &
          .and. (index(wannier_plot_spinor_mode, 'down') .eq. 0)) &
        call io_error('Error: wannier_plot_spinor_mode not recognised')
      if (wannier_plot_radius < 0.0_dp) call io_error('Error: wannier_plot_radius must be positive')
      if (wannier_plot_scale < 0.0_dp) call io_error('Error: wannier_plot_scale must be positive')
    endif

    write_u_matrices = .false.
    call param_get_keyword('write_u_matrices', found, l_value=write_u_matrices)

    bands_plot = .false.
    call param_get_keyword('bands_plot', found, l_value=bands_plot)

    write_bvec = .false.
    call param_get_keyword('write_bvec', found, l_value=write_bvec)

    bands_num_points = 100
    call param_get_keyword('bands_num_points', found, i_value=bands_num_points)

    bands_plot_format = 'gnuplot'
    call param_get_keyword('bands_plot_format', found, c_value=bands_plot_format)

    bands_plot_mode = 's-k'
    call param_get_keyword('bands_plot_mode', found, c_value=bands_plot_mode)

    bands_plot_dim = 3
    call param_get_keyword('bands_plot_dim', found, i_value=bands_plot_dim)

    num_bands_project = 0
    call param_get_range_vector('bands_plot_project', found, num_bands_project, lcount=.true.)
    if (found) then
      if (num_bands_project < 1) call io_error('Error: problem reading bands_plot_project')
      if (allocated(bands_plot_project)) deallocate (bands_plot_project)
      allocate (bands_plot_project(num_bands_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating bands_plot_project in param_read')
      call param_get_range_vector('bands_plot_project', found, num_bands_project, .false., bands_plot_project)
      if (any(bands_plot_project < 1) .or. any(bands_plot_project > num_wann)) &
        call io_error('Error: bands_plot_project asks for a non-valid wannier function to be projected')
    endif

    bands_num_spec_points = 0
    call param_get_block_length('kpoint_path', found, i_temp)
    if (found) then
      bands_num_spec_points = i_temp*2
      if (allocated(bands_label)) deallocate (bands_label)
      allocate (bands_label(bands_num_spec_points), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating bands_label in param_read')
      if (allocated(bands_spec_points)) deallocate (bands_spec_points)
      allocate (bands_spec_points(3, bands_num_spec_points), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating bands_spec_points in param_read')
      call param_get_keyword_kpath
    end if
    if (.not. found .and. bands_plot) &
      call io_error('A bandstructure plot has been requested but there is no kpoint_path block')

    ! checks
    if (bands_plot) then
      if ((index(bands_plot_format, 'gnu') .eq. 0) .and. (index(bands_plot_format, 'xmgr') .eq. 0)) &
        call io_error('Error: bands_plot_format not recognised')
      if ((index(bands_plot_mode, 's-k') .eq. 0) .and. (index(bands_plot_mode, 'cut') .eq. 0)) &
        call io_error('Error: bands_plot_mode not recognised')
      if (bands_num_points < 0) call io_error('Error: bands_num_points must be positive')
    endif

    fermi_surface_plot = .false.
    call param_get_keyword('fermi_surface_plot', found, l_value=fermi_surface_plot)

    fermi_surface_num_points = 50
    call param_get_keyword('fermi_surface_num_points', found, i_value=fermi_surface_num_points)

    fermi_surface_plot_format = 'xcrysden'
    call param_get_keyword('fermi_surface_plot_format', &
                           found, c_value=fermi_surface_plot_format)

    nfermi = 0
    found_fermi_energy = .false.
    call param_get_keyword('fermi_energy', found, r_value=fermi_energy)
    if (found) then
      found_fermi_energy = .true.
      nfermi = 1
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
      nfermi = nint(abs((fermi_energy_max - fermi_energy_min)/fermi_energy_step)) + 1
    endif
    !
    if (found_fermi_energy) then
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(1), stat=ierr)
      fermi_energy_list(1) = fermi_energy
    elseif (fermi_energy_scan) then
      if (nfermi .eq. 1) then
        fermi_energy_step = 0.0_dp
      else
        fermi_energy_step = (fermi_energy_max - fermi_energy_min)/real(nfermi - 1, dp)
      endif
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(nfermi), stat=ierr)
      do i = 1, nfermi
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
      'Error allocating fermi_energy_list in param_read')

    ! checks
    if (fermi_surface_plot) then
      if ((index(fermi_surface_plot_format, 'xcrys') .eq. 0)) &
        call io_error('Error: fermi_surface_plot_format not recognised')
      if (fermi_surface_num_points < 0) &
        call io_error('Error: fermi_surface_num_points must be positive')
    endif

    kslice = .false.
    call param_get_keyword('kslice', found, l_value=kslice)

    kslice_task = 'fermi_lines'
    call param_get_keyword('kslice_task', found, c_value=kslice_task)
    if (kslice .and. index(kslice_task, 'fermi_lines') == 0 .and. &
        index(kslice_task, 'curv') == 0 .and. &
        index(kslice_task, 'morb') == 0 .and. &
        index(kslice_task, 'shc') == 0) call io_error &
      ('Error: value of kslice_task not recognised in param_read')
    if (kslice .and. index(kslice_task, 'curv') > 0 .and. &
        index(kslice_task, 'morb') > 0) call io_error &
      ("Error: kslice_task cannot include both 'curv' and 'morb'")
    if (kslice .and. index(kslice_task, 'shc') > 0 .and. &
        index(kslice_task, 'morb') > 0) call io_error &
      ("Error: kslice_task cannot include both 'shc' and 'morb'")
    if (kslice .and. index(kslice_task, 'shc') > 0 .and. &
        index(kslice_task, 'curv') > 0) call io_error &
      ("Error: kslice_task cannot include both 'shc' and 'curv'")

    kslice_2dkmesh(1:2) = 50
    call param_get_vector_length('kslice_2dkmesh', found, length=i)
    if (found) then
      if (i == 1) then
        call param_get_keyword_vector('kslice_2dkmesh', found, 1, &
                                      i_value=kslice_2dkmesh)
        kslice_2dkmesh(2) = kslice_2dkmesh(1)
      elseif (i == 2) then
        call param_get_keyword_vector('kslice_2dkmesh', found, 2, &
                                      i_value=kslice_2dkmesh)
      else
        call io_error('Error: kslice_2dkmesh must be provided as either' &
                      //' one integer or a vector of two integers')
      endif
      if (any(kslice_2dkmesh <= 0)) &
        call io_error('Error: kslice_2dkmesh elements must be' &
                      //' greater than zero')
    endif

    kslice_corner = 0.0_dp
    call param_get_keyword_vector('kslice_corner', found, 3, r_value=kslice_corner)

    kslice_b1(1) = 1.0_dp
    kslice_b1(2) = 0.0_dp
    kslice_b1(3) = 0.0_dp
    call param_get_keyword_vector('kslice_b1', found, 3, r_value=kslice_b1)

    kslice_b2(1) = 0.0_dp
    kslice_b2(2) = 1.0_dp
    kslice_b2(3) = 0.0_dp
    call param_get_keyword_vector('kslice_b2', found, 3, r_value=kslice_b2)

    kslice_fermi_lines_colour = 'none'
    call param_get_keyword('kslice_fermi_lines_colour', found, &
                           c_value=kslice_fermi_lines_colour)
    if (kslice .and. index(kslice_fermi_lines_colour, 'none') == 0 .and. &
        index(kslice_fermi_lines_colour, 'spin') == 0) call io_error &
      ('Error: value of kslice_fermi_lines_colour not recognised ' &
       //'in param_read')

!    slice_plot_format         = 'plotmv'
!    call param_get_keyword('slice_plot_format',found,c_value=slice_plot_format)

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

    !IVO

    dos = .false.
    call param_get_keyword('dos', found, l_value=dos)

    berry = .false.
    call param_get_keyword('berry', found, l_value=berry)

    transl_inv = .false.
    call param_get_keyword('transl_inv', found, l_value=transl_inv)

    berry_task = ' '
    call param_get_keyword('berry_task', found, c_value=berry_task)
    if (berry .and. .not. found) call io_error &
      ('Error: berry=T and berry_task is not set')
    if (berry .and. index(berry_task, 'ahc') == 0 .and. index(berry_task, 'morb') == 0 &
        .and. index(berry_task, 'kubo') == 0 .and. index(berry_task, 'sc') == 0 &
        .and. index(berry_task, 'shc') == 0) call io_error &
      ('Error: value of berry_task not recognised in param_read')

    ! Stepan
    gyrotropic = .false.
    call param_get_keyword('gyrotropic', found, l_value=gyrotropic)
    gyrotropic_task = 'all'
    call param_get_keyword('gyrotropic_task', found, c_value=gyrotropic_task)
    gyrotropic_box(:, :) = 0.0
    gyrotropic_degen_thresh = 0.0_dp
    call param_get_keyword('gyrotropic_degen_thresh', found, r_value=gyrotropic_degen_thresh)

    do i = 1, 3
      gyrotropic_box(i, i) = 1.0_dp
      gyrotropic_box_tmp(:) = 0.0_dp
      call param_get_keyword_vector('gyrotropic_box_b'//achar(48 + i), found, 3, r_value=gyrotropic_box_tmp)
      if (found) gyrotropic_box(i, :) = gyrotropic_box_tmp(:)
    enddo
    gyrotropic_box_corner(:) = 0.0_dp
    call param_get_keyword_vector('gyrotropic_box_center', found, 3, r_value=gyrotropic_box_tmp)
    if (found) gyrotropic_box_corner(:) = &
      gyrotropic_box_tmp(:) - 0.5*(gyrotropic_box(1, :) + gyrotropic_box(2, :) + gyrotropic_box(3, :))

    call param_get_range_vector('gyrotropic_band_list', found, gyrotropic_num_bands, lcount=.true.)
    if (found) then
      if (gyrotropic_num_bands < 1) call io_error('Error: problem reading gyrotropic_band_list')
      if (allocated(gyrotropic_band_list)) deallocate (gyrotropic_band_list)
      allocate (gyrotropic_band_list(gyrotropic_num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating gyrotropic_band_list in param_read')
      call param_get_range_vector('gyrotropic_band_list', found, gyrotropic_num_bands, .false., gyrotropic_band_list)
      if (any(gyrotropic_band_list < 1) .or. any(gyrotropic_band_list > num_wann)) &
        call io_error('Error: gyrotropic_band_list asks for a non-valid bands')
    else
      ! include all bands in the calculation
      gyrotropic_num_bands = num_wann
      if (allocated(gyrotropic_band_list)) deallocate (gyrotropic_band_list)
      allocate (gyrotropic_band_list(gyrotropic_num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating gyrotropic_band_list in param_read')
      do loop = 1, num_wann
        gyrotropic_band_list(loop) = loop
      end do
    end if

    smr_max_arg = 5.0
    call param_get_keyword('smr_max_arg', found, r_value=smr_max_arg)
    if (found .and. (smr_max_arg <= 0._dp)) &
      call io_error('Error: smr_max_arg must be greater than zero')

    gyrotropic_smr_max_arg = smr_max_arg
    call param_get_keyword('gyrotropic_smr_max_arg', found, &
                           r_value=gyrotropic_smr_max_arg)
    if (found .and. (gyrotropic_smr_max_arg <= 0._dp)) call io_error &
      ('Error: gyrotropic_smr_max_arg must be greater than zero')

!-------------------------------------------------------
!    alpha=0
!    call param_get_keyword('alpha',found,i_value=alpha)

!    beta=0
!    call param_get_keyword('beta',found,i_value=beta)

!    gamma=0
!    call param_get_keyword('gamma',found,i_value=gamma)
!-------------------------------------------------------

    berry_curv_adpt_kmesh = 1
    call param_get_keyword('berry_curv_adpt_kmesh', found, &
                           i_value=berry_curv_adpt_kmesh)
    if (berry_curv_adpt_kmesh < 1) &
      call io_error( &
      'Error:  berry_curv_adpt_kmesh must be a positive integer')

    berry_curv_adpt_kmesh_thresh = 100.0_dp
    call param_get_keyword('berry_curv_adpt_kmesh_thresh', found, &
                           r_value=berry_curv_adpt_kmesh_thresh)

    berry_curv_unit = 'ang2'
    call param_get_keyword('berry_curv_unit', found, c_value=berry_curv_unit)
    if (berry_curv_unit .ne. 'ang2' .and. berry_curv_unit .ne. 'bohr2') &
      call io_error &
      ('Error: value of berry_curv_unit not recognised in param_read')

    wanint_kpoint_file = .false.
    call param_get_keyword('wanint_kpoint_file', found, &
                           l_value=wanint_kpoint_file)

!    smear_temp = -1.0_dp
!    call param_get_keyword('smear_temp',found,r_value=smear_temp)

    kubo_adpt_smr = adpt_smr
    call param_get_keyword('kubo_adpt_smr', found, l_value=kubo_adpt_smr)

    kubo_adpt_smr_fac = adpt_smr_fac
    call param_get_keyword('kubo_adpt_smr_fac', found, &
                           r_value=kubo_adpt_smr_fac)
    if (found .and. (kubo_adpt_smr_fac <= 0._dp)) call io_error &
      ('Error: kubo_adpt_smr_fac must be greater than zero')

    kubo_adpt_smr_max = adpt_smr_max
    call param_get_keyword('kubo_adpt_smr_max', found, &
                           r_value=kubo_adpt_smr_max)
    if (kubo_adpt_smr_max <= 0._dp) call io_error &
      ('Error: kubo_adpt_smr_max must be greater than zero')

    kubo_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('kubo_smr_fixed_en_width', found, &
                           r_value=kubo_smr_fixed_en_width)
    if (found .and. (kubo_smr_fixed_en_width < 0._dp)) call io_error &
      ('Error: kubo_smr_fixed_en_width must be greater than or equal to zero')

    gyrotropic_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('gyrotropic_smr_fixed_en_width', found, &
                           r_value=gyrotropic_smr_fixed_en_width)
    if (found .and. (gyrotropic_smr_fixed_en_width < 0._dp)) call io_error &
      ('Error: gyrotropic_smr_fixed_en_width must be greater than or equal to zero')

    sc_phase_conv = 1
    call param_get_keyword('sc_phase_conv', found, i_value=sc_phase_conv)
    if ((sc_phase_conv .ne. 1) .and. ((sc_phase_conv .ne. 2))) call io_error('Error: sc_phase_conv must be either 1 or 2')

    scissors_shift = 0.0_dp
    call param_get_keyword('scissors_shift', found, &
                           r_value=scissors_shift)

    shc_freq_scan = .false.
    call param_get_keyword('shc_freq_scan', found, l_value=shc_freq_scan)

    shc_alpha = 1
    call param_get_keyword('shc_alpha', found, i_value=shc_alpha)
    if (found .and. (shc_alpha < 1 .or. shc_alpha > 3)) call io_error &
      ('Error:  shc_alpha must be 1, 2 or 3')

    shc_beta = 2
    call param_get_keyword('shc_beta', found, i_value=shc_beta)
    if (found .and. (shc_beta < 1 .or. shc_beta > 3)) call io_error &
      ('Error:  shc_beta must be 1, 2 or 3')

    shc_gamma = 3
    call param_get_keyword('shc_gamma', found, i_value=shc_gamma)
    if (found .and. (shc_gamma < 1 .or. shc_gamma > 3)) call io_error &
      ('Error:  shc_gamma must be 1, 2 or 3')

    shc_bandshift = .false.
    call param_get_keyword('shc_bandshift', found, l_value=shc_bandshift)
    shc_bandshift = shc_bandshift .and. berry .and. .not. (index(berry_task, 'shc') == 0)
    if ((abs(scissors_shift) > 1.0e-7_dp) .and. shc_bandshift) &
      call io_error('Error: shc_bandshift and scissors_shift cannot be used simultaneously')

    shc_bandshift_firstband = 0
    call param_get_keyword('shc_bandshift_firstband', found, i_value=shc_bandshift_firstband)
    if (shc_bandshift .and. (.not. found)) &
      call io_error('Error: shc_bandshift required but no shc_bandshift_firstband provided')
    if ((shc_bandshift_firstband < 1) .and. found) &
      call io_error('Error: shc_bandshift_firstband must >= 1')

    shc_bandshift_energyshift = 0._dp
    call param_get_keyword('shc_bandshift_energyshift', found, r_value=shc_bandshift_energyshift)
    if (shc_bandshift .and. (.not. found)) &
      call io_error('Error: shc_bandshift required but no shc_bandshift_energyshift provided')

    spin_moment = .false.
    call param_get_keyword('spin_moment', found, &
                           l_value=spin_moment)

    spin_axis_polar = 0.0_dp
    call param_get_keyword('spin_axis_polar', found, &
                           r_value=spin_axis_polar)

    spin_axis_azimuth = 0.0_dp
    call param_get_keyword('spin_axis_azimuth', found, &
                           r_value=spin_axis_azimuth)

    spin_decomp = .false.
    call param_get_keyword('spin_decomp', found, l_value=spin_decomp)

    if (spin_decomp .and. (num_elec_per_state .ne. 1)) then
      call io_error('spin_decomp can be true only if num_elec_per_state is 1')
    end if

    use_degen_pert = .false.
    call param_get_keyword('use_degen_pert', found, &
                           l_value=use_degen_pert)

    degen_thr = 1.0d-4
    call param_get_keyword('degen_thr', found, r_value=degen_thr)

    kpath = .false.
    call param_get_keyword('kpath', found, l_value=kpath)

    kpath_task = 'bands'
    call param_get_keyword('kpath_task', found, c_value=kpath_task)
    if (kpath .and. index(kpath_task, 'bands') == 0 .and. &
        index(kpath_task, 'curv') == 0 .and. &
        index(kpath_task, 'morb') == 0 .and. &
        index(kpath_task, 'shc') == 0) call io_error &
      ('Error: value of kpath_task not recognised in param_read')
    if (bands_num_spec_points == 0 .and. kpath) &
      call io_error('Error: a kpath plot has been requested but there is no kpoint_path block')

    kpath_num_points = 100
    call param_get_keyword('kpath_num_points', found, &
                           i_value=kpath_num_points)
    if (kpath_num_points < 0) &
      call io_error('Error: kpath_num_points must be positive')

    kpath_bands_colour = 'none'
    call param_get_keyword('kpath_bands_colour', found, &
                           c_value=kpath_bands_colour)
    if (kpath .and. index(kpath_bands_colour, 'none') == 0 .and. &
        index(kpath_bands_colour, 'spin') == 0 .and. &
        index(kpath_bands_colour, 'shc') == 0) call io_error &
      ('Error: value of kpath_bands_colour not recognised in param_read')
    if (kpath .and. index(kpath_task, 'shc') > 0 .and. &
        index(kpath_task, 'spin') > 0) call io_error &
      ("Error: kpath_task cannot include both 'shc' and 'spin'")

    ! set to a negative default value
    num_valence_bands = -99
    call param_get_keyword('num_valence_bands', found, i_value=num_valence_bands)
    if (found .and. (num_valence_bands .le. 0)) &
      call io_error('Error: num_valence_bands should be greater than zero')
    ! there is a check on this parameter later

    dos_task = 'dos_plot'
    if (dos) then
      dos_plot = .true.
    else
      dos_plot = .false.
    endif
    call param_get_keyword('dos_task', found, c_value=dos_task)
    if (dos) then
      if (index(dos_task, 'dos_plot') == 0 .and. &
          index(dos_task, 'find_fermi_energy') == 0) call io_error &
        ('Error: value of dos_task not recognised in param_read')
      if (index(dos_task, 'dos_plot') > 0) dos_plot = .true.
      if (index(dos_task, 'find_fermi_energy') > 0 .and. found_fermi_energy) &
        call io_error &
        ('Error: Cannot set "dos_task = find_fermi_energy" and give a value to "fermi_energy"')
    end if

!    sigma_abc_onlyorb=.false.
!    call param_get_keyword('sigma_abc_onlyorb',found,l_value=sigma_abc_onlyorb)

! -------------------------------------------------------------------

    !IVO_END

    dos_energy_step = 0.01_dp
    call param_get_keyword('dos_energy_step', found, r_value=dos_energy_step)

    dos_adpt_smr = adpt_smr
    call param_get_keyword('dos_adpt_smr', found, l_value=dos_adpt_smr)

    dos_adpt_smr_fac = adpt_smr_fac
    call param_get_keyword('dos_adpt_smr_fac', found, r_value=dos_adpt_smr_fac)
    if (found .and. (dos_adpt_smr_fac <= 0._dp)) &
      call io_error('Error: dos_adpt_smr_fac must be greater than zero')

    dos_adpt_smr_max = adpt_smr_max
    call param_get_keyword('dos_adpt_smr_max', found, r_value=dos_adpt_smr_max)
    if (dos_adpt_smr_max <= 0._dp) call io_error &
      ('Error: dos_adpt_smr_max must be greater than zero')

    dos_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('dos_smr_fixed_en_width', found, r_value=dos_smr_fixed_en_width)
    if (found .and. (dos_smr_fixed_en_width < 0._dp)) &
      call io_error('Error: dos_smr_fixed_en_width must be greater than or equal to zero')

!    dos_gaussian_width        = 0.1_dp
!    call param_get_keyword('dos_gaussian_width',found,r_value=dos_gaussian_width)

!    dos_plot_format           = 'gnuplot'
!    call param_get_keyword('dos_plot_format',found,c_value=dos_plot_format)

    call param_get_range_vector('dos_project', found, num_dos_project, &
                                lcount=.true.)
    if (found) then
      if (num_dos_project < 1) call io_error('Error: problem reading dos_project')
      if (allocated(dos_project)) deallocate (dos_project)
      allocate (dos_project(num_dos_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dos_project in param_read')
      call param_get_range_vector('dos_project', found, num_dos_project, &
                                  .false., dos_project)
      if (any(dos_project < 1) .or. any(dos_project > num_wann)) call io_error &
        ('Error: dos_project asks for out-of-range Wannier functions')
    else
      ! by default plot all
      num_dos_project = num_wann
      if (allocated(dos_project)) deallocate (dos_project)
      allocate (dos_project(num_dos_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dos_project in param_read')
      do i = 1, num_dos_project
        dos_project(i) = i
      end do
    endif

    hr_plot = .false.
    call param_get_keyword('hr_plot', found, l_value=hr_plot)
    if (found) call io_error('Input parameter hr_plot is no longer used. Please use write_hr instead.')
    write_hr = .false.
    call param_get_keyword('write_hr', found, l_value=write_hr)

    write_rmn = .false.
    call param_get_keyword('write_rmn', found, l_value=write_rmn)

    write_tb = .false.
    call param_get_keyword('write_tb', found, l_value=write_tb)

    hr_cutoff = 0.0_dp
    call param_get_keyword('hr_cutoff', found, r_value=hr_cutoff)

    dist_cutoff_mode = 'three_dim'
    call param_get_keyword('dist_cutoff_mode', found, c_value=dist_cutoff_mode)
    if ((index(dist_cutoff_mode, 'three_dim') .eq. 0) &
        .and. (index(dist_cutoff_mode, 'two_dim') .eq. 0) &
        .and. (index(dist_cutoff_mode, 'one_dim') .eq. 0)) &
      call io_error('Error: dist_cutoff_mode not recognised')

! aam_2012-04-13: moved later
!    dist_cutoff                 = 1000.0_dp
!    call param_get_keyword('dist_cutoff',found,r_value=dist_cutoff)

    one_dim_axis = 'none'
    call param_get_keyword('one_dim_axis', found, c_value=one_dim_axis)
    one_dim_dir = 0
    if (index(one_dim_axis, 'x') > 0) one_dim_dir = 1
    if (index(one_dim_axis, 'y') > 0) one_dim_dir = 2
    if (index(one_dim_axis, 'z') > 0) one_dim_dir = 3
    if (transport .and. .not. tran_read_ht .and. (one_dim_dir .eq. 0)) call io_error('Error: one_dim_axis not recognised')
    if (bands_plot .and. (index(bands_plot_mode, 'cut') .ne. 0)&
       & .and. ((bands_plot_dim .ne. 3) .or. (index(dist_cutoff_mode, 'three_dim') .eq. 0))&
       & .and. (one_dim_dir .eq. 0)) &
         call io_error('Error: one_dim_axis not recognised')

301 continue

    use_ws_distance = .true.
    call param_get_keyword('use_ws_distance', found, l_value=use_ws_distance)

    ws_distance_tol = 1.e-5_dp
    call param_get_keyword('ws_distance_tol', found, r_value=ws_distance_tol)

    ws_search_size = 2

    call param_get_vector_length('ws_search_size', found, length=i)
    if (found) then
      if (i .eq. 1) then
        call param_get_keyword_vector('ws_search_size', found, 1, &
                                      i_value=ws_search_size)
        ws_search_size(2) = ws_search_size(1)
        ws_search_size(3) = ws_search_size(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector('ws_search_size', found, 3, &
                                      i_value=ws_search_size)
      else
        call io_error('Error: ws_search_size must be provided as either one integer or a vector of three integers')
      end if
      if (any(ws_search_size <= 0)) &
        call io_error('Error: ws_search_size elements must be greater than zero')
    end if

    !%%%%%%%%%%%%%%%%
    ! Transport
    !%%%%%%%%%%%%%%%%

    transport_mode = 'bulk'
    call param_get_keyword('transport_mode', found, c_value=transport_mode)

!    if ( .not.tran_read_ht  .and. (index(transport_mode,'lcr').ne.0) ) &
!       call io_error('Error: transport_mode.eq.lcr not compatible with tran_read_ht.eq.false')

    tran_win_min = -3.0_dp
    call param_get_keyword('tran_win_min', found, r_value=tran_win_min)

    tran_win_max = 3.0_dp
    call param_get_keyword('tran_win_max', found, r_value=tran_win_max)

    tran_energy_step = 0.01_dp
    call param_get_keyword('tran_energy_step', found, r_value=tran_energy_step)

    tran_num_bb = 0
    call param_get_keyword('tran_num_bb', found, i_value=tran_num_bb)

    tran_num_ll = 0
    call param_get_keyword('tran_num_ll', found, i_value=tran_num_ll)

    tran_num_rr = 0
    call param_get_keyword('tran_num_rr', found, i_value=tran_num_rr)

    tran_num_cc = 0
    call param_get_keyword('tran_num_cc', found, i_value=tran_num_cc)

    tran_num_lc = 0
    call param_get_keyword('tran_num_lc', found, i_value=tran_num_lc)

    tran_num_cr = 0
    call param_get_keyword('tran_num_cr', found, i_value=tran_num_cr)

    tran_num_bandc = 0
    call param_get_keyword('tran_num_bandc', found, i_value=tran_num_bandc)

    tran_write_ht = .false.
    call param_get_keyword('tran_write_ht', found, l_value=tran_write_ht)

    tran_use_same_lead = .true.
    call param_get_keyword('tran_use_same_lead', found, l_value=tran_use_same_lead)

    tran_num_cell_ll = 0
    call param_get_keyword('tran_num_cell_ll', found, i_value=tran_num_cell_ll)

    tran_num_cell_rr = 0
    call param_get_keyword('tran_num_cell_rr', found, i_value=tran_num_cell_rr)

    tran_group_threshold = 0.15_dp
    call param_get_keyword('tran_group_threshold', found, r_value=tran_group_threshold)

    dist_cutoff = 1000.0_dp
    call param_get_keyword('dist_cutoff', found, r_value=dist_cutoff)

    dist_cutoff_hc = dist_cutoff
    call param_get_keyword('dist_cutoff_hc', found, r_value=dist_cutoff_hc)

    ! checks
    if (transport) then
      if ((index(transport_mode, 'bulk') .eq. 0) .and. (index(transport_mode, 'lcr') .eq. 0)) &
        call io_error('Error: transport_mode not recognised')
      if (tran_num_bb < 0) call io_error('Error: tran_num_bb < 0')
      if (tran_num_ll < 0) call io_error('Error: tran_num_ll < 0')
      if (tran_num_rr < 0) call io_error('Error: tran_num_rr < 0')
      if (tran_num_cc < 0) call io_error('Error: tran_num_cc < 0')
      if (tran_num_lc < 0) call io_error('Error: tran_num_lc < 0')
      if (tran_num_cr < 0) call io_error('Error: tran_num_cr < 0')
      if (tran_num_bandc < 0) call io_error('Error: tran_num_bandc < 0')
      if (tran_num_cell_ll < 0) call io_error('Error: tran_num_cell_ll < 0')
      if (tran_num_cell_rr < 0) call io_error('Error: tran_num_cell_rr < 0')
      if (tran_group_threshold < 0.0_dp) call io_error('Error: tran_group_threshold < 0')
    endif

    if (transport .and. tran_read_ht) goto 302

    !%%%%%%%%%%%%%%%%
    ! Disentanglement
    !%%%%%%%%%%%%%%%%

    disentanglement = .false.
    if (num_bands > num_wann) disentanglement = .true.

    ! These must be read here, before the check on the existence of the .eig file!
    geninterp = .false.
    call param_get_keyword('geninterp', found, l_value=geninterp)
    boltzwann = .false.
    call param_get_keyword('boltzwann', found, l_value=boltzwann)

    ! Read the eigenvalues from wannier.eig
    eig_found = .false.
    if (.not. library .and. .not. effective_model) then

      if (.not. postproc_setup) then
        inquire (file=trim(seedname)//'.eig', exist=eig_found)
        if (.not. eig_found) then
          if (disentanglement) then
            call io_error('No '//trim(seedname)//'.eig file found. Needed for disentanglement')
          else if ((bands_plot .or. dos_plot .or. fermi_surface_plot .or. write_hr .or. boltzwann &
                    .or. geninterp)) then
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

    dis_win_min = -1.0_dp; dis_win_max = 0.0_dp
    if (eig_found) dis_win_min = minval(eigval)
    call param_get_keyword('dis_win_min', found, r_value=dis_win_min)

    if (eig_found) dis_win_max = maxval(eigval)
    call param_get_keyword('dis_win_max', found, r_value=dis_win_max)
    if (eig_found .and. (dis_win_max .lt. dis_win_min)) &
      call io_error('Error: param_read: check disentanglement windows')

    dis_froz_min = -1.0_dp; dis_froz_max = 0.0_dp
    ! no default for dis_froz_max
    frozen_states = .false.
    call param_get_keyword('dis_froz_max', found, r_value=dis_froz_max)
    if (found) then
      frozen_states = .true.
      dis_froz_min = dis_win_min ! default value for the bottom of frozen window
    end if
    call param_get_keyword('dis_froz_min', found2, r_value=dis_froz_min)
    if (eig_found) then
      if (dis_froz_max .lt. dis_froz_min) &
        call io_error('Error: param_read: check disentanglement frozen windows')
      if (found2 .and. .not. found) &
        call io_error('Error: param_read: found dis_froz_min but not dis_froz_max')
    endif

    dis_num_iter = 200
    call param_get_keyword('dis_num_iter', found, i_value=dis_num_iter)
    if (dis_num_iter < 0) call io_error('Error: dis_num_iter must be positive')

    dis_mix_ratio = 0.5_dp
    call param_get_keyword('dis_mix_ratio', found, r_value=dis_mix_ratio)
    if (dis_mix_ratio <= 0.0_dp .or. dis_mix_ratio > 1.0_dp) &
      call io_error('Error: dis_mix_ratio must be greater than 0.0 but not greater than 1.0')

    dis_conv_tol = 1.0e-10_dp
    call param_get_keyword('dis_conv_tol', found, r_value=dis_conv_tol)
    if (dis_conv_tol < 0.0_dp) call io_error('Error: dis_conv_tol must be positive')

    dis_conv_window = 3
    call param_get_keyword('dis_conv_window', found, i_value=dis_conv_window)
    if (dis_conv_window < 0) call io_error('Error: dis_conv_window must be positive')

    ! GS-start
    dis_spheres_first_wann = 1
    call param_get_keyword('dis_spheres_first_wann', found, i_value=dis_spheres_first_wann)
    if (dis_spheres_first_wann < 1) call io_error('Error: dis_spheres_first_wann must be greater than 0')
    if (dis_spheres_first_wann > num_bands - num_wann + 1) &
      call io_error('Error: dis_spheres_first_wann is larger than num_bands-num_wann+1')
    dis_spheres_num = 0
    call param_get_keyword('dis_spheres_num', found, i_value=dis_spheres_num)
    if (dis_spheres_num < 0) call io_error('Error: dis_spheres_num cannot be negative')
    if (dis_spheres_num > 0) then
      allocate (dis_spheres(4, dis_spheres_num), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dis_spheres in param_read')
      call param_get_keyword_block('dis_spheres', found, dis_spheres_num, 4, r_value=dis_spheres)
      if (.not. found) call io_error('Error: Did not find dis_spheres in the input file')
      do nkp = 1, dis_spheres_num
        if (dis_spheres(4, nkp) < 1.0e-15_dp) &
          call io_error('Error: radius for dis_spheres must be > 0')
      enddo
    endif
    ! GS-end

    ! [gp-begin, Jun 1, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! General band interpolator (geninterp)
    !%%%%%%%%%%%%%%%%%%%%
    geninterp_alsofirstder = .false.
    call param_get_keyword('geninterp_alsofirstder', found, l_value=geninterp_alsofirstder)
    geninterp_single_file = .true.
    call param_get_keyword('geninterp_single_file', found, l_value=geninterp_single_file)
    ! [gp-end, Jun 1, 2012]

    ! [gp-begin, Apr 12, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! Boltzmann transport
    !%%%%%%%%%%%%%%%%%%%%
    ! Note: to be put AFTER the disentanglement routines!

    boltz_calc_also_dos = .false.
    call param_get_keyword('boltz_calc_also_dos', found, l_value=boltz_calc_also_dos)

    boltz_calc_also_dos = boltz_calc_also_dos .and. boltzwann

    ! 0 means the normal 3d case for the calculation of the Seebeck coefficient
    ! The other valid possibilities are 1,2,3 for x,y,z respectively
    boltz_2d_dir_num = 0
    call param_get_keyword('boltz_2d_dir', found, c_value=boltz_2d_dir)
    if (found) then
      if (trim(boltz_2d_dir) == 'no') then
        boltz_2d_dir_num = 0
      elseif (trim(boltz_2d_dir) == 'x') then
        boltz_2d_dir_num = 1
      elseif (trim(boltz_2d_dir) == 'y') then
        boltz_2d_dir_num = 2
      elseif (trim(boltz_2d_dir) == 'z') then
        boltz_2d_dir_num = 3
      else
        call io_error('Error: boltz_2d_dir can only be "no", "x", "y" or "z".')
      end if
    end if

    boltz_dos_energy_step = 0.001_dp
    call param_get_keyword('boltz_dos_energy_step', found, r_value=boltz_dos_energy_step)
    if (found .and. (boltz_dos_energy_step <= 0._dp)) &
      call io_error('Error: boltz_dos_energy_step must be positive')

    if (allocated(eigval)) then
      boltz_dos_energy_min = minval(eigval) - 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      boltz_dos_energy_min = -1.0_dp
    end if
    call param_get_keyword('boltz_dos_energy_min', found, r_value=boltz_dos_energy_min)
    if (allocated(eigval)) then
      boltz_dos_energy_max = maxval(eigval) + 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      boltz_dos_energy_max = 0.0_dp
    end if
    call param_get_keyword('boltz_dos_energy_max', found, r_value=boltz_dos_energy_max)
    if (boltz_dos_energy_max <= boltz_dos_energy_min) &
      call io_error('Error: boltz_dos_energy_max must be greater than boltz_dos_energy_min')

    boltz_dos_adpt_smr = adpt_smr
    call param_get_keyword('boltz_dos_adpt_smr', found, l_value=boltz_dos_adpt_smr)

    boltz_dos_adpt_smr_fac = adpt_smr_fac
    call param_get_keyword('boltz_dos_adpt_smr_fac', found, r_value=boltz_dos_adpt_smr_fac)
    if (found .and. (boltz_dos_adpt_smr_fac <= 0._dp)) &
      call io_error('Error: boltz_dos_adpt_smr_fac must be greater than zero')

    boltz_dos_adpt_smr_max = adpt_smr_max
    call param_get_keyword('boltz_dos_adpt_smr_max', found, r_value=boltz_dos_adpt_smr_max)
    if (boltz_dos_adpt_smr_max <= 0._dp) call io_error &
      ('Error: boltz_dos_adpt_smr_max must be greater than zero')

    boltz_dos_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('boltz_dos_smr_fixed_en_width', found, r_value=boltz_dos_smr_fixed_en_width)
    if (found .and. (boltz_dos_smr_fixed_en_width < 0._dp)) &
      call io_error('Error: boltz_dos_smr_fixed_en_width must be greater than or equal to zero')

    boltz_mu_min = -999._dp
    call param_get_keyword('boltz_mu_min', found, r_value=boltz_mu_min)
    if ((.not. found) .and. boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_min provided')
    boltz_mu_max = -999._dp
    call param_get_keyword('boltz_mu_max', found2, r_value=boltz_mu_max)
    if ((.not. found2) .and. boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_max provided')
    if (found .and. found2 .and. (boltz_mu_max < boltz_mu_min)) &
      call io_error('Error: boltz_mu_max must be greater than boltz_mu_min')
    boltz_mu_step = 0._dp
    call param_get_keyword('boltz_mu_step', found, r_value=boltz_mu_step)
    if ((.not. found) .and. boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_step provided')
    if (found .and. (boltz_mu_step <= 0._dp)) &
      call io_error('Error: boltz_mu_step must be greater than zero')

    boltz_temp_min = -999._dp
    call param_get_keyword('boltz_temp_min', found, r_value=boltz_temp_min)
    if ((.not. found) .and. boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_min provided')
    boltz_temp_max = -999._dp
    call param_get_keyword('boltz_temp_max', found2, r_value=boltz_temp_max)
    if ((.not. found2) .and. boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_max provided')
    if (found .and. found2 .and. (boltz_temp_max < boltz_temp_min)) &
      call io_error('Error: boltz_temp_max must be greater than boltz_temp_min')
    if (found .and. (boltz_temp_min <= 0._dp)) &
      call io_error('Error: boltz_temp_min must be greater than zero')
    boltz_temp_step = 0._dp
    call param_get_keyword('boltz_temp_step', found, r_value=boltz_temp_step)
    if ((.not. found) .and. boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_step provided')
    if (found .and. (boltz_temp_step <= 0._dp)) &
      call io_error('Error: boltz_temp_step must be greater than zero')

    ! The interpolation mesh is read later on

    ! By default, the energy step for the TDF is 1 meV
    boltz_tdf_energy_step = 0.001_dp
    call param_get_keyword('boltz_tdf_energy_step', found, r_value=boltz_tdf_energy_step)
    if (boltz_tdf_energy_step <= 0._dp) &
      call io_error('Error: boltz_tdf_energy_step must be greater than zero')

    ! For TDF: TDF smeared in a NON-adaptive way; value in eV, default = 0._dp
    ! (i.e., no smearing)
    boltz_TDF_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('boltz_tdf_smr_fixed_en_width', found, r_value=boltz_TDF_smr_fixed_en_width)
    if (found .and. (boltz_TDF_smr_fixed_en_width < 0._dp)) &
      call io_error('Error: boltz_TDF_smr_fixed_en_width must be greater than or equal to zero')

    ! By default: use the "global" smearing index
    boltz_TDF_smr_index = smr_index
    call param_get_keyword('boltz_tdf_smr_type', found, c_value=ctmp)
    if (found) boltz_TDF_smr_index = get_smearing_index(ctmp, 'boltz_tdf_smr_type')

    ! By default: use the "global" smearing index
    boltz_dos_smr_index = smr_index
    call param_get_keyword('boltz_dos_smr_type', found, c_value=ctmp)
    if (found) boltz_dos_smr_index = get_smearing_index(ctmp, 'boltz_dos_smr_type')

    ! By default: use the "global" smearing index
    dos_smr_index = smr_index
    call param_get_keyword('dos_smr_type', found, c_value=ctmp)
    if (found) dos_smr_index = get_smearing_index(ctmp, 'dos_smr_type')

    ! By default: use the "global" smearing index
    kubo_smr_index = smr_index
    call param_get_keyword('kubo_smr_type', found, c_value=ctmp)
    if (found) kubo_smr_index = get_smearing_index(ctmp, 'kubo_smr_type')

    ! By default: use the "global" smearing index
    gyrotropic_smr_index = smr_index
    call param_get_keyword('gyrotropic_smr_type', found, c_value=ctmp)
    if (found) gyrotropic_smr_index = get_smearing_index(ctmp, 'gyrotropic_smr_type')

    ! By default: 10 fs relaxation time
    boltz_relax_time = 10._dp
    call param_get_keyword('boltz_relax_time', found, r_value=boltz_relax_time)

    boltz_bandshift = .false.
    call param_get_keyword('boltz_bandshift', found, l_value=boltz_bandshift)
    boltz_bandshift = boltz_bandshift .and. boltzwann

    boltz_bandshift_firstband = 0
    call param_get_keyword('boltz_bandshift_firstband', found, i_value=boltz_bandshift_firstband)
    if (boltz_bandshift .and. (.not. found)) &
      call io_error('Error: boltz_bandshift required but no boltz_bandshift_firstband provided')
    boltz_bandshift_energyshift = 0._dp
    call param_get_keyword('boltz_bandshift_energyshift', found, r_value=boltz_bandshift_energyshift)
    if (boltz_bandshift .and. (.not. found)) &
      call io_error('Error: boltz_bandshift required but no boltz_bandshift_energyshift provided')
    ! [gp-end, Apr 12, 2012]

    !%%%%%%%%%%%%%%%%
    !  Other Stuff
    !%%%%%%%%%%%%%%%%

    ! aam: vdW
    write_vdw_data = .false.
    call param_get_keyword('write_vdw_data', found, l_value=write_vdw_data)
    if (write_vdw_data) then
      if ((.not. gamma_only) .or. (num_kpts .ne. 1)) &
        call io_error('Error: write_vdw_data may only be used with a single k-point at Gamma')
    endif
    if (write_vdw_data .and. disentanglement .and. num_valence_bands .le. 0) &
      call io_error('If writing vdw data and disentangling then num_valence_bands must be defined')

    if (frozen_states) then
      dos_energy_max = dis_froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      dos_energy_max = maxval(eigval) + 0.6667_dp
    else
      dos_energy_max = dis_win_max + 0.6667_dp
    end if
    call param_get_keyword('dos_energy_max', found, r_value=dos_energy_max)

    if (allocated(eigval)) then
      dos_energy_min = minval(eigval) - 0.6667_dp
    else
      dos_energy_min = dis_win_min - 0.6667_dp
    end if
    call param_get_keyword('dos_energy_min', found, r_value=dos_energy_min)

    kubo_freq_min = 0.0_dp
    gyrotropic_freq_min = kubo_freq_min
    call param_get_keyword('kubo_freq_min', found, r_value=kubo_freq_min)
    !
    if (frozen_states) then
      kubo_freq_max = dis_froz_max - fermi_energy_list(1) + 0.6667_dp
    elseif (allocated(eigval)) then
      kubo_freq_max = maxval(eigval) - minval(eigval) + 0.6667_dp
    else
      kubo_freq_max = dis_win_max - dis_win_min + 0.6667_dp
    end if
    gyrotropic_freq_max = kubo_freq_max
    call param_get_keyword('kubo_freq_max', found, r_value=kubo_freq_max)

    !
    kubo_freq_step = 0.01_dp
    call param_get_keyword('kubo_freq_step', found, r_value=kubo_freq_step)
    if (found .and. kubo_freq_step < 0.0_dp) call io_error( &
      'Error: kubo_freq_step must be positive')
    !
    kubo_nfreq = nint((kubo_freq_max - kubo_freq_min)/kubo_freq_step) + 1
    if (kubo_nfreq <= 1) kubo_nfreq = 2
    kubo_freq_step = (kubo_freq_max - kubo_freq_min)/(kubo_nfreq - 1)
    !
    if (allocated(kubo_freq_list)) deallocate (kubo_freq_list)
    allocate (kubo_freq_list(kubo_nfreq), stat=ierr)
    if (ierr /= 0) &
      call io_error('Error allocating kubo_freq_list in param_read')
    do i = 1, kubo_nfreq
      kubo_freq_list(i) = kubo_freq_min &
                          + (i - 1)*(kubo_freq_max - kubo_freq_min)/(kubo_nfreq - 1)
    enddo
    !
    ! TODO: Alternatively, read list of (complex) frequencies; kubo_nfreq is
    !       the length of the list

    gyrotropic_freq_step = 0.01_dp
    call param_get_keyword('gyrotropic_freq_min', found, r_value=gyrotropic_freq_min)
    call param_get_keyword('gyrotropic_freq_max', found, r_value=gyrotropic_freq_max)
    call param_get_keyword('gyrotropic_freq_step', found, r_value=gyrotropic_freq_step)
    gyrotropic_nfreq = nint((gyrotropic_freq_max - gyrotropic_freq_min)/gyrotropic_freq_step) + 1
    if (gyrotropic_nfreq <= 1) gyrotropic_nfreq = 2
    gyrotropic_freq_step = (gyrotropic_freq_max - gyrotropic_freq_min)/(gyrotropic_nfreq - 1)
    if (allocated(gyrotropic_freq_list)) deallocate (gyrotropic_freq_list)
    allocate (gyrotropic_freq_list(gyrotropic_nfreq), stat=ierr)
    if (ierr /= 0) &
      call io_error('Error allocating gyrotropic_freq_list in param_read')
    do i = 1, gyrotropic_nfreq
      gyrotropic_freq_list(i) = gyrotropic_freq_min &
                                + (i - 1)*(gyrotropic_freq_max - gyrotropic_freq_min)/(gyrotropic_nfreq - 1) &
                                + cmplx_i*gyrotropic_smr_fixed_en_width
    enddo

    if (frozen_states) then
      kubo_eigval_max = dis_froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      kubo_eigval_max = maxval(eigval) + 0.6667_dp
    else
      kubo_eigval_max = dis_win_max + 0.6667_dp
    end if
    gyrotropic_eigval_max = kubo_eigval_max

    call param_get_keyword('kubo_eigval_max', found, r_value=kubo_eigval_max)
    call param_get_keyword('gyrotropic_eigval_max', found, r_value=gyrotropic_eigval_max)

    automatic_translation = .true.
    translation_centre_frac = 0.0_dp
    call param_get_keyword_vector('translation_centre_frac', found, 3, r_value=rv_temp)
    if (found) then
      translation_centre_frac = rv_temp
      automatic_translation = .false.
    endif

    sc_eta = 0.04
    call param_get_keyword('sc_eta', found, r_value=sc_eta)

    sc_w_thr = 5.0d0
    call param_get_keyword('sc_w_thr', found, r_value=sc_w_thr)

    use_bloch_phases = .false.
    call param_get_keyword('use_bloch_phases', found, l_value=use_bloch_phases)
    if (disentanglement .and. use_bloch_phases) &
      call io_error('Error: Cannot use bloch phases for disentanglement')

    search_shells = 36
    call param_get_keyword('search_shells', found, i_value=search_shells)
    if (search_shells < 0) call io_error('Error: search_shells must be positive')

    kmesh_tol = 0.000001_dp
    call param_get_keyword('kmesh_tol', found, r_value=kmesh_tol)
    if (kmesh_tol < 0.0_dp) call io_error('Error: kmesh_tol must be positive')

    num_shells = 0
    call param_get_range_vector('shell_list', found, num_shells, lcount=.true.)
    if (found) then
      if (num_shells < 0 .or. num_shells > max_shells) &
        call io_error('Error: number of shell in shell_list must be between zero and six')
      if (allocated(shell_list)) deallocate (shell_list)
      allocate (shell_list(num_shells), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating shell_list in param_read')
      call param_get_range_vector('shell_list', found, num_shells, .false., shell_list)
      if (any(shell_list < 1)) &
        call io_error('Error: shell_list must contain positive numbers')
    else
      if (allocated(shell_list)) deallocate (shell_list)
      allocate (shell_list(max_shells), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating shell_list in param_read')
    end if

    call param_get_keyword('num_shells', found, i_value=itmp)
    if (found .and. (itmp /= num_shells)) &
      call io_error('Error: Found obsolete keyword num_shells. Its value does not agree with shell_list')

    ! If .true., does not perform the check of B1 of
    ! Marzari, Vanderbild, PRB 56, 12847 (1997)
    ! in kmesh.F90
    ! mainly needed for the interaction with Z2PACK
    ! By default: .false. (perform the tests)
    skip_B1_tests = .false.
    call param_get_keyword('skip_b1_tests', found, l_value=skip_B1_tests)

    call param_get_keyword_block('unit_cell_cart', found, 3, 3, r_value=real_lattice_tmp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <unit_cell_cart> in input file'
    if (.not. library) then
      real_lattice = transpose(real_lattice_tmp)
      if (.not. found) call io_error('Error: Did not find the cell information in the input file')
    end if

    if (.not. library) &
      call utility_recip_lattice(real_lattice, recip_lattice, cell_volume)
    call utility_metric(real_lattice, recip_lattice, real_metric, recip_metric)

    if (.not. effective_model) allocate (kpt_cart(3, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating kpt_cart in param_read')
    if (.not. library .and. .not. effective_model) then
      allocate (kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpt_latt in param_read')
    end if

    call param_get_keyword_block('kpoints', found, num_kpts, 3, r_value=kpt_cart)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <kpoints> in input file'
    if (.not. library .and. .not. effective_model) then
      kpt_latt = kpt_cart
      if (.not. found) call io_error('Error: Did not find the kpoint information in the input file')
    end if

    ! Calculate the kpoints in cartesian coordinates
    if (.not. effective_model) then
      do nkp = 1, num_kpts
        kpt_cart(:, nkp) = matmul(kpt_latt(:, nkp), recip_lattice(:, :))
      end do
    endif

    ! get the nnkpts block -- this is allowed only in postproc-setup mode
    call param_get_block_length('nnkpts', explicit_nnkpts, rows)
    if (explicit_nnkpts) then
      nntot = rows/num_kpts
      if (modulo(rows, num_kpts) /= 0) then
        call io_error('The number of rows in nnkpts must be a multiple of num_kpts')
      end if
      if (allocated(nnkpts_block)) deallocate (nnkpts_block)
      allocate (nnkpts_block(5, rows), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nnkpts_block in param_read')
      call param_get_keyword_block('nnkpts', found, rows, 5, i_value=nnkpts_block)
      ! check that postproc_setup is true
      if (.not. postproc_setup) &
        call io_error('Input parameter nnkpts_block is allowed only if postproc_setup = .true.')
      ! assign the values in nnkpts_block to nnlist and nncell
      ! this keeps track of how many neighbours have been seen for each k-point
      if (allocated(nnkpts_idx)) deallocate (nnkpts_idx)
      allocate (nnkpts_idx(num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nnkpts_idx in param_read')
      nnkpts_idx = 1
      ! allocating "global" nnlist & nncell
      ! These are deallocated in kmesh_dealloc
      if (allocated(nnlist)) deallocate (nnlist)
      allocate (nnlist(num_kpts, nntot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nnlist in param_read')
      if (allocated(nncell)) deallocate (nncell)
      allocate (nncell(3, num_kpts, nntot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nncell in param_read')
      do i = 1, num_kpts*nntot
        k = nnkpts_block(1, i)
        nnlist(k, nnkpts_idx(k)) = nnkpts_block(2, i)
        nncell(:, k, nnkpts_idx(k)) = nnkpts_block(3:, i)
        nnkpts_idx(k) = nnkpts_idx(k) + 1
      end do
      ! check that all k-points have the same number of neighbours
      if (any(nnkpts_idx /= (/(nntot + 1, i=1, num_kpts)/))) then
        call io_error('Inconsistent number of nearest neighbours.')
      end if
      deallocate (nnkpts_idx, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating nnkpts_idx in param_read')
      deallocate (nnkpts_block, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating nnkpts_block in param_read')
    end if

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

    ! To be called after having read the global flag
    call get_module_kmesh(moduleprefix='boltz', &
                          should_be_defined=boltzwann, &
                          module_kmesh=boltz_kmesh, &
                          module_kmesh_spacing=boltz_kmesh_spacing)

    call get_module_kmesh(moduleprefix='berry', &
                          should_be_defined=berry, &
                          module_kmesh=berry_kmesh, &
                          module_kmesh_spacing=berry_kmesh_spacing)

    call get_module_kmesh(moduleprefix='gyrotropic', &
                          should_be_defined=gyrotropic, &
                          module_kmesh=gyrotropic_kmesh, &
                          module_kmesh_spacing=gyrotropic_kmesh_spacing)

    call get_module_kmesh(moduleprefix='spin', &
                          should_be_defined=spin_moment, &
                          module_kmesh=spin_kmesh, &
                          module_kmesh_spacing=spin_kmesh_spacing)

    call get_module_kmesh(moduleprefix='dos', &
                          should_be_defined=dos, &
                          module_kmesh=dos_kmesh, &
                          module_kmesh_spacing=dos_kmesh_spacing)

    ! Atoms
    if (.not. library) num_atoms = 0
    call param_get_block_length('atoms_frac', found, i_temp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <atoms_frac> in input file'
    call param_get_block_length('atoms_cart', found2, i_temp2, lunits)
    if (found2 .and. library) write (stdout, '(a)') ' Ignoring <atoms_cart> in input file'
    if (.not. library) then
      if (found .and. found2) call io_error('Error: Cannot specify both atoms_frac and atoms_cart')
      if (found .and. i_temp > 0) then
        lunits = .false.
        num_atoms = i_temp
      elseif (found2 .and. i_temp2 > 0) then
        num_atoms = i_temp2
        if (lunits) num_atoms = num_atoms - 1
      end if
      if (num_atoms > 0) then
        call param_get_atoms(lunits)
      end if
    endif

    ! Projections
    auto_projections = .false.
    call param_get_keyword('auto_projections', found, l_value=auto_projections)
    num_proj = 0
    call param_get_block_length('projections', found, i_temp)
    ! check to see that there are no unrecognised keywords
    if (found) then
      if (auto_projections) call io_error('Error: Cannot specify both auto_projections and projections block')
      lhasproj = .true.
      call param_get_projections(num_proj, lcount=.true.)
    else
      if (guiding_centres .and. .not. (gamma_only .and. use_bloch_phases)) &
        call io_error('param_read: Guiding centres requested, but no projection block found')
      lhasproj = .false.
      num_proj = num_wann
    end if

    lselproj = .false.
    num_select_projections = 0
    call param_get_range_vector('select_projections', found, num_select_projections, lcount=.true.)
    if (found) then
      if (num_select_projections < 1) call io_error('Error: problem reading select_projections')
      if (allocated(select_projections)) deallocate (select_projections)
      allocate (select_projections(num_select_projections), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating select_projections in param_read')
      call param_get_range_vector('select_projections', found, num_select_projections, .false., select_projections)
      if (any(select_projections < 1)) &
        call io_error('Error: select_projections must contain positive numbers')
      if (num_select_projections < num_wann) &
        call io_error('Error: too few projections selected')
      if (num_select_projections > num_wann) &
        call io_error('Error: too many projections selected')
      if (.not. lhasproj) &
        call io_error('Error: select_projections cannot be used without defining the projections')
      if (maxval(select_projections(:)) > num_proj) &
        call io_error('Error: select_projections contains a number greater than num_proj')
      lselproj = .true.
    end if

    if (allocated(proj2wann_map)) deallocate (proj2wann_map)
    allocate (proj2wann_map(num_proj), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating proj2wann_map in param_read')
    proj2wann_map = -1

    if (lselproj) then
      do i = 1, num_proj
        do j = 1, num_wann
          if (select_projections(j) == i) proj2wann_map(i) = j
        enddo
      enddo
    else
      do i = 1, num_wann
        proj2wann_map(i) = i
      enddo
    endif

    if (lhasproj) call param_get_projections(num_proj, lcount=.false.)

    ! Constrained centres
    call param_get_block_length('slwf_centres', found, i_temp)
    if (found) then
      if (slwf_constrain) then
        ! Allocate array for constrained centres
        call param_get_centre_constraints
      else
        write (stdout, '(a)') ' slwf_constrain set to false. Ignoring <slwf_centres> block '
      end if
      ! Check that either projections or constrained centres are specified if slwf_constrain=.true.
    elseif (.not. found) then
      if (slwf_constrain) then
        if (.not. allocated(proj_site)) then
          call io_error('Error: slwf_constrain = true, but neither &
               & <slwf_centre> block  nor &
               & <projection_block> are specified.')
        else
          ! Allocate array for constrained centres
          call param_get_centre_constraints
        end if
      end if
    end if
    ! Warning
    if (slwf_constrain .and. allocated(proj_site) .and. .not. found) &
         & write (stdout, '(a)') ' Warning: No <slwf_centres> block found, but slwf_constrain set to true. &
           & Desired centres for SLWF same as projection centres.'

302 continue

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

    if (transport .and. tran_read_ht) goto 303

    ! For aesthetic purposes, convert some things to uppercase
    call param_uppercase()

303 continue

    deallocate (in_data, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating in_data in param_read')

    if (transport .and. tran_read_ht) return

    ! =============================== !
    ! Some checks and initialisations !
    ! =============================== !

!    if (restart.ne.' ') disentanglement=.false.

    if (disentanglement) then
      if (allocated(ndimwin)) deallocate (ndimwin)
      allocate (ndimwin(num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating ndimwin in param_read')
      if (allocated(lwindow)) deallocate (lwindow)
      allocate (lwindow(num_bands, num_kpts), stat=ierr)
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
    omega_total = -999.0_dp
    omega_tilde = -999.0_dp
    omega_invariant = -999.0_dp
    have_disentangled = .false.

    if (allocated(wannier_centres)) deallocate (wannier_centres)
    allocate (wannier_centres(3, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wannier_centres in param_read')
    wannier_centres = 0.0_dp
    if (allocated(wannier_spreads)) deallocate (wannier_spreads)
    allocate (wannier_spreads(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating wannier_spreads in param_read')
    wannier_spreads = 0.0_dp

    return

105 call io_error('Error: Problem opening eigenvalue file '//trim(seedname)//'.eig')
106 call io_error('Error: Problem reading eigenvalue file '//trim(seedname)//'.eig')

  end subroutine param_read

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

  subroutine get_module_kmesh(moduleprefix, should_be_defined, module_kmesh, module_kmesh_spacing)
    !! This function reads and sets the interpolation mesh variables needed by a given module
    !>
    !!  This function MUST be called after having read the global kmesh and kmesh_spacing!!
    !!  if the user didn't provide an interpolation_mesh_spacing, it is set to -1, so that
    !!       one can check in the code what the user asked for
    !!  The function takes care also of setting the default value to the global one if no local
    !!       keyword is defined
    use w90_io, only: io_error
    character(len=*), intent(in)       :: moduleprefix
    !!The prefix that is appended before the name of the variables. In particular,
    !!if the prefix is for instance XXX, the two variables that are read from the
    !!input file are XXX_kmesh and XXX_kmesh_spacing.
    logical, intent(in)                :: should_be_defined
    !! A logical flag: if it is true, at the end the code stops if no value is specified.
    !! Define it to .false. if no check should be performed.
    !! Often, you can simply pass the flag that activates the module itself.
    integer, dimension(3), intent(out) :: module_kmesh
    !! the integer array (length 3) where the interpolation mesh will be saved
    real(kind=dp), intent(out)         :: module_kmesh_spacing
    !! the real number on which the min mesh spacing is saved. -1 if it the
    !!user specifies in input the mesh and not the mesh_spacing

    logical :: found, found2
    integer :: i

    ! Default values
    module_kmesh_spacing = -1._dp
    module_kmesh = 0
    call param_get_keyword(trim(moduleprefix)//'_kmesh_spacing', found, r_value=module_kmesh_spacing)
    if (found) then
      if (module_kmesh_spacing .le. 0._dp) &
        call io_error('Error: '//trim(moduleprefix)//'_kmesh_spacing must be greater than zero')

      call internal_set_kmesh(module_kmesh_spacing, recip_lattice, module_kmesh)
    end if
    call param_get_vector_length(trim(moduleprefix)//'_kmesh', found2, length=i)
    if (found2) then
      if (found) &
        call io_error('Error: cannot set both '//trim(moduleprefix)//'_kmesh and ' &
                      //trim(moduleprefix)//'_kmesh_spacing')
      if (i .eq. 1) then
        call param_get_keyword_vector(trim(moduleprefix)//'_kmesh', found2, 1, i_value=module_kmesh)
        module_kmesh(2) = module_kmesh(1)
        module_kmesh(3) = module_kmesh(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector(trim(moduleprefix)//'_kmesh', found2, 3, i_value=module_kmesh)
      else
        call io_error('Error: '//trim(moduleprefix)// &
                      '_kmesh must be provided as either one integer or a vector of 3 integers')
      end if
      if (any(module_kmesh <= 0)) &
        call io_error('Error: '//trim(moduleprefix)//'_kmesh elements must be greater than zero')
    end if

    if ((found .eqv. .false.) .and. (found2 .eqv. .false.)) then
      ! This is the case where no  "local" interpolation k-mesh is provided in the input
      if (global_kmesh_set) then
        module_kmesh = kmesh
        ! I set also boltz_kmesh_spacing so that I can check if it is < 0 or not, and if it is
        ! > 0 I can print on output the mesh spacing that was chosen
        module_kmesh_spacing = kmesh_spacing
      else
        if (should_be_defined) &
          call io_error('Error: '//trim(moduleprefix)//' module required, but no interpolation mesh given.')
      end if
    end if
  end subroutine get_module_kmesh

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

    character(len=4)   :: orderstr

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
    do nsp = 1, num_species
      ic = ichar(atoms_label(nsp) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atoms_label(nsp) (1:1) = char(ic + ichar('Z') - ichar('z'))
    enddo

    do nsp = 1, num_species
      ic = ichar(atoms_symbol(nsp) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atoms_symbol(nsp) (1:1) = char(ic + ichar('Z') - ichar('z'))
    enddo

    ! Bands labels (eg, x --> X)
    do loop = 1, bands_num_spec_points
      do inner_loop = 1, len(bands_label(loop))
        ic = ichar(bands_label(loop) (inner_loop:inner_loop))
        if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
          bands_label(loop) (inner_loop:inner_loop) = char(ic + ichar('Z') - ichar('z'))
      enddo
    enddo

    ! Length unit (ang --> Ang, bohr --> Bohr)
    ic = ichar(length_unit(1:1))
    if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
      length_unit(1:1) = char(ic + ichar('Z') - ichar('z'))

    return

  end subroutine param_uppercase

!===================================================================
  subroutine param_write
    !==================================================================!
    !                                                                  !
    !! write wannier90 parameters to stdout
    !                                                                  !
    !===================================================================

    implicit none

    integer :: i, nkp, loop, nat, nsp

    if (transport .and. tran_read_ht) goto 401

    ! System
    write (stdout, *)
    write (stdout, '(36x,a6)') '------'
    write (stdout, '(36x,a6)') 'SYSTEM'
    write (stdout, '(36x,a6)') '------'
    write (stdout, *)
    if (lenconfac .eq. 1.0_dp) then
      write (stdout, '(30x,a21)') 'Lattice Vectors (Ang)'
    else
      write (stdout, '(28x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write (stdout, 101) 'a_1', (real_lattice(1, I)*lenconfac, i=1, 3)
    write (stdout, 101) 'a_2', (real_lattice(2, I)*lenconfac, i=1, 3)
    write (stdout, 101) 'a_3', (real_lattice(3, I)*lenconfac, i=1, 3)
    write (stdout, *)
    write (stdout, '(19x,a17,3x,f11.5)', advance='no') &
      'Unit Cell Volume:', cell_volume*lenconfac**3
    if (lenconfac .eq. 1.0_dp) then
      write (stdout, '(2x,a7)') '(Ang^3)'
    else
      write (stdout, '(2x,a8)') '(Bohr^3)'
    endif
    write (stdout, *)
    if (lenconfac .eq. 1.0_dp) then
      write (stdout, '(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
      write (stdout, '(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write (stdout, 101) 'b_1', (recip_lattice(1, I)/lenconfac, i=1, 3)
    write (stdout, 101) 'b_2', (recip_lattice(2, I)/lenconfac, i=1, 3)
    write (stdout, 101) 'b_3', (recip_lattice(3, I)/lenconfac, i=1, 3)
    write (stdout, *) ' '
    ! Atoms
    if (num_atoms > 0) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
      else
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, num_species
        do nat = 1, atoms_species_num(nsp)
          write (stdout, '(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                 '|', atoms_symbol(nsp), nat, atoms_pos_frac(:, nat, nsp),&
  &                 '|', atoms_pos_cart(:, nat, nsp)*lenconfac, '|'
        end do
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
    else
      write (stdout, '(25x,a)') 'No atom positions specified'
    end if
    ! Constrained centres
    if (selective_loc .and. slwf_constrain) then
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '| Wannier#        Original Centres              Constrained centres          |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do i = 1, slwf_num
        write (stdout, '(1x,a1,2x,i3,2x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                    '|', i, ccentres_frac(i, :), '|', wannier_centres(:, i), '|'
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
    end if
    ! Projections
    if (iprint > 1 .and. allocated(input_proj_site)) then
      write (stdout, '(32x,a)') '-----------'
      write (stdout, '(32x,a)') 'PROJECTIONS'
      write (stdout, '(32x,a)') '-----------'
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a)') '|     Frac. Coord.   l mr  r        z-axis               x-axis          Z/a |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, num_proj
        write (stdout, '(1x,a1,3(1x,f5.2),1x,i2,1x,i2,1x,i2,3(1x,f6.3),3(1x,f6.3),2x,f4.1,1x,a1)')&
  &              '|', input_proj_site(1, nsp), input_proj_site(2, nsp), &
             input_proj_site(3, nsp), input_proj_l(nsp), input_proj_m(nsp), input_proj_radial(nsp), &
             input_proj_z(1, nsp), input_proj_z(2, nsp), input_proj_z(3, nsp), input_proj_x(1, nsp), &
             input_proj_x(2, nsp), input_proj_x(3, nsp), input_proj_zona(nsp), '|'
      end do
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, *) ' '
    end if

    if (iprint > 1 .and. lselproj .and. allocated(proj_site)) then
      write (stdout, '(30x,a)') '--------------------'
      write (stdout, '(30x,a)') 'SELECTED PROJECTIONS'
      write (stdout, '(30x,a)') '--------------------'
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a)') '|     Frac. Coord.   l mr  r        z-axis               x-axis          Z/a |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, num_wann
        if (proj2wann_map(nsp) < 0) cycle
        write (stdout, '(1x,a1,3(1x,f5.2),1x,i2,1x,i2,1x,i2,3(1x,f6.3),3(1x,f6.3),2x,f4.1,1x,a1)')&
  &              '|', proj_site(1, nsp), proj_site(2, nsp), &
             proj_site(3, nsp), proj_l(nsp), proj_m(nsp), proj_radial(nsp), &
             proj_z(1, nsp), proj_z(2, nsp), proj_z(3, nsp), proj_x(1, nsp), &
             proj_x(2, nsp), proj_x(3, nsp), proj_zona(nsp), '|'
      end do
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, *) ' '
    end if

    ! K-points
    write (stdout, '(32x,a)') '------------'
    write (stdout, '(32x,a)') 'K-POINT GRID'
    write (stdout, '(32x,a)') '------------'
    write (stdout, *) ' '
    write (stdout, '(13x,a,i3,1x,a1,i3,1x,a1,i3,6x,a,i5)') 'Grid size =', mp_grid(1), 'x', mp_grid(2), 'x', mp_grid(3), &
      'Total points =', num_kpts
    write (stdout, *) ' '
    if (iprint > 1) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |'
      else
        write (stdout, '(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nkp = 1, num_kpts
        write (stdout, '(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|', nkp, kpt_latt(:, nkp), '|', &
          kpt_cart(:, nkp)/lenconfac, '|'
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, *) ' '
    end if
    ! Main
    write (stdout, *) ' '
    write (stdout, '(1x,a78)') '*---------------------------------- MAIN ------------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Wannier Functions               :', num_wann, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Objective Wannier Functions     :', slwf_num, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of input Bloch states              :', num_bands, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :', iprint, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :', timing_level, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Optimisation (0=memory, 3=speed)          :', optimisation, '|'
    write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :', trim(length_unit), '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Post-processing setup (write *.nnkp)      :', postproc_setup, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using Gamma-only branch of algorithms     :', gamma_only, '|'
    !YN: RS:
    if (lsitesymmetry) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using symmetry-adapted WF mode            :', lsitesymmetry, '|'
      write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Tolerance for symmetry condition on U     :', symmetrize_eps, '|'
    endif

    if (cp_pp .or. iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  CP code post-processing                   :', cp_pp, '|'
    if (wannier_plot .or. iprint > 2) then
      if (wvfn_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  Wavefunction (UNK) file-type              :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  Wavefunction (UNK) file-type              :', 'unformatted', '|'
      endif
      if (spin == 1) then
        write (stdout, '(1x,a46,16x,a2,13x,a1)') '|  Wavefunction spin channel                 :', 'up', '|'
      else
        write (stdout, '(1x,a46,14x,a4,13x,a1)') '|  Wavefunction spin channel                 :', 'down', '|'
      endif
    endif

    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    ! Wannierise
    write (stdout, '(1x,a78)') '*------------------------------- WANNIERISE ---------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Total number of iterations                :', num_iter, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of CG steps before reset           :', num_cg_steps, '|'
    if (lfixstep) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fixed step length for minimisation        :', fixed_step, '|'
    else
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Trial step length for line search         :', trial_step, '|'
    endif
    write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Convergence tolerence                     :', conv_tol, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Convergence window                        :', conv_window, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between writing output         :', num_print_cycles, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between backing up to disk     :', num_dump_cycles, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write r^2_nm to file                      :', write_r2mn, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write xyz WF centres to file              :', write_xyz, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write on-site energies <0n|H|0n> to file  :', write_hr_diag, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use guiding centre to control phases      :', guiding_centres, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use phases for initial projections        :', use_bloch_phases, '|'
    if (guiding_centres .or. iprint > 2) then
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations before starting guiding centres:', num_no_guide_iter, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between using guiding centres  :', num_guide_cycles, '|'
    end if
    if (selective_loc .or. iprint > 2) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Perform selective localization            :', selective_loc, '|'
    end if
    if (slwf_constrain .or. iprint > 2) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use constrains in selective localization  :', slwf_constrain, '|'
      write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Value of the Lagrange multiplier          :',&
           &slwf_lambda, '|'
    end if
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    !
    ! Disentanglement
    !
    if (disentanglement .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------------- DISENTANGLE --------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using band disentanglement                :', disentanglement, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Total number of iterations                :', dis_num_iter, '|'
      write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|  Mixing ratio                              :', dis_mix_ratio, '|'
      write (stdout, '(1x,a46,8x,ES10.3,13x,a1)') '|  Convergence tolerence                     :', dis_conv_tol, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Convergence window                        :', dis_conv_window, '|'
      ! GS-start
      if (dis_spheres_num .gt. 0) then
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of spheres in k-space              :', dis_spheres_num, '|'
        do nkp = 1, dis_spheres_num
          write (stdout, '(1x,a13,I4,a2,2x,3F8.3,a15,F8.3,9x,a1)') &
            '|   center n.', nkp, ' :', dis_spheres(1:3, nkp), ',    radius   =', dis_spheres(4, nkp), '|'
        enddo
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Index of first Wannier band               :', dis_spheres_first_wann, '|'
      endif
      ! GS-end
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    end if
    !
    ! Plotting
    !
    if (wannier_plot .or. bands_plot .or. fermi_surface_plot .or. kslice &
        .or. dos_plot .or. write_hr .or. iprint > 2) then
      !
      write (stdout, '(1x,a78)') '*-------------------------------- PLOTTING ----------------------------------*'
      !
      if (wannier_plot .or. iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Wannier functions                :', wannier_plot, '|'
        write (stdout, '(1x,a46,1x,I5,a1,I5,a1,I5,13x,a1)') &
          '|   Size of supercell for plotting           :', &
          wannier_plot_supercell(1), 'x', wannier_plot_supercell(2), 'x', &
          wannier_plot_supercell(3), '|'

        if (translate_home_cell) then
          write (stdout, '(1x,a46,10x,L8,13x,a1)') &
            '|  Translating WFs to home cell              :', translate_home_cell, '|'
        end if

        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting mode (molecule or crystal)      :', trim(wannier_plot_mode), '|'
        if (spinors) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting mode for spinor WFs             :', &
            trim(wannier_plot_spinor_mode), '|'
          write (stdout, '(1x,a46,10x,L8,13x,a1)') '|   Include phase for spinor WFs             :', &
            wannier_plot_spinor_phase, '|'
        end if
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :', trim(wannier_plot_format), '|'
        if (index(wannier_plot_format, 'cub') > 0 .or. iprint > 2) then
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Plot radius                              :', wannier_plot_radius, '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Plot scale                               :', wannier_plot_scale, '|'
        endif
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (fermi_surface_plot .or. iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Fermi surface                    :', fermi_surface_plot, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Number of plotting points (along b_1)    :', fermi_surface_num_points, '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :' &
          , trim(fermi_surface_plot_format), '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (bands_plot .or. iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting interpolated bandstructure       :', bands_plot, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Number of K-path sections                :', bands_num_spec_points/2, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Divisions along first K-path section     :', bands_num_points, '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Output format                            :', trim(bands_plot_format), '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Output mode                              :', trim(bands_plot_mode), '|'
        if (index(bands_plot_mode, 'cut') .ne. 0) then
          write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Dimension of the system                  :', bands_plot_dim, '|'
          if (bands_plot_dim .eq. 1) &
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System extended in                       :', trim(one_dim_axis), '|'
          if (bands_plot_dim .eq. 2) &
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System confined in                       :', trim(one_dim_axis), '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off value                :', hr_cutoff, '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off distance             :', dist_cutoff, '|'
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian cut-off distance mode        :', trim(dist_cutoff_mode), '|'
        endif
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
        write (stdout, '(1x,a78)') '|   K-space path sections:                                                   |'
        if (bands_num_spec_points == 0) then
          write (stdout, '(1x,a78)') '|     None defined                                                           |'
        else
          do loop = 1, bands_num_spec_points, 2
            write (stdout, '(1x,a10,1x,a5,1x,3F7.3,5x,a3,1x,a5,1x,3F7.3,3x,a1)') '|    From:', bands_label(loop), &
              (bands_spec_points(i, loop), i=1, 3), 'To:', bands_label(loop + 1), (bands_spec_points(i, loop + 1), i=1, 3), '|'
          end do
        end if
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (write_hr .or. iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Hamiltonian in WF basis          :', write_hr, '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      endif
      if (write_vdw_data .or. iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Writing data for Van der Waals post-proc  :', write_vdw_data, '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      endif
      !
    endif

401 continue
    !
    ! Transport
    !
    if (transport .or. iprint > 2) then
      !
      write (stdout, '(1x,a78)') '*------------------------------- TRANSPORT ----------------------------------*'
      !
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Transport mode                            :', trim(transport_mode), '|'
      !
      if (tran_read_ht) then
        !
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian from external files          :', 'T', '|'
        !
      else
        !
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian from external files          :', 'F', '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System extended in                       :', trim(one_dim_axis), '|'
        !
      end if

      write (stdout, '(1x,a78)') '|   Centre of the unit cell to which WF are translated (fract. coords):      |'
      write (stdout, '(1x,a1,35x,F12.6,a1,F12.6,a1,F12.6,3x,a1)') '|', translation_centre_frac(1), ',', &
        translation_centre_frac(2), ',', &
        translation_centre_frac(3), '|'

      if (size(fermi_energy_list) == 1) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi energy (eV)                         :', fermi_energy_list(1), '|'
      else
        write (stdout, '(1x,a21,I8,a12,f8.3,a4,f8.3,a3,13x,a1)') '|  Fermi energy     :', size(fermi_energy_list), &
          ' steps from ', fermi_energy_list(1), ' to ', &
          fermi_energy_list(size(fermi_energy_list)), ' eV', '|'
      end if
      !
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      !
    endif

101 format(20x, a3, 2x, 3F11.6)

  end subroutine param_write

!===================================================================
  subroutine param_postw90_write
    !==================================================================!
    !                                                                  !
    !! write postw90 parameters to stdout
    !                                                                  !
    !===================================================================

    implicit none

    integer :: i, loop, nat, nsp

    ! System
    write (stdout, *)
    write (stdout, '(36x,a6)') '------'
    write (stdout, '(36x,a6)') 'SYSTEM'
    write (stdout, '(36x,a6)') '------'
    write (stdout, *)
    if (lenconfac .eq. 1.0_dp) then
      write (stdout, '(30x,a21)') 'Lattice Vectors (Ang)'
    else
      write (stdout, '(28x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write (stdout, 101) 'a_1', (real_lattice(1, I)*lenconfac, i=1, 3)
    write (stdout, 101) 'a_2', (real_lattice(2, I)*lenconfac, i=1, 3)
    write (stdout, 101) 'a_3', (real_lattice(3, I)*lenconfac, i=1, 3)
    write (stdout, *)
    write (stdout, '(19x,a17,3x,f11.5)', advance='no') &
      'Unit Cell Volume:', cell_volume*lenconfac**3
    if (lenconfac .eq. 1.0_dp) then
      write (stdout, '(2x,a7)') '(Ang^3)'
    else
      write (stdout, '(2x,a8)') '(Bohr^3)'
    endif
    write (stdout, *)
    if (lenconfac .eq. 1.0_dp) then
      write (stdout, '(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
      write (stdout, '(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write (stdout, 101) 'b_1', (recip_lattice(1, I)/lenconfac, i=1, 3)
    write (stdout, 101) 'b_2', (recip_lattice(2, I)/lenconfac, i=1, 3)
    write (stdout, 101) 'b_3', (recip_lattice(3, I)/lenconfac, i=1, 3)
    write (stdout, *) ' '
    ! Atoms
    if (num_atoms > 0) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
      else
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, num_species
        do nat = 1, atoms_species_num(nsp)
          write (stdout, '(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                 '|', atoms_symbol(nsp), nat, atoms_pos_frac(:, nat, nsp),&
  &                 '|', atoms_pos_cart(:, nat, nsp)*lenconfac, '|'
        end do
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
    else
      write (stdout, '(25x,a)') 'No atom positions specified'
    end if
    write (stdout, *) ' '
    ! Main
    write (stdout, *) ' '

    write (stdout, '(1x,a78)') '*-------------------------------- POSTW90 -----------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Wannier Functions               :', num_wann, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of electrons per state             :', num_elec_per_state, '|'
    if (abs(scissors_shift) > 1.0e-7_dp .or. iprint > 0) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Scissor shift applied to conduction bands :', scissors_shift, '|'
      if (num_valence_bands > 0) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Number of valence bands                   :', num_valence_bands, '|'
      else
        write (stdout, '(1x,a78)') '|  Number of valence bands                   :       not defined             |'
      endif
    endif
    if (spin_decomp .or. iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Spin decomposition                        :', spin_decomp, '|'
    if (spin_moment .or. iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Spin moment                       :', spin_moment, '|'
    if (spin_decomp .or. spin_moment .or. iprint > 2) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Polar angle of spin quantisation axis     :', spin_axis_polar, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Azimuthal angle of spin quantisation axis :', spin_axis_azimuth, '|'
      if (spn_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  Spn file-type                   :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  Spn file-type                   :', 'unformatted', '|'
      endif
      if (uHu_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  uHu file-type                   :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  uHu file-type                   :', 'unformatted', '|'
      endif
    end if

    if (size(fermi_energy_list) == 1) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi energy (eV)                         :', fermi_energy_list(1), '|'
    else
      write (stdout, '(1x,a21,I8,a12,f8.3,a4,f8.3,a3,13x,a1)') '|  Fermi energy     :', size(fermi_energy_list), &
        ' steps from ', fermi_energy_list(1), ' to ', &
        fermi_energy_list(size(fermi_energy_list)), ' eV', '|'
    end if

    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :', iprint, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :', timing_level, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Optimisation (0=memory, 3=speed)          :', optimisation, '|'
    write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :', trim(length_unit), '|'
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    write (stdout, '(1x,a78)') '*------------------------ Global Smearing Parameters ------------------------*'
    if (adpt_smr) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', adpt_smr_fac, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width (eV)       :', adpt_smr_max, '|'

    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', smr_fixed_en_width, '|'
    endif
    write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', trim(param_get_smearing_type(smr_index)), '|'
    if (global_kmesh_set) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       T', '|'
      if (kmesh_spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          kmesh(1), 'x', kmesh(2), 'x', kmesh(3), ' Spacing = ', kmesh_spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , kmesh(1), 'x', kmesh(2), 'x', kmesh(3), '|'
      endif
    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       F', '|'
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    ! DOS
    if (dos .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*---------------------------------- DOS -------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Density of States                :', dos, '|'
      if (num_dos_project > 1) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', dos_energy_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', dos_energy_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', dos_energy_step, '|'
      if (dos_adpt_smr .eqv. adpt_smr .and. dos_adpt_smr_fac == adpt_smr_fac .and. dos_adpt_smr_max == adpt_smr_max &
          .and. dos_smr_fixed_en_width == smr_fixed_en_width .and. smr_index == dos_smr_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (dos_adpt_smr) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', dos_adpt_smr_fac, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', dos_adpt_smr_max, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', dos_smr_fixed_en_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', trim(param_get_smearing_type(dos_smr_index)), '|'
      endif
      if (kmesh(1) == dos_kmesh(1) .and. kmesh(2) == dos_kmesh(2) .and. kmesh(3) == dos_kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (dos_kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            dos_kmesh(1), 'x', dos_kmesh(2), 'x', dos_kmesh(3), ' Spacing = ', dos_kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , dos_kmesh(1), 'x', dos_kmesh(2), 'x', dos_kmesh(3), '|'
        endif
      endif
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    if (kpath .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KPATH ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a path in k-space   :', kpath, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', kpath_num_points, '|'
      if (index(kpath_task, 'bands') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       F', '|'
      endif
      if (index(kpath_task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       F', '|'
      endif
      if (index(kpath_task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(kpath_task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the bands    :', trim(kpath_bands_colour), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      write (stdout, '(1x,a78)') '|   K-space path sections:                                                   |'
      if (bands_num_spec_points == 0) then
        write (stdout, '(1x,a78)') '|     None defined                                                           |'
      else
        do loop = 1, bands_num_spec_points, 2
          write (stdout, '(1x,a10,2x,a1,2x,3F7.3,5x,a3,2x,a1,2x,3F7.3,7x,a1)') '|    From:', bands_label(loop), &
            (bands_spec_points(i, loop), i=1, 3), 'To:', bands_label(loop + 1), (bands_spec_points(i, loop + 1), i=1, 3), '|'
        end do
      end if
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (kslice .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KSLICE -----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a slice in k-space  :', kslice, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi level used for slice                :', fermi_energy_list(1), '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', kpath_num_points, '|'
      if (index(kslice_task, 'fermi_lines') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       F', '|'
      endif
      if (index(kslice_task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       F', '|'
      endif
      if (index(kslice_task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(kslice_task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the lines    :', &
        trim(kslice_fermi_lines_colour), '|'
      write (stdout, '(1x,a78)') '|  2D slice parameters (in reduced coordinates):                             |'
      write (stdout, '(1x,a14,2x,3F8.3,37x,a1)') '|     Corner: ', (kslice_corner(i), i=1, 3), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector1: ', (kslice_b1(i), i=1, 3), ' Divisions:', kslice_2dkmesh(1), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector2: ', (kslice_b2(i), i=1, 3), ' Divisions:', kslice_2dkmesh(1), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (berry .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- BERRY ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Berry Phase related properties    :', berry, '|'
      if (index(berry_task, 'kubo') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       F', '|'
      endif
      if (index(berry_task, 'ahc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       F', '|'
      endif
      if (index(berry_task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       F', '|'
      endif
      if (index(berry_task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       F', '|'
      endif
      if (index(berry_task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for optical responses     :', kubo_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency for optical responses     :', kubo_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for optical responses           :', kubo_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue for optical responses    :', kubo_eigval_max, '|'
      if (index(berry_task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing factor for shift current         :', sc_eta, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Frequency theshold for shift current      :', sc_w_thr, '|'
        write (stdout, '(1x,a46,1x,a27,3x,a1)') '|  Bloch sums                                :', &
          trim(param_get_convention_type(sc_phase_conv)), '|'
      end if
      if (kubo_adpt_smr .eqv. adpt_smr .and. kubo_adpt_smr_fac == adpt_smr_fac .and. kubo_adpt_smr_max == adpt_smr_max &
          .and. kubo_smr_fixed_en_width == smr_fixed_en_width .and. smr_index == kubo_smr_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (kubo_adpt_smr) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', kubo_adpt_smr_fac, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', kubo_adpt_smr_max, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', kubo_smr_fixed_en_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', trim(param_get_smearing_type(kubo_smr_index)), '|'
      endif
      if (kmesh(1) == berry_kmesh(1) .and. kmesh(2) == berry_kmesh(2) .and. kmesh(3) == berry_kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (berry_kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            berry_kmesh(1), 'x', berry_kmesh(2), 'x', berry_kmesh(3), ' Spacing = ', berry_kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , berry_kmesh(1), 'x', berry_kmesh(2), 'x', berry_kmesh(3), '|'
        endif
      endif
      if (berry_curv_adpt_kmesh > 1) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Using an adaptive refinement mesh of size :', berry_curv_adpt_kmesh, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Threshold for adaptive refinement         :', &
          berry_curv_adpt_kmesh_thresh, '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    none', '|'
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (gyrotropic .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- GYROTROPIC   ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '| Compute Gyrotropic properties              :', gyrotropic, '|'
      write (stdout, '(1x,a46,10x,a20,1x,a1)') '| gyrotropic_task                            :', gyrotropic_task, '|'
      call parameters_gyro_write_task(gyrotropic_task, '-d0', 'calculate the D tensor')
      call parameters_gyro_write_task(gyrotropic_task, '-dw', 'calculate the tildeD tensor')
      call parameters_gyro_write_task(gyrotropic_task, '-c', 'calculate the C tensor')
      call parameters_gyro_write_task(gyrotropic_task, '-k', 'calculate the K tensor')
      call parameters_gyro_write_task(gyrotropic_task, '-noa', 'calculate the interbad natural optical activity')
      call parameters_gyro_write_task(gyrotropic_task, '-dos', 'calculate the density of states')

      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for tildeD,NOA            :', gyrotropic_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency                           :', gyrotropic_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for frequency                   :', gyrotropic_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue                          :', gyrotropic_eigval_max, '|'
      if (gyrotropic_smr_fixed_en_width == smr_fixed_en_width .and. smr_index == gyrotropic_smr_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        write (stdout, '(1x,a78)') '|  Using local  smearing parameters                                          |'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
        gyrotropic_smr_fixed_en_width, '|'
      write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function                         :', &
        trim(param_get_smearing_type(gyrotropic_smr_index)), '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  degen_thresh                              :', gyrotropic_degen_thresh, '|'

      if (kmesh(1) == gyrotropic_kmesh(1) .and. kmesh(2) == gyrotropic_kmesh(2) .and. kmesh(3) == gyrotropic_kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      elseif (gyrotropic_kmesh_spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          gyrotropic_kmesh(1), 'x', gyrotropic_kmesh(2), 'x', gyrotropic_kmesh(3), ' Spacing = ', gyrotropic_kmesh_spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , gyrotropic_kmesh(1), 'x', gyrotropic_kmesh(2), 'x', gyrotropic_kmesh(3), '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    not implemented', '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (boltzwann .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------------- BOLTZWANN ----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Boltzmann transport properties    :', boltzwann, '|'
      if (boltz_2d_dir_num > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  2d structure: non-periodic dimension  :', trim(boltz_2d_dir), '|'
      else
        write (stdout, '(1x,a78)') '|  3d Structure                              :                 T             |'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Relaxation Time (fs)                      :', boltz_relax_time, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Chemical Potential (eV)  :', boltz_mu_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Chemical Potential (eV)  :', boltz_mu_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Chemical Potential (eV)     :', boltz_mu_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Temperature (K)          :', boltz_temp_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Temperature (K)          :', boltz_temp_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Temperature (K)             :', boltz_temp_step, '|'

      if (kmesh(1) == boltz_kmesh(1) .and. kmesh(2) == boltz_kmesh(2) .and. kmesh(3) == boltz_kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (boltz_kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            boltz_kmesh(1), 'x', boltz_kmesh(2), 'x', boltz_kmesh(3), ' Spacing = ', boltz_kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , boltz_kmesh(1), 'x', boltz_kmesh(2), 'x', boltz_kmesh(3), '|'
        endif
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for TDF (eV)                    :', boltz_tdf_energy_step, '|'
      write (stdout, '(1x,a25,5x,a43,4x,a1)') '|  TDF Smearing Function ', trim(param_get_smearing_type(boltz_tdf_smr_index)), '|'
      if (boltz_tdf_smr_fixed_en_width > 0.0_dp) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
          '|  TDF fixed Smearing width (eV)             :', boltz_tdf_smr_fixed_en_width, '|'
      else
        write (stdout, '(1x,a78)') '|  TDF fixed Smearing width                  :         unsmeared             |'
      endif
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute DOS at same time                  :', boltz_calc_also_dos, '|'
      if (boltz_calc_also_dos .and. iprint > 2) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', boltz_dos_energy_min, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', boltz_dos_energy_max, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', boltz_dos_energy_step, '|'
        if (boltz_dos_adpt_smr .eqv. adpt_smr .and. boltz_dos_adpt_smr_fac == adpt_smr_fac &
            .and. boltz_dos_adpt_smr_max == adpt_smr_max &
            .and. boltz_dos_smr_fixed_en_width == smr_fixed_en_width .and. smr_index == boltz_dos_smr_index) then
          write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
        else
          if (boltz_dos_adpt_smr) then
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Adaptive width smearing               :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Adaptive smearing factor              :', boltz_dos_adpt_smr_fac, '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Maximum allowed smearing width        :', boltz_dos_adpt_smr_max, '|'
          else
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Fixed width smearing                  :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  DOS Smearing width                         :', &
              boltz_dos_smr_fixed_en_width, '|'
          endif
          write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', trim(param_get_smearing_type(boltz_dos_smr_index)), '|'
        endif
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (geninterp .or. iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------Generic Band Interpolation--------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Properties at given k-points      :', geninterp, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Calculate band gradients                  :', geninterp_alsofirstder, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write data into a single file             :', geninterp_single_file, '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

101 format(20x, a3, 2x, 3F11.6)

  end subroutine param_postw90_write

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
    write (stdout, *) '            |        G. Pizzi et al., arXiv:1907:09788 (2019)   |'
    write (stdout, *) '            |        https://arxiv.org/abs/1907.09788           |'
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
    write (stdout, *) '            | Copyright (c) 1996-2019                           |'
    write (stdout, *) '            |        The Wannier90 Developer Group and          |'
    write (stdout, *) '            |        individual contributors                    |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |      Release: ', adjustl(w90_version), '  27th February 2019      |'
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

    if (allocated(ndimwin)) then
      deallocate (ndimwin, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ndimwin in param_dealloc')
    end if
    if (allocated(lwindow)) then
      deallocate (lwindow, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating lwindow in param_dealloc')
    end if
    if (allocated(eigval)) then
      deallocate (eigval, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating eigval in param_dealloc')
    endif
    if (allocated(shell_list)) then
      deallocate (shell_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating shell_list in param_dealloc')
    endif
    if (allocated(kpt_latt)) then
      deallocate (kpt_latt, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kpt_latt in param_dealloc')
    endif
    if (allocated(kpt_cart)) then
      deallocate (kpt_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kpt_cart in param_dealloc')
    endif
    if (allocated(bands_label)) then
      deallocate (bands_label, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating bands_label in param_dealloc')
    end if
    if (allocated(bands_spec_points)) then
      deallocate (bands_spec_points, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating bands_spec_points in param_dealloc')
    end if
    if (allocated(atoms_label)) then
      deallocate (atoms_label, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_label in param_dealloc')
    end if
    if (allocated(atoms_symbol)) then
      deallocate (atoms_symbol, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_symbol in param_dealloc')
    end if
    if (allocated(atoms_pos_frac)) then
      deallocate (atoms_pos_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atom_pos_frac in param_dealloc')
    end if
    if (allocated(atoms_pos_cart)) then
      deallocate (atoms_pos_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_pos_cart in param_dealloc')
    end if
    if (allocated(atoms_species_num)) then
      deallocate (atoms_species_num, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating atoms_species_num in param_dealloc')
    end if
    if (allocated(input_proj_site)) then
      deallocate (input_proj_site, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_site in param_dealloc')
    end if
    if (allocated(input_proj_l)) then
      deallocate (input_proj_l, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_l in param_dealloc')
    end if
    if (allocated(input_proj_m)) then
      deallocate (input_proj_m, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_m in param_dealloc')
    end if
    if (allocated(input_proj_s)) then
      deallocate (input_proj_s, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_s in param_dealloc')
    end if
    if (allocated(input_proj_s_qaxis)) then
      deallocate (input_proj_s_qaxis, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_s_qaxis in param_dealloc')
    end if
    if (allocated(input_proj_z)) then
      deallocate (input_proj_z, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_z in param_dealloc')
    end if
    if (allocated(input_proj_x)) then
      deallocate (input_proj_x, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_x in param_dealloc')
    end if
    if (allocated(input_proj_radial)) then
      deallocate (input_proj_radial, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_radial in param_dealloc')
    end if
    if (allocated(input_proj_zona)) then
      deallocate (input_proj_zona, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating input_proj_zona in param_dealloc')
    end if
    if (allocated(proj_site)) then
      deallocate (proj_site, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_site in param_dealloc')
    end if
    if (allocated(proj_l)) then
      deallocate (proj_l, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_l in param_dealloc')
    end if
    if (allocated(proj_m)) then
      deallocate (proj_m, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_m in param_dealloc')
    end if
    if (allocated(proj_s)) then
      deallocate (proj_s, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_s in param_dealloc')
    end if
    if (allocated(proj_s_qaxis)) then
      deallocate (proj_s_qaxis, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_s_qaxis in param_dealloc')
    end if
    if (allocated(proj_z)) then
      deallocate (proj_z, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_z in param_dealloc')
    end if
    if (allocated(proj_x)) then
      deallocate (proj_x, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_x in param_dealloc')
    end if
    if (allocated(proj_radial)) then
      deallocate (proj_radial, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_radial in param_dealloc')
    end if
    if (allocated(proj_zona)) then
      deallocate (proj_zona, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_zona in param_dealloc')
    end if
    if (allocated(wannier_plot_list)) then
      deallocate (wannier_plot_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_plot_list in param_dealloc')
    end if
    if (allocated(exclude_bands)) then
      deallocate (exclude_bands, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating exclude_bands in param_dealloc')
    end if
    if (allocated(wannier_centres)) then
      deallocate (wannier_centres, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_centres in param_dealloc')
    end if
    if (allocated(wannier_spreads)) then
      deallocate (wannier_spreads, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_spreads in param_dealloc')
    endif
    if (allocated(bands_plot_project)) then
      deallocate (bands_plot_project, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating bands_plot_project in param_dealloc')
    endif
    if (allocated(dos_project)) then
      deallocate (dos_project, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating dos_project in param_dealloc')
    endif
    if (allocated(fermi_energy_list)) then
      deallocate (fermi_energy_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating fermi_energy_list in param_dealloc')
    endif
    if (allocated(kubo_freq_list)) then
      deallocate (kubo_freq_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kubo_freq_list in param_dealloc')
    endif
    if (allocated(dis_spheres)) then
      deallocate (dis_spheres, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating dis_spheres in param_dealloc')
    endif
    if (allocated(ccentres_frac)) then
      deallocate (ccentres_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating ccentres_frac in param_dealloc')
    endif
    if (allocated(ccentres_cart)) then
      deallocate (ccentres_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating ccentres_cart in param_dealloc')
    end if
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
  subroutine param_write_chkpt(chkpt)
    !=================================================!
    !! Write checkpoint file
    !! IMPORTANT! If you change the chkpt format, adapt
    !! accordingly also the w90chk2chk.x utility!
    !=================================================!

    use w90_io, only: io_file_unit, io_date, seedname

    implicit none

    character(len=*), intent(in) :: chkpt

    integer :: chk_unit, nkp, i, j, k, l
    character(len=9) :: cdate, ctime
    character(len=33) :: header
    character(len=20) :: chkpt1

    write (stdout, '(/1x,3a)', advance='no') 'Writing checkpoint file ', trim(seedname), '.chk...'

    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk', form='unformatted')

    write (chk_unit) header                                   ! Date and time
    write (chk_unit) num_bands                                ! Number of bands
    write (chk_unit) num_exclude_bands                        ! Number of excluded bands
    write (chk_unit) (exclude_bands(i), i=1, num_exclude_bands) ! Excluded bands
    write (chk_unit) ((real_lattice(i, j), i=1, 3), j=1, 3)        ! Real lattice
    write (chk_unit) ((recip_lattice(i, j), i=1, 3), j=1, 3)       ! Reciprocal lattice
    write (chk_unit) num_kpts                                 ! Number of k-points
    write (chk_unit) (mp_grid(i), i=1, 3)                       ! M-P grid
    write (chk_unit) ((kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts) ! K-points
    write (chk_unit) nntot                  ! Number of nearest k-point neighbours
    write (chk_unit) num_wann               ! Number of wannier functions
    chkpt1 = adjustl(trim(chkpt))
    write (chk_unit) chkpt1                 ! Position of checkpoint
    write (chk_unit) have_disentangled      ! Whether a disentanglement has been performed
    if (have_disentangled) then
      write (chk_unit) omega_invariant     ! Omega invariant
      ! lwindow, ndimwin and U_matrix_opt
      write (chk_unit) ((lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
      write (chk_unit) (ndimwin(nkp), nkp=1, num_kpts)
      write (chk_unit) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
    endif
    write (chk_unit) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)               ! U_matrix
    write (chk_unit) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, nntot), l=1, num_kpts) ! M_matrix
    write (chk_unit) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)
    write (chk_unit) (wannier_spreads(i), i=1, num_wann)
    close (chk_unit)

    write (stdout, '(a/)') ' done'

    return

  end subroutine param_write_chkpt

!=================================================!
  subroutine param_read_chkpt()
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
    if (ntmp .ne. num_bands) call io_error('param_read_chk: Mismatch in num_bands')
    read (chk_unit) ntmp                           ! Number of excluded bands
    if (ntmp .ne. num_exclude_bands) &
      call io_error('param_read_chk: Mismatch in num_exclude_bands')
    read (chk_unit) (tmp_excl_bands(i), i=1, num_exclude_bands) ! Excluded bands
    do i = 1, num_exclude_bands
      if (tmp_excl_bands(i) .ne. exclude_bands(i)) &
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
        if (abs(tmp_kpt_latt(i, nkp) - kpt_latt(i, nkp)) .gt. eps6) &
          call io_error('param_read_chk: Mismatch in kpt_latt')
      enddo
    enddo
    read (chk_unit) ntmp                ! nntot
    if (ntmp .ne. nntot) &
      call io_error('param_read_chk: Mismatch in nntot')
    read (chk_unit) ntmp                ! num_wann
    if (ntmp .ne. num_wann) &
      call io_error('param_read_chk: Mismatch in num_wann')
    ! End of consistency checks

    read (chk_unit) checkpoint             ! checkpoint
    checkpoint = adjustl(trim(checkpoint))

    read (chk_unit) have_disentangled      ! whether a disentanglement has been performed

    if (have_disentangled) then

      read (chk_unit) omega_invariant     ! omega invariant

      ! lwindow
      if (.not. allocated(lwindow)) then
        allocate (lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating lwindow in param_read_chkpt')
      endif
      read (chk_unit, err=122) ((lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)

      ! ndimwin
      if (.not. allocated(ndimwin)) then
        allocate (ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ndimwin in param_read_chkpt')
      endif
      read (chk_unit, err=123) (ndimwin(nkp), nkp=1, num_kpts)

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
    if (.not. allocated(m_matrix)) then
      allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating m_matrix in param_read_chkpt')
    endif
    read (chk_unit, err=126) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, nntot), l=1, num_kpts)

    ! wannier_centres
    read (chk_unit, err=127) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)

    ! wannier spreads
    read (chk_unit, err=128) (wannier_spreads(i), i=1, num_wann)

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
  subroutine param_chkpt_dist
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

    integer :: ierr, loop_kpt, m, i, j

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

    call comms_bcast(have_disentangled, 1)

    if (have_disentangled) then
      if (.not. on_root) then

        if (.not. allocated(u_matrix_opt)) then
          allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating u_matrix_opt in param_chkpt_dist')
        endif

        if (.not. allocated(lwindow)) then
          allocate (lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating lwindow in param_chkpt_dist')
        endif

        if (.not. allocated(ndimwin)) then
          allocate (ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating ndimwin in param_chkpt_dist')
        endif

      end if

      call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts)
      call comms_bcast(lwindow(1, 1), num_bands*num_kpts)
      call comms_bcast(ndimwin(1), num_kpts)
      call comms_bcast(omega_invariant, 1)
    end if
    call comms_bcast(wannier_centres(1, 1), 3*num_wann)
    call comms_bcast(wannier_spreads(1), num_wann)

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
  subroutine param_get_block_length(keyword, found, rows, lunits)
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
  subroutine param_get_atoms(lunits)
    !===================================!
    !                                   !
    !!   Fills the atom data block
    !                                   !
    !===================================!

    use w90_constants, only: bohr
    use w90_utility, only: utility_frac_to_cart, utility_cart_to_frac
    use w90_io, only: io_error
    implicit none

    logical, intent(in) :: lunits
    !! Do we expect a first line with the units

    real(kind=dp)     :: atoms_pos_frac_tmp(3, num_atoms)
    real(kind=dp)     :: atoms_pos_cart_tmp(3, num_atoms)
    character(len=20) :: keyword
    integer           :: in, ins, ine, loop, i, line_e, line_s, counter
    integer           :: i_temp, loop2, max_sites, ierr, ic
    logical           :: found_e, found_s, found, frac
    character(len=maxlen) :: dummy, end_st, start_st
    character(len=maxlen) :: ctemp(num_atoms)
    character(len=maxlen) :: atoms_label_tmp(num_atoms)
    logical           :: lconvert

    keyword = "atoms_cart"
    frac = .false.
    call param_get_block_length("atoms_frac", found, i_temp)
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
      do loop = 1, num_atoms
        call utility_frac_to_cart(atoms_pos_frac_tmp(:, loop), atoms_pos_cart_tmp(:, loop), real_lattice)
      end do
    else
      do loop = 1, num_atoms
        call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), atoms_pos_frac_tmp(:, loop), recip_lattice)
      end do
    end if

    ! Now we sort the data into the proper structures
    num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          num_species = num_species + 1
          ctemp(num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atoms_species_num(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in param_get_atoms')
    allocate (atoms_label(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in param_get_atoms')
    allocate (atoms_symbol(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in param_get_atoms')
    atoms_species_num(:) = 0

    do loop = 1, num_species
      atoms_label(loop) = ctemp(loop)
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          atoms_species_num(loop) = atoms_species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atoms_species_num)
    allocate (atoms_pos_frac(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in param_get_atoms')
    allocate (atoms_pos_cart(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in param_get_atoms')

    do loop = 1, num_species
      counter = 0
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          atoms_pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atoms_pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, num_species
      atoms_symbol(loop) (1:2) = atoms_label(loop) (1:2)
      ic = ichar(atoms_symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atoms_symbol(loop) (2:2) = ' '
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

    character(len=*), intent(in) :: atoms_label_tmp(num_atoms)
    !! Atom labels
    real(kind=dp), intent(in)      :: atoms_pos_cart_tmp(3, num_atoms)
    !! Atom positions

    real(kind=dp)     :: atoms_pos_frac_tmp(3, num_atoms)
    integer           :: loop2, max_sites, ierr, ic, loop, counter
    character(len=maxlen) :: ctemp(num_atoms)
    character(len=maxlen) :: tmp_string

    do loop = 1, num_atoms
      call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), &
                                atoms_pos_frac_tmp(:, loop), recip_lattice)
    enddo

    ! Now we sort the data into the proper structures
    num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          num_species = num_species + 1
          ctemp(num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atoms_species_num(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in param_lib_set_atoms')
    allocate (atoms_label(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in param_lib_set_atoms')
    allocate (atoms_symbol(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in param_lib_set_atoms')
    atoms_species_num(:) = 0

    do loop = 1, num_species
      atoms_label(loop) = ctemp(loop)
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          atoms_species_num(loop) = atoms_species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atoms_species_num)
    allocate (atoms_pos_frac(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in param_lib_set_atoms')
    allocate (atoms_pos_cart(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in param_lib_set_atoms')

    do loop = 1, num_species
      counter = 0
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          atoms_pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atoms_pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, num_species
      atoms_symbol(loop) (1:2) = atoms_label(loop) (1:2)
      ic = ichar(atoms_symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atoms_symbol(loop) (2:2) = ' '
      tmp_string = trim(adjustl(utility_lowercase(atoms_symbol(loop))))
      atoms_symbol(loop) (1:2) = tmp_string(1:2)
      tmp_string = trim(adjustl(utility_lowercase(atoms_label(loop))))
      atoms_label(loop) (1:2) = tmp_string(1:2)
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

  subroutine param_get_centre_constraints
    !=============================================================================!
    !                                                                             !
    !!  assigns projection centres as default centre constraints and global
    !!  Lagrange multiplier as individual Lagrange multipliers then reads
    !!  the centre_constraints block for individual centre constraint parameters
    !                                                                             !
    !=============================================================================!
    use w90_io, only: io_error
    use w90_utility, only: utility_frac_to_cart
    integer           :: loop1, index1, constraint_num, index2, loop2
    integer           :: column, start, finish, wann, ierr
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
              call param_get_centre_constraint_from_column(column, start, finish, wann, dummy)
              start = loop2 + 1
              finish = start
            end if
          end if
          if (loop2 == len_trim(dummy) .and. dummy(loop2:loop2) /= ' ') then
            finish = loop2
            call param_get_centre_constraint_from_column(column, start, finish, wann, dummy)
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
      call utility_frac_to_cart(ccentres_frac(loop1, :), ccentres_cart(loop1, :), real_lattice)
    end do
  end subroutine param_get_centre_constraints

  subroutine param_get_centre_constraint_from_column(column, start, finish, wann, dummy)
    !===================================!
    !                                   !
    !!  assigns value read to constraint
    !!  parameters based on column
    !                                   !
    !===================================!
    use w90_io, only: io_error
    integer, intent(inout):: column, start, finish, wann
    character(len=maxlen), intent(inout):: dummy
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
  subroutine param_get_projections(num_proj, lcount)
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
      allocate (input_proj_site(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_site in param_get_projections')
      allocate (input_proj_l(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_l in param_get_projections')
      allocate (input_proj_m(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_m in param_get_projections')
      allocate (input_proj_z(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_z in param_get_projections')
      allocate (input_proj_x(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_x in param_get_projections')
      allocate (input_proj_radial(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_radial in param_get_projections')
      allocate (input_proj_zona(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_zona in param_get_projections')
      if (spinors) then
        allocate (input_proj_s(num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s in param_get_projections')
        allocate (input_proj_s_qaxis(3, num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s_qaxis in param_get_projections')
      endif

      allocate (proj_site(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_site in param_get_projections')
      allocate (proj_l(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_l in param_get_projections')
      allocate (proj_m(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_m in param_get_projections')
      allocate (proj_z(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_z in param_get_projections')
      allocate (proj_x(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_x in param_get_projections')
      allocate (proj_radial(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_radial in param_get_projections')
      allocate (proj_zona(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_zona in param_get_projections')
      if (spinors) then
        allocate (proj_s(num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_s in param_get_projections')
        allocate (proj_s_qaxis(3, num_wann), stat=ierr)
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
          if (num_species == 0) &
            call io_error('param_get_projection: Atom centred projection requested but no atoms defined')
          do loop = 1, num_species
            if (trim(ctemp) == atoms_label(loop)) then
              species = loop
              sites = atoms_species_num(loop)
              exit
            end if
            if (loop == num_species) call io_error('param_get_projection: Atom site not recognised '//trim(ctemp))
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
        if (spinors) then
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
                  input_proj_site(:, counter) = pos_frac
                  input_proj_l(counter) = loop_l
                  input_proj_m(counter) = loop_m
                  input_proj_z(:, counter) = proj_z_tmp
                  input_proj_x(:, counter) = proj_x_tmp
                  input_proj_radial(counter) = proj_radial_tmp
                  input_proj_zona(counter) = proj_zona_tmp
                  if (spinors) then
                    if (spn_counter == 1) then
                      if (proj_u_tmp) input_proj_s(counter) = 1
                      if (proj_d_tmp) input_proj_s(counter) = -1
                    else
                      if (loop_s == 1) input_proj_s(counter) = 1
                      if (loop_s == 2) input_proj_s(counter) = -1
                    endif
                    input_proj_s_qaxis(:, counter) = proj_s_qaxis_tmp
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
                    input_proj_site(:, counter) = atoms_pos_frac(:, loop_sites, species)
                    input_proj_l(counter) = loop_l
                    input_proj_m(counter) = loop_m
                    input_proj_z(:, counter) = proj_z_tmp
                    input_proj_x(:, counter) = proj_x_tmp
                    input_proj_radial(counter) = proj_radial_tmp
                    input_proj_zona(counter) = proj_zona_tmp
                    if (spinors) then
                      if (spn_counter == 1) then
                        if (proj_u_tmp) input_proj_s(counter) = 1
                        if (proj_d_tmp) input_proj_s(counter) = -1
                      else
                        if (loop_s == 1) input_proj_s(counter) = 1
                        if (loop_s == 2) input_proj_s(counter) = -1
                      endif
                      input_proj_s_qaxis(:, counter) = proj_s_qaxis_tmp
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
        call random_number(input_proj_site(:, loop))
        input_proj_l(loop) = 0
        input_proj_m(loop) = 1
        input_proj_z(:, loop) = proj_z_def
        input_proj_x(:, loop) = proj_x_def
        input_proj_zona(loop) = proj_zona_def
        input_proj_radial(loop) = proj_radial_def
        if (spinors) then
          if (modulo(loop, 2) == 1) then
            input_proj_s(loop) = 1
          else
            input_proj_s(loop) = -1
          end if
          input_proj_s_qaxis(1, loop) = 0.
          input_proj_s_qaxis(2, loop) = 0.
          input_proj_s_qaxis(3, loop) = 1.
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

      znorm = sqrt(sum(input_proj_z(:, loop)*input_proj_z(:, loop)))
      xnorm = sqrt(sum(input_proj_x(:, loop)*input_proj_x(:, loop)))
      input_proj_z(:, loop) = input_proj_z(:, loop)/znorm             ! normalise z
      input_proj_x(:, loop) = input_proj_x(:, loop)/xnorm             ! normalise x
      cosphi = sum(input_proj_z(:, loop)*input_proj_x(:, loop))

      ! Check whether z-axis and z-axis are orthogonal
      if (abs(cosphi) .gt. eps6) then

        ! Special case of circularly symmetric projections (pz, dz2, fz3)
        ! just choose an x-axis that is perpendicular to the given z-axis
        if ((input_proj_l(loop) .ge. 0) .and. (input_proj_m(loop) .eq. 1)) then
          proj_x_tmp(:) = input_proj_x(:, loop)            ! copy of original x-axis
          call random_seed()
          call random_number(proj_z_tmp(:))         ! random vector
          ! calculate new x-axis as the cross (vector) product of random vector with z-axis
          input_proj_x(1, loop) = proj_z_tmp(2)*input_proj_z(3, loop) - proj_z_tmp(3)*input_proj_z(2, loop)
          input_proj_x(2, loop) = proj_z_tmp(3)*input_proj_z(1, loop) - proj_z_tmp(1)*input_proj_z(3, loop)
          input_proj_x(3, loop) = proj_z_tmp(1)*input_proj_z(2, loop) - proj_z_tmp(2)*input_proj_z(1, loop)
          xnorm_new = sqrt(sum(input_proj_x(:, loop)*input_proj_x(:, loop)))
          input_proj_x(:, loop) = input_proj_x(:, loop)/xnorm_new   ! normalise
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
        proj_x_tmp(:) = input_proj_x(:, loop)               ! copy of original x-axis
        ! calculate new x-axis:
        ! x = z \cross (x_tmp \cross z) / sinphi = ( x_tmp - z(z.x_tmp) ) / sinphi
        input_proj_x(:, loop) = (proj_x_tmp(:) - cosphi*input_proj_z(:, loop))/sinphi

        ! Final check
555     cosphi_new = sum(input_proj_z(:, loop)*input_proj_x(:, loop))
        if (abs(cosphi_new) .gt. eps6) then
          write (stdout, *) ' Projection:', loop
          call io_error(' Error: z and x axes are still not orthogonal after projection')
        endif

      endif

    enddo

    do loop = 1, num_proj
      if (proj2wann_map(loop) < 0) cycle
      proj_site(:, proj2wann_map(loop)) = input_proj_site(:, loop)
      proj_l(proj2wann_map(loop)) = input_proj_l(loop)
      proj_m(proj2wann_map(loop)) = input_proj_m(loop)
      proj_z(:, proj2wann_map(loop)) = input_proj_z(:, loop)
      proj_x(:, proj2wann_map(loop)) = input_proj_x(:, loop)
      proj_radial(proj2wann_map(loop)) = input_proj_radial(loop)
      proj_zona(proj2wann_map(loop)) = input_proj_zona(loop)
    enddo

    if (spinors) then
      do loop = 1, num_proj
        if (proj2wann_map(loop) < 0) cycle
        proj_s(proj2wann_map(loop)) = input_proj_s(loop)
        proj_s_qaxis(:, proj2wann_map(loop)) = input_proj_s_qaxis(:, loop)
      enddo
    endif

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
      read (dummy, *, err=240, end=240) bands_label(counter - 1), (bands_spec_points(i, counter - 1), i=1, 3) &
        , bands_label(counter), (bands_spec_points(i, counter), i=1, 3)
    end do

    in_data(line_s:line_e) (1:maxlen) = ' '

    return

240 call io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy))

  end subroutine param_get_keyword_kpath

!===========================================!
  subroutine param_memory_estimate
    !===========================================!
    !                                           !
    !! Estimate how much memory we will allocate
    !                                           !
    !===========================================!

    use w90_comms, only: on_root

    implicit none

    real(kind=dp), parameter :: size_log = 1.0_dp
    real(kind=dp), parameter :: size_int = 4.0_dp
    real(kind=dp), parameter :: size_real = 8.0_dp
    real(kind=dp), parameter :: size_cmplx = 16.0_dp
    real(kind=dp) :: mem_wan, mem_wan1, mem_param, mem_dis, mem_dis2, mem_dis1
    real(kind=dp) :: mem_bw
    integer :: NumPoints1, NumPoints2, NumPoints3, ndim
    real(kind=dp) :: TDF_exceeding_energy

    mem_param = 0
    mem_dis = 0
    mem_dis1 = 0
    mem_dis2 = 0
    mem_wan = 0
    mem_wan1 = 0
    mem_bw = 0

    ! First the data stored in the parameters module
    mem_param = mem_param + num_wann*num_wann*num_kpts*size_cmplx                   !u_matrix
    if (.not. disentanglement) &
      mem_param = mem_param + num_wann*num_wann*nntot*num_kpts*size_cmplx       !m_matrix

    if (disentanglement) then
      mem_param = mem_param + num_bands*num_wann*num_kpts*size_cmplx             ! u_matrix_opt
    endif

    if (allocated(atoms_species_num)) then
      mem_param = mem_param + (num_species)*size_int                               !atoms_species_num
      mem_param = mem_param + (num_species)*size_real                              !atoms_label
      mem_param = mem_param + (num_species)*size_real                              !atoms_symbol
      mem_param = mem_param + (3*maxval(atoms_species_num)*num_species)*size_real  !atoms_pos_frac
      mem_param = mem_param + (3*maxval(atoms_species_num)*num_species)*size_real  !atoms_pos_cart
    endif

    if (allocated(input_proj_site)) then
      mem_param = mem_param + (3*num_proj)*size_real              !input_proj_site
      mem_param = mem_param + (num_proj)*size_int                !input_proj_l
      mem_param = mem_param + (num_proj)*size_int                 !input_proj_m
      mem_param = mem_param + (3*num_proj)*size_real             !input_proj_z
      mem_param = mem_param + (3*num_proj)*size_real             !input_proj_x
      mem_param = mem_param + (num_proj)*size_real                !input_proj_radial
      mem_param = mem_param + (num_proj)*size_real                !input_proj_zona
    endif

    if (allocated(proj_site)) then
      mem_param = mem_param + (3*num_wann)*size_real              !proj_site
      mem_param = mem_param + (num_wann)*size_int                !proj_l
      mem_param = mem_param + (num_wann)*size_int                 !proj_m
      mem_param = mem_param + (3*num_wann)*size_real             !proj_z
      mem_param = mem_param + (3*num_wann)*size_real             !proj_x
      mem_param = mem_param + (num_wann)*size_real                !proj_radial
      mem_param = mem_param + (num_wann)*size_real                !proj_zona
    endif

    mem_param = mem_param + num_kpts*nntot*size_int                  !nnlist
    mem_param = mem_param + num_kpts*nntot/2*size_int                !neigh
    mem_param = mem_param + 3*num_kpts*nntot*size_int                !nncell
    mem_param = mem_param + nntot*size_real                          !wb
    mem_param = mem_param + 3*nntot/2*size_real                      !bka
    mem_param = mem_param + 3*nntot*num_kpts*size_real               !bk

    mem_param = mem_param + num_bands*num_kpts*size_real             !eigval
    mem_param = mem_param + 3*num_kpts*size_real                     !kpt_cart
    mem_param = mem_param + 3*num_kpts*size_real                     !kpt_latt
    if (disentanglement) then
      mem_param = mem_param + num_kpts*size_int                     !ndimwin
      mem_param = mem_param + num_bands*num_kpts*size_log           !lwindow
    endif
    mem_param = mem_param + 3*num_wann*size_real                     !wannier_centres
    mem_param = mem_param + num_wann*size_real                       !wannier_spreads

    if (disentanglement) then
      ! Module vars
      mem_dis = mem_dis + num_bands*num_kpts*size_real              !eigval_opt
      mem_dis = mem_dis + num_kpts*size_int                         !nfirstwin
      mem_dis = mem_dis + num_kpts*size_int                         !ndimfroz
      mem_dis = mem_dis + num_bands*num_kpts*size_int               !indxfroz
      mem_dis = mem_dis + num_bands*num_kpts*size_int               !indxnfroz
      mem_dis = mem_dis + num_bands*num_kpts*size_log               !lfrozen

      !the memory high-water wiil occur in dis_extract or when we allocate m_matrix

      mem_dis1 = mem_dis1 + num_wann*num_bands*size_cmplx              !cwb
      mem_dis1 = mem_dis1 + num_wann*num_wann*size_cmplx               !cww
      mem_dis1 = mem_dis1 + num_bands*num_wann*size_cmplx              !cbw
      mem_dis1 = mem_dis1 + 5*num_bands*size_int                       !iwork
      mem_dis1 = mem_dis1 + num_bands*size_int                         !ifail
      mem_dis1 = mem_dis1 + num_bands*size_real                        !w
      if (gamma_only) then
        mem_dis1 = mem_dis1 + (num_bands*(num_bands + 1))/2*size_real    !cap_r
        mem_dis1 = mem_dis1 + 8*num_bands*size_real                    !work
        mem_dis1 = mem_dis1 + num_bands*num_bands*size_real            !rz
      else
        mem_dis1 = mem_dis1 + 7*num_bands*size_real                    !rwork
        mem_dis1 = mem_dis1 + (num_bands*(num_bands + 1))/2*size_cmplx   !cap
        mem_dis1 = mem_dis1 + 2*num_bands*size_cmplx                   !cwork
        mem_dis1 = mem_dis1 + num_bands*num_bands*size_cmplx           !cz
      end if
      mem_dis1 = mem_dis1 + num_kpts*size_real                         !wkomegai1
      mem_dis1 = mem_dis1 + num_bands*num_bands*num_kpts*size_cmplx    !ceamp
      mem_dis1 = mem_dis1 + num_bands*num_bands*num_kpts*size_cmplx    !cham
      mem_dis2 = mem_dis2 + num_wann*num_wann*nntot*num_kpts*size_cmplx!m_matrix

      if (optimisation <= 0) then
        mem_dis = mem_dis + mem_dis1
      else
        mem_dis = mem_dis + max(mem_dis1, mem_dis2)
      endif

      mem_dis = mem_dis + num_bands*num_bands*nntot*num_kpts*size_cmplx      ! m_matrix_orig
      mem_dis = mem_dis + num_bands*num_wann*num_kpts*size_cmplx             ! a_matrix

    endif

    !Wannierise

    mem_wan1 = mem_wan1 + (num_wann*num_wann*nntot*num_kpts)*size_cmplx     !  'm0'
    if (optimisation > 0) then
      mem_wan = mem_wan + mem_wan1
    endif
    mem_wan = mem_wan + (num_wann*num_wann*num_kpts)*size_cmplx           !  'u0'
    mem_wan = mem_wan + (num_wann*nntot*num_kpts)*size_real               !  'rnkb'
    mem_wan = mem_wan + (num_wann*nntot*num_kpts)*size_real               !  'ln_tmp'
    mem_wan = mem_wan + (num_wann*nntot*num_kpts)*size_cmplx              !  'csheet'
    mem_wan = mem_wan + (num_wann*nntot*num_kpts)*size_real               !  'sheet'
    mem_wan = mem_wan + (3*num_wann)*size_real                             !  'rave'
    mem_wan = mem_wan + (num_wann)*size_real                              !  'r2ave'
    mem_wan = mem_wan + (num_wann)*size_real                               !  'rave2'
    mem_wan = mem_wan + (3*num_wann)*size_real                            !  'rguide'
    mem_wan = mem_wan + (num_wann*num_wann)*size_cmplx                  !  'cz'
    if (gamma_only) then
      mem_wan = mem_wan + num_wann*num_wann*nntot*2*size_cmplx    ! m_w
      mem_wan = mem_wan + num_wann*num_wann*size_cmplx            ! uc_rot
      mem_wan = mem_wan + num_wann*num_wann*size_real             ! ur_rot
      !internal_svd_omega_i
      mem_wan = mem_wan + 10*num_wann*size_cmplx                   ! cw1
      mem_wan = mem_wan + 10*num_wann*size_cmplx                   ! cw2
      mem_wan = mem_wan + num_wann*num_wann*size_cmplx             ! cv1
      mem_wan = mem_wan + num_wann*num_wann*size_cmplx             ! cv2
      mem_wan = mem_wan + num_wann*num_wann*size_real              ! cpad1
      mem_wan = mem_wan + num_wann*size_cmplx                      ! singvd
    else
      mem_wan = mem_wan + (num_wann)*size_cmplx                              !  'cwschur1'
      mem_wan = mem_wan + (10*num_wann)*size_cmplx                        !  'cwschur2'
      mem_wan = mem_wan + (num_wann)*size_cmplx                              !  'cwschur3'
      mem_wan = mem_wan + (num_wann)*size_cmplx                             !  'cwschur4'
      mem_wan = mem_wan + (num_wann*num_wann*num_kpts)*size_cmplx         !  'cdq'
      mem_wan = mem_wan + (num_wann*num_wann)*size_cmplx                   !  'cmtmp'
      mem_wan = mem_wan + (num_wann*num_wann*num_kpts)*size_cmplx         !  'cdqkeep'
      mem_wan = mem_wan + (num_wann*num_wann)*size_cmplx                      !  'tmp_cdq'
      mem_wan = mem_wan + (num_wann)*size_real                               !  'evals'
      mem_wan = mem_wan + (4*num_wann)*size_cmplx                            !  'cwork'
      mem_wan = mem_wan + (3*num_wann - 2)*size_real                           !  'rwork'
      !d_omega
      mem_wan = mem_wan + (num_wann*num_wann)*size_cmplx   !  'cr'
      mem_wan = mem_wan + (num_wann*num_wann)*size_cmplx    !  'crt'
    end if

    if (ispostw90) then
      if (boltzwann) then
        if (spin_decomp) then
          ndim = 3
        else
          ndim = 1
        end if

        ! I set a big value to have a rough estimate
        TDF_exceeding_energy = 2._dp
        NumPoints1 = int(floor((boltz_temp_max - boltz_temp_min)/boltz_temp_step)) + 1 ! temperature array
        NumPoints2 = int(floor((boltz_mu_max - boltz_mu_min)/boltz_mu_step)) + 1  ! mu array
        NumPoints3 = int(floor((dis_win_max - dis_win_min + 2._dp*TDF_exceeding_energy)/boltz_tdf_energy_step)) + 1 ! tdfenergyarray
        mem_bw = mem_bw + NumPoints1*size_real                         !TempArray
        mem_bw = mem_bw + NumPoints1*size_real                         !KTArray
        mem_bw = mem_bw + NumPoints2*size_real                         !MuArray
        mem_bw = mem_bw + NumPoints3*size_real                         !TDFEnergyArray
        mem_bw = mem_bw + 6*NumPoints3*ndim*size_real                  !TDFArray
        mem_bw = mem_bw + 6*NumPoints3*size_real                       !IntegrandArray
        mem_bw = mem_bw + (9*4 + 6)*size_real
        !ElCondTimesSeebeckFP,ThisElCond,ElCondInverse,ThisSeebeck,ElCondTimesSeebeck
        mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !ElCond
        mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !Seebeck
        mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !ThermCond
        ! I put a upper bound here below (as if there was only 1 node), because I do not have any knowledge at this point
        ! of the number of processors, so I cannot have a correct estimate
        mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !LocalElCond
        mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !LocalSeebeck
        mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !LocalThermCond

        mem_bw = mem_bw + num_wann*num_wann*size_cmplx                 !HH
        mem_bw = mem_bw + 3*num_wann*num_wann*size_cmplx               !delHH
        mem_bw = mem_bw + num_wann*num_wann*size_cmplx                 !UU
        mem_bw = mem_bw + 3*num_wann*size_real                         !del_eig
        mem_bw = mem_bw + num_wann*size_real                           !eig
        mem_bw = mem_bw + num_wann*size_real                           !levelspacing_k

        NumPoints1 = int(floor((boltz_dos_energy_max - boltz_dos_energy_min)/boltz_dos_energy_step)) + 1!dosnumpoints
        mem_bw = mem_bw + NumPoints1*size_real                         !DOS_EnergyArray
        mem_bw = mem_bw + 6*ndim*NumPoints3*size_real                  !TDF_k
        mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_k
        mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_all
      end if
    end if

    if (disentanglement) &
      mem_wan = mem_wan + num_wann*num_wann*nntot*num_kpts*size_cmplx       !m_matrix

    if (on_root) then
      write (stdout, '(1x,a)') '*============================================================================*'
      write (stdout, '(1x,a)') '|                              MEMORY ESTIMATE                               |'
      write (stdout, '(1x,a)') '|         Maximum RAM allocated during each phase of the calculation         |'
      write (stdout, '(1x,a)') '*============================================================================*'
      if (disentanglement) &
        write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Disentanglement:', (mem_param + mem_dis)/(1024**2), ' Mb'
      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:', (mem_param + mem_wan)/(1024**2), ' Mb'
      if (optimisation > 0 .and. iprint > 1) then
        write (stdout, '(1x,a)') '|                                                                            |'
        write (stdout, '(1x,a)') '|   N.B. by setting optimisation=0 memory usage will be reduced to:          |'
        if (disentanglement) &
          write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Disentanglement:', &
          (mem_param + mem_dis - max(mem_dis1, mem_dis2) + mem_dis1)/(1024**2), ' Mb'
        if (gamma_only) then
          write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:', (mem_param + mem_wan)/(1024**2), ' Mb'
        else
          write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:', &
            (mem_param + mem_wan - mem_wan1)/(1024**2), ' Mb'
        end if
        write (stdout, '(1x,a)') '|   However, this will result in more i/o and slow down the calculation      |'
      endif

      if (ispostw90) then
        if (boltzwann) &
          write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'BoltzWann:', (mem_param + mem_bw)/(1024**2), ' Mb'
      end if

      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'plot_wannier:', (mem_param + mem_wan)/(1024**2), ' Mb'
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, *) ' '
    endif

!    if(disentanglement) then
!       write(*,'(a12,f12.4,a)') 'Disentangle',(mem_param+mem_dis)/(1024**2),' Mb'
!    end if
!    write(*,'(a12,f12.4,a)') 'Wannierise ',(mem_wan+mem_param)/(1024**2),' Mb'
!    write(*,'(a12,f12.4,a)') 'Module',(mem_param)/(1024**2),' Mb'

    return
  end subroutine param_memory_estimate

!===========================================================!
  subroutine param_dist
    !===========================================================!
    !                                                           !
    !! distribute the parameters across processors              !
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, io_date, io_time, &
      io_stopwatch
    use w90_comms, only: comms_bcast, on_root

    integer :: ierr

    call comms_bcast(effective_model, 1)
    call comms_bcast(eig_found, 1)
    call comms_bcast(postproc_setup, 1)
    if (.not. effective_model) then
      call comms_bcast(mp_grid(1), 3)
      call comms_bcast(num_kpts, 1)
      call comms_bcast(num_bands, 1)
    endif
    call comms_bcast(num_wann, 1)
    call comms_bcast(timing_level, 1)
    call comms_bcast(iprint, 1)
    call comms_bcast(energy_unit, 1)
    call comms_bcast(length_unit, 1)
    call comms_bcast(wvfn_formatted, 1)
    call comms_bcast(spn_formatted, 1)
    call comms_bcast(uHu_formatted, 1)
    call comms_bcast(berry_uHu_formatted, 1)
    call comms_bcast(spin, 1)
    call comms_bcast(num_dump_cycles, 1)
    call comms_bcast(num_print_cycles, 1)
    call comms_bcast(num_atoms, 1)   ! Ivo: not used in postw90, right?
    call comms_bcast(num_species, 1) ! Ivo: not used in postw90, right?
    call comms_bcast(real_lattice(1, 1), 9)
    call comms_bcast(recip_lattice(1, 1), 9)
    call comms_bcast(real_metric(1, 1), 9)
    call comms_bcast(recip_metric(1, 1), 9)
    call comms_bcast(cell_volume, 1)
    call comms_bcast(dos_energy_step, 1)
    call comms_bcast(dos_adpt_smr, 1)
    call comms_bcast(dos_smr_index, 1)
    call comms_bcast(dos_kmesh_spacing, 1)
    call comms_bcast(dos_kmesh(1), 3)
    call comms_bcast(dos_adpt_smr_max, 1)
    call comms_bcast(dos_smr_fixed_en_width, 1)
    call comms_bcast(dos_adpt_smr_fac, 1)
    call comms_bcast(num_dos_project, 1)
    call comms_bcast(num_exclude_bands, 1)
    if (num_exclude_bands > 0) then
      if (.not. on_root) then
        allocate (exclude_bands(num_exclude_bands), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating exclude_bands in param_dist')
      endif
      call comms_bcast(exclude_bands(1), num_exclude_bands)
    end if

    call comms_bcast(gamma_only, 1)
    call comms_bcast(dis_win_min, 1)
    call comms_bcast(dis_win_max, 1)
    call comms_bcast(dis_froz_min, 1)
    call comms_bcast(dis_froz_max, 1)
    call comms_bcast(dis_num_iter, 1)
    call comms_bcast(dis_mix_ratio, 1)
    call comms_bcast(dis_conv_tol, 1)
    call comms_bcast(dis_conv_window, 1)
    call comms_bcast(dis_spheres_first_wann, 1)
    call comms_bcast(dis_spheres_num, 1)
    if (dis_spheres_num > 0) then
      if (.not. on_root) then
        allocate (dis_spheres(4, dis_spheres_num), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating dis_spheres in param_dist')
      endif
      call comms_bcast(dis_spheres(1, 1), 4*dis_spheres_num)
    end if
    call comms_bcast(num_iter, 1)
    call comms_bcast(num_cg_steps, 1)
    call comms_bcast(conv_tol, 1)
    call comms_bcast(conv_window, 1)
    call comms_bcast(guiding_centres, 1)
    call comms_bcast(wannier_plot, 1)
    call comms_bcast(num_wannier_plot, 1)
    if (num_wannier_plot > 0) then
      if (.not. on_root) then
        allocate (wannier_plot_list(num_wannier_plot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating wannier_plot_list in param_dist')
      endif
      call comms_bcast(wannier_plot_list(1), num_wannier_plot)
    end if
    call comms_bcast(wannier_plot_supercell(1), 3)
    call comms_bcast(wannier_plot_format, len(wannier_plot_format))
    call comms_bcast(wannier_plot_mode, len(wannier_plot_mode))
    call comms_bcast(wannier_plot_spinor_mode, len(wannier_plot_spinor_mode))
    call comms_bcast(write_u_matrices, 1)
    call comms_bcast(bands_plot, 1)
    call comms_bcast(write_bvec, 1)
    call comms_bcast(bands_num_points, 1)
    call comms_bcast(bands_plot_format, len(bands_plot_format))
    call comms_bcast(bands_plot_mode, len(bands_plot_mode))
    call comms_bcast(num_bands_project, 1)

    if (num_bands_project > 0) then
      if (.not. on_root) then
        allocate (bands_plot_project(num_bands_project), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bands_plot_project in param_dist')
      endif
      call comms_bcast(bands_plot_project(1), num_bands_project)
    end if
    call comms_bcast(bands_plot_dim, 1)
    call comms_bcast(write_hr, 1)
    call comms_bcast(write_rmn, 1)
    call comms_bcast(write_tb, 1)
    call comms_bcast(hr_cutoff, 1)
    call comms_bcast(dist_cutoff, 1)
    call comms_bcast(dist_cutoff_mode, len(dist_cutoff_mode))
    call comms_bcast(dist_cutoff_hc, 1)
    call comms_bcast(one_dim_axis, len(one_dim_axis))
    call comms_bcast(use_ws_distance, 1)
    call comms_bcast(ws_distance_tol, 1)
    call comms_bcast(ws_search_size(1), 3)
    call comms_bcast(fermi_surface_plot, 1)
    call comms_bcast(fermi_surface_num_points, 1)
    call comms_bcast(fermi_surface_plot_format, len(fermi_surface_plot_format))
    call comms_bcast(fermi_energy, 1) !! used?

    call comms_bcast(berry, 1)
    call comms_bcast(berry_task, len(berry_task))
    call comms_bcast(berry_kmesh_spacing, 1)
    call comms_bcast(berry_kmesh(1), 3)
    call comms_bcast(berry_curv_adpt_kmesh, 1)
    call comms_bcast(berry_curv_adpt_kmesh_thresh, 1)
    call comms_bcast(berry_curv_unit, len(berry_curv_unit))
!  Stepan Tsirkin
    call comms_bcast(gyrotropic, 1)
    call comms_bcast(gyrotropic_task, len(gyrotropic_task))
    call comms_bcast(gyrotropic_kmesh_spacing, 1)
    call comms_bcast(gyrotropic_kmesh(1), 3)
    call comms_bcast(gyrotropic_smr_fixed_en_width, 1)
    call comms_bcast(gyrotropic_smr_index, 1)
    call comms_bcast(gyrotropic_eigval_max, 1)
    call comms_bcast(gyrotropic_nfreq, 1)
    call comms_bcast(gyrotropic_degen_thresh, 1)
    call comms_bcast(gyrotropic_num_bands, 1)
    call comms_bcast(gyrotropic_box(1, 1), 9)
    call comms_bcast(gyrotropic_box_corner(1), 3)
    call comms_bcast(gyrotropic_smr_max_arg, 1)
    call comms_bcast(gyrotropic_smr_fixed_en_width, 1)
    call comms_bcast(gyrotropic_smr_index, 1)

    call comms_bcast(kubo_adpt_smr, 1)
    call comms_bcast(kubo_adpt_smr_fac, 1)
    call comms_bcast(kubo_adpt_smr_max, 1)
    call comms_bcast(kubo_smr_fixed_en_width, 1)
    call comms_bcast(kubo_smr_index, 1)
    call comms_bcast(kubo_eigval_max, 1)
    call comms_bcast(kubo_nfreq, 1)
    call comms_bcast(nfermi, 1)
    call comms_bcast(dos_energy_min, 1)
    call comms_bcast(dos_energy_max, 1)
    call comms_bcast(spin_kmesh_spacing, 1)
    call comms_bcast(spin_kmesh(1), 3)
    call comms_bcast(wanint_kpoint_file, 1)
! Junfeng Qiao
    call comms_bcast(shc_freq_scan, 1)
    call comms_bcast(shc_alpha, 1)
    call comms_bcast(shc_beta, 1)
    call comms_bcast(shc_gamma, 1)
    call comms_bcast(shc_bandshift, 1)
    call comms_bcast(shc_bandshift_firstband, 1)
    call comms_bcast(shc_bandshift_energyshift, 1)

    call comms_bcast(devel_flag, len(devel_flag))
    call comms_bcast(spin_moment, 1)
    call comms_bcast(spin_axis_polar, 1)
    call comms_bcast(spin_axis_azimuth, 1)
    call comms_bcast(spin_decomp, 1)
    call comms_bcast(use_degen_pert, 1)
    call comms_bcast(degen_thr, 1)
    call comms_bcast(num_valence_bands, 1)
    call comms_bcast(dos, 1)
    call comms_bcast(dos_task, len(dos_task))
    call comms_bcast(kpath, 1)
    call comms_bcast(kpath_task, len(kpath_task))
    call comms_bcast(kpath_bands_colour, len(kpath_bands_colour))
    call comms_bcast(kslice, 1)
    call comms_bcast(kslice_task, len(kslice_task))
    call comms_bcast(transl_inv, 1)
    call comms_bcast(num_elec_per_state, 1)
    call comms_bcast(scissors_shift, 1)
    !

! ----------------------------------------------
    call comms_bcast(geninterp, 1)
    call comms_bcast(geninterp_alsofirstder, 1)
    call comms_bcast(geninterp_single_file, 1)
    ! [gp-begin, Apr 12, 2012]
    ! BoltzWann variables
    call comms_bcast(boltzwann, 1)
    call comms_bcast(boltz_calc_also_dos, 1)
    call comms_bcast(boltz_2d_dir_num, 1)
    call comms_bcast(boltz_dos_energy_step, 1)
    call comms_bcast(boltz_dos_energy_min, 1)
    call comms_bcast(boltz_dos_energy_max, 1)
    call comms_bcast(boltz_dos_adpt_smr, 1)
    call comms_bcast(boltz_dos_smr_fixed_en_width, 1)
    call comms_bcast(boltz_dos_adpt_smr_fac, 1)
    call comms_bcast(boltz_dos_adpt_smr_max, 1)
    call comms_bcast(boltz_mu_min, 1)
    call comms_bcast(boltz_mu_max, 1)
    call comms_bcast(boltz_mu_step, 1)
    call comms_bcast(boltz_temp_min, 1)
    call comms_bcast(boltz_temp_max, 1)
    call comms_bcast(boltz_temp_step, 1)
    call comms_bcast(boltz_kmesh_spacing, 1)
    call comms_bcast(boltz_kmesh(1), 3)
    call comms_bcast(boltz_tdf_energy_step, 1)
    call comms_bcast(boltz_relax_time, 1)
    call comms_bcast(boltz_TDF_smr_fixed_en_width, 1)
    call comms_bcast(boltz_TDF_smr_index, 1)
    call comms_bcast(boltz_dos_smr_index, 1)
    call comms_bcast(boltz_bandshift, 1)
    call comms_bcast(boltz_bandshift_firstband, 1)
    call comms_bcast(boltz_bandshift_energyshift, 1)
    ! [gp-end]
    call comms_bcast(use_ws_distance, 1)
    call comms_bcast(disentanglement, 1)

    call comms_bcast(transport, 1)
    call comms_bcast(tran_easy_fix, 1)
    call comms_bcast(transport_mode, len(transport_mode))
    call comms_bcast(tran_win_min, 1)
    call comms_bcast(tran_win_max, 1)
    call comms_bcast(tran_energy_step, 1)
    call comms_bcast(tran_num_bb, 1)
    call comms_bcast(tran_num_ll, 1)
    call comms_bcast(tran_num_rr, 1)
    call comms_bcast(tran_num_cc, 1)
    call comms_bcast(tran_num_lc, 1)
    call comms_bcast(tran_num_cr, 1)
    call comms_bcast(tran_num_bandc, 1)
    call comms_bcast(tran_write_ht, 1)
    call comms_bcast(tran_read_ht, 1)
    call comms_bcast(tran_use_same_lead, 1)
    call comms_bcast(tran_num_cell_ll, 1)
    call comms_bcast(tran_num_cell_rr, 1)
    call comms_bcast(tran_group_threshold, 1)
    call comms_bcast(translation_centre_frac(1), 3)
    call comms_bcast(num_shells, 1)
    call comms_bcast(skip_B1_tests, 1)
    call comms_bcast(explicit_nnkpts, 1)

    call comms_bcast(calc_only_A, 1)
    call comms_bcast(use_bloch_phases, 1)
    call comms_bcast(restart, len(restart))
    call comms_bcast(write_r2mn, 1)
    call comms_bcast(num_guide_cycles, 1)
    call comms_bcast(num_no_guide_iter, 1)
    call comms_bcast(fixed_step, 1)
    call comms_bcast(trial_step, 1)
    call comms_bcast(precond, 1)
    call comms_bcast(write_proj, 1)
    call comms_bcast(timing_level, 1)
    call comms_bcast(spinors, 1)
    call comms_bcast(num_elec_per_state, 1)
    call comms_bcast(translate_home_cell, 1)
    call comms_bcast(write_xyz, 1)
    call comms_bcast(write_hr_diag, 1)
    call comms_bcast(conv_noise_amp, 1)
    call comms_bcast(conv_noise_num, 1)
    call comms_bcast(wannier_plot_radius, 1)
    call comms_bcast(wannier_plot_scale, 1)
    call comms_bcast(kmesh_tol, 1)
    call comms_bcast(optimisation, 1)
    call comms_bcast(write_vdw_data, 1)
    call comms_bcast(lenconfac, 1)
    call comms_bcast(lfixstep, 1)
    call comms_bcast(lsitesymmetry, 1)
    call comms_bcast(frozen_states, 1)

    !vv: Constrained centres
    call comms_bcast(slwf_num, 1)
    call comms_bcast(slwf_constrain, 1)
    call comms_bcast(slwf_lambda, 1)
    call comms_bcast(selective_loc, 1)
    if (selective_loc .and. slwf_constrain) then
      if (.not. on_root) then
        allocate (ccentres_frac(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_frac in param_get_centre_constraints')
        allocate (ccentres_cart(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_cart in param_get_centre_constraints')
      endif
      call comms_bcast(ccentres_frac(1, 1), 3*num_wann)
      call comms_bcast(ccentres_cart(1, 1), 3*num_wann)
    end if

    ! vv: automatic projections
    call comms_bcast(auto_projections, 1)

    call comms_bcast(num_proj, 1)
    call comms_bcast(lhasproj, 1)
    if (lhasproj) then
      if (.not. on_root) then
        allocate (input_proj_site(3, num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_site in param_dist')
        allocate (proj_site(3, num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_site in param_dist')
      endif
      call comms_bcast(input_proj_site(1, 1), 3*num_proj)
      call comms_bcast(proj_site(1, 1), 3*num_wann)
    endif

    ! These variables are different from the ones above in that they are
    ! allocatable, and in param_read they were allocated on the root node only
    !
    if (.not. on_root) then
      allocate (fermi_energy_list(nfermi), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating fermi_energy_read in postw90_param_dist')
      allocate (kubo_freq_list(kubo_nfreq), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating kubo_freq_list in postw90_param_dist')
      allocate (dos_project(num_dos_project), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating dos_project in postw90_param_dist')
      if (.not. effective_model) then
        if (eig_found) then
          allocate (eigval(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating eigval in postw90_param_dist')
        end if
        allocate (kpt_latt(3, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error allocating kpt_latt in postw90_param_dist')
      endif
      allocate (gyrotropic_band_list(gyrotropic_num_bands), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating gyrotropic_num_bands in postw90_param_dist')
      allocate (gyrotropic_freq_list(gyrotropic_nfreq), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating gyrotropic_freq_list in postw90_param_dist')
    end if

    if (nfermi > 0) call comms_bcast(fermi_energy_list(1), nfermi)
    if (kubo_nfreq > 0) call comms_bcast(kubo_freq_list(1), kubo_nfreq)
    call comms_bcast(gyrotropic_freq_list(1), gyrotropic_nfreq)
    call comms_bcast(gyrotropic_band_list(1), gyrotropic_num_bands)
    if (num_dos_project > 0) call comms_bcast(dos_project(1), num_dos_project)
    if (.not. effective_model) then
      if (eig_found) then
        call comms_bcast(eigval(1, 1), num_bands*num_kpts)
      end if
      call comms_bcast(kpt_latt(1, 1), 3*num_kpts)
    endif

    if (.not. effective_model .and. .not. explicit_nnkpts) then

      call comms_bcast(nnh, 1)
      call comms_bcast(nntot, 1)
      call comms_bcast(wbtot, 1)

      if (.not. on_root) then
        allocate (nnlist(num_kpts, nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nnlist in param_dist')
        allocate (neigh(num_kpts, nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating neigh in param_dist')
        allocate (nncell(3, num_kpts, nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nncell in param_dist')
        allocate (wb(nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating wb in param_dist')
        allocate (bka(3, nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bka in param_dist')
        allocate (bk(3, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bk in param_dist')
      end if

      call comms_bcast(nnlist(1, 1), num_kpts*nntot)
      call comms_bcast(neigh(1, 1), num_kpts*nntot/2)
      call comms_bcast(nncell(1, 1, 1), 3*num_kpts*nntot)
      call comms_bcast(wb(1), nntot)
      call comms_bcast(bka(1, 1), 3*nntot/2)
      call comms_bcast(bk(1, 1, 1), 3*nntot*num_kpts)

    endif

    call comms_bcast(omega_total, 1)
    call comms_bcast(omega_tilde, 1)
    call comms_bcast(omega_invariant, 1)
    call comms_bcast(have_disentangled, 1)

    if (.not. on_root) then
      allocate (wannier_centres(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in param_dist')
      wannier_centres = 0.0_dp
      allocate (wannier_spreads(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating wannier_spreads in param_dist')
      wannier_spreads = 0.0_dp
      if (disentanglement) then
        allocate (ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ndimwin in param_dist')
        allocate (lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating lwindow in param_dist')
      endif
    endif

  end subroutine param_dist

  subroutine parameters_gyro_write_task(task, key, comment)
    use w90_io, only: stdout

    character(len=*), intent(in) :: task, key, comment
    character(len=42) :: comment1

    comment1 = comment
    if ((index(task, key) > 0) .or. (index(task, 'all') > 0)) then
      write (stdout, '(1x,a2,a42,a2,10x,a8,13x,a1)') '| ', comment1, ' :', '       T', '|'
    else
      write (stdout, '(1x,a2,a42,a2,10x,a8,13x,a1)') '| ', comment1, ' :', '       F', '|'
    endif
  end subroutine parameters_gyro_write_task

end module w90_parameters
