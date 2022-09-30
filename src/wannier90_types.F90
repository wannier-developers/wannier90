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
!  w90_wannier90_types: data types specific to wannier90.x   !
!                                                            !
!------------------------------------------------------------!

module w90_wannier90_types

  ! Definition of types encapsulating various quantities, data and parameters.
  ! Variables are grouped according to physical meaning and their use in the Wannier90 project.
  !
  !! Here are defined types specific to wannier90.x (not used by postw90.x).
  !! Types used by both wannier90.x and postw90.x are defined in types.F90.
  !! Types specific to postw90.x (not used by wannier90.x) are defined in postw90/postw90_types.F90.

  use w90_constants, only: dp, maxlen

  implicit none

  public

  type w90_calculation_type
    !!==================================================
    !! Contains variables to control the execution path of the program.
    !!==================================================
    logical :: postproc_setup = .false.
    character(len=20) :: restart = ' '
    logical :: bands_plot = .false. !hamiltonian (setup only), plot, wannier_lib
    logical :: wannier_plot = .false. !plot, wannier_lib
    logical :: fermi_surface_plot = .false. ! plot, wannier_lib!
    logical :: transport = .false. ! also hamiltonian, wannier_prog, wannier_lib
  end type w90_calculation_type

  type output_file_type
    logical :: write_hr = .false. !plot, transport and wannier_lib
    logical :: write_r2mn = .false.
    logical :: write_proj = .false.
    logical :: write_hr_diag = .false.
    logical :: write_vdw_data = .false.
    ! aam: for WF-based calculation of vdW C6 coefficients
    logical :: write_u_matrices = .false.
    logical :: write_bvec = .false.
    logical :: write_rmn = .false.
    logical :: write_tb = .false.
    logical :: write_xyz = .false. !wannierise and transport
    logical :: svd_omega = .false. !print svd omega invormation (if iprint >= 2)
  end type output_file_type

  type real_space_ham_type
    !!==================================================
    !! Contains information to control the structure of the real-space Hamiltonian
    !! and how it is calculated.
    !!==================================================
    real(kind=dp) :: hr_cutoff = 0.0_dp !plot and transport
    ! dist_cutoff - only plot and transport
    real(kind=dp) :: dist_cutoff = 1000.0_dp !plot and transport
    character(len=20) :: dist_cutoff_mode = 'three_dim' !plot and transport
    real(kind=dp) :: dist_cutoff_hc = 1000.0_dp !plot and transport
    integer :: one_dim_dir = 0 ! transport and plot
    ! REVIEW_2021-07-22: system_dim is really providing information about the dimensionality
    ! REVIEW_2021-07-22: of the system. Whilst currently it is only used for plotting, its
    ! REVIEW_2021-07-22: use may be generalised in the future. Therefore it makes more sense
    ! REVIEW_2021-07-22: to put it here.
    integer :: system_dim = 3
    ! REVIEW_2021-07-22: There's been some discussion in the past about generalising
    ! REVIEW_2021-07-22: the use of translate_home_cell to also take into account of
    ! REVIEW_2021-07-22: changes in H(R) when WFs are translated. As this is something
    ! REVIEW_2021-07-22: we plan to do, translate_home_cell should probably be here
    logical :: translate_home_cell = .false. ! currently used by wann_write_xyz when write_xyz=.true.
    ! REVIEW_2021-08-09: future plan is that these variables (translation_centre_frac and
    ! REVIEW_2021-08-09: automatic_translation will also result in the hamiltonian being
    ! REVIEW_2021-08-09: modified to be consistent with the translated Wannier centres.
    ! REVIEW_2021-08-09: This is related to Issue 39 in the main repo.
    real(kind=dp) :: translation_centre_frac(3) = 0.0_dp
    ! For Hamiltonian matrix in WF representation
    logical              :: automatic_translation = .true.
  end type real_space_ham_type

  type band_plot_type
    !!==================================================
    !! Contains information to control how the bandstructure plotting is performed and formatted.
    !!==================================================
    character(len=20) :: mode = 's-k' !hamiltonian (setup only), plot
    character(len=20) :: format = 'gnuplot'
    integer, allocatable :: project(:)
  end type band_plot_type

  type wannier_plot_type
    !!==================================================
    !! Contains information for how to plot the wannier functions.
    !!==================================================
    integer, allocatable :: list(:)
    integer :: supercell(3) = 2
    real(kind=dp) :: radius = 3.5_dp
    real(kind=dp) :: scale = 1.0_dp
    character(len=20) :: format = 'xcrysden'
    character(len=20) :: mode = 'crystal'
    character(len=20) :: spinor_mode = 'total'
    logical :: spinor_phase = .true.
  end type wannier_plot_type

  type wvfn_read_type ! only in plot.F90
    !!==================================================
    !! Contains information for how to read the wavefunction files
    !!==================================================
    logical :: formatted = .false.
    !! Read the wvfn from fortran formatted file
    integer :: spin_channel = 1
    !! Spin up=1 down=2
  end type wvfn_read_type

  ! parameters used to control the minimisation of the disentanglement process
  type dis_control_type
    !!==================================================
    !! Contains parameters that control the disentanglement minimisation procedure
    !!==================================================
    integer :: num_iter = 200
    !! number of disentanglement iteration steps
    real(kind=dp) :: mix_ratio = 0.5_dp
    !! Mixing ratio for the disentanglement routine
    real(kind=dp) :: conv_tol = 1.0e-10_dp
    !! Convergence tolerance for the disentanglement
    integer :: conv_window = 3
    !! Size of the convergence window for disentanglement
  end type dis_control_type

  type dis_spheres_type
    ! GS-start
    integer :: first_wann = 1
    integer :: num = 0
    real(kind=dp), allocatable :: spheres(:, :)
    ! GS-end
  end type dis_spheres_type

  type wann_slwf_type
    !!==================================================
    !! Contains parameters that control the selective localisation and constrained centres algorithm
    !!==================================================
    integer :: slwf_num !num_wann
    !! Number of objective Wannier functions (others excluded from spread functional)
    logical :: selective_loc = .false.
    !! Selective localization
    logical :: constrain = .false.
    !! Constrained centres in Cartesian coordinates in angstrom.
    real(kind=dp), allocatable :: centres(:, :)
    real(kind=dp) :: lambda = 1.0_dp
    !! Centre constraints for each Wannier function. Co-ordinates of centre constraint defaults
    !! to centre of trial orbital. Individual Lagrange multipliers, lambdas, default to global Lagrange multiplier.
  end type wann_slwf_type

  type guiding_centres_type
    logical :: enable = .false.
    integer :: num_guide_cycles = 1
    integer :: num_no_guide_iter = 0
    real(kind=dp), allocatable :: centres(:, :)
  end type guiding_centres_type

  type wann_control_type ! only in wannierise.F90
    !!==================================================
    !! Contains parameters that control the wannierisation minimisation procedure
    !!==================================================
    integer :: num_dump_cycles = 100
    !! Number of steps before writing checkpoint
    integer :: num_print_cycles = 1
    !! Number of steps between writing output
    integer :: num_iter = 100
    !! Number of wannierisation iterations
    integer :: num_cg_steps = 5
    !! Number of Conjugate Gradient steps
    real(kind=dp) :: conv_tol = 1.0e-10_dp
    integer :: conv_window
    type(guiding_centres_type) :: guiding_centres
    logical :: use_ss_functional = .false.
    real(kind=dp) :: fixed_step = -999.0_dp
    real(kind=dp) :: trial_step = 2.0_dp
    logical :: precond = .false.
    logical :: lfixstep = .false. ! derived from input
    real(kind=dp) :: conv_noise_amp = -1.0_dp
    integer :: conv_noise_num = 3
    type(wann_slwf_type) :: constrain
  end type wann_control_type

  type wann_omega_type
    !!==================================================
    !! Contains the total spread and its decomposition into gauge invariant and gauge dependent parts.
    !!==================================================
    ! REVIEW_2021-08-04: This type is mainly used for the library mode to be returned back to the external code.
    ! REVIEW_2021-08-04: Internally the code mostly uses the localisation_vars type in wannierise.F90.
    !==================================================
    real(kind=dp) :: invariant = -999.0_dp !wannierise, disentangle and chk2chk
    real(kind=dp) :: total = -999.0_dp
    real(kind=dp) :: tilde = -999.0_dp
  end type wann_omega_type

  ! REVIEW_2021-08-09: We are thinking that this functionality should be probably moved to postw90 at some point.
  type fermi_surface_plot_type
    !!==================================================
    !! Contains variables to control Fermi surface plotting in the main wannier code.
    !!==================================================
    integer :: num_points = 50
    character(len=20) :: plot_format = 'xcrysden'
  end type fermi_surface_plot_type

  ! REVIEW_2021-08-09: This functionality should be moved to postw90 at some point.
  ! REVIEW_2021-08-09: See Issue 31 in the main repo.
  type transport_type ! transport.F90
    !!==================================================
    !! Contains variables to control the calculation of quantum (Landauer-Buttiker) transport
    !!==================================================
    logical :: easy_fix = .false.
    character(len=20) :: mode = 'bulk' ! also hamiltonian
    real(kind=dp) :: win_min = -3.0_dp
    real(kind=dp) :: win_max = 3.0_dp
    real(kind=dp) :: energy_step = 0.01_dp
    integer :: num_bb = 0
    integer :: num_ll = 0
    integer :: num_rr = 0
    integer :: num_cc = 0
    integer :: num_lc = 0
    integer :: num_cr = 0
    integer :: num_bandc = 0
    logical :: write_ht = .false.
    logical :: read_ht = .false. ! also wannier_prog
    logical :: use_same_lead = .true.
    integer :: num_cell_ll = 0
    integer :: num_cell_rr = 0
    real(kind=dp) :: group_threshold = 0.15_dp
  end type transport_type

  ! projections selection - overlap.F90
  ! REVIEW_2021-08-09: At first sight it might appear that select_projections should go in
  ! REVIEW_2021-08-09: the proj_input_type container; but the way the code is structured
  ! REVIEW_2021-08-09: makes this less appealing because there are two proj_input_type variables
  ! REVIEW_2021-08-09: proj_input and proj, the latter containing the subset of selected projections.
  ! REVIEW_2021-08-09: Perhaps best to keep it as currently coded (for now at least).
  type select_projection_type
    !!==================================================
    !! Contains variables relevant to selecting a subset of the projections for the calculation.
    !!==================================================
    logical :: lselproj = .false.
    !integer, save :: num_select_projections
    !integer, allocatable, save :: select_projections(:)
    integer, allocatable :: proj2wann_map(:)
    logical :: auto_projections = .false.
    !! a u t o m a t i c   p r o j e c t i o
  end type select_projection_type

  ! from sitesym
  type sitesym_type
    ! Variables and parameters needed by other modules
    integer :: nkptirr = 9999
    integer :: nsymmetry = 9999
    integer, allocatable :: kptsym(:, :), ir2ik(:), ik2ir(:)
    real(kind=dp) :: symmetrize_eps = 1.d-3
    complex(kind=dp), allocatable :: d_matrix_band(:, :, :, :)
    complex(kind=dp), allocatable :: d_matrix_wann(:, :, :, :)
  end type sitesym_type

  ! from hamiltonian
  type ham_logical_type
    logical :: ham_have_setup = .false.
    logical :: have_ham_k = .false.
    logical :: have_ham_r = .false.
    logical :: have_translated = .false.
    logical :: use_translation = .false.
  end type ham_logical_type

end module w90_wannier90_types
