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

module wannier_param_types
  !! This module contains parameters to control the actions of wannier90.
  !! Also routines to read the parameters and write them out again.

  use w90_constants, only: dp
  use w90_io, only: maxlen

  implicit none

  public

  type w90_calculation_type
    !! =========================
    !! Contains variables to control the execution path of the program.
    !! =========================
    logical :: postproc_setup
    character(len=20) :: restart
    logical :: bands_plot !hamiltonian (setup only), plot, wannier_lib
    logical :: wannier_plot !plot, wannier_lib
    logical :: fermi_surface_plot ! plot, wannier_lib!
    logical :: transport ! also hamiltonian, wannier_prog, wannier_lib
  end type w90_calculation_type

  type output_file_type
    logical :: write_hr !plot, transport and wannier_lib
    logical :: write_r2mn
    logical :: write_proj
    logical :: write_hr_diag
    logical :: write_vdw_data
    ! aam: for WF-based calculation of vdW C6 coefficients
    logical :: write_u_matrices
    logical :: write_bvec
    logical :: write_rmn
    logical :: write_tb
    logical :: write_xyz !wannierise and transport
    ! REVIEW_2021-07-22: There's been some discussion in the past about generalising
    ! REVIEW_2021-07-22: the use of translate_home_cell to also take into account of
    ! REVIEW_2021-07-22: changes in H(R) when WFs are translated. As this is something
    ! REVIEW_2021-07-22: we plan to do, translate_home_cell probably should not be in
    ! REVIEW_2021-07-22: this type as it will end up being more general.
    logical :: translate_home_cell ! MAYBE used by wann_write_xyz when write_xyz=.true.
  end type output_file_type

  type real_space_ham_type
    !! ===============================
    !! Contains information to control the structure of the real-space Hamiltonian
    !! and how it is calculated.
    !! ===============================
    real(kind=dp) :: hr_cutoff !plot and transport
    ! dist_cutoff - only plot and transport
    real(kind=dp) :: dist_cutoff !plot and transport
    character(len=20) :: dist_cutoff_mode !plot and transport
    real(kind=dp) :: dist_cutoff_hc !plot and transport
    integer :: one_dim_dir ! transport and plot
    ! REVIEW_2021-07-22: system_dim is really providing information about the dimensionality
    ! REVIEW_2021-07-22: of the system. Whilst currently it is only used for plotting, its
    ! REVIEW_2021-07-22: use may be generalised in the future. Therefore it makes more sense
    ! REVIEW_2021-07-22: to put it here.
    integer :: system_dim
    ! REVIEW_2021-08-09: future plan is that these variables (translation_centre_frac and
    ! REVIEW_2021-08-09: automatic_translation will also result in the hamiltonian being
    ! REVIEW_2021-08-09: modified to be consistent with the translated Wannier centres.
    ! REVIEW_2021-08-09: This is related to Issue 39 in the main repo.
    real(kind=dp) :: translation_centre_frac(3)
    ! For Hamiltonian matrix in WF representation
    logical              :: automatic_translation
  end type real_space_ham_type

  type band_plot_type
    !! ========================
    !! Contains information to control how the bandstructure plotting is performed and formatted.
    !! ========================
    character(len=20) :: plot_mode !hamiltonian (setup only), plot
    character(len=20) :: plot_format
    integer, allocatable :: plot_project(:)
    ! REVIEW_2021-07-22: num_project can be removed (similar to num_exclude_bands -- see discussion there)
    integer :: num_project
  end type band_plot_type

  type wannier_plot_type
    !! ======================
    !! Contains information for how to plot the wannier functions.
    !! ======================
    integer :: num_plot
    integer, allocatable :: plot_list(:)
    integer :: plot_supercell(3)
    real(kind=dp) :: plot_radius
    real(kind=dp) :: plot_scale
    character(len=20) :: plot_format
    character(len=20) :: plot_mode
    character(len=20) :: plot_spinor_mode
    logical :: plot_spinor_phase
  end type wannier_plot_type

  type wvfn_read_type ! only in plot.F90
    !! ========================
    !! Contains information for how to read the wavefunction files
    !! ========================
    logical :: formatted
    !! Read the wvfn from fortran formatted file
    integer :: spin_channel
    !! Spin up=1 down=2
  end type wvfn_read_type

  ! parameters used to control the minimisation of the disentanglement process
  type dis_control_type
    !! ========================
    !! Contains parameters that control the disentanglement minimisation procedure
    !! ========================
    integer :: num_iter
    !! number of disentanglement iteration steps
    real(kind=dp) :: mix_ratio
    !! Mixing ratio for the disentanglement routine
    real(kind=dp) :: conv_tol
    !! Convergence tolerance for the disentanglement
    integer :: conv_window
    !! Size of the convergence window for disentanglement
    ! REVIEW_2021-08-04: these spheres variables belong in disentangle_manifold_type
    ! GS-start
    integer :: spheres_first_wann
    integer :: spheres_num
    real(kind=dp), allocatable :: spheres(:, :)
    ! GS-end
  end type dis_control_type

  ! REVIEW_2021-08-04: Nest this type in wann_control_type
  type wann_slwf_type
    !! ===================
    !! Contains parameters that control the selective localisation and constrained centres algorithm
    !! ===================
    integer :: slwf_num
    !! Number of objective Wannier functions (others excluded from spread functional)
    logical :: selective_loc
    !! Selective localization
    logical :: constrain
    !! Constrained centres in Cartesian coordinates in angstrom.
    real(kind=dp), allocatable :: centres(:, :)
    real(kind=dp) :: lambda
    !! Centre constraints for each Wannier function. Co-ordinates of centre constraint defaults
    !! to centre of trial orbital. Individual Lagrange multipliers, lambdas, default to global Lagrange multiplier.
  end type wann_slwf_type

  type wann_control_type ! only in wannierise.F90
    !! ========================
    !! Contains parameters that control the wannierisation minimisation procedure
    !! ========================
    integer :: num_dump_cycles
    !! Number of steps before writing checkpoint
    integer :: num_print_cycles
    !! Number of steps between writing output
    integer :: num_iter
    !! Number of wannierisation iterations
    integer :: num_cg_steps
    !! Number of Conjugate Gradient steps
    real(kind=dp) :: conv_tol
    integer :: conv_window
    ! REVIEW_2021-08-04: move these three guiding centres parameters to a new type (see below) which can be
    ! REVIEW_2021-08-04: nested inside this type.
    logical :: guiding_centres
    integer :: num_guide_cycles
    integer :: num_no_guide_iter
    ! REVIEW_2021-08-04: make a new guiding_centres_type that includes guiding_centres(logical) [which we can rename "enable"]
    ! REVIEW_2021-08-04: num_guide_cycles, num_no_guide_iter and guiding_centres(real array) [which we can rename "centres"]
    real(kind=dp) :: fixed_step
    real(kind=dp) :: trial_step
    logical :: precond
    logical :: lfixstep ! derived from input
    real(kind=dp) :: conv_noise_amp
    integer :: conv_noise_num
    ! REVIEW_2021-08-04: Nest wann_slwf_type in here (currently called wann_localise_type)
  end type wann_control_type

  type wann_omega_type
    !! ===================
    !! Contains the total spread and its decomposition into gauge invariant and gauge dependent parts.
    !! ===================
    ! REVIEW_2021-08-04: This type is mainly used for the library mode to be returned back to the external code.
    ! REVIEW_2021-08-04: Internally the code mostly uses the localisation_vars type in wannierise.F90.
    ! ===================
    real(kind=dp) :: invariant !wannierise, disentangle and chk2chk
    real(kind=dp) :: total
    real(kind=dp) :: tilde
  end type wann_omega_type

  ! REVIEW_2021-08-04: We don't think there is a need to have this type.
  type wannierise_type
    type(wann_control_type) :: control
    ! REVIEW_2021-08-04: Nest this type inside wann_control_type (see comments above)
    type(wann_slwf_type) :: constrain
    type(wann_omega_type) :: omega
    ! Projections, used when guiding_centres=.true. or for constrained if slwf_constrain=.true.
    ! REVIEW_2021-08-04: See whether it is possible to avoid using this variable and instead use
    ! REVIEW_2021-08-04: guiding_centres_type%centres or wann_slwf_type%centres as appropriate and
    ! REVIEW_2021-08-04: remove proj_site altogether.
    real(kind=dp), allocatable :: proj_site(:, :) ! MAYBE
  end type wannierise_type

  ! REVIEW_2021-08-09: We are thinking that this functionality should be probably moved to postw90 at some point.
  type fermi_surface_plot_type
    !! ==================
    !! Contains variables to control Fermi surface plotting in the main wannier code.
    !! ==================
    integer :: num_points
    character(len=20) :: plot_format
  end type fermi_surface_plot_type

  ! REVIEW_2021-08-09: This functionality should be moved to postw90 at some point.
  ! REVIEW_2021-08-09: See Issue 31 in the main repo.
  type transport_type ! transport.F90
    !! =================
    !! Contains variables to control the calculation of quantum (Landauer-Buttiker) transport
    !! =================
    logical :: easy_fix
    character(len=20) :: mode ! also hamiltonian
    real(kind=dp) :: win_min
    real(kind=dp) :: win_max
    real(kind=dp) :: energy_step
    integer :: num_bb
    integer :: num_ll
    integer :: num_rr
    integer :: num_cc
    integer :: num_lc
    integer :: num_cr
    integer :: num_bandc
    logical :: write_ht
    logical :: read_ht ! also wannier_prog
    logical :: use_same_lead
    integer :: num_cell_ll
    integer :: num_cell_rr
    real(kind=dp) :: group_threshold
  end type transport_type

  ! projections selection - overlap.F90
  ! REVIEW_2021-08-09: Check that select_projections(:) is deallocated
  ! REVIEW_2021-08-09: At first sight it might appear that select_projections should go in
  ! REVIEW_2021-08-09: the proj_input_type container; but the way the code is structured
  ! REVIEW_2021-08-09: makes this less appealing because there are two proj_input_type variables
  ! REVIEW_2021-08-09: proj_input and proj, the latter containing the subset of selected projections.
  ! REVIEW_2021-08-09: Perhaps best to keep it as currently coded (for now at least).
  type select_projection_type
    !! ================
    !! Contains variables relevant to selecting a subset of the projections for the calculation.
    !! ================
    logical :: lselproj
    !integer, save :: num_select_projections
    !integer, allocatable, save :: select_projections(:)
    integer, allocatable :: proj2wann_map(:)
  end type select_projection_type

end module wannier_param_types
