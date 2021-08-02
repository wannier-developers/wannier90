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

  type band_plot_type
    !! ========================
    !! Contains information to control how the bandstructure plotting is performed and formatted.
    !! ========================
    character(len=20) :: plot_mode !hamiltonian (setup only), plot
    ! REVIEW_2021-07-22: move num_points to kpoint_path_type (see discussion in parameters.F90)
    integer :: num_points
    character(len=20) :: plot_format
    integer, allocatable :: plot_project(:)
    ! REVIEW_2021-07-22: num_project can be removed (similar to num_exclude_bands -- see discussion there)
    integer :: num_project
    ! REVIEW_2021-07-22: plot_dim is really providing information about the dimensionality
    ! REVIEW_2021-07-22: of the system. Whilst currently it is only used for plotting, its
    ! REVIEW_2021-07-22: use may be generalised in the future. Therefore it makes more sense
    ! REVIEW_2021-07-22: to put it in real_space_ham_type in parameters.F90 and to call it
    ! REVIEW_2021-07-22: something more general, such as system_dim.
    integer :: plot_dim
    ! others from param_input? MAYBE
  end type band_plot_type

  ! REVIEW_2021-07-22: REVIEWED UP TO HERE SO FAR. WE ARE MEETING AGAIN TO COMPLETE
  ! REVIEW_2021-07-22: REVIEW ON 4 AUG, BEFORE THE NEXT MEETING WITH STFC

  type wannier_plot_type
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

  ! MAYBE - seems almost redundant
  type param_plot_type ! only in plot.F90
    logical :: wvfn_formatted
    !! Read the wvfn from fortran formatted file
    integer :: spin
    !! Spin up=1 down=2
  end type param_plot_type

  ! parameters used to control the minimisation of the disentanglement process
  type disentangle_type
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
  end type disentangle_type

  type wann_control_type ! only in wannierise.F90
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
    logical :: guiding_centres
    integer :: num_guide_cycles
    integer :: num_no_guide_iter
    real(kind=dp) :: fixed_step
    real(kind=dp) :: trial_step
    logical :: precond
    logical :: lfixstep ! derived from input
    real(kind=dp) :: conv_noise_amp
    integer :: conv_noise_num
  end type wann_control_type

  type wann_omega_type
    real(kind=dp) :: invariant !wannierise, disentangle and chk2chk
    real(kind=dp) :: total
    real(kind=dp) :: tilde
  end type wann_omega_type

  type wann_localise_type
    integer :: slwf_num
    !! Number of objective Wannier functions (others excluded from spread functional)
    logical :: selective_loc
    !! Selective localization
    logical :: slwf_constrain ! MAYBE or in main wannier_type
    !! Constrained centres
    real(kind=dp), allocatable :: ccentres_cart(:, :)
    real(kind=dp) :: slwf_lambda
    !! Centre constraints for each Wannier function. Co-ordinates of centre constraint defaults
    !! to centre of trial orbital. Individual Lagrange multipliers, lambdas, default to global Lagrange multiplier.
  end type wann_localise_type

  type param_wannierise_type
    type(wann_control_type) :: control
    type(wann_localise_type) :: constrain
    type(wann_omega_type) :: omega
    ! Projections, used when guiding_centres=.true. or for constrained if slwf_constrain=.true.
    real(kind=dp), allocatable :: proj_site(:, :) ! MAYBE
  end type param_wannierise_type

  ! MAYBE which vars did Arash mean from parameter_input?
  type param_hamiltonian_type
    real(kind=dp) :: translation_centre_frac(3)
    ! For Hamiltonian matrix in WF representation
    logical              :: automatic_translation
  end type param_hamiltonian_type

  type fermi_surface_type
    integer :: num_points
    character(len=20) :: plot_format
  end type fermi_surface_type

  type transport_type ! transport.F90
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
  type select_projection_type
    logical :: lselproj
    !integer, save :: num_select_projections
    !integer, allocatable, save :: select_projections(:)
    integer, allocatable :: proj2wann_map(:)
  end type select_projection_type

end module wannier_param_types
