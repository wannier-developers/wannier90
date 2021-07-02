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

  !BGS probably can be w90 only
  type param_driver_type
    logical :: explicit_nnkpts
    !! nnkpts block is in the input file (allowed only for post-proc setup)
    logical :: postproc_setup
    character(len=20) :: restart
    character(len=20) :: checkpoint
    ! Are we running as a library
    !logical :: library = .false.
    ! Projection data in wannier_lib
    !type(projection_type) :: proj
  end type param_driver_type

  type w90_calculation_type
    logical :: disentanglement !disentangle, overlap, wannier_prog, wannier_lib
    logical :: bands_plot !hamiltonian (setup only), plot, wannier_lib
    logical :: wannier_plot !plot, wannier_lib
    logical :: fermi_surface_plot ! plot, wannier_lib!
    logical :: transport ! also hamiltonian, wannier_prog, wannier_lib
    logical :: write_hr !plot, transport and wannier_lib
    logical :: write_r2mn
    logical :: write_proj
    logical :: write_hr_diag
    logical :: translate_home_cell ! MAYBE used by wann_write_xyz when write_xyz=.true.
    logical :: write_vdw_data
    ! aam: for WF-based calculation of vdW C6 coefficients
    logical :: write_u_matrices
    logical :: write_bvec
    logical :: write_rmn
    logical :: write_tb
    logical :: write_xyz !wannierise and transport
  end type w90_calculation_type

  type param_plot_type ! only in plot.F90
    logical :: wvfn_formatted
    !! Read the wvfn from fortran formatted file
    integer :: spin
    !! Spin up=1 down=2
    integer :: num_wannier_plot
    integer, allocatable :: wannier_plot_list(:)
    integer :: wannier_plot_supercell(3)
    real(kind=dp) :: wannier_plot_radius
    real(kind=dp) :: wannier_plot_scale
    character(len=20) :: wannier_plot_format
    character(len=20) :: wannier_plot_mode
    character(len=20) :: wannier_plot_spinor_mode
    logical :: wannier_plot_spinor_phase
    integer :: bands_num_points
    character(len=20) :: bands_plot_format
    integer, allocatable :: bands_plot_project(:)
    integer :: num_bands_project
    integer :: bands_plot_dim
    !BGS bands_ in band_plot_type? (with ones in param_input)
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

  type param_wannierise_type ! only in wannierise.F90
    integer :: num_dump_cycles
    !! Number of steps before writing checkpoint
    integer :: num_print_cycles
    !! Number of steps between writing output
    integer :: slwf_num
    !! Number of objective Wannier functions (others excluded from spread functional)
    logical :: selective_loc
    !! Selective localization
    logical :: slwf_constrain
    !! Constrained centres
    real(kind=dp), allocatable :: ccentres_cart(:, :)
    real(kind=dp) :: slwf_lambda
    !! Centre constraints for each Wannier function. Co-ordinates of centre constraint defaults
    !! to centre of trial orbital. Individual Lagrange multipliers, lambdas, default to global Lagrange multiplier.
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
    logical :: translate_home_cell
    real(kind=dp) :: omega_total
    real(kind=dp) :: omega_tilde
    ! Projections
    real(kind=dp), allocatable :: proj_site(:, :)
  end type param_wannierise_type

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
