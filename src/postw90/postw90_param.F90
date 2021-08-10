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

module pw90_parameters

  use w90_constants, only: dp
  use w90_comms, only: w90commtype

  implicit none

  public

  type pw90_calculation_type !postw90.F90
    !! =================
    !! Contains information about the high-level task that pw90 is being asked to do,
    !! including any global variables that affect all calculation branches.
    !! =================
    logical :: kpath
    logical :: kslice
    logical :: dos
    logical :: berry
    logical :: gyrotropic
    logical :: geninterp
    logical :: boltzwann
    logical :: spin_moment !postw90_common and postw90
    logical :: spin_decomp !postw90_common, berry, dos and boltzwann
  end type pw90_calculation_type

  type pw90_oper_read_type ! only in postw90/get_oper.F90
    !! ==============
    !! Contains variables for determining whether formatted or unformatted input is read by get_oper.F90
    !! ==============
    logical :: spn_formatted
    !! Read the spin from fortran formatted file
    logical :: uHu_formatted
    !! Read the uHu from fortran formatted file
  end type pw90_oper_read_type

  ! consider removing this
  ! REVIEW_2021-08-09: See Issue 34 in the main repo for a discussion about scissors_shift
  ! REVIEW_2021-08-09: which has been deprecated and should be removed in a future version of
  ! REVIEW_2021-08-09: the code. For now, please keep it as a standalone variable.
  ! REVIEW_2021-08-09: effective_model should also be a standalone variable.
  type postw90_common_type
    real(kind=dp) :: scissors_shift ! get_oper and berry
    ! IVO
    ! Are we running postw90 starting from an effective model?
    logical :: effective_model = .false.
  end type postw90_common_type

  ! Module  s p i n
  type pw90_spin_mod_type
    !! ===============
    !! Contains variables used in the spin module of postw90
    !! ===============
    real(kind=dp) :: axis_polar
    real(kind=dp) :: axis_azimuth
    real(kind=dp) :: kmesh_spacing
    integer :: kmesh(3)
  end type pw90_spin_mod_type

  ! REVIEW_2021-08-09: A general comment, which we record here so we don't forget.
  ! REVIEW_2021-08-09: It would be good to try to keep type and variable names
  ! REVIEW_2021-08-09: as consistent as possible. As an example, at the top of berry.F90, we have
  ! REVIEW_2021-08-09: type(kmesh_info_type), intent(in) :: kmesh_info   ! consistent
  ! REVIEW_2021-08-09: type(kpoint_dist_type), intent(in) :: kdist       ! not consistent
  ! REVIEW_2021-08-09: type(k_point_type), intent(in) :: k_points        ! fixed
  ! REVIEW_2021-08-09: This would make reading and debugging the code easier in future.

  type pw90_band_deriv_degen_type
    !! ================
    !! Contains variables for doing degenerate perturbation theory when bands are degenerate
    !! and band derivatives are needed.
    !! ================
    logical :: use_degen_pert
    real(kind=dp) :: degen_thr
  end type pw90_band_deriv_degen_type

  ! module  k p a t h (used by postw90/kpath)
  type pw90_kpath_mod_type
    !! ================
    !! Contains control variables for the kpath module of postw90
    !! ================
    character(len=20) :: task
    integer :: num_points
    character(len=20) :: bands_colour
  end type pw90_kpath_mod_type

  ! module  k s l i c e (postw90/kslice)
  type pw90_kslice_mod_type
    !! ===============
    !! Contains control variables for the kslice module of postw90
    !! ===============
    character(len=20) :: task
    real(kind=dp) :: corner(3)
    real(kind=dp) :: b1(3)
    real(kind=dp) :: b2(3)
    integer :: kmesh2d(2)
    character(len=20) :: fermi_lines_colour
  end type pw90_kslice_mod_type

  type pw90_smearing_type
    !! =============
    !! Contains variables for controlling the smearing.
    !! =============
    logical    :: use_adaptive
    real(kind=dp)    :: adaptive_prefactor
    integer    :: type_index
    real(kind=dp)    :: fixed_width
    real(kind=dp)    :: adaptive_max_width
  end type pw90_smearing_type

  type pw90_dos_mod_type
    !! ===============
    !! Contains variables for the dos module of postw90
    !! ===============
    character(len=20)    :: task
    type(pw90_smearing_type) :: smearing
    real(kind=dp)    :: energy_max
    real(kind=dp)    :: energy_min
    real(kind=dp)    :: energy_step
    integer    :: num_project
    integer, allocatable :: project(:)
    !  character(len=20)    :: plot_format
    real(kind=dp)    :: kmesh_spacing
    integer    :: kmesh(3)
    !  real(kind=dp) :: gaussian_width
  end type pw90_dos_mod_type

  ! Module  b e r r y (mainly postw90/berry)
  type pw90_berry_mod_type
    !! =============
    !! Contains variables for the berry module of postw90
    !! =============
    character(len=120) :: task
    real(kind=dp) :: kmesh_spacing
    integer :: kmesh(3)
    integer :: curv_adpt_kmesh
    real(kind=dp) :: curv_adpt_kmesh_thresh
    character(len=20) :: curv_unit ! postw90/kpath, kslice as well
    type(pw90_smearing_type) :: kubo_smearing
    integer :: sc_phase_conv
    real(kind=dp) :: sc_eta ! also postw90/wan_ham
    real(kind=dp) :: sc_w_thr
    logical :: wanint_kpoint_file ! also postw90/spin, postw90/dos, postw90.F90
    logical :: transl_inv !also used in postw90/get_oper, postw90/gyrotropic
    integer :: kubo_nfreq
    complex(kind=dp), allocatable :: kubo_freq_list(:)
    real(kind=dp) :: kubo_eigval_max
  end type pw90_berry_mod_type

  ! spin Hall conductivity (postw90 - common, get_oper, berry, kpath)
  type pw90_spin_hall_type
    !! =============
    !! Contains variables controlling the calculation of spin hall conductivity in postw90
    !! =============
    logical :: freq_scan
    integer :: alpha
    integer :: beta
    integer :: gamma
    logical :: bandshift
    integer :: bandshift_firstband
    real(kind=dp) :: bandshift_energyshift
  end type pw90_spin_hall_type

  type pw90_gyrotropic_type ! postw90 - common, gyrotropic
    !! =============
    !! Contains variables for the gyrotropic module of postw90
    !! =============
    character(len=120) :: task
    integer :: kmesh(3)
    real(kind=dp) :: kmesh_spacing
    ! REVIEW_2021-08-09: Should this use pw90_smearing_type?
    integer :: smr_index
    real(kind=dp) :: smr_fixed_en_width
    integer :: nfreq
    complex(kind=dp), allocatable :: freq_list(:)
    real(kind=dp) :: box_corner(3), box(3, 3)
    real(kind=dp) :: degen_thresh
    integer, allocatable :: band_list(:)
    integer :: num_bands
    ! REVIEW_2021-08-09: Should this use pw90_smearing_type?
    ! REVIEW_2021-08-09: Is this a speed-up that could be applied more generally?
    real(kind=dp) :: smr_max_arg
    real(kind=dp) :: eigval_max
  end type pw90_gyrotropic_type

  ! [gp-begin, Jun 1, 2012]
  ! GeneralInterpolator variables - postw90/geninterp
  type pw90_geninterp_mod_type
    !! =============
    !! Contains variables for the geninterp module of postw90
    !! =============
    logical :: alsofirstder
    logical :: single_file
  end type pw90_geninterp_mod_type
  ! [gp-end, Jun 1, 2012]

  ! [gp-begin, Apr 12, 2012]
  ! BoltzWann variables (postw90/boltzwann.F90)
  type pw90_boltzwann_type
    !! =============
    !! Contains variables for the boltzwann module of postw90
    !! =============
    logical :: calc_also_dos
    integer :: dir_num_2d
    real(kind=dp) :: dos_energy_step
    real(kind=dp) :: dos_energy_min
    real(kind=dp) :: dos_energy_max
    type(pw90_smearing_type) :: dos_smearing
    real(kind=dp) :: mu_min
    real(kind=dp) :: mu_max
    real(kind=dp) :: mu_step
    real(kind=dp) :: temp_min
    real(kind=dp) :: temp_max
    real(kind=dp) :: temp_step
    real(kind=dp) :: kmesh_spacing
    integer :: kmesh(3)
    real(kind=dp) :: tdf_energy_step
    ! REVIEW_2021-08-09: Should this use pw90_smearing_type?
    ! REVIEW_2021-08-09: If we use the smearing type, then rename tdf_smearing
    integer :: TDF_smr_index
    real(kind=dp) :: TDF_smr_fixed_en_width
    real(kind=dp) :: relax_time
    logical :: bandshift
    integer :: bandshift_firstband
    real(kind=dp) :: bandshift_energyshift
  end type pw90_boltzwann_type
  ! [gp-end, Apr 12, 2012]

end module pw90_parameters

module pw90_param_methods

  use w90_constants, only: dp
  use w90_io, only: maxlen
  use w90_param_types, only: print_output_type, print_output_type, wannier_data_type, &
    kmesh_input_type, kmesh_info_type, k_points_type, dis_manifold_type, &
    fermi_data_type, atom_data_type, kpoint_path_type, proj_input_type, w90_system_type, &
    exclude_bands_type, ws_region_type
  use w90_param_methods
  use pw90_parameters

  implicit none

  private

  ! These could be local to param_read if they weren't also used by param_write
  type pw90_extra_io_type
    ! from gyrotropic section
    real(kind=dp)                               :: gyrotropic_freq_min
    real(kind=dp)                               :: gyrotropic_freq_max
    real(kind=dp)                               :: gyrotropic_freq_step
    real(kind=dp)                               :: kubo_freq_min
    real(kind=dp)                               :: kubo_freq_max
    real(kind=dp)                               :: kubo_freq_step
    ! Adaptive vs. fixed smearing stuff [GP, Jul 12, 2012]
    ! Only internal, always use the local variables defined by each module
    ! that take this value as default
    type(pw90_smearing_type) :: smear
    ! [gp-begin, Apr 13, 2012]
    ! Global interpolation k mesh variables
    ! These don't need to be public, since their values are copied in the variables of the
    ! local interpolation meshes. JRY: added save attribute
    real(kind=dp)             :: kmesh_spacing
    integer                   :: kmesh(3)
    logical                   :: global_kmesh_set
    ! [gp-end]
    character(len=4) :: boltz_2d_dir ! this could be local to read_boltzwann
  end type pw90_extra_io_type

  public :: pw90_extra_io_type
  public :: param_postw90_read
  public :: param_postw90_write

contains

  subroutine param_postw90_read(rs_region, system, excluded_bands, verbose, wann_data, &
                                kmesh_data, k_points, num_kpts, dis_window, fermi, atoms, &
                                num_bands, num_wann, eigval, mp_grid, real_lattice, &
                                recip_lattice, spec_points, pw90_calcs, &
                                postw90_oper, pw90_common, pw90_spin, pw90_ham, &
                                kpath, kslice, dos_data, berry, spin_hall, &
                                gyrotropic, geninterp, boltz, eig_found, write_data, &
                                gamma_only, bohr, stdout, seedname)
    !==================================================================!
    !                                                                  !
    !! Read parameters and calculate derived values
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !                                                                  !
    !===================================================================

    implicit none

    ! arguments
    type(atom_data_type), intent(inout) :: atoms
    type(pw90_berry_mod_type), intent(inout) :: berry
    type(pw90_boltzwann_type), intent(inout) :: boltz
    type(dis_manifold_type), intent(inout) :: dis_window
    type(pw90_dos_mod_type), intent(inout) :: dos_data
    type(exclude_bands_type), intent(inout) :: excluded_bands
    type(fermi_data_type), intent(inout) :: fermi
    type(pw90_geninterp_mod_type), intent(inout) :: geninterp
    type(pw90_gyrotropic_type), intent(inout) :: gyrotropic
    type(pw90_kpath_mod_type), intent(inout) :: kpath
    type(k_points_type), intent(inout) :: k_points
    type(pw90_kslice_mod_type), intent(inout) :: kslice
    type(kmesh_input_type), intent(inout) :: kmesh_data
    type(postw90_common_type), intent(inout) :: pw90_common
    type(pw90_band_deriv_degen_type), intent(inout) :: pw90_ham
    type(pw90_oper_read_type), intent(inout) :: postw90_oper
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    type(print_output_type), intent(inout) :: verbose
    type(pw90_calculation_type), intent(inout) :: pw90_calcs
    type(pw90_extra_io_type), intent(inout) :: write_data
    type(ws_region_type), intent(inout) :: rs_region
    type(kpoint_path_type), intent(inout) :: spec_points
    type(pw90_spin_hall_type), intent(inout) :: spin_hall
    type(w90_system_type), intent(inout) :: system
    type(wannier_data_type), intent(inout) :: wann_data

    integer, intent(inout) :: mp_grid(3)
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_kpts
    integer, intent(inout) :: num_wann
    integer, intent(in) :: stdout

    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: recip_lattice(3, 3)

    character(len=50), intent(in)  :: seedname

    logical, intent(inout) :: eig_found
    logical, intent(inout) :: gamma_only

    ! local variables
    logical :: dos_plot
    logical :: found_fermi_energy
    logical :: disentanglement, library, ok
    character(len=20) :: energy_unit

    library = .false.
    call param_in_file(seedname, stdout)
    call param_read_verbosity(verbose, stdout, seedname)
    call param_read_pw90_calcs(pw90_calcs, stdout, seedname)
    call param_read_effective_model(pw90_common%effective_model, stdout, seedname)
    call param_read_units(verbose%lenconfac, verbose%length_unit, energy_unit, bohr, &
                          stdout, seedname)
    call param_read_oper(postw90_oper, stdout, seedname)
    call param_read_num_wann(num_wann, stdout, seedname)
    call param_read_exclude_bands(excluded_bands, stdout, seedname) !for read_chkpt
    call param_read_num_bands(pw90_common%effective_model, library, &
                              excluded_bands%num_exclude_bands, num_bands, num_wann, .false., &
                              stdout, seedname)
    disentanglement = (num_bands > num_wann)
    call param_read_devel(verbose%devel_flag, stdout, seedname)
    call param_read_mp_grid(pw90_common%effective_model, library, mp_grid, num_kpts, &
                            stdout, seedname)
    call param_read_gamma_only(gamma_only, num_kpts, library, stdout, seedname)
    call param_read_system(library, system, stdout, seedname)
    call param_read_kpath(library, spec_points, ok, .false., stdout, seedname)
    call param_read_fermi_energy(found_fermi_energy, fermi, stdout, seedname)
    call param_read_kslice(pw90_calcs%kslice, kslice, stdout, seedname)
    call param_read_smearing(write_data%smear, stdout, seedname)
    call param_read_scissors_shift(pw90_common, stdout, seedname)
    call param_read_pw90spin(pw90_calcs%spin_moment, pw90_calcs%spin_decomp, pw90_spin, &
                             system%num_elec_per_state, stdout, seedname)
    call param_read_gyrotropic(gyrotropic, num_wann, write_data%smear%fixed_width, &
                               write_data%smear%type_index, stdout, seedname)
    call param_read_berry(pw90_calcs, berry, write_data%smear, stdout, seedname)
    call param_read_spin_hall(pw90_calcs, pw90_common, spin_hall, berry%task, stdout, seedname)
    call param_read_pw90ham(pw90_ham, stdout, seedname)
    call param_read_pw90_kpath(pw90_calcs, kpath, spec_points, stdout, seedname)
    call param_read_dos(pw90_calcs, dos_data, found_fermi_energy, num_wann, write_data%smear, &
                        dos_plot, stdout, seedname)
    call param_read_ws_data(rs_region, stdout, seedname)
    call param_read_eigvals(pw90_common%effective_model, pw90_calcs%boltzwann, &
                            pw90_calcs%geninterp, dos_plot, disentanglement, eig_found, eigval, &
                            library, .false., num_bands, num_kpts, stdout, seedname)
    dis_window%win_min = -1.0_dp
    dis_window%win_max = 0.0_dp
    if (eig_found) dis_window%win_min = minval(eigval)
    if (eig_found) dis_window%win_max = maxval(eigval)
    call param_read_dis_manifold(eig_found, dis_window, stdout, seedname)
    call param_read_geninterp(geninterp, stdout, seedname)
    call param_read_boltzwann(boltz, eigval, write_data%smear, pw90_calcs%boltzwann, &
                              write_data%boltz_2d_dir, stdout, seedname)
    call param_read_energy_range(berry, dos_data, gyrotropic, dis_window, fermi, eigval, write_data, &
                                 stdout, seedname)
    call param_read_lattice(library, real_lattice, recip_lattice, bohr, stdout, seedname)
    call param_read_kmesh_data(kmesh_data, stdout, seedname)
    call param_read_kpoints(pw90_common%effective_model, library, k_points, num_kpts, &
                            recip_lattice, bohr, stdout, seedname)
    call param_read_global_kmesh(write_data%global_kmesh_set, write_data%kmesh_spacing, &
                                 write_data%kmesh, recip_lattice, stdout, seedname)
    call param_read_local_kmesh(pw90_calcs, berry, dos_data, pw90_spin, gyrotropic, &
                                boltz, recip_lattice, write_data%global_kmesh_set, &
                                write_data%kmesh, write_data%kmesh_spacing, stdout, seedname)
    call param_read_atoms(library, atoms, real_lattice, recip_lattice, bohr, stdout, seedname) !pw90_write
    call param_clean_infile(stdout, seedname)
    ! For aesthetic purposes, convert some things to uppercase
    call param_uppercase(atoms, spec_points, verbose%length_unit)
    call param_read_final_alloc(disentanglement, dis_window, wann_data, num_wann, num_bands, num_kpts, stdout, seedname)
  end subroutine param_postw90_read

  subroutine param_read_pw90_calcs(pw90_calcs, stdout, seedname)
    implicit none
    integer, intent(in) :: stdout
    type(pw90_calculation_type), intent(out) :: pw90_calcs
    character(len=50), intent(in)  :: seedname
    logical :: found

    pw90_calcs%dos = .false.
    call param_get_keyword(stdout, seedname, 'dos', found, l_value=pw90_calcs%dos)

    pw90_calcs%berry = .false.
    call param_get_keyword(stdout, seedname, 'berry', found, l_value=pw90_calcs%berry)

    pw90_calcs%kpath = .false.
    call param_get_keyword(stdout, seedname, 'kpath', found, l_value=pw90_calcs%kpath)

    pw90_calcs%kslice = .false.
    call param_get_keyword(stdout, seedname, 'kslice', found, l_value=pw90_calcs%kslice)

    pw90_calcs%gyrotropic = .false.
    call param_get_keyword(stdout, seedname, 'gyrotropic', found, l_value=pw90_calcs%gyrotropic)

    pw90_calcs%geninterp = .false.
    call param_get_keyword(stdout, seedname, 'geninterp', found, l_value=pw90_calcs%geninterp)
    pw90_calcs%boltzwann = .false.
    call param_get_keyword(stdout, seedname, 'boltzwann', found, l_value=pw90_calcs%boltzwann)

  end subroutine param_read_pw90_calcs

  subroutine param_read_effective_model(effective_model, stdout, seedname)
    implicit none
    integer, intent(in) :: stdout
    logical, intent(inout) :: effective_model
    character(len=50), intent(in)  :: seedname

    logical :: found

    !ivo
    call param_get_keyword(stdout, seedname, 'effective_model', found, l_value=effective_model)
  end subroutine param_read_effective_model

  subroutine param_read_oper(postw90_oper, stdout, seedname)
    implicit none
    integer, intent(in) :: stdout
    type(pw90_oper_read_type), intent(inout) :: postw90_oper
    character(len=50), intent(in)  :: seedname

    logical :: found

    postw90_oper%spn_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword(stdout, seedname, 'spn_formatted', found, l_value=postw90_oper%spn_formatted)

    postw90_oper%uHu_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword(stdout, seedname, 'uhu_formatted', found, l_value=postw90_oper%uHu_formatted)
  end subroutine param_read_oper

  subroutine param_read_kslice(pw90_kslice, kslice, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: pw90_kslice
    type(pw90_kslice_mod_type), intent(inout) :: kslice
    character(len=50), intent(in)  :: seedname

    integer :: i
    logical :: found

    kslice%task = 'fermi_lines'
    call param_get_keyword(stdout, seedname, 'kslice_task', found, c_value=kslice%task)
    if (pw90_kslice .and. index(kslice%task, 'fermi_lines') == 0 .and. &
        index(kslice%task, 'curv') == 0 .and. &
        index(kslice%task, 'morb') == 0 .and. &
        index(kslice%task, 'shc') == 0) call io_error &
      ('Error: value of kslice_task not recognised in param_read', stdout, seedname)
    if (pw90_kslice .and. index(kslice%task, 'curv') > 0 .and. &
        index(kslice%task, 'morb') > 0) call io_error &
      ("Error: kslice_task cannot include both 'curv' and 'morb'", stdout, seedname)
    if (pw90_kslice .and. index(kslice%task, 'shc') > 0 .and. &
        index(kslice%task, 'morb') > 0) call io_error &
      ("Error: kslice_task cannot include both 'shc' and 'morb'", stdout, seedname)
    if (pw90_kslice .and. index(kslice%task, 'shc') > 0 .and. &
        index(kslice%task, 'curv') > 0) call io_error &
      ("Error: kslice_task cannot include both 'shc' and 'curv'", stdout, seedname)

    kslice%kmesh2d(1:2) = 50
    call param_get_vector_length(stdout, seedname, 'kslice_2dkmesh', found, length=i)
    if (found) then
      if (i == 1) then
        call param_get_keyword_vector(stdout, seedname, 'kslice_2dkmesh', found, 1, &
                                      i_value=kslice%kmesh2d)
        kslice%kmesh2d(2) = kslice%kmesh2d(1)
      elseif (i == 2) then
        call param_get_keyword_vector(stdout, seedname, 'kslice_2dkmesh', found, 2, &
                                      i_value=kslice%kmesh2d)
      else
        call io_error('Error: kslice_2dkmesh must be provided as either' &
                      //' one integer or a vector of two integers', stdout, seedname)
      endif
      if (any(kslice%kmesh2d <= 0)) &
        call io_error('Error: kslice_2dkmesh elements must be' &
                      //' greater than zero', stdout, seedname)
    endif

    kslice%corner = 0.0_dp
    call param_get_keyword_vector(stdout, seedname, 'kslice_corner', found, 3, r_value=kslice%corner)

    kslice%b1(1) = 1.0_dp
    kslice%b1(2) = 0.0_dp
    kslice%b1(3) = 0.0_dp
    call param_get_keyword_vector(stdout, seedname, 'kslice_b1', found, 3, r_value=kslice%b1)

    kslice%b2(1) = 0.0_dp
    kslice%b2(2) = 1.0_dp
    kslice%b2(3) = 0.0_dp
    call param_get_keyword_vector(stdout, seedname, 'kslice_b2', found, 3, r_value=kslice%b2)

    kslice%fermi_lines_colour = 'none'
    call param_get_keyword(stdout, seedname, 'kslice_fermi_lines_colour', found, &
                           c_value=kslice%fermi_lines_colour)
    if (pw90_kslice .and. index(kslice%fermi_lines_colour, 'none') == 0 .and. &
        index(kslice%fermi_lines_colour, 'spin') == 0) call io_error &
      ('Error: value of kslice_fermi_lines_colour not recognised ' &
       //'in param_read', stdout, seedname)

!    slice_plot_format         = 'plotmv'
!    call param_get_keyword('slice_plot_format',found,c_value=slice_plot_format)
  end subroutine param_read_kslice

  subroutine param_read_smearing(smearing, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    type(pw90_smearing_type), intent(out) :: smearing
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    logical :: found
    character(len=maxlen)              :: ctmp
    ! [gp-begin, Apr 20, 2012]

    ! By default: Gaussian
    smearing%type_index = 0
    call param_get_keyword(stdout, seedname, 'smr_type', found, c_value=ctmp)
    if (found) smearing%type_index = get_smearing_index(ctmp, 'smr_type', stdout, seedname)

    ! By default: adaptive smearing
    smearing%use_adaptive = .true.
    call param_get_keyword(stdout, seedname, 'adpt_smr', found, l_value=smearing%use_adaptive)

    ! By default: a=sqrt(2)
    smearing%adaptive_prefactor = sqrt(2.0_dp)
    call param_get_keyword(stdout, seedname, 'adpt_smr_fac', found, r_value=smearing%adaptive_prefactor)
    if (found .and. (smearing%adaptive_prefactor <= 0._dp)) &
      call io_error('Error: adpt_smr_fac must be greater than zero', stdout, seedname)

    ! By default: 1 eV
    smearing%adaptive_max_width = 1.0_dp
    call param_get_keyword(stdout, seedname, 'adpt_smr_max', found, r_value=smearing%adaptive_max_width)
    if (smearing%adaptive_max_width <= 0._dp) &
      call io_error('Error: adpt_smr_max must be greater than zero', stdout, seedname)

    ! By default: if adpt_smr is manually set to false by the user, but he/she doesn't
    ! define smr_fixed_en_width: NO smearing, i.e. just the histogram
    smearing%fixed_width = 0.0_dp
    call param_get_keyword(stdout, seedname, 'smr_fixed_en_width', found, &
                           r_value=smearing%fixed_width)
    if (found .and. (smearing%fixed_width < 0._dp)) &
      call io_error('Error: smr_fixed_en_width must be greater than or equal to zero', stdout, &
                    seedname)
    ! [gp-end]
  end subroutine param_read_smearing

  subroutine param_read_scissors_shift(pw90_common, stdout, seedname)
    implicit none
    integer, intent(in) :: stdout
    type(postw90_common_type), intent(inout) :: pw90_common
    character(len=50), intent(in)  :: seedname

    logical :: found

    pw90_common%scissors_shift = 0.0_dp
    call param_get_keyword(stdout, seedname, 'scissors_shift', found, r_value=pw90_common%scissors_shift)

  end subroutine param_read_scissors_shift

  subroutine param_read_pw90spin(spin_moment, spin_decomp, pw90_spin, num_elec_per_state, &
                                 stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(out) :: spin_moment ! from pw90_calcs
    logical, intent(out) :: spin_decomp ! from pw90_common
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    integer, intent(in) :: num_elec_per_state
    character(len=50), intent(in)  :: seedname

    logical :: found

    spin_moment = .false.
    call param_get_keyword(stdout, seedname, 'spin_moment', found, l_value=spin_moment)

    pw90_spin%axis_polar = 0.0_dp
    call param_get_keyword(stdout, seedname, 'spin_axis_polar', found, r_value=pw90_spin%axis_polar)

    pw90_spin%axis_azimuth = 0.0_dp
    call param_get_keyword(stdout, seedname, 'spin_axis_azimuth', found, r_value=pw90_spin%axis_azimuth)

    spin_decomp = .false.
    call param_get_keyword(stdout, seedname, 'spin_decomp', found, l_value=spin_decomp)

    if (spin_decomp .and. (num_elec_per_state .ne. 1)) then
      call io_error('spin_decomp can be true only if num_elec_per_state is 1', stdout, seedname)
    end if

  end subroutine param_read_pw90spin

  subroutine param_read_gyrotropic(gyrotropic, num_wann, smr_fixed_en_width, smr_index, &
                                   stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_gyrotropic_type), intent(out) :: gyrotropic
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: smr_fixed_en_width
    integer, intent(in) :: smr_index
    character(len=50), intent(in)  :: seedname

    real(kind=dp) :: smr_max_arg
    real(kind=dp)                   :: gyrotropic_box_tmp(3)
    integer :: i, ierr, loop
    logical :: found
    character(len=maxlen)              :: ctmp

    ! Stepan
    gyrotropic%task = 'all'
    call param_get_keyword(stdout, seedname, 'gyrotropic_task', found, c_value=gyrotropic%task)
    gyrotropic%box(:, :) = 0.0
    gyrotropic%degen_thresh = 0.0_dp
    call param_get_keyword(stdout, seedname, 'gyrotropic_degen_thresh', found, r_value=gyrotropic%degen_thresh)

    do i = 1, 3
      gyrotropic%box(i, i) = 1.0_dp
      gyrotropic_box_tmp(:) = 0.0_dp
      call param_get_keyword_vector(stdout, seedname, 'gyrotropic_box_b'//achar(48 + i), found, 3, r_value=gyrotropic_box_tmp)
      if (found) gyrotropic%box(i, :) = gyrotropic_box_tmp(:)
    enddo
    gyrotropic%box_corner(:) = 0.0_dp
    call param_get_keyword_vector(stdout, seedname, 'gyrotropic_box_center', found, 3, r_value=gyrotropic_box_tmp)
    if (found) gyrotropic%box_corner(:) = &
      gyrotropic_box_tmp(:) - 0.5*(gyrotropic%box(1, :) + gyrotropic%box(2, :) + gyrotropic%box(3, :))

    call param_get_range_vector(stdout, seedname, 'gyrotropic_band_list', found, gyrotropic%num_bands, lcount=.true.)
    if (found) then
      if (gyrotropic%num_bands < 1) call io_error('Error: problem reading gyrotropic_band_list', stdout, seedname)
      if (allocated(gyrotropic%band_list)) deallocate (gyrotropic%band_list)
      allocate (gyrotropic%band_list(gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating gyrotropic_band_list in param_read', stdout, seedname)
      call param_get_range_vector(stdout, seedname, 'gyrotropic_band_list', found, &
                                  gyrotropic%num_bands, .false., gyrotropic%band_list)
      if (any(gyrotropic%band_list < 1) .or. any(gyrotropic%band_list > num_wann)) &
        call io_error('Error: gyrotropic_band_list asks for a non-valid bands', stdout, seedname)
    else
      ! include all bands in the calculation
      gyrotropic%num_bands = num_wann
      if (allocated(gyrotropic%band_list)) deallocate (gyrotropic%band_list)
      allocate (gyrotropic%band_list(gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating gyrotropic_band_list in param_read', stdout, seedname)
      do loop = 1, num_wann
        gyrotropic%band_list(loop) = loop
      end do
    end if

    smr_max_arg = 5.0
    call param_get_keyword(stdout, seedname, 'smr_max_arg', found, r_value=smr_max_arg)
    if (found .and. (smr_max_arg <= 0._dp)) &
      call io_error('Error: smr_max_arg must be greater than zero', stdout, seedname)

    gyrotropic%smr_max_arg = smr_max_arg
    call param_get_keyword(stdout, seedname, 'gyrotropic_smr_max_arg', found, &
                           r_value=gyrotropic%smr_max_arg)
    if (found .and. (gyrotropic%smr_max_arg <= 0._dp)) call io_error &
      ('Error: gyrotropic_smr_max_arg must be greater than zero', stdout, seedname)

    gyrotropic%smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword(stdout, seedname, 'gyrotropic_smr_fixed_en_width', found, &
                           r_value=gyrotropic%smr_fixed_en_width)
    if (found .and. (gyrotropic%smr_fixed_en_width < 0._dp)) call io_error &
      ('Error: gyrotropic_smr_fixed_en_width must be greater than or equal to zero', stdout, seedname)

    ! By default: use the "global" smearing index
    gyrotropic%smr_index = smr_index
    call param_get_keyword(stdout, seedname, 'gyrotropic_smr_type', found, c_value=ctmp)
    if (found) gyrotropic%smr_index = get_smearing_index(ctmp, 'gyrotropic_smr_type', stdout, seedname)

  end subroutine param_read_gyrotropic

  subroutine param_read_berry(pw90_calcs, berry, smearing, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_calculation_type), intent(in) :: pw90_calcs
    type(pw90_berry_mod_type), intent(out) :: berry
    type(pw90_smearing_type), intent(in) :: smearing
    character(len=50), intent(in)  :: seedname

    logical :: found
    character(len=maxlen)              :: ctmp

!-------------------------------------------------------
!    alpha=0
!    call param_get_keyword('alpha',found,i_value=alpha)

!    beta=0
!    call param_get_keyword('beta',found,i_value=beta)

!    gamma=0
!    call param_get_keyword('gamma',found,i_value=gamma)
!-------------------------------------------------------

    berry%transl_inv = .false.
    call param_get_keyword(stdout, seedname, 'transl_inv', found, l_value=berry%transl_inv)

    berry%task = ' '
    call param_get_keyword(stdout, seedname, 'berry_task', found, c_value=berry%task)
    if (pw90_calcs%berry .and. .not. found) call io_error &
      ('Error: berry=T and berry_task is not set', stdout, seedname)
    if (pw90_calcs%berry .and. index(berry%task, 'ahc') == 0 .and. index(berry%task, 'morb') == 0 &
        .and. index(berry%task, 'kubo') == 0 .and. index(berry%task, 'sc') == 0 &
        .and. index(berry%task, 'shc') == 0) call io_error &
      ('Error: value of berry_task not recognised in param_read', stdout, seedname)

    berry%curv_adpt_kmesh = 1
    call param_get_keyword(stdout, seedname, 'berry_curv_adpt_kmesh', found, &
                           i_value=berry%curv_adpt_kmesh)
    if (berry%curv_adpt_kmesh < 1) &
      call io_error( &
      'Error:  berry_curv_adpt_kmesh must be a positive integer', stdout, seedname)

    berry%curv_adpt_kmesh_thresh = 100.0_dp
    call param_get_keyword(stdout, seedname, 'berry_curv_adpt_kmesh_thresh', found, &
                           r_value=berry%curv_adpt_kmesh_thresh)

    berry%curv_unit = 'ang2'
    call param_get_keyword(stdout, seedname, 'berry_curv_unit', found, c_value=berry%curv_unit)
    if (berry%curv_unit .ne. 'ang2' .and. berry%curv_unit .ne. 'bohr2') &
      call io_error &
      ('Error: value of berry_curv_unit not recognised in param_read', stdout, seedname)

    berry%wanint_kpoint_file = .false.
    call param_get_keyword(stdout, seedname, 'wanint_kpoint_file', found, &
                           l_value=berry%wanint_kpoint_file)

!    smear_temp = -1.0_dp
!    call param_get_keyword('smear_temp',found,r_value=smear_temp)

    berry%kubo_smearing%use_adaptive = smearing%use_adaptive
    call param_get_keyword(stdout, seedname, 'kubo_adpt_smr', found, l_value=berry%kubo_smearing%use_adaptive)

    berry%kubo_smearing%adaptive_prefactor = smearing%adaptive_prefactor
    call param_get_keyword(stdout, seedname, 'kubo_adpt_smr_fac', found, &
                           r_value=berry%kubo_smearing%adaptive_prefactor)
    if (found .and. (berry%kubo_smearing%adaptive_prefactor <= 0._dp)) call io_error &
      ('Error: kubo_adpt_smr_fac must be greater than zero', stdout, seedname)

    berry%kubo_smearing%adaptive_max_width = smearing%adaptive_max_width
    call param_get_keyword(stdout, seedname, 'kubo_adpt_smr_max', found, &
                           r_value=berry%kubo_smearing%adaptive_max_width)
    if (berry%kubo_smearing%adaptive_max_width <= 0._dp) call io_error &
      ('Error: kubo_adpt_smr_max must be greater than zero', stdout, seedname)

    berry%kubo_smearing%fixed_width = smearing%fixed_width
    call param_get_keyword(stdout, seedname, 'kubo_smr_fixed_en_width', found, &
                           r_value=berry%kubo_smearing%fixed_width)
    if (found .and. (berry%kubo_smearing%fixed_width < 0._dp)) call io_error &
      ('Error: kubo_smr_fixed_en_width must be greater than or equal to zero', stdout, seedname)

    berry%sc_phase_conv = 1
    call param_get_keyword(stdout, seedname, 'sc_phase_conv', found, i_value=berry%sc_phase_conv)
    if ((berry%sc_phase_conv .ne. 1) .and. ((berry%sc_phase_conv .ne. 2))) &
      call io_error('Error: sc_phase_conv must be either 1 or 2', stdout, seedname)

    ! By default: use the "global" smearing index
    berry%kubo_smearing%type_index = smearing%type_index
    call param_get_keyword(stdout, seedname, 'kubo_smr_type', found, c_value=ctmp)
    if (found) berry%kubo_smearing%type_index = get_smearing_index(ctmp, 'kubo_smr_type', stdout, seedname)

    berry%sc_eta = 0.04
    call param_get_keyword(stdout, seedname, 'sc_eta', found, r_value=berry%sc_eta)

    berry%sc_w_thr = 5.0d0
    call param_get_keyword(stdout, seedname, 'sc_w_thr', found, r_value=berry%sc_w_thr)

  end subroutine param_read_berry

  subroutine param_read_spin_hall(pw90_calcs, pw90_common, spin_hall, berry_task, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_calculation_type), intent(in) :: pw90_calcs
    type(postw90_common_type), intent(in) :: pw90_common
    type(pw90_spin_hall_type), intent(out) :: spin_hall
    character(len=*), intent(in) :: berry_task
    character(len=50), intent(in)  :: seedname

    logical :: found

    spin_hall%freq_scan = .false.
    call param_get_keyword(stdout, seedname, 'shc_freq_scan', found, l_value=spin_hall%freq_scan)

    spin_hall%alpha = 1
    call param_get_keyword(stdout, seedname, 'shc_alpha', found, i_value=spin_hall%alpha)
    if (found .and. (spin_hall%alpha < 1 .or. spin_hall%alpha > 3)) call io_error &
      ('Error:  shc_alpha must be 1, 2 or 3', stdout, seedname)

    spin_hall%beta = 2
    call param_get_keyword(stdout, seedname, 'shc_beta', found, i_value=spin_hall%beta)
    if (found .and. (spin_hall%beta < 1 .or. spin_hall%beta > 3)) call io_error &
      ('Error:  shc_beta must be 1, 2 or 3', stdout, seedname)

    spin_hall%gamma = 3
    call param_get_keyword(stdout, seedname, 'shc_gamma', found, i_value=spin_hall%gamma)
    if (found .and. (spin_hall%gamma < 1 .or. spin_hall%gamma > 3)) call io_error &
      ('Error:  shc_gamma must be 1, 2 or 3', stdout, seedname)

    spin_hall%bandshift = .false.
    call param_get_keyword(stdout, seedname, 'shc_bandshift', found, l_value=spin_hall%bandshift)
    spin_hall%bandshift = spin_hall%bandshift .and. pw90_calcs%berry .and. .not. (index(berry_task, 'shc') == 0)
    if ((abs(pw90_common%scissors_shift) > 1.0e-7_dp) .and. spin_hall%bandshift) &
      call io_error('Error: shc_bandshift and scissors_shift cannot be used simultaneously', stdout, seedname)

    spin_hall%bandshift_firstband = 0
    call param_get_keyword(stdout, seedname, 'shc_bandshift_firstband', found, i_value=spin_hall%bandshift_firstband)
    if (spin_hall%bandshift .and. (.not. found)) &
      call io_error('Error: shc_bandshift required but no shc_bandshift_firstband provided', stdout, seedname)
    if ((spin_hall%bandshift_firstband < 1) .and. found) &
      call io_error('Error: shc_bandshift_firstband must >= 1', stdout, seedname)

    spin_hall%bandshift_energyshift = 0._dp
    call param_get_keyword(stdout, seedname, 'shc_bandshift_energyshift', found, r_value=spin_hall%bandshift_energyshift)
    if (spin_hall%bandshift .and. (.not. found)) &
      call io_error('Error: shc_bandshift required but no shc_bandshift_energyshift provided', stdout, seedname)

  end subroutine param_read_spin_hall

  subroutine param_read_pw90ham(pw90_ham, stdout, seedname)
!   use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_band_deriv_degen_type), intent(out) :: pw90_ham
    character(len=50), intent(in)  :: seedname

    logical :: found

    pw90_ham%use_degen_pert = .false.
    call param_get_keyword(stdout, seedname, 'use_degen_pert', found, l_value=pw90_ham%use_degen_pert)

    pw90_ham%degen_thr = 1.0d-4
    call param_get_keyword(stdout, seedname, 'degen_thr', found, r_value=pw90_ham%degen_thr)

  end subroutine param_read_pw90ham

  subroutine param_read_pw90_kpath(pw90_calcs, kpath, spec_points, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_calculation_type), intent(in) :: pw90_calcs
    type(pw90_kpath_mod_type), intent(out) :: kpath
    type(kpoint_path_type), intent(in) :: spec_points
    character(len=50), intent(in)  :: seedname

    logical :: found

    kpath%task = 'bands'
    call param_get_keyword(stdout, seedname, 'kpath_task', found, c_value=kpath%task)
    if (pw90_calcs%kpath .and. index(kpath%task, 'bands') == 0 .and. &
        index(kpath%task, 'curv') == 0 .and. &
        index(kpath%task, 'morb') == 0 .and. &
        index(kpath%task, 'shc') == 0) call io_error &
      ('Error: value of kpath_task not recognised in param_read', stdout, seedname)
    if (spec_points%bands_num_spec_points == 0 .and. pw90_calcs%kpath) &
      call io_error('Error: a kpath plot has been requested but there is no kpoint_path block', stdout, seedname)

    kpath%num_points = 100
    call param_get_keyword(stdout, seedname, 'kpath_num_points', found, &
                           i_value=kpath%num_points)
    if (kpath%num_points < 0) &
      call io_error('Error: kpath_num_points must be positive', stdout, seedname)

    kpath%bands_colour = 'none'
    call param_get_keyword(stdout, seedname, 'kpath_bands_colour', found, &
                           c_value=kpath%bands_colour)
    if (pw90_calcs%kpath .and. index(kpath%bands_colour, 'none') == 0 .and. &
        index(kpath%bands_colour, 'spin') == 0 .and. &
        index(kpath%bands_colour, 'shc') == 0) call io_error &
      ('Error: value of kpath_bands_colour not recognised in param_read', stdout, seedname)
    if (pw90_calcs%kpath .and. index(kpath%task, 'shc') > 0 .and. &
        index(kpath%task, 'spin') > 0) call io_error &
      ("Error: kpath_task cannot include both 'shc' and 'spin'", stdout, seedname)

  end subroutine param_read_pw90_kpath

  subroutine param_read_dos(pw90_calcs, dos_data, found_fermi_energy, num_wann, smearing, &
                            dos_plot, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_calculation_type), intent(in) :: pw90_calcs
    type(pw90_dos_mod_type), intent(out) :: dos_data
    logical, intent(in) :: found_fermi_energy
    integer, intent(in) :: num_wann
    type(pw90_smearing_type), intent(in) :: smearing
    logical, intent(out) :: dos_plot
    character(len=50), intent(in)  :: seedname

    integer :: i, ierr
    logical :: found
    character(len=maxlen)              :: ctmp

    dos_data%task = 'dos_plot'
    if (pw90_calcs%dos) then
      dos_plot = .true.
    else
      dos_plot = .false.
    endif
    call param_get_keyword(stdout, seedname, 'dos_task', found, c_value=dos_data%task)
    if (pw90_calcs%dos) then
      if (index(dos_data%task, 'dos_plot') == 0 .and. &
          index(dos_data%task, 'find_fermi_energy') == 0) call io_error &
        ('Error: value of dos_task not recognised in param_read', stdout, seedname)
      if (index(dos_data%task, 'dos_plot') > 0) dos_plot = .true.
      if (index(dos_data%task, 'find_fermi_energy') > 0 .and. found_fermi_energy) &
        call io_error &
        ('Error: Cannot set "dos_task = find_fermi_energy" and give a value to "fermi_energy"', stdout, seedname)
    end if

!    sigma_abc_onlyorb=.false.
!    call param_get_keyword('sigma_abc_onlyorb',found,l_value=sigma_abc_onlyorb)

! -------------------------------------------------------------------

    !IVO_END

    dos_data%energy_step = 0.01_dp
    call param_get_keyword(stdout, seedname, 'dos_energy_step', found, &
                           r_value=dos_data%energy_step)

    dos_data%smearing%use_adaptive = smearing%use_adaptive
    call param_get_keyword(stdout, seedname, 'dos_adpt_smr', found, &
                           l_value=dos_data%smearing%use_adaptive)

    dos_data%smearing%adaptive_prefactor = smearing%adaptive_prefactor
    call param_get_keyword(stdout, seedname, 'dos_adpt_smr_fac', found, &
                           r_value=dos_data%smearing%adaptive_prefactor)
    if (found .and. (dos_data%smearing%adaptive_prefactor <= 0._dp)) &
      call io_error('Error: dos_adpt_smr_fac must be greater than zero', stdout, seedname)

    dos_data%smearing%adaptive_max_width = smearing%adaptive_max_width
    call param_get_keyword(stdout, seedname, 'dos_adpt_smr_max', found, &
                           r_value=dos_data%smearing%adaptive_max_width)
    if (dos_data%smearing%adaptive_max_width <= 0._dp) call io_error &
      ('Error: dos_adpt_smr_max must be greater than zero', stdout, seedname)

    dos_data%smearing%fixed_width = smearing%fixed_width
    call param_get_keyword(stdout, seedname, 'dos_smr_fixed_en_width', found, &
                           r_value=dos_data%smearing%fixed_width)
    if (found .and. (dos_data%smearing%fixed_width < 0._dp)) &
      call io_error('Error: dos_smr_fixed_en_width must be greater than or equal to zero', stdout, seedname)

!    dos_gaussian_width        = 0.1_dp
!    call param_get_keyword('dos_gaussian_width',found,r_value=dos_gaussian_width)

!    dos_plot_format           = 'gnuplot'
!    call param_get_keyword('dos_plot_format',found,c_value=dos_plot_format)

    call param_get_range_vector(stdout, seedname, 'dos_project', found, &
                                dos_data%num_project, lcount=.true.)
    if (found) then
      if (dos_data%num_project < 1) call io_error('Error: problem reading dos_project', stdout, seedname)
      if (allocated(dos_data%project)) deallocate (dos_data%project)
      allocate (dos_data%project(dos_data%num_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dos_project in param_read', stdout, seedname)
      call param_get_range_vector(stdout, seedname, 'dos_project', found, &
                                  dos_data%num_project, .false., &
                                  dos_data%project)
      if (any(dos_data%project < 1) .or. &
          any(dos_data%project > num_wann)) &
        call io_error('Error: dos_project asks for out-of-range Wannier functions', stdout, seedname)
    else
      ! by default plot all
      dos_data%num_project = num_wann
      if (allocated(dos_data%project)) deallocate (dos_data%project)
      allocate (dos_data%project(dos_data%num_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dos_project in param_read', stdout, seedname)
      do i = 1, dos_data%num_project
        dos_data%project(i) = i
      end do
    endif

    ! By default: use the "global" smearing index
    dos_data%smearing%type_index = smearing%type_index
    call param_get_keyword(stdout, seedname, 'dos_smr_type', found, c_value=ctmp)
    if (found) dos_data%smearing%type_index = get_smearing_index(ctmp, 'dos_smr_type', stdout, seedname)

  end subroutine param_read_dos

  subroutine param_read_geninterp(geninterp, stdout, seedname)
    ! [gp-begin, Jun 1, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! General band interpolator (geninterp)
    !%%%%%%%%%%%%%%%%%%%%
!   use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_geninterp_mod_type), intent(out) :: geninterp
    character(len=50), intent(in)  :: seedname

    logical :: found

    geninterp%alsofirstder = .false.
    call param_get_keyword(stdout, seedname, 'geninterp_alsofirstder', found, l_value=geninterp%alsofirstder)
    geninterp%single_file = .true.
    call param_get_keyword(stdout, seedname, 'geninterp_single_file', found, l_value=geninterp%single_file)
    ! [gp-end, Jun 1, 2012]

  end subroutine param_read_geninterp

  subroutine param_read_boltzwann(boltz, eigval, smearing, do_boltzwann, boltz_2d_dir, &
                                  stdout, seedname)
    ! [gp-begin, Jun 1, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! General band interpolator (geninterp)
    !%%%%%%%%%%%%%%%%%%%%
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_boltzwann_type), intent(inout) :: boltz
    real(kind=dp), allocatable, intent(in) :: eigval(:, :)
    type(pw90_smearing_type), intent(in) :: smearing
    logical, intent(in) :: do_boltzwann
    character(len=4), intent(out) :: boltz_2d_dir
    character(len=50), intent(in)  :: seedname

    logical :: found, found2
    character(len=maxlen)              :: ctmp

    ! [gp-begin, Apr 12, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! Boltzmann transport
    !%%%%%%%%%%%%%%%%%%%%
    ! Note: to be put AFTER the disentanglement routines!

    boltz%calc_also_dos = .false.
    call param_get_keyword(stdout, seedname, 'boltz_calc_also_dos', found, l_value=boltz%calc_also_dos)

    boltz%calc_also_dos = boltz%calc_also_dos .and. do_boltzwann

    ! 0 means the normal 3d case for the calculation of the Seebeck coefficient
    ! The other valid possibilities are 1,2,3 for x,y,z respectively
    boltz%dir_num_2d = 0
    call param_get_keyword(stdout, seedname, 'boltz_2d_dir', found, c_value=boltz_2d_dir)
    if (found) then
      if (trim(boltz_2d_dir) == 'no') then
        boltz%dir_num_2d = 0
      elseif (trim(boltz_2d_dir) == 'x') then
        boltz%dir_num_2d = 1
      elseif (trim(boltz_2d_dir) == 'y') then
        boltz%dir_num_2d = 2
      elseif (trim(boltz_2d_dir) == 'z') then
        boltz%dir_num_2d = 3
      else
        call io_error('Error: boltz_2d_dir can only be "no", "x", "y" or "z".', stdout, seedname)
      end if
    end if

    boltz%dos_energy_step = 0.001_dp
    call param_get_keyword(stdout, seedname, 'boltz_dos_energy_step', found, r_value=boltz%dos_energy_step)
    if (found .and. (boltz%dos_energy_step <= 0._dp)) &
      call io_error('Error: boltz_dos_energy_step must be positive', stdout, seedname)

    if (allocated(eigval)) then
      boltz%dos_energy_min = minval(eigval) - 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      boltz%dos_energy_min = -1.0_dp
    end if
    call param_get_keyword(stdout, seedname, 'boltz_dos_energy_min', found, r_value=boltz%dos_energy_min)
    if (allocated(eigval)) then
      boltz%dos_energy_max = maxval(eigval) + 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      boltz%dos_energy_max = 0.0_dp
    end if
    call param_get_keyword(stdout, seedname, 'boltz_dos_energy_max', found, r_value=boltz%dos_energy_max)
    if (boltz%dos_energy_max <= boltz%dos_energy_min) &
      call io_error('Error: boltz_dos_energy_max must be greater than boltz_dos_energy_min', stdout, seedname)

    boltz%dos_smearing%use_adaptive = smearing%use_adaptive
    call param_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr', found, &
                           l_value=boltz%dos_smearing%use_adaptive)

    boltz%dos_smearing%adaptive_prefactor = smearing%adaptive_prefactor
    call param_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr_fac', found, &
                           r_value=boltz%dos_smearing%adaptive_prefactor)
    if (found .and. (boltz%dos_smearing%adaptive_prefactor <= 0._dp)) &
      call io_error('Error: boltz_dos_adpt_smr_fac must be greater than zero', stdout, seedname)

    boltz%dos_smearing%adaptive_max_width = smearing%adaptive_max_width
    call param_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr_max', found, r_value=boltz%dos_smearing%adaptive_max_width)
    if (boltz%dos_smearing%adaptive_max_width <= 0._dp) call io_error &
      ('Error: boltz_dos_adpt_smr_max must be greater than zero', stdout, seedname)

    boltz%dos_smearing%fixed_width = smearing%fixed_width
    call param_get_keyword(stdout, seedname, 'boltz_dos_smr_fixed_en_width', found, r_value=boltz%dos_smearing%fixed_width)
    if (found .and. (boltz%dos_smearing%fixed_width < 0._dp)) &
      call io_error('Error: boltz_dos_smr_fixed_en_width must be greater than or equal to zero', stdout, seedname)

    boltz%mu_min = -999._dp
    call param_get_keyword(stdout, seedname, 'boltz_mu_min', found, r_value=boltz%mu_min)
    if ((.not. found) .and. do_boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_min provided', stdout, seedname)
    boltz%mu_max = -999._dp
    call param_get_keyword(stdout, seedname, 'boltz_mu_max', found2, r_value=boltz%mu_max)
    if ((.not. found2) .and. do_boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_max provided', stdout, seedname)
    if (found .and. found2 .and. (boltz%mu_max < boltz%mu_min)) &
      call io_error('Error: boltz_mu_max must be greater than boltz_mu_min', stdout, seedname)
    boltz%mu_step = 0._dp
    call param_get_keyword(stdout, seedname, 'boltz_mu_step', found, r_value=boltz%mu_step)
    if ((.not. found) .and. do_boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_step provided', stdout, seedname)
    if (found .and. (boltz%mu_step <= 0._dp)) &
      call io_error('Error: boltz_mu_step must be greater than zero', stdout, seedname)

    boltz%temp_min = -999._dp
    call param_get_keyword(stdout, seedname, 'boltz_temp_min', found, r_value=boltz%temp_min)
    if ((.not. found) .and. do_boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_min provided', stdout, seedname)
    boltz%temp_max = -999._dp
    call param_get_keyword(stdout, seedname, 'boltz_temp_max', found2, r_value=boltz%temp_max)
    if ((.not. found2) .and. do_boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_max provided', stdout, seedname)
    if (found .and. found2 .and. (boltz%temp_max < boltz%temp_min)) &
      call io_error('Error: boltz_temp_max must be greater than boltz_temp_min', stdout, seedname)
    if (found .and. (boltz%temp_min <= 0._dp)) &
      call io_error('Error: boltz_temp_min must be greater than zero', stdout, seedname)
    boltz%temp_step = 0._dp
    call param_get_keyword(stdout, seedname, 'boltz_temp_step', found, r_value=boltz%temp_step)
    if ((.not. found) .and. do_boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_step provided', stdout, seedname)
    if (found .and. (boltz%temp_step <= 0._dp)) &
      call io_error('Error: boltz_temp_step must be greater than zero', stdout, seedname)

    ! The interpolation mesh is read later on

    ! By default, the energy step for the TDF is 1 meV
    boltz%tdf_energy_step = 0.001_dp
    call param_get_keyword(stdout, seedname, 'boltz_tdf_energy_step', found, r_value=boltz%tdf_energy_step)
    if (boltz%tdf_energy_step <= 0._dp) &
      call io_error('Error: boltz_tdf_energy_step must be greater than zero', stdout, seedname)

    ! For TDF: TDF smeared in a NON-adaptive way; value in eV, default = 0._dp
    ! (i.e., no smearing)
    boltz%TDF_smr_fixed_en_width = smearing%fixed_width
    call param_get_keyword(stdout, seedname, 'boltz_tdf_smr_fixed_en_width', found, r_value=boltz%TDF_smr_fixed_en_width)
    if (found .and. (boltz%TDF_smr_fixed_en_width < 0._dp)) &
      call io_error('Error: boltz_TDF_smr_fixed_en_width must be greater than or equal to zero', stdout, seedname)

    ! By default: use the "global" smearing index
    boltz%TDF_smr_index = smearing%type_index
    call param_get_keyword(stdout, seedname, 'boltz_tdf_smr_type', found, c_value=ctmp)
    if (found) boltz%TDF_smr_index = get_smearing_index(ctmp, 'boltz_tdf_smr_type', stdout, seedname)

    ! By default: use the "global" smearing index
    boltz%dos_smearing%type_index = smearing%type_index
    call param_get_keyword(stdout, seedname, 'boltz_dos_smr_type', found, c_value=ctmp)
    if (found) boltz%dos_smearing%type_index = get_smearing_index(ctmp, 'boltz_dos_smr_type', stdout, seedname)

    ! By default: 10 fs relaxation time
    boltz%relax_time = 10._dp
    call param_get_keyword(stdout, seedname, 'boltz_relax_time', found, r_value=boltz%relax_time)

    boltz%bandshift = .false.
    call param_get_keyword(stdout, seedname, 'boltz_bandshift', found, l_value=boltz%bandshift)
    boltz%bandshift = boltz%bandshift .and. do_boltzwann

    boltz%bandshift_firstband = 0
    call param_get_keyword(stdout, seedname, 'boltz_bandshift_firstband', found, i_value=boltz%bandshift_firstband)
    if (boltz%bandshift .and. (.not. found)) &
      call io_error('Error: boltz_bandshift required but no boltz_bandshift_firstband provided', stdout, seedname)
    boltz%bandshift_energyshift = 0._dp
    call param_get_keyword(stdout, seedname, 'boltz_bandshift_energyshift', found, r_value=boltz%bandshift_energyshift)
    if (boltz%bandshift .and. (.not. found)) &
      call io_error('Error: boltz_bandshift required but no boltz_bandshift_energyshift provided', stdout, seedname)
    ! [gp-end, Apr 12, 2012]
  end subroutine param_read_boltzwann

  subroutine param_read_energy_range(berry, dos_data, gyrotropic, dis_window, fermi, eigval, &
                                     write_data, stdout, seedname)
    use w90_constants, only: cmplx_i
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(pw90_berry_mod_type), intent(inout) :: berry
    type(pw90_dos_mod_type), intent(inout) :: dos_data
    type(pw90_gyrotropic_type), intent(inout) :: gyrotropic
    type(dis_manifold_type), intent(in) :: dis_window
    type(fermi_data_type), intent(in) :: fermi
    real(kind=dp), allocatable, intent(in) :: eigval(:, :)
    type(pw90_extra_io_type), intent(inout) :: write_data
    character(len=50), intent(in)  :: seedname

    integer :: i, ierr
    logical :: found

    if (dis_window%frozen_states) then
      dos_data%energy_max = dis_window%froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      dos_data%energy_max = maxval(eigval) + 0.6667_dp
    else
      dos_data%energy_max = dis_window%win_max + 0.6667_dp
    end if
    call param_get_keyword(stdout, seedname, 'dos_energy_max', found, &
                           r_value=dos_data%energy_max)

    if (allocated(eigval)) then
      dos_data%energy_min = minval(eigval) - 0.6667_dp
    else
      dos_data%energy_min = dis_window%win_min - 0.6667_dp
    end if
    call param_get_keyword(stdout, seedname, 'dos_energy_min', found, &
                           r_value=dos_data%energy_min)

    write_data%kubo_freq_min = 0.0_dp
    write_data%gyrotropic_freq_min = write_data%kubo_freq_min
    call param_get_keyword(stdout, seedname, 'kubo_freq_min', found, &
                           r_value=write_data%kubo_freq_min)
    !
    if (dis_window%frozen_states) then
      write_data%kubo_freq_max = dis_window%froz_max - fermi%energy_list(1) + 0.6667_dp
    elseif (allocated(eigval)) then
      write_data%kubo_freq_max = maxval(eigval) - minval(eigval) + 0.6667_dp
    else
      write_data%kubo_freq_max = dis_window%win_max - dis_window%win_min + 0.6667_dp
    end if
    write_data%gyrotropic_freq_max = write_data%kubo_freq_max
    call param_get_keyword(stdout, seedname, 'kubo_freq_max', found, &
                           r_value=write_data%kubo_freq_max)

    !
    write_data%kubo_freq_step = 0.01_dp
    call param_get_keyword(stdout, seedname, 'kubo_freq_step', found, &
                           r_value=write_data%kubo_freq_step)
    if (found .and. write_data%kubo_freq_step < 0.0_dp) call io_error( &
      'Error: kubo_freq_step must be positive', stdout, seedname)
    !
    berry%kubo_nfreq = nint((write_data%kubo_freq_max - write_data%kubo_freq_min) &
                            /write_data%kubo_freq_step) + 1
    if (berry%kubo_nfreq <= 1) berry%kubo_nfreq = 2
    write_data%kubo_freq_step = (write_data%kubo_freq_max - write_data%kubo_freq_min) &
                                /(berry%kubo_nfreq - 1)
    !
    if (allocated(berry%kubo_freq_list)) deallocate (berry%kubo_freq_list)
    allocate (berry%kubo_freq_list(berry%kubo_nfreq), stat=ierr)
    if (ierr /= 0) &
      call io_error('Error allocating kubo_freq_list in param_read', stdout, seedname)
    do i = 1, berry%kubo_nfreq
      berry%kubo_freq_list(i) = write_data%kubo_freq_min &
                                + (i - 1)*(write_data%kubo_freq_max - write_data%kubo_freq_min) &
                                /(berry%kubo_nfreq - 1)
    enddo
    !
    ! TODO: Alternatively, read list of (complex) frequencies; kubo_nfreq is
    !       the length of the list

    write_data%gyrotropic_freq_step = 0.01_dp
    call param_get_keyword(stdout, seedname, 'gyrotropic_freq_min', found, &
                           r_value=write_data%gyrotropic_freq_min)
    call param_get_keyword(stdout, seedname, 'gyrotropic_freq_max', found, &
                           r_value=write_data%gyrotropic_freq_max)
    call param_get_keyword(stdout, seedname, 'gyrotropic_freq_step', found, &
                           r_value=write_data%gyrotropic_freq_step)
    gyrotropic%nfreq = nint((write_data%gyrotropic_freq_max - write_data%gyrotropic_freq_min) &
                            /write_data%gyrotropic_freq_step) + 1
    if (gyrotropic%nfreq <= 1) gyrotropic%nfreq = 2
    write_data%gyrotropic_freq_step = (write_data%gyrotropic_freq_max &
                                       - write_data%gyrotropic_freq_min)/(gyrotropic%nfreq - 1)
    if (allocated(gyrotropic%freq_list)) deallocate (gyrotropic%freq_list)
    allocate (gyrotropic%freq_list(gyrotropic%nfreq), stat=ierr)
    if (ierr /= 0) &
      call io_error('Error allocating gyrotropic_freq_list in param_read', stdout, seedname)
    do i = 1, gyrotropic%nfreq
      gyrotropic%freq_list(i) = write_data%gyrotropic_freq_min &
                                + (i - 1)*(write_data%gyrotropic_freq_max &
                                           - write_data%gyrotropic_freq_min)/(gyrotropic%nfreq - 1) &
                                + cmplx_i*gyrotropic%smr_fixed_en_width
    enddo

    if (dis_window%frozen_states) then
      berry%kubo_eigval_max = dis_window%froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      berry%kubo_eigval_max = maxval(eigval) + 0.6667_dp
    else
      berry%kubo_eigval_max = dis_window%win_max + 0.6667_dp
    end if
    gyrotropic%eigval_max = berry%kubo_eigval_max

    call param_get_keyword(stdout, seedname, 'kubo_eigval_max', found, r_value=berry%kubo_eigval_max)
    call param_get_keyword(stdout, seedname, 'gyrotropic_eigval_max', found, r_value=gyrotropic%eigval_max)
  end subroutine param_read_energy_range

  subroutine param_read_local_kmesh(pw90_calcs, berry, dos_data, pw90_spin, &
                                    gyrotropic, boltz, recip_lattice, global_kmesh_set, kmesh, &
                                    kmesh_spacing, stdout, seedname)
    implicit none
    integer, intent(in) :: stdout
    type(pw90_calculation_type), intent(in) :: pw90_calcs
    type(pw90_berry_mod_type), intent(inout) :: berry
    type(pw90_dos_mod_type), intent(inout) :: dos_data
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    type(pw90_gyrotropic_type), intent(inout) :: gyrotropic
    type(pw90_boltzwann_type), intent(inout) :: boltz
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    logical, intent(in) :: global_kmesh_set
    real(kind=dp), intent(in) :: kmesh_spacing
    integer, intent(in) :: kmesh(3)
    character(len=50), intent(in)  :: seedname

    ! To be called after having read the global flag
    call get_module_kmesh(stdout, seedname, recip_lattice, global_kmesh_set, kmesh, kmesh_spacing, &
                          moduleprefix='boltz', should_be_defined=pw90_calcs%boltzwann, &
                          module_kmesh=boltz%kmesh, &
                          module_kmesh_spacing=boltz%kmesh_spacing)

    call get_module_kmesh(stdout, seedname, recip_lattice, global_kmesh_set, kmesh, kmesh_spacing, &
                          moduleprefix='berry', should_be_defined=pw90_calcs%berry, &
                          module_kmesh=berry%kmesh, &
                          module_kmesh_spacing=berry%kmesh_spacing)

    call get_module_kmesh(stdout, seedname, recip_lattice, global_kmesh_set, kmesh, kmesh_spacing, &
                          moduleprefix='gyrotropic', should_be_defined=pw90_calcs%gyrotropic, &
                          module_kmesh=gyrotropic%kmesh, &
                          module_kmesh_spacing=gyrotropic%kmesh_spacing)

    call get_module_kmesh(stdout, seedname, recip_lattice, global_kmesh_set, kmesh, kmesh_spacing, &
                          moduleprefix='spin', should_be_defined=pw90_calcs%spin_moment, &
                          module_kmesh=pw90_spin%kmesh, &
                          module_kmesh_spacing=pw90_spin%kmesh_spacing)

    call get_module_kmesh(stdout, seedname, recip_lattice, global_kmesh_set, kmesh, kmesh_spacing, &
                          moduleprefix='dos', should_be_defined=pw90_calcs%dos, &
                          module_kmesh=dos_data%kmesh, &
                          module_kmesh_spacing=dos_data%kmesh_spacing)
  end subroutine param_read_local_kmesh

  subroutine get_module_kmesh(stdout, seedname, recip_lattice, global_kmesh_set, kmesh, &
                              kmesh_spacing, moduleprefix, should_be_defined, module_kmesh, &
                              module_kmesh_spacing)
    !! This function reads and sets the interpolation mesh variables needed by a given module
    !>
    !!  This function MUST be called after having read the global kmesh and kmesh_spacing!!
    !!  if the user didn't provide an interpolation_mesh_spacing, it is set to -1, so that
    !!       one can check in the code what the user asked for
    !!  The function takes care also of setting the default value to the global one if no local
    !!       keyword is defined
    use w90_io, only: io_error
    integer, intent(in) :: stdout
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
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
    logical, intent(in) :: global_kmesh_set
    real(kind=dp), intent(in) :: kmesh_spacing
    integer, intent(in) :: kmesh(3)
    character(len=50), intent(in)  :: seedname

    logical :: found, found2
    integer :: i

    ! Default values
    module_kmesh_spacing = -1._dp
    module_kmesh = 0
    call param_get_keyword(stdout, seedname, trim(moduleprefix)//'_kmesh_spacing', found, r_value=module_kmesh_spacing)
    if (found) then
      if (module_kmesh_spacing .le. 0._dp) &
        call io_error('Error: '//trim(moduleprefix)//'_kmesh_spacing must be greater than zero', stdout, seedname)

      call internal_set_kmesh(module_kmesh_spacing, recip_lattice, module_kmesh)
    end if
    call param_get_vector_length(stdout, seedname, trim(moduleprefix)//'_kmesh', found2, length=i)
    if (found2) then
      if (found) &
        call io_error('Error: cannot set both '//trim(moduleprefix)//'_kmesh and ' &
                      //trim(moduleprefix)//'_kmesh_spacing', stdout, seedname)
      if (i .eq. 1) then
        call param_get_keyword_vector(stdout, seedname, trim(moduleprefix)//'_kmesh', found2, 1, i_value=module_kmesh)
        module_kmesh(2) = module_kmesh(1)
        module_kmesh(3) = module_kmesh(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector(stdout, seedname, trim(moduleprefix)//'_kmesh', found2, 3, i_value=module_kmesh)
      else
        call io_error('Error: '//trim(moduleprefix)// &
                      '_kmesh must be provided as either one integer or a vector of 3 integers', stdout, seedname)
      end if
      if (any(module_kmesh <= 0)) &
        call io_error('Error: '//trim(moduleprefix)//'_kmesh elements must be greater than zero', stdout, seedname)
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
          call io_error('Error: '//trim(moduleprefix)//' module required, but no interpolation mesh given.', stdout, seedname)
      end if
    end if
  end subroutine get_module_kmesh

!===================================================================
  subroutine param_postw90_write(param_input, system, fermi, atoms, num_wann, &
                                 real_lattice, recip_lattice, spec_points, &
                                 pw90_calcs, postw90_oper, pw90_common, &
                                 pw90_spin, kpath, kslice, dos_data, berry, &
                                 gyrotropic, geninterp, boltz, write_data, stdout)
    !==================================================================!
    !                                                                  !
    !! write postw90 parameters to stdout
    !                                                                  !
    !===================================================================

    implicit none

    !data from parameters module
    type(print_output_type), intent(in) :: param_input
    type(w90_system_type), intent(in) :: system
    type(fermi_data_type), intent(in) :: fermi
    type(atom_data_type), intent(in) :: atoms
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    type(kpoint_path_type), intent(in) :: spec_points
    type(pw90_calculation_type), intent(in) :: pw90_calcs
    type(pw90_oper_read_type), intent(in) :: postw90_oper
    type(postw90_common_type), intent(in) :: pw90_common
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(pw90_kpath_mod_type), intent(in) :: kpath
    type(pw90_kslice_mod_type), intent(in) :: kslice
    type(pw90_dos_mod_type), intent(in) :: dos_data
    type(pw90_berry_mod_type), intent(in) :: berry
    type(pw90_gyrotropic_type), intent(in) :: gyrotropic
    type(pw90_geninterp_mod_type), intent(in) :: geninterp
    type(pw90_boltzwann_type), intent(in) :: boltz
    type(pw90_extra_io_type), intent(in) :: write_data

    integer :: i, loop, nat, nsp
    real(kind=dp) :: cell_volume

    ! System
    write (stdout, *)
    write (stdout, '(36x,a6)') '------'
    write (stdout, '(36x,a6)') 'SYSTEM'
    write (stdout, '(36x,a6)') '------'
    write (stdout, *)
    if (param_input%lenconfac .eq. 1.0_dp) then
      write (stdout, '(30x,a21)') 'Lattice Vectors (Ang)'
    else
      write (stdout, '(28x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write (stdout, 101) 'a_1', (real_lattice(1, I)*param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'a_2', (real_lattice(2, I)*param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'a_3', (real_lattice(3, I)*param_input%lenconfac, i=1, 3)
    write (stdout, *)
    cell_volume = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - real_lattice(3, 2)*real_lattice(2, 3)) + &
                  real_lattice(1, 2)*(real_lattice(2, 3)*real_lattice(3, 1) - real_lattice(3, 3)*real_lattice(2, 1)) + &
                  real_lattice(1, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - real_lattice(3, 1)*real_lattice(2, 2))
    write (stdout, '(19x,a17,3x,f11.5)', advance='no') &
      'Unit Cell Volume:', cell_volume*param_input%lenconfac**3
    if (param_input%lenconfac .eq. 1.0_dp) then
      write (stdout, '(2x,a7)') '(Ang^3)'
    else
      write (stdout, '(2x,a8)') '(Bohr^3)'
    endif
    write (stdout, *)
    if (param_input%lenconfac .eq. 1.0_dp) then
      write (stdout, '(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
      write (stdout, '(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write (stdout, 101) 'b_1', (recip_lattice(1, I)/param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'b_2', (recip_lattice(2, I)/param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'b_3', (recip_lattice(3, I)/param_input%lenconfac, i=1, 3)
    write (stdout, *) ' '
    ! Atoms
    if (atoms%num_atoms > 0) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (param_input%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
      else
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, atoms%num_species
        do nat = 1, atoms%species_num(nsp)
          write (stdout, '(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                 '|', atoms%symbol(nsp), nat, atoms%pos_frac(:, nat, nsp),&
  &                 '|', atoms%pos_cart(:, nat, nsp)*param_input%lenconfac, '|'
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
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of electrons per state             :', &
      system%num_elec_per_state, '|'
    if (abs(pw90_common%scissors_shift) > 1.0e-7_dp .or. param_input%iprint > 0) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Scissor shift applied to conduction bands :', pw90_common%scissors_shift, '|'
      if (system%num_valence_bands > 0) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Number of valence bands                   :', &
          system%num_valence_bands, '|'
      else
        write (stdout, '(1x,a78)') '|  Number of valence bands                   :       not defined             |'
      endif
    endif
    if (pw90_calcs%spin_decomp .or. param_input%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Spin decomposition                        :', pw90_calcs%spin_decomp, '|'
    if (pw90_calcs%spin_moment .or. param_input%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Spin moment                       :', pw90_calcs%spin_moment, '|'
    if (pw90_calcs%spin_decomp .or. pw90_calcs%spin_moment .or. param_input%iprint > 2) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Polar angle of spin quantisation axis     :', pw90_spin%axis_polar, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Azimuthal angle of spin quantisation axis :', pw90_spin%axis_azimuth, '|'
      if (postw90_oper%spn_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  Spn file-type                   :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  Spn file-type                   :', 'unformatted', '|'
      endif
      if (postw90_oper%uHu_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  uHu file-type                   :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  uHu file-type                   :', 'unformatted', '|'
      endif
    end if

    if (size(fermi%energy_list) == 1) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi energy (eV)                         :', fermi%energy_list(1), '|'
    else
      write (stdout, '(1x,a21,I8,a12,f8.3,a4,f8.3,a3,13x,a1)') '|  Fermi energy     :', size(fermi%energy_list), &
        ' steps from ', fermi%energy_list(1), ' to ', &
        fermi%energy_list(size(fermi%energy_list)), ' eV', '|'
    end if

    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :', param_input%iprint, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :', param_input%timing_level, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Optimisation (0=memory, 3=speed)          :', param_input%optimisation, '|'
    write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :', trim(param_input%length_unit), '|'
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    write (stdout, '(1x,a78)') '*------------------------ Global Smearing Parameters ------------------------*'
    if (write_data%smear%use_adaptive) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', &
        write_data%smear%adaptive_prefactor, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width (eV)       :', &
        write_data%smear%adaptive_max_width, '|'

    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
        write_data%smear%fixed_width, '|'
    endif
    write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
      trim(param_get_smearing_type(write_data%smear%type_index)), '|'
    if (write_data%global_kmesh_set) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       T', '|'
      if (write_data%kmesh_spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          write_data%kmesh(1), 'x', write_data%kmesh(2), 'x', write_data%kmesh(3), ' Spacing = ', &
          write_data%kmesh_spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , write_data%kmesh(1), 'x', write_data%kmesh(2), 'x', write_data%kmesh(3), '|'
      endif
    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       F', '|'
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    ! DOS
    if (pw90_calcs%dos .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*---------------------------------- DOS -------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Density of States                :', pw90_calcs%dos, '|'
      if (dos_data%num_project > 1) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', dos_data%energy_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', dos_data%energy_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', dos_data%energy_step, '|'
      if (dos_data%smearing%use_adaptive .eqv. write_data%smear%use_adaptive .and. &
          dos_data%smearing%adaptive_prefactor == write_data%smear%adaptive_prefactor .and. &
          dos_data%smearing%adaptive_max_width == write_data%smear%adaptive_max_width .and. &
          dos_data%smearing%fixed_width == write_data%smear%fixed_width .and. &
          write_data%smear%type_index == dos_data%smearing%type_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (dos_data%smearing%use_adaptive) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', &
            dos_data%smearing%adaptive_prefactor, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', &
            dos_data%smearing%adaptive_max_width, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
            dos_data%smearing%fixed_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
          trim(param_get_smearing_type(dos_data%smearing%type_index)), '|'
      endif
      if (write_data%kmesh(1) == dos_data%kmesh(1) .and. &
          write_data%kmesh(2) == dos_data%kmesh(2) .and. &
          write_data%kmesh(3) == dos_data%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (dos_data%kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            dos_data%kmesh(1), 'x', dos_data%kmesh(2), 'x', &
            dos_data%kmesh(3), ' Spacing = ', dos_data%kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :', &
            dos_data%kmesh(1), 'x', dos_data%kmesh(2), 'x', &
            dos_data%kmesh(3), '|'
        endif
      endif
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    if (pw90_calcs%kpath .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KPATH ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a path in k-space   :', pw90_calcs%kpath, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', kpath%num_points, '|'
      if (index(kpath%task, 'bands') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       F', '|'
      endif
      if (index(kpath%task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       F', '|'
      endif
      if (index(kpath%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(kpath%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the bands    :', trim(kpath%bands_colour), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      write (stdout, '(1x,a78)') '|   K-space path sections:                                                   |'
      if (spec_points%bands_num_spec_points == 0) then
        write (stdout, '(1x,a78)') '|     None defined                                                           |'
      else
        do loop = 1, spec_points%bands_num_spec_points, 2
          write (stdout, '(1x,a10,2x,a1,2x,3F7.3,5x,a3,2x,a1,2x,3F7.3,7x,a1)') '|    From:', &
            spec_points%labels(loop), (spec_points%points(i, loop), i=1, 3), &
            'To:', spec_points%labels(loop + 1), (spec_points%points(i, loop + 1), i=1, 3), '|'
        end do
      end if
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%kslice .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KSLICE -----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a slice in k-space  :', pw90_calcs%kslice, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi level used for slice                :', fermi%energy_list(1), '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', kpath%num_points, '|'
      if (index(kslice%task, 'fermi_lines') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       F', '|'
      endif
      if (index(kslice%task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       F', '|'
      endif
      if (index(kslice%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(kslice%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the lines    :', &
        trim(kslice%fermi_lines_colour), '|'
      write (stdout, '(1x,a78)') '|  2D slice parameters (in reduced coordinates):                             |'
      write (stdout, '(1x,a14,2x,3F8.3,37x,a1)') '|     Corner: ', (kslice%corner(i), i=1, 3), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector1: ', (kslice%b1(i), i=1, 3), ' Divisions:', kslice%kmesh2d(1), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector2: ', (kslice%b2(i), i=1, 3), ' Divisions:', kslice%kmesh2d(1), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%berry .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- BERRY ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Berry Phase related properties    :', pw90_calcs%berry, '|'
      if (index(berry%task, 'kubo') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       F', '|'
      endif
      if (index(berry%task, 'ahc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       F', '|'
      endif
      if (index(berry%task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       F', '|'
      endif
      if (index(berry%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       F', '|'
      endif
      if (index(berry%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for optical responses     :', &
        write_data%kubo_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency for optical responses     :', &
        write_data%kubo_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for optical responses           :', &
        write_data%kubo_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue for optical responses    :', berry%kubo_eigval_max, '|'
      if (index(berry%task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing factor for shift current         :', berry%sc_eta, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Frequency theshold for shift current      :', berry%sc_w_thr, '|'
        write (stdout, '(1x,a46,1x,a27,3x,a1)') '|  Bloch sums                                :', &
          trim(param_get_convention_type(berry%sc_phase_conv)), '|'
      end if
      if (berry%kubo_smearing%use_adaptive .eqv. write_data%smear%use_adaptive .and. &
          berry%kubo_smearing%adaptive_prefactor == write_data%smear%adaptive_prefactor .and. &
          berry%kubo_smearing%adaptive_max_width == write_data%smear%adaptive_max_width &
          .and. berry%kubo_smearing%fixed_width == write_data%smear%fixed_width .and. &
          write_data%smear%type_index == berry%kubo_smearing%type_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (berry%kubo_smearing%use_adaptive) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', &
            berry%kubo_smearing%adaptive_prefactor, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', &
            berry%kubo_smearing%adaptive_max_width, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
            berry%kubo_smearing%fixed_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
          trim(param_get_smearing_type(berry%kubo_smearing%type_index)), '|'
      endif
      if (write_data%kmesh(1) == berry%kmesh(1) .and. write_data%kmesh(2) == berry%kmesh(2) .and. &
          write_data%kmesh(3) == berry%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (berry%kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            berry%kmesh(1), 'x', berry%kmesh(2), 'x', berry%kmesh(3), ' Spacing = ', berry%kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , berry%kmesh(1), 'x', berry%kmesh(2), 'x', berry%kmesh(3), '|'
        endif
      endif
      if (berry%curv_adpt_kmesh > 1) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Using an adaptive refinement mesh of size :', berry%curv_adpt_kmesh, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Threshold for adaptive refinement         :', &
          berry%curv_adpt_kmesh_thresh, '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    none', '|'
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%gyrotropic .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- GYROTROPIC   ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '| Compute Gyrotropic properties              :', pw90_calcs%gyrotropic, '|'
      write (stdout, '(1x,a46,10x,a20,1x,a1)') '| gyrotropic_task                            :', gyrotropic%task, '|'
      call parameters_gyro_write_task(gyrotropic%task, '-d0', 'calculate the D tensor', stdout)
      call parameters_gyro_write_task(gyrotropic%task, '-dw', 'calculate the tildeD tensor', stdout)
      call parameters_gyro_write_task(gyrotropic%task, '-c', 'calculate the C tensor', stdout)
      call parameters_gyro_write_task(gyrotropic%task, '-k', 'calculate the K tensor', stdout)
      call parameters_gyro_write_task(gyrotropic%task, '-noa', 'calculate the interbad natural optical activity', stdout)
      call parameters_gyro_write_task(gyrotropic%task, '-dos', 'calculate the density of states', stdout)

      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for tildeD,NOA            :', &
        write_data%gyrotropic_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency                           :', &
        write_data%gyrotropic_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for frequency                   :', &
        write_data%gyrotropic_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue                          :', &
        gyrotropic%eigval_max, '|'
      if (gyrotropic%smr_fixed_en_width == write_data%smear%fixed_width &
          .and. write_data%smear%type_index == gyrotropic%smr_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        write (stdout, '(1x,a78)') '|  Using local  smearing parameters                                          |'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
        gyrotropic%smr_fixed_en_width, '|'
      write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function                         :', &
        trim(param_get_smearing_type(gyrotropic%smr_index)), '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  degen_thresh                              :', gyrotropic%degen_thresh, '|'

      if (write_data%kmesh(1) == gyrotropic%kmesh(1) .and. &
          write_data%kmesh(2) == gyrotropic%kmesh(2) .and. &
          write_data%kmesh(3) == gyrotropic%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      elseif (gyrotropic%kmesh_spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          gyrotropic%kmesh(1), 'x', gyrotropic%kmesh(2), 'x', gyrotropic%kmesh(3), ' Spacing = ', gyrotropic%kmesh_spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , gyrotropic%kmesh(1), 'x', gyrotropic%kmesh(2), 'x', gyrotropic%kmesh(3), '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    not implemented', '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%boltzwann .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------------- BOLTZWANN ----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Boltzmann transport properties    :', &
        pw90_calcs%boltzwann, '|'
      if (boltz%dir_num_2d > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  2d structure: non-periodic dimension  :', &
          trim(write_data%boltz_2d_dir), '|'
      else
        write (stdout, '(1x,a78)') '|  3d Structure                              :                 T             |'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Relaxation Time (fs)                      :', boltz%relax_time, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Chemical Potential (eV)  :', boltz%mu_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Chemical Potential (eV)  :', boltz%mu_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Chemical Potential (eV)     :', boltz%mu_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Temperature (K)          :', boltz%temp_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Temperature (K)          :', boltz%temp_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Temperature (K)             :', boltz%temp_step, '|'

      if (write_data%kmesh(1) == boltz%kmesh(1) .and. write_data%kmesh(2) == boltz%kmesh(2) &
          .and. write_data%kmesh(3) == boltz%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (boltz%kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            boltz%kmesh(1), 'x', boltz%kmesh(2), 'x', boltz%kmesh(3), ' Spacing = ', boltz%kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , boltz%kmesh(1), 'x', boltz%kmesh(2), 'x', boltz%kmesh(3), '|'
        endif
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for TDF (eV)                    :', boltz%tdf_energy_step, '|'
      write (stdout, '(1x,a25,5x,a43,4x,a1)') '|  TDF Smearing Function ', trim(param_get_smearing_type(boltz%tdf_smr_index)), '|'
      if (boltz%tdf_smr_fixed_en_width > 0.0_dp) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
          '|  TDF fixed Smearing width (eV)             :', boltz%tdf_smr_fixed_en_width, '|'
      else
        write (stdout, '(1x,a78)') '|  TDF fixed Smearing width                  :         unsmeared             |'
      endif
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute DOS at same time                  :', boltz%calc_also_dos, '|'
      if (boltz%calc_also_dos .and. param_input%iprint > 2) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', boltz%dos_energy_min, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', boltz%dos_energy_max, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', boltz%dos_energy_step, '|'
        if (boltz%dos_smearing%use_adaptive .eqv. write_data%smear%use_adaptive .and. &
            boltz%dos_smearing%adaptive_prefactor == write_data%smear%adaptive_prefactor &
            .and. boltz%dos_smearing%adaptive_max_width == write_data%smear%adaptive_max_width &
            .and. boltz%dos_smearing%fixed_width == write_data%smear%fixed_width .and. &
            write_data%smear%type_index == boltz%dos_smearing%type_index) then
          write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
        else
          if (boltz%dos_smearing%use_adaptive) then
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Adaptive width smearing               :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Adaptive smearing factor              :', boltz%dos_smearing%adaptive_prefactor, '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Maximum allowed smearing width        :', boltz%dos_smearing%adaptive_max_width, '|'
          else
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Fixed width smearing                  :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  DOS Smearing width                         :', &
              boltz%dos_smearing%fixed_width, '|'
          endif
          write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
            trim(param_get_smearing_type(boltz%dos_smearing%type_index)), '|'
        endif
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%geninterp .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------Generic Band Interpolation--------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Properties at given k-points      :', pw90_calcs%geninterp, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Calculate band gradients                  :', geninterp%alsofirstder, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write data into a single file             :', geninterp%single_file, '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

101 format(20x, a3, 2x, 3F11.6)

  end subroutine param_postw90_write

  subroutine param_pw90_dealloc(excluded_bands, wann_data, kmesh_data, k_points, dis_window, &
                                fermi, atoms, eigval, spec_points, dos_data, berry, proj_input, &
                                stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    !data from parameters module
    type(exclude_bands_type), intent(inout) :: excluded_bands
    type(wannier_data_type), intent(inout) :: wann_data
    type(kmesh_input_type), intent(inout) :: kmesh_data
    type(proj_input_type), intent(inout) :: proj_input
    type(k_points_type), intent(inout) :: k_points
    type(dis_manifold_type), intent(inout) :: dis_window
    type(fermi_data_type), intent(inout) :: fermi
    type(atom_data_type), intent(inout) :: atoms
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    type(kpoint_path_type), intent(inout) :: spec_points
    type(pw90_dos_mod_type), intent(inout) :: dos_data
    type(pw90_berry_mod_type), intent(inout) :: berry
    character(len=50), intent(in)  :: seedname

    integer :: ierr

    call param_dealloc(excluded_bands, wann_data, proj_input, kmesh_data, k_points, &
                       dis_window, atoms, eigval, spec_points, stdout, seedname)
    if (allocated(dos_data%project)) then
      deallocate (dos_data%project, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating dos_project in param_pw90_dealloc', stdout, seedname)
    endif
    if (allocated(fermi%energy_list)) then
      deallocate (fermi%energy_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating fermi_energy_list in param_pw90_dealloc', stdout, seedname)
    endif
    if (allocated(berry%kubo_freq_list)) then
      deallocate (berry%kubo_freq_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kubo_freq_list in param_pw90_dealloc', stdout, seedname)
    endif
  end subroutine param_pw90_dealloc

  ! extra postw90 memory
  subroutine param_pw90_mem_estimate(mem_param, mem_bw, dis_window, do_boltzwann, boltz, &
                                     spin_decomp, num_wann, stdout)

    ! JJ, should only be called from root node

    implicit none

    type(dis_manifold_type), intent(in) :: dis_window
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    !real(kind=dp), parameter :: size_log = 1.0_dp
    !real(kind=dp), parameter :: size_int = 4.0_dp
    real(kind=dp), parameter :: size_real = 8.0_dp
    real(kind=dp), parameter :: size_cmplx = 16.0_dp
    real(kind=dp), intent(in) :: mem_param
    real(kind=dp), intent(inout) :: mem_bw
    logical, intent(in) :: do_boltzwann, spin_decomp
    type(pw90_boltzwann_type) :: boltz
    integer :: NumPoints1, NumPoints2, NumPoints3, ndim
    real(kind=dp) :: TDF_exceeding_energy

    if (do_boltzwann) then
      if (spin_decomp) then
        ndim = 3
      else
        ndim = 1
      end if

      ! I set a big value to have a rough estimate
      TDF_exceeding_energy = 2._dp
      NumPoints1 = int(floor((boltz%temp_max - boltz%temp_min)/boltz%temp_step)) + 1 ! temperature array
      NumPoints2 = int(floor((boltz%mu_max - boltz%mu_min)/boltz%mu_step)) + 1  ! mu array
      NumPoints3 = int(floor((dis_window%win_max - dis_window%win_min &
                              + 2._dp*TDF_exceeding_energy)/boltz%tdf_energy_step)) + 1 ! tdfenergyarray
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

      NumPoints1 = int(floor((boltz%dos_energy_max - boltz%dos_energy_min)/boltz%dos_energy_step)) + 1!dosnumpoints
      mem_bw = mem_bw + NumPoints1*size_real                         !DOS_EnergyArray
      mem_bw = mem_bw + 6*ndim*NumPoints3*size_real                  !TDF_k
      mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_k
      mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_all
    end if

    if (do_boltzwann) &
      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'BoltzWann:', (mem_param + mem_bw)/(1024**2), ' Mb'

  end subroutine param_pw90_mem_estimate

  subroutine parameters_gyro_write_task(task, key, comment, stdout)
!   use w90_io, only: stdout

    integer, intent(in) :: stdout
    character(len=*), intent(in) :: task, key, comment
    character(len=42) :: comment1

    comment1 = comment
    if ((index(task, key) > 0) .or. (index(task, 'all') > 0)) then
      write (stdout, '(1x,a2,a42,a2,10x,a8,13x,a1)') '| ', comment1, ' :', '       T', '|'
    else
      write (stdout, '(1x,a2,a42,a2,10x,a8,13x,a1)') '| ', comment1, ' :', '       F', '|'
    endif
  end subroutine parameters_gyro_write_task

end module pw90_param_methods
