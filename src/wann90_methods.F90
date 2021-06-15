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

module wannier_methods

  use w90_constants, only: dp
  use w90_types
  use w90_param_types
  use w90_param_methods
  use wannier_param_types

  implicit none

  private

  type w90_extra_io_type
    character(len=20) :: one_dim_axis
    !! Constrained centres
    real(kind=dp), allocatable :: ccentres_frac(:, :)
  end type w90_extra_io_type

  public :: w90_extra_io_type
  public :: param_read
  public :: param_write
  public :: param_w90_dealloc
  public :: param_write_chkpt
  public :: param_memory_estimate
  public :: param_dist

contains

  !==================================================================!
  subroutine param_read(driver, w90_calcs, pp_calc, param_input, param_plot, param_wannierise, &
                        lsitesymmetry, symmetrize_eps, wann_data, param_hamil, kmesh_data, &
                        kmesh_info, k_points, num_kpts, dis_data, fermi_surface_data, fermi, &
                        tran, atoms, num_bands, num_wann, eigval, mp_grid, num_proj, select_proj, &
                        real_lattice, recip_lattice, spec_points, eig_found, library, &
                        library_param_read_first_pass, bohr, stdout, seedname, write_data, &
                        proj, lhasproj)
    !==================================================================!
    !                                                                  !
    !! Read parameters and calculate derived values
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !                                                                  !
    !===================================================================
    use w90_constants, only: w90_physical_constants
    implicit none

    !data from parameters module
    type(param_driver_type), intent(inout) :: driver
    type(w90_calculation_type), intent(inout) :: w90_calcs
    type(postproc_type), intent(inout) :: pp_calc
    type(parameter_input_type), intent(inout) :: param_input
    type(param_plot_type), intent(inout) :: param_plot
    type(param_wannierise_type), intent(inout) :: param_wannierise
    ! RS: symmetry-adapted Wannier functions
    logical, intent(inout) :: lsitesymmetry
    real(kind=dp), intent(inout) :: symmetrize_eps
    type(wannier_data_type), intent(inout) :: wann_data
    type(param_hamiltonian_type), intent(inout) :: param_hamil
    type(param_kmesh_type), intent(inout) :: kmesh_data
    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(k_point_type), intent(inout) :: k_points
    integer, intent(inout) :: num_kpts
    type(disentangle_type), intent(inout) :: dis_data
    type(fermi_surface_type), intent(inout) :: fermi_surface_data
    type(fermi_data_type), intent(inout) :: fermi
    type(transport_type), intent(inout) :: tran
    type(atom_data_type), intent(inout) :: atoms
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_wann
    integer, intent(in) :: stdout
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    integer, intent(inout) :: mp_grid(3)
    integer, intent(inout) :: num_proj
    type(select_projection_type), intent(inout) :: select_proj
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: recip_lattice(3, 3)
    type(special_kpoints_type), intent(inout) :: spec_points
    logical, intent(inout) :: eig_found
    logical, intent(in) :: library
    logical, intent(in) :: library_param_read_first_pass
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname
    type(w90_extra_io_type), intent(inout) :: write_data
    ! was in driver, only used by wannier_lib
    type(projection_type), intent(inout) :: proj
    !Projections
    logical, intent(out) :: lhasproj

    !local variables
    character(len=20) :: energy_unit
    !! Units for energy
    logical                   :: found_fermi_energy
    real(kind=dp)             :: kmesh_spacing
    integer                   :: kmesh(3)
    logical                   :: global_kmesh_set
    logical :: has_kpath

    w90_calcs%disentanglement = .false.
    call param_in_file(stdout, seedname)
    call param_read_sym(lsitesymmetry, symmetrize_eps, stdout, seedname)
    call param_read_verbosity(param_input, stdout, seedname)
    call param_read_w90_calcs(w90_calcs, stdout, seedname)
    call param_read_transport(w90_calcs%transport, tran, driver%restart, stdout, seedname)
    call param_read_dist_cutoff(param_input, stdout, seedname)
    if (.not. (w90_calcs%transport .and. tran%read_ht)) then
      call param_read_units(param_input, energy_unit, bohr, stdout, seedname)
      call param_read_num_wann(num_wann, stdout, seedname)
      call param_read_exclude_bands(param_input, stdout, seedname)
      call param_read_num_bands(.false., library, param_input, num_bands, num_wann, &
                                library_param_read_first_pass, stdout, seedname)
      w90_calcs%disentanglement = (num_bands > num_wann)
      call param_read_lattice(library, real_lattice, recip_lattice, bohr, stdout, seedname)
      call param_read_wannierise(param_wannierise, num_wann, write_data%ccentres_frac, &
                                 stdout, seedname)
      call param_read_devel(param_input%devel_flag, stdout, seedname)
      call param_read_mp_grid(.false., library, mp_grid, num_kpts, stdout, seedname)
      call param_read_gamma_only(param_input%gamma_only, num_kpts, library, stdout, seedname)
      call param_read_post_proc(w90_calcs%cp_pp, pp_calc%only_A, driver%postproc_setup, stdout, seedname)
      call param_read_restart(driver, stdout, seedname)
      call param_read_system(library, param_input, stdout, seedname)
      call param_read_kpath(library, spec_points, has_kpath, stdout, seedname)
      call param_read_plot(w90_calcs, param_plot, param_input%bands_plot_mode, num_wann, has_kpath, stdout, seedname)
      call param_read_fermi_surface(fermi_surface_data, w90_calcs%fermi_surface_plot, stdout, seedname)
      call param_read_fermi_energy(found_fermi_energy, fermi, stdout, seedname)
      call param_read_outfiles(w90_calcs, param_input, param_wannierise, param_plot, num_kpts, stdout, seedname)
    endif
    ! BGS tran/plot related stuff...
    call param_read_one_dim(w90_calcs, param_plot, param_input, write_data%one_dim_axis, &
                            tran%read_ht, stdout, seedname)
    call param_read_ws_data(param_input, stdout, seedname) !ws_search etc
    if (.not. (w90_calcs%transport .and. tran%read_ht)) then
      call param_read_eigvals(.false., .false., .false., &
                              w90_calcs%bands_plot .or. w90_calcs%fermi_surface_plot .or. &
                              w90_calcs%write_hr, w90_calcs%disentanglement, eig_found, &
                              eigval, library, driver%postproc_setup, num_bands, num_kpts, stdout, seedname)
      dis_data%win_min = -1.0_dp
      dis_data%win_max = 0.0_dp
      if (eig_found) dis_data%win_min = minval(eigval)
      if (eig_found) dis_data%win_max = maxval(eigval)
      call param_read_disentangle_all(eig_found, dis_data, stdout, seedname)
      call param_read_disentangle_w90(dis_data, num_bands, num_wann, bohr, stdout, seedname)
      call param_read_hamil(param_hamil, stdout, seedname)
      call param_read_bloch_phase(w90_calcs%use_bloch_phases, w90_calcs%disentanglement, stdout, seedname)
      call param_read_kmesh_data(kmesh_data, stdout, seedname)
      call param_read_kpoints(.false., library, k_points, num_kpts, recip_lattice, bohr, stdout, seedname)
      call param_read_explicit_kpts(library, driver, kmesh_info, num_kpts, bohr, stdout, seedname)
      call param_read_global_kmesh(global_kmesh_set, kmesh_spacing, kmesh, recip_lattice, stdout, seedname)
      call param_read_atoms(library, atoms, real_lattice, recip_lattice, bohr, stdout, seedname)
      call param_read_projections(proj, w90_calcs%use_bloch_phases, lhasproj, &
                                  param_wannierise%guiding_centres, param_wannierise%proj_site, &
                                  kmesh_data, select_proj, num_proj, param_input, atoms, &
                                  recip_lattice, num_wann, library, bohr, stdout, seedname)
      ! projections needs to be allocated before reading constrained centres
      if (param_wannierise%slwf_constrain) then
        call param_read_constrained_centres(write_data%ccentres_frac, param_wannierise, &
                                            real_lattice, num_wann, library, stdout, seedname)
      endif
    endif
    call param_clean_infile(stdout, seedname)
    if (.not. (w90_calcs%transport .and. tran%read_ht)) then
      ! For aesthetic purposes, convert some things to uppercase
      call param_uppercase(param_input, atoms, spec_points)

      param_wannierise%omega_total = -999.0_dp
      param_wannierise%omega_tilde = -999.0_dp
      ! Initialise
      param_wannierise%omega_total = -999.0_dp
      param_wannierise%omega_tilde = -999.0_dp
      param_input%omega_invariant = -999.0_dp
      param_input%have_disentangled = .false.
      call param_read_final_alloc(w90_calcs%disentanglement, dis_data, &
                                  wann_data, num_wann, num_bands, num_kpts, stdout, seedname)
    endif
  end subroutine param_read

  !==================================================================!
  subroutine param_read_sym(lsitesymmetry, symmetrize_eps, stdout, seedname)
    !%%%%%%%%%%%%%%%%
    ! Site symmetry
    !%%%%%%%%%%%%%%%%
    implicit none
    logical, intent(inout) :: lsitesymmetry
    real(kind=dp), intent(inout) :: symmetrize_eps
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    logical :: found

    ! default value is lsitesymmetry=.false.
    call param_get_keyword(stdout, seedname, 'site_symmetry', found, l_value=lsitesymmetry)!YN:

    ! default value is symmetrize_eps=0.001
    call param_get_keyword(stdout, seedname, 'symmetrize_eps', found, r_value=symmetrize_eps)!YN:
  end subroutine param_read_sym

  subroutine param_read_w90_calcs(w90_calcs, stdout, seedname)
    !%%%%%%%%%%%%%%%%
    ! Transport
    !%%%%%%%%%%%%%%%%
!   use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(w90_calculation_type), intent(out) :: w90_calcs
    character(len=50), intent(in)  :: seedname

    logical :: found

    w90_calcs%transport = .false.
    call param_get_keyword(stdout, seedname, 'transport', found, l_value=w90_calcs%transport)

    w90_calcs%wannier_plot = .false.
    call param_get_keyword(stdout, seedname, 'wannier_plot', found, l_value=w90_calcs%wannier_plot)

    w90_calcs%bands_plot = .false.
    call param_get_keyword(stdout, seedname, 'bands_plot', found, l_value=w90_calcs%bands_plot)

    w90_calcs%fermi_surface_plot = .false.
    call param_get_keyword(stdout, seedname, 'fermi_surface_plot', found, l_value=w90_calcs%fermi_surface_plot)

  end subroutine param_read_w90_calcs

  subroutine param_read_transport(transport, tran, restart, stdout, seedname)
    !%%%%%%%%%%%%%%%%
    ! Transport
    !%%%%%%%%%%%%%%%%
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: transport
    type(transport_type), intent(out) :: tran
    character(len=*), intent(inout) :: restart
    character(len=50), intent(in)  :: seedname

    logical :: found

    tran%read_ht = .false.
    call param_get_keyword(stdout, seedname, 'tran_read_ht', found, l_value=tran%read_ht)

    tran%easy_fix = .false.
    call param_get_keyword(stdout, seedname, 'tran_easy_fix', found, l_value=tran%easy_fix)

    if (transport .and. tran%read_ht) restart = ' '

    tran%mode = 'bulk'
    call param_get_keyword(stdout, seedname, 'transport_mode', found, c_value=tran%mode)

!    if ( .not.tran_read_ht  .and. (index(transport_mode,'lcr').ne.0) ) &
!       call io_error('Error: transport_mode.eq.lcr not compatible with tran_read_ht.eq.false')

    tran%win_min = -3.0_dp
    call param_get_keyword(stdout, seedname, 'tran_win_min', found, r_value=tran%win_min)

    tran%win_max = 3.0_dp
    call param_get_keyword(stdout, seedname, 'tran_win_max', found, r_value=tran%win_max)

    tran%energy_step = 0.01_dp
    call param_get_keyword(stdout, seedname, 'tran_energy_step', found, r_value=tran%energy_step)

    tran%num_bb = 0
    call param_get_keyword(stdout, seedname, 'tran_num_bb', found, i_value=tran%num_bb)

    tran%num_ll = 0
    call param_get_keyword(stdout, seedname, 'tran_num_ll', found, i_value=tran%num_ll)

    tran%num_rr = 0
    call param_get_keyword(stdout, seedname, 'tran_num_rr', found, i_value=tran%num_rr)

    tran%num_cc = 0
    call param_get_keyword(stdout, seedname, 'tran_num_cc', found, i_value=tran%num_cc)

    tran%num_lc = 0
    call param_get_keyword(stdout, seedname, 'tran_num_lc', found, i_value=tran%num_lc)

    tran%num_cr = 0
    call param_get_keyword(stdout, seedname, 'tran_num_cr', found, i_value=tran%num_cr)

    tran%num_bandc = 0
    call param_get_keyword(stdout, seedname, 'tran_num_bandc', found, i_value=tran%num_bandc)

    tran%write_ht = .false.
    call param_get_keyword(stdout, seedname, 'tran_write_ht', found, l_value=tran%write_ht)

    tran%use_same_lead = .true.
    call param_get_keyword(stdout, seedname, 'tran_use_same_lead', found, l_value=tran%use_same_lead)

    tran%num_cell_ll = 0
    call param_get_keyword(stdout, seedname, 'tran_num_cell_ll', found, i_value=tran%num_cell_ll)

    tran%num_cell_rr = 0
    call param_get_keyword(stdout, seedname, 'tran_num_cell_rr', found, i_value=tran%num_cell_rr)

    tran%group_threshold = 0.15_dp
    call param_get_keyword(stdout, seedname, 'tran_group_threshold', found, r_value=tran%group_threshold)

    ! checks
    if (transport) then
      if ((index(tran%mode, 'bulk') .eq. 0) .and. (index(tran%mode, 'lcr') .eq. 0)) &
        call io_error('Error: transport_mode not recognised', stdout, seedname)
      if (tran%num_bb < 0) call io_error('Error: tran_num_bb < 0', stdout, seedname)
      if (tran%num_ll < 0) call io_error('Error: tran_num_ll < 0', stdout, seedname)
      if (tran%num_rr < 0) call io_error('Error: tran_num_rr < 0', stdout, seedname)
      if (tran%num_cc < 0) call io_error('Error: tran_num_cc < 0', stdout, seedname)
      if (tran%num_lc < 0) call io_error('Error: tran_num_lc < 0', stdout, seedname)
      if (tran%num_cr < 0) call io_error('Error: tran_num_cr < 0', stdout, seedname)
      if (tran%num_bandc < 0) call io_error('Error: tran_num_bandc < 0', stdout, seedname)
      if (tran%num_cell_ll < 0) call io_error('Error: tran_num_cell_ll < 0', stdout, seedname)
      if (tran%num_cell_rr < 0) call io_error('Error: tran_num_cell_rr < 0', stdout, seedname)
      if (tran%group_threshold < 0.0_dp) call io_error('Error: tran_group_threshold < 0', stdout, seedname)
    endif

  end subroutine param_read_transport

  subroutine param_read_dist_cutoff(param_input, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(parameter_input_type), intent(inout) :: param_input
    character(len=50), intent(in)  :: seedname

    logical :: found

    param_input%dist_cutoff_mode = 'three_dim'
    call param_get_keyword(stdout, seedname, 'dist_cutoff_mode', found, c_value=param_input%dist_cutoff_mode)
    if ((index(param_input%dist_cutoff_mode, 'three_dim') .eq. 0) &
        .and. (index(param_input%dist_cutoff_mode, 'two_dim') .eq. 0) &
        .and. (index(param_input%dist_cutoff_mode, 'one_dim') .eq. 0)) &
      call io_error('Error: dist_cutoff_mode not recognised', stdout, seedname)

    param_input%dist_cutoff = 1000.0_dp
    call param_get_keyword(stdout, seedname, 'dist_cutoff', found, r_value=param_input%dist_cutoff)

    param_input%dist_cutoff_hc = param_input%dist_cutoff
    call param_get_keyword(stdout, seedname, 'dist_cutoff_hc', found, r_value=param_input%dist_cutoff_hc)

    param_input%hr_cutoff = 0.0_dp
    call param_get_keyword(stdout, seedname, 'hr_cutoff', found, r_value=param_input%hr_cutoff)

  end subroutine param_read_dist_cutoff

  subroutine param_read_wannierise(param_wannierise, num_wann, ccentres_frac, stdout, seedname)
    !%%%%%%%%%%%
    ! Wannierise
    !%%%%%%%%%%%
    use w90_io, only: io_error
    implicit none
    type(param_wannierise_type), intent(out) :: param_wannierise
    integer, intent(in) :: num_wann
    real(kind=dp), allocatable, intent(inout) :: ccentres_frac(:, :)
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: ierr
    logical :: found

    param_wannierise%num_dump_cycles = 100     ! frequency to write backups at
    call param_get_keyword(stdout, seedname, 'num_dump_cycles', found, i_value=param_wannierise%num_dump_cycles)
    if (param_wannierise%num_dump_cycles < 0) call io_error('Error: num_dump_cycles must be positive', stdout, seedname)

    param_wannierise%num_print_cycles = 1          ! frequency to write at
    call param_get_keyword(stdout, seedname, 'num_print_cycles', found, i_value=param_wannierise%num_print_cycles)
    if (param_wannierise%num_print_cycles < 0) call io_error('Error: num_print_cycles must be positive', stdout, seedname)

    param_wannierise%num_iter = 100
    call param_get_keyword(stdout, seedname, 'num_iter', found, i_value=param_wannierise%num_iter)
    if (param_wannierise%num_iter < 0) call io_error('Error: num_iter must be positive', stdout, seedname)

    param_wannierise%num_cg_steps = 5
    call param_get_keyword(stdout, seedname, 'num_cg_steps', found, i_value=param_wannierise%num_cg_steps)
    if (param_wannierise%num_cg_steps < 0) call io_error('Error: num_cg_steps must be positive', stdout, seedname)

    param_wannierise%conv_tol = 1.0e-10_dp
    call param_get_keyword(stdout, seedname, 'conv_tol', found, r_value=param_wannierise%conv_tol)
    if (param_wannierise%conv_tol < 0.0_dp) call io_error('Error: conv_tol must be positive', stdout, seedname)

    param_wannierise%conv_noise_amp = -1.0_dp
    call param_get_keyword(stdout, seedname, 'conv_noise_amp', found, r_value=param_wannierise%conv_noise_amp)

    ! JJ why is this -1 by default?  it implies that no checking is made for convergence
    param_wannierise%conv_window = -1
    if (param_wannierise%conv_noise_amp > 0.0_dp) param_wannierise%conv_window = 5
    call param_get_keyword(stdout, seedname, 'conv_window', found, i_value=param_wannierise%conv_window)

    param_wannierise%conv_noise_num = 3
    call param_get_keyword(stdout, seedname, 'conv_noise_num', found, i_value=param_wannierise%conv_noise_num)
    if (param_wannierise%conv_noise_num < 0) call io_error('Error: conv_noise_num must be positive', stdout, seedname)

    param_wannierise%guiding_centres = .false.
    call param_get_keyword(stdout, seedname, 'guiding_centres', found, l_value=param_wannierise%guiding_centres)

    param_wannierise%num_guide_cycles = 1
    call param_get_keyword(stdout, seedname, 'num_guide_cycles', found, i_value=param_wannierise%num_guide_cycles)
    if (param_wannierise%num_guide_cycles < 0) call io_error('Error: num_guide_cycles must be >= 0', stdout, seedname)

    param_wannierise%num_no_guide_iter = 0
    call param_get_keyword(stdout, seedname, 'num_no_guide_iter', found, i_value=param_wannierise%num_no_guide_iter)
    if (param_wannierise%num_no_guide_iter < 0) call io_error('Error: num_no_guide_iter must be >= 0', stdout, seedname)

    param_wannierise%fixed_step = -999.0_dp; 
    param_wannierise%lfixstep = .false.
    call param_get_keyword(stdout, seedname, 'fixed_step', found, r_value=param_wannierise%fixed_step)
    if (found .and. (param_wannierise%fixed_step < 0.0_dp)) call io_error('Error: fixed_step must be > 0', stdout, seedname)
    if (param_wannierise%fixed_step > 0.0_dp) param_wannierise%lfixstep = .true.

    param_wannierise%trial_step = 2.0_dp
    call param_get_keyword(stdout, seedname, 'trial_step', found, r_value=param_wannierise%trial_step)
    if (found .and. param_wannierise%lfixstep) then
      call io_error('Error: cannot specify both fixed_step and trial_step', stdout, seedname)
    endif

    param_wannierise%precond = .false.
    call param_get_keyword(stdout, seedname, 'precond', found, l_value=param_wannierise%precond)

    param_wannierise%slwf_num = num_wann
    param_wannierise%selective_loc = .false.
    call param_get_keyword(stdout, seedname, 'slwf_num', found, i_value=param_wannierise%slwf_num)
    if (found) then
      if (param_wannierise%slwf_num .gt. num_wann .or. param_wannierise%slwf_num .lt. 1) then
        call io_error('Error: slwf_num must be an integer between 1 and num_wann', stdout, seedname)
      end if
      if (param_wannierise%slwf_num .lt. num_wann) param_wannierise%selective_loc = .true.
    end if

    param_wannierise%slwf_constrain = .false.
    call param_get_keyword(stdout, seedname, 'slwf_constrain', found, l_value=param_wannierise%slwf_constrain)
    if (found .and. param_wannierise%slwf_constrain) then
      if (param_wannierise%selective_loc) then
        allocate (ccentres_frac(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_frac in param_get_centre_constraints', stdout, seedname)
        allocate (param_wannierise%ccentres_cart(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_cart in param_get_centre_constraints', stdout, seedname)
      else
        write (stdout, *) ' No selective localisation requested. Ignoring constraints on centres'
        param_wannierise%slwf_constrain = .false.
      end if
    end if

    param_wannierise%slwf_lambda = 1.0_dp
    call param_get_keyword(stdout, seedname, 'slwf_lambda', found, r_value=param_wannierise%slwf_lambda)
    if (found) then
      if (param_wannierise%slwf_lambda < 0.0_dp) call io_error('Error: slwf_lambda  must be positive.', stdout, seedname)
    endif

    param_wannierise%translate_home_cell = .false.
    call param_get_keyword(stdout, seedname, 'translate_home_cell', found, l_value=param_wannierise%translate_home_cell)
  end subroutine param_read_wannierise

  subroutine param_read_disentangle_w90(dis_data, num_bands, num_wann, bohr, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    !logical, intent(in) :: eig_found
    !real(kind=dp), intent(in) :: eigval(:, :)
    integer, intent(in) :: stdout
    type(disentangle_type), intent(inout) :: dis_data
    integer, intent(in) :: num_bands, num_wann
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname

    integer :: nkp, ierr
    logical :: found

    dis_data%num_iter = 200
    call param_get_keyword(stdout, seedname, 'dis_num_iter', found, i_value=dis_data%num_iter)
    if (dis_data%num_iter < 0) call io_error('Error: dis_num_iter must be positive', stdout, seedname)

    dis_data%mix_ratio = 0.5_dp
    call param_get_keyword(stdout, seedname, 'dis_mix_ratio', found, r_value=dis_data%mix_ratio)
    if (dis_data%mix_ratio <= 0.0_dp .or. dis_data%mix_ratio > 1.0_dp) &
      call io_error('Error: dis_mix_ratio must be greater than 0.0 but not greater than 1.0', stdout, seedname)

    dis_data%conv_tol = 1.0e-10_dp
    call param_get_keyword(stdout, seedname, 'dis_conv_tol', found, r_value=dis_data%conv_tol)
    if (dis_data%conv_tol < 0.0_dp) call io_error('Error: dis_conv_tol must be positive', stdout, seedname)

    dis_data%conv_window = 3
    call param_get_keyword(stdout, seedname, 'dis_conv_window', found, i_value=dis_data%conv_window)
    if (dis_data%conv_window < 0) call io_error('Error: dis_conv_window must be positive', stdout, seedname)

    ! GS-start
    dis_data%spheres_first_wann = 1
    call param_get_keyword(stdout, seedname, 'dis_spheres_first_wann', found, i_value=dis_data%spheres_first_wann)
    if (dis_data%spheres_first_wann < 1) call io_error('Error: dis_spheres_first_wann must be greater than 0', stdout, seedname)
    if (dis_data%spheres_first_wann > num_bands - num_wann + 1) &
      call io_error('Error: dis_spheres_first_wann is larger than num_bands-num_wann+1', stdout, seedname)
    dis_data%spheres_num = 0
    call param_get_keyword(stdout, seedname, 'dis_spheres_num', found, i_value=dis_data%spheres_num)
    if (dis_data%spheres_num < 0) call io_error('Error: dis_spheres_num cannot be negative', stdout, seedname)
    if (dis_data%spheres_num > 0) then
      allocate (dis_data%spheres(4, dis_data%spheres_num), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dis_spheres in param_read', stdout, seedname)
      call param_get_keyword_block(stdout, seedname, 'dis_spheres', found, dis_data%spheres_num, 4, bohr, &
                                   r_value=dis_data%spheres)
      if (.not. found) call io_error('Error: Did not find dis_spheres in the input file', stdout, seedname)
      do nkp = 1, dis_data%spheres_num
        if (dis_data%spheres(4, nkp) < 1.0e-15_dp) &
          call io_error('Error: radius for dis_spheres must be > 0', stdout, seedname)
      enddo
    endif
    ! GS-end
  end subroutine param_read_disentangle_w90

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

  subroutine param_read_post_proc(cp_pp, pp_only_A, postproc_setup, stdout, seedname)
!   use w90_io, only: post_proc_flag, io_error
    use w90_io, only: post_proc_flag
    implicit none
    integer, intent(in) :: stdout
    logical, intent(out) :: cp_pp, pp_only_A, postproc_setup
    character(len=50), intent(in)  :: seedname

    logical :: found

    postproc_setup = .false.            ! set to true to write .nnkp file and exit
    call param_get_keyword(stdout, seedname, 'postproc_setup', found, l_value=postproc_setup)
    ! We allow this keyword to be overriden by a command line arg -pp
    if (post_proc_flag) postproc_setup = .true.

    cp_pp = .false.         ! set to true if doing CP post-processing
    call param_get_keyword(stdout, seedname, 'cp_pp', found, l_value=cp_pp)

    pp_only_A = .false.
    call param_get_keyword(stdout, seedname, 'calc_only_A', found, l_value=pp_only_A)
  end subroutine param_read_post_proc

  subroutine param_read_restart(driver, stdout, seedname)
!   use w90_io, only: seedname, io_error
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(param_driver_type), intent(inout) :: driver
    character(len=50), intent(in)  :: seedname

    logical :: found, chk_found

    driver%restart = ' '
    call param_get_keyword(stdout, seedname, 'restart', found, c_value=driver%restart)
    if (found) then
      if ((driver%restart .ne. 'default') .and. (driver%restart .ne. 'wannierise') &
          .and. (driver%restart .ne. 'plot') .and. (driver%restart .ne. 'transport')) then
        call io_error('Error in input file: value of restart not recognised', stdout, seedname)
      else
        inquire (file=trim(seedname)//'.chk', exist=chk_found)
        if (.not. chk_found) &
          call io_error('Error: restart requested but '//trim(seedname)//'.chk file not found', stdout, seedname)
      endif
    endif
    !post processing takes priority (user is not warned of this)
    if (driver%postproc_setup) driver%restart = ' '
  end subroutine param_read_restart

  subroutine param_read_outfiles(w90_calcs, param_input, param_wannierise, param_plot, num_kpts, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    type(w90_calculation_type), intent(inout) :: w90_calcs
    type(parameter_input_type), intent(inout) :: param_input ! write_xyz
    type(param_wannierise_type), intent(inout) :: param_wannierise
    type(param_plot_type), intent(inout) :: param_plot
    integer, intent(in) :: stdout
    integer, intent(in) :: num_kpts
    character(len=50), intent(in)  :: seedname

    logical :: found, hr_plot

    param_input%write_xyz = .false.
    call param_get_keyword(stdout, seedname, 'write_xyz', found, l_value=param_input%write_xyz)

    param_wannierise%write_r2mn = .false.
    call param_get_keyword(stdout, seedname, 'write_r2mn', found, l_value=param_wannierise%write_r2mn)

    param_wannierise%write_proj = .false.
    call param_get_keyword(stdout, seedname, 'write_proj', found, l_value=param_wannierise%write_proj)

    param_wannierise%write_hr_diag = .false.
    call param_get_keyword(stdout, seedname, 'write_hr_diag', found, l_value=param_wannierise%write_hr_diag)

    hr_plot = .false.
    call param_get_keyword(stdout, seedname, 'hr_plot', found, l_value=hr_plot)
    if (found) call io_error('Input parameter hr_plot is no longer used. Please use write_hr instead.', stdout, seedname)
    w90_calcs%write_hr = .false.
    call param_get_keyword(stdout, seedname, 'write_hr', found, l_value=w90_calcs%write_hr)

    param_plot%write_rmn = .false.
    call param_get_keyword(stdout, seedname, 'write_rmn', found, l_value=param_plot%write_rmn)

    param_plot%write_tb = .false.
    call param_get_keyword(stdout, seedname, 'write_tb', found, l_value=param_plot%write_tb)

    !%%%%%%%%%%%%%%%%
    !  Other Stuff
    !%%%%%%%%%%%%%%%%

    ! aam: vdW
    param_wannierise%write_vdw_data = .false.
    call param_get_keyword(stdout, seedname, 'write_vdw_data', found, l_value=param_wannierise%write_vdw_data)
    if (param_wannierise%write_vdw_data) then
      if ((.not. param_input%gamma_only) .or. (num_kpts .ne. 1)) &
        call io_error('Error: write_vdw_data may only be used with a single k-point at Gamma', stdout, seedname)
    endif
    if (param_wannierise%write_vdw_data .and. w90_calcs%disentanglement .and. &
        param_input%num_valence_bands <= 0) &
      call io_error('If writing vdw data and disentangling then num_valence_bands must be defined', stdout, seedname)

  end subroutine param_read_outfiles

  subroutine param_read_plot(w90_calcs, param_plot, bands_plot_mode, num_wann, has_kpath, stdout, seedname)
    !%%%%%%%%%
    ! Plotting
    !%%%%%%%%%
    use w90_io, only: io_error
    implicit none
    type(w90_calculation_type), intent(in) :: w90_calcs
    type(param_plot_type), intent(out) :: param_plot
    character(len=*), intent(out) :: bands_plot_mode
    integer, intent(in) :: stdout
    integer, intent(in) :: num_wann
    logical, intent(in) :: has_kpath
    character(len=50), intent(in)  :: seedname

    integer :: i, loop, ierr
    logical :: found
    character(len=6) :: spin_str

    param_plot%wvfn_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword(stdout, seedname, 'wvfn_formatted', found, l_value=param_plot%wvfn_formatted)

    param_plot%spin = 1
    call param_get_keyword(stdout, seedname, 'spin', found, c_value=spin_str)
    if (found) then
      if (index(spin_str, 'up') > 0) then
        param_plot%spin = 1
      elseif (index(spin_str, 'down') > 0) then
        param_plot%spin = 2
      else
        call io_error('Error: unrecognised value of spin found: '//trim(spin_str), stdout, seedname)
      end if
    end if

    param_plot%wannier_plot_supercell = 2

    call param_get_vector_length(stdout, seedname, 'wannier_plot_supercell', found, length=i)
    if (found) then
      if (i .eq. 1) then
        call param_get_keyword_vector(stdout, seedname, 'wannier_plot_supercell', found, 1, &
                                      i_value=param_plot%wannier_plot_supercell)
        param_plot%wannier_plot_supercell(2) = param_plot%wannier_plot_supercell(1)
        param_plot%wannier_plot_supercell(3) = param_plot%wannier_plot_supercell(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector(stdout, seedname, 'wannier_plot_supercell', found, 3, &
                                      i_value=param_plot%wannier_plot_supercell)
      else
        call io_error('Error: wannier_plot_supercell must be provided as either one integer or a vector of three integers', &
                      stdout, seedname)
      end if
      if (any(param_plot%wannier_plot_supercell <= 0)) &
        call io_error('Error: wannier_plot_supercell elements must be greater than zero', stdout, seedname)
    end if

    param_plot%wannier_plot_format = 'xcrysden'
    call param_get_keyword(stdout, seedname, 'wannier_plot_format', found, c_value=param_plot%wannier_plot_format)

    param_plot%wannier_plot_mode = 'crystal'
    call param_get_keyword(stdout, seedname, 'wannier_plot_mode', found, c_value=param_plot%wannier_plot_mode)

    param_plot%wannier_plot_spinor_mode = 'total'
    call param_get_keyword(stdout, seedname, 'wannier_plot_spinor_mode', found, c_value=param_plot%wannier_plot_spinor_mode)
    param_plot%wannier_plot_spinor_phase = .true.
    call param_get_keyword(stdout, seedname, 'wannier_plot_spinor_phase', found, l_value=param_plot%wannier_plot_spinor_phase)

    call param_get_range_vector(stdout, seedname, 'wannier_plot_list', found, param_plot%num_wannier_plot, lcount=.true.)
    if (found) then
      if (param_plot%num_wannier_plot < 1) call io_error('Error: problem reading wannier_plot_list', stdout, seedname)
      if (allocated(param_plot%wannier_plot_list)) deallocate (param_plot%wannier_plot_list)
      allocate (param_plot%wannier_plot_list(param_plot%num_wannier_plot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_plot_list in param_read', stdout, seedname)
      call param_get_range_vector(stdout, seedname, 'wannier_plot_list', found, &
                                  param_plot%num_wannier_plot, .false., &
                                  param_plot%wannier_plot_list)
      if (any(param_plot%wannier_plot_list < 1) .or. any(param_plot%wannier_plot_list > num_wann)) &
        call io_error('Error: wannier_plot_list asks for a non-valid wannier function to be plotted', stdout, seedname)
    else
      ! we plot all wannier functions
      param_plot%num_wannier_plot = num_wann
      if (allocated(param_plot%wannier_plot_list)) deallocate (param_plot%wannier_plot_list)
      allocate (param_plot%wannier_plot_list(param_plot%num_wannier_plot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_plot_list in param_read', stdout, seedname)
      do loop = 1, num_wann
        param_plot%wannier_plot_list(loop) = loop
      end do
    end if

    param_plot%wannier_plot_radius = 3.5_dp
    call param_get_keyword(stdout, seedname, 'wannier_plot_radius', found, r_value=param_plot%wannier_plot_radius)

    param_plot%wannier_plot_scale = 1.0_dp
    call param_get_keyword(stdout, seedname, 'wannier_plot_scale', found, r_value=param_plot%wannier_plot_scale)

    ! checks
    if (w90_calcs%wannier_plot) then
      if ((index(param_plot%wannier_plot_format, 'xcrys') .eq. 0) .and. (index(param_plot%wannier_plot_format, 'cub') .eq. 0)) &
        call io_error('Error: wannier_plot_format not recognised', stdout, seedname)
      if ((index(param_plot%wannier_plot_mode, 'crys') .eq. 0) .and. (index(param_plot%wannier_plot_mode, 'mol') .eq. 0)) &
        call io_error('Error: wannier_plot_mode not recognised', stdout, seedname)
      if ((index(param_plot%wannier_plot_spinor_mode, 'total') .eq. 0) &
          .and. (index(param_plot%wannier_plot_spinor_mode, 'up') .eq. 0) &
          .and. (index(param_plot%wannier_plot_spinor_mode, 'down') .eq. 0)) &
        call io_error('Error: wannier_plot_spinor_mode not recognised', stdout, seedname)
      if (param_plot%wannier_plot_radius < 0.0_dp) call io_error('Error: wannier_plot_radius must be positive', stdout, seedname)
      if (param_plot%wannier_plot_scale < 0.0_dp) call io_error('Error: wannier_plot_scale must be positive', stdout, seedname)
    endif

    param_plot%write_u_matrices = .false.
    call param_get_keyword(stdout, seedname, 'write_u_matrices', found, l_value=param_plot%write_u_matrices)

    param_plot%write_bvec = .false.
    call param_get_keyword(stdout, seedname, 'write_bvec', found, l_value=param_plot%write_bvec)

    param_plot%bands_num_points = 100
    call param_get_keyword(stdout, seedname, 'bands_num_points', found, i_value=param_plot%bands_num_points)

    param_plot%bands_plot_format = 'gnuplot'
    call param_get_keyword(stdout, seedname, 'bands_plot_format', found, c_value=param_plot%bands_plot_format)

    bands_plot_mode = 's-k'
    call param_get_keyword(stdout, seedname, 'bands_plot_mode', found, c_value=bands_plot_mode)

    param_plot%bands_plot_dim = 3
    call param_get_keyword(stdout, seedname, 'bands_plot_dim', found, i_value=param_plot%bands_plot_dim)

    param_plot%num_bands_project = 0
    call param_get_range_vector(stdout, seedname, 'bands_plot_project', found, param_plot%num_bands_project, lcount=.true.)
    if (found) then
      if (param_plot%num_bands_project < 1) call io_error('Error: problem reading bands_plot_project', stdout, seedname)
      if (allocated(param_plot%bands_plot_project)) deallocate (param_plot%bands_plot_project)
      allocate (param_plot%bands_plot_project(param_plot%num_bands_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating bands_plot_project in param_read', stdout, seedname)
      call param_get_range_vector(stdout, seedname, 'bands_plot_project', found, &
                                  param_plot%num_bands_project, .false., &
                                  param_plot%bands_plot_project)
      if (any(param_plot%bands_plot_project < 1) .or. any(param_plot%bands_plot_project > num_wann)) &
        call io_error('Error: bands_plot_project asks for a non-valid wannier function to be projected', stdout, seedname)
    endif

    if (.not. has_kpath .and. w90_calcs%bands_plot) &
      call io_error('A bandstructure plot has been requested but there is no kpoint_path block', stdout, seedname)

    ! checks
    if (w90_calcs%bands_plot) then
      if ((index(param_plot%bands_plot_format, 'gnu') .eq. 0) .and. &
          (index(param_plot%bands_plot_format, 'xmgr') .eq. 0)) &
        call io_error('Error: bands_plot_format not recognised', stdout, seedname)
      if ((index(bands_plot_mode, 's-k') .eq. 0) .and. (index(bands_plot_mode, 'cut') .eq. 0)) &
        call io_error('Error: bands_plot_mode not recognised', stdout, seedname)
      if (param_plot%bands_num_points < 0) call io_error('Error: bands_num_points must be positive', stdout, seedname)
    endif

  end subroutine param_read_plot

  subroutine param_read_fermi_surface(fermi_surface_data, fermi_surface_plot, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(fermi_surface_type), intent(out) :: fermi_surface_data
    logical, intent(in) :: fermi_surface_plot
    character(len=50), intent(in)  :: seedname

    logical :: found

    fermi_surface_data%num_points = 50
    call param_get_keyword(stdout, seedname, 'fermi_surface_num_points', found, i_value=fermi_surface_data%num_points)

    fermi_surface_data%plot_format = 'xcrysden'
    call param_get_keyword(stdout, seedname, 'fermi_surface_plot_format', &
                           found, c_value=fermi_surface_data%plot_format)

    if (fermi_surface_plot) then
      if ((index(fermi_surface_data%plot_format, 'xcrys') .eq. 0)) &
        call io_error('Error: fermi_surface_plot_format not recognised', stdout, seedname)
      if (fermi_surface_data%num_points < 0) &
        call io_error('Error: fermi_surface_num_points must be positive', stdout, seedname)
    endif
  end subroutine param_read_fermi_surface

  subroutine param_read_one_dim(w90_calcs, param_plot, param_input, one_dim_axis, tran_read_ht, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    type(w90_calculation_type), intent(in) :: w90_calcs
    type(param_plot_type), intent(in) :: param_plot
    type(parameter_input_type), intent(inout) :: param_input
    character(len=*), intent(out) :: one_dim_axis
    logical, intent(in) :: tran_read_ht
    character(len=50), intent(in)  :: seedname

    logical :: found

    one_dim_axis = 'none'
    call param_get_keyword(stdout, seedname, 'one_dim_axis', found, c_value=one_dim_axis)
    param_input%one_dim_dir = 0
    if (index(one_dim_axis, 'x') > 0) param_input%one_dim_dir = 1
    if (index(one_dim_axis, 'y') > 0) param_input%one_dim_dir = 2
    if (index(one_dim_axis, 'z') > 0) param_input%one_dim_dir = 3
    if (w90_calcs%transport .and. .not. tran_read_ht .and. &
        (param_input%one_dim_dir .eq. 0)) call io_error('Error: one_dim_axis not recognised', stdout, seedname)
    if (w90_calcs%bands_plot .and. (index(param_input%bands_plot_mode, 'cut') .ne. 0) .and. &
        ((param_plot%bands_plot_dim .ne. 3) .or. &
         (index(param_input%dist_cutoff_mode, 'three_dim') .eq. 0)) .and. &
        (param_input%one_dim_dir .eq. 0)) &
      call io_error('Error: one_dim_axis not recognised', stdout, seedname)

  end subroutine param_read_one_dim

  subroutine param_read_hamil(param_hamil, stdout, seedname)
    implicit none
    integer, intent(in) :: stdout
    type(param_hamiltonian_type), intent(out) :: param_hamil
    real(kind=dp) :: rv_temp(3)
    character(len=50), intent(in)  :: seedname

    logical :: found

    param_hamil%automatic_translation = .true.
    param_hamil%translation_centre_frac = 0.0_dp
    call param_get_keyword_vector(stdout, seedname, 'translation_centre_frac', found, 3, r_value=rv_temp)
    if (found) then
      param_hamil%translation_centre_frac = rv_temp
      param_hamil%automatic_translation = .false.
    endif
  end subroutine param_read_hamil

  subroutine param_read_bloch_phase(use_bloch_phases, disentanglement, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: stdout
    logical, intent(out) :: use_bloch_phases
    logical, intent(in) :: disentanglement
    character(len=50), intent(in)  :: seedname

    logical :: found

    use_bloch_phases = .false.
    call param_get_keyword(stdout, seedname, 'use_bloch_phases', found, l_value=use_bloch_phases)
    if (disentanglement .and. use_bloch_phases) &
      call io_error('Error: Cannot use bloch phases for disentanglement', stdout, seedname)
  end subroutine param_read_bloch_phase

  subroutine param_read_explicit_kpts(library, driver, kmesh_info, num_kpts, bohr, stdout, seedname)
    use w90_io, only: io_error
    use w90_utility, only: utility_recip_lattice
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: library
    type(param_driver_type), intent(inout) :: driver
    type(kmesh_info_type), intent(inout) :: kmesh_info
    integer, intent(in) :: num_kpts
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname

    integer :: i, k, ierr, rows
    logical :: found
    integer, allocatable, dimension(:, :) :: nnkpts_block
    integer, allocatable, dimension(:) :: nnkpts_idx

    ! get the nnkpts block -- this is allowed only in postproc-setup mode
    call param_get_block_length(stdout, seedname, 'nnkpts', driver%explicit_nnkpts, rows, library)
    if (driver%explicit_nnkpts) then
      kmesh_info%nntot = rows/num_kpts
      if (modulo(rows, num_kpts) /= 0) then
        call io_error('The number of rows in nnkpts must be a multiple of num_kpts', stdout, seedname)
      end if
      if (allocated(nnkpts_block)) deallocate (nnkpts_block)
      allocate (nnkpts_block(5, rows), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nnkpts_block in param_read', stdout, seedname)
      call param_get_keyword_block(stdout, seedname, 'nnkpts', found, rows, 5, bohr, i_value=nnkpts_block)
      ! check that postproc_setup is true
      if (.not. driver%postproc_setup) &
        call io_error('Input parameter nnkpts_block is allowed only if postproc_setup = .true.', stdout, seedname)
      ! assign the values in nnkpts_block to nnlist and nncell
      ! this keeps track of how many neighbours have been seen for each k-point
      if (allocated(nnkpts_idx)) deallocate (nnkpts_idx)
      allocate (nnkpts_idx(num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nnkpts_idx in param_read', stdout, seedname)
      nnkpts_idx = 1
      ! allocating "global" nnlist & nncell
      ! These are deallocated in kmesh_dealloc
      if (allocated(kmesh_info%nnlist)) deallocate (kmesh_info%nnlist)
      allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nnlist in param_read', stdout, seedname)
      if (allocated(kmesh_info%nncell)) deallocate (kmesh_info%nncell)
      allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating nncell in param_read', stdout, seedname)
      do i = 1, num_kpts*kmesh_info%nntot
        k = nnkpts_block(1, i)
        kmesh_info%nnlist(k, nnkpts_idx(k)) = nnkpts_block(2, i)
        kmesh_info%nncell(:, k, nnkpts_idx(k)) = nnkpts_block(3:, i)
        nnkpts_idx(k) = nnkpts_idx(k) + 1
      end do
      ! check that all k-points have the same number of neighbours
      if (any(nnkpts_idx /= (/(kmesh_info%nntot + 1, i=1, num_kpts)/))) then
        call io_error('Inconsistent number of nearest neighbours.', stdout, seedname)
      end if
      deallocate (nnkpts_idx, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating nnkpts_idx in param_read', stdout, seedname)
      deallocate (nnkpts_block, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating nnkpts_block in param_read', stdout, seedname)
    end if

  end subroutine param_read_explicit_kpts

  subroutine param_read_projections(proj, use_bloch_phases, lhasproj, guiding_centres, &
                                    proj_site, kmesh_data, select_proj, num_proj, &
                                    param_input, atoms, recip_lattice, num_wann, library, &
                                    bohr, stdout, seedname)
    use w90_io, only: io_error
    implicit none
    type(projection_type), intent(inout) :: proj
    logical, intent(in) :: use_bloch_phases, guiding_centres, library
    logical, intent(out) :: lhasproj
    real(kind=dp), allocatable, dimension(:, :), intent(out) :: proj_site
    type(param_kmesh_type), intent(inout) :: kmesh_data
    type(select_projection_type), intent(inout) :: select_proj
    integer, intent(inout) :: num_proj
    type(parameter_input_type), intent(in) :: param_input
    type(atom_data_type), intent(in) :: atoms
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname

    integer, intent(in) :: stdout
    integer :: i, j, i_temp, loop, ierr
    logical :: found
    ! projections selection
    integer :: num_select_projections
    integer, allocatable :: select_projections(:)

    ! Projections
    kmesh_data%auto_projections = .false.
    call param_get_keyword(stdout, seedname, 'auto_projections', found, l_value=kmesh_data%auto_projections)
    num_proj = 0
    call param_get_block_length(stdout, seedname, 'projections', found, i_temp, library)
    ! check to see that there are no unrecognised keywords
    if (found) then
      if (kmesh_data%auto_projections) call io_error('Error: Cannot specify both auto_projections and projections block', &
                                                     stdout, seedname)
      lhasproj = .true.
      call param_get_projections(num_proj, atoms, kmesh_data, param_input, &
                                 num_wann, proj_site, proj, recip_lattice, &
                                 .true., bohr, stdout, seedname)
    else
      if (guiding_centres .and. .not. (param_input%gamma_only .and. use_bloch_phases)) &
        call io_error('param_read: Guiding centres requested, but no projection block found', stdout, seedname)
      lhasproj = .false.
      num_proj = num_wann
    end if

    select_proj%lselproj = .false.
    num_select_projections = 0
    call param_get_range_vector(stdout, seedname, 'select_projections', found, num_select_projections, lcount=.true.)
    if (found) then
      if (num_select_projections < 1) call io_error('Error: problem reading select_projections', stdout, seedname)
      if (allocated(select_projections)) deallocate (select_projections)
      allocate (select_projections(num_select_projections), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating select_projections in param_read', stdout, seedname)
      call param_get_range_vector(stdout, seedname, 'select_projections', found, &
                                  num_select_projections, .false., select_projections)
      if (any(select_projections < 1)) &
        call io_error('Error: select_projections must contain positive numbers', stdout, seedname)
      if (num_select_projections < num_wann) &
        call io_error('Error: too few projections selected', stdout, seedname)
      if (num_select_projections > num_wann) &
        call io_error('Error: too many projections selected', stdout, seedname)
      if (.not. lhasproj) &
        call io_error('Error: select_projections cannot be used without defining the projections', stdout, seedname)
      if (maxval(select_projections(:)) > num_proj) &
        call io_error('Error: select_projections contains a number greater than num_proj', stdout, seedname)
      select_proj%lselproj = .true.
    end if

    if (allocated(select_proj%proj2wann_map)) deallocate (select_proj%proj2wann_map)
    allocate (select_proj%proj2wann_map(num_proj), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating proj2wann_map in param_read', stdout, seedname)
    select_proj%proj2wann_map = -1

    if (select_proj%lselproj) then
      do i = 1, num_proj
        do j = 1, num_wann
          if (select_projections(j) == i) select_proj%proj2wann_map(i) = j
        enddo
      enddo
    else
      do i = 1, num_wann
        select_proj%proj2wann_map(i) = i
      enddo
    endif

    if (lhasproj) then
      call param_get_projections(num_proj, atoms, kmesh_data, param_input, &
                                 num_wann, proj_site, proj, &
                                 recip_lattice, .false., bohr, stdout, seedname)
      do loop = 1, num_proj
        if (select_proj%proj2wann_map(loop) < 0) cycle
        proj_site(:, select_proj%proj2wann_map(loop)) = kmesh_data%input_proj_site(:, loop)
        proj%l(select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%l(loop)
        proj%m(select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%m(loop)
        proj%z(:, select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%z(:, loop)
        proj%x(:, select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%x(:, loop)
        proj%radial(select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%radial(loop)
        proj%zona(select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%zona(loop)
      enddo

      if (param_input%spinors) then
        do loop = 1, num_proj
          if (select_proj%proj2wann_map(loop) < 0) cycle
          proj%s(select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%s(loop)
          proj%s_qaxis(:, select_proj%proj2wann_map(loop)) = kmesh_data%input_proj%s_qaxis(:, loop)
        enddo
      endif
    endif

  end subroutine param_read_projections

  subroutine param_read_constrained_centres(ccentres_frac, param_wannierise, real_lattice, &
                                            num_wann, library, stdout, seedname)
!   use w90_io, only: io_error, stdout
    use w90_io, only: io_error
    implicit none
    real(kind=dp), intent(inout) :: ccentres_frac(:, :)
    type(param_wannierise_type), intent(inout) :: param_wannierise
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    logical, intent(in) :: library
    character(len=50), intent(in)  :: seedname

    integer :: i_temp
    logical :: found

    ! Constrained centres
    call param_get_block_length(stdout, seedname, 'slwf_centres', found, i_temp, library)
    if (found) then
      if (param_wannierise%slwf_constrain) then
        ! Allocate array for constrained centres
        call param_get_centre_constraints(ccentres_frac, &
                                          param_wannierise%ccentres_cart, &
                                          param_wannierise%proj_site, &
                                          num_wann, real_lattice, stdout, seedname)
      else
        write (stdout, '(a)') ' slwf_constrain set to false. Ignoring <slwf_centres> block '
      end if
      ! Check that either projections or constrained centres are specified if slwf_constrain=.true.
    elseif (.not. found) then
      if (param_wannierise%slwf_constrain) then
        if (.not. allocated(param_wannierise%proj_site)) then
          call io_error('Error: slwf_constrain = true, but neither &
               & <slwf_centre> block  nor &
               & <projection_block> are specified.', stdout, seedname)
        else
          ! Allocate array for constrained centres
          call param_get_centre_constraints(ccentres_frac, &
                                            param_wannierise%ccentres_cart, &
                                            param_wannierise%proj_site, &
                                            num_wann, real_lattice, stdout, seedname)
        end if
      end if
    end if
    ! Warning
    if (param_wannierise%slwf_constrain .and. allocated(param_wannierise%proj_site) .and. .not. found) &
         & write (stdout, '(a)') ' Warning: No <slwf_centres> block found, but slwf_constrain set to true. &
           & Desired centres for SLWF same as projection centres.'
  end subroutine param_read_constrained_centres

!===================================================================
  subroutine param_write(driver, w90_calcs, param_input, param_plot, param_wannierise, &
                         lsitesymmetry, symmetrize_eps, wann_data, param_hamil, kmesh_data, &
                         k_points, num_kpts, dis_data, fermi_surface_data, fermi, tran, atoms, &
                         num_bands, num_wann, mp_grid, num_proj, select_proj, real_lattice, &
                         recip_lattice, spec_points, stdout, write_data, proj)
    !==================================================================!
    !                                                                  !
    !! write wannier90 parameters to stdout
    !                                                                  !
    !===================================================================

    implicit none

    !data from parameters module
    type(param_driver_type), intent(in) :: driver
    type(w90_calculation_type), intent(in) :: w90_calcs
    type(parameter_input_type), intent(in) :: param_input
    type(param_plot_type), intent(in) :: param_plot
    type(param_wannierise_type), intent(in) :: param_wannierise
    ! RS: symmetry-adapted Wannier functions
    logical, intent(in) :: lsitesymmetry
    real(kind=dp), intent(in) :: symmetrize_eps
    type(wannier_data_type), intent(in) :: wann_data
    type(param_hamiltonian_type), intent(in) :: param_hamil
    type(param_kmesh_type), intent(in) :: kmesh_data
    type(k_point_type), intent(in) :: k_points
    integer, intent(in) :: num_kpts
    type(disentangle_type), intent(in) :: dis_data
    type(fermi_surface_type), intent(in) :: fermi_surface_data
    type(fermi_data_type), intent(in) :: fermi
    type(transport_type), intent(in) :: tran
    type(atom_data_type), intent(in) :: atoms
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_proj
    type(select_projection_type), intent(in) :: select_proj
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    type(special_kpoints_type), intent(in) :: spec_points
    !type(pw90_calculation_type), intent(in) :: pw90_calcs
    type(w90_extra_io_type), intent(in) :: write_data
    type(projection_type), intent(in) :: proj

    integer :: i, nkp, loop, nat, nsp
    real(kind=dp) :: cell_volume

    if (w90_calcs%transport .and. tran%read_ht) goto 401

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
    ! Constrained centres
    if (param_wannierise%selective_loc .and. param_wannierise%slwf_constrain) then
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '| Wannier#        Original Centres              Constrained centres          |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do i = 1, param_wannierise%slwf_num
        write (stdout, '(1x,a1,2x,i3,2x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                    '|', i, write_data%ccentres_frac(i, :), '|', wann_data%centres(:, i), '|'
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
    end if
    ! Projections
    if (param_input%iprint > 1 .and. allocated(kmesh_data%input_proj_site)) then
      write (stdout, '(32x,a)') '-----------'
      write (stdout, '(32x,a)') 'PROJECTIONS'
      write (stdout, '(32x,a)') '-----------'
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a)') '|     Frac. Coord.   l mr  r        z-axis               x-axis          Z/a |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, num_proj
        write (stdout, '(1x,a1,3(1x,f5.2),1x,i2,1x,i2,1x,i2,3(1x,f6.3),3(1x,f6.3),2x,f4.1,1x,a1)') &
          '|', kmesh_data%input_proj_site(1, nsp), kmesh_data%input_proj_site(2, nsp), &
          kmesh_data%input_proj_site(3, nsp), kmesh_data%input_proj%l(nsp), &
          kmesh_data%input_proj%m(nsp), kmesh_data%input_proj%radial(nsp), &
          kmesh_data%input_proj%z(1, nsp), kmesh_data%input_proj%z(2, nsp), &
          kmesh_data%input_proj%z(3, nsp), kmesh_data%input_proj%x(1, nsp), &
          kmesh_data%input_proj%x(2, nsp), kmesh_data%input_proj%x(3, nsp), &
          kmesh_data%input_proj%zona(nsp), '|'
      end do
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, *) ' '
    end if

    if (param_input%iprint > 1 .and. select_proj%lselproj .and. allocated(param_wannierise%proj_site)) then
      write (stdout, '(30x,a)') '--------------------'
      write (stdout, '(30x,a)') 'SELECTED PROJECTIONS'
      write (stdout, '(30x,a)') '--------------------'
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a)') '|     Frac. Coord.   l mr  r        z-axis               x-axis          Z/a |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, num_wann
        if (select_proj%proj2wann_map(nsp) < 0) cycle
        write (stdout, '(1x,a1,3(1x,f5.2),1x,i2,1x,i2,1x,i2,3(1x,f6.3),3(1x,f6.3),2x,f4.1,1x,a1)')&
  &              '|', param_wannierise%proj_site(1, nsp), param_wannierise%proj_site(2, nsp), &
             param_wannierise%proj_site(3, nsp), proj%l(nsp), proj%m(nsp), proj%radial(nsp), &
             proj%z(1, nsp), proj%z(2, nsp), proj%z(3, nsp), proj%x(1, nsp), &
             proj%x(2, nsp), proj%x(3, nsp), proj%zona(nsp), '|'
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
    if (param_input%iprint > 1) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (param_input%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |'
      else
        write (stdout, '(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nkp = 1, num_kpts
        write (stdout, '(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|', nkp, k_points%kpt_latt(:, nkp), '|', &
          k_points%kpt_cart(:, nkp)/param_input%lenconfac, '|'
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, *) ' '
    end if
    ! Main
    write (stdout, *) ' '
    write (stdout, '(1x,a78)') '*---------------------------------- MAIN ------------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Wannier Functions               :', num_wann, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Objective Wannier Functions     :', param_wannierise%slwf_num, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of input Bloch states              :', num_bands, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :', param_input%iprint, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :', param_input%timing_level, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Optimisation (0=memory, 3=speed)          :', param_input%optimisation, '|'
    write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :', trim(param_input%length_unit), '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Post-processing setup (write *.nnkp)      :', driver%postproc_setup, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using Gamma-only branch of algorithms     :', param_input%gamma_only, '|'
    !YN: RS:
    if (lsitesymmetry) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using symmetry-adapted WF mode            :', lsitesymmetry, '|'
      write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Tolerance for symmetry condition on U     :', symmetrize_eps, '|'
    endif

    if (w90_calcs%cp_pp .or. param_input%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  CP code post-processing                   :', w90_calcs%cp_pp, '|'
    if (w90_calcs%wannier_plot .or. param_input%iprint > 2) then
      if (param_plot%wvfn_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  Wavefunction (UNK) file-type              :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  Wavefunction (UNK) file-type              :', 'unformatted', '|'
      endif
      if (param_plot%spin == 1) then
        write (stdout, '(1x,a46,16x,a2,13x,a1)') '|  Wavefunction spin channel                 :', 'up', '|'
      else
        write (stdout, '(1x,a46,14x,a4,13x,a1)') '|  Wavefunction spin channel                 :', 'down', '|'
      endif
    endif

    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    ! Wannierise
    write (stdout, '(1x,a78)') '*------------------------------- WANNIERISE ---------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Total number of iterations                :', param_wannierise%num_iter, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of CG steps before reset           :', param_wannierise%num_cg_steps, '|'
    if (param_wannierise%lfixstep) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fixed step length for minimisation        :', &
        param_wannierise%fixed_step, '|'
    else
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Trial step length for line search         :', &
        param_wannierise%trial_step, '|'
    endif
    write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Convergence tolerence                     :', param_wannierise%conv_tol, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Convergence window                        :', param_wannierise%conv_window, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between writing output         :', &
      param_wannierise%num_print_cycles, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between backing up to disk     :', &
      param_wannierise%num_dump_cycles, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write r^2_nm to file                      :', param_wannierise%write_r2mn, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write xyz WF centres to file              :', param_input%write_xyz, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write on-site energies <0n|H|0n> to file  :', param_wannierise%write_hr_diag, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use guiding centre to control phases      :', param_wannierise%guiding_centres, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use phases for initial projections        :', w90_calcs%use_bloch_phases, '|'
    if (param_wannierise%guiding_centres .or. param_input%iprint > 2) then
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations before starting guiding centres:', &
        param_wannierise%num_no_guide_iter, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between using guiding centres  :', &
        param_wannierise%num_guide_cycles, '|'
    end if
    if (param_wannierise%selective_loc .or. param_input%iprint > 2) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Perform selective localization            :', param_wannierise%selective_loc, '|'
    end if
    if (param_wannierise%slwf_constrain .or. param_input%iprint > 2) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use constrains in selective localization  :', &
        param_wannierise%slwf_constrain, '|'
      write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Value of the Lagrange multiplier          :',&
           &param_wannierise%slwf_lambda, '|'
    end if
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    !
    ! Disentanglement
    !
    if (w90_calcs%disentanglement .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------------- DISENTANGLE --------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using band disentanglement                :', w90_calcs%disentanglement, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Total number of iterations                :', dis_data%num_iter, '|'
      write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|  Mixing ratio                              :', dis_data%mix_ratio, '|'
      write (stdout, '(1x,a46,8x,ES10.3,13x,a1)') '|  Convergence tolerence                     :', dis_data%conv_tol, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Convergence window                        :', dis_data%conv_window, '|'
      ! GS-start
      if (dis_data%spheres_num .gt. 0) then
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of spheres in k-space              :', dis_data%spheres_num, '|'
        do nkp = 1, dis_data%spheres_num
          write (stdout, '(1x,a13,I4,a2,2x,3F8.3,a15,F8.3,9x,a1)') &
            '|   center n.', nkp, ' :', dis_data%spheres(1:3, nkp), ',    radius   =', dis_data%spheres(4, nkp), '|'
        enddo
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Index of first Wannier band               :', &
          dis_data%spheres_first_wann, '|'
      endif
      ! GS-end
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    end if
    !
    ! Plotting
    !
    if (w90_calcs%wannier_plot .or. w90_calcs%bands_plot .or. w90_calcs%fermi_surface_plot &
        .or. w90_calcs%write_hr .or. param_input%iprint > 2) then
      !
      write (stdout, '(1x,a78)') '*-------------------------------- PLOTTING ----------------------------------*'
      !
      if (w90_calcs%wannier_plot .or. param_input%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Wannier functions                :', w90_calcs%wannier_plot, '|'
        write (stdout, '(1x,a46,1x,I5,a1,I5,a1,I5,13x,a1)') &
          '|   Size of supercell for plotting           :', &
          param_plot%wannier_plot_supercell(1), 'x', param_plot%wannier_plot_supercell(2), 'x', &
          param_plot%wannier_plot_supercell(3), '|'

        if (param_wannierise%translate_home_cell) then
          write (stdout, '(1x,a46,10x,L8,13x,a1)') &
            '|  Translating WFs to home cell              :', param_wannierise%translate_home_cell, '|'
        end if

        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting mode (molecule or crystal)      :', &
          trim(param_plot%wannier_plot_mode), '|'
        if (param_input%spinors) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting mode for spinor WFs             :', &
            trim(param_plot%wannier_plot_spinor_mode), '|'
          write (stdout, '(1x,a46,10x,L8,13x,a1)') '|   Include phase for spinor WFs             :', &
            param_plot%wannier_plot_spinor_phase, '|'
        end if
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :', &
          trim(param_plot%wannier_plot_format), '|'
        if (index(param_plot%wannier_plot_format, 'cub') > 0 .or. param_input%iprint > 2) then
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Plot radius                              :', &
            param_plot%wannier_plot_radius, '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Plot scale                               :', &
            param_plot%wannier_plot_scale, '|'
        endif
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (w90_calcs%fermi_surface_plot .or. param_input%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Fermi surface                    :', w90_calcs%fermi_surface_plot, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Number of plotting points (along b_1)    :', &
          fermi_surface_data%num_points, '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :' &
          , trim(fermi_surface_data%plot_format), '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (w90_calcs%bands_plot .or. param_input%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting interpolated bandstructure       :', w90_calcs%bands_plot, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Number of K-path sections                :', &
          spec_points%bands_num_spec_points/2, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Divisions along first K-path section     :', &
          param_plot%bands_num_points, '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Output format                            :', &
          trim(param_plot%bands_plot_format), '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Output mode                              :', &
          trim(param_input%bands_plot_mode), '|'
        if (index(param_input%bands_plot_mode, 'cut') .ne. 0) then
          write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Dimension of the system                  :', &
            param_plot%bands_plot_dim, '|'
          if (param_plot%bands_plot_dim .eq. 1) &
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System extended in                       :', &
            trim(write_data%one_dim_axis), '|'
          if (param_plot%bands_plot_dim .eq. 2) &
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System confined in                       :', &
            trim(write_data%one_dim_axis), '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off value                :', &
            param_input%hr_cutoff, '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off distance             :', &
            param_input%dist_cutoff, '|'
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian cut-off distance mode        :', &
            trim(param_input%dist_cutoff_mode), '|'
        endif
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
        write (stdout, '(1x,a78)') '|   K-space path sections:                                                   |'
        if (spec_points%bands_num_spec_points == 0) then
          write (stdout, '(1x,a78)') '|     None defined                                                           |'
        else
          do loop = 1, spec_points%bands_num_spec_points, 2
            write (stdout, '(1x,a10,1x,a5,1x,3F7.3,5x,a3,1x,a5,1x,3F7.3,3x,a1)') '|    From:', &
              spec_points%bands_label(loop), (spec_points%bands_spec_points(i, loop), i=1, 3), &
              'To:', spec_points%bands_label(loop + 1), (spec_points%bands_spec_points(i, loop + 1), i=1, 3), '|'
          end do
        end if
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (w90_calcs%write_hr .or. param_input%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Hamiltonian in WF basis          :', w90_calcs%write_hr, '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      endif
      if (param_wannierise%write_vdw_data .or. param_input%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Writing data for Van der Waals post-proc  :', &
          param_wannierise%write_vdw_data, '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      endif
      !
    endif

401 continue
    !
    ! Transport
    !
    if (w90_calcs%transport .or. param_input%iprint > 2) then
      !
      write (stdout, '(1x,a78)') '*------------------------------- TRANSPORT ----------------------------------*'
      !
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Transport mode                            :', trim(tran%mode), '|'
      !
      if (tran%read_ht) then
        !
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian from external files          :', 'T', '|'
        !
      else
        !
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian from external files          :', 'F', '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System extended in                       :', &
          trim(write_data%one_dim_axis), '|'
        !
      end if

      write (stdout, '(1x,a78)') '|   Centre of the unit cell to which WF are translated (fract. coords):      |'
      write (stdout, '(1x,a1,35x,F12.6,a1,F12.6,a1,F12.6,3x,a1)') '|', param_hamil%translation_centre_frac(1), ',', &
        param_hamil%translation_centre_frac(2), ',', &
        param_hamil%translation_centre_frac(3), '|'

      if (size(fermi%energy_list) == 1) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi energy (eV)                         :', fermi%energy_list(1), '|'
      else
        write (stdout, '(1x,a21,I8,a12,f8.3,a4,f8.3,a3,13x,a1)') '|  Fermi energy     :', size(fermi%energy_list), &
          ' steps from ', fermi%energy_list(1), ' to ', &
          fermi%energy_list(size(fermi%energy_list)), ' eV', '|'
      end if
      !
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      !
    endif

101 format(20x, a3, 2x, 3F11.6)

  end subroutine param_write

  subroutine param_w90_dealloc(param_input, param_plot, param_wannierise, &
                               wann_data, kmesh_data, k_points, dis_data, &
                               atoms, eigval, spec_points, stdout, seedname, write_data, proj)
    use w90_io, only: io_error
    implicit none
    !data from parameters module
    !type(param_driver_type), intent(inout) :: driver
    integer, intent(in) :: stdout
    type(parameter_input_type), intent(inout) :: param_input
    type(param_plot_type), intent(inout) :: param_plot
    type(param_wannierise_type), intent(inout) :: param_wannierise
    type(wannier_data_type), intent(inout) :: wann_data
    type(param_kmesh_type), intent(inout) :: kmesh_data
    type(k_point_type), intent(inout) :: k_points
    type(disentangle_type), intent(inout) :: dis_data
    !type(fermi_data_type), intent(inout) :: fermi
    type(atom_data_type), intent(inout) :: atoms
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    type(special_kpoints_type), intent(inout) :: spec_points
    character(len=50), intent(in)  :: seedname
    !type(dos_plot_type), intent(inout) :: dos_data
    !type(berry_type), intent(inout) :: berry
    type(w90_extra_io_type), intent(inout) :: write_data
    type(projection_type), intent(inout) :: proj

    integer :: ierr

    call param_dealloc(param_input, wann_data, kmesh_data, k_points, &
                       dis_data, atoms, eigval, spec_points, stdout, seedname)
    if (allocated(param_plot%wannier_plot_list)) then
      deallocate (param_plot%wannier_plot_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_plot_list in param_dealloc', stdout, seedname)
    end if
    if (allocated(param_plot%bands_plot_project)) then
      deallocate (param_plot%bands_plot_project, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating bands_plot_project in param_dealloc', stdout, seedname)
    endif
    if (allocated(write_data%ccentres_frac)) then
      deallocate (write_data%ccentres_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating ccentres_frac in param_w90_dealloc', stdout, seedname)
    endif
    if (allocated(param_wannierise%proj_site)) then
      deallocate (param_wannierise%proj_site, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_site in param_dealloc', stdout, seedname)
    end if
    if (allocated(param_wannierise%ccentres_cart)) then
      deallocate (param_wannierise%ccentres_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating ccentres_cart in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%l)) then
      deallocate (proj%l, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_l in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%m)) then
      deallocate (proj%m, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_m in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%s)) then
      deallocate (proj%s, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_s in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%s_qaxis)) then
      deallocate (proj%s_qaxis, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_s_qaxis in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%z)) then
      deallocate (proj%z, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_z in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%x)) then
      deallocate (proj%x, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_x in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%radial)) then
      deallocate (proj%radial, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_radial in param_dealloc', stdout, seedname)
    end if
    if (allocated(proj%zona)) then
      deallocate (proj%zona, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating proj_zona in param_dealloc', stdout, seedname)
    end if
  end subroutine param_w90_dealloc

!=================================================!
  subroutine param_write_chkpt(chkpt, param_input, wann_data, kmesh_info, &
                               k_points, num_kpts, dis_data, num_bands, &
                               num_wann, u_matrix, u_matrix_opt, m_matrix, &
                               mp_grid, real_lattice, recip_lattice, stdout, seedname)
    !=================================================!
    !! Write checkpoint file
    !! IMPORTANT! If you change the chkpt format, adapt
    !! accordingly also the w90chk2chk.x utility!
    !! Also, note that this routine writes the u_matrix and the m_matrix - in parallel
    !! mode these are however stored in distributed form in, e.g., u_matrix_loc only, so
    !! if you are changing the u_matrix, remember to gather it from u_matrix_loc first!
    !=================================================!

!   use w90_io, only: io_file_unit, io_date, seedname
    use w90_io, only: io_file_unit, io_date

    implicit none

    character(len=*), intent(in) :: chkpt
    !data from parameters module
    type(parameter_input_type), intent(in) :: param_input
    type(wannier_data_type), intent(in) :: wann_data
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(k_point_type), intent(in) :: k_points
    integer, intent(in) :: num_kpts
    type(disentangle_type), intent(in) :: dis_data
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    integer, intent(in) :: mp_grid(3)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    character(len=50), intent(in)  :: seedname

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
    write (chk_unit) param_input%num_exclude_bands            ! Number of excluded bands
    write (chk_unit) (param_input%exclude_bands(i), i=1, param_input%num_exclude_bands) ! Excluded bands
    write (chk_unit) ((real_lattice(i, j), i=1, 3), j=1, 3)        ! Real lattice
    write (chk_unit) ((recip_lattice(i, j), i=1, 3), j=1, 3)       ! Reciprocal lattice
    write (chk_unit) num_kpts                                 ! Number of k-points
    write (chk_unit) (mp_grid(i), i=1, 3)                       ! M-P grid
    write (chk_unit) ((k_points%kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts) ! K-points
    write (chk_unit) kmesh_info%nntot                  ! Number of nearest k-point neighbours
    write (chk_unit) num_wann               ! Number of wannier functions
    chkpt1 = adjustl(trim(chkpt))
    write (chk_unit) chkpt1                 ! Position of checkpoint
    write (chk_unit) param_input%have_disentangled      ! Whether a disentanglement has been performed
    if (param_input%have_disentangled) then
      write (chk_unit) param_input%omega_invariant     ! Omega invariant
      ! lwindow, ndimwin and U_matrix_opt
      write (chk_unit) ((dis_data%lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
      write (chk_unit) (dis_data%ndimwin(nkp), nkp=1, num_kpts)
      write (chk_unit) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
    endif
    write (chk_unit) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)               ! U_matrix
    write (chk_unit) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, kmesh_info%nntot), l=1, num_kpts) ! M_matrix
    write (chk_unit) ((wann_data%centres(i, j), i=1, 3), j=1, num_wann)
    write (chk_unit) (wann_data%spreads(i), i=1, num_wann)
    close (chk_unit)

    write (stdout, '(a/)') ' done'

    return

  end subroutine param_write_chkpt

!===========================================!
  subroutine param_memory_estimate(w90_calcs, param_input, param_wannierise, &
                                   kmesh_data, kmesh_info, num_kpts, &
                                   atoms, num_bands, num_wann, num_proj, stdout)
    !===========================================!
    !                                           !
    !! Estimate how much memory we will allocate
    !                                           !
    !===========================================!

    !use w90_comms, only:

    implicit none

    !data from parameters module
    type(w90_calculation_type), intent(in) :: w90_calcs
    type(parameter_input_type), intent(in) :: param_input
    type(param_wannierise_type), intent(in) :: param_wannierise
    type(param_kmesh_type), intent(in) :: kmesh_data
    type(kmesh_info_type), intent(in) :: kmesh_info
    integer, intent(in) :: num_kpts
    !type(disentangle_type), intent(in) :: dis_data
    type(atom_data_type), intent(in) :: atoms
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: num_proj
    !type(pw90_calculation_type), intent(in) :: pw90_calcs
    !type(postw90_common_type), intent(in) :: pw90_common
    !type(boltzwann_type), intent(in) :: boltz
    !logical, intent(in) :: ispostw90 ! Are we running postw90?

    real(kind=dp), parameter :: size_log = 1.0_dp
    real(kind=dp), parameter :: size_int = 4.0_dp
    real(kind=dp), parameter :: size_real = 8.0_dp
    real(kind=dp), parameter :: size_cmplx = 16.0_dp
    real(kind=dp) :: mem_wan, mem_wan1, mem_param, mem_dis, mem_dis2, mem_dis1
    real(kind=dp) :: mem_bw
    !integer :: NumPoints1, NumPoints2, NumPoints3, ndim
    !real(kind=dp) :: TDF_exceeding_energy

    mem_param = 0
    mem_dis = 0
    mem_dis1 = 0
    mem_dis2 = 0
    mem_wan = 0
    mem_wan1 = 0
    mem_bw = 0

    ! First the data stored in the parameters module
    mem_param = mem_param + num_wann*num_wann*num_kpts*size_cmplx                   !u_matrix
    if (.not. w90_calcs%disentanglement) &
      mem_param = mem_param + num_wann*num_wann*kmesh_info%nntot*num_kpts*size_cmplx       !m_matrix

    if (w90_calcs%disentanglement) then
      mem_param = mem_param + num_bands*num_wann*num_kpts*size_cmplx             ! u_matrix_opt
    endif

    if (allocated(atoms%species_num)) then
      mem_param = mem_param + (atoms%num_species)*size_int                               !atoms_species_num
      mem_param = mem_param + (atoms%num_species)*size_real                              !atoms_label
      mem_param = mem_param + (atoms%num_species)*size_real                              !atoms_symbol
      mem_param = mem_param + (3*maxval(atoms%species_num)*atoms%num_species)*size_real  !atoms_pos_frac
      mem_param = mem_param + (3*maxval(atoms%species_num)*atoms%num_species)*size_real  !atoms_pos_cart
    endif

    if (allocated(kmesh_data%input_proj_site)) then
      mem_param = mem_param + (3*num_proj)*size_real              !input_proj_site
      mem_param = mem_param + (num_proj)*size_int                !input_proj_l
      mem_param = mem_param + (num_proj)*size_int                 !input_proj_m
      mem_param = mem_param + (3*num_proj)*size_real             !input_proj_z
      mem_param = mem_param + (3*num_proj)*size_real             !input_proj_x
      mem_param = mem_param + (num_proj)*size_real                !input_proj_radial
      mem_param = mem_param + (num_proj)*size_real                !input_proj_zona
    endif

    if (allocated(param_wannierise%proj_site)) then
      mem_param = mem_param + (3*num_wann)*size_real              !proj_site
      mem_param = mem_param + (num_wann)*size_int                !proj_l
      mem_param = mem_param + (num_wann)*size_int                 !proj_m
      mem_param = mem_param + (3*num_wann)*size_real             !proj_z
      mem_param = mem_param + (3*num_wann)*size_real             !proj_x
      mem_param = mem_param + (num_wann)*size_real                !proj_radial
      mem_param = mem_param + (num_wann)*size_real                !proj_zona
    endif

    mem_param = mem_param + num_kpts*kmesh_info%nntot*size_int    !nnlist
    mem_param = mem_param + num_kpts*kmesh_info%nntot/2*size_int  !neigh
    mem_param = mem_param + 3*num_kpts*kmesh_info%nntot*size_int  !nncell
    mem_param = mem_param + kmesh_info%nntot*size_real            !wb
    mem_param = mem_param + 3*kmesh_info%nntot/2*size_real        !bka
    mem_param = mem_param + 3*kmesh_info%nntot*num_kpts*size_real !bk

    mem_param = mem_param + num_bands*num_kpts*size_real             !eigval
    mem_param = mem_param + 3*num_kpts*size_real                     !kpt_cart
    mem_param = mem_param + 3*num_kpts*size_real                     !kpt_latt
    if (w90_calcs%disentanglement) then
      mem_param = mem_param + num_kpts*size_int                     !ndimwin
      mem_param = mem_param + num_bands*num_kpts*size_log           !lwindow
    endif
    mem_param = mem_param + 3*num_wann*size_real                     !wannier_centres
    mem_param = mem_param + num_wann*size_real                       !wannier_spreads

    if (w90_calcs%disentanglement) then
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
      if (param_input%gamma_only) then
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
      mem_dis2 = mem_dis2 + num_wann*num_wann*kmesh_info%nntot*num_kpts*size_cmplx!m_matrix

      if (param_input%optimisation <= 0) then
        mem_dis = mem_dis + mem_dis1
      else
        mem_dis = mem_dis + max(mem_dis1, mem_dis2)
      endif

      mem_dis = mem_dis + num_bands*num_bands*kmesh_info%nntot*num_kpts*size_cmplx      ! m_matrix_orig
      mem_dis = mem_dis + num_bands*num_wann*num_kpts*size_cmplx             ! a_matrix

    endif

    !Wannierise

    mem_wan1 = mem_wan1 + (num_wann*num_wann*kmesh_info%nntot*num_kpts)*size_cmplx     !  'm0'
    if (param_input%optimisation > 0) then
      mem_wan = mem_wan + mem_wan1
    endif
    mem_wan = mem_wan + (num_wann*num_wann*num_kpts)*size_cmplx           !  'u0'
    mem_wan = mem_wan + (num_wann*kmesh_info%nntot*num_kpts)*size_real    !  'rnkb'
    mem_wan = mem_wan + (num_wann*kmesh_info%nntot*num_kpts)*size_real     !  'ln_tmp'
    mem_wan = mem_wan + (num_wann*kmesh_info%nntot*num_kpts)*size_cmplx    !  'csheet'
    mem_wan = mem_wan + (num_wann*kmesh_info%nntot*num_kpts)*size_real     !  'sheet'
    mem_wan = mem_wan + (3*num_wann)*size_real                             !  'rave'
    mem_wan = mem_wan + (num_wann)*size_real                              !  'r2ave'
    mem_wan = mem_wan + (num_wann)*size_real                               !  'rave2'
    mem_wan = mem_wan + (3*num_wann)*size_real                            !  'rguide'
    mem_wan = mem_wan + (num_wann*num_wann)*size_cmplx                  !  'cz'
    if (param_input%gamma_only) then
      mem_wan = mem_wan + num_wann*num_wann*kmesh_info%nntot*2*size_cmplx    ! m_w
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

    if (w90_calcs%disentanglement) &
      mem_wan = mem_wan + num_wann*num_wann*kmesh_info%nntot*num_kpts*size_cmplx       !m_matrix

    if (param_input%iprint > 0) then
      write (stdout, '(1x,a)') '*============================================================================*'
      write (stdout, '(1x,a)') '|                              MEMORY ESTIMATE                               |'
      write (stdout, '(1x,a)') '|         Maximum RAM allocated during each phase of the calculation         |'
      write (stdout, '(1x,a)') '*============================================================================*'
      if (w90_calcs%disentanglement) &
        write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Disentanglement:', (mem_param + mem_dis)/(1024**2), ' Mb'
      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:', (mem_param + mem_wan)/(1024**2), ' Mb'
      if (param_input%optimisation > 0 .and. param_input%iprint > 1) then
        write (stdout, '(1x,a)') '|                                                                            |'
        write (stdout, '(1x,a)') '|   N.B. by setting optimisation=0 memory usage will be reduced to:          |'
        if (w90_calcs%disentanglement) &
          write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Disentanglement:', &
          (mem_param + mem_dis - max(mem_dis1, mem_dis2) + mem_dis1)/(1024**2), ' Mb'
        if (param_input%gamma_only) then
          write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:', (mem_param + mem_wan)/(1024**2), ' Mb'
        else
          write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:', &
            (mem_param + mem_wan - mem_wan1)/(1024**2), ' Mb'
        end if
        write (stdout, '(1x,a)') '|   However, this will result in more i/o and slow down the calculation      |'
      endif

      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'plot_wannier:', (mem_param + mem_wan)/(1024**2), ' Mb'
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, *) ' '
    endif

!    if(w90_calcs%disentanglement) then
!       write(*,'(a12,f12.4,a)') 'Disentangle',(mem_param+mem_dis)/(1024**2),' Mb'
!    end if
!    write(*,'(a12,f12.4,a)') 'Wannierise ',(mem_wan+mem_param)/(1024**2),' Mb'
!    write(*,'(a12,f12.4,a)') 'Module',(mem_param)/(1024**2),' Mb'

    return
  end subroutine param_memory_estimate

!===========================================================!
  subroutine param_dist(driver, w90_calcs, pp_calc, param_input, param_plot, param_wannierise, &
                        lsitesymmetry, symmetrize_eps, wann_data, param_hamil, kmesh_data, &
                        kmesh_info, k_points, num_kpts, dis_data, fermi_surface_data, fermi, &
                        tran, atoms, num_bands, num_wann, eigval, mp_grid, num_proj, real_lattice, &
                        recip_lattice, eig_found, lhasproj, stdout, seedname, comm)
    !===========================================================!
    !                                                           !
    !! distribute the parameters across processors              !
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp !, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, io_date, io_time
    use w90_comms, only: comms_bcast, w90commtype, mpirank

    implicit none
    !data from parameters module
    type(param_driver_type), intent(inout) :: driver
    type(w90_calculation_type), intent(inout) :: w90_calcs
    type(postproc_type), intent(inout) :: pp_calc
    type(parameter_input_type), intent(inout) :: param_input
    type(param_plot_type), intent(inout) :: param_plot
    type(param_wannierise_type), intent(inout) :: param_wannierise
    ! RS: symmetry-adapted Wannier functions
    logical, intent(inout) :: lsitesymmetry
    real(kind=dp), intent(inout) :: symmetrize_eps
    type(wannier_data_type), intent(inout) :: wann_data
    type(param_hamiltonian_type), intent(inout) :: param_hamil
    type(param_kmesh_type), intent(inout) :: kmesh_data
    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(k_point_type), intent(inout) :: k_points
    integer, intent(inout) :: num_kpts
    type(disentangle_type), intent(inout) :: dis_data
    type(fermi_surface_type), intent(inout) :: fermi_surface_data
    type(fermi_data_type), intent(inout) :: fermi
    type(transport_type), intent(inout) :: tran
    type(atom_data_type), intent(inout) :: atoms
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_wann
    integer, intent(inout) :: stdout
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    integer, intent(inout) :: mp_grid(3)
    integer, intent(inout) :: num_proj
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: recip_lattice(3, 3)
    !type(pw90_calculation_type), intent(inout) :: pw90_calcs
    !type(postw90_oper_type), intent(inout) :: postw90_oper
    !type(postw90_common_type), intent(inout) :: pw90_common
    !type(postw90_spin_type), intent(inout) :: pw90_spin
    !type(postw90_ham_type), intent(inout) :: pw90_ham
    !type(kpath_type), intent(inout) :: kpath
    !type(kslice_type), intent(inout) :: kslice
    !type(dos_plot_type), intent(inout) :: dos_data
    !type(berry_type), intent(inout) :: berry
    !type(spin_hall_type), intent(inout) :: spin_hall
    !type(gyrotropic_type), intent(inout) :: gyrotropic
    !type(geninterp_type), intent(inout) :: geninterp
    !type(boltzwann_type), intent(inout) :: boltz
    logical, intent(inout) :: eig_found
    logical, intent(inout) :: lhasproj
    type(w90commtype), intent(in) :: comm
    character(len=50), intent(in)  :: seedname
    logical :: on_root = .false.
    integer :: ierr
    integer :: iprintroot !JJ

    if (mpirank(comm) == 0) on_root = .true.

    !call comms_bcast(pw90_common%effective_model, 1)
    call comms_bcast(eig_found, 1, stdout, seedname, comm)
    call comms_bcast(driver%postproc_setup, 1, stdout, seedname, comm)
    call comms_bcast(w90_calcs%cp_pp, 1, stdout, seedname, comm)
    !if (.not. pw90_common%effective_model) then
    call comms_bcast(mp_grid(1), 3, stdout, seedname, comm)
    call comms_bcast(num_kpts, 1, stdout, seedname, comm)
    call comms_bcast(num_bands, 1, stdout, seedname, comm)
    !endif
    call comms_bcast(num_wann, 1, stdout, seedname, comm)
    call comms_bcast(param_input%timing_level, 1, stdout, seedname, comm)

    !______________________________________
    !JJ fixme maybe? not so pretty solution to setting iprint to zero on non-root processes
    iprintroot = param_input%iprint
    param_input%iprint = 0
    call comms_bcast(param_input%iprint, 1, stdout, seedname, comm)
    if (on_root) param_input%iprint = iprintroot
    !______________________________________

    !call comms_bcast(energy_unit, 1, stdout, seedname, comm)
    call comms_bcast(param_input%length_unit, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%wvfn_formatted, 1, stdout, seedname, comm)
    !call comms_bcast(postw90_oper%spn_formatted, 1)
    !call comms_bcast(postw90_oper%uHu_formatted, 1)
    !call comms_bcast(berry_uHu_formatted, 1)
    call comms_bcast(param_plot%spin, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%num_dump_cycles, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%num_print_cycles, 1, stdout, seedname, comm)
    call comms_bcast(atoms%num_atoms, 1, stdout, seedname, comm)   ! Ivo: not used in postw90, right?
    call comms_bcast(atoms%num_species, 1, stdout, seedname, comm) ! Ivo: not used in postw90, right?
    call comms_bcast(real_lattice(1, 1), 9, stdout, seedname, comm)
    call comms_bcast(recip_lattice(1, 1), 9, stdout, seedname, comm)
    !call comms_bcast(real_metric(1, 1), 9)
    !call comms_bcast(recip_metric(1, 1), 9)
    !call comms_bcast(cell_volume, 1)
    !call comms_bcast(dos_data%energy_step, 1)
    !call comms_bcast(dos_data%adpt_smr, 1)
    !call comms_bcast(dos_data%smr_index, 1)
    !call comms_bcast(dos_data%kmesh_spacing, 1)
    !call comms_bcast(dos_data%kmesh(1), 3)
    !call comms_bcast(dos_data%adpt_smr_max, 1)
    !call comms_bcast(dos_data%smr_fixed_en_width, 1)
    !call comms_bcast(dos_data%adpt_smr_fac, 1)
    !call comms_bcast(dos_data%num_project, 1)
    call comms_bcast(param_input%num_exclude_bands, 1, stdout, seedname, comm)
    if (param_input%num_exclude_bands > 0) then
      if (.not. on_root) then
        allocate (param_input%exclude_bands(param_input%num_exclude_bands), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating exclude_bands in param_dist', stdout, seedname)
      endif

      call comms_bcast(param_input%exclude_bands(1), param_input%num_exclude_bands, stdout, &
                       seedname, comm)
    end if

    call comms_bcast(param_input%gamma_only, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%win_min, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%win_max, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%froz_min, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%froz_max, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%num_iter, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%mix_ratio, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%conv_tol, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%conv_window, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%spheres_first_wann, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%spheres_num, 1, stdout, seedname, comm)
    if (dis_data%spheres_num > 0) then
      if (.not. on_root) then
        allocate (dis_data%spheres(4, dis_data%spheres_num), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating dis_spheres in param_dist', stdout, seedname)
      endif
      call comms_bcast(dis_data%spheres(1, 1), 4*dis_data%spheres_num, stdout, seedname, comm)
    end if
    call comms_bcast(param_wannierise%num_iter, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%num_cg_steps, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%conv_tol, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%conv_window, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%guiding_centres, 1, stdout, seedname, comm)
    call comms_bcast(w90_calcs%wannier_plot, 1, stdout, seedname, comm) !!BGS!!
    call comms_bcast(param_plot%num_wannier_plot, 1, stdout, seedname, comm)
    if (param_plot%num_wannier_plot > 0) then
      if (.not. on_root) then
        allocate (param_plot%wannier_plot_list(param_plot%num_wannier_plot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating wannier_plot_list in param_dist', stdout, seedname)
      endif
      call comms_bcast(param_plot%wannier_plot_list(1), param_plot%num_wannier_plot, stdout, &
                       seedname, comm)
    end if
    call comms_bcast(param_plot%wannier_plot_supercell(1), 3, stdout, seedname, comm)
    call comms_bcast(param_plot%wannier_plot_format, len(param_plot%wannier_plot_format), stdout, &
                     seedname, comm)
    call comms_bcast(param_plot%wannier_plot_mode, len(param_plot%wannier_plot_mode), stdout, &
                     seedname, comm)
    call comms_bcast(param_plot%wannier_plot_spinor_mode, len(param_plot%wannier_plot_spinor_mode), &
                     stdout, seedname, comm)
    call comms_bcast(param_plot%write_u_matrices, 1, stdout, seedname, comm)
    call comms_bcast(w90_calcs%bands_plot, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%write_bvec, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%bands_num_points, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%bands_plot_format, len(param_plot%bands_plot_format), stdout, &
                     seedname, comm)
    call comms_bcast(param_input%bands_plot_mode, len(param_input%bands_plot_mode), stdout, &
                     seedname, comm)
    call comms_bcast(param_plot%num_bands_project, 1, stdout, seedname, comm)

    if (param_plot%num_bands_project > 0) then
      if (.not. on_root) then
        allocate (param_plot%bands_plot_project(param_plot%num_bands_project), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bands_plot_project in param_dist', stdout, seedname)
      endif
      call comms_bcast(param_plot%bands_plot_project(1), param_plot%num_bands_project, stdout, &
                       seedname, comm)
    end if
    call comms_bcast(param_plot%bands_plot_dim, 1, stdout, seedname, comm)
    call comms_bcast(w90_calcs%write_hr, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%write_rmn, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%write_tb, 1, stdout, seedname, comm)
    call comms_bcast(param_input%hr_cutoff, 1, stdout, seedname, comm)
    call comms_bcast(param_input%dist_cutoff, 1, stdout, seedname, comm)
    call comms_bcast(param_input%dist_cutoff_mode, len(param_input%dist_cutoff_mode), stdout, &
                     seedname, comm)
    call comms_bcast(param_input%dist_cutoff_hc, 1, stdout, seedname, comm)
    !call comms_bcast(one_dim_axis, len(one_dim_axis), stdout, seedname, comm)
    call comms_bcast(param_input%use_ws_distance, 1, stdout, seedname, comm)
    call comms_bcast(param_input%ws_distance_tol, 1, stdout, seedname, comm)
    call comms_bcast(param_input%ws_search_size(1), 3, stdout, seedname, comm)
    call comms_bcast(w90_calcs%fermi_surface_plot, 1, stdout, seedname, comm)
    call comms_bcast(fermi_surface_data%num_points, 1, stdout, seedname, comm)
    call comms_bcast(fermi_surface_data%plot_format, len(fermi_surface_data%plot_format), stdout, &
                     seedname, comm)
    !call comms_bcast(fermi_energy, 1) !! used?

    !call comms_bcast(pw90_calcs%berry, 1)
    !call comms_bcast(berry%task, len(berry%task))
    !call comms_bcast(berry%kmesh_spacing, 1)
    !call comms_bcast(berry%kmesh(1), 3)
    !call comms_bcast(berry%curv_adpt_kmesh, 1)
    !call comms_bcast(berry%curv_adpt_kmesh_thresh, 1)
    !call comms_bcast(berry%curv_unit, len(berry%curv_unit))
!  Stepan Tsirkin
    !call comms_bcast(pw90_calcs%gyrotropic, 1)
    !call comms_bcast(gyrotropic%task, len(gyrotropic%task))
    !call comms_bcast(gyrotropic%kmesh_spacing, 1)
    !call comms_bcast(gyrotropic%kmesh(1), 3)
    !call comms_bcast(gyrotropic%smr_fixed_en_width, 1)
    !call comms_bcast(gyrotropic%smr_index, 1)
    !call comms_bcast(gyrotropic%eigval_max, 1)
    !call comms_bcast(gyrotropic%nfreq, 1)
    !call comms_bcast(gyrotropic%degen_thresh, 1)
    !call comms_bcast(gyrotropic%num_bands, 1)
    !call comms_bcast(gyrotropic%box(1, 1), 9)
    !call comms_bcast(gyrotropic%box_corner(1), 3)
    !call comms_bcast(gyrotropic%smr_max_arg, 1)
    !call comms_bcast(gyrotropic%smr_fixed_en_width, 1)
    !call comms_bcast(gyrotropic%smr_index, 1)

    !call comms_bcast(berry%kubo_adpt_smr, 1)
    !call comms_bcast(berry%kubo_adpt_smr_fac, 1)
    !call comms_bcast(berry%kubo_adpt_smr_max, 1)
    !call comms_bcast(berry%kubo_smr_fixed_en_width, 1)
    !call comms_bcast(berry%kubo_smr_index, 1)
    !call comms_bcast(berry%kubo_eigval_max, 1)
    !call comms_bcast(berry%kubo_nfreq, 1)
    call comms_bcast(fermi%n, 1, stdout, seedname, comm)
    !call comms_bcast(dos_data%energy_min, 1)
    !call comms_bcast(dos_data%energy_max, 1)
    !call comms_bcast(pw90_spin%spin_kmesh_spacing, 1)
    !call comms_bcast(pw90_spin%spin_kmesh(1), 3)
    !call comms_bcast(berry%wanint_kpoint_file, 1)
! Junfeng Qiao
    !call comms_bcast(spin_hall%freq_scan, 1)
    !call comms_bcast(spin_hall%alpha, 1)
    !call comms_bcast(spin_hall%beta, 1)
    !call comms_bcast(spin_hall%gamma, 1)
    !call comms_bcast(spin_hall%bandshift, 1)
    !call comms_bcast(spin_hall%bandshift_firstband, 1)
    !call comms_bcast(spin_hall%bandshift_energyshift, 1)

    call comms_bcast(param_input%devel_flag, len(param_input%devel_flag), stdout, seedname, comm)
    !call comms_bcast(pw90_common%spin_moment, 1)
    !call comms_bcast(pw90_spin%spin_axis_polar, 1)
    !call comms_bcast(pw90_spin%spin_axis_azimuth, 1)
    !call comms_bcast(pw90_common%spin_decomp, 1)
    !call comms_bcast(pw90_ham%use_degen_pert, 1)
    !call comms_bcast(pw90_ham%degen_thr, 1)
    call comms_bcast(param_input%num_valence_bands, 1, stdout, seedname, comm)
    !call comms_bcast(pw90_calcs%dos, 1)
    !call comms_bcast(dos_data%task, len(dos_data%task))
    !call comms_bcast(pw90_calcs%kpath, 1)
    !call comms_bcast(kpath%task, len(kpath%task))
    !call comms_bcast(kpath%bands_colour, len(kpath%bands_colour))
    !call comms_bcast(pw90_calcs%kslice, 1)
    !call comms_bcast(kslice%task, len(kslice%task))
    !call comms_bcast(berry%transl_inv, 1)
    call comms_bcast(param_input%num_elec_per_state, 1, stdout, seedname, comm)
    !call comms_bcast(pw90_common%scissors_shift, 1)
    !

! ----------------------------------------------
    !call comms_bcast(pw90_calcs%geninterp, 1)
    !call comms_bcast(geninterp%alsofirstder, 1)
    !call comms_bcast(geninterp%single_file, 1)
    ! [gp-begin, Apr 12, 2012]
    ! BoltzWann variables
    !call comms_bcast(pw90_calcs%boltzwann, 1)
    !call comms_bcast(boltz%calc_also_dos, 1)
    !call comms_bcast(boltz%dir_num_2d, 1)
    !call comms_bcast(boltz%dos_energy_step, 1)
    !call comms_bcast(boltz%dos_energy_min, 1)
    !call comms_bcast(boltz%dos_energy_max, 1)
    !call comms_bcast(boltz%dos_adpt_smr, 1)
    !call comms_bcast(boltz%dos_smr_fixed_en_width, 1)
    !call comms_bcast(boltz%dos_adpt_smr_fac, 1)
    !call comms_bcast(boltz%dos_adpt_smr_max, 1)
    !call comms_bcast(boltz%mu_min, 1)
    !call comms_bcast(boltz%mu_max, 1)
    !call comms_bcast(boltz%mu_step, 1)
    !call comms_bcast(boltz%temp_min, 1)
    !call comms_bcast(boltz%temp_max, 1)
    !call comms_bcast(boltz%temp_step, 1)
    !call comms_bcast(boltz%kmesh_spacing, 1)
    !call comms_bcast(boltz%kmesh(1), 3)
    !call comms_bcast(boltz%tdf_energy_step, 1)
    !call comms_bcast(boltz%relax_time, 1)
    !call comms_bcast(boltz%TDF_smr_fixed_en_width, 1)
    !call comms_bcast(boltz%TDF_smr_index, 1)
    !call comms_bcast(boltz%dos_smr_index, 1)
    !call comms_bcast(boltz%bandshift, 1)
    !call comms_bcast(boltz%bandshift_firstband, 1)
    !call comms_bcast(boltz%bandshift_energyshift, 1)
    ! [gp-end]

    call comms_bcast(param_input%use_ws_distance, 1, stdout, seedname, comm)
    call comms_bcast(w90_calcs%disentanglement, 1, stdout, seedname, comm)

    call comms_bcast(w90_calcs%transport, 1, stdout, seedname, comm)
    call comms_bcast(tran%easy_fix, 1, stdout, seedname, comm)
    call comms_bcast(tran%mode, len(tran%mode), stdout, seedname, comm)
    call comms_bcast(tran%win_min, 1, stdout, seedname, comm)
    call comms_bcast(tran%win_max, 1, stdout, seedname, comm)
    call comms_bcast(tran%energy_step, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_bb, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_ll, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_rr, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_cc, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_lc, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_cr, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_bandc, 1, stdout, seedname, comm)
    call comms_bcast(tran%write_ht, 1, stdout, seedname, comm)
    call comms_bcast(tran%read_ht, 1, stdout, seedname, comm)
    call comms_bcast(tran%use_same_lead, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_cell_ll, 1, stdout, seedname, comm)
    call comms_bcast(tran%num_cell_rr, 1, stdout, seedname, comm)
    call comms_bcast(tran%group_threshold, 1, stdout, seedname, comm)
    call comms_bcast(param_hamil%translation_centre_frac(1), 3, stdout, seedname, comm)
    call comms_bcast(kmesh_data%num_shells, 1, stdout, seedname, comm)
    call comms_bcast(kmesh_data%skip_B1_tests, 1, stdout, seedname, comm)
    call comms_bcast(driver%explicit_nnkpts, 1, stdout, seedname, comm)

    call comms_bcast(pp_calc%only_A, 1, stdout, seedname, comm)
    call comms_bcast(w90_calcs%use_bloch_phases, 1, stdout, seedname, comm)
    call comms_bcast(driver%restart, len(driver%restart), stdout, seedname, comm)
    call comms_bcast(param_wannierise%write_r2mn, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%num_guide_cycles, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%num_no_guide_iter, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%fixed_step, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%trial_step, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%precond, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%write_proj, 1, stdout, seedname, comm)
    call comms_bcast(param_input%timing_level, 1, stdout, seedname, comm)
    call comms_bcast(param_input%spinors, 1, stdout, seedname, comm)
    call comms_bcast(param_input%num_elec_per_state, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%translate_home_cell, 1, stdout, seedname, comm)
    call comms_bcast(param_input%write_xyz, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%write_hr_diag, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%conv_noise_amp, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%conv_noise_num, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%wannier_plot_radius, 1, stdout, seedname, comm)
    call comms_bcast(param_plot%wannier_plot_scale, 1, stdout, seedname, comm)
    call comms_bcast(kmesh_data%tol, 1, stdout, seedname, comm)
    call comms_bcast(param_input%optimisation, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%write_vdw_data, 1, stdout, seedname, comm)
    call comms_bcast(param_input%lenconfac, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%lfixstep, 1, stdout, seedname, comm)
    call comms_bcast(lsitesymmetry, 1, stdout, seedname, comm)
    call comms_bcast(dis_data%frozen_states, 1, stdout, seedname, comm)
    call comms_bcast(symmetrize_eps, 1, stdout, seedname, comm)

    !vv: Constrained centres
    call comms_bcast(param_wannierise%slwf_num, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%slwf_constrain, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%slwf_lambda, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%selective_loc, 1, stdout, seedname, comm)
    if (param_wannierise%selective_loc .and. param_wannierise%slwf_constrain) then
      if (.not. on_root) then
        !allocate (ccentres_frac(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_frac in param_get_centre_constraints', stdout, seedname)
        allocate (param_wannierise%ccentres_cart(num_wann, 3), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ccentres_cart in param_get_centre_constraints', stdout, seedname)
      endif
      !call comms_bcast(ccentres_frac(1, 1), 3*num_wann, stdout, seedname, comm)
      call comms_bcast(param_wannierise%ccentres_cart(1, 1), 3*num_wann, stdout, seedname, comm)
    end if

    ! vv: automatic projections
    call comms_bcast(kmesh_data%auto_projections, 1, stdout, seedname, comm)

    call comms_bcast(num_proj, 1, stdout, seedname, comm)
    call comms_bcast(lhasproj, 1, stdout, seedname, comm)
    if (lhasproj) then
      if (.not. on_root) then
        allocate (kmesh_data%input_proj_site(3, num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_site in param_dist', stdout, seedname)
        allocate (param_wannierise%proj_site(3, num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_site in param_dist', stdout, seedname)
      endif
      call comms_bcast(kmesh_data%input_proj_site(1, 1), 3*num_proj, stdout, seedname, comm)
      call comms_bcast(param_wannierise%proj_site(1, 1), 3*num_wann, stdout, seedname, comm)
    endif

    ! These variables are different from the ones above in that they are
    ! allocatable, and in param_read they were allocated on the root node only
    !
    if (.not. on_root) then
      allocate (fermi%energy_list(fermi%n), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating fermi_energy_read in postw90_param_dist', stdout, seedname)
      !allocate (berry%kubo_freq_list(berry%kubo_nfreq), stat=ierr)
      !if (ierr /= 0) call io_error( &
      !  'Error allocating kubo_freq_list in postw90_param_dist')
      !allocate (dos_data%project(dos_data%num_project), stat=ierr)
      !if (ierr /= 0) &
      !  call io_error('Error allocating dos_project in postw90_param_dist')
      !if (.not. pw90_common%effective_model) then
      if (eig_found) then
        allocate (eigval(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error allocating eigval in postw90_param_dist', stdout, seedname)
      end if
      allocate (k_points%kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating kpt_latt in postw90_param_dist', stdout, seedname)
      !endif
      !allocate (gyrotropic%band_list(gyrotropic%num_bands), stat=ierr)
      !if (ierr /= 0) call io_error( &
      !  'Error allocating gyrotropic_num_bands in postw90_param_dist')
      !allocate (gyrotropic%freq_list(gyrotropic%nfreq), stat=ierr)
      !if (ierr /= 0) call io_error( &
      !  'Error allocating gyrotropic_freq_list in postw90_param_dist')
    end if

    if (fermi%n > 0) call comms_bcast(fermi%energy_list(1), fermi%n, stdout, seedname, comm)
    !if (berry%kubo_nfreq > 0) call comms_bcast(berry%kubo_freq_list(1), berry%kubo_nfreq)
    !call comms_bcast(gyrotropic%freq_list(1), gyrotropic%nfreq)
    !call comms_bcast(gyrotropic%band_list(1), gyrotropic%num_bands)
    !if (dos_data%num_project > 0) &
    !  call comms_bcast(dos_data%project(1), dos_data%num_project)
    !if (.not. pw90_common%effective_model) then
    if (eig_found) then
      call comms_bcast(eigval(1, 1), num_bands*num_kpts, stdout, seedname, comm)
    end if
    call comms_bcast(k_points%kpt_latt(1, 1), 3*num_kpts, stdout, seedname, comm)
    !endif

    !if (.not. pw90_common%effective_model .and. .not. driver%explicit_nnkpts) then
    if (.not. driver%explicit_nnkpts) then

      call comms_bcast(kmesh_info%nnh, 1, stdout, seedname, comm)
      call comms_bcast(kmesh_info%nntot, 1, stdout, seedname, comm)
      call comms_bcast(kmesh_info%wbtot, 1, stdout, seedname, comm)

      if (.not. on_root) then
        allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nnlist in param_dist', stdout, seedname)
        allocate (kmesh_info%neigh(num_kpts, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating neigh in param_dist', stdout, seedname)
        allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nncell in param_dist', stdout, seedname)
        allocate (kmesh_info%wb(kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating wb in param_dist', stdout, seedname)
        allocate (kmesh_info%bka(3, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bka in param_dist', stdout, seedname)
        allocate (kmesh_info%bk(3, kmesh_info%nntot, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bk in param_dist', stdout, seedname)
      end if

      call comms_bcast(kmesh_info%nnlist(1, 1), num_kpts*kmesh_info%nntot, stdout, seedname, comm)
      call comms_bcast(kmesh_info%neigh(1, 1), num_kpts*kmesh_info%nntot/2, stdout, seedname, comm)
      call comms_bcast(kmesh_info%nncell(1, 1, 1), 3*num_kpts*kmesh_info%nntot, stdout, seedname, &
                       comm)
      call comms_bcast(kmesh_info%wb(1), kmesh_info%nntot, stdout, seedname, comm)
      call comms_bcast(kmesh_info%bka(1, 1), 3*kmesh_info%nntot/2, stdout, seedname, comm)
      call comms_bcast(kmesh_info%bk(1, 1, 1), 3*kmesh_info%nntot*num_kpts, stdout, seedname, comm)

    endif

    call comms_bcast(param_wannierise%omega_total, 1, stdout, seedname, comm)
    call comms_bcast(param_wannierise%omega_tilde, 1, stdout, seedname, comm)
    call comms_bcast(param_input%omega_invariant, 1, stdout, seedname, comm)
    call comms_bcast(param_input%have_disentangled, 1, stdout, seedname, comm)

    if (.not. on_root) then
      allocate (wann_data%centres(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in param_dist', stdout, seedname)
      wann_data%centres = 0.0_dp
      allocate (wann_data%spreads(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating wannier_spreads in param_dist', stdout, seedname)
      wann_data%spreads = 0.0_dp
      if (w90_calcs%disentanglement) then
        allocate (dis_data%ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ndimwin in param_dist', stdout, seedname)
        allocate (dis_data%lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating lwindow in param_dist', stdout, seedname)
      endif
    endif

  end subroutine param_dist

end module wannier_methods
