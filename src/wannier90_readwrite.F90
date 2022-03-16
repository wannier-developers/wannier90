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
!  w90_wannier90_readwrite: input/output routines            !
!     specific to wannier90.x                                !
!                                                            !
!------------------------------------------------------------!

module w90_wannier90_readwrite

  !! Read/write routines specific to wannier90.x data types

  use w90_constants, only: dp
  use w90_types
  use w90_error
  use w90_readwrite
  use w90_wannier90_types

  implicit none

  private

  type w90_extra_io_type
    character(len=20) :: one_dim_axis = 'none'
    !! Constrained centres
    real(kind=dp), allocatable :: ccentres_frac(:, :)
  end type w90_extra_io_type

  public :: w90_extra_io_type
  public :: w90_wannier90_readwrite_dist
  public :: w90_wannier90_readwrite_memory_estimate
  public :: w90_wannier90_readwrite_read
  public :: w90_wannier90_readwrite_w90_dealloc
  public :: w90_wannier90_readwrite_write
  public :: w90_wannier90_readwrite_write_chkpt

contains

  !================================================!
  subroutine w90_wannier90_readwrite_read(atom_data, band_plot, dis_control, dis_spheres, &
                                          dis_manifold, exclude_bands, fermi_energy_list, &
                                          fermi_surface_data, kmesh_input, kmesh_info, kpt_latt, &
                                          output_file, wvfn_read, wann_control, wann_omega, proj, &
                                          proj_input, real_space_ham, select_proj, kpoint_path, &
                                          w90_system, tran, print_output, wannier_data, wann_plot, &
                                          w90_extra_io, ws_region, w90_calculation, eigval, &
                                          real_lattice, bohr, symmetrize_eps, mp_grid, num_bands, &
                                          num_kpts, num_proj, num_wann, optimisation, eig_found, &
                                          calc_only_A, cp_pp, gamma_only, lhasproj, library, &
                                          library_param_read_first_pass, lsitesymmetry, &
                                          use_bloch_phases, seedname, stdout, error, comm)
    !================================================!
    !
    !! Read parameters and calculate derived values
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !
    !================================================

    use w90_constants, only: w90_physical_constants_type
    use w90_utility, only: utility_recip_lattice, utility_inverse_mat

    implicit none

    ! arguments
    type(atom_data_type), intent(inout) :: atom_data
    type(band_plot_type), intent(inout) :: band_plot
    type(dis_control_type), intent(inout) :: dis_control
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(dis_spheres_type), intent(inout) :: dis_spheres
    type(fermi_surface_plot_type), intent(inout) :: fermi_surface_data
    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(kpoint_path_type), intent(inout) :: kpoint_path
    type(output_file_type), intent(inout) :: output_file
    type(print_output_type), intent(inout) :: print_output
    type(proj_input_type), intent(inout) :: proj
    type(proj_input_type), intent(inout) :: proj_input
    type(real_space_ham_type), intent(inout) :: real_space_ham
    type(select_projection_type), intent(inout) :: select_proj
    type(transport_type), intent(inout) :: tran
    type(w90_calculation_type), intent(inout) :: w90_calculation
    type(w90_extra_io_type), intent(inout) :: w90_extra_io
    type(w90_system_type), intent(inout) :: w90_system
    type(wann_control_type), intent(inout) :: wann_control
    type(wannier_data_type), intent(inout) :: wannier_data
    type(wannier_plot_type), intent(inout) :: wann_plot
    type(wann_omega_type), intent(inout) :: wann_omega
    type(ws_region_type), intent(inout) :: ws_region
    type(wvfn_read_type), intent(inout) :: wvfn_read
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer, allocatable, intent(inout) :: exclude_bands(:)
    integer, intent(inout) :: mp_grid(3)
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_kpts
    integer, intent(inout) :: num_proj
    integer, intent(inout) :: num_wann
    integer, intent(inout) :: optimisation
    integer, intent(in) :: stdout

    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    real(kind=dp), allocatable, intent(inout) :: fermi_energy_list(:)
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: symmetrize_eps

    character(len=*), intent(in)  :: seedname

    real(kind=dp) :: recip_lattice(3, 3), volume, inv_lattice(3, 3)
    logical, intent(inout) :: eig_found
    logical, intent(in) :: library
    logical, intent(in) :: library_param_read_first_pass
    !Projections
    logical, intent(out) :: lhasproj
    ! RS: symmetry-adapted Wannier functions
    logical, intent(inout) :: lsitesymmetry
    logical, intent(out) :: use_bloch_phases, cp_pp, calc_only_A
    logical, intent(inout) :: gamma_only

    ! local variables
    integer :: num_exclude_bands
    !! Units for energy
    logical :: found_fermi_energy
    logical :: has_kpath
    logical :: disentanglement
    character(len=20) :: energy_unit

    disentanglement = .false.
    call w90_readwrite_in_file(seedname, error, comm)
    if (allocated(error)) return

    call w90_wannier90_readwrite_read_sym(symmetrize_eps, lsitesymmetry, error, comm)
    if (allocated(error)) return

    call w90_readwrite_read_verbosity(print_output, error, comm)
    if (allocated(error)) return

    call w90_readwrite_read_algorithm_control(optimisation, error, comm)
    if (allocated(error)) return

    call w90_wannier90_readwrite_read_w90_calcs(w90_calculation, error, comm)
    if (allocated(error)) return

    call w90_wannier90_readwrite_read_transport(w90_calculation%transport, tran, w90_calculation%restart, &
                                                error, comm)
    if (allocated(error)) return

    call w90_wannier90_readwrite_read_dist_cutoff(real_space_ham, error, comm)
    if (allocated(error)) return

    if (.not. (w90_calculation%transport .and. tran%read_ht)) then
      call w90_readwrite_read_units(print_output%lenconfac, print_output%length_unit, energy_unit, &
                                    bohr, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_num_wann(num_wann, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_exclude_bands(exclude_bands, num_exclude_bands, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_num_bands(.false., library, num_exclude_bands, num_bands, &
                                        num_wann, library_param_read_first_pass, stdout, error, comm)
      if (allocated(error)) return
      disentanglement = (num_bands > num_wann)

      call w90_readwrite_read_lattice(library, real_lattice, bohr, stdout, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_wannierise(wann_control, num_wann, w90_extra_io%ccentres_frac, &
                                                   stdout, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_mp_grid(.false., library, mp_grid, num_kpts, stdout, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_gamma_only(gamma_only, num_kpts, library, stdout, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_post_proc(cp_pp, calc_only_A, w90_calculation%postproc_setup, &
                                                  error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_restart(w90_calculation, seedname, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_system(library, w90_system, stdout, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_kpath(library, kpoint_path, has_kpath, w90_calculation%bands_plot, &
                                    error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_plot_info(wvfn_read, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_band_plot(band_plot, num_wann, has_kpath, &
                                                  w90_calculation%bands_plot, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_wann_plot(wann_plot, num_wann, &
                                                  w90_calculation%wannier_plot, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_fermi_surface(fermi_surface_data, w90_calculation%fermi_surface_plot, &
                                                      error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_fermi_energy(found_fermi_energy, fermi_energy_list, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_outfiles(output_file, num_kpts, w90_system%num_valence_bands, &
                                                 disentanglement, gamma_only, error, comm)
      if (allocated(error)) return

    endif
    ! BGS tran/plot related stuff...
    call w90_wannier90_readwrite_read_one_dim(w90_calculation, band_plot, real_space_ham, &
                                              w90_extra_io%one_dim_axis, tran%read_ht, error, comm)
    if (allocated(error)) return

    call w90_readwrite_read_ws_data(ws_region, error, comm) !ws_search etc
    if (allocated(error)) return

    if (.not. (w90_calculation%transport .and. tran%read_ht)) then
      call w90_readwrite_read_eigvals(.false., .false., .false., &
                                      w90_calculation%bands_plot .or. w90_calculation%fermi_surface_plot .or. &
                                      output_file%write_hr, disentanglement, eig_found, &
                                      eigval, library, w90_calculation%postproc_setup, num_bands, &
                                      num_kpts, stdout, seedname, error, comm)
      if (allocated(error)) return

      dis_manifold%win_min = -1.0_dp
      dis_manifold%win_max = 0.0_dp
      if (eig_found) dis_manifold%win_min = minval(eigval)
      if (eig_found) dis_manifold%win_max = maxval(eigval)
      call w90_readwrite_read_dis_manifold(eig_found, dis_manifold, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_disentangle(dis_control, dis_spheres, num_bands, num_wann, &
                                                    bohr, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_hamil(real_space_ham, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_bloch_phase(use_bloch_phases, disentanglement, error, comm)
      if (allocated(error)) return

      call w90_readwrite_read_kmesh_data(kmesh_input, error, comm)
      if (allocated(error)) return

      call utility_recip_lattice(real_lattice, recip_lattice, volume, error, comm)
      if (allocated(error)) return

      call utility_inverse_mat(real_lattice, inv_lattice)

      call w90_readwrite_read_kpoints(.false., library, kpt_latt, num_kpts, bohr, stdout, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_explicit_kpts(library, w90_calculation, kmesh_info, &
                                                      num_kpts, bohr, error, comm)
      if (allocated(error)) return

      !call w90_wannier90_readwrite_read_global_kmesh(global_kmesh_set, kmesh_spacing, kmesh, recip_lattice, &
      !                             stdout, seedname)
      call w90_readwrite_read_atoms(library, atom_data, real_lattice, bohr, stdout, error, comm)
      if (allocated(error)) return

      call w90_wannier90_readwrite_read_projections(proj, use_bloch_phases, lhasproj, &
                                                    wann_control%guiding_centres%enable, &
                                                    proj_input, select_proj, num_proj, &
                                                    atom_data, inv_lattice, num_wann, gamma_only, &
                                                    w90_system%spinors, library, bohr, stdout, &
                                                    error, comm)
      if (allocated(error)) return

      if (allocated(proj%site)) then
        if (allocated(wann_control%guiding_centres%centres)) &
          deallocate (wann_control%guiding_centres%centres)
        allocate (wann_control%guiding_centres%centres(3, num_wann))
        wann_control%guiding_centres%centres(:, :) = proj%site(:, :)
      endif
      ! projections needs to be allocated before reading constrained centres
      if (wann_control%constrain%constrain) then
        call w90_wannier90_readwrite_read_constrained_centres(w90_extra_io%ccentres_frac, &
                                                              wann_control, real_lattice, &
                                                              num_wann, library, stdout, &
                                                              error, comm)
        if (allocated(error)) return

      endif
    endif
    call w90_readwrite_clean_infile(stdout, seedname, error, comm)
    if (allocated(error)) return

    if (.not. (w90_calculation%transport .and. tran%read_ht)) then
      ! For aesthetic purposes, convert some things to uppercase
      call w90_readwrite_uppercase(atom_data, kpoint_path, print_output%length_unit)
      ! Initialise
      call w90_readwrite_read_final_alloc(disentanglement, dis_manifold, wannier_data, num_wann, &
                                          num_bands, num_kpts, error, comm)
      if (allocated(error)) return
    endif
  end subroutine w90_wannier90_readwrite_read

  !================================================!
  subroutine w90_wannier90_readwrite_read_sym(symmetrize_eps, lsitesymmetry, error, comm)
    !================================================!
    ! Site symmetry
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    logical, intent(inout) :: lsitesymmetry
    real(kind=dp), intent(inout) :: symmetrize_eps
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    ! default value is lsitesymmetry=.false.
    call w90_readwrite_get_keyword('site_symmetry', found, error, comm, l_value=lsitesymmetry)!YN:
    if (allocated(error)) return

    ! default value is symmetrize_eps=0.001
    call w90_readwrite_get_keyword('symmetrize_eps', found, error, comm, r_value=symmetrize_eps)!YN:
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_sym

  !================================================!
  subroutine w90_wannier90_readwrite_read_w90_calcs(w90_calculation, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(w90_calculation_type), intent(inout) :: w90_calculation
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('transport', found, error, comm, l_value=w90_calculation%transport)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('wannier_plot', found, error, comm, l_value=w90_calculation%wannier_plot)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('bands_plot', found, error, comm, l_value=w90_calculation%bands_plot)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('fermi_surface_plot', found, error, comm, l_value=w90_calculation%fermi_surface_plot)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_w90_calcs

  !================================================!
  subroutine w90_wannier90_readwrite_read_transport(transport, tran, restart, error, comm)
    !================================================!
    ! Transport
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    logical, intent(in) :: transport
    type(transport_type), intent(inout) :: tran
    character(len=*), intent(inout) :: restart
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('tran_read_ht', found, error, comm, l_value=tran%read_ht)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_easy_fix', found, error, comm, l_value=tran%easy_fix)
    if (allocated(error)) return

    if (transport .and. tran%read_ht) restart = ' '

    call w90_readwrite_get_keyword('transport_mode', found, error, comm, c_value=tran%mode)
    if (allocated(error)) return

!    if ( .not.tran_read_ht  .and. (index(transport_mode,'lcr').ne.0) ) then
!       call set_error_input(error, 'Error: transport_mode.eq.lcr not compatible with tran_read_ht.eq.false', comm)
!       return
!    endif

    call w90_readwrite_get_keyword('tran_win_min', found, error, comm, r_value=tran%win_min)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_win_max', found, error, comm, r_value=tran%win_max)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_energy_step', found, error, comm, r_value=tran%energy_step)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_bb', found, error, comm, i_value=tran%num_bb)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_ll', found, error, comm, i_value=tran%num_ll)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_rr', found, error, comm, i_value=tran%num_rr)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_cc', found, error, comm, i_value=tran%num_cc)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_lc', found, error, comm, i_value=tran%num_lc)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_cr', found, error, comm, i_value=tran%num_cr)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_bandc', found, error, comm, i_value=tran%num_bandc)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_write_ht', found, error, comm, l_value=tran%write_ht)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_use_same_lead', found, error, comm, l_value=tran%use_same_lead)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_cell_ll', found, error, comm, i_value=tran%num_cell_ll)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_num_cell_rr', found, error, comm, i_value=tran%num_cell_rr)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('tran_group_threshold', found, error, comm, r_value=tran%group_threshold)
    if (allocated(error)) return

    ! checks
    if (transport) then
      if ((index(tran%mode, 'bulk') .eq. 0) .and. (index(tran%mode, 'lcr') .eq. 0)) then
        call set_error_input(error, 'Error: transport_mode not recognised', comm)
        return
      endif
      if (tran%num_bb < 0) then
        call set_error_input(error, 'Error: tran_num_bb < 0', comm)
        return
      endif
      if (tran%num_ll < 0) then
        call set_error_input(error, 'Error: tran_num_ll < 0', comm)
        return
      endif
      if (tran%num_rr < 0) then
        call set_error_input(error, 'Error: tran_num_rr < 0', comm)
        return
      endif
      if (tran%num_cc < 0) then
        call set_error_input(error, 'Error: tran_num_cc < 0', comm)
        return
      endif
      if (tran%num_lc < 0) then
        call set_error_input(error, 'Error: tran_num_lc < 0', comm)
        return
      endif
      if (tran%num_cr < 0) then
        call set_error_input(error, 'Error: tran_num_cr < 0', comm)
        return
      endif
      if (tran%num_bandc < 0) then
        call set_error_input(error, 'Error: tran_num_bandc < 0', comm)
        return
      endif
      if (tran%num_cell_ll < 0) then
        call set_error_input(error, 'Error: tran_num_cell_ll < 0', comm)
        return
      endif
      if (tran%num_cell_rr < 0) then
        call set_error_input(error, 'Error: tran_num_cell_rr < 0', comm)
        return
      endif
      if (tran%group_threshold < 0.0_dp) then
        call set_error_input(error, 'Error: tran_group_threshold < 0', comm)
        return
      endif
    endif

  end subroutine w90_wannier90_readwrite_read_transport

  !================================================!
  subroutine w90_wannier90_readwrite_read_dist_cutoff(real_space_ham, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(real_space_ham_type), intent(inout) :: real_space_ham
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('dist_cutoff_mode', found, error, comm, &
                                   c_value=real_space_ham%dist_cutoff_mode)
    if (allocated(error)) return
    if ((index(real_space_ham%dist_cutoff_mode, 'three_dim') .eq. 0) &
        .and. (index(real_space_ham%dist_cutoff_mode, 'two_dim') .eq. 0) &
        .and. (index(real_space_ham%dist_cutoff_mode, 'one_dim') .eq. 0)) then
      call set_error_input(error, 'Error: dist_cutoff_mode not recognised', comm)
      return
    endif

    call w90_readwrite_get_keyword('dist_cutoff', found, error, comm, &
                                   r_value=real_space_ham%dist_cutoff)
    if (allocated(error)) return

    real_space_ham%dist_cutoff_hc = real_space_ham%dist_cutoff
    call w90_readwrite_get_keyword('dist_cutoff_hc', found, error, comm, &
                                   r_value=real_space_ham%dist_cutoff_hc)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('hr_cutoff', found, error, comm, r_value=real_space_ham%hr_cutoff)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('bands_plot_dim', found, error, comm, i_value=real_space_ham%system_dim)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_dist_cutoff

  !================================================!
  subroutine w90_wannier90_readwrite_read_wannierise(wann_control, num_wann, ccentres_frac, stdout, error, comm)
    !================================================!
    ! Wannierise
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(wann_control_type), intent(inout) :: wann_control
    integer, intent(in) :: num_wann
    real(kind=dp), allocatable, intent(inout) :: ccentres_frac(:, :)
    integer, intent(in) :: stdout
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer :: ierr
    logical :: found

    call w90_readwrite_get_keyword('num_dump_cycles', found, error, comm, &
                                   i_value=wann_control%num_dump_cycles)
    if (allocated(error)) return

    if (wann_control%num_dump_cycles < 0) then
      call set_error_input(error, 'Error: num_dump_cycles must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('num_print_cycles', found, error, comm, &
                                   i_value=wann_control%num_print_cycles)
    if (allocated(error)) return

    if (wann_control%num_print_cycles < 0) then
      call set_error_input(error, 'Error: num_print_cycles must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('num_iter', found, error, comm, &
                                   i_value=wann_control%num_iter)
    if (allocated(error)) return

    if (wann_control%num_iter < 0) then
      call set_error_input(error, 'Error: num_iter must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('num_cg_steps', found, error, comm, &
                                   i_value=wann_control%num_cg_steps)
    if (allocated(error)) return

    if (wann_control%num_cg_steps < 0) then
      call set_error_input(error, 'Error: num_cg_steps must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('conv_tol', found, error, comm, &
                                   r_value=wann_control%conv_tol)
    if (allocated(error)) return

    if (wann_control%conv_tol < 0.0_dp) then
      call set_error_input(error, 'Error: conv_tol must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('conv_noise_amp', found, error, comm, &
                                   r_value=wann_control%conv_noise_amp)
    if (allocated(error)) return

    ! note that the default here is not to check convergence
    wann_control%conv_window = -1
    if (wann_control%conv_noise_amp > 0.0_dp) wann_control%conv_window = 5
    call w90_readwrite_get_keyword('conv_window', found, error, comm, &
                                   i_value=wann_control%conv_window)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('conv_noise_num', found, error, comm, &
                                   i_value=wann_control%conv_noise_num)
    if (allocated(error)) return

    if (wann_control%conv_noise_num < 0) then
      call set_error_input(error, 'Error: conv_noise_num must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('guiding_centres', found, error, comm, &
                                   l_value=wann_control%guiding_centres%enable)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('num_guide_cycles', found, error, comm, &
                                   i_value=wann_control%guiding_centres%num_guide_cycles)
    if (allocated(error)) return

    if (wann_control%guiding_centres%num_guide_cycles < 0) then
      call set_error_input(error, 'Error: num_guide_cycles must be >= 0', comm)
      return
    endif

    call w90_readwrite_get_keyword('num_no_guide_iter', found, error, comm, &
                                   i_value=wann_control%guiding_centres%num_no_guide_iter)
    if (allocated(error)) return

    if (wann_control%guiding_centres%num_no_guide_iter < 0) then
      call set_error_input(error, 'Error: num_no_guide_iter must be >= 0', comm)
      return
    endif

    call w90_readwrite_get_keyword('fixed_step', found, error, comm, &
                                   r_value=wann_control%fixed_step)
    if (allocated(error)) return

    if (found .and. (wann_control%fixed_step < 0.0_dp)) then
      call set_error_input(error, 'Error: fixed_step must be > 0', comm)
      return
    endif
    if (wann_control%fixed_step > 0.0_dp) wann_control%lfixstep = .true.

    call w90_readwrite_get_keyword('trial_step', found, error, comm, &
                                   r_value=wann_control%trial_step)
    if (allocated(error)) return

    if (found .and. wann_control%lfixstep) then
      call set_error_input(error, 'Error: cannot specify both fixed_step and trial_step', comm)
      return
    endif

    call w90_readwrite_get_keyword('precond', found, error, comm, &
                                   l_value=wann_control%precond)
    if (allocated(error)) return

    wann_control%constrain%slwf_num = num_wann
    call w90_readwrite_get_keyword('slwf_num', found, error, comm, &
                                   i_value=wann_control%constrain%slwf_num)
    if (allocated(error)) return

    if (found) then
      if (wann_control%constrain%slwf_num .gt. num_wann .or. &
          wann_control%constrain%slwf_num .lt. 1) then
        call set_error_input(error, 'Error: slwf_num must be an integer between 1 and num_wann', comm)
        return
      end if
      if (wann_control%constrain%slwf_num .lt. num_wann) &
        wann_control%constrain%selective_loc = .true.
    end if

    call w90_readwrite_get_keyword('slwf_constrain', found, error, comm, &
                                   l_value=wann_control%constrain%constrain)
    if (allocated(error)) return

    if (found .and. wann_control%constrain%constrain) then
      if (wann_control%constrain%selective_loc) then
        allocate (ccentres_frac(num_wann, 3), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating ccentres_frac in w90_readwrite_get_centre_constraints', comm)
          return
        endif
        allocate (wann_control%constrain%centres(num_wann, 3), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating ccentres_cart in w90_readwrite_get_centre_constraints', comm)
          return
        endif
      else
        write (stdout, *) ' No selective localisation requested. Ignoring constraints on centres'
        wann_control%constrain%constrain = .false.
      end if
    end if

    call w90_readwrite_get_keyword('slwf_lambda', found, error, comm, &
                                   r_value=wann_control%constrain%lambda)
    if (allocated(error)) return

    if (found) then
      if (wann_control%constrain%lambda < 0.0_dp) then
        call set_error_input(error, 'Error: slwf_lambda  must be positive.', comm)
        return
      endif
    endif
  end subroutine w90_wannier90_readwrite_read_wannierise

  !================================================!
  subroutine w90_wannier90_readwrite_read_disentangle(dis_control, dis_spheres, num_bands, &
                                                      num_wann, bohr, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(dis_control_type), intent(inout) :: dis_control
    type(dis_spheres_type), intent(inout) :: dis_spheres
    integer, intent(in) :: num_bands, num_wann
    real(kind=dp), intent(in) :: bohr
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer :: nkp, ierr
    logical :: found

    call w90_readwrite_get_keyword('dis_num_iter', found, error, comm, i_value=dis_control%num_iter)
    if (allocated(error)) return
    if (dis_control%num_iter < 0) then
      call set_error_input(error, 'Error: dis_num_iter must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('dis_mix_ratio', found, error, comm, r_value=dis_control%mix_ratio)
    if (allocated(error)) return
    if (dis_control%mix_ratio <= 0.0_dp .or. dis_control%mix_ratio > 1.0_dp) then
      call set_error_input(error, 'Error: dis_mix_ratio must be greater than 0.0 but not greater than 1.0', comm)
      return
    endif

    call w90_readwrite_get_keyword('dis_conv_tol', found, error, comm, r_value=dis_control%conv_tol)
    if (allocated(error)) return
    if (dis_control%conv_tol < 0.0_dp) then
      call set_error_input(error, 'Error: dis_conv_tol must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('dis_conv_window', found, error, comm, i_value=dis_control%conv_window)
    if (allocated(error)) return
    if (dis_control%conv_window < 0) then
      call set_error_input(error, 'Error: dis_conv_window must be positive', comm)
      return
    endif

    ! GS-start
    call w90_readwrite_get_keyword('dis_spheres_first_wann', found, error, comm, i_value=dis_spheres%first_wann)
    if (allocated(error)) return
    if (dis_spheres%first_wann < 1) then
      call set_error_input(error, 'Error: dis_spheres_first_wann must be greater than 0', comm)
      return
    endif
    if (dis_spheres%first_wann > num_bands - num_wann + 1) then
      call set_error_input(error, 'Error: dis_spheres_first_wann is larger than num_bands-num_wann+1', comm)
      return
    endif
    call w90_readwrite_get_keyword('dis_spheres_num', found, error, comm, i_value=dis_spheres%num)
    if (allocated(error)) return
    if (dis_spheres%num < 0) then
      call set_error_input(error, 'Error: dis_spheres_num cannot be negative', comm)
      return
    endif
    if (dis_spheres%num > 0) then
      allocate (dis_spheres%spheres(4, dis_spheres%num), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating dis_spheres in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_keyword_block('dis_spheres', found, dis_spheres%num, 4, &
                                           bohr, error, comm, r_value=dis_spheres%spheres)
      if (allocated(error)) return
      if (.not. found) then
        call set_error_input(error, 'Error: Did not find dis_spheres in the input file', comm)
        return
      endif
      do nkp = 1, dis_spheres%num
        if (dis_spheres%spheres(4, nkp) < 1.0e-15_dp) then
          call set_error_input(error, 'Error: radius for dis_spheres must be > 0', comm)
          return
        endif
      enddo
    endif
    ! GS-end
  end subroutine w90_wannier90_readwrite_read_disentangle

  !================================================!
  subroutine w90_wannier90_readwrite_read_post_proc(cp_pp, pp_only_A, postproc_setup, error, comm)
    !================================================!
    use w90_io, only: post_proc_flag
    use w90_error, only: w90_error_type
    implicit none
    logical, intent(inout) :: cp_pp, pp_only_A
    logical, intent(inout) :: postproc_setup ! alreadt init in w90_calc
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('postproc_setup', found, error, comm, l_value=postproc_setup)
    if (allocated(error)) return
    ! We allow this keyword to be overriden by a command line arg -pp
    if (post_proc_flag) postproc_setup = .true.

    call w90_readwrite_get_keyword('cp_pp', found, error, comm, l_value=cp_pp)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('calc_only_A', found, error, comm, l_value=pp_only_A)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_post_proc

  !================================================!
  subroutine w90_wannier90_readwrite_read_restart(w90_calculation, seedname, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(w90_calculation_type), intent(inout) :: w90_calculation
    character(len=*), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found, chk_found

    call w90_readwrite_get_keyword('restart', found, error, comm, c_value=w90_calculation%restart)
    if (allocated(error)) return
    if (found) then
      if ((w90_calculation%restart .ne. 'default') .and. (w90_calculation%restart .ne. 'wannierise') &
          .and. (w90_calculation%restart .ne. 'plot') .and. (w90_calculation%restart .ne. 'transport')) then
        call set_error_input(error, 'Error in input file: value of restart not recognised', comm)
        return
      else
        inquire (file=trim(seedname)//'.chk', exist=chk_found)
        if (.not. chk_found) then
          call set_error_file(error, 'Error: restart requested but '//trim(seedname)//'.chk file not found', comm)
          return
        endif
      endif
    endif
    !post processing takes priority (user is not warned of this)
    if (w90_calculation%postproc_setup) w90_calculation%restart = ' '
  end subroutine w90_wannier90_readwrite_read_restart

  !================================================!
  subroutine w90_wannier90_readwrite_read_outfiles(output_file, num_kpts, num_valence_bands, disentanglement, &
                                                   gamma_only, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(output_file_type), intent(inout) :: output_file
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_valence_bands
    logical, intent(in) :: disentanglement, gamma_only
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found, hr_plot

    call w90_readwrite_get_keyword('write_xyz', found, error, comm, l_value=output_file%write_xyz)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('write_r2mn', found, error, comm, l_value=output_file%write_r2mn)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('write_proj', found, error, comm, l_value=output_file%write_proj)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('write_hr_diag', found, error, comm, &
                                   l_value=output_file%write_hr_diag)
    if (allocated(error)) return

    hr_plot = .false.
    call w90_readwrite_get_keyword('hr_plot', found, error, comm, l_value=hr_plot)
    if (allocated(error)) return
    if (found) then
      call set_error_input(error, 'Input parameter hr_plot is no longer used. Please use write_hr instead.', comm)
      return
    endif
    call w90_readwrite_get_keyword('write_hr', found, error, comm, l_value=output_file%write_hr)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('write_rmn', found, error, comm, l_value=output_file%write_rmn)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('write_tb', found, error, comm, l_value=output_file%write_tb)
    if (allocated(error)) return

    !%%%%%%%%%%%%%%%%
    !  Other Stuff
    !%%%%%%%%%%%%%%%%

    ! aam: vdW
    call w90_readwrite_get_keyword('write_vdw_data', found, error, comm, &
                                   l_value=output_file%write_vdw_data)
    if (allocated(error)) return
    if (output_file%write_vdw_data) then
      if ((.not. gamma_only) .or. (num_kpts .ne. 1)) then
        call set_error_input(error, 'Error: write_vdw_data may only be used with a single k-point at Gamma', comm)
        return
      endif
    endif
    if (output_file%write_vdw_data .and. disentanglement .and. num_valence_bands <= 0) then
      call set_error_input(error, 'If writing vdw data and disentangling then num_valence_bands must be defined', comm)
      return
    endif

    call w90_readwrite_get_keyword('write_u_matrices', found, error, comm, &
                                   l_value=output_file%write_u_matrices)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('write_bvec', found, error, comm, l_value=output_file%write_bvec)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_outfiles

  !================================================!
  subroutine w90_wannier90_readwrite_read_plot_info(wvfn_read, error, comm)
    !================================================!
    ! Plotting
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(wvfn_read_type), intent(inout) :: wvfn_read
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found
    character(len=6) :: spin_str

    call w90_readwrite_get_keyword('wvfn_formatted', found, error, comm, l_value=wvfn_read%formatted)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('spin', found, error, comm, c_value=spin_str)
    if (allocated(error)) return
    if (found) then
      if (index(spin_str, 'up') > 0) then
        wvfn_read%spin_channel = 1
      elseif (index(spin_str, 'down') > 0) then
        wvfn_read%spin_channel = 2
      else
        call set_error_input(error, 'Error: unrecognised value of spin found: '//trim(spin_str), comm)
        return
      end if
    end if

  end subroutine w90_wannier90_readwrite_read_plot_info

  !================================================!
  subroutine w90_wannier90_readwrite_read_band_plot(band_plot, num_wann, has_kpath, bands_plot, error, comm)
    !================================================!
    ! Plotting
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(band_plot_type), intent(inout) :: band_plot
    integer, intent(in) :: num_wann
    logical, intent(in) :: has_kpath
    logical, intent(in) :: bands_plot
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer :: ierr, num_project
    logical :: found

    call w90_readwrite_get_keyword('bands_plot_format', found, error, comm, c_value=band_plot%format)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('bands_plot_mode', found, error, comm, c_value=band_plot%mode)
    if (allocated(error)) return

    num_project = 0
    call w90_readwrite_get_range_vector('bands_plot_project', found, num_project, .true., error, comm)
    if (allocated(error)) return
    if (found) then
      if (num_project < 1) then
        call set_error_input(error, 'Error: problem reading bands_plot_project', comm)
        return
      endif
      if (allocated(band_plot%project)) deallocate (band_plot%project)
      allocate (band_plot%project(num_project), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating bands_plot_project in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_range_vector('bands_plot_project', found, &
                                          num_project, .false., error, comm, band_plot%project)
      if (allocated(error)) return
      if (any(band_plot%project < 1) .or. any(band_plot%project > num_wann)) then
        call set_error_input(error, 'Error: bands_plot_project asks for a non-valid wannier function to be projected', comm)
        return
      endif
    endif

    if (.not. has_kpath .and. bands_plot) then
      call set_error_input(error, 'A bandstructure plot has been requested but there is no kpoint_path block', comm)
      return
    endif

    ! checks
    if (bands_plot) then
      if ((index(band_plot%format, 'gnu') .eq. 0) .and. &
          (index(band_plot%format, 'xmgr') .eq. 0)) then
        call set_error_input(error, 'Error: bands_plot_format not recognised', comm)
        return
      endif
      if ((index(band_plot%mode, 's-k') .eq. 0) .and. (index(band_plot%mode, 'cut') .eq. 0)) then
        call set_error_input(error, 'Error: bands_plot_mode not recognised', comm)
        return
      endif
    endif

  end subroutine w90_wannier90_readwrite_read_band_plot

  !================================================!
  subroutine w90_wannier90_readwrite_read_wann_plot(wann_plot, num_wann, wannier_plot, error, comm)
    !================================================!
    ! Plotting
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(wannier_plot_type), intent(inout) :: wann_plot
    integer, intent(in) :: num_wann
    logical, intent(in) :: wannier_plot
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer :: i, loop, ierr, wann_plot_num
    logical :: found

    call w90_readwrite_get_vector_length('wannier_plot_supercell', found, i, error, comm)
    if (allocated(error)) return

    if (found) then
      if (i .eq. 1) then
        call w90_readwrite_get_keyword_vector('wannier_plot_supercell', found, 1, error, comm, &
                                              i_value=wann_plot%supercell)
        if (allocated(error)) return

        wann_plot%supercell(2) = wann_plot%supercell(1)
        wann_plot%supercell(3) = wann_plot%supercell(1)
      elseif (i .eq. 3) then
        call w90_readwrite_get_keyword_vector('wannier_plot_supercell', found, 3, error, comm, &
                                              i_value=wann_plot%supercell)
        if (allocated(error)) return

      else
        call set_error_input(error, 'Error: wannier_plot_supercell must be provided as either' &
                             //'one integer or a vector of three integers', comm)
        return
      end if
      if (any(wann_plot%supercell <= 0)) then
        call set_error_input(error, 'Error: wannier_plot_supercell elements must be greater than zero', comm)
        return
      endif
    end if

    call w90_readwrite_get_keyword('wannier_plot_format', found, error, comm, c_value=wann_plot%format)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('wannier_plot_mode', found, error, comm, c_value=wann_plot%mode)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('wannier_plot_spinor_mode', found, error, comm, &
                                   c_value=wann_plot%spinor_mode)
    if (allocated(error)) return
    call w90_readwrite_get_keyword('wannier_plot_spinor_phase', found, error, comm, &
                                   l_value=wann_plot%spinor_phase)
    if (allocated(error)) return

    wann_plot_num = 0
    call w90_readwrite_get_range_vector('wannier_plot_list', found, wann_plot_num, .true., error, comm)
    if (allocated(error)) return
    if (found) then
      if (wann_plot_num < 1) then
        call set_error_input(error, 'Error: problem reading wannier_plot_list', comm)
        return
      endif
      if (allocated(wann_plot%list)) deallocate (wann_plot%list)
      allocate (wann_plot%list(wann_plot_num), stat=ierr)
      if (ierr /= 0) then
        call set_error_input(error, 'Error allocating wannier_plot_list in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_range_vector('wannier_plot_list', found, wann_plot_num, .false., &
                                          error, comm, wann_plot%list)
      if (allocated(error)) return
      if (any(wann_plot%list < 1) .or. any(wann_plot%list > num_wann)) then
        call set_error_input(error, 'Error: wannier_plot_list asks for a non-valid wannier function to be plotted', comm)
        return
      endif
    else
      ! we plot all wannier functions
      wann_plot_num = num_wann
      if (allocated(wann_plot%list)) deallocate (wann_plot%list)
      allocate (wann_plot%list(wann_plot_num), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating wannier_plot_list in w90_wannier90_readwrite_read', comm)
        return
      endif
      do loop = 1, num_wann
        wann_plot%list(loop) = loop
      end do
    end if

    call w90_readwrite_get_keyword('wannier_plot_radius', found, error, comm, r_value=wann_plot%radius)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('wannier_plot_scale', found, error, comm, r_value=wann_plot%scale)
    if (allocated(error)) return

    ! checks
    if (wannier_plot) then
      if ((index(wann_plot%format, 'xcrys') .eq. 0) .and. (index(wann_plot%format, 'cub') .eq. 0)) then
        call set_error_input(error, 'Error: wannier_plot_format not recognised', comm)
        return
      endif
      if ((index(wann_plot%mode, 'crys') .eq. 0) .and. (index(wann_plot%mode, 'mol') .eq. 0)) then
        call set_error_input(error, 'Error: wannier_plot_mode not recognised', comm)
        return
      endif
      if ((index(wann_plot%spinor_mode, 'total') .eq. 0) &
          .and. (index(wann_plot%spinor_mode, 'up') .eq. 0) &
          .and. (index(wann_plot%spinor_mode, 'down') .eq. 0)) then
        call set_error_input(error, 'Error: wannier_plot_spinor_mode not recognised', comm)
        return
      endif
      if (wann_plot%radius < 0.0_dp) then
        call set_error_input(error, 'Error: wannier_plot_radius must be positive', comm)
        return
      endif
      if (wann_plot%scale < 0.0_dp) then
        call set_error_input(error, 'Error: wannier_plot_scale must be positive', comm)
        return
      endif
    endif

  end subroutine w90_wannier90_readwrite_read_wann_plot

  !================================================!
  subroutine w90_wannier90_readwrite_read_fermi_surface(fermi_surface_data, fermi_surface_plot, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(fermi_surface_plot_type), intent(inout) :: fermi_surface_data
    logical, intent(in) :: fermi_surface_plot
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('fermi_surface_num_points', found, error, comm, &
                                   i_value=fermi_surface_data%num_points)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('fermi_surface_plot_format', found, error, comm, &
                                   c_value=fermi_surface_data%plot_format)
    if (allocated(error)) return

    if (fermi_surface_plot) then
      if ((index(fermi_surface_data%plot_format, 'xcrys') .eq. 0)) then
        call set_error_input(error, 'Error: fermi_surface_plot_format not recognised', comm)
        return
      endif
      if (fermi_surface_data%num_points < 0) then
        call set_error_input(error, 'Error: fermi_surface_num_points must be positive', comm)
        return
      endif
    endif
  end subroutine w90_wannier90_readwrite_read_fermi_surface

  !================================================!
  subroutine w90_wannier90_readwrite_read_one_dim(w90_calculation, band_plot, real_space_ham, one_dim_axis, &
                                                  tran_read_ht, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(w90_calculation_type), intent(in) :: w90_calculation
    type(band_plot_type), intent(in) :: band_plot
    type(real_space_ham_type), intent(inout) :: real_space_ham
    character(len=*), intent(inout) :: one_dim_axis
    logical, intent(in) :: tran_read_ht
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('one_dim_axis', found, error, comm, c_value=one_dim_axis)
    if (allocated(error)) return

    if (index(one_dim_axis, 'x') > 0) real_space_ham%one_dim_dir = 1
    if (index(one_dim_axis, 'y') > 0) real_space_ham%one_dim_dir = 2
    if (index(one_dim_axis, 'z') > 0) real_space_ham%one_dim_dir = 3
    if (w90_calculation%transport .and. .not. tran_read_ht .and. &
        (real_space_ham%one_dim_dir .eq. 0)) then
      call set_error_input(error, 'Error: one_dim_axis not recognised', comm)
      return
    endif
    if (w90_calculation%bands_plot .and. (index(band_plot%mode, 'cut') .ne. 0) .and. &
        ((real_space_ham%system_dim .ne. 3) .or. &
         (index(real_space_ham%dist_cutoff_mode, 'three_dim') .eq. 0)) .and. &
        (real_space_ham%one_dim_dir .eq. 0)) then
      call set_error_input(error, 'Error: one_dim_axis not recognised', comm)
      return
    endif

  end subroutine w90_wannier90_readwrite_read_one_dim

  !================================================!
  subroutine w90_wannier90_readwrite_read_hamil(hamiltonian, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    type(real_space_ham_type), intent(inout) :: hamiltonian
    real(kind=dp) :: rv_temp(3)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('translate_home_cell', found, error, comm, &
                                   l_value=hamiltonian%translate_home_cell)
    if (allocated(error)) return

    call w90_readwrite_get_keyword_vector('translation_centre_frac', found, 3, error, comm, &
                                          r_value=rv_temp)
    if (allocated(error)) return

    if (found) then
      hamiltonian%translation_centre_frac = rv_temp
      hamiltonian%automatic_translation = .false.
    endif
  end subroutine w90_wannier90_readwrite_read_hamil

  !================================================!
  subroutine w90_wannier90_readwrite_read_bloch_phase(use_bloch_phases, disentanglement, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    implicit none
    logical, intent(inout) :: use_bloch_phases
    logical, intent(in) :: disentanglement
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('use_bloch_phases', found, error, comm, l_value=use_bloch_phases)
    if (allocated(error)) return

    if (disentanglement .and. use_bloch_phases) then
      call set_error_input(error, 'Error: Cannot use bloch phases for disentanglement', comm)
      return
    endif
  end subroutine w90_wannier90_readwrite_read_bloch_phase

  !================================================!
  subroutine w90_wannier90_readwrite_read_explicit_kpts(library, w90_calculation, kmesh_info, &
                                                        num_kpts, bohr, error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_utility, only: utility_recip_lattice

    implicit none

    ! arguments
    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(w90_calculation_type), intent(in) :: w90_calculation
    integer, intent(in) :: num_kpts
    real(kind=dp), intent(in) :: bohr
    logical, intent(in) :: library
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    ! local variables
    integer, allocatable :: nnkpts_block(:, :)
    integer, allocatable :: nnkpts_idx(:)
    integer :: i, k, ierr, rows
    logical :: found

    ! get the nnkpts block -- this is allowed only in postproc-setup mode
    call w90_readwrite_get_block_length('nnkpts', kmesh_info%explicit_nnkpts, rows, library, error, comm)
    if (allocated(error)) return

    if (kmesh_info%explicit_nnkpts) then
      kmesh_info%nntot = rows/num_kpts
      if (modulo(rows, num_kpts) /= 0) then
        call set_error_input(error, 'The number of rows in nnkpts must be a multiple of num_kpts', comm)
        return
      end if
      if (allocated(nnkpts_block)) deallocate (nnkpts_block)
      allocate (nnkpts_block(5, rows), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nnkpts_block in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_keyword_block('nnkpts', found, rows, 5, bohr, error, comm, &
                                           i_value=nnkpts_block)
      if (allocated(error)) return

      ! check that postproc_setup is true
      if (.not. w90_calculation%postproc_setup) then
        call set_error_input(error, 'Input parameter nnkpts_block is allowed only if postproc_setup = .true.', comm)
        return
      endif

      ! assign the values in nnkpts_block to nnlist and nncell
      ! this keeps track of how many neighbours have been seen for each k-point
      if (allocated(nnkpts_idx)) deallocate (nnkpts_idx)
      allocate (nnkpts_idx(num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nnkpts_idx in w90_wannier90_readwrite_read', comm)
        return
      endif
      nnkpts_idx = 1
      ! allocating "global" nnlist & nncell
      ! These are deallocated in kmesh_dealloc
      if (allocated(kmesh_info%nnlist)) deallocate (kmesh_info%nnlist)
      allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nnlist in w90_wannier90_readwrite_read', comm)
        return
      endif
      if (allocated(kmesh_info%nncell)) deallocate (kmesh_info%nncell)
      allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nncell in w90_wannier90_readwrite_read', comm)
        return
      endif
      do i = 1, num_kpts*kmesh_info%nntot
        k = nnkpts_block(1, i)
        kmesh_info%nnlist(k, nnkpts_idx(k)) = nnkpts_block(2, i)
        kmesh_info%nncell(:, k, nnkpts_idx(k)) = nnkpts_block(3:, i)
        nnkpts_idx(k) = nnkpts_idx(k) + 1
      end do
      ! check that all k-points have the same number of neighbours
      if (any(nnkpts_idx /= (/(kmesh_info%nntot + 1, i=1, num_kpts)/))) then
        call set_error_input(error, 'Inconsistent number of nearest neighbours.', comm)
        return
      end if
      deallocate (nnkpts_idx, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating nnkpts_idx in w90_wannier90_readwrite_read', comm)
        return
      endif
      deallocate (nnkpts_block, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating nnkpts_block in w90_wannier90_readwrite_read', comm)
        return
      endif
    end if

  end subroutine w90_wannier90_readwrite_read_explicit_kpts

  !================================================!
  subroutine w90_wannier90_readwrite_read_projections(proj, use_bloch_phases, lhasproj, &
                                                      guiding_centres, proj_input, select_proj, &
                                                      num_proj, atom_data, recip_lattice, &
                                                      num_wann, gamma_only, spinors, library, &
                                                      bohr, stdout, error, comm)
    !================================================!
    use w90_error, only: w90_error_type

    implicit none

    type(atom_data_type), intent(in) :: atom_data
    type(proj_input_type), intent(inout) :: proj
    type(proj_input_type), intent(inout) :: proj_input
    type(select_projection_type), intent(inout) :: select_proj
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm
    integer, intent(in) :: num_wann
    integer, intent(inout) :: num_proj
    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    logical, intent(in) :: gamma_only
    logical, intent(in) :: spinors
    logical, intent(in) :: use_bloch_phases, guiding_centres, library
    logical, intent(out) :: lhasproj

    integer, intent(in) :: stdout
    integer :: i, j, i_temp, loop, ierr
    logical :: found
    ! projections selection
    integer :: num_select_projections
    integer, allocatable :: select_projections(:)

    ! Projections
    call w90_readwrite_get_keyword('auto_projections', found, error, comm, &
                                   l_value=proj_input%auto_projections)
    if (allocated(error)) return

    call w90_readwrite_get_block_length('projections', found, i_temp, library, error, comm)
    if (allocated(error)) return
    ! check to see that there are no unrecognised keywords
    if (found) then
      if (proj_input%auto_projections) then
        call set_error_input(error, 'Error: Cannot specify both auto_projections and projections block', comm)
        return
      endif
      lhasproj = .true.
      call w90_readwrite_get_projections(num_proj, atom_data, num_wann, proj_input, proj, &
                                         recip_lattice, .true., spinors, bohr, stdout, error, comm)
      if (allocated(error)) return
    else
      if (guiding_centres .and. .not. (gamma_only .and. use_bloch_phases)) then
        call set_error_input(error, 'w90_wannier90_readwrite_read: Guiding centres requested, but no projection block found', comm)
        return
      endif
      lhasproj = .false.
      num_proj = num_wann
    end if

    num_select_projections = 0
    call w90_readwrite_get_range_vector('select_projections', found, num_select_projections, &
                                        .true., error, comm)
    if (allocated(error)) return

    if (found) then
      if (num_select_projections < 1) then
        call set_error_input(error, 'Error: problem reading select_projections', comm)
        return
      endif
      if (allocated(select_projections)) deallocate (select_projections)
      allocate (select_projections(num_select_projections), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating select_projections in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_range_vector('select_projections', found, num_select_projections, &
                                          .false., error, comm, select_projections)
      if (allocated(error)) return

      if (any(select_projections < 1)) then
        call set_error_input(error, 'Error: select_projections must contain positive numbers', comm)
        return
      endif
      if (num_select_projections < num_wann) then
        call set_error_input(error, 'Error: too few projections selected', comm)
        return
      endif
      if (num_select_projections > num_wann) then
        call set_error_input(error, 'Error: too many projections selected', comm)
        return
      endif
      if (.not. lhasproj) then
        call set_error_input(error, 'Error: select_projections cannot be used without defining the projections', comm)
        return
      endif
      if (maxval(select_projections(:)) > num_proj) then
        call set_error_input(error, 'Error: select_projections contains a number greater than num_proj', comm)
        return
      endif
      select_proj%lselproj = .true.
    end if

    if (allocated(select_proj%proj2wann_map)) deallocate (select_proj%proj2wann_map)
    allocate (select_proj%proj2wann_map(num_proj), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating proj2wann_map in w90_wannier90_readwrite_read', comm)
      return
    endif
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
      call w90_readwrite_get_projections(num_proj, atom_data, num_wann, proj_input, proj, &
                                         recip_lattice, .false., spinors, bohr, stdout, error, comm)
      if (allocated(error)) return

      do loop = 1, num_proj
        if (select_proj%proj2wann_map(loop) < 0) cycle
        proj%site(:, select_proj%proj2wann_map(loop)) = proj_input%site(:, loop)
        proj%l(select_proj%proj2wann_map(loop)) = proj_input%l(loop)
        proj%m(select_proj%proj2wann_map(loop)) = proj_input%m(loop)
        proj%z(:, select_proj%proj2wann_map(loop)) = proj_input%z(:, loop)
        proj%x(:, select_proj%proj2wann_map(loop)) = proj_input%x(:, loop)
        proj%radial(select_proj%proj2wann_map(loop)) = proj_input%radial(loop)
        proj%zona(select_proj%proj2wann_map(loop)) = proj_input%zona(loop)
      enddo

      if (spinors) then
        do loop = 1, num_proj
          if (select_proj%proj2wann_map(loop) < 0) cycle
          proj%s(select_proj%proj2wann_map(loop)) = proj_input%s(loop)
          proj%s_qaxis(:, select_proj%proj2wann_map(loop)) = proj_input%s_qaxis(:, loop)
        enddo
      endif
    endif

  end subroutine w90_wannier90_readwrite_read_projections

  !================================================!
  subroutine w90_wannier90_readwrite_read_constrained_centres(ccentres_frac, wann_control, &
                                                              real_lattice, num_wann, library, &
                                                              stdout, error, comm)
    !================================================!
    implicit none
    real(kind=dp), intent(inout) :: ccentres_frac(:, :)
    type(wann_control_type), intent(inout) :: wann_control
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    logical, intent(in) :: library

    integer :: i_temp
    logical :: found

    ! Constrained centres
    call w90_readwrite_get_block_length('slwf_centres', found, i_temp, library, error, comm)
    if (allocated(error)) return

    if (found) then
      if (wann_control%constrain%constrain) then
        ! Allocate array for constrained centres
        call w90_readwrite_get_centre_constraints(ccentres_frac, &
                                                  wann_control%constrain%centres, &
                                                  wann_control%guiding_centres%centres, &
                                                  num_wann, real_lattice, error, comm)
        if (allocated(error)) return

      else
        write (stdout, '(a)') ' slwf_constrain set to false. Ignoring <slwf_centres> block '
      end if
      ! Check that either projections or constrained centres are specified if slwf_constrain=.true.
    elseif (.not. found) then
      if (wann_control%constrain%constrain) then
        if (.not. allocated(wann_control%guiding_centres%centres)) then
          call set_error_input(error, 'Error: slwf_constrain = true, but neither <slwf_centre> block  nor &
               & <projection_block> are specified.', comm)
          return
        else
          ! Allocate array for constrained centres
          call w90_readwrite_get_centre_constraints(ccentres_frac, &
                                                    wann_control%constrain%centres, &
                                                    wann_control%guiding_centres%centres, &
                                                    num_wann, real_lattice, error, comm)
          if (allocated(error)) return

        end if
      end if
    end if
    ! Warning
!jj fixme
    if (wann_control%constrain%constrain .and. allocated(wann_control%guiding_centres%centres) &
        .and. .not. found) &
         & write (stdout, '(a)') ' Warning: No <slwf_centres> block found, but slwf_constrain set to true. &
           & Desired centres for SLWF same as projection centres.'
  end subroutine w90_wannier90_readwrite_read_constrained_centres

  !================================================!
  subroutine w90_wannier90_readwrite_write(atom_data, band_plot, dis_control, dis_spheres, &
                                           fermi_energy_list, fermi_surface_data, kpt_latt, &
                                           output_file, wvfn_read, wann_control, proj, proj_input, &
                                           real_space_ham, select_proj, kpoint_path, tran, &
                                           print_output, wannier_data, wann_plot, w90_extra_io, &
                                           w90_calculation, real_lattice, symmetrize_eps, mp_grid, &
                                           num_bands, num_kpts, num_proj, num_wann, optimisation, &
                                           cp_pp, gamma_only, lsitesymmetry, spinors, &
                                           use_bloch_phases, stdout)
    !================================================!
    !
    !! write wannier90 parameters to stdout
    !
    !================================================
    use w90_utility, only: utility_recip_lattice_base, utility_inverse_mat, utility_cart_to_frac, &
      utility_frac_to_cart

    implicit none

    !passed vaiables
    type(w90_calculation_type), intent(in) :: w90_calculation
    type(output_file_type), intent(in) :: output_file
    type(real_space_ham_type), intent(in) :: real_space_ham
    type(wvfn_read_type), intent(in) :: wvfn_read
    type(print_output_type), intent(in) :: print_output
    type(band_plot_type), intent(in) :: band_plot
    type(wann_control_type), intent(in) :: wann_control
    type(wannier_data_type), intent(in) :: wannier_data
    type(dis_control_type), intent(in) :: dis_control
    type(dis_spheres_type), intent(in) :: dis_spheres
    type(fermi_surface_plot_type), intent(in) :: fermi_surface_data
    type(transport_type), intent(in) :: tran
    type(atom_data_type), intent(in) :: atom_data
    type(select_projection_type), intent(in) :: select_proj
    type(proj_input_type), intent(in) :: proj_input
    type(kpoint_path_type), intent(in) :: kpoint_path
    type(w90_extra_io_type), intent(in) :: w90_extra_io
    type(wannier_plot_type), intent(in) :: wann_plot
    type(proj_input_type), intent(in) :: proj

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_proj
    integer, intent(in) :: num_kpts
    integer, intent(in) :: optimisation

    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: symmetrize_eps
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)

    ! RS: symmetry-adapted Wannier functions
    logical, intent(in) :: lsitesymmetry
    logical, intent(in) :: cp_pp, use_bloch_phases
    logical, intent(in) :: gamma_only
    logical, intent(in) :: spinors

    ! local variables
    real(kind=dp) :: recip_lattice(3, 3), inv_lattice(3, 3), pos_frac(3), kpt_cart(3), volume
    integer :: i, nkp, loop, nat, nsp, bands_num_spec_points
    real(kind=dp) :: cell_volume
    logical :: disentanglement

    disentanglement = (num_bands > num_wann)
    if (w90_calculation%transport .and. tran%read_ht) goto 401

    ! System
    write (stdout, *)
    write (stdout, '(36x,a6)') '------'
    write (stdout, '(36x,a6)') 'SYSTEM'
    write (stdout, '(36x,a6)') '------'
    write (stdout, *)
    if (print_output%lenconfac .eq. 1.0_dp) then
      write (stdout, '(30x,a21)') 'Lattice Vectors (Ang)'
    else
      write (stdout, '(28x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write (stdout, 101) 'a_1', (real_lattice(1, I)*print_output%lenconfac, i=1, 3)
    write (stdout, 101) 'a_2', (real_lattice(2, I)*print_output%lenconfac, i=1, 3)
    write (stdout, 101) 'a_3', (real_lattice(3, I)*print_output%lenconfac, i=1, 3)
    write (stdout, *)
    cell_volume = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - real_lattice(3, 2)*real_lattice(2, 3)) + &
                  real_lattice(1, 2)*(real_lattice(2, 3)*real_lattice(3, 1) - real_lattice(3, 3)*real_lattice(2, 1)) + &
                  real_lattice(1, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - real_lattice(3, 1)*real_lattice(2, 2))
    write (stdout, '(19x,a17,3x,f11.5)', advance='no') &
      'Unit Cell Volume:', cell_volume*print_output%lenconfac**3
    if (print_output%lenconfac .eq. 1.0_dp) then
      write (stdout, '(2x,a7)') '(Ang^3)'
    else
      write (stdout, '(2x,a8)') '(Bohr^3)'
    endif
    write (stdout, *)
    if (print_output%lenconfac .eq. 1.0_dp) then
      write (stdout, '(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
      write (stdout, '(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
    write (stdout, 101) 'b_1', (recip_lattice(1, I)/print_output%lenconfac, i=1, 3)
    write (stdout, 101) 'b_2', (recip_lattice(2, I)/print_output%lenconfac, i=1, 3)
    write (stdout, 101) 'b_3', (recip_lattice(3, I)/print_output%lenconfac, i=1, 3)
    write (stdout, *) ' '
    ! Atoms
    if (atom_data%num_atoms > 0) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (print_output%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
      else
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      call utility_inverse_mat(real_lattice, inv_lattice)
      do nsp = 1, atom_data%num_species
        do nat = 1, atom_data%species_num(nsp)
          call utility_cart_to_frac(atom_data%pos_cart(:, nat, nsp), pos_frac, inv_lattice)
          write (stdout, '(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                 '|', atom_data%symbol(nsp), nat, pos_frac(:),&
  &                 '|', atom_data%pos_cart(:, nat, nsp)*print_output%lenconfac, '|'
        end do
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
    else
      write (stdout, '(25x,a)') 'No atom positions specified'
    end if
    ! Constrained centres
    if (wann_control%constrain%selective_loc .and. &
        wann_control%constrain%constrain) then
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '| Wannier#        Original Centres              Constrained centres          |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do i = 1, wann_control%constrain%slwf_num
        write (stdout, '(1x,a1,2x,i3,2x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                    '|', i, w90_extra_io%ccentres_frac(i, :), '|', wannier_data%centres(:, i), '|'
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
    end if
    ! Projections
    if (print_output%iprint > 1 .and. allocated(proj_input%site)) then
      write (stdout, '(32x,a)') '-----------'
      write (stdout, '(32x,a)') 'PROJECTIONS'
      write (stdout, '(32x,a)') '-----------'
      write (stdout, *) ' '
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a)') '|     Frac. Coord.   l mr  r        z-axis               x-axis          Z/a |'
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, num_proj
        write (stdout, '(1x,a1,3(1x,f5.2),1x,i2,1x,i2,1x,i2,3(1x,f6.3),3(1x,f6.3),2x,f4.1,1x,a1)') &
          '|', proj_input%site(1, nsp), proj_input%site(2, nsp), &
          proj_input%site(3, nsp), proj_input%l(nsp), &
          proj_input%m(nsp), proj_input%radial(nsp), &
          proj_input%z(1, nsp), proj_input%z(2, nsp), &
          proj_input%z(3, nsp), proj_input%x(1, nsp), &
          proj_input%x(2, nsp), proj_input%x(3, nsp), &
          proj_input%zona(nsp), '|'
      end do
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, *) ' '
    end if

    if (print_output%iprint > 1 .and. select_proj%lselproj .and. &
        allocated(wann_control%guiding_centres%centres)) then
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
            &              '|', wann_control%guiding_centres%centres(1, nsp), &
            wann_control%guiding_centres%centres(2, nsp), &
            wann_control%guiding_centres%centres(3, nsp), proj%l(nsp), &
            proj%m(nsp), proj%radial(nsp), &
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
    if (print_output%iprint > 1) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (print_output%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |'
      else
        write (stdout, '(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nkp = 1, num_kpts
        call utility_frac_to_cart(kpt_latt(:, nkp), kpt_cart, recip_lattice)
        write (stdout, '(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|', nkp, kpt_latt(:, nkp), '|', &
          kpt_cart(:)/print_output%lenconfac, '|'
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, *) ' '
    end if
    ! Main
    write (stdout, *) ' '
    write (stdout, '(1x,a78)') '*---------------------------------- MAIN ------------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Wannier Functions               :', num_wann, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Objective Wannier Functions     :', &
      wann_control%constrain%slwf_num, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of input Bloch states              :', num_bands, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :', print_output%iprint, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :', print_output%timing_level, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Optimisation (0=memory, 3=speed)          :', optimisation, '|'
    write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :', trim(print_output%length_unit), '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Post-processing setup (write *.nnkp)      :', &
      w90_calculation%postproc_setup, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using Gamma-only branch of algorithms     :', gamma_only, '|'
    !YN: RS:
    if (lsitesymmetry) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using symmetry-adapted WF mode            :', lsitesymmetry, '|'
      write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Tolerance for symmetry condition on U     :', symmetrize_eps, '|'
    endif

    if (cp_pp .or. print_output%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  CP code post-processing                   :', &
      cp_pp, '|'
    if (w90_calculation%wannier_plot .or. print_output%iprint > 2) then
      if (wvfn_read%formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  Wavefunction (UNK) file-type              :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  Wavefunction (UNK) file-type              :', 'unformatted', '|'
      endif
      if (wvfn_read%spin_channel == 1) then
        write (stdout, '(1x,a46,16x,a2,13x,a1)') '|  Wavefunction spin channel                 :', 'up', '|'
      else
        write (stdout, '(1x,a46,14x,a4,13x,a1)') '|  Wavefunction spin channel                 :', 'down', '|'
      endif
    endif

    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    ! Wannierise
    write (stdout, '(1x,a78)') '*------------------------------- WANNIERISE ---------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Total number of iterations                :', &
      wann_control%num_iter, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of CG steps before reset           :', &
      wann_control%num_cg_steps, '|'
    if (wann_control%lfixstep) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fixed step length for minimisation        :', &
        wann_control%fixed_step, '|'
    else
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Trial step length for line search         :', &
        wann_control%trial_step, '|'
    endif
    write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Convergence tolerence                     :', &
      wann_control%conv_tol, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Convergence window                        :', &
      wann_control%conv_window, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between writing output         :', &
      wann_control%num_print_cycles, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between backing up to disk     :', &
      wann_control%num_dump_cycles, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write r^2_nm to file                      :', &
      output_file%write_r2mn, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write xyz WF centres to file              :', &
      output_file%write_xyz, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write on-site energies <0n|H|0n> to file  :', &
      output_file%write_hr_diag, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use guiding centre to control phases      :', &
      wann_control%guiding_centres%enable, '|'
    write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use phases for initial projections        :', &
      use_bloch_phases, '|'
    if (wann_control%guiding_centres%enable .or. print_output%iprint > 2) then
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations before starting guiding centres:', &
        wann_control%guiding_centres%num_no_guide_iter, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Iterations between using guiding centres  :', &
        wann_control%guiding_centres%num_guide_cycles, '|'
    end if
    if (wann_control%constrain%selective_loc .or. print_output%iprint > 2) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Perform selective localization            :', &
        wann_control%constrain%selective_loc, '|'
    end if
    if (wann_control%constrain%constrain .or. print_output%iprint > 2) then
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Use constrains in selective localization  :', &
        wann_control%constrain%constrain, '|'
      write (stdout, '(1x,a46,8x,E10.3,13x,a1)') '|  Value of the Lagrange multiplier          :',&
           &wann_control%constrain%lambda, '|'
    end if
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    !
    ! Disentanglement
    !
    if (disentanglement .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------------- DISENTANGLE --------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Using band disentanglement                :', &
        disentanglement, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Total number of iterations                :', dis_control%num_iter, '|'
      write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|  Mixing ratio                              :', dis_control%mix_ratio, '|'
      write (stdout, '(1x,a46,8x,ES10.3,13x,a1)') '|  Convergence tolerence                     :', dis_control%conv_tol, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Convergence window                        :', dis_control%conv_window, '|'
      ! GS-start
      if (dis_spheres%num .gt. 0) then
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of spheres in k-space              :', dis_spheres%num, '|'
        do nkp = 1, dis_spheres%num
          write (stdout, '(1x,a13,I4,a2,2x,3F8.3,a15,F8.3,9x,a1)') &
            '|   center n.', nkp, ' :', dis_spheres%spheres(1:3, nkp), ',    radius   =', dis_spheres%spheres(4, nkp), '|'
        enddo
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Index of first Wannier band               :', &
          dis_spheres%first_wann, '|'
      endif
      ! GS-end
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    end if
    !
    ! Plotting
    !
    if (w90_calculation%wannier_plot .or. w90_calculation%bands_plot .or. w90_calculation%fermi_surface_plot &
        .or. output_file%write_hr .or. print_output%iprint > 2) then
      !
      write (stdout, '(1x,a78)') '*-------------------------------- PLOTTING ----------------------------------*'
      !
      if (w90_calculation%wannier_plot .or. print_output%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Wannier functions                :', w90_calculation%wannier_plot, '|'
        write (stdout, '(1x,a46,1x,I5,a1,I5,a1,I5,13x,a1)') &
          '|   Size of supercell for plotting           :', &
          wann_plot%supercell(1), 'x', wann_plot%supercell(2), 'x', wann_plot%supercell(3), '|'

        if (real_space_ham%translate_home_cell) then
          write (stdout, '(1x,a46,10x,L8,13x,a1)') &
            '|  Translating WFs to home cell              :', real_space_ham%translate_home_cell, '|'
        end if

        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting mode (molecule or crystal)      :', &
          trim(wann_plot%mode), '|'
        if (spinors) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting mode for spinor WFs             :', &
            trim(wann_plot%spinor_mode), '|'
          write (stdout, '(1x,a46,10x,L8,13x,a1)') '|   Include phase for spinor WFs             :', &
            wann_plot%spinor_phase, '|'
        end if
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :', &
          trim(wann_plot%format), '|'
        if (index(wann_plot%format, 'cub') > 0 .or. print_output%iprint > 2) then
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Plot radius                              :', &
            wann_plot%radius, '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Plot scale                               :', &
            wann_plot%scale, '|'
        endif
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (w90_calculation%fermi_surface_plot .or. print_output%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Fermi surface                    :', &
          w90_calculation%fermi_surface_plot, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Number of plotting points (along b_1)    :', &
          fermi_surface_data%num_points, '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :', &
          trim(fermi_surface_data%plot_format), '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (w90_calculation%bands_plot .or. print_output%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting interpolated bandstructure       :', w90_calculation%bands_plot, '|'
        bands_num_spec_points = 0
        if (allocated(kpoint_path%labels)) bands_num_spec_points = size(kpoint_path%labels)
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Number of K-path sections                :', &
          bands_num_spec_points/2, '|'
        write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Divisions along first K-path section     :', &
          kpoint_path%num_points_first_segment, '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Output format                            :', &
          trim(band_plot%format), '|'
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Output mode                              :', &
          trim(band_plot%mode), '|'
        if (index(band_plot%mode, 'cut') .ne. 0) then
          write (stdout, '(1x,a46,10x,I8,13x,a1)') '|   Dimension of the system                  :', &
            real_space_ham%system_dim, '|'
          if (real_space_ham%system_dim .eq. 1) &
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System extended in                       :', &
            trim(w90_extra_io%one_dim_axis), '|'
          if (real_space_ham%system_dim .eq. 2) &
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   System confined in                       :', &
            trim(w90_extra_io%one_dim_axis), '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off value                :', &
            real_space_ham%hr_cutoff, '|'
          write (stdout, '(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off distance             :', &
            real_space_ham%dist_cutoff, '|'
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian cut-off distance mode        :', &
            trim(real_space_ham%dist_cutoff_mode), '|'
        endif
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
        write (stdout, '(1x,a78)') '|   K-space path sections:                                                   |'
        if (bands_num_spec_points == 0) then
          write (stdout, '(1x,a78)') '|     None defined                                                           |'
        else
          do loop = 1, bands_num_spec_points, 2
            write (stdout, '(1x,a10,1x,a5,1x,3F7.3,5x,a3,1x,a5,1x,3F7.3,3x,a1)') '|    From:', &
              kpoint_path%labels(loop), (kpoint_path%points(i, loop), i=1, 3), &
              'To:', kpoint_path%labels(loop + 1), (kpoint_path%points(i, loop + 1), i=1, 3), '|'
          end do
        end if
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      end if
      !
      if (output_file%write_hr .or. print_output%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Hamiltonian in WF basis          :', output_file%write_hr, '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      endif
      if (output_file%write_vdw_data .or. print_output%iprint > 2) then
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Writing data for Van der Waals post-proc  :', &
          output_file%write_vdw_data, '|'
        write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      endif
      !
    endif

401 continue
    !
    ! Transport
    !
    if (w90_calculation%transport .or. print_output%iprint > 2) then
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
          trim(w90_extra_io%one_dim_axis), '|'
        !
      end if

      write (stdout, '(1x,a78)') '|   Centre of the unit cell to which WF are translated (fract. coords):      |'
      write (stdout, '(1x,a1,35x,F12.6,a1,F12.6,a1,F12.6,3x,a1)') '|', real_space_ham%translation_centre_frac(1), ',', &
        real_space_ham%translation_centre_frac(2), ',', &
        real_space_ham%translation_centre_frac(3), '|'

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

  end subroutine w90_wannier90_readwrite_write

  !================================================!
  subroutine w90_wannier90_readwrite_w90_dealloc(atom_data, band_plot, dis_spheres, dis_manifold, &
                                                 exclude_bands, kmesh_input, kpt_latt, &
                                                 wann_control, proj, proj_input, select_proj, &
                                                 kpoint_path, wannier_data, wann_plot, &
                                                 w90_extra_io, eigval, error, comm)
    !================================================!
    use w90_error, only: w90_error_type

    implicit none

    ! arguments
    type(band_plot_type), intent(inout) :: band_plot
    type(wann_control_type), intent(inout) :: wann_control
    type(wannier_data_type), intent(inout) :: wannier_data
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(dis_spheres_type), intent(inout) :: dis_spheres
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(atom_data_type), intent(inout) :: atom_data
    type(kpoint_path_type), intent(inout) :: kpoint_path
    type(select_projection_type), intent(inout) :: select_proj
    type(w90_extra_io_type), intent(inout) :: w90_extra_io
    type(wannier_plot_type), intent(inout) :: wann_plot
    type(proj_input_type), intent(inout) :: proj
    type(proj_input_type), intent(inout) :: proj_input
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer, allocatable, intent(inout) :: exclude_bands(:)

    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)

    ! local variables
    integer :: ierr

    call w90_readwrite_dealloc(exclude_bands, wannier_data, proj_input, kmesh_input, kpt_latt, &
                               dis_manifold, atom_data, eigval, kpoint_path, error, comm)
    if (allocated(error)) return

    if (allocated(wann_plot%list)) then
      deallocate (wann_plot%list, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wannier_plot_list in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(band_plot%project)) then
      deallocate (band_plot%project, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating bands_plot_project in w90_readwrite_dealloc', comm)
        return
      endif
    endif
    if (allocated(w90_extra_io%ccentres_frac)) then
      deallocate (w90_extra_io%ccentres_frac, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating ccentres_frac in w90_wannier90_readwrite_w90_dealloc', comm)
        return
      endif
    endif
    if (allocated(wann_control%guiding_centres%centres)) then
      deallocate (wann_control%guiding_centres%centres, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wannier proj_site in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(wann_control%constrain%centres)) then
      deallocate (wann_control%constrain%centres, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating ccentres_cart in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%l)) then
      deallocate (proj%l, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_l in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%site)) then
      deallocate (proj%site, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_site in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%m)) then
      deallocate (proj%m, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_m in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%s)) then
      deallocate (proj%s, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_s in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%s_qaxis)) then
      deallocate (proj%s_qaxis, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_s_qaxis in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%z)) then
      deallocate (proj%z, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_z in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%x)) then
      deallocate (proj%x, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_x in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%radial)) then
      deallocate (proj%radial, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_radial in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(proj%zona)) then
      deallocate (proj%zona, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating proj_zona in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(dis_spheres%spheres)) then
      deallocate (dis_spheres%spheres, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating dis_spheres in w90_readwrite_dealloc', comm)
        return
      endif
    endif
    if (allocated(select_proj%proj2wann_map)) then
      deallocate (select_proj%proj2wann_map, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating select_projections in w90_readwrite_dealloc', comm)
        return
      endif
    endif
  end subroutine w90_wannier90_readwrite_w90_dealloc

  !================================================!
  subroutine w90_wannier90_readwrite_write_chkpt(chkpt, exclude_bands, wannier_data, kmesh_info, &
                                                 kpt_latt, num_kpts, dis_manifold, num_bands, &
                                                 num_wann, u_matrix, u_matrix_opt, m_matrix, &
                                                 mp_grid, real_lattice, omega_invariant, &
                                                 have_disentangled, stdout, seedname)
    !================================================!
    !! Write checkpoint file
    !! IMPORTANT! If you change the chkpt format, adapt
    !! accordingly also the w90chk2chk.x utility!
    !! Also, note that this routine writes the u_matrix and the m_matrix - in parallel
    !! mode these are however stored in distributed form in, e.g., u_matrix_loc only, so
    !! if you are changing the u_matrix, remember to gather it from u_matrix_loc first!
    !================================================!

    use w90_io, only: io_file_unit, io_date
    use w90_utility, only: utility_recip_lattice_base

    implicit none

    ! arguments
    type(wannier_data_type), intent(in) :: wannier_data
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(dis_manifold_type), intent(in) :: dis_manifold

    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)

    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: omega_invariant
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    integer, allocatable, intent(in) :: exclude_bands(:)
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    character(len=*), intent(in)  :: seedname
    character(len=*), intent(in) :: chkpt

    logical, intent(in) :: have_disentangled

    ! local variables
    integer :: chk_unit, nkp, i, j, k, l, num_exclude_bands
    real(kind=dp) :: recip_lattice(3, 3), volume
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
    if (allocated(exclude_bands)) then
      num_exclude_bands = size(exclude_bands)
    else
      num_exclude_bands = 0
    endif
    write (chk_unit) num_exclude_bands         ! Number of excluded bands
    write (chk_unit) (exclude_bands(i), i=1, num_exclude_bands) ! Excluded bands
    write (chk_unit) ((real_lattice(i, j), i=1, 3), j=1, 3)        ! Real lattice
    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
    write (chk_unit) ((recip_lattice(i, j), i=1, 3), j=1, 3)       ! Reciprocal lattice
    write (chk_unit) num_kpts                                 ! Number of k-points
    write (chk_unit) (mp_grid(i), i=1, 3)                       ! M-P grid
    write (chk_unit) ((kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts) ! K-points
    write (chk_unit) kmesh_info%nntot                  ! Number of nearest k-point neighbours
    write (chk_unit) num_wann               ! Number of wannier functions
    chkpt1 = adjustl(trim(chkpt))
    write (chk_unit) chkpt1                 ! Position of checkpoint
    write (chk_unit) have_disentangled      ! Whether a disentanglement has been performed
    if (have_disentangled) then
      write (chk_unit) omega_invariant     ! Omega invariant
      ! lwindow, ndimwin and U_matrix_opt
      write (chk_unit) ((dis_manifold%lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
      write (chk_unit) (dis_manifold%ndimwin(nkp), nkp=1, num_kpts)
      write (chk_unit) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
    endif
    write (chk_unit) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)               ! U_matrix
    write (chk_unit) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, kmesh_info%nntot), l=1, num_kpts) ! M_matrix
    write (chk_unit) ((wannier_data%centres(i, j), i=1, 3), j=1, num_wann)
    write (chk_unit) (wannier_data%spreads(i), i=1, num_wann)
    close (chk_unit)

    write (stdout, '(a/)') ' done'

    return

  end subroutine w90_wannier90_readwrite_write_chkpt

  !================================================!
  subroutine w90_wannier90_readwrite_memory_estimate(atom_data, kmesh_info, wann_control, proj_input, print_output, &
                                                     num_bands, num_kpts, num_proj, num_wann, optimisation, &
                                                     gamma_only, stdout)
    !================================================!
    !
    !! Estimate how much memory we will allocate
    !
    !================================================!

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(wann_control_type), intent(in) :: wann_control
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(proj_input_type), intent(in) :: proj_input
    type(atom_data_type), intent(in) :: atom_data

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: num_proj
    integer, intent(in) :: num_kpts
    integer, intent(in) :: optimisation
    logical, intent(in) :: gamma_only

    ! local variables
    real(kind=dp), parameter :: size_log = 1.0_dp
    real(kind=dp), parameter :: size_int = 4.0_dp
    real(kind=dp), parameter :: size_real = 8.0_dp
    real(kind=dp), parameter :: size_cmplx = 16.0_dp
    real(kind=dp) :: mem_wan, mem_wan1, mem_param, mem_dis, mem_dis2, mem_dis1
    real(kind=dp) :: mem_bw
    logical :: disentanglement

    disentanglement = (num_bands > num_wann)
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
      mem_param = mem_param + num_wann*num_wann*kmesh_info%nntot*num_kpts*size_cmplx       !m_matrix

    if (disentanglement) then
      mem_param = mem_param + num_bands*num_wann*num_kpts*size_cmplx             ! u_matrix_opt
    endif

    if (allocated(atom_data%species_num)) then
      mem_param = mem_param + (atom_data%num_species)*size_int                               !atoms_species_num
      mem_param = mem_param + (atom_data%num_species)*size_real                              !atoms_label
      mem_param = mem_param + (atom_data%num_species)*size_real                              !atoms_symbol
      !mem_param = mem_param + (3*maxval(atom_data%species_num)*atom_data%num_species)*size_real  !atoms_pos_frac
      mem_param = mem_param + (3*maxval(atom_data%species_num)*atom_data%num_species)*size_real  !atoms_pos_cart
    endif

    if (allocated(proj_input%site)) then
      mem_param = mem_param + (3*num_proj)*size_real              !input_proj_site
      mem_param = mem_param + (num_proj)*size_int                !input_proj_l
      mem_param = mem_param + (num_proj)*size_int                 !input_proj_m
      mem_param = mem_param + (3*num_proj)*size_real             !input_proj_z
      mem_param = mem_param + (3*num_proj)*size_real             !input_proj_x
      mem_param = mem_param + (num_proj)*size_real                !input_proj_radial
      mem_param = mem_param + (num_proj)*size_real                !input_proj_zona
    endif

    if (allocated(wann_control%guiding_centres%centres)) then
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
    !mem_param = mem_param + 3*num_kpts*size_real                     !kpt_cart
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
      mem_dis2 = mem_dis2 + num_wann*num_wann*kmesh_info%nntot*num_kpts*size_cmplx!m_matrix

      if (optimisation <= 0) then
        mem_dis = mem_dis + mem_dis1
      else
        mem_dis = mem_dis + max(mem_dis1, mem_dis2)
      endif

      mem_dis = mem_dis + num_bands*num_bands*kmesh_info%nntot*num_kpts*size_cmplx      ! m_matrix_orig
      mem_dis = mem_dis + num_bands*num_wann*num_kpts*size_cmplx             ! a_matrix

    endif

    !Wannierise

    mem_wan1 = mem_wan1 + (num_wann*num_wann*kmesh_info%nntot*num_kpts)*size_cmplx     !  'm0'
    if (optimisation > 0) then
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
    if (gamma_only) then
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

    if (disentanglement) &
      mem_wan = mem_wan + num_wann*num_wann*kmesh_info%nntot*num_kpts*size_cmplx       !m_matrix

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a)') '*============================================================================*'
      write (stdout, '(1x,a)') '|                              MEMORY ESTIMATE                               |'
      write (stdout, '(1x,a)') '|         Maximum RAM allocated during each phase of the calculation         |'
      write (stdout, '(1x,a)') '*============================================================================*'
      if (disentanglement) &
        write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Disentanglement:', (mem_param + mem_dis)/(1024**2), ' Mb'
      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:', (mem_param + mem_wan)/(1024**2), ' Mb'
      if (optimisation > 0 .and. print_output%iprint > 1) then
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

      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'plot_wannier:', (mem_param + mem_wan)/(1024**2), ' Mb'
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      write (stdout, *) ' '
    endif

!    if(w90_calculation%disentanglement) then
!       write(*,'(a12,f12.4,a)') 'Disentangle',(mem_param+mem_dis)/(1024**2),' Mb'
!    end if
!    write(*,'(a12,f12.4,a)') 'Wannierise ',(mem_wan+mem_param)/(1024**2),' Mb'
!    write(*,'(a12,f12.4,a)') 'Module',(mem_param)/(1024**2),' Mb'

    return
  end subroutine w90_wannier90_readwrite_memory_estimate

  !================================================!
  subroutine w90_wannier90_readwrite_dist(atom_data, band_plot, dis_control, dis_spheres, &
                                          dis_manifold, exclude_bands, fermi_energy_list, &
                                          fermi_surface_data, kmesh_input, kmesh_info, kpt_latt, &
                                          output_file, wvfn_read, wann_control, wann_omega, &
                                          proj_input, real_space_ham, w90_system, tran, &
                                          print_output, wannier_data, wann_plot, ws_region, &
                                          w90_calculation, eigval, real_lattice, symmetrize_eps, &
                                          mp_grid, first_segment, num_bands, num_kpts, num_proj, &
                                          num_wann, optimisation, eig_found, cp_pp, gamma_only, &
                                          have_disentangled, lhasproj, lsitesymmetry, &
                                          use_bloch_phases, error, comm)
    !================================================!
    !
    !! Distribute the various parameters across processors
    !
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_file_unit, io_date, io_time
    use w90_comms, only: comms_bcast, w90comm_type, mpirank

    implicit none

    ! arguments
    type(atom_data_type), intent(inout) :: atom_data
    type(band_plot_type), intent(inout) :: band_plot
    type(dis_control_type), intent(inout) :: dis_control
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(dis_spheres_type), intent(inout) :: dis_spheres
    type(fermi_surface_plot_type), intent(inout) :: fermi_surface_data
    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(output_file_type), intent(inout) :: output_file
    type(print_output_type), intent(inout) :: print_output
    type(proj_input_type), intent(inout) :: proj_input
    type(real_space_ham_type), intent(inout) :: real_space_ham
    type(transport_type), intent(inout) :: tran
    type(w90_calculation_type), intent(inout) :: w90_calculation
    type(w90comm_type), intent(in) :: comm
    type(w90_system_type), intent(inout) :: w90_system
    type(wann_control_type), intent(inout) :: wann_control
    type(wannier_data_type), intent(inout) :: wannier_data
    type(wannier_plot_type), intent(inout) :: wann_plot
    type(wann_omega_type), intent(inout) :: wann_omega
    type(ws_region_type), intent(inout) :: ws_region
    type(wvfn_read_type), intent(inout) :: wvfn_read
    type(w90_error_type), allocatable, intent(out) :: error

    integer, allocatable, intent(inout) :: exclude_bands(:)
    integer, intent(inout) :: first_segment
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_wann
    integer, intent(inout) :: mp_grid(3)
    integer, intent(inout) :: num_proj
    integer, intent(inout) :: num_kpts
    integer, intent(inout) :: optimisation

    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    real(kind=dp), allocatable, intent(inout) :: fermi_energy_list(:)
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: symmetrize_eps

    logical, intent(inout) :: eig_found
    logical, intent(inout) :: cp_pp
    logical, intent(inout) :: gamma_only
    logical, intent(inout) :: have_disentangled
    logical, intent(inout) :: lhasproj
    logical, intent(inout) :: lsitesymmetry ! RS: symmetry-adapted Wannier functions
    logical, intent(inout) :: use_bloch_phases

    ! local variables
    logical :: on_root = .false.
    integer :: ierr
    integer :: iprintroot
    integer :: num_project, wann_plot_num, num_exclude_bands, fermi_n
    logical :: disentanglement

    !character(len=128) :: errmesg

    if (mpirank(comm) == 0) on_root = .true.

    call comms_bcast(eig_found, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(w90_calculation%postproc_setup, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(cp_pp, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(mp_grid(1), 3, error, comm)
    if (allocated(error)) return

    call comms_bcast(num_kpts, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(num_bands, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(num_wann, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(print_output%timing_level, 1, error, comm)
    if (allocated(error)) return

    disentanglement = (num_bands > num_wann)

    !______________________________________
    !JJ fixme maybe? not so pretty solution to setting iprint to zero on non-root processes
    iprintroot = print_output%iprint
    print_output%iprint = 0
    call comms_bcast(print_output%iprint, 1, error, comm)
    if (allocated(error)) return
    if (on_root) print_output%iprint = iprintroot
    !______________________________________

    !call comms_bcast(energy_unit, 1, error, comm)
    call comms_bcast(print_output%length_unit, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wvfn_read%formatted, 1, error, comm)
    if (allocated(error)) return

    !call comms_bcast(postw90_oper%spn_formatted, 1)
    !if (allocated(error)) return

    !call comms_bcast(postw90_oper%uHu_formatted, 1)
    !if (allocated(error)) return

    !call comms_bcast(berry_uHu_formatted, 1)
    !if (allocated(error)) return

    call comms_bcast(wvfn_read%spin_channel, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%num_dump_cycles, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%num_print_cycles, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(atom_data%num_atoms, 1, error, comm)   ! Ivo: not used in postw90, right?
    if (allocated(error)) return

    call comms_bcast(atom_data%num_species, 1, error, comm) ! Ivo: not used in postw90, right?
    if (allocated(error)) return

    call comms_bcast(real_lattice(1, 1), 9, error, comm)
    if (allocated(error)) return

    num_exclude_bands = 0
    if (on_root) then
      if (allocated(exclude_bands)) num_exclude_bands = size(exclude_bands)
    endif
    call comms_bcast(num_exclude_bands, 1, error, comm)
    if (allocated(error)) return

    if (num_exclude_bands > 0) then
      if (.not. on_root) then
        allocate (exclude_bands(num_exclude_bands), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating exclude_bands in w90_wannier90_readwrite_dist', comm)
          return
        endif
      endif

      call comms_bcast(exclude_bands(1), num_exclude_bands, error, comm)
      if (allocated(error)) return
    end if

    call comms_bcast(gamma_only, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_manifold%win_min, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_manifold%win_max, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_manifold%froz_min, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_manifold%froz_max, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_control%num_iter, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_control%mix_ratio, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_control%conv_tol, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_control%conv_window, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_spheres%first_wann, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_spheres%num, 1, error, comm)
    if (allocated(error)) return

    if (dis_spheres%num > 0) then
      if (.not. on_root) then
        allocate (dis_spheres%spheres(4, dis_spheres%num), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating dis_spheres in w90_wannier90_readwrite_dist', comm)
          return
        endif
      endif
      call comms_bcast(dis_spheres%spheres(1, 1), 4*dis_spheres%num, error, comm)
      if (allocated(error)) return

    end if
    call comms_bcast(wann_control%num_iter, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%num_cg_steps, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%conv_tol, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%conv_window, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%guiding_centres%enable, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(w90_calculation%wannier_plot, 1, error, comm)
    if (allocated(error)) return

    wann_plot_num = 0
    if (on_root) then
      if (allocated(wann_plot%list)) wann_plot_num = size(wann_plot%list)
    endif
    call comms_bcast(wann_plot_num, 1, error, comm)
    if (allocated(error)) return

    if (wann_plot_num > 0) then
      if (.not. on_root) then
        allocate (wann_plot%list(wann_plot_num), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating wannier_plot_list in w90_wannier90_readwrite_dist', comm)
          return
        endif
      endif
      call comms_bcast(wann_plot%list(1), wann_plot_num, error, comm)
      if (allocated(error)) return

    end if
    call comms_bcast(wann_plot%supercell(1), 3, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_plot%format, len(wann_plot%format), error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_plot%mode, len(wann_plot%mode), error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_plot%spinor_mode, len(wann_plot%spinor_mode), error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_u_matrices, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(w90_calculation%bands_plot, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_bvec, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(first_segment, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(band_plot%format, len(band_plot%format), error, comm)
    if (allocated(error)) return

    call comms_bcast(band_plot%mode, len(band_plot%mode), error, comm)
    if (allocated(error)) return

    num_project = 0
    if (on_root) then
      if (allocated(band_plot%project)) num_project = size(band_plot%project)
    endif
    call comms_bcast(num_project, 1, error, comm)
    if (allocated(error)) return

    if (num_project > 0) then
      if (.not. on_root) then
        allocate (band_plot%project(num_project), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating bands_plot_project in w90_wannier90_readwrite_dist', comm)
          return
        endif
      endif
      call comms_bcast(band_plot%project(1), num_project, error, comm)
      if (allocated(error)) return

    end if

    call comms_bcast(real_space_ham%system_dim, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_hr, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_rmn, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_tb, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(real_space_ham%hr_cutoff, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(real_space_ham%dist_cutoff, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(real_space_ham%dist_cutoff_mode, len(real_space_ham%dist_cutoff_mode), error, comm)
    if (allocated(error)) return

    call comms_bcast(real_space_ham%dist_cutoff_hc, 1, error, comm)
    if (allocated(error)) return

    !call comms_bcast(one_dim_axis, len(one_dim_axis), error, comm)
    !if (allocated(error)) return

    call comms_bcast(ws_region%use_ws_distance, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(ws_region%ws_distance_tol, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(ws_region%ws_search_size(1), 3, error, comm)
    if (allocated(error)) return

    call comms_bcast(w90_calculation%fermi_surface_plot, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(fermi_surface_data%num_points, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(fermi_surface_data%plot_format, len(fermi_surface_data%plot_format), error, comm)
    if (allocated(error)) return

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
    fermi_n = 0
    if (on_root) then
      if (allocated(fermi_energy_list)) fermi_n = size(fermi_energy_list)
    endif
    call comms_bcast(fermi_n, 1, error, comm)
    if (allocated(error)) return

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

    !call comms_bcast(print_output%devel_flag, len(print_output%devel_flag), error, comm)
    !call comms_bcast(pw90_common%spin_moment, 1)
    !call comms_bcast(pw90_spin%spin_axis_polar, 1)
    !call comms_bcast(pw90_spin%spin_axis_azimuth, 1)
    !call comms_bcast(pw90_common%spin_decomp, 1)
    !call comms_bcast(pw90_ham%use_degen_pert, 1)
    !call comms_bcast(pw90_ham%degen_thr, 1)
    call comms_bcast(w90_system%num_valence_bands, 1, error, comm)
    if (allocated(error)) return

    !call comms_bcast(pw90_calcs%dos, 1)
    !call comms_bcast(dos_data%task, len(dos_data%task))
    !call comms_bcast(pw90_calcs%kpath, 1)
    !call comms_bcast(kpath%task, len(kpath%task))
    !call comms_bcast(kpath%bands_colour, len(kpath%bands_colour))
    !call comms_bcast(pw90_calcs%kslice, 1)
    !call comms_bcast(kslice%task, len(kslice%task))
    !call comms_bcast(berry%transl_inv, 1)
    call comms_bcast(w90_system%num_elec_per_state, 1, error, comm)
    if (allocated(error)) return
    !call comms_bcast(pw90_common%scissors_shift, 1)

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

    call comms_bcast(ws_region%use_ws_distance, 1, error, comm)
    if (allocated(error)) return

    !call comms_bcast(w90_calculation%disentanglement, 1, error, comm)
    !if (allocated(error)) return

    call comms_bcast(w90_calculation%transport, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%easy_fix, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%mode, len(tran%mode), error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%win_min, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%win_max, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%energy_step, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_bb, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_ll, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_rr, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_cc, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_lc, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_cr, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_bandc, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%write_ht, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%read_ht, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%use_same_lead, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_cell_ll, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%num_cell_rr, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(tran%group_threshold, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(real_space_ham%translation_centre_frac(1), 3, error, comm)
    if (allocated(error)) return

    call comms_bcast(kmesh_input%num_shells, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(kmesh_input%skip_B1_tests, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(kmesh_info%explicit_nnkpts, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(use_bloch_phases, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(w90_calculation%restart, len(w90_calculation%restart), error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_r2mn, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%guiding_centres%num_guide_cycles, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%guiding_centres%num_no_guide_iter, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%fixed_step, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%trial_step, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%precond, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_proj, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(print_output%timing_level, 1, error, comm)
    if (allocated(error)) return

    ! test mpi error handling using "unlucky" input token
    if (print_output%timing_level < 0 .and. mpirank(comm) == abs(print_output%timing_level)) then
      call set_error_input(error, 'received unlucky_rank', comm)
      return
    endif

    call comms_bcast(w90_system%spinors, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(w90_system%num_elec_per_state, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(real_space_ham%translate_home_cell, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_xyz, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_hr_diag, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%conv_noise_amp, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%conv_noise_num, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_plot%radius, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_plot%scale, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(kmesh_input%tol, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(optimisation, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(output_file%write_vdw_data, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(print_output%lenconfac, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%lfixstep, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(lsitesymmetry, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(dis_manifold%frozen_states, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(symmetrize_eps, 1, error, comm)
    if (allocated(error)) return

    !vv: Constrained centres
    call comms_bcast(wann_control%constrain%slwf_num, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%constrain%constrain, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%constrain%lambda, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_control%constrain%selective_loc, 1, error, comm)
    if (allocated(error)) return

    if (wann_control%constrain%selective_loc .and. wann_control%constrain%constrain) then
      if (.not. on_root) then
        !allocate (ccentres_frac(num_wann, 3), stat=ierr)
        !if (ierr /= 0) then
        !  call set_error_alloc(error, 'Error allocating ccentres_frac in w90_readwrite_get_centre_constraints', comm)
        !  return
        !endif
        allocate (wann_control%constrain%centres(num_wann, 3), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating ccentres_cart in w90_readwrite_get_centre_constraints', comm)
          return
        endif
      endif
      call comms_bcast(wann_control%constrain%centres(1, 1), 3*num_wann, error, comm)
      if (allocated(error)) return
    end if

    ! vv: automatic projections
    call comms_bcast(proj_input%auto_projections, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(num_proj, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(lhasproj, 1, error, comm)
    if (allocated(error)) return

    if (lhasproj) then
      if (.not. on_root) then
        allocate (proj_input%site(3, num_proj), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating input_proj_site in w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (wann_control%guiding_centres%centres(3, num_wann), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating proj_site in w90_wannier90_readwrite_dist', comm)
          return
        endif
      endif
      call comms_bcast(proj_input%site(1, 1), 3*num_proj, error, comm)
      if (allocated(error)) return

      call comms_bcast(wann_control%guiding_centres%centres(1, 1), 3*num_wann, error, comm)
      if (allocated(error)) return
    endif

    ! These variables are different from the ones above in that they are
    ! allocatable, and in w90_wannier90_readwrite_read they were allocated on the root node only
    !
    if (.not. on_root) then
      allocate (fermi_energy_list(fermi_n), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating fermi_energy_read in postw90_w90_wannier90_readwrite_dist', comm)
        return
      endif
      !allocate (berry%kubo_freq_list(berry%kubo_nfreq), stat=ierr)
      !if (ierr /= 0) then
      !  call set_error_alloc(error, 'Error allocating kubo_freq_list in postw90_w90_wannier90_readwrite_dist', comm)
      !  return
      !endif
      !allocate (dos_data%project(dos_data%num_project), stat=ierr)
      !if (ierr /= 0) then
      !  call set_error_alloc(error, 'Error allocating dos_project in postw90_w90_wannier90_readwrite_dist', comm)
      !  return
      !endif
      !if (.not. pw90_common%effective_model) then
      if (eig_found) then
        allocate (eigval(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating eigval in postw90_w90_wannier90_readwrite_dist', comm)
          return
        endif
      end if
      allocate (kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating kpt_latt in postw90_w90_wannier90_readwrite_dist', comm)
        return
      endif
      !endif
      !allocate (gyrotropic%band_list(gyrotropic%num_bands), stat=ierr)
      !if (ierr /= 0) call set_error_alloc( &
      !  'Error allocating gyrotropic_num_bands in postw90_w90_wannier90_readwrite_dist')
      !allocate (gyrotropic%freq_list(gyrotropic%nfreq), stat=ierr)
      !if (ierr /= 0) call set_error_alloc( &
      !  'Error allocating gyrotropic_freq_list in postw90_w90_wannier90_readwrite_dist')
    end if

    if (fermi_n > 0) then
      call comms_bcast(fermi_energy_list(1), fermi_n, error, comm)
      if (allocated(error)) return
    endif
    !if (berry%kubo_nfreq > 0) call comms_bcast(berry%kubo_freq_list(1), berry%kubo_nfreq)
    !call comms_bcast(gyrotropic%freq_list(1), gyrotropic%nfreq)
    !call comms_bcast(gyrotropic%band_list(1), gyrotropic%num_bands)
    !if (dos_data%num_project > 0) &
    !  call comms_bcast(dos_data%project(1), dos_data%num_project)
    !if (.not. pw90_common%effective_model) then
    if (eig_found) then
      call comms_bcast(eigval(1, 1), num_bands*num_kpts, error, comm)
      if (allocated(error)) return
    end if
    call comms_bcast(kpt_latt(1, 1), 3*num_kpts, error, comm)
    if (allocated(error)) return

    !if (.not. pw90_common%effective_model .and. .not. w90_calculation%explicit_nnkpts) then
    if (.not. kmesh_info%explicit_nnkpts) then

      call comms_bcast(kmesh_info%nnh, 1, error, comm)
      if (allocated(error)) return

      call comms_bcast(kmesh_info%nntot, 1, error, comm)
      if (allocated(error)) return

      call comms_bcast(kmesh_info%wbtot, 1, error, comm)
      if (allocated(error)) return

      if (.not. on_root) then
        allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating nnlist in w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%neigh(num_kpts, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating neigh in w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating nncell in w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%wb(kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating wb in w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%bka(3, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating bka in w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (kmesh_info%bk(3, kmesh_info%nntot, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating bk in w90_wannier90_readwrite_dist', comm)
          return
        endif
      end if

      call comms_bcast(kmesh_info%nnlist(1, 1), num_kpts*kmesh_info%nntot, error, comm)
      if (allocated(error)) return

      call comms_bcast(kmesh_info%neigh(1, 1), num_kpts*kmesh_info%nntot/2, error, comm)
      if (allocated(error)) return

      call comms_bcast(kmesh_info%nncell(1, 1, 1), 3*num_kpts*kmesh_info%nntot, error, comm)
      if (allocated(error)) return

      call comms_bcast(kmesh_info%wb(1), kmesh_info%nntot, error, comm)
      if (allocated(error)) return

      call comms_bcast(kmesh_info%bka(1, 1), 3*kmesh_info%nntot/2, error, comm)
      if (allocated(error)) return

      call comms_bcast(kmesh_info%bk(1, 1, 1), 3*kmesh_info%nntot*num_kpts, error, comm)
      if (allocated(error)) return
    endif

    call comms_bcast(wann_omega%total, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_omega%tilde, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(wann_omega%invariant, 1, error, comm)
    if (allocated(error)) return

    call comms_bcast(have_disentangled, 1, error, comm)
    if (allocated(error)) return

    if (.not. on_root) then
      allocate (wannier_data%centres(3, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating wannier_centres in w90_wannier90_readwrite_dist', comm)
        return
      endif
      wannier_data%centres = 0.0_dp
      allocate (wannier_data%spreads(num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating wannier_spreads in w90_wannier90_readwrite_dist', comm)
        return
      endif
      wannier_data%spreads = 0.0_dp
      if (disentanglement) then
        allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating ndimwin in w90_wannier90_readwrite_dist', comm)
          return
        endif
        allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating lwindow in w90_wannier90_readwrite_dist', comm)
          return
        endif
      endif
    endif

  end subroutine w90_wannier90_readwrite_dist

end module w90_wannier90_readwrite
