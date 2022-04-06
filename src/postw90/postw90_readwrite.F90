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
!  w90_postw90_readwrite: input and output routines          !
!     specific to postw90.x                                  !
!                                                            !
!------------------------------------------------------------!

module w90_postw90_readwrite

  !! Read/write routines specific to postw90.x data types

  use w90_constants, only: dp, maxlen
  use w90_types, only: print_output_type, print_output_type, wannier_data_type, &
    kmesh_input_type, kmesh_info_type, dis_manifold_type, atom_data_type, kpoint_path_type, &
    proj_input_type, w90_system_type, ws_region_type
  use w90_readwrite
  use w90_postw90_types
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
    set_error_input, set_error_fatal, set_error_file

  implicit none

  private

  ! These could be local to w90_wannier90_readwrite_read if they weren't also used by w90_wannier90_readwrite_write
  type pw90_extra_io_type
    ! from gyrotropic section
    real(kind=dp) :: gyrotropic_freq_min = 0.0_dp
    real(kind=dp) :: gyrotropic_freq_max
    real(kind=dp) :: gyrotropic_freq_step = 0.01_dp
    real(kind=dp) :: kubo_freq_min = 0.0_dp
    real(kind=dp) :: kubo_freq_max
    real(kind=dp) :: kubo_freq_step = 0.01_dp
    ! Adaptive vs. fixed smearing stuff [GP, Jul 12, 2012]
    ! Only internal, always use the local variables defined by each module
    ! that take this value as default
    type(pw90_smearing_type) :: smear
    ! [gp-begin, Apr 13, 2012]
    ! Global interpolation k mesh variables
    ! These don't need to be public, since their values are copied in the variables of the
    ! local interpolation meshes. JRY: added save attribute
    type(kmesh_spacing_type) :: global_kmesh
    logical :: global_kmesh_set = .false.
    ! [gp-end]
    character(len=4) :: boltz_2d_dir ! this could be local to read_boltzwann
  end type pw90_extra_io_type

  public :: pw90_extra_io_type
  public :: w90_postw90_readwrite_read
  public :: w90_postw90_readwrite_write

contains

  !================================================!

  subroutine w90_postw90_readwrite_read(ws_region, w90_system, exclude_bands, print_output, &
                                        wannier_data, kmesh_input, kpt_latt, num_kpts, &
                                        dis_manifold, fermi_energy_list, atom_data, num_bands, &
                                        num_wann, eigval, mp_grid, real_lattice, kpoint_path, &
                                        pw90_calculation, pw90_oper_read, scissors_shift, &
                                        effective_model, pw90_spin, pw90_band_deriv_degen, &
                                        pw90_kpath, pw90_kslice, pw90_dos, pw90_berry, &
                                        pw90_spin_hall, pw90_gyrotropic, pw90_geninterp, &
                                        pw90_boltzwann, eig_found, pw90_extra_io, gamma_only, &
                                        bohr, optimisation, stdout, seedname, error, comm)
    !================================================!
    !
    !! Read parameters and calculate derived values
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !
    !================================================
    use w90_utility, only: utility_recip_lattice
    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none

    ! arguments
    type(atom_data_type), intent(inout) :: atom_data
    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(pw90_boltzwann_type), intent(inout) :: pw90_boltzwann
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(pw90_dos_mod_type), intent(inout) :: pw90_dos
    type(pw90_geninterp_mod_type), intent(inout) :: pw90_geninterp
    type(pw90_gyrotropic_type), intent(inout) :: pw90_gyrotropic
    type(pw90_kpath_mod_type), intent(inout) :: pw90_kpath
    type(pw90_kslice_mod_type), intent(inout) :: pw90_kslice
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(pw90_band_deriv_degen_type), intent(inout) :: pw90_band_deriv_degen
    type(pw90_oper_read_type), intent(inout) :: pw90_oper_read
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    type(print_output_type), intent(inout) :: print_output
    type(pw90_calculation_type), intent(inout) :: pw90_calculation
    type(pw90_extra_io_type), intent(inout) :: pw90_extra_io
    type(ws_region_type), intent(inout) :: ws_region
    type(kpoint_path_type), intent(inout) :: kpoint_path
    type(pw90_spin_hall_type), intent(inout) :: pw90_spin_hall
    type(w90_system_type), intent(inout) :: w90_system
    type(wannier_data_type), intent(inout) :: wannier_data
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer, intent(inout) :: mp_grid(3)
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_kpts
    integer, intent(inout) :: num_wann
    integer, intent(inout) :: optimisation
    integer, intent(in) :: stdout
    integer, allocatable, intent(inout) :: exclude_bands(:)

    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: scissors_shift
    real(kind=dp), allocatable, intent(inout) :: fermi_energy_list(:)
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)

    character(len=50), intent(in)  :: seedname

    logical, intent(inout) :: eig_found
    logical, intent(inout) :: gamma_only
    logical, intent(inout) :: effective_model

    ! local variables
    real(kind=dp) :: recip_lattice(3, 3), volume
    integer :: num_exclude_bands
    logical :: dos_plot
    logical :: found_fermi_energy
    logical :: disentanglement, library, ok
    character(len=20) :: energy_unit

    library = .false.
    pw90_kslice%corner = 0.0_dp
    pw90_kslice%b1 = [1.0_dp, 0.0_dp, 0.0_dp]
    pw90_kslice%b2 = [0.0_dp, 1.0_dp, 0.0_dp]
    pw90_kslice%kmesh2d = 50

    call w90_readwrite_read_verbosity(print_output, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_algorithm_control(optimisation, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_pw90_calcs(pw90_calculation, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_effective_model(effective_model, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_units(print_output%lenconfac, print_output%length_unit, energy_unit, &
                                  bohr, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_oper(pw90_oper_read, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_num_wann(num_wann, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_exclude_bands(exclude_bands, num_exclude_bands, error, comm) !for read_chkpt
    if (allocated(error)) return
    call w90_readwrite_read_num_bands(effective_model, library, num_exclude_bands, num_bands, &
                                      num_wann, .false., stdout, error, comm)
    if (allocated(error)) return
    disentanglement = (num_bands > num_wann)
    !call w90_readwrite_read_devel(print_output%devel_flag, stdout, seedname)
    call w90_readwrite_read_mp_grid(effective_model, library, mp_grid, num_kpts, stdout, &
                                    error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_gamma_only(gamma_only, num_kpts, library, stdout, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_system(library, w90_system, stdout, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_kpath(library, kpoint_path, ok, .false., error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_fermi_energy(found_fermi_energy, fermi_energy_list, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_kslice(pw90_calculation%kslice, pw90_kslice, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_smearing(pw90_extra_io%smear, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_scissors_shift(scissors_shift, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_pw90spin(pw90_calculation%spin_moment, &
                                               pw90_calculation%spin_decomp, pw90_spin, &
                                               w90_system%num_elec_per_state, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_gyrotropic(pw90_gyrotropic, num_wann, &
                                                 pw90_extra_io%smear%fixed_width, &
                                                 pw90_extra_io%smear%type_index, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_berry(pw90_calculation, pw90_berry, pw90_extra_io%smear, &
                                            error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_spin_hall(pw90_calculation, scissors_shift, pw90_spin_hall, &
                                                pw90_berry%task, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_pw90ham(pw90_band_deriv_degen, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_pw90_kpath(pw90_calculation, pw90_kpath, kpoint_path, &
                                                 error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_dos(pw90_calculation, pw90_dos, found_fermi_energy, &
                                          num_wann, pw90_extra_io%smear, dos_plot, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_ws_data(ws_region, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_eigvals(effective_model, pw90_calculation%boltzwann, &
                                    pw90_calculation%geninterp, dos_plot, disentanglement, &
                                    eig_found, eigval, library, .false., num_bands, num_kpts, &
                                    stdout, seedname, error, comm)
    if (allocated(error)) return
    dis_manifold%win_min = -1.0_dp
    dis_manifold%win_max = 0.0_dp
    if (eig_found) dis_manifold%win_min = minval(eigval)
    if (eig_found) dis_manifold%win_max = maxval(eigval)
    call w90_readwrite_read_dis_manifold(eig_found, dis_manifold, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_geninterp(pw90_geninterp, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_boltzwann(pw90_boltzwann, eigval, pw90_extra_io%smear, &
                                                pw90_calculation%boltzwann, &
                                                pw90_extra_io%boltz_2d_dir, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_energy_range(pw90_berry, pw90_dos, pw90_gyrotropic, &
                                                   dis_manifold, fermi_energy_list, eigval, &
                                                   pw90_extra_io, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_lattice(library, real_lattice, bohr, stdout, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_kmesh_data(kmesh_input, error, comm)
    if (allocated(error)) return
    call utility_recip_lattice(real_lattice, recip_lattice, volume, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_kpoints(effective_model, library, kpt_latt, num_kpts, &
                                    bohr, stdout, error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_global_kmesh(pw90_extra_io%global_kmesh_set, &
                                                   pw90_extra_io%global_kmesh, recip_lattice, &
                                                   error, comm)
    if (allocated(error)) return
    call w90_wannier90_readwrite_read_local_kmesh(pw90_calculation, pw90_berry, pw90_dos, &
                                                  pw90_spin, pw90_gyrotropic, pw90_boltzwann, &
                                                  recip_lattice, pw90_extra_io%global_kmesh_set, &
                                                  pw90_extra_io%global_kmesh, error, comm)
    if (allocated(error)) return
    call w90_readwrite_read_atoms(library, atom_data, real_lattice, bohr, stdout, error, comm)
    if (allocated(error)) return
  end subroutine w90_postw90_readwrite_read

  !================================================!
  subroutine w90_wannier90_readwrite_read_pw90_calcs(pw90_calculation, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(pw90_calculation_type), intent(inout) :: pw90_calculation
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm
    logical :: found

    call w90_readwrite_get_keyword('dos', found, error, comm, l_value=pw90_calculation%dos)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('berry', found, error, comm, l_value=pw90_calculation%berry)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('kpath', found, error, comm, l_value=pw90_calculation%kpath)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('kslice', found, error, comm, l_value=pw90_calculation%kslice)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('gyrotropic', found, error, comm, &
                                   l_value=pw90_calculation%gyrotropic)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('geninterp', found, error, comm, &
                                   l_value=pw90_calculation%geninterp)
    if (allocated(error)) return
    call w90_readwrite_get_keyword('boltzwann', found, error, comm, &
                                   l_value=pw90_calculation%boltzwann)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_pw90_calcs

  !================================================!
  subroutine w90_wannier90_readwrite_read_effective_model(effective_model, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    logical, intent(inout) :: effective_model
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('effective_model', found, error, comm, l_value=effective_model)
    if (allocated(error)) return
  end subroutine w90_wannier90_readwrite_read_effective_model

  !================================================!
  subroutine w90_wannier90_readwrite_read_oper(pw90_oper_read, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(pw90_oper_read_type), intent(inout) :: pw90_oper_read
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    ! formatted or "binary" file
    call w90_readwrite_get_keyword('spn_formatted', found, error, comm, &
                                   l_value=pw90_oper_read%spn_formatted)
    if (allocated(error)) return

    ! formatted or "binary" file
    call w90_readwrite_get_keyword('uhu_formatted', found, error, comm, &
                                   l_value=pw90_oper_read%uHu_formatted)
    if (allocated(error)) return
  end subroutine w90_wannier90_readwrite_read_oper

  !================================================!
  subroutine w90_wannier90_readwrite_read_kslice(kslicel, pw90_kslice, error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none
    logical, intent(in) :: kslicel
    type(pw90_kslice_mod_type), intent(inout) :: pw90_kslice
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer :: i
    logical :: found

    call w90_readwrite_get_keyword('kslice_task', found, error, comm, c_value=pw90_kslice%task)
    if (allocated(error)) return
    if (kslicel .and. index(pw90_kslice%task, 'fermi_lines') == 0 .and. &
        index(pw90_kslice%task, 'curv') == 0 .and. &
        index(pw90_kslice%task, 'morb') == 0 .and. &
        index(pw90_kslice%task, 'shc') == 0) then
      call set_error_input(error, 'Error: value of kslice_task not recognised in w90_wannier90_readwrite_read', comm)
      return
    endif
    if (kslicel .and. index(pw90_kslice%task, 'curv') > 0 .and. &
        index(pw90_kslice%task, 'morb') > 0) then
      call set_error_input(error, "Error: kslice_task cannot include both 'curv' and 'morb'", comm)
      return
    endif
    if (kslicel .and. index(pw90_kslice%task, 'shc') > 0 .and. &
        index(pw90_kslice%task, 'morb') > 0) then
      call set_error_input(error, "Error: kslice_task cannot include both 'shc' and 'morb'", comm)
      return
    endif
    if (kslicel .and. index(pw90_kslice%task, 'shc') > 0 .and. &
        index(pw90_kslice%task, 'curv') > 0) then
      call set_error_input(error, "Error: kslice_task cannot include both 'shc' and 'curv'", comm)
      return
    endif

    call w90_readwrite_get_vector_length('kslice_2dkmesh', found, i, error, comm)
    if (allocated(error)) return
    if (found) then
      if (i == 1) then
        call w90_readwrite_get_keyword_vector('kslice_2dkmesh', found, 1, error, comm, &
                                              i_value=pw90_kslice%kmesh2d)
        if (allocated(error)) return
        pw90_kslice%kmesh2d(2) = pw90_kslice%kmesh2d(1)
      elseif (i == 2) then
        call w90_readwrite_get_keyword_vector('kslice_2dkmesh', found, 2, error, comm, &
                                              i_value=pw90_kslice%kmesh2d)
        if (allocated(error)) return
      else
        call set_error_input(error, 'Error: kslice_2dkmesh must be provided as either' &
                             //' one integer or a vector of two integers', comm)
        return
      endif
      if (any(pw90_kslice%kmesh2d <= 0)) then
        call set_error_input(error, 'Error: kslice_2dkmesh elements must be greater than zero', &
                             comm)
        return
      endif
    endif

    call w90_readwrite_get_keyword_vector('kslice_corner', found, 3, error, comm, &
                                          r_value=pw90_kslice%corner)
    if (allocated(error)) return

    call w90_readwrite_get_keyword_vector('kslice_b1', found, 3, error, comm, &
                                          r_value=pw90_kslice%b1)
    if (allocated(error)) return

    call w90_readwrite_get_keyword_vector('kslice_b2', found, 3, error, comm, &
                                          r_value=pw90_kslice%b2)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('kslice_fermi_lines_colour', found, error, comm, &
                                   c_value=pw90_kslice%fermi_lines_colour)
    if (allocated(error)) return
    if (kslicel .and. index(pw90_kslice%fermi_lines_colour, 'none') == 0 .and. &
        index(pw90_kslice%fermi_lines_colour, 'spin') == 0) then
      call set_error_input(error, 'Error: value of kslice_fermi_lines_colour not recognised ' &
                           //'in w90_wannier90_readwrite_read', comm)
      return
    endif
!    slice_plot_format         = 'plotmv'
!    call w90_readwrite_get_keyword('slice_plot_format',found,c_value=slice_plot_format)
  end subroutine w90_wannier90_readwrite_read_kslice

  !================================================!
  subroutine w90_wannier90_readwrite_read_smearing(pw90_smearing, error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none
    type(pw90_smearing_type), intent(inout) :: pw90_smearing
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found
    character(len=maxlen)              :: ctmp
    ! [gp-begin, Apr 20, 2012]

    ! By default: Gaussian
    call w90_readwrite_get_keyword('smr_type', found, error, comm, c_value=ctmp)
    if (allocated(error)) return
    if (found) then
      pw90_smearing%type_index = w90_readwrite_get_smearing_index(ctmp, 'smr_type', error, comm)
      if (allocated(error)) return
    endif

    ! By default: adaptive smearing
    call w90_readwrite_get_keyword('adpt_smr', found, error, comm, l_value=pw90_smearing%use_adaptive)
    if (allocated(error)) return

    ! By default: a=sqrt(2)
    call w90_readwrite_get_keyword('adpt_smr_fac', found, error, comm, &
                                   r_value=pw90_smearing%adaptive_prefactor)
    if (allocated(error)) return
    if (found .and. (pw90_smearing%adaptive_prefactor <= 0._dp)) then
      call set_error_input(error, 'Error: adpt_smr_fac must be greater than zero', comm)
      return
    endif

    ! By default: 1 eV
    call w90_readwrite_get_keyword('adpt_smr_max', found, error, comm, &
                                   r_value=pw90_smearing%adaptive_max_width)
    if (allocated(error)) return
    if (pw90_smearing%adaptive_max_width <= 0._dp) then
      call set_error_input(error, 'Error: adpt_smr_max must be greater than zero', comm)
      return
    endif

    ! By default: if adpt_smr is manually set to false by the user, but he/she doesn't
    ! define smr_fixed_en_width: NO smearing, i.e. just the histogram
    call w90_readwrite_get_keyword('smr_fixed_en_width', found, error, comm, &
                                   r_value=pw90_smearing%fixed_width)
    if (allocated(error)) return
    if (found .and. (pw90_smearing%fixed_width < 0._dp)) then
      call set_error_input(error, 'Error: smr_fixed_en_width must be greater than or equal to zero', comm)
      return
    endif
  end subroutine w90_wannier90_readwrite_read_smearing

  !================================================!
  subroutine w90_wannier90_readwrite_read_scissors_shift(scissors_shift, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    real(kind=dp), intent(inout) :: scissors_shift
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('scissors_shift', found, error, comm, r_value=scissors_shift)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_scissors_shift

  !================================================!
  subroutine w90_wannier90_readwrite_read_pw90spin(spin_moment, spin_decomp, pw90_spin, &
                                                   num_elec_per_state, error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none
    logical, intent(inout) :: spin_moment ! from pw90_calculation
    logical, intent(inout) :: spin_decomp ! from pw90_common
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    integer, intent(in) :: num_elec_per_state
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('spin_moment', found, error, comm, l_value=spin_moment)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('spin_axis_polar', found, error, comm, r_value=pw90_spin%axis_polar)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('spin_axis_azimuth', found, error, comm, &
                                   r_value=pw90_spin%axis_azimuth)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('spin_decomp', found, error, comm, l_value=spin_decomp)
    if (allocated(error)) return

    if (spin_decomp .and. (num_elec_per_state .ne. 1)) then
      call set_error_input(error, 'spin_decomp can be true only if num_elec_per_state is 1', comm)
      return
    end if

  end subroutine w90_wannier90_readwrite_read_pw90spin

  !================================================!
  subroutine w90_wannier90_readwrite_read_gyrotropic(pw90_gyrotropic, num_wann, &
                                                     smr_fixed_en_width, smr_index, error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none
    type(pw90_gyrotropic_type), intent(out) :: pw90_gyrotropic
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: smr_fixed_en_width
    integer, intent(in) :: smr_index
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    real(kind=dp) :: smr_max_arg
    real(kind=dp)                   :: gyrotropic_box_tmp(3)
    integer :: i, ierr, loop
    logical :: found
    character(len=maxlen)              :: ctmp

    ! Stepan
    call w90_readwrite_get_keyword('gyrotropic_task', found, error, comm, c_value=pw90_gyrotropic%task)
    if (allocated(error)) return
    call w90_readwrite_get_keyword('gyrotropic_degen_thresh', found, error, comm, &
                                   r_value=pw90_gyrotropic%degen_thresh)
    if (allocated(error)) return

    do i = 1, 3
      pw90_gyrotropic%box(i, i) = 1.0_dp
      gyrotropic_box_tmp(:) = 0.0_dp
      call w90_readwrite_get_keyword_vector('gyrotropic_box_b'//achar(48 + i), found, 3, error, &
                                            comm, r_value=gyrotropic_box_tmp)
      if (allocated(error)) return
      if (found) pw90_gyrotropic%box(i, :) = gyrotropic_box_tmp(:)
    enddo
    call w90_readwrite_get_keyword_vector('gyrotropic_box_center', found, 3, error, comm, &
                                          r_value=gyrotropic_box_tmp)
    if (allocated(error)) return
    if (found) pw90_gyrotropic%box_corner(:) = &
      gyrotropic_box_tmp(:) - 0.5*(pw90_gyrotropic%box(1, :) + pw90_gyrotropic%box(2, :) + &
                                   pw90_gyrotropic%box(3, :))

    call w90_readwrite_get_range_vector('gyrotropic_band_list', found, pw90_gyrotropic%num_bands, &
                                        .true., error, comm)
    if (allocated(error)) return
    if (found) then
      if (pw90_gyrotropic%num_bands < 1) then
        call set_error_input(error, 'Error: problem reading gyrotropic_band_list', comm)
        return
      endif
      if (allocated(pw90_gyrotropic%band_list)) deallocate (pw90_gyrotropic%band_list)
      allocate (pw90_gyrotropic%band_list(pw90_gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating gyrotropic_band_list in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_range_vector('gyrotropic_band_list', found, &
                                          pw90_gyrotropic%num_bands, .false., error, comm, &
                                          pw90_gyrotropic%band_list)
      if (allocated(error)) return
      if (any(pw90_gyrotropic%band_list < 1) .or. any(pw90_gyrotropic%band_list > num_wann)) then
        call set_error_input(error, 'Error: gyrotropic_band_list asks for a non-valid bands', comm)
        return
      endif
    else
      ! include all bands in the calculation
      pw90_gyrotropic%num_bands = num_wann
      if (allocated(pw90_gyrotropic%band_list)) deallocate (pw90_gyrotropic%band_list)
      allocate (pw90_gyrotropic%band_list(pw90_gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating gyrotropic_band_list in w90_wannier90_readwrite_read', comm)
        return
      endif
      do loop = 1, num_wann
        pw90_gyrotropic%band_list(loop) = loop
      end do
    end if

    pw90_gyrotropic%smearing%use_adaptive = .false.
    smr_max_arg = 5.0
    call w90_readwrite_get_keyword('smr_max_arg', found, error, comm, r_value=smr_max_arg)
    if (allocated(error)) return
    if (found .and. (smr_max_arg <= 0._dp)) then
      call set_error_input(error, 'Error: smr_max_arg must be greater than zero', comm)
      return
    endif

    pw90_gyrotropic%smearing%max_arg = smr_max_arg
    call w90_readwrite_get_keyword('gyrotropic_smr_max_arg', found, error, comm, &
                                   r_value=pw90_gyrotropic%smearing%max_arg)
    if (allocated(error)) return
    if (found .and. (pw90_gyrotropic%smearing%max_arg <= 0._dp)) then
      call set_error_input(error, 'Error: gyrotropic_smr_max_arg must be greater than zero', comm)
      return
    endif

    pw90_gyrotropic%smearing%fixed_width = smr_fixed_en_width
    call w90_readwrite_get_keyword('gyrotropic_smr_fixed_en_width', found, error, comm, &
                                   r_value=pw90_gyrotropic%smearing%fixed_width)
    if (allocated(error)) return
    if (found .and. (pw90_gyrotropic%smearing%fixed_width < 0._dp)) then
      call set_error_input(error, 'Error: gyrotropic_smr_fixed_en_width must be greater than or equal to zero', comm)
      return
    endif

    ! By default: use the "global" smearing index
    pw90_gyrotropic%smearing%type_index = smr_index
    call w90_readwrite_get_keyword('gyrotropic_smr_type', found, error, comm, c_value=ctmp)
    if (allocated(error)) return
    if (found) then
      pw90_gyrotropic%smearing%type_index = w90_readwrite_get_smearing_index(ctmp, &
                                                                             'gyrotropic_smr_type', &
                                                                             error, comm)
      if (allocated(error)) return
    endif

  end subroutine w90_wannier90_readwrite_read_gyrotropic

  !================================================!
  subroutine w90_wannier90_readwrite_read_berry(pw90_calculation, pw90_berry, pw90_smearing, &
                                                error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none
    type(pw90_calculation_type), intent(in) :: pw90_calculation
    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(pw90_smearing_type), intent(in) :: pw90_smearing
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found
    integer :: kdotp_num_bands, ierr
    character(len=maxlen)              :: ctmp

!-------------------------------------------------------
!    alpha=0
!    call w90_readwrite_get_keyword('alpha',found,i_value=alpha)

!    beta=0
!    call w90_readwrite_get_keyword('beta',found,i_value=beta)

!    gamma=0
!    call w90_readwrite_get_keyword('gamma',found,i_value=gamma)
!-------------------------------------------------------

    call w90_readwrite_get_keyword('transl_inv', found, error, comm, l_value=pw90_berry%transl_inv)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('berry_task', found, error, comm, c_value=pw90_berry%task)
    if (allocated(error)) return
    if (pw90_calculation%berry .and. .not. found) then
      call set_error_input(error, 'Error: berry=T and berry_task is not set', comm)
      return
    endif
    if (pw90_calculation%berry .and. index(pw90_berry%task, 'ahc') == 0 &
        .and. index(pw90_berry%task, 'morb') == 0 &
        .and. index(pw90_berry%task, 'kubo') == 0 .and. index(pw90_berry%task, 'sc') == 0 &
        .and. index(pw90_berry%task, 'shc') == 0 .and. index(pw90_berry%task, 'kdotp') == 0) then

      call set_error_input(error, 'Error: value of berry_task not recognised in w90_wannier90_readwrite_read', comm)
      return
    endif

    call w90_readwrite_get_keyword('berry_curv_adpt_kmesh', found, error, comm, &
                                   i_value=pw90_berry%curv_adpt_kmesh)
    if (allocated(error)) return
    if (pw90_berry%curv_adpt_kmesh < 1) then
      call set_error_input(error, 'Error:  berry_curv_adpt_kmesh must be a positive integer', comm)
      return
    endif

    call w90_readwrite_get_keyword('berry_curv_adpt_kmesh_thresh', found, error, comm, &
                                   r_value=pw90_berry%curv_adpt_kmesh_thresh)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('berry_curv_unit', found, error, comm, c_value=pw90_berry%curv_unit)
    if (allocated(error)) return
    if (pw90_berry%curv_unit .ne. 'ang2' .and. pw90_berry%curv_unit .ne. 'bohr2') then
      call set_error_input(error, 'Error: value of berry_curv_unit not recognised in w90_wannier90_readwrite_read', comm)
      return
    endif

    call w90_readwrite_get_keyword('wanint_kpoint_file', found, error, comm, &
                                   l_value=pw90_berry%wanint_kpoint_file)
    if (allocated(error)) return

!    smear_temp = -1.0_dp
!    call w90_readwrite_get_keyword('smear_temp',found,r_value=smear_temp)

    pw90_berry%kubo_smearing%use_adaptive = pw90_smearing%use_adaptive
    call w90_readwrite_get_keyword('kubo_adpt_smr', found, error, comm, &
                                   l_value=pw90_berry%kubo_smearing%use_adaptive)
    if (allocated(error)) return

    pw90_berry%kubo_smearing%adaptive_prefactor = pw90_smearing%adaptive_prefactor
    call w90_readwrite_get_keyword('kubo_adpt_smr_fac', found, error, comm, &
                                   r_value=pw90_berry%kubo_smearing%adaptive_prefactor)
    if (allocated(error)) return
    if (found .and. (pw90_berry%kubo_smearing%adaptive_prefactor <= 0._dp)) then
      call set_error_input(error, 'Error: kubo_adpt_smr_fac must be greater than zero', comm)
      return
    endif

    pw90_berry%kubo_smearing%adaptive_max_width = pw90_smearing%adaptive_max_width
    call w90_readwrite_get_keyword('kubo_adpt_smr_max', found, error, comm, &
                                   r_value=pw90_berry%kubo_smearing%adaptive_max_width)
    if (allocated(error)) return
    if (pw90_berry%kubo_smearing%adaptive_max_width <= 0._dp) then
      call set_error_input(error, 'Error: kubo_adpt_smr_max must be greater than zero', comm)
      return
    endif

    pw90_berry%kubo_smearing%fixed_width = pw90_smearing%fixed_width
    call w90_readwrite_get_keyword('kubo_smr_fixed_en_width', found, error, comm, &
                                   r_value=pw90_berry%kubo_smearing%fixed_width)
    if (allocated(error)) return
    if (found .and. (pw90_berry%kubo_smearing%fixed_width < 0._dp)) then
      call set_error_input(error, 'Error: kubo_smr_fixed_en_width must be greater than or equal to zero', comm)
      return
    endif

    call w90_readwrite_get_keyword('sc_phase_conv', found, error, comm, &
                                   i_value=pw90_berry%sc_phase_conv)
    if (allocated(error)) return
    if ((pw90_berry%sc_phase_conv .ne. 1) .and. ((pw90_berry%sc_phase_conv .ne. 2))) then
      call set_error_input(error, 'Error: sc_phase_conv must be either 1 or 2', comm)
      return
    endif

    call w90_readwrite_get_keyword('sc_use_eta_corr', found, error, comm, &
                                   l_value=pw90_berry%sc_use_eta_corr)
    if (allocated(error)) return

    ! By default: use the "global" smearing index
    pw90_berry%kubo_smearing%type_index = pw90_smearing%type_index
    call w90_readwrite_get_keyword('kubo_smr_type', found, error, comm, c_value=ctmp)
    if (allocated(error)) return
    if (found) then
      pw90_berry%kubo_smearing%type_index = w90_readwrite_get_smearing_index(ctmp, 'kubo_smr_type', &
                                                                             error, comm)
      if (allocated(error)) return
    endif

    call w90_readwrite_get_keyword('sc_eta', found, error, comm, r_value=pw90_berry%sc_eta)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('sc_w_thr', found, error, comm, r_value=pw90_berry%sc_w_thr)
    if (allocated(error)) return

    call w90_readwrite_get_keyword_vector('kdotp_kpoint', found, 3, error, comm, &
                                          r_value=pw90_berry%kdotp_kpoint)
    if (allocated(error)) return

    kdotp_num_bands = 0
    call w90_readwrite_get_keyword('kdotp_num_bands', found, error, comm, i_value=kdotp_num_bands)
    if (allocated(error)) return
    if (found) then
      if (kdotp_num_bands < 1) then
        call set_error_input(error, 'Error: problem reading kdotp_num_bands', comm)
        return
      endif
      allocate (pw90_berry%kdotp_bands(kdotp_num_bands), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating kdotp_num_bands in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_range_vector('kdotp_bands', found, kdotp_num_bands, &
                                          .false., error, comm, pw90_berry%kdotp_bands)
      if (allocated(error)) return
      if (any(pw90_berry%kdotp_bands < 1)) then
        call set_error_input(error, 'Error: kdotp_bands must contain positive numbers', comm)
        return
      endif
    end if

  end subroutine w90_wannier90_readwrite_read_berry

  !================================================!
  subroutine w90_wannier90_readwrite_read_spin_hall(pw90_calculation, scissors_shift, &
                                                    pw90_spin_hall, berry_task, error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none

    type(pw90_calculation_type), intent(in) :: pw90_calculation
    type(pw90_spin_hall_type), intent(inout) :: pw90_spin_hall
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    real(kind=dp), intent(in) :: scissors_shift

    character(len=*), intent(in) :: berry_task

    logical :: found

    call w90_readwrite_get_keyword('shc_freq_scan', found, error, comm, &
                                   l_value=pw90_spin_hall%freq_scan)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('shc_alpha', found, error, comm, i_value=pw90_spin_hall%alpha)
    if (allocated(error)) return
    if (found .and. (pw90_spin_hall%alpha < 1 .or. pw90_spin_hall%alpha > 3)) then
      call set_error_input(error, 'Error:  shc_alpha must be 1, 2 or 3', comm)
      return
    endif

    call w90_readwrite_get_keyword('shc_beta', found, error, comm, i_value=pw90_spin_hall%beta)
    if (allocated(error)) return
    if (found .and. (pw90_spin_hall%beta < 1 .or. pw90_spin_hall%beta > 3)) then
      call set_error_input(error, 'Error:  shc_beta must be 1, 2 or 3', comm)
      return
    endif

    call w90_readwrite_get_keyword('shc_gamma', found, error, comm, i_value=pw90_spin_hall%gamma)
    if (allocated(error)) return
    if (found .and. (pw90_spin_hall%gamma < 1 .or. pw90_spin_hall%gamma > 3)) then
      call set_error_input(error, 'Error:  shc_gamma must be 1, 2 or 3', comm)
      return
    endif

    call w90_readwrite_get_keyword('shc_bandshift', found, error, comm, l_value=pw90_spin_hall%bandshift)
    if (allocated(error)) return
    pw90_spin_hall%bandshift = pw90_spin_hall%bandshift .and. pw90_calculation%berry .and. &
                               .not. (index(berry_task, 'shc') == 0)
    if ((abs(scissors_shift) > 1.0e-7_dp) .and. pw90_spin_hall%bandshift) then
      call set_error_input(error, 'Error: shc_bandshift and scissors_shift cannot be used simultaneously', comm)
      return
    endif

    call w90_readwrite_get_keyword('shc_bandshift_firstband', found, error, comm, &
                                   i_value=pw90_spin_hall%bandshift_firstband)
    if (allocated(error)) return
    if (pw90_spin_hall%bandshift .and. (.not. found)) then
      call set_error_input(error, 'Error: shc_bandshift required but no shc_bandshift_firstband provided', comm)
      return
    endif
    if ((pw90_spin_hall%bandshift_firstband < 1) .and. found) then
      call set_error_input(error, 'Error: shc_bandshift_firstband must >= 1', comm)
      return
    endif

    call w90_readwrite_get_keyword('shc_bandshift_energyshift', found, error, comm, &
                                   r_value=pw90_spin_hall%bandshift_energyshift)
    if (allocated(error)) return
    if (pw90_spin_hall%bandshift .and. (.not. found)) then
      call set_error_input(error, 'Error: shc_bandshift required but no shc_bandshift_energyshift provided', comm)
      return
    endif

    call w90_readwrite_get_keyword('shc_method', found, error, comm, c_value=pw90_spin_hall%method)
    if (allocated(error)) return
    if (index(berry_task, 'shc') > 0 .and. .not. found) then
      call set_error_input(error, 'Error: berry_task=shc and shc_method is not set', comm)
      return
    endif
    if (index(berry_task, 'shc') > 0 .and. index(pw90_spin_hall%method, 'qiao') == 0 &
        .and. index(pw90_spin_hall%method, 'ryoo') == 0) then
      call set_error_input(error, 'Error: value of shc_method not recognised in w90_wannier90_readwrite_read', comm)
      return
    endif

  end subroutine w90_wannier90_readwrite_read_spin_hall

  !================================================!
  subroutine w90_wannier90_readwrite_read_pw90ham(pw90_band_deriv_degen, error, comm)
    !================================================!
    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(pw90_band_deriv_degen_type), intent(inout) :: pw90_band_deriv_degen
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('use_degen_pert', found, error, comm, &
                                   l_value=pw90_band_deriv_degen%use_degen_pert)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('degen_thr', found, error, comm, &
                                   r_value=pw90_band_deriv_degen%degen_thr)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_pw90ham

  !================================================!
  subroutine w90_wannier90_readwrite_read_pw90_kpath(pw90_calculation, pw90_kpath, kpoint_path, &
                                                     error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none

    type(pw90_calculation_type), intent(in) :: pw90_calculation
    type(pw90_kpath_mod_type), intent(inout) :: pw90_kpath
    type(kpoint_path_type), intent(in) :: kpoint_path
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('kpath_task', found, error, comm, c_value=pw90_kpath%task)
    if (allocated(error)) return
    if (pw90_calculation%kpath .and. index(pw90_kpath%task, 'bands') == 0 .and. &
        index(pw90_kpath%task, 'curv') == 0 .and. &
        index(pw90_kpath%task, 'morb') == 0 .and. &
        index(pw90_kpath%task, 'shc') == 0) then
      call set_error_input(error, 'Error: value of kpath_task not recognised in w90_wannier90_readwrite_read', comm)
      return
    endif
    if (.not. allocated(kpoint_path%labels) .and. pw90_calculation%kpath) then
      call set_error_input(error, 'Error: a kpath plot has been requested but there is no kpoint_path block', comm)
      return
    endif

    call w90_readwrite_get_keyword('kpath_num_points', found, error, comm, &
                                   i_value=pw90_kpath%num_points)
    if (allocated(error)) return
    if (pw90_kpath%num_points < 0) then
      call set_error_input(error, 'Error: kpath_num_points must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword('kpath_bands_colour', found, error, comm, &
                                   c_value=pw90_kpath%bands_colour)
    if (allocated(error)) return
    if (pw90_calculation%kpath .and. index(pw90_kpath%bands_colour, 'none') == 0 .and. &
        index(pw90_kpath%bands_colour, 'spin') == 0 .and. &
        index(pw90_kpath%bands_colour, 'shc') == 0) then
      call set_error_input(error, 'Error: value of kpath_bands_colour not recognised in w90_wannier90_readwrite_read', comm)
      return
    endif
    if (pw90_calculation%kpath .and. index(pw90_kpath%task, 'shc') > 0 .and. &
        index(pw90_kpath%task, 'spin') > 0) then
      call set_error_input(error, "Error: kpath_task cannot include both 'shc' and 'spin'", comm)
      return
    endif

  end subroutine w90_wannier90_readwrite_read_pw90_kpath

  !================================================!
  subroutine w90_wannier90_readwrite_read_dos(pw90_calculation, pw90_dos, found_fermi_energy, &
                                              num_wann, pw90_smearing, dos_plot, error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none

    type(pw90_calculation_type), intent(in) :: pw90_calculation
    type(pw90_dos_mod_type), intent(inout) :: pw90_dos
    type(pw90_smearing_type), intent(in) :: pw90_smearing
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    logical, intent(out) :: dos_plot
    logical, intent(in) :: found_fermi_energy

    integer :: i, ierr
    logical :: found
    character(len=maxlen)              :: ctmp

    if (pw90_calculation%dos) then
      dos_plot = .true.
    else
      dos_plot = .false.
    endif
    call w90_readwrite_get_keyword('dos_task', found, error, comm, c_value=pw90_dos%task)
    if (allocated(error)) return
    if (pw90_calculation%dos) then
      if (index(pw90_dos%task, 'dos_plot') == 0 .and. &
          index(pw90_dos%task, 'find_fermi_energy') == 0) then
        call set_error_input(error, 'Error: value of dos_task not recognised in w90_wannier90_readwrite_read', comm)
        return
      endif
      if (index(pw90_dos%task, 'dos_plot') > 0) dos_plot = .true.
      if (index(pw90_dos%task, 'find_fermi_energy') > 0 .and. found_fermi_energy) then
        call set_error_input(error, 'Error: Cannot set "dos_task = find_fermi_energy" and give a value to "fermi_energy"', comm)
        return
      endif
    end if

!    sigma_abc_onlyorb=.false.
!    call w90_readwrite_get_keyword('sigma_abc_onlyorb',found,l_value=sigma_abc_onlyorb)

! -------------------------------------------------------------------

    !IVO_END

    call w90_readwrite_get_keyword('dos_energy_step', found, error, comm, &
                                   r_value=pw90_dos%energy_step)
    if (allocated(error)) return

    pw90_dos%smearing%use_adaptive = pw90_smearing%use_adaptive
    call w90_readwrite_get_keyword('dos_adpt_smr', found, error, comm, &
                                   l_value=pw90_dos%smearing%use_adaptive)
    if (allocated(error)) return

    pw90_dos%smearing%adaptive_prefactor = pw90_smearing%adaptive_prefactor
    call w90_readwrite_get_keyword('dos_adpt_smr_fac', found, error, comm, &
                                   r_value=pw90_dos%smearing%adaptive_prefactor)
    if (allocated(error)) return
    if (found .and. (pw90_dos%smearing%adaptive_prefactor <= 0._dp)) then
      call set_error_input(error, 'Error: dos_adpt_smr_fac must be greater than zero', comm)
      return
    endif

    pw90_dos%smearing%adaptive_max_width = pw90_smearing%adaptive_max_width
    call w90_readwrite_get_keyword('dos_adpt_smr_max', found, error, comm, &
                                   r_value=pw90_dos%smearing%adaptive_max_width)
    if (allocated(error)) return
    if (pw90_dos%smearing%adaptive_max_width <= 0._dp) then
      call set_error_input(error, 'Error: dos_adpt_smr_max must be greater than zero', comm)
      return
    endif

    pw90_dos%smearing%fixed_width = pw90_smearing%fixed_width
    call w90_readwrite_get_keyword('dos_smr_fixed_en_width', found, error, comm, &
                                   r_value=pw90_dos%smearing%fixed_width)
    if (allocated(error)) return
    if (found .and. (pw90_dos%smearing%fixed_width < 0._dp)) then
      call set_error_input(error, 'Error: dos_smr_fixed_en_width must be greater than or equal to zero', comm)
      return
    endif

!    dos_gaussian_width        = 0.1_dp
!    call w90_readwrite_get_keyword('dos_gaussian_width',found,r_value=dos_gaussian_width)

!    dos_plot_format           = 'gnuplot'
!    call w90_readwrite_get_keyword('dos_plot_format',found,c_value=dos_plot_format)

    call w90_readwrite_get_range_vector('dos_project', found, pw90_dos%num_project, .true., &
                                        error, comm)
    if (allocated(error)) return
    if (found) then
      if (pw90_dos%num_project < 1) then
        call set_error_input(error, 'Error: problem reading dos_project', comm)
        return
      endif
      if (allocated(pw90_dos%project)) deallocate (pw90_dos%project)
      allocate (pw90_dos%project(pw90_dos%num_project), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating dos_project in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_range_vector('dos_project', found, pw90_dos%num_project, .false., &
                                          error, comm, pw90_dos%project)
      if (allocated(error)) return
      if (any(pw90_dos%project < 1) .or. &
          any(pw90_dos%project > num_wann)) then
        call set_error_input(error, 'Error: dos_project asks for out-of-range Wannier functions', comm)
        return
      endif
    else
      ! by default plot all
      pw90_dos%num_project = num_wann
      if (allocated(pw90_dos%project)) deallocate (pw90_dos%project)
      allocate (pw90_dos%project(pw90_dos%num_project), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating dos_project in w90_wannier90_readwrite_read', comm)
        return
      endif
      do i = 1, pw90_dos%num_project
        pw90_dos%project(i) = i
      end do
    endif

    ! By default: use the "global" smearing index
    pw90_dos%smearing%type_index = pw90_smearing%type_index
    call w90_readwrite_get_keyword('dos_smr_type', found, error, comm, c_value=ctmp)
    if (allocated(error)) return
    if (found) then
      pw90_dos%smearing%type_index = w90_readwrite_get_smearing_index(ctmp, 'dos_smr_type', &
                                                                      error, comm)
      if (allocated(error)) return
    endif

  end subroutine w90_wannier90_readwrite_read_dos

  !================================================!
  subroutine w90_wannier90_readwrite_read_geninterp(pw90_geninterp, error, comm)
    !================================================!
    ! [gp-begin, Jun 1, 2012]
    ! General band interpolator (pw90_geninterp)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none

    type(pw90_geninterp_mod_type), intent(inout) :: pw90_geninterp
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found

    call w90_readwrite_get_keyword('geninterp_alsofirstder', found, error, comm, &
                                   l_value=pw90_geninterp%alsofirstder)
    if (allocated(error)) return
    call w90_readwrite_get_keyword('geninterp_single_file', found, error, comm, &
                                   l_value=pw90_geninterp%single_file)
    if (allocated(error)) return
    ! [gp-end, Jun 1, 2012]

  end subroutine w90_wannier90_readwrite_read_geninterp

  !================================================!
  subroutine w90_wannier90_readwrite_read_boltzwann(pw90_boltzwann, eigval, pw90_smearing, &
                                                    do_boltzwann, boltz_2d_dir, error, comm)
    !================================================!
    ! [gp-begin, Jun 1, 2012]
    ! General band interpolator (pw90_geninterp)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none
    type(pw90_boltzwann_type), intent(inout) :: pw90_boltzwann
    type(pw90_smearing_type), intent(in) :: pw90_smearing
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    real(kind=dp), allocatable, intent(in) :: eigval(:, :)
    logical, intent(in) :: do_boltzwann
    character(len=4), intent(out) :: boltz_2d_dir

    logical :: found, found2
    character(len=maxlen)              :: ctmp

    ! [gp-begin, Apr 12, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! Boltzmann transport
    !%%%%%%%%%%%%%%%%%%%%
    ! Note: to be put AFTER the disentanglement routines!
    pw90_boltzwann%TDF_smearing%use_adaptive = .false.

    call w90_readwrite_get_keyword('boltz_calc_also_dos', found, error, comm, &
                                   l_value=pw90_boltzwann%calc_also_dos)
    if (allocated(error)) return

    pw90_boltzwann%calc_also_dos = pw90_boltzwann%calc_also_dos .and. do_boltzwann

    ! 0 means the normal 3d case for the calculation of the Seebeck coefficient
    ! The other valid possibilities are 1,2,3 for x,y,z respectively
    call w90_readwrite_get_keyword('boltz_2d_dir', found, error, comm, c_value=boltz_2d_dir)
    if (allocated(error)) return
    if (found) then
      if (trim(boltz_2d_dir) == 'no') then
        pw90_boltzwann%dir_num_2d = 0
      elseif (trim(boltz_2d_dir) == 'x') then
        pw90_boltzwann%dir_num_2d = 1
      elseif (trim(boltz_2d_dir) == 'y') then
        pw90_boltzwann%dir_num_2d = 2
      elseif (trim(boltz_2d_dir) == 'z') then
        pw90_boltzwann%dir_num_2d = 3
      else
        call set_error_input(error, 'Error: boltz_2d_dir can only be "no", "x", "y" or "z".', comm)
        return
      end if
    end if

    call w90_readwrite_get_keyword('boltz_dos_energy_step', found, error, comm, &
                                   r_value=pw90_boltzwann%dos_energy_step)
    if (allocated(error)) return
    if (found .and. (pw90_boltzwann%dos_energy_step <= 0._dp)) then
      call set_error_input(error, 'Error: boltz_dos_energy_step must be positive', comm)
      return
    endif

    if (allocated(eigval)) then
      pw90_boltzwann%dos_energy_min = minval(eigval) - 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      pw90_boltzwann%dos_energy_min = -1.0_dp
    end if
    call w90_readwrite_get_keyword('boltz_dos_energy_min', found, error, comm, &
                                   r_value=pw90_boltzwann%dos_energy_min)
    if (allocated(error)) return
    if (allocated(eigval)) then
      pw90_boltzwann%dos_energy_max = maxval(eigval) + 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      pw90_boltzwann%dos_energy_max = 0.0_dp
    end if
    call w90_readwrite_get_keyword('boltz_dos_energy_max', found, error, comm, &
                                   r_value=pw90_boltzwann%dos_energy_max)
    if (allocated(error)) return
    if (pw90_boltzwann%dos_energy_max <= pw90_boltzwann%dos_energy_min) then
      call set_error_input(error, 'Error: boltz_dos_energy_max must be greater than boltz_dos_energy_min', comm)
      return
    endif

    pw90_boltzwann%dos_smearing%use_adaptive = pw90_smearing%use_adaptive
    call w90_readwrite_get_keyword('boltz_dos_adpt_smr', found, error, comm, &
                                   l_value=pw90_boltzwann%dos_smearing%use_adaptive)
    if (allocated(error)) return

    pw90_boltzwann%dos_smearing%adaptive_prefactor = pw90_smearing%adaptive_prefactor
    call w90_readwrite_get_keyword('boltz_dos_adpt_smr_fac', found, error, comm, &
                                   r_value=pw90_boltzwann%dos_smearing%adaptive_prefactor)
    if (allocated(error)) return
    if (found .and. (pw90_boltzwann%dos_smearing%adaptive_prefactor <= 0._dp)) then
      call set_error_input(error, 'Error: boltz_dos_adpt_smr_fac must be greater than zero', comm)
      return
    endif

    pw90_boltzwann%dos_smearing%adaptive_max_width = pw90_smearing%adaptive_max_width
    call w90_readwrite_get_keyword('boltz_dos_adpt_smr_max', found, error, comm, &
                                   r_value=pw90_boltzwann%dos_smearing%adaptive_max_width)
    if (allocated(error)) return
    if (pw90_boltzwann%dos_smearing%adaptive_max_width <= 0._dp) then
      call set_error_input(error, 'Error: boltz_dos_adpt_smr_max must be greater than zero', comm)
      return
    endif

    pw90_boltzwann%dos_smearing%fixed_width = pw90_smearing%fixed_width
    call w90_readwrite_get_keyword('boltz_dos_smr_fixed_en_width', found, error, comm, &
                                   r_value=pw90_boltzwann%dos_smearing%fixed_width)
    if (allocated(error)) return
    if (found .and. (pw90_boltzwann%dos_smearing%fixed_width < 0._dp)) then
      call set_error_input(error, 'Error: boltz_dos_smr_fixed_en_width must be greater than or equal to zero', comm)
      return
    endif

    call w90_readwrite_get_keyword('boltz_mu_min', found, error, comm, r_value=pw90_boltzwann%mu_min)
    if (allocated(error)) return
    if ((.not. found) .and. do_boltzwann) then
      call set_error_input(error, 'Error: BoltzWann required but no boltz_mu_min provided', comm)
      return
    endif
    call w90_readwrite_get_keyword('boltz_mu_max', found2, error, comm, &
                                   r_value=pw90_boltzwann%mu_max)
    if (allocated(error)) return
    if ((.not. found2) .and. do_boltzwann) then
      call set_error_input(error, 'Error: BoltzWann required but no boltz_mu_max provided', comm)
      return
    endif
    if (found .and. found2 .and. (pw90_boltzwann%mu_max < pw90_boltzwann%mu_min)) then
      call set_error_input(error, 'Error: boltz_mu_max must be greater than boltz_mu_min', comm)
      return
    endif
    call w90_readwrite_get_keyword('boltz_mu_step', found, error, comm, r_value=pw90_boltzwann%mu_step)
    if (allocated(error)) return
    if ((.not. found) .and. do_boltzwann) then
      call set_error_input(error, 'Error: BoltzWann required but no boltz_mu_step provided', comm)
      return
    endif
    if (found .and. (pw90_boltzwann%mu_step <= 0._dp)) then
      call set_error_input(error, 'Error: boltz_mu_step must be greater than zero', comm)
      return
    endif

    call w90_readwrite_get_keyword('boltz_temp_min', found, error, comm, &
                                   r_value=pw90_boltzwann%temp_min)
    if (allocated(error)) return
    if ((.not. found) .and. do_boltzwann) then
      call set_error_input(error, 'Error: BoltzWann required but no boltz_temp_min provided', comm)
      return
    endif
    call w90_readwrite_get_keyword('boltz_temp_max', found2, error, comm, &
                                   r_value=pw90_boltzwann%temp_max)
    if (allocated(error)) return
    if ((.not. found2) .and. do_boltzwann) then
      call set_error_input(error, 'Error: BoltzWann required but no boltz_temp_max provided', comm)
      return
    endif
    if (found .and. found2 .and. (pw90_boltzwann%temp_max < pw90_boltzwann%temp_min)) then
      call set_error_input(error, 'Error: boltz_temp_max must be greater than boltz_temp_min', comm)
      return
    endif
    if (found .and. (pw90_boltzwann%temp_min <= 0._dp)) then
      call set_error_input(error, 'Error: boltz_temp_min must be greater than zero', comm)
      return
    endif
    call w90_readwrite_get_keyword('boltz_temp_step', found, error, comm, &
                                   r_value=pw90_boltzwann%temp_step)
    if (allocated(error)) return
    if ((.not. found) .and. do_boltzwann) then
      call set_error_input(error, 'Error: BoltzWann required but no boltz_temp_step provided', comm)
      return
    endif
    if (found .and. (pw90_boltzwann%temp_step <= 0._dp)) then
      call set_error_input(error, 'Error: boltz_temp_step must be greater than zero', comm)
      return
    endif

    ! The interpolation mesh is read later on

    ! By default, the energy step for the TDF is 1 meV
    call w90_readwrite_get_keyword('boltz_tdf_energy_step', found, error, comm, &
                                   r_value=pw90_boltzwann%tdf_energy_step)
    if (allocated(error)) return
    if (pw90_boltzwann%tdf_energy_step <= 0._dp) then
      call set_error_input(error, 'Error: boltz_tdf_energy_step must be greater than zero', comm)
      return
    endif

    ! For TDF: TDF smeared in a NON-adaptive way; value in eV, default = 0._dp
    ! (i.e., no smearing)
    pw90_boltzwann%tdf_smearing%fixed_width = pw90_smearing%fixed_width
    call w90_readwrite_get_keyword('boltz_tdf_smr_fixed_en_width', found, error, comm, &
                                   r_value=pw90_boltzwann%tdf_smearing%fixed_width)
    if (allocated(error)) return
    if (found .and. (pw90_boltzwann%tdf_smearing%fixed_width < 0._dp)) then
      call set_error_input(error, 'Error: boltz_TDF_smr_fixed_en_width must be greater than or equal to zero', comm)
      return
    endif

    ! By default: use the "global" smearing index
    pw90_boltzwann%tdf_smearing%type_index = pw90_smearing%type_index
    call w90_readwrite_get_keyword('boltz_tdf_smr_type', found, error, comm, c_value=ctmp)
    if (allocated(error)) return
    if (found) then
      pw90_boltzwann%tdf_smearing%type_index = &
        w90_readwrite_get_smearing_index(ctmp, 'boltz_tdf_smr_type', error, comm)
      if (allocated(error)) return
    endif

    ! By default: use the "global" smearing index
    pw90_boltzwann%dos_smearing%type_index = pw90_smearing%type_index
    call w90_readwrite_get_keyword('boltz_dos_smr_type', found, error, comm, c_value=ctmp)
    if (allocated(error)) return
    if (found) then
      pw90_boltzwann%dos_smearing%type_index = &
        w90_readwrite_get_smearing_index(ctmp, 'boltz_dos_smr_type', error, comm)
    endif

    ! By default: 10 fs relaxation time
    call w90_readwrite_get_keyword('boltz_relax_time', found, error, comm, &
                                   r_value=pw90_boltzwann%relax_time)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('boltz_bandshift', found, error, comm, &
                                   l_value=pw90_boltzwann%bandshift)
    if (allocated(error)) return
    pw90_boltzwann%bandshift = pw90_boltzwann%bandshift .and. do_boltzwann

    call w90_readwrite_get_keyword('boltz_bandshift_firstband', found, error, comm, &
                                   i_value=pw90_boltzwann%bandshift_firstband)
    if (allocated(error)) return
    if (pw90_boltzwann%bandshift .and. (.not. found)) then
      call set_error_input(error, 'Error: boltz_bandshift required but no boltz_bandshift_firstband provided', comm)
      return
    endif
    call w90_readwrite_get_keyword('boltz_bandshift_energyshift', found, error, comm, &
                                   r_value=pw90_boltzwann%bandshift_energyshift)
    if (allocated(error)) return
    if (pw90_boltzwann%bandshift .and. (.not. found)) then
      call set_error_input(error, 'Error: boltz_bandshift required but no boltz_bandshift_energyshift provided', comm)
      return
    endif
  end subroutine w90_wannier90_readwrite_read_boltzwann

  !================================================!
  subroutine w90_wannier90_readwrite_read_energy_range(pw90_berry, pw90_dos, pw90_gyrotropic, &
                                                       dis_manifold, fermi_energy_list, eigval, &
                                                       pw90_extra_io, error, comm)
    !================================================!

    use w90_constants, only: cmplx_i
    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none

    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(pw90_dos_mod_type), intent(inout) :: pw90_dos
    type(pw90_gyrotropic_type), intent(inout) :: pw90_gyrotropic
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), allocatable, intent(in) :: eigval(:, :)
    type(pw90_extra_io_type), intent(inout) :: pw90_extra_io

    integer :: i, ierr
    logical :: found

    if (dis_manifold%frozen_states) then
      pw90_dos%energy_max = dis_manifold%froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      pw90_dos%energy_max = maxval(eigval) + 0.6667_dp
    else
      pw90_dos%energy_max = dis_manifold%win_max + 0.6667_dp
    end if
    call w90_readwrite_get_keyword('dos_energy_max', found, error, comm, &
                                   r_value=pw90_dos%energy_max)
    if (allocated(error)) return

    if (allocated(eigval)) then
      pw90_dos%energy_min = minval(eigval) - 0.6667_dp
    else
      pw90_dos%energy_min = dis_manifold%win_min - 0.6667_dp
    end if
    call w90_readwrite_get_keyword('dos_energy_min', found, error, comm, &
                                   r_value=pw90_dos%energy_min)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('kubo_freq_min', found, error, comm, &
                                   r_value=pw90_extra_io%kubo_freq_min)
    if (allocated(error)) return

    if (dis_manifold%frozen_states) then
      pw90_extra_io%kubo_freq_max = dis_manifold%froz_max - fermi_energy_list(1) + 0.6667_dp
    elseif (allocated(eigval)) then
      pw90_extra_io%kubo_freq_max = maxval(eigval) - minval(eigval) + 0.6667_dp
    else
      pw90_extra_io%kubo_freq_max = dis_manifold%win_max - dis_manifold%win_min + 0.6667_dp
    end if
    pw90_extra_io%gyrotropic_freq_max = pw90_extra_io%kubo_freq_max
    call w90_readwrite_get_keyword('kubo_freq_max', found, error, comm, &
                                   r_value=pw90_extra_io%kubo_freq_max)
    if (allocated(error)) return

    call w90_readwrite_get_keyword('kubo_freq_step', found, error, comm, &
                                   r_value=pw90_extra_io%kubo_freq_step)
    if (allocated(error)) return
    if (found .and. pw90_extra_io%kubo_freq_step < 0.0_dp) then
      call set_error_input(error, 'Error: kubo_freq_step must be positive', comm)
      return
    endif

    pw90_berry%kubo_nfreq = nint((pw90_extra_io%kubo_freq_max - pw90_extra_io%kubo_freq_min) &
                                 /pw90_extra_io%kubo_freq_step) + 1
    if (pw90_berry%kubo_nfreq <= 1) pw90_berry%kubo_nfreq = 2
    pw90_extra_io%kubo_freq_step = (pw90_extra_io%kubo_freq_max - pw90_extra_io%kubo_freq_min) &
                                   /(pw90_berry%kubo_nfreq - 1)

    if (allocated(pw90_berry%kubo_freq_list)) deallocate (pw90_berry%kubo_freq_list)
    allocate (pw90_berry%kubo_freq_list(pw90_berry%kubo_nfreq), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating kubo_freq_list in w90_wannier90_readwrite_read', comm)
      return
    endif
    do i = 1, pw90_berry%kubo_nfreq
      pw90_berry%kubo_freq_list(i) = pw90_extra_io%kubo_freq_min + &
                                     (i - 1)*(pw90_extra_io%kubo_freq_max - &
                                              pw90_extra_io%kubo_freq_min)/(pw90_berry%kubo_nfreq - 1)
    enddo

    ! TODO: Alternatively, read list of (complex) frequencies; kubo_nfreq is
    !       the length of the list

    call w90_readwrite_get_keyword('gyrotropic_freq_min', found, error, comm, &
                                   r_value=pw90_extra_io%gyrotropic_freq_min)
    if (allocated(error)) return
    call w90_readwrite_get_keyword('gyrotropic_freq_max', found, error, comm, &
                                   r_value=pw90_extra_io%gyrotropic_freq_max)
    if (allocated(error)) return
    call w90_readwrite_get_keyword('gyrotropic_freq_step', found, error, comm, &
                                   r_value=pw90_extra_io%gyrotropic_freq_step)
    if (allocated(error)) return
    pw90_gyrotropic%nfreq = nint((pw90_extra_io%gyrotropic_freq_max - &
                                  pw90_extra_io%gyrotropic_freq_min)/ &
                                 pw90_extra_io%gyrotropic_freq_step) + 1
    if (pw90_gyrotropic%nfreq <= 1) pw90_gyrotropic%nfreq = 2
    pw90_extra_io%gyrotropic_freq_step = (pw90_extra_io%gyrotropic_freq_max &
                                          - pw90_extra_io%gyrotropic_freq_min)/ &
                                         (pw90_gyrotropic%nfreq - 1)
    if (allocated(pw90_gyrotropic%freq_list)) deallocate (pw90_gyrotropic%freq_list)
    allocate (pw90_gyrotropic%freq_list(pw90_gyrotropic%nfreq), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating gyrotropic_freq_list in w90_wannier90_readwrite_read', comm)
      return
    endif
    do i = 1, pw90_gyrotropic%nfreq
      pw90_gyrotropic%freq_list(i) = pw90_extra_io%gyrotropic_freq_min &
                                     + (i - 1)*(pw90_extra_io%gyrotropic_freq_max &
                                                - pw90_extra_io%gyrotropic_freq_min)/(pw90_gyrotropic%nfreq - 1) &
                                     + cmplx_i*pw90_gyrotropic%smearing%fixed_width
    enddo

    if (dis_manifold%frozen_states) then
      pw90_berry%kubo_eigval_max = dis_manifold%froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      pw90_berry%kubo_eigval_max = maxval(eigval) + 0.6667_dp
    else
      pw90_berry%kubo_eigval_max = dis_manifold%win_max + 0.6667_dp
    end if
    pw90_gyrotropic%eigval_max = pw90_berry%kubo_eigval_max

    call w90_readwrite_get_keyword('kubo_eigval_max', found, error, comm, &
                                   r_value=pw90_berry%kubo_eigval_max)
    if (allocated(error)) return
    call w90_readwrite_get_keyword('gyrotropic_eigval_max', found, error, comm, &
                                   r_value=pw90_gyrotropic%eigval_max)
    if (allocated(error)) return
  end subroutine w90_wannier90_readwrite_read_energy_range

  !================================================!
  subroutine w90_wannier90_readwrite_read_global_kmesh(global_kmesh_set, kmesh, recip_lattice, &
                                                       error, comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none

    type(kmesh_spacing_type), intent(out) :: kmesh

    logical, intent(inout) :: global_kmesh_set
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer :: i
    logical :: found

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! k meshes
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! [GP-begin, Apr13, 2012]
    ! Global interpolation k-mesh; this is overridden by "local" meshes of a given submodule
    ! This bit of code must appear *before* all other codes for the local interpolation meshes,
    ! BUT *after* having calculated the reciprocal-space vectors.
    global_kmesh_set = .false.
    kmesh%spacing = -1._dp
    kmesh%mesh = 0
    call w90_readwrite_get_keyword('kmesh_spacing', found, error, comm, r_value=kmesh%spacing)
    if (allocated(error)) return
    if (found) then
      if (kmesh%spacing .le. 0._dp) then
        call set_error_input(error, 'Error: kmesh_spacing must be greater than zero', comm)
        return
      endif
      global_kmesh_set = .true.

      call w90_readwrite_set_kmesh(kmesh%spacing, recip_lattice, kmesh%mesh)
    end if
    call w90_readwrite_get_vector_length('kmesh', found, i, error, comm)
    if (allocated(error)) return
    if (found) then
      if (global_kmesh_set) then
        call set_error_input(error, 'Error: cannot set both kmesh and kmesh_spacing', comm)
        return
      endif
      if (i .eq. 1) then
        global_kmesh_set = .true.
        call w90_readwrite_get_keyword_vector('kmesh', found, 1, error, comm, i_value=kmesh%mesh)
        if (allocated(error)) return
        kmesh%mesh(2) = kmesh%mesh(1)
        kmesh%mesh(3) = kmesh%mesh(1)
      elseif (i .eq. 3) then
        global_kmesh_set = .true.
        call w90_readwrite_get_keyword_vector('kmesh', found, 3, error, comm, i_value=kmesh%mesh)
        if (allocated(error)) return
      else
        call set_error_input(error, 'Error: kmesh must be provided as either one integer or a vector of three integers', comm)
        return
      end if
      if (any(kmesh%mesh <= 0)) then
        call set_error_input(error, 'Error: kmesh elements must be greater than zero', comm)
        return
      endif
    end if
    ! [GP-end]
  end subroutine w90_wannier90_readwrite_read_global_kmesh

  !================================================!
  subroutine w90_wannier90_readwrite_read_local_kmesh(pw90_calculation, pw90_berry, pw90_dos, &
                                                      pw90_spin, pw90_gyrotropic, pw90_boltzwann, &
                                                      recip_lattice, global_kmesh_set, &
                                                      global_kmesh, error, comm)
    !================================================!
    use w90_comms, only: w90comm_type
    implicit none

    type(pw90_calculation_type), intent(in) :: pw90_calculation
    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(pw90_dos_mod_type), intent(inout) :: pw90_dos
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    type(pw90_gyrotropic_type), intent(inout) :: pw90_gyrotropic
    type(pw90_boltzwann_type), intent(inout) :: pw90_boltzwann
    type(kmesh_spacing_type), intent(in) :: global_kmesh
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    logical, intent(in) :: global_kmesh_set

    ! To be called after having read the global flag
    call get_module_kmesh(recip_lattice, global_kmesh_set, global_kmesh, error, comm, &
                          moduleprefix='boltz', should_be_defined=pw90_calculation%boltzwann, &
                          module_kmesh=pw90_boltzwann%kmesh)
    if (allocated(error)) return

    call get_module_kmesh(recip_lattice, global_kmesh_set, global_kmesh, error, comm, &
                          moduleprefix='berry', should_be_defined=pw90_calculation%berry, &
                          module_kmesh=pw90_berry%kmesh)
    if (allocated(error)) return

    call get_module_kmesh(recip_lattice, global_kmesh_set, global_kmesh, error, comm, &
                          moduleprefix='gyrotropic', &
                          should_be_defined=pw90_calculation%gyrotropic, &
                          module_kmesh=pw90_gyrotropic%kmesh)
    if (allocated(error)) return

    call get_module_kmesh(recip_lattice, global_kmesh_set, global_kmesh, error, comm, &
                          moduleprefix='spin', should_be_defined=pw90_calculation%spin_moment, &
                          module_kmesh=pw90_spin%kmesh)
    if (allocated(error)) return

    call get_module_kmesh(recip_lattice, global_kmesh_set, global_kmesh, error, comm, &
                          moduleprefix='dos', should_be_defined=pw90_calculation%dos, &
                          module_kmesh=pw90_dos%kmesh)
    if (allocated(error)) return

  end subroutine w90_wannier90_readwrite_read_local_kmesh

  !================================================!
  subroutine get_module_kmesh(recip_lattice, global_kmesh_set, global_kmesh, error, comm, &
                              moduleprefix, should_be_defined, module_kmesh)
    !================================================!
    !! This function reads and sets the interpolation mesh variables needed by a given module
    !>
    !!  This function MUST be called after having read the global kmesh and kmesh_spacing!!
    !!  if the user didn't provide an interpolation_mesh_spacing, it is set to -1, so that
    !!       one can check in the code what the user asked for
    !!  The function takes care also of setting the default value to the global one if no local
    !!       keyword is defined
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    character(len=*), intent(in)       :: moduleprefix
    !!The prefix that is appended before the name of the variables. In particular,
    !!if the prefix is for instance XXX, the two variables that are read from the
    !!input file are XXX_kmesh and XXX_kmesh_spacing.
    logical, intent(in)                :: should_be_defined
    !! A logical flag: if it is true, at the end the code stops if no value is specified.
    !! Define it to .false. if no check should be performed.
    !! Often, you can simply pass the flag that activates the module itself.
    type(kmesh_spacing_type), intent(inout) :: module_kmesh
    !! the integer array (length 3) where the interpolation mesh will be saved, and
    !! the real number on which the min mesh spacing is saved. -1 if it the
    !!user specifies in input the mesh and not the mesh_spacing
    logical, intent(in) :: global_kmesh_set
    type(kmesh_spacing_type), intent(in) :: global_kmesh
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    logical :: found, found2
    integer :: i

    ! Default values
    module_kmesh%spacing = -1._dp
    module_kmesh%mesh = 0
    call w90_readwrite_get_keyword(trim(moduleprefix)//'_kmesh_spacing', found, error, comm, &
                                   r_value=module_kmesh%spacing)
    if (allocated(error)) return
    if (found) then
      if (module_kmesh%spacing .le. 0._dp) then
        call set_error_input(error, 'Error: '//trim(moduleprefix)//'_kmesh_spacing must be greater than zero', comm)
        return
      endif

      call w90_readwrite_set_kmesh(module_kmesh%spacing, recip_lattice, module_kmesh%mesh)
    end if
    call w90_readwrite_get_vector_length(trim(moduleprefix)//'_kmesh', found2, i, error, comm)
    if (allocated(error)) return
    if (found2) then
      if (found) then
        call set_error_input(error, 'Error: cannot set both '//trim(moduleprefix)//'_kmesh and ' &
                             //trim(moduleprefix)//'_kmesh_spacing', comm)
        return
      endif
      if (i .eq. 1) then
        call w90_readwrite_get_keyword_vector(trim(moduleprefix)//'_kmesh', found2, &
                                              1, error, comm, i_value=module_kmesh%mesh)
        if (allocated(error)) return
        module_kmesh%mesh(2) = module_kmesh%mesh(1)
        module_kmesh%mesh(3) = module_kmesh%mesh(1)
      elseif (i .eq. 3) then
        call w90_readwrite_get_keyword_vector(trim(moduleprefix)//'_kmesh', found2, &
                                              3, error, comm, i_value=module_kmesh%mesh)
        if (allocated(error)) return
      else
        call set_error_input(error, 'Error: '//trim(moduleprefix)// &
                             '_kmesh must be provided as either one integer or a vector of 3 integers', &
                             comm)
        return
      end if
      if (any(module_kmesh%mesh <= 0)) then
        call set_error_input(error, 'Error: '//trim(moduleprefix)//'_kmesh elements must be greater than zero', comm)
        return
      endif
    end if

    if ((found .eqv. .false.) .and. (found2 .eqv. .false.)) then
      ! This is the case where no  "local" interpolation k-mesh is provided in the input
      if (global_kmesh_set) then
        module_kmesh%mesh = global_kmesh%mesh
        ! I set also boltz_kmesh_spacing so that I can check if it is < 0 or not, and if it is
        ! > 0 I can print on output the mesh spacing that was chosen
        module_kmesh%spacing = global_kmesh%spacing
      else
        if (should_be_defined) then
          call set_error_input(error, 'Error: '//trim(moduleprefix)//' module required, but no interpolation mesh given.', comm)
          return
        endif
      end if
    end if
  end subroutine get_module_kmesh

  !================================================
  subroutine w90_postw90_readwrite_write(print_output, w90_system, fermi_energy_list, atom_data, &
                                         num_wann, real_lattice, kpoint_path, pw90_calculation, &
                                         pw90_oper_read, scissors_shift, pw90_spin, pw90_kpath, &
                                         pw90_kslice, pw90_dos, pw90_berry, pw90_gyrotropic, &
                                         pw90_geninterp, pw90_boltzwann, pw90_extra_io, &
                                         optimisation, stdout)
    !================================================!
    !
    !! write postw90 parameters to stdout
    !
    !================================================
    use w90_utility, only: utility_recip_lattice_base, utility_inverse_mat, utility_cart_to_frac

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(w90_system_type), intent(in) :: w90_system
    type(atom_data_type), intent(in) :: atom_data
    type(kpoint_path_type), intent(in) :: kpoint_path
    type(pw90_calculation_type), intent(in) :: pw90_calculation
    type(pw90_oper_read_type), intent(in) :: pw90_oper_read
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(pw90_kpath_mod_type), intent(in) :: pw90_kpath
    type(pw90_kslice_mod_type), intent(in) :: pw90_kslice
    type(pw90_dos_mod_type), intent(in) :: pw90_dos
    type(pw90_berry_mod_type), intent(in) :: pw90_berry
    type(pw90_gyrotropic_type), intent(in) :: pw90_gyrotropic
    type(pw90_geninterp_mod_type), intent(in) :: pw90_geninterp
    type(pw90_boltzwann_type), intent(in) :: pw90_boltzwann
    type(pw90_extra_io_type), intent(in) :: pw90_extra_io

    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift
    integer, intent(in) :: num_wann
    integer, intent(in) :: optimisation
    integer, intent(in) :: stdout

    ! local variables
    real(kind=dp) :: recip_lattice(3, 3), inv_lattice(3, 3), pos_frac(3), volume
    real(kind=dp) :: cell_volume
    integer :: i, loop, nat, nsp

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
    cell_volume = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - &
                                      real_lattice(3, 2)*real_lattice(2, 3)) + &
                  real_lattice(1, 2)*(real_lattice(2, 3)*real_lattice(3, 1) - &
                                      real_lattice(3, 3)*real_lattice(2, 1)) + &
                  real_lattice(1, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - &
                                      real_lattice(3, 1)*real_lattice(2, 2))
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
    write (stdout, *) ' '
    ! Main
    write (stdout, *) ' '

    write (stdout, '(1x,a78)') '*-------------------------------- POSTW90 -----------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Wannier Functions               :', num_wann, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of electrons per state             :', &
      w90_system%num_elec_per_state, '|'
    if (abs(scissors_shift) > 1.0e-7_dp .or. print_output%iprint > 0) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Scissor shift applied to conduction bands :', scissors_shift, '|'
      if (w90_system%num_valence_bands > 0) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Number of valence bands                   :', &
          w90_system%num_valence_bands, '|'
      else
        write (stdout, '(1x,a78)') '|  Number of valence bands                   :       not defined             |'
      endif
    endif
    if (pw90_calculation%spin_decomp .or. print_output%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Spin decomposition                        :', pw90_calculation%spin_decomp, '|'
    if (pw90_calculation%spin_moment .or. print_output%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Spin moment                       :', pw90_calculation%spin_moment, '|'
    if (pw90_calculation%spin_decomp .or. pw90_calculation%spin_moment .or. print_output%iprint > 2) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Polar angle of spin quantisation axis     :', pw90_spin%axis_polar, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Azimuthal angle of spin quantisation axis :', pw90_spin%axis_azimuth, '|'
      if (pw90_oper_read%spn_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  Spn file-type                   :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  Spn file-type                   :', 'unformatted', '|'
      endif
      if (pw90_oper_read%uHu_formatted) then
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

    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :', print_output%iprint, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :', print_output%timing_level, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Optimisation (0=memory, 3=speed)          :', optimisation, '|'
    write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :', trim(print_output%length_unit), '|'
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    write (stdout, '(1x,a78)') '*------------------------ Global Smearing Parameters ------------------------*'
    if (pw90_extra_io%smear%use_adaptive) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', &
        pw90_extra_io%smear%adaptive_prefactor, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width (eV)       :', &
        pw90_extra_io%smear%adaptive_max_width, '|'

    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
        pw90_extra_io%smear%fixed_width, '|'
    endif
    write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
      trim(w90_readwrite_get_smearing_type(pw90_extra_io%smear%type_index)), '|'
    if (pw90_extra_io%global_kmesh_set) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       T', '|'
      if (pw90_extra_io%global_kmesh%spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          pw90_extra_io%global_kmesh%mesh(1), 'x', pw90_extra_io%global_kmesh%mesh(2), 'x', &
          pw90_extra_io%global_kmesh%mesh(3), ' Spacing = ', pw90_extra_io%global_kmesh%spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , pw90_extra_io%global_kmesh%mesh(1), 'x', pw90_extra_io%global_kmesh%mesh(2), 'x', &
          pw90_extra_io%global_kmesh%mesh(3), '|'
      endif
    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       F', '|'
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    ! DOS
    if (pw90_calculation%dos .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*---------------------------------- DOS -------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Density of States                :', pw90_calculation%dos, '|'
      if (pw90_dos%num_project > 1) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', pw90_dos%energy_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', pw90_dos%energy_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', pw90_dos%energy_step, '|'
      if (pw90_dos%smearing%use_adaptive .eqv. pw90_extra_io%smear%use_adaptive .and. &
          pw90_dos%smearing%adaptive_prefactor == pw90_extra_io%smear%adaptive_prefactor .and. &
          pw90_dos%smearing%adaptive_max_width == pw90_extra_io%smear%adaptive_max_width .and. &
          pw90_dos%smearing%fixed_width == pw90_extra_io%smear%fixed_width .and. &
          pw90_extra_io%smear%type_index == pw90_dos%smearing%type_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (pw90_dos%smearing%use_adaptive) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', &
            pw90_dos%smearing%adaptive_prefactor, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', &
            pw90_dos%smearing%adaptive_max_width, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
            pw90_dos%smearing%fixed_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
          trim(w90_readwrite_get_smearing_type(pw90_dos%smearing%type_index)), '|'
      endif
      if (pw90_extra_io%global_kmesh%mesh(1) == pw90_dos%kmesh%mesh(1) .and. &
          pw90_extra_io%global_kmesh%mesh(2) == pw90_dos%kmesh%mesh(2) .and. &
          pw90_extra_io%global_kmesh%mesh(3) == pw90_dos%kmesh%mesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (pw90_dos%kmesh%spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            pw90_dos%kmesh%mesh(1), 'x', pw90_dos%kmesh%mesh(2), 'x', &
            pw90_dos%kmesh%mesh(3), ' Spacing = ', pw90_dos%kmesh%spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :', &
            pw90_dos%kmesh%mesh(1), 'x', pw90_dos%kmesh%mesh(2), 'x', &
            pw90_dos%kmesh%mesh(3), '|'
        endif
      endif
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    if (pw90_calculation%kpath .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KPATH ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a path in k-space   :', pw90_calculation%kpath, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', pw90_kpath%num_points, '|'
      if (index(pw90_kpath%task, 'bands') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       F', '|'
      endif
      if (index(pw90_kpath%task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       F', '|'
      endif
      if (index(pw90_kpath%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(pw90_kpath%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the bands    :', trim(pw90_kpath%bands_colour), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      write (stdout, '(1x,a78)') '|   K-space path sections:                                                   |'
      if (.not. allocated(kpoint_path%labels)) then
        write (stdout, '(1x,a78)') '|     None defined                                                           |'
      else
        do loop = 1, size(kpoint_path%labels), 2
          write (stdout, '(1x,a10,2x,a1,2x,3F7.3,5x,a3,2x,a1,2x,3F7.3,7x,a1)') '|    From:', &
            kpoint_path%labels(loop), (kpoint_path%points(i, loop), i=1, 3), &
            'To:', kpoint_path%labels(loop + 1), (kpoint_path%points(i, loop + 1), i=1, 3), '|'
        end do
      end if
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calculation%kslice .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KSLICE -----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a slice in k-space  :', pw90_calculation%kslice, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi level used for slice                :', fermi_energy_list(1), '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', pw90_kpath%num_points, '|'
      if (index(pw90_kslice%task, 'fermi_lines') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       F', '|'
      endif
      if (index(pw90_kslice%task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       F', '|'
      endif
      if (index(pw90_kslice%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(pw90_kslice%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the lines    :', &
        trim(pw90_kslice%fermi_lines_colour), '|'
      write (stdout, '(1x,a78)') '|  2D slice parameters (in reduced coordinates):                             |'
      write (stdout, '(1x,a14,2x,3F8.3,37x,a1)') '|     Corner: ', (pw90_kslice%corner(i), i=1, 3), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector1: ', (pw90_kslice%b1(i), i=1, 3), ' Divisions:', pw90_kslice%kmesh2d(1), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector2: ', (pw90_kslice%b2(i), i=1, 3), ' Divisions:', pw90_kslice%kmesh2d(1), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calculation%berry .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- BERRY ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Berry Phase related properties    :', pw90_calculation%berry, '|'
      if (index(pw90_berry%task, 'kubo') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       F', '|'
      endif
      if (index(pw90_berry%task, 'ahc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       F', '|'
      endif
      if (index(pw90_berry%task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       F', '|'
      endif
      if (index(pw90_berry%task, 'kdotp') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute k.p expansion coefficients        :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute k.p expansion coefficients        :', '       F', '|'
      endif
      if (index(pw90_berry%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       F', '|'
      endif
      if (index(pw90_berry%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for optical responses     :', &
        pw90_extra_io%kubo_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency for optical responses     :', &
        pw90_extra_io%kubo_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for optical responses           :', &
        pw90_extra_io%kubo_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue for optical responses    :', pw90_berry%kubo_eigval_max, '|'
      if (index(pw90_berry%task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing factor for shift current         :', pw90_berry%sc_eta, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Frequency theshold for shift current      :', pw90_berry%sc_w_thr, '|'
        write (stdout, '(1x,a46,1x,a27,3x,a1)') '|  Bloch sums                                :', &
          trim(w90_readwrite_get_convention_type(pw90_berry%sc_phase_conv)), '|'
        write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Finite eta correction for shift current   :', &
          pw90_berry%sc_use_eta_corr, '|'
      end if
      if (index(pw90_berry%task, 'kdotp') > 0) then
        write (stdout, '(1x,a46,10x,f8.3,1x,f8.3,1x,f8.3,1x,13x,a1)') '|  Chosen k-point kdotp_kpoint                 :', &
          pw90_berry%kdotp_kpoint(1), pw90_berry%kdotp_kpoint(2), pw90_berry%kdotp_kpoint(3), '|'
        write (stdout, '(1x,a46,10x,i4,13x,a1)') '|  kdotp_num_bands                             :', &
          size(pw90_berry%kdotp_bands), '|'
        write (stdout, '(1x,a46,10x,*(i4))') '|  kdotp_bands                                 :', &
          pw90_berry%kdotp_bands(:)
      end if
      if (pw90_berry%kubo_smearing%use_adaptive .eqv. pw90_extra_io%smear%use_adaptive .and. &
          pw90_berry%kubo_smearing%adaptive_prefactor == pw90_extra_io%smear%adaptive_prefactor .and. &
          pw90_berry%kubo_smearing%adaptive_max_width == pw90_extra_io%smear%adaptive_max_width &
          .and. pw90_berry%kubo_smearing%fixed_width == pw90_extra_io%smear%fixed_width .and. &
          pw90_extra_io%smear%type_index == pw90_berry%kubo_smearing%type_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (pw90_berry%kubo_smearing%use_adaptive) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', &
            pw90_berry%kubo_smearing%adaptive_prefactor, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', &
            pw90_berry%kubo_smearing%adaptive_max_width, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
            pw90_berry%kubo_smearing%fixed_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
          trim(w90_readwrite_get_smearing_type(pw90_berry%kubo_smearing%type_index)), '|'
      endif
      if (pw90_extra_io%global_kmesh%mesh(1) == pw90_berry%kmesh%mesh(1) .and. &
          pw90_extra_io%global_kmesh%mesh(2) == pw90_berry%kmesh%mesh(2) .and. &
          pw90_extra_io%global_kmesh%mesh(3) == pw90_berry%kmesh%mesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (pw90_berry%kmesh%spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            pw90_berry%kmesh%mesh(1), 'x', pw90_berry%kmesh%mesh(2), 'x', pw90_berry%kmesh%mesh(3), &
            ' Spacing = ', pw90_berry%kmesh%spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , pw90_berry%kmesh%mesh(1), 'x', pw90_berry%kmesh%mesh(2), 'x', pw90_berry%kmesh%mesh(3), '|'
        endif
      endif
      if (pw90_berry%curv_adpt_kmesh > 1) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Using an adaptive refinement mesh of size :', pw90_berry%curv_adpt_kmesh, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Threshold for adaptive refinement         :', &
          pw90_berry%curv_adpt_kmesh_thresh, '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    none', '|'
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calculation%gyrotropic .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- GYROTROPIC   ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '| Compute Gyrotropic properties              :', pw90_calculation%gyrotropic, '|'
      write (stdout, '(1x,a46,10x,a20,1x,a1)') '| gyrotropic_task                            :', pw90_gyrotropic%task, '|'
      call parameters_gyro_write_task(pw90_gyrotropic%task, '-d0', 'calculate the D tensor', stdout)
      call parameters_gyro_write_task(pw90_gyrotropic%task, '-dw', 'calculate the tildeD tensor', stdout)
      call parameters_gyro_write_task(pw90_gyrotropic%task, '-c', 'calculate the C tensor', stdout)
      call parameters_gyro_write_task(pw90_gyrotropic%task, '-k', 'calculate the K tensor', stdout)
      call parameters_gyro_write_task(pw90_gyrotropic%task, '-noa', 'calculate the interbad natural optical activity', stdout)
      call parameters_gyro_write_task(pw90_gyrotropic%task, '-dos', 'calculate the density of states', stdout)

      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for tildeD,NOA            :', &
        pw90_extra_io%gyrotropic_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency                           :', &
        pw90_extra_io%gyrotropic_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for frequency                   :', &
        pw90_extra_io%gyrotropic_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue                          :', &
        pw90_gyrotropic%eigval_max, '|'
      if (pw90_gyrotropic%smearing%fixed_width == pw90_extra_io%smear%fixed_width &
          .and. pw90_extra_io%smear%type_index == pw90_gyrotropic%smearing%type_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        write (stdout, '(1x,a78)') '|  Using local  smearing parameters                                          |'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
        pw90_gyrotropic%smearing%fixed_width, '|'
      write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function                         :', &
        trim(w90_readwrite_get_smearing_type(pw90_gyrotropic%smearing%type_index)), '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  degen_thresh                              :', &
        pw90_gyrotropic%degen_thresh, '|'

      if (pw90_extra_io%global_kmesh%mesh(1) == pw90_gyrotropic%kmesh%mesh(1) .and. &
          pw90_extra_io%global_kmesh%mesh(2) == pw90_gyrotropic%kmesh%mesh(2) .and. &
          pw90_extra_io%global_kmesh%mesh(3) == pw90_gyrotropic%kmesh%mesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      elseif (pw90_gyrotropic%kmesh%spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          pw90_gyrotropic%kmesh%mesh(1), 'x', pw90_gyrotropic%kmesh%mesh(2), 'x', pw90_gyrotropic%kmesh%mesh(3), &
          ' Spacing = ', pw90_gyrotropic%kmesh%spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , pw90_gyrotropic%kmesh%mesh(1), 'x', pw90_gyrotropic%kmesh%mesh(2), 'x', pw90_gyrotropic%kmesh%mesh(3), '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    not implemented', '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calculation%boltzwann .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------------- BOLTZWANN ----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Boltzmann transport properties    :', &
        pw90_calculation%boltzwann, '|'
      if (pw90_boltzwann%dir_num_2d > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  2d structure: non-periodic dimension  :', &
          trim(pw90_extra_io%boltz_2d_dir), '|'
      else
        write (stdout, '(1x,a78)') '|  3d Structure                              :                 T             |'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Relaxation Time (fs)                      :', pw90_boltzwann%relax_time, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Chemical Potential (eV)  :', pw90_boltzwann%mu_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Chemical Potential (eV)  :', pw90_boltzwann%mu_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Chemical Potential (eV)     :', pw90_boltzwann%mu_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Temperature (K)          :', pw90_boltzwann%temp_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Temperature (K)          :', pw90_boltzwann%temp_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Temperature (K)             :', pw90_boltzwann%temp_step, '|'

      if (pw90_extra_io%global_kmesh%mesh(1) == pw90_boltzwann%kmesh%mesh(1) .and. &
          pw90_extra_io%global_kmesh%mesh(2) == pw90_boltzwann%kmesh%mesh(2) .and. &
          pw90_extra_io%global_kmesh%mesh(3) == pw90_boltzwann%kmesh%mesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (pw90_boltzwann%kmesh%spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            pw90_boltzwann%kmesh%mesh(1), 'x', pw90_boltzwann%kmesh%mesh(2), 'x', pw90_boltzwann%kmesh%mesh(3), &
            ' Spacing = ', pw90_boltzwann%kmesh%spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , pw90_boltzwann%kmesh%mesh(1), 'x', pw90_boltzwann%kmesh%mesh(2), 'x', pw90_boltzwann%kmesh%mesh(3), '|'
        endif
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for TDF (eV)                    :', &
        pw90_boltzwann%tdf_energy_step, '|'
      write (stdout, '(1x,a25,5x,a43,4x,a1)') '|  TDF Smearing Function ', &
        trim(w90_readwrite_get_smearing_type(pw90_boltzwann%tdf_smearing%type_index)), '|'
      if (pw90_boltzwann%tdf_smearing%fixed_width > 0.0_dp) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
          '|  TDF fixed Smearing width (eV)             :', pw90_boltzwann%tdf_smearing%fixed_width, '|'
      else
        write (stdout, '(1x,a78)') '|  TDF fixed Smearing width                  :         unsmeared             |'
      endif
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute DOS at same time                  :', pw90_boltzwann%calc_also_dos, '|'
      if (pw90_boltzwann%calc_also_dos .and. print_output%iprint > 2) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', &
          pw90_boltzwann%dos_energy_min, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', &
          pw90_boltzwann%dos_energy_max, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', &
          pw90_boltzwann%dos_energy_step, '|'
        if (pw90_boltzwann%dos_smearing%use_adaptive .eqv. pw90_extra_io%smear%use_adaptive .and. &
            pw90_boltzwann%dos_smearing%adaptive_prefactor == pw90_extra_io%smear%adaptive_prefactor &
            .and. pw90_boltzwann%dos_smearing%adaptive_max_width == pw90_extra_io%smear%adaptive_max_width &
            .and. pw90_boltzwann%dos_smearing%fixed_width == pw90_extra_io%smear%fixed_width .and. &
            pw90_extra_io%smear%type_index == pw90_boltzwann%dos_smearing%type_index) then
          write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
        else
          if (pw90_boltzwann%dos_smearing%use_adaptive) then
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Adaptive width smearing               :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Adaptive smearing factor              :', pw90_boltzwann%dos_smearing%adaptive_prefactor, '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Maximum allowed smearing width        :', pw90_boltzwann%dos_smearing%adaptive_max_width, '|'
          else
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Fixed width smearing                  :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  DOS Smearing width                         :', &
              pw90_boltzwann%dos_smearing%fixed_width, '|'
          endif
          write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
            trim(w90_readwrite_get_smearing_type(pw90_boltzwann%dos_smearing%type_index)), '|'
        endif
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calculation%geninterp .or. print_output%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------Generic Band Interpolation--------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Properties at given k-points      :', pw90_calculation%geninterp, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Calculate band gradients                  :', pw90_geninterp%alsofirstder, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write data into a single file             :', pw90_geninterp%single_file, '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

101 format(20x, a3, 2x, 3F11.6)

  end subroutine w90_postw90_readwrite_write

  !================================================!
  subroutine w90_postw90_readwrite_dealloc(exclude_bands, wannier_data, kmesh_input, kpt_latt, &
                                           dis_manifold, fermi_energy_list, atom_data, eigval, &
                                           kpoint_path, pw90_dos, pw90_berry, proj_input, error, &
                                           comm)
    !================================================!

    use w90_error, only: w90_error_type
    use w90_comms, only: w90comm_type

    implicit none

    type(wannier_data_type), intent(inout) :: wannier_data
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(proj_input_type), intent(inout) :: proj_input
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(atom_data_type), intent(inout) :: atom_data
    type(kpoint_path_type), intent(inout) :: kpoint_path
    type(pw90_dos_mod_type), intent(inout) :: pw90_dos
    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer, allocatable, intent(inout) :: exclude_bands(:)
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    real(kind=dp), allocatable, intent(inout) :: fermi_energy_list(:)
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)

    integer :: ierr

    call w90_readwrite_dealloc(exclude_bands, wannier_data, proj_input, kmesh_input, kpt_latt, &
                               dis_manifold, atom_data, eigval, kpoint_path, error, comm)
    if (allocated(error)) return
    if (allocated(pw90_dos%project)) then
      deallocate (pw90_dos%project, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating dos_project in w90_postw90_readwrite_dealloc', comm)
        return
      endif
    endif
    if (allocated(fermi_energy_list)) then
      deallocate (fermi_energy_list, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating fermi_energy_list in w90_postw90_readwrite_dealloc', comm)
        return
      endif
    endif
    if (allocated(pw90_berry%kubo_freq_list)) then
      deallocate (pw90_berry%kubo_freq_list, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating kubo_freq_list in w90_postw90_readwrite_dealloc', comm)
        return
      endif
    endif
  end subroutine w90_postw90_readwrite_dealloc

  ! extra postw90 memory
  !================================================!
  subroutine w90_postw90_readwrite_mem_estimate(mem_param, mem_bw, dis_manifold, do_boltzwann, &
                                                pw90_boltzwann, spin_decomp, num_wann, stdout)
    !================================================!
    ! note, should only be called from root node
    !================================================!

    implicit none

    type(dis_manifold_type), intent(in) :: dis_manifold
    type(pw90_boltzwann_type) :: pw90_boltzwann

    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    !real(kind=dp), parameter :: size_log = 1.0_dp
    !real(kind=dp), parameter :: size_int = 4.0_dp
    real(kind=dp), parameter :: size_real = 8.0_dp
    real(kind=dp), parameter :: size_cmplx = 16.0_dp
    real(kind=dp), intent(in) :: mem_param
    real(kind=dp), intent(inout) :: mem_bw
    logical, intent(in) :: do_boltzwann, spin_decomp
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
      NumPoints1 = int(floor((pw90_boltzwann%temp_max - pw90_boltzwann%temp_min)/ &
                             pw90_boltzwann%temp_step)) + 1 ! temperature array
      NumPoints2 = int(floor((pw90_boltzwann%mu_max - pw90_boltzwann%mu_min)/ &
                             pw90_boltzwann%mu_step)) + 1  ! mu array
      NumPoints3 = int(floor((dis_manifold%win_max - dis_manifold%win_min &
                              + 2._dp*TDF_exceeding_energy)/ &
                             pw90_boltzwann%tdf_energy_step)) + 1 ! tdfenergyarray
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

      NumPoints1 = int(floor((pw90_boltzwann%dos_energy_max - pw90_boltzwann%dos_energy_min)/ &
                             pw90_boltzwann%dos_energy_step)) + 1!dosnumpoints
      mem_bw = mem_bw + NumPoints1*size_real                         !DOS_EnergyArray
      mem_bw = mem_bw + 6*ndim*NumPoints3*size_real                  !TDF_k
      mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_k
      mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_all
    end if

    if (do_boltzwann) &
      write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'BoltzWann:', &
      (mem_param + mem_bw)/(1024**2), ' Mb'

  end subroutine w90_postw90_readwrite_mem_estimate

  !================================================!
  subroutine parameters_gyro_write_task(task, key, comment, stdout)
    !================================================!
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

end module w90_postw90_readwrite
