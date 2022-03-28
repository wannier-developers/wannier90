
module w90_helper_types

  use w90_constants
  use w90_types
  use w90_wannier90_types

  implicit none

  ! Todo - disentangle, restarts, distribute data, transport
  ! stdout!!!

  type lib_global_type
    type(w90_calculation_type) :: w90_calculation ! separate this?

    type(atom_data_type) :: atom_data
    type(dis_control_type) :: dis_control
    type(dis_manifold_type) :: dis_manifold
    type(dis_spheres_type) :: dis_spheres
    type(kmesh_info_type) :: kmesh_info
    type(kmesh_input_type) :: kmesh_input
    type(kpoint_path_type) :: kpoint_path
    type(output_file_type) :: output_file
    type(print_output_type) :: print_output
    type(proj_input_type) :: proj
    type(proj_input_type) :: proj_input
    type(real_space_ham_type) :: real_space_ham
    type(select_projection_type) :: select_proj
    !type(transport_type) :: tran
    !type(w90_extra_io_type) :: w90_extra_io
    type(w90_system_type) :: w90_system
    type(wann_control_type) :: wann_control
    type(wannier_data_type) :: wannier_data
    type(wann_omega_type) :: wann_omega
    type(ws_region_type) :: ws_region
    type(wvfn_read_type) :: wvfn_read
    !type(w90comm_type) :: comm

    type(timer_list_type) :: timer

    integer, allocatable :: exclude_bands(:)
    integer :: mp_grid(3)
    integer :: num_bands
    integer :: num_kpts
    integer :: num_proj = 0
    integer :: num_wann = -99
    integer :: optimisation = 3
    !integer :: stdout

    real(kind=dp), allocatable :: eigval(:, :)
    real(kind=dp), allocatable :: fermi_energy_list(:)
    real(kind=dp), allocatable :: kpt_latt(:, :)
    real(kind=dp) :: real_lattice(3, 3)
    !real(kind=dp) :: symmetrize_eps

    logical :: eig_found = .false.
    !Projections
    logical :: lhasproj = .false.
    ! RS: symmetry-adapted Wannier functions
    logical :: lsitesymmetry = .false.
    logical :: use_bloch_phases = .false.
    logical :: calc_only_A = .false.
    logical :: gamma_only

    ! added for wannierise
    type(wann_omega_type) :: omega
    type(ham_logical_type) :: ham_logical
    type(sitesym_type) :: sitesym
    complex(kind=dp), allocatable :: ham_k(:, :, :)
    complex(kind=dp), allocatable :: ham_r(:, :, :)
    real(kind=dp), allocatable :: wannier_centres_translated(:, :)
    integer, allocatable :: irvec(:, :)
    integer, allocatable :: shift_vec(:, :)
    integer, allocatable :: ndegen(:)
    integer :: rpt_origin
    integer :: nrpts
    logical :: have_disentangled = .false.
    character(len=128) :: seedname
  end type lib_global_type

  type lib_plot_type
    type(band_plot_type) :: band_plot
    type(wannier_plot_type) :: wann_plot
    type(fermi_surface_plot_type) :: fermi_surface_data
  end type lib_plot_type

  type lib_transport_type
    type(transport_type) :: tran
  end type lib_transport_type

  public:: input_reader, create_kmesh, overlaps, wannierise

contains

  subroutine input_reader(helper, plot, transport, seedname, comm) !, stdout, comm)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    type(lib_plot_type), intent(inout) :: plot
    type(lib_transport_type), intent(inout) :: transport
    !integer, intent(in) :: stdout
    character(len=*), intent(in) :: seedname
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    integer :: stdout
    logical :: cp_pp

    stdout = 6
    call w90_wannier90_readwrite_read(helper%atom_data, plot%band_plot, helper%dis_control, &
                                      helper%dis_spheres, helper%dis_manifold, helper%exclude_bands, &
                                      helper%fermi_energy_list, plot%fermi_surface_data, &
                                      helper%kmesh_input, helper%kmesh_info, helper%kpt_latt, &
                                      helper%output_file, helper%wvfn_read, helper%wann_control, &
                                      helper%wann_omega, helper%proj, helper%proj_input, &
                                      helper%real_space_ham, helper%select_proj, helper%kpoint_path, &
                                      helper%w90_system, transport%tran, helper%print_output, &
                                      helper%wannier_data, plot%wann_plot, io_params, &
                                      helper%ws_region, helper%w90_calculation, helper%eigval, &
                                      helper%real_lattice, physics%bohr, helper%sitesym%symmetrize_eps, &
                                      helper%mp_grid, helper%num_bands, helper%num_kpts, &
                                      helper%num_proj, helper%num_wann, helper%optimisation, &
                                      helper%eig_found, helper%calc_only_A, cp_pp, helper%gamma_only, &
                                      helper%lhasproj, .false., .false., helper%lsitesymmetry, &
                                      helper%use_bloch_phases, seedname, stdout, error, comm)
    if (allocated(error)) then
      write (0, *) 'Error in reader', error%code, error%message
      deallocate (error)
    endif
    helper%seedname = seedname
  end subroutine input_reader

  subroutine checkpoint(helper, label, m_matrix, u_matrix, u_opt, comm)
    !use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    implicit none
    type(lib_global_type), intent(in) :: helper
    character(len=*), intent(in) :: label
    complex(kind=dp), intent(in) :: u_opt(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    !integer, intent(in) :: stdout
    type(w90comm_type), intent(in) :: comm
    !
    integer :: stdout

    if (mpirank(comm) == 0) then
      stdout = 6
      ! e.g. label = 'postwann' after wannierisation
      call w90_wannier90_readwrite_write_chkpt(label, helper%exclude_bands, helper%wannier_data, &
                                               helper%kmesh_info, helper%kpt_latt, helper%num_kpts, helper%dis_manifold, &
                                               helper%num_bands, helper%num_wann, u_matrix, u_opt, m_matrix, helper%mp_grid, &
                                               helper%real_lattice, helper%omega%invariant, helper%have_disentangled, stdout, &
                                               helper%seedname)
    endif
  end subroutine checkpoint

  subroutine create_kmesh(helper, comm)
    use w90_kmesh, only: kmesh_get
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    !integer, intent(in) :: stdout
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_error_type), allocatable :: error
    integer :: stdout

    stdout = 6
    call kmesh_get(helper%kmesh_input, helper%kmesh_info, helper%print_output, helper%kpt_latt, &
                   helper%real_lattice, helper%num_kpts, helper%gamma_only, stdout, helper%timer, &
                   error, comm)
    if (allocated(error)) then
      write (0, *) 'Error in kmesh_get', error%code, error%message
      deallocate (error)
    endif
  end subroutine create_kmesh

  subroutine overlaps(helper, a_matrix, m_matrix, u_matrix, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    use w90_overlap, only: overlap_read
    implicit none
    type(lib_global_type), intent(inout) :: helper
    !integer, intent(in) :: stdout
    type(w90comm_type), intent(in) :: comm
    complex(kind=dp), intent(inout) :: a_matrix(:, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_local(:, :, :, :)
    complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    integer :: stdout
    logical :: cp_pp

    !if disentangle
    !disentanglement = (num_bands > num_wann)
    !allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
    !allocate (m_matrix_orig(num_bands, num_bands, nntot, num_kpts), stat=ierr)
    !allocate (m_matrix_orig_local(num_bands, num_bands, nntot, counts(my_node_id)), stat=ierr)
    !allocate (a_matrix(num_bands, num_wann, num_kpts), stat=ierr)
    !allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
    !else
    !allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
    !allocate (m_matrix_local(num_wann, num_wann, nntot, counts(my_node_id)), stat=ierr)
    cp_pp = .false.
    stdout = 6
    ! should be distributed if MPI
    allocate (m_matrix_local(helper%num_wann, helper%num_wann, helper%kmesh_info%nntot, helper%num_kpts))
    call overlap_read(helper%kmesh_info, helper%select_proj, helper%sitesym, a_matrix, &
                      m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, u_matrix, u_matrix_opt, &
                      helper%num_bands, helper%num_kpts, helper%num_proj, helper%num_wann, &
                      helper%print_output%timing_level, cp_pp, helper%gamma_only, helper%lsitesymmetry, &
                      helper%use_bloch_phases, helper%seedname, stdout, helper%timer, error, comm)
    deallocate (m_matrix_local)
    if (allocated(error)) then
      write (0, *) 'Error in overlaps', error%code, error%message
      deallocate (error)
    endif
  end subroutine overlaps

  subroutine wannierise(helper, plot, transport, m_matrix, u_matrix, u_opt, comm)
    use w90_wannierise, only: wann_main, wann_main_gamma
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    type(lib_plot_type), intent(in) :: plot
    type(lib_transport_type), intent(in) :: transport
    !integer, intent(in) :: stdout
    type(w90comm_type), intent(in) :: comm
    complex(kind=dp), intent(in) :: u_opt(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    integer :: stdout

    stdout = 6
    if (helper%gamma_only) then
      call wann_main_gamma(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                           helper%kmesh_info, helper%kpt_latt, helper%output_file, helper%wann_control, &
                           helper%omega, helper%w90_system, helper%print_output, helper%wannier_data, m_matrix, &
                           u_matrix, u_opt, helper%eigval, helper%real_lattice, helper%mp_grid, &
                           helper%num_bands, helper%num_kpts, helper%num_wann, helper%have_disentangled, &
                           helper%real_space_ham%translate_home_cell, helper%seedname, stdout, &
                           helper%timer, error, comm)
    else
      call wann_main(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                     helper%ham_logical, helper%kmesh_info, helper%kpt_latt, helper%output_file, &
                     helper%real_space_ham, helper%wann_control, helper%omega, helper%sitesym, &
                     helper%w90_system, helper%print_output, helper%wannier_data, helper%ws_region, &
                     helper%w90_calculation, helper%ham_k, helper%ham_r, m_matrix, u_matrix, u_opt, &
                     helper%eigval, helper%real_lattice, helper%wannier_centres_translated, helper%irvec, &
                     helper%mp_grid, helper%ndegen, helper%shift_vec, helper%nrpts, helper%num_bands, &
                     helper%num_kpts, helper%num_proj, helper%num_wann, helper%optimisation, &
                     helper%rpt_origin, plot%band_plot%mode, transport%tran%mode, &
                     helper%have_disentangled, helper%lsitesymmetry, helper%seedname, stdout, &
                     helper%timer, error, comm)
    endif
    if (allocated(error)) then
      write (0, *) 'Error in wannierise', error%code, error%message
      deallocate (error)
    endif
  end subroutine wannierise

  subroutine plot_files(helper, plot, transport, m_matrix, u_matrix, u_matrix_opt, comm)
    use w90_plot, only: plot_main
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper ! inout due to ham_logical
    type(lib_plot_type), intent(in) :: plot
    type(lib_transport_type), intent(in) :: transport
    !integer, intent(in) :: stdout
    type(w90comm_type), intent(in) :: comm
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    integer :: stdout

    stdout = 6
    call plot_main(helper%atom_data, plot%band_plot, helper%dis_manifold, helper%fermi_energy_list, &
                   plot%fermi_surface_data, helper%ham_logical, helper%kmesh_info, helper%kpt_latt, &
                   helper%output_file, helper%wvfn_read, helper%real_space_ham, helper%kpoint_path, &
                   helper%print_output, helper%wannier_data, plot%wann_plot, helper%ws_region, &
                   helper%w90_calculation, helper%ham_k, helper%ham_r, m_matrix, u_matrix, &
                   u_matrix_opt, helper%eigval, helper%real_lattice, helper%wannier_centres_translated, &
                   physics%bohr, helper%irvec, helper%mp_grid, helper%ndegen, helper%shift_vec, helper%nrpts, &
                   helper%num_bands, helper%num_kpts, helper%num_wann, helper%rpt_origin, &
                   transport%tran%mode, helper%have_disentangled, helper%lsitesymmetry, helper%w90_system%spinors, &
                   helper%seedname, stdout, helper%timer, error, comm)
    if (allocated(error)) then
      write (0, *) 'Error in plotting', error%code, error%message
      deallocate (error)
    endif
  end subroutine plot_files

  subroutine print_times(helper)
    use w90_io, only: io_print_timings
    implicit none
    type(lib_global_type), intent(in) :: helper
    integer stdout

    stdout = 6
    call io_print_timings(helper%timer, stdout)
  end subroutine print_times

end module w90_helper_types
