
module w90_helper_types

  use w90_constants
  use w90_types
  use w90_wannier90_types
  use iso_fortran_env, only: error_unit

  implicit none

  ! Todo - disentangle, restarts, distribute data
  ! should we have a lib_wannierise_type?

  type lib_global_type
    type(w90_calculation_type) :: w90_calculation ! separate this?
    ! matrices
    complex(kind=dp), pointer :: u_opt(:, :, :) => null()
    complex(kind=dp), pointer :: a_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: u_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: m_matrix(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_orig(:, :, :, :) => null()

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

  public:: input_reader, create_kmesh, checkpoint, overlaps, wannierise, plot_files, &
           transport, print_times, get_fortran_stdout

contains

  subroutine get_fortran_stdout(output)
    use iso_fortran_env, only: output_unit
    implicit none
    integer, intent(out) :: output

    output = output_unit
  end subroutine get_fortran_stdout

  subroutine input_reader(helper, plot, transport, seedname, output, comm)
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_uppercase, &
      w90_readwrite_clean_infile, w90_readwrite_read_final_alloc
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    implicit none
    type(lib_global_type), intent(inout) :: helper
    type(lib_plot_type), intent(inout) :: plot
    type(lib_transport_type), intent(inout) :: transport
    integer, intent(in) :: output
    character(len=*), intent(in) :: seedname
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    logical :: cp_pp, disentanglement

    call w90_readwrite_in_file(seedname, error, comm)
    if (allocated(error)) then
      write (error_unit, *) 'Error in input file access', error%code, error%message
      deallocate (error)
    else
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
                                        helper%use_bloch_phases, seedname, output, error, comm)
      if (allocated(error)) then
        write (error_unit, *) 'Error in reader', error%code, error%message
        deallocate (error)
      else
        ! For aesthetic purposes, convert some things to uppercase
        call w90_readwrite_uppercase(helper%atom_data, helper%kpoint_path, &
                                     helper%print_output%length_unit)
        disentanglement = (helper%num_bands > helper%num_wann)
        call w90_readwrite_read_final_alloc(disentanglement, helper%dis_manifold, &
                                            helper%wannier_data, helper%num_wann, &
                                            helper%num_bands, helper%num_kpts, error, comm)
        if (allocated(error)) then
          write (error_unit, *) 'Error in read alloc', error%code, error%message
          deallocate (error)
        endif
      endif
      call w90_readwrite_clean_infile(output, seedname, error, comm)
      if (allocated(error)) then
        write (error_unit, *) 'Error in input close', error%code, error%message
        deallocate (error)
      endif
      helper%seedname = seedname
      if (mpirank(comm) /= 0) helper%print_output%iprint = 0
    endif
  end subroutine input_reader

  subroutine checkpoint(helper, label, output, comm)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write_chkpt
    !use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    implicit none
    type(lib_global_type), intent(in) :: helper
    character(len=*), intent(in) :: label
    !complex(kind=dp), intent(in) :: u_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), allocatable :: m_matrix(:, :, :, :)
    integer, intent(in) :: output
    type(w90comm_type), intent(in) :: comm

    if (mpirank(comm) == 0) then
      if ((.not. associated(helper%u_opt)) .or. (.not. associated(helper%u_matrix))) then
        write (error_unit, *) 'Matrices not set for checkpoint read'
        return
      endif
      ! might need to read this in to the helper for a plot/transport restart, but not for pw90
      allocate (m_matrix(helper%num_wann, helper%num_wann, helper%kmesh_info%nntot, helper%num_kpts))
      ! e.g. label = 'postwann' after wannierisation
      call w90_wannier90_readwrite_write_chkpt(label, helper%exclude_bands, helper%wannier_data, &
                                               helper%kmesh_info, helper%kpt_latt, helper%num_kpts, helper%dis_manifold, &
                                               helper%num_bands, helper%num_wann, helper%u_matrix, &
                                               helper%u_opt, m_matrix, helper%mp_grid, &
                                               helper%real_lattice, helper%omega%invariant, helper%have_disentangled, output, &
                                               helper%seedname)
      deallocate (m_matrix)
    endif
  end subroutine checkpoint

  subroutine create_kmesh(helper, output, comm)
    use w90_kmesh, only: kmesh_get
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    integer, intent(in) :: output
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_error_type), allocatable :: error

    call kmesh_get(helper%kmesh_input, helper%kmesh_info, helper%print_output, helper%kpt_latt, &
                   helper%real_lattice, helper%num_kpts, helper%gamma_only, output, helper%timer, &
                   error, comm)
    if (allocated(error)) then
      write (error_unit, *) 'Error in kmesh_get', error%code, error%message
      deallocate (error)
    endif
  end subroutine create_kmesh

  subroutine overlaps(helper, output, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    use w90_overlap, only: overlap_read
    implicit none
    type(lib_global_type), intent(inout) :: helper
    integer, intent(in) :: output
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(inout) :: a_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: m_orig(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_local(:, :, :, :)
    complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    logical :: cp_pp

    if ((.not. associated(helper%m_matrix)) .or. (.not. associated(helper%m_orig)) .or. &
        (.not. associated(helper%a_matrix)) .or. (.not. associated(helper%u_matrix))) then
      write (error_unit, *) 'Matrices not set for overlap call'
      return
    endif
    if (helper%num_bands > helper%num_wann) then
      ! disentnglement
      !allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      !allocate (m_matrix_orig(helper%num_bands, helper%num_bands, helper%kmesh_info%nntot, &
      !                        helper%num_kpts))
      !allocate (m_matrix_orig_local(helper%num_bands, helper%num_bands, nntot, counts(my_node_id)), stat=ierr)
      allocate (m_matrix_orig_local(helper%num_bands, helper%num_bands, helper%kmesh_info%nntot, &
                                    helper%num_kpts))
      !allocate (a_matrix(num_bands, num_wann, num_kpts), stat=ierr)
      !allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
      !else
      !allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
      !allocate (m_matrix_local(num_wann, num_wann, nntot, counts(my_node_id)), stat=ierr)
    endif
    cp_pp = .false.
    ! should be distributed if MPI
    allocate (m_matrix_local(helper%num_wann, helper%num_wann, helper%kmesh_info%nntot, helper%num_kpts))
    call overlap_read(helper%kmesh_info, helper%select_proj, helper%sitesym, helper%a_matrix, &
                      helper%m_matrix, m_matrix_local, helper%m_orig, m_matrix_orig_local, helper%u_matrix, u_matrix_opt, &
                      helper%num_bands, helper%num_kpts, helper%num_proj, helper%num_wann, &
                      helper%print_output%timing_level, cp_pp, helper%gamma_only, helper%lsitesymmetry, &
                      helper%use_bloch_phases, helper%seedname, output, helper%timer, error, comm)
    deallocate (m_matrix_local)
    if (helper%num_bands > helper%num_wann) then
      deallocate (m_matrix_orig_local)
      !deallocate (m_matrix_orig)
    endif
    if (allocated(u_matrix_opt)) deallocate (u_matrix_opt)
    if (allocated(error)) then
      write (error_unit, *) 'Error in overlaps', error%code, error%message
      deallocate (error)
    endif
  end subroutine overlaps

  subroutine disentangle(helper, output, comm)
    use w90_disentangle, only: dis_main
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    !type(lib_plot_type), intent(in) :: plot
    !type(lib_transport_type), intent(in) :: transport
    integer, intent(in) :: output
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(inout) :: u_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: a_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: m_matrix_local(:, :, :, :)
    !complex(kind=dp), intent(inout) :: m_orig(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)

    if ((.not. associated(helper%m_matrix)) .or. (.not. associated(helper%m_orig)) .or. &
        (.not. associated(helper%a_matrix)) .or. (.not. associated(helper%u_matrix)) .or. &
        (.not. associated(helper%u_opt))) then
      write (error_unit, *) 'Matrices not set for disentangle call'
      return
    endif
    !allocate (m_matrix_orig(helper%num_bands, helper%num_bands, helper%kmesh_info%nntot, &
    !                        helper%num_kpts))
    ! MPI issue
    allocate (m_matrix_orig_local(helper%num_bands, helper%num_bands, helper%kmesh_info%nntot, &
                                  helper%num_kpts))
    m_matrix_orig_local(1:helper%num_bands, 1:helper%num_bands, 1:helper%kmesh_info%nntot, &
                        1:helper%num_kpts) = helper%m_orig(1:helper%num_bands, 1:helper%num_bands, &
                                                           1:helper%kmesh_info%nntot, 1:helper%num_kpts)
    ! MPI issue
    allocate (m_matrix_local(helper%num_wann, helper%num_wann, helper%kmesh_info%nntot, &
                             helper%num_kpts))
    call dis_main(helper%dis_control, helper%dis_spheres, helper%dis_manifold, helper%kmesh_info, &
                  helper%kpt_latt, helper%sitesym, helper%print_output, helper%a_matrix, helper%m_matrix, &
                  m_matrix_local, helper%m_orig, m_matrix_orig_local, helper%u_matrix, helper%u_opt, &
                  helper%eigval, helper%real_lattice, helper%omega%invariant, helper%num_bands, &
                  helper%num_kpts, helper%num_wann, helper%optimisation, helper%gamma_only, &
                  helper%lsitesymmetry, output, helper%timer, error, comm)
    helper%have_disentangled = .true.
    if (allocated(m_matrix_local)) deallocate (m_matrix_local)
    if (allocated(m_matrix_orig_local)) deallocate (m_matrix_orig_local)
    !if (allocated(m_matrix_orig)) deallocate (m_matrix_orig)
    if (allocated(error)) then
      write (error_unit, *) 'Error in disentangle', error%code, error%message
      deallocate (error)
    endif
  end subroutine disentangle

  subroutine wannierise(helper, plot, transport, output, comm)
    use w90_wannierise, only: wann_main, wann_main_gamma
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    type(lib_plot_type), intent(in) :: plot
    type(lib_transport_type), intent(in) :: transport
    integer, intent(in) :: output
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(in) :: u_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error

    if ((.not. associated(helper%m_matrix)) .or. (.not. associated(helper%u_opt)) .or. &
        (.not. associated(helper%u_matrix))) then
      write (error_unit, *) 'Matrices not set for wannierise call'
      return
    endif
    if (helper%gamma_only) then
      call wann_main_gamma(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                           helper%kmesh_info, helper%kpt_latt, helper%output_file, helper%wann_control, &
                           helper%omega, helper%w90_system, helper%print_output, helper%wannier_data, helper%m_matrix, &
                           helper%u_matrix, helper%u_opt, helper%eigval, helper%real_lattice, helper%mp_grid, &
                           helper%num_bands, helper%num_kpts, helper%num_wann, helper%have_disentangled, &
                           helper%real_space_ham%translate_home_cell, helper%seedname, output, &
                           helper%timer, error, comm)
    else
      call wann_main(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                     helper%ham_logical, helper%kmesh_info, helper%kpt_latt, helper%output_file, &
                     helper%real_space_ham, helper%wann_control, helper%omega, helper%sitesym, &
                     helper%w90_system, helper%print_output, helper%wannier_data, helper%ws_region, &
                     helper%w90_calculation, helper%ham_k, helper%ham_r, helper%m_matrix, helper%u_matrix, helper%u_opt, &
                     helper%eigval, helper%real_lattice, helper%wannier_centres_translated, helper%irvec, &
                     helper%mp_grid, helper%ndegen, helper%shift_vec, helper%nrpts, helper%num_bands, &
                     helper%num_kpts, helper%num_proj, helper%num_wann, helper%optimisation, &
                     helper%rpt_origin, plot%band_plot%mode, transport%tran%mode, &
                     helper%have_disentangled, helper%lsitesymmetry, helper%seedname, output, &
                     helper%timer, error, comm)
    endif
    if (allocated(error)) then
      write (error_unit, *) 'Error in wannierise', error%code, error%message
      deallocate (error)
    endif
  end subroutine wannierise

  subroutine plot_files(helper, plot, transport, output, comm)
    use w90_plot, only: plot_main
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper ! inout due to ham_logical
    type(lib_plot_type), intent(in) :: plot
    type(lib_transport_type), intent(in) :: transport
    integer, intent(in) :: output
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    !complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error

    call plot_main(helper%atom_data, plot%band_plot, helper%dis_manifold, helper%fermi_energy_list, &
                   plot%fermi_surface_data, helper%ham_logical, helper%kmesh_info, helper%kpt_latt, &
                   helper%output_file, helper%wvfn_read, helper%real_space_ham, helper%kpoint_path, &
                   helper%print_output, helper%wannier_data, plot%wann_plot, helper%ws_region, &
                   helper%w90_calculation, helper%ham_k, helper%ham_r, helper%m_matrix, helper%u_matrix, &
                   helper%u_opt, helper%eigval, helper%real_lattice, helper%wannier_centres_translated, &
                   physics%bohr, helper%irvec, helper%mp_grid, helper%ndegen, helper%shift_vec, helper%nrpts, &
                   helper%num_bands, helper%num_kpts, helper%num_wann, helper%rpt_origin, &
                   transport%tran%mode, helper%have_disentangled, helper%lsitesymmetry, helper%w90_system%spinors, &
                   helper%seedname, output, helper%timer, error, comm)
    if (allocated(error)) then
      write (error_unit, *) 'Error in plotting', error%code, error%message
      deallocate (error)
    endif
  end subroutine plot_files

  subroutine transport(helper, plot, tran, output, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    use w90_transport, only: tran_main
    implicit none
    type(lib_global_type), intent(inout) :: helper ! because of ham_logical
    type(lib_plot_type), intent(in) :: plot
    type(lib_transport_type), intent(inout) :: tran
    integer, intent(in) :: output
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(in) :: u_opt(:, :, :)
    !complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error

    ! should these tests be done outside?
    if (mpirank(comm) == 0) then
      if (helper%w90_calculation%transport) then
        call tran_main(helper%atom_data, helper%dis_manifold, helper%fermi_energy_list, &
                       helper%ham_logical, helper%kpt_latt, helper%output_file, helper%real_space_ham, &
                       tran%tran, helper%print_output, helper%wannier_data, helper%ws_region, &
                       helper%w90_calculation, helper%ham_k, helper%ham_r, helper%u_matrix, helper%u_opt, &
                       helper%eigval, helper%real_lattice, helper%wannier_centres_translated, helper%irvec, &
                       helper%mp_grid, helper%ndegen, helper%shift_vec, helper%nrpts, helper%num_bands, &
                       helper%num_kpts, helper%num_wann, helper%rpt_origin, plot%band_plot%mode, &
                       helper%have_disentangled, helper%lsitesymmetry, helper%seedname, output, &
                       helper%timer, error, comm)
        if (allocated(error)) then
          write (error_unit, *) 'Error in transport', error%code, error%message
          deallocate (error)
        endif
      endif
    endif
  end subroutine transport

  subroutine print_times(helper, output)
    use w90_io, only: io_print_timings
    implicit none
    type(lib_global_type), intent(in) :: helper
    integer, intent(in) :: output

    call io_print_timings(helper%timer, output)
  end subroutine print_times

  subroutine set_a_matrix(helper, a_matrix)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: a_matrix(:, :, :)

    helper%a_matrix => a_matrix
  end subroutine set_a_matrix

  subroutine set_m_matrix(helper, m_matrix)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: m_matrix(:, :, :, :)

    helper%m_matrix => m_matrix
  end subroutine set_m_matrix

  subroutine set_m_orig(helper, m_orig)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: m_orig(:, :, :, :)

    helper%m_orig => m_orig
  end subroutine set_m_orig

  subroutine set_u_matrix(helper, u_matrix)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: u_matrix(:, :, :)

    helper%u_matrix => u_matrix
  end subroutine set_u_matrix

  subroutine set_u_opt(helper, u_opt)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: u_opt(:, :, :)

    helper%u_opt => u_opt
  end subroutine set_u_opt

end module w90_helper_types
