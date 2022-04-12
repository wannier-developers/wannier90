
module w90_lib_all

  use w90_constants
  use w90_types
  use w90_wannier90_types
  use w90_helper_types
  use w90_postw90_types

  ! Todo - initialisation issues that we had to fix
  ! Todo - read_chkpt and allocatable
  ! Todo - num_valence_bands init

  implicit none

  type lib_postw90_type
    type(pw90_calculation_type) :: calculation
    type(pw90_berry_mod_type) :: berry
    type(pw90_boltzwann_type) :: boltzwann
    type(pw90_dos_mod_type) :: dos
    type(pw90_geninterp_mod_type) :: geninterp
    type(pw90_gyrotropic_type) :: gyrotropic
    type(pw90_kpath_mod_type) :: kpath
    type(pw90_kslice_mod_type) :: kslice
    type(pw90_band_deriv_degen_type) :: band_deriv_degen
    type(pw90_oper_read_type) :: oper_read
    type(pw90_spin_mod_type) :: spin
    type(pw90_spin_hall_type) :: spin_hall
    real(kind=dp) :: scissors_shift
    logical :: effective_model

    character(len=20) :: checkpoint

    ! for dos - make them local to dos call?
    type(kpoint_dist_type) :: kpt_dist
    type(ws_distance_type) :: ws_distance
    type(wigner_seitz_type) :: ws_vec
  end type lib_postw90_type

  public :: read_all_input, read_checkpoint, calc_dos

contains

  subroutine read_all_input(wann90, pw90, plot, transport, seedname, output, comm)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_readall, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    use w90_postw90_readwrite, only: w90_postw90_readwrite_readall, pw90_extra_io_type
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_uppercase, &
      w90_readwrite_clean_infile, w90_readwrite_read_final_alloc

    implicit none
    type(lib_global_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    type(lib_plot_type), intent(inout) :: plot
    type(lib_transport_type), intent(inout) :: transport
    integer, intent(in) :: output
    character(len=*), intent(in) :: seedname
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    type(pw90_extra_io_type) :: pw90_params
    logical :: cp_pp
    logical :: disentanglement

    call w90_readwrite_in_file(seedname, error, comm)
    if (allocated(error)) then
      write (0, *) 'Error in input file access', error%code, error%message
      deallocate (error)
    else
      call w90_wannier90_readwrite_readall(wann90%atom_data, plot%band_plot, wann90%dis_control, &
                                           wann90%dis_spheres, wann90%dis_manifold, &
                                           wann90%exclude_bands, wann90%fermi_energy_list, &
                                           plot%fermi_surface_data, wann90%kmesh_input, &
                                           wann90%kmesh_info, wann90%kpt_latt, wann90%output_file, &
                                           wann90%wvfn_read, wann90%wann_control, &
                                           wann90%wann_omega, wann90%proj, wann90%proj_input, &
                                           wann90%real_space_ham, wann90%select_proj, &
                                           wann90%kpoint_path, wann90%w90_system, transport%tran, &
                                           wann90%print_output, wann90%wannier_data, plot%wann_plot, &
                                           io_params, wann90%ws_region, wann90%w90_calculation, &
                                           wann90%eigval, wann90%real_lattice, physics%bohr, &
                                           wann90%sitesym%symmetrize_eps, wann90%mp_grid, &
                                           wann90%num_bands, wann90%num_kpts, wann90%num_proj, &
                                           wann90%num_wann, wann90%optimisation, wann90%eig_found, &
                                           wann90%calc_only_A, cp_pp, wann90%gamma_only, &
                                           wann90%lhasproj, wann90%lsitesymmetry, &
                                           wann90%use_bloch_phases, seedname, output, error, comm)
      wann90%seedname = seedname
      if (mpirank(comm) /= 0) wann90%print_output%iprint = 0
      if (allocated(error)) then
        write (0, *) 'Error in wannier90 read', error%code, error%message
        deallocate (error)
      else
        call w90_postw90_readwrite_readall(wann90%w90_system, wann90%dis_manifold, &
                                           wann90%fermi_energy_list, wann90%num_bands, &
                                           wann90%num_wann, wann90%eigval, wann90%real_lattice, &
                                           wann90%kpoint_path, pw90%calculation, pw90%oper_read, &
                                           pw90%scissors_shift, pw90%effective_model, pw90%spin, &
                                           pw90%band_deriv_degen, pw90%kpath, pw90%kslice, &
                                           pw90%dos, pw90%berry, pw90%spin_hall, pw90%gyrotropic, &
                                           pw90%geninterp, pw90%boltzwann, pw90_params, error, comm)
        if (allocated(error)) then
          write (0, *) 'Error in postw90 read', error%code, error%message
          deallocate (error)
        else
          ! For aesthetic purposes, convert some things to uppercase
          call w90_readwrite_uppercase(wann90%atom_data, wann90%kpoint_path, &
                                       wann90%print_output%length_unit)
          disentanglement = (wann90%num_bands > wann90%num_wann)
          call w90_readwrite_read_final_alloc(disentanglement, wann90%dis_manifold, &
                                              wann90%wannier_data, wann90%num_wann, &
                                              wann90%num_bands, wann90%num_kpts, error, comm)
          if (allocated(error)) then
            write (0, *) 'Error in read alloc', error%code, error%message
            deallocate (error)
          endif
        endif
        call w90_readwrite_clean_infile(output, seedname, error, comm)
        if (allocated(error)) then
          write (0, *) 'Error in input close', error%code, error%message
          deallocate (error)
        endif
      endif
    endif
  end subroutine read_all_input

  subroutine read_checkpoint(wann90, pw90, m_matrix, u_matrix, u_opt, output, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    use w90_readwrite, only: w90_readwrite_read_chkpt_header, w90_readwrite_read_chkpt_matrices

    implicit none
    type(lib_global_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: output
    complex(kind=dp), intent(inout) :: u_opt(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_error_type), allocatable :: error
    integer :: chk_unit
    integer :: num_exclude_bands

    ! Todo - fix ispostw90 flag
    num_exclude_bands = 0
    if (allocated(wann90%exclude_bands)) num_exclude_bands = size(wann90%exclude_bands)
    call w90_readwrite_read_chkpt_header(wann90%exclude_bands, wann90%kmesh_info, wann90%kpt_latt, &
                                         wann90%real_lattice, wann90%mp_grid, wann90%num_bands, &
                                         num_exclude_bands, wann90%num_kpts, wann90%num_wann, &
                                         pw90%checkpoint, wann90%have_disentangled, .true., &
                                         wann90%seedname, chk_unit, output, error, comm)
    if (allocated(error)) then
      write (0, *) 'Error in reading checkpoint header', error%code, error%message
      deallocate (error)
    else
      call w90_readwrite_read_chkpt_matrices(wann90%dis_manifold, wann90%kmesh_info, &
                                             wann90%wannier_data, m_matrix, u_matrix, u_opt, &
                                             wann90%omega%invariant, wann90%num_bands, &
                                             wann90%num_kpts, wann90%num_wann, &
                                             wann90%have_disentangled, wann90%seedname, chk_unit, &
                                             output, error, comm)
      if (allocated(error)) then
        write (0, *) 'Error in reading checkpoint matrices', error%code, error%message
        deallocate (error)
      endif
    endif
  end subroutine read_checkpoint

  subroutine calc_dos(wann90, pw90, u_matrix, u_opt, output, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    use w90_dos, only: dos_main
    use w90_postw90_common, only: pw90common_wanint_setup

    implicit none
    type(lib_global_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    !type(lib_plot_type), intent(inout) :: plot
    !type(lib_transport_type), intent(inout) :: transport
    integer, intent(in) :: output
    complex(kind=dp), intent(inout) :: u_opt(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !character(len=*), intent(in) :: seedname
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_error_type), allocatable :: error
    !type(w90_extra_io_type) :: io_params
    !type(pw90_extra_io_type) :: pw90_params
    !logical :: cp_pp
    !logical :: disentanglement
    !integer :: num_exclude_bands
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)
    complex(kind=dp), allocatable :: v_matrix(:, :, :)
    integer :: i, j, m, loop_kpt, ierr

    ! put this in a separate setup? (since may be coming from wannierise rather than checkpoint
    call pw90common_wanint_setup(wann90%num_wann, wann90%print_output, wann90%real_lattice, &
                                 wann90%mp_grid, pw90%effective_model, pw90%ws_vec, output, &
                                 wann90%seedname, wann90%timer, error, comm)
    if (pw90%calculation%dos .and. index(pw90%dos%task, 'dos_plot') > 0) then
      ! build v_matrix, shouldn't really be here
      allocate (v_matrix(wann90%num_bands, wann90%num_wann, wann90%num_kpts), stat=ierr)
      if (ierr /= 0) then
        write (0, *) 'Error allocating v_matrix in calc_dos'
      else
        ! u_matrix and u_opt are stored on root only
        if (.not. wann90%have_disentangled) then
          v_matrix(1:wann90%num_wann, :, :) = u_matrix(1:wann90%num_wann, :, :)
        else
          v_matrix = cmplx_0
          do loop_kpt = 1, wann90%num_kpts
            do j = 1, wann90%num_wann
              do m = 1, wann90%dis_manifold%ndimwin(loop_kpt)
                do i = 1, wann90%num_wann
                  v_matrix(m, j, loop_kpt) = v_matrix(m, j, loop_kpt) &
                                             + u_opt(m, i, loop_kpt)*u_matrix(i, j, loop_kpt)
                enddo
              enddo
            enddo
          enddo
        endif
        call dos_main(pw90%berry, wann90%dis_manifold, pw90%dos, pw90%kpt_dist, wann90%kpt_latt, &
                      pw90%oper_read, pw90%band_deriv_degen, pw90%spin, wann90%ws_region, &
                      wann90%w90_system, wann90%print_output, wann90%wannier_data, pw90%ws_distance, &
                      pw90%ws_vec, HH_R, SS_R, u_matrix, v_matrix, wann90%eigval, wann90%real_lattice, &
                      pw90%scissors_shift, wann90%mp_grid, wann90%num_bands, wann90%num_kpts, &
                      wann90%num_wann, pw90%effective_model, wann90%have_disentangled, &
                      pw90%calculation%spin_decomp, wann90%seedname, output, wann90%timer, &
                      error, comm)
        if (allocated(HH_R)) deallocate (HH_R)
        if (allocated(SS_R)) deallocate (SS_R)
        if (allocated(v_matrix)) deallocate (v_matrix)
        if (allocated(error)) then
          write (0, *) 'Error in dos', error%code, error%message
          deallocate (error)
        endif
      endif
    else
      write (output, *) ' No dos calculation requested'
    endif
  end subroutine calc_dos

end module w90_lib_all
