
module w90_lib_all

  use w90_constants
  use w90_types
  use w90_wannier90_types
  use w90_library
  use w90_postw90_types

  ! Todo - initialisation issues that we had to fix
  ! Todo - num_valence_bands init
  ! Todo - need to test MPI (data dist routines...?)
  ! have AA_R etc as numpy matrices that are passed in so calculate once separately?
  ! or have 2 versions, one which takes those matrices and a hihger level wrapper that
  ! creates them and then calls the other, so that for single one-off calls you don't need
  ! to keep them in the SCF or other code. BUT m_matrix is part of the checkpoint isn't it
  ! so it shouldn't be re-read!

  implicit none

  type lib_postw90_type
    complex(kind=dp), pointer :: v_matrix(:, :, :)

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
    real(kind=dp) :: scissors_shift = 0.0_dp
    logical :: effective_model = .false.

    character(len=20) :: checkpoint

    ! for dos - make them local to dos call?
    type(kpoint_dist_type) :: kpt_dist
    type(ws_distance_type) :: ws_distance
    type(wigner_seitz_type) :: ws_vec

    ! put eigenvalues here for the moment
    !real(kind=dp), allocatable :: eigval(:, :)
    logical :: eig_found = .false.
  end type lib_postw90_type

  public :: read_checkpoint, calc_dos, boltzwann, gyrotropic, berry, kpath, kslice, &
            read_all_input_has_eigs, read_all_input_and_eigs
  private :: read_all_input

contains

  subroutine read_all_input_and_eigs(wann90, w90only, pw90, seedname, istdout, istderr, ierr)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_postw90_readwrite, only: w90_postw90_readwrite_readall, pw90_extra_io_type
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_clean_infile, &
      w90_readwrite_read_final_alloc, w90_readwrite_read_eigvals

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_wannier_type), intent(inout) :: w90only
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    character(len=*), intent(in) :: seedname
    integer, intent(out) :: ierr
    real(kind=dp) :: eigval(1, 1)

    call read_all_input(wann90, w90only, pw90, eigval, .false., seedname, istdout, istderr, ierr)
  end subroutine read_all_input_and_eigs

  subroutine read_all_input_has_eigs(wann90, w90only, pw90, eigval, seedname, istdout, istderr, &
                                     ierr)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_postw90_readwrite, only: w90_postw90_readwrite_readall, pw90_extra_io_type
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_clean_infile, &
      w90_readwrite_read_final_alloc, w90_readwrite_read_eigvals

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_wannier_type), intent(inout) :: w90only
    type(lib_postw90_type), intent(inout) :: pw90
    real(kind=dp), intent(inout) :: eigval(:, :)
    integer, intent(in) :: istdout, istderr
    character(len=*), intent(in) :: seedname
    integer, intent(out) :: ierr

    call read_all_input(wann90, w90only, pw90, eigval, .true., seedname, istdout, istderr, ierr)
  end subroutine read_all_input_has_eigs

  subroutine read_all_input(wann90, w90only, pw90, eigval, eig_ok, seedname, istdout, istderr, ierr)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_postw90_readwrite, only: w90_postw90_readwrite_readall, pw90_extra_io_type
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_clean_infile, &
      w90_readwrite_read_final_alloc, w90_readwrite_read_eigvals

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_wannier_type), intent(inout) :: w90only
    type(lib_postw90_type), intent(inout) :: pw90
    real(kind=dp), intent(in) :: eigval(:, :)
    integer, intent(in) :: istdout, istderr
    character(len=*), intent(in) :: seedname
    integer, intent(out) :: ierr
    logical, intent(in) :: eig_ok
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    real(kind=dp), pointer :: read_eigs(:, :)
    type(pw90_extra_io_type) :: pw90_params
    logical :: cp_pp
    logical :: disentanglement

    ierr = 0
    call w90_readwrite_in_file(wann90%settings, seedname, error, wann90%comm)
    if (allocated(error)) then
      write (istderr, *) 'Error in input file access', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    else
      call w90_wannier90_readwrite_read(wann90%settings, wann90%atom_data, w90only%band_plot, &
                                        w90only%dis_control, w90only%dis_spheres, &
                                        wann90%dis_manifold, wann90%exclude_bands, &
                                        wann90%fermi_energy_list, w90only%fermi_surface_data, &
                                        wann90%kmesh_input, wann90%kmesh_info, wann90%kpt_latt, &
                                        w90only%output_file, w90only%wvfn_read, &
                                        w90only%wann_control, w90only%proj, w90only%proj_input, &
                                        w90only%real_space_ham, w90only%select_proj, &
                                        wann90%kpoint_path, wann90%w90_system, w90only%tran, &
                                        wann90%print_output, w90only%wann_plot, &
                                        io_params, wann90%ws_region, w90only%w90_calculation, &
                                        wann90%real_lattice, physics%bohr, &
                                        w90only%sitesym%symmetrize_eps, wann90%mp_grid, &
                                        wann90%num_bands, wann90%num_kpts, w90only%num_proj, &
                                        wann90%num_wann, w90only%optimisation, &
                                        w90only%calc_only_A, cp_pp, wann90%gamma_only, &
                                        w90only%lhasproj, w90only%lsitesymmetry, &
                                        w90only%use_bloch_phases, seedname, istdout, error, wann90%comm)
      wann90%seedname = seedname
      if (mpirank(wann90%comm) /= 0) wann90%print_output%iprint = 0
      if (allocated(error)) then
        write (istderr, *) 'Error in wannier90 read', error%code, error%message
        ierr = sign(1, error%code)
        deallocate (error)
      else
        disentanglement = (wann90%num_bands > wann90%num_wann)
        ! If the user is setting values from outside then the eigvals should be associated,
        ! if they're not then we assume its a 'driver' program that needs to read then
        if (eig_ok) then
          call set_eigval(wann90, eigval)
          pw90%eig_found = .true.
        else
          allocate (read_eigs(wann90%num_bands, wann90%num_kpts))
          call w90_readwrite_read_eigvals(pw90%eig_found, read_eigs, wann90%num_bands, &
                                          wann90%num_kpts, istdout, seedname, error, wann90%comm)
          if (.not. allocated(error)) then
            call set_eigval(wann90, read_eigs)
          endif
        endif
        if (allocated(error)) then
          write (istderr, *) 'Error in wannier90 eigenvalues', error%code, error%message
          ierr = sign(1, error%code)
          deallocate (error)
        else
          call w90_postw90_readwrite_readall(wann90%settings, wann90%w90_system, &
                                             wann90%dis_manifold, wann90%fermi_energy_list, &
                                             wann90%num_bands, wann90%num_wann, wann90%eigval, &
                                             wann90%real_lattice, wann90%kpoint_path, &
                                             pw90%calculation, pw90%oper_read, &
                                             pw90%scissors_shift, pw90%effective_model, pw90%spin, &
                                             pw90%band_deriv_degen, pw90%kpath, pw90%kslice, &
                                             pw90%dos, pw90%berry, pw90%spin_hall, &
                                             pw90%gyrotropic, pw90%geninterp, pw90%boltzwann, &
                                             pw90_params, error, wann90%comm)
          if (allocated(error)) then
            write (istderr, *) 'Error in postw90 read', error%code, error%message
            ierr = sign(1, error%code)
            deallocate (error)
          else
            call w90_readwrite_read_final_alloc(disentanglement, wann90%dis_manifold, &
                                                wann90%wannier_data, wann90%num_wann, &
                                                wann90%num_bands, wann90%num_kpts, error, wann90%comm)
            if (allocated(error)) then
              write (istderr, *) 'Error in read alloc', error%code, error%message
              ierr = sign(1, error%code)
              deallocate (error)
            endif
          endif
          call w90_readwrite_clean_infile(wann90%settings, istdout, seedname, error, wann90%comm)
          if (allocated(error)) then
            write (istderr, *) 'Error in input close', error%code, error%message
            ierr = sign(1, error%code)
            deallocate (error)
          endif
        endif
      endif
    endif
  end subroutine read_all_input

  subroutine read_checkpoint(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_readwrite, only: w90_readwrite_read_chkpt_header, w90_readwrite_read_chkpt_matrices
    use w90_postw90_common, only: pw90common_wanint_setup

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: m_matrix(:, :, :, :)
    real(kind=dp) :: omega_invariant ! assuming postw90 in doing this
    integer :: chk_unit
    integer :: num_exclude_bands

    ierr = 0
    ! Todo - fix ispostw90 flag
    num_exclude_bands = 0
    if (allocated(wann90%exclude_bands)) num_exclude_bands = size(wann90%exclude_bands)
    call w90_readwrite_read_chkpt_header(wann90%exclude_bands, wann90%kmesh_info, wann90%kpt_latt, &
                                         wann90%real_lattice, wann90%mp_grid, wann90%num_bands, &
                                         num_exclude_bands, wann90%num_kpts, wann90%num_wann, &
                                         pw90%checkpoint, wann90%have_disentangled, .true., &
                                         wann90%seedname, chk_unit, istdout, error, wann90%comm)
    if (allocated(error)) then
      write (istderr, *) 'Error in reading checkpoint header', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    else
      allocate (m_matrix(wann90%num_wann, wann90%num_wann, wann90%kmesh_info%nntot, &
                         wann90%num_kpts))
      call w90_readwrite_read_chkpt_matrices(wann90%dis_manifold, wann90%kmesh_info, &
                                             wann90%wannier_data, m_matrix, wann90%u_matrix, &
                                             wann90%u_opt, omega_invariant, wann90%num_bands, &
                                             wann90%num_kpts, wann90%num_wann, &
                                             wann90%have_disentangled, wann90%seedname, chk_unit, &
                                             istdout, error, wann90%comm)
      deallocate (m_matrix)
      if (allocated(error)) then
        write (istderr, *) 'Error in reading checkpoint matrices', error%code, error%message
        ierr = sign(1, error%code)
        deallocate (error)
      endif
    endif
  end subroutine read_checkpoint

  subroutine pw_setup(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_postw90_common, only: pw90common_wanint_setup

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    !
    type(w90_error_type), allocatable :: error

    ierr = 0
    call pw90common_wanint_setup(wann90%num_wann, wann90%print_output, wann90%real_lattice, &
                                 wann90%mp_grid, pw90%effective_model, wann90%ws_region, &
                                 pw90%ws_vec, istdout, wann90%seedname, wann90%timer, error, &
                                 wann90%comm)
    if (allocated(error)) then
      write (istderr, *) 'Error in post setup', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine pw_setup

  subroutine calc_v_matrix(wann90, pw90, v_matrix)
    !use w90_error_base, only: w90_error_type
    !use w90_comms, only: w90_comm_type, mpirank

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    !integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout), target :: v_matrix(:, :, :)
    !
    integer :: i, j, m, loop_kpt

    !allocate (v_matrix(wann90%num_bands, wann90%num_wann, wann90%num_kpts), stat=ierr)
    ! u_matrix and u_opt are stored on root only
    if (.not. wann90%have_disentangled) then
      v_matrix(1:wann90%num_wann, :, :) = wann90%u_matrix(1:wann90%num_wann, :, :)
    else
      !this should be initialised by the caller really
      v_matrix(1:wann90%num_bands, 1:wann90%num_wann, 1:wann90%num_kpts) = cmplx_0
      do loop_kpt = 1, wann90%num_kpts
        do j = 1, wann90%num_wann
          do m = 1, wann90%dis_manifold%ndimwin(loop_kpt)
            do i = 1, wann90%num_wann
              v_matrix(m, j, loop_kpt) = v_matrix(m, j, loop_kpt) &
                                         + wann90%u_opt(m, i, loop_kpt)*wann90%u_matrix(i, j, loop_kpt)
            enddo
          enddo
        enddo
      enddo
    endif
    pw90%v_matrix => v_matrix
  end subroutine calc_v_matrix

  subroutine calc_dos(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_dos, only: dos_main

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    !type(lib_plot_type), intent(inout) :: plot
    !type(lib_transport_type), intent(inout) :: transport
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    !character(len=*), intent(in) :: seedname
    integer, intent(out) :: ierr
    !
    type(w90_error_type), allocatable :: error
    !type(w90_extra_io_type) :: io_params
    !type(pw90_extra_io_type) :: pw90_params
    !logical :: cp_pp
    !logical :: disentanglement
    !integer :: num_exclude_bands
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)

    ierr = 0
    if (pw90%calculation%dos .and. index(pw90%dos%task, 'dos_plot') > 0) then
      call dos_main(pw90%berry, wann90%dis_manifold, pw90%dos, pw90%kpt_dist, wann90%kpt_latt, &
                    pw90%oper_read, pw90%band_deriv_degen, pw90%spin, wann90%ws_region, &
                    wann90%w90_system, wann90%print_output, wann90%wannier_data, pw90%ws_distance, &
                    pw90%ws_vec, HH_R, SS_R, wann90%u_matrix, pw90%v_matrix, wann90%eigval, &
                    wann90%real_lattice, pw90%scissors_shift, wann90%mp_grid, wann90%num_bands, &
                    wann90%num_kpts, wann90%num_wann, pw90%effective_model, &
                    wann90%have_disentangled, pw90%calculation%spin_decomp, wann90%seedname, &
                    istdout, wann90%timer, error, wann90%comm)
      if (allocated(HH_R)) deallocate (HH_R)
      if (allocated(SS_R)) deallocate (SS_R)
      if (allocated(error)) then
        write (istderr, *) 'Error in dos', error%code, error%message
        ierr = sign(1, error%code)
        deallocate (error)
      endif
    else
      write (istdout, *) ' No dos calculation requested'
    endif
  end subroutine calc_dos

  subroutine boltzwann(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_boltzwann, only: boltzwann_main

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(pw90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)

    ierr = 0
    call boltzwann_main(pw90%boltzwann, wann90%dis_manifold, pw90%dos, wann90%kpt_latt, &
                        pw90%band_deriv_degen, pw90%oper_read, pw90%spin, physics, &
                        wann90%ws_region, wann90%w90_system, wann90%wannier_data, &
                        pw90%ws_distance, pw90%ws_vec, wann90%print_output, HH_R, SS_R, &
                        pw90%v_matrix, wann90%u_matrix, wann90%eigval, wann90%real_lattice, &
                        pw90%scissors_shift, wann90%mp_grid, wann90%num_wann, wann90%num_bands, &
                        wann90%num_kpts, pw90%effective_model, wann90%have_disentangled, &
                        pw90%calculation%spin_decomp, wann90%seedname, istdout, wann90%timer, &
                        error, wann90%comm)
    if (allocated(SS_R)) deallocate (SS_R)
    if (allocated(HH_R)) deallocate (HH_R)
    if (allocated(error)) then
      write (istderr, *) 'Error in boltzwann', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine boltzwann

  subroutine gyrotropic(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_gyrotropic, only: gyrotropic_main

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(pw90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: AA_R(:, :, :, :)
    complex(kind=dp), allocatable :: BB_R(:, :, :, :)
    complex(kind=dp), allocatable :: CC_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)

    ierr = 0
    call gyrotropic_main(pw90%berry, wann90%dis_manifold, wann90%fermi_energy_list, &
                         pw90%gyrotropic, wann90%kmesh_info, wann90%kpt_latt, physics, &
                         pw90%oper_read, pw90%band_deriv_degen, wann90%ws_region, &
                         wann90%w90_system, wann90%print_output, wann90%wannier_data, &
                         pw90%ws_vec, pw90%ws_distance, AA_R, BB_R, CC_R, HH_R, SS_R, &
                         wann90%u_matrix, pw90%v_matrix, wann90%eigval, wann90%real_lattice, &
                         pw90%scissors_shift, wann90%mp_grid, wann90%num_bands, wann90%num_kpts, &
                         wann90%num_wann, pw90%effective_model, wann90%have_disentangled, &
                         wann90%seedname, istdout, wann90%timer, error, wann90%comm)
    if (allocated(SS_R)) deallocate (SS_R)
    if (allocated(HH_R)) deallocate (HH_R)
    if (allocated(CC_R)) deallocate (CC_R)
    if (allocated(BB_R)) deallocate (BB_R)
    if (allocated(AA_R)) deallocate (AA_R)
    if (allocated(error)) then
      write (istderr, *) 'Error in gyrotropic', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine gyrotropic

  subroutine berry(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_berry, only: berry_main

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(pw90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: AA_R(:, :, :, :)
    complex(kind=dp), allocatable :: BB_R(:, :, :, :)
    complex(kind=dp), allocatable :: CC_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)
    complex(kind=dp), allocatable :: SR_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SHR_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SH_R(:, :, :, :)
    complex(kind=dp), allocatable :: SAA_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SBB_R(:, :, :, :, :)
    integer :: fermi_n

    ierr = 0
    fermi_n = 0
    if (allocated(wann90%fermi_energy_list)) fermi_n = size(wann90%fermi_energy_list)
    call berry_main(pw90%berry, wann90%dis_manifold, wann90%fermi_energy_list, wann90%kmesh_info, &
                    pw90%kpt_dist, wann90%kpt_latt, pw90%band_deriv_degen, pw90%oper_read, &
                    pw90%spin, physics, wann90%ws_region, pw90%spin_hall, wann90%wannier_data, &
                    pw90%ws_distance, pw90%ws_vec, wann90%print_output, AA_R, BB_R, CC_R, HH_R, &
                    SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, wann90%u_matrix, pw90%v_matrix, &
                    wann90%eigval, wann90%real_lattice, pw90%scissors_shift, wann90%mp_grid, &
                    fermi_n, wann90%num_wann, wann90%num_kpts, wann90%num_bands, &
                    wann90%w90_system%num_valence_bands, pw90%effective_model, &
                    wann90%have_disentangled, pw90%calculation%spin_decomp, &
                    wann90%seedname, istdout, wann90%timer, error, wann90%comm)
    if (allocated(SBB_R)) deallocate (SBB_R)
    if (allocated(SAA_R)) deallocate (SAA_R)
    if (allocated(SHR_R)) deallocate (SHR_R)
    if (allocated(SH_R)) deallocate (SH_R)
    if (allocated(SR_R)) deallocate (SR_R)
    if (allocated(SS_R)) deallocate (SS_R)
    if (allocated(HH_R)) deallocate (HH_R)
    if (allocated(CC_R)) deallocate (CC_R)
    if (allocated(BB_R)) deallocate (BB_R)
    if (allocated(AA_R)) deallocate (AA_R)
    if (allocated(error)) then
      write (istderr, *) 'Error in berry', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine berry

  subroutine kpath(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_kpath, only: k_path

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(pw90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: AA_R(:, :, :, :)
    complex(kind=dp), allocatable :: BB_R(:, :, :, :)
    complex(kind=dp), allocatable :: CC_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)
    complex(kind=dp), allocatable :: SR_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SHR_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SH_R(:, :, :, :)
    complex(kind=dp), allocatable :: SAA_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SBB_R(:, :, :, :, :)
    integer :: fermi_n

    ierr = 0
    fermi_n = 0
    if (allocated(wann90%fermi_energy_list)) fermi_n = size(wann90%fermi_energy_list)
    call k_path(pw90%berry, wann90%dis_manifold, wann90%fermi_energy_list, wann90%kmesh_info, &
                pw90%kpath, wann90%kpt_latt, pw90%oper_read, pw90%band_deriv_degen, pw90%spin, &
                wann90%ws_region, wann90%kpoint_path, pw90%spin_hall, wann90%print_output, &
                wann90%wannier_data, pw90%ws_distance, pw90%ws_vec, AA_R, BB_R, CC_R, HH_R, SH_R, &
                SHR_R, SR_R, SS_R, SAA_R, SBB_R, pw90%v_matrix, wann90%u_matrix, physics%bohr, &
                wann90%eigval, wann90%real_lattice, pw90%scissors_shift, wann90%mp_grid, fermi_n, &
                wann90%num_wann, wann90%num_bands, wann90%num_kpts, &
                wann90%w90_system%num_valence_bands, pw90%effective_model, &
                wann90%have_disentangled, wann90%seedname, istdout, wann90%timer, error, wann90%comm)
    if (allocated(SBB_R)) deallocate (SBB_R)
    if (allocated(SAA_R)) deallocate (SAA_R)
    if (allocated(SHR_R)) deallocate (SHR_R)
    if (allocated(SH_R)) deallocate (SH_R)
    if (allocated(SR_R)) deallocate (SR_R)
    if (allocated(SS_R)) deallocate (SS_R)
    if (allocated(HH_R)) deallocate (HH_R)
    if (allocated(CC_R)) deallocate (CC_R)
    if (allocated(BB_R)) deallocate (BB_R)
    if (allocated(AA_R)) deallocate (AA_R)
    if (allocated(error)) then
      write (istderr, *) 'Error in kpath', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine kpath

  subroutine kslice(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_kslice, only: k_slice

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(pw90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: AA_R(:, :, :, :)
    complex(kind=dp), allocatable :: BB_R(:, :, :, :)
    complex(kind=dp), allocatable :: CC_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)
    complex(kind=dp), allocatable :: SR_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SHR_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SH_R(:, :, :, :)
    complex(kind=dp), allocatable :: SAA_R(:, :, :, :, :)
    complex(kind=dp), allocatable :: SBB_R(:, :, :, :, :)
    integer :: fermi_n

    ierr = 0
    fermi_n = 0
    if (allocated(wann90%fermi_energy_list)) fermi_n = size(wann90%fermi_energy_list)
    call k_slice(pw90%berry, wann90%dis_manifold, wann90%fermi_energy_list, wann90%kmesh_info, &
                 wann90%kpt_latt, pw90%kslice, pw90%oper_read, pw90%band_deriv_degen, pw90%spin, &
                 wann90%ws_region, pw90%spin_hall, wann90%print_output, &
                 wann90%wannier_data, pw90%ws_distance, pw90%ws_vec, AA_R, BB_R, CC_R, HH_R, SH_R, &
                 SHR_R, SR_R, SS_R, SAA_R, SBB_R, pw90%v_matrix, wann90%u_matrix, physics%bohr, &
                 wann90%eigval, wann90%real_lattice, pw90%scissors_shift, wann90%mp_grid, fermi_n, &
                 wann90%num_bands, wann90%num_kpts, wann90%num_wann, &
                 wann90%w90_system%num_valence_bands, pw90%effective_model, &
                 wann90%have_disentangled, wann90%seedname, istdout, wann90%timer, error, wann90%comm)
    if (allocated(SBB_R)) deallocate (SBB_R)
    if (allocated(SAA_R)) deallocate (SAA_R)
    if (allocated(SHR_R)) deallocate (SHR_R)
    if (allocated(SH_R)) deallocate (SH_R)
    if (allocated(SR_R)) deallocate (SR_R)
    if (allocated(SS_R)) deallocate (SS_R)
    if (allocated(HH_R)) deallocate (HH_R)
    if (allocated(CC_R)) deallocate (CC_R)
    if (allocated(BB_R)) deallocate (BB_R)
    if (allocated(AA_R)) deallocate (AA_R)
    if (allocated(error)) then
      write (istderr, *) 'Error in kslice', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine kslice

  subroutine spin_moment(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_spin, only: spin_get_moment

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: HH_R(:, :, :)
    complex(kind=dp), allocatable :: SS_R(:, :, :, :)

    ierr = 0
    call spin_get_moment(wann90%dis_manifold, wann90%fermi_energy_list, pw90%kpt_dist, &
                         wann90%kpt_latt, pw90%oper_read, pw90%spin, wann90%ws_region, &
                         wann90%print_output, wann90%wannier_data, pw90%ws_distance, pw90%ws_vec, &
                         HH_R, SS_R, wann90%u_matrix, pw90%v_matrix, wann90%eigval, &
                         wann90%real_lattice, pw90%scissors_shift, wann90%mp_grid, &
                         wann90%num_wann, wann90%num_bands, wann90%num_kpts, &
                         wann90%w90_system%num_valence_bands, pw90%effective_model, &
                         wann90%have_disentangled, pw90%berry%wanint_kpoint_file, wann90%seedname, &
                         istdout, wann90%timer, error, wann90%comm)
    if (allocated(SS_R)) deallocate (SS_R)
    if (allocated(HH_R)) deallocate (HH_R)
    if (allocated(error)) then
      write (istderr, *) 'Error in spin_moment', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine spin_moment

  subroutine geninterp(wann90, pw90, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_geninterp, only: geninterp_main

    implicit none
    type(lib_common_type), intent(inout) :: wann90
    type(lib_postw90_type), intent(inout) :: pw90
    integer, intent(in) :: istdout, istderr
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: v_matrix(:, :, :)
    integer, intent(out) :: ierr
    !
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: HH_R(:, :, :)

    ierr = 0
    call geninterp_main(wann90%dis_manifold, pw90%geninterp, wann90%kpt_latt, &
                        pw90%band_deriv_degen, wann90%ws_region, wann90%print_output, &
                        wann90%wannier_data, pw90%ws_distance, pw90%ws_vec, HH_R, pw90%v_matrix, &
                        wann90%u_matrix, wann90%eigval, wann90%real_lattice, pw90%scissors_shift, &
                        wann90%mp_grid, wann90%num_bands, wann90%num_kpts, wann90%num_wann, &
                        wann90%w90_system%num_valence_bands, pw90%effective_model, &
                        wann90%have_disentangled, wann90%seedname, istdout, wann90%timer, &
                        error, wann90%comm)
    if (allocated(HH_R)) deallocate (HH_R)
    if (allocated(error)) then
      write (istderr, *) 'Error in geninterp', error%code, error%message
      ierr = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine geninterp

end module w90_lib_all
