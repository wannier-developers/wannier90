
module w90_helper_types

  ! as with fortran routines like allocate, the status variable indicates an error if non-zero
  ! positive is an error, negative is a warning (such as non-convergence) which is recoverable.
  use w90_constants
  use w90_types
  use w90_wannier90_types
  use w90_readwrite, only: update_settings
  use iso_fortran_env, only: error_unit

  implicit none

  ! should we have a lib_wannierise_type?

  type lib_global_type
    ! matrices
    complex(kind=dp), pointer :: u_opt(:, :, :) => null()
    complex(kind=dp), pointer :: u_matrix(:, :, :) => null()

    ! array of num_kpts containing the nodeid that the kpoint is on
    integer, pointer :: dist_kpoints(:) => null()
    ! num kpoints per rank (size of contiguous blcok)
    integer, pointer :: counts(:) => null()
    ! kpoint block offset
    integer, pointer :: displs(:) => null()

    type(atom_data_type) :: atom_data
    !type(dis_control_type) :: dis_control
    type(dis_manifold_type) :: dis_manifold
    !type(dis_spheres_type) :: dis_spheres
    type(kmesh_info_type) :: kmesh_info
    type(kmesh_input_type) :: kmesh_input
    type(kpoint_path_type) :: kpoint_path
    !type(output_file_type) :: output_file
    type(print_output_type) :: print_output
    !type(proj_input_type) :: proj
    !type(proj_input_type) :: proj_input
    !type(real_space_ham_type) :: real_space_ham
    !type(select_projection_type) :: select_proj
    type(w90_system_type) :: w90_system
    !type(wann_control_type) :: wann_control
    type(wannier_data_type) :: wannier_data
    !type(wann_omega_type) :: wann_omega
    type(ws_region_type) :: ws_region
    !type(wvfn_read_type) :: wvfn_read

    type(timer_list_type) :: timer

    integer, allocatable :: exclude_bands(:)
    integer :: mp_grid(3)
    integer :: num_bands
    integer :: num_kpts
    !integer :: num_proj = 0
    integer :: num_wann = -99
    !integer :: optimisation = 3

    real(kind=dp), allocatable :: eigval(:, :)
    real(kind=dp), allocatable :: fermi_energy_list(:)
    real(kind=dp), allocatable :: kpt_latt(:, :)
    real(kind=dp) :: real_lattice(3, 3)
    !!real(kind=dp) :: symmetrize_eps

    !logical :: eig_found = .false.
    !Projections
    !logical :: lhasproj = .false.
    ! RS: symmetry-adapted Wannier functions
    !logical :: lsitesymmetry = .false.
    !logical :: use_bloch_phases = .false.
    !logical :: calc_only_A = .false.
    logical :: gamma_only

    ! added for wannierise
    !type(wann_omega_type) :: omega
    !type(ham_logical_type) :: ham_logical
    !type(sitesym_type) :: sitesym
    !complex(kind=dp), allocatable :: ham_k(:, :, :)
    !complex(kind=dp), allocatable :: ham_r(:, :, :)
    !real(kind=dp), allocatable :: wannier_centres_translated(:, :)
    !integer, allocatable :: irvec(:, :)
    !integer, allocatable :: shift_vec(:, :)
    !integer, allocatable :: ndegen(:)
    !integer :: rpt_origin
    !integer :: nrpts
    logical :: have_disentangled = .false.
    character(len=128) :: seedname
  end type lib_global_type

  type lib_w90_type
    type(w90_calculation_type) :: w90_calculation ! separate this?
    ! matrices
    complex(kind=dp), pointer :: a_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: m_matrix(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_matrix_local(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_orig(:, :, :, :) => null()

    type(dis_control_type) :: dis_control
    type(dis_spheres_type) :: dis_spheres
    type(output_file_type) :: output_file
    type(proj_input_type) :: proj
    type(proj_input_type) :: proj_input
    type(real_space_ham_type) :: real_space_ham
    type(select_projection_type) :: select_proj
    type(wann_control_type) :: wann_control
    type(wann_omega_type) :: wann_omega
    type(wvfn_read_type) :: wvfn_read

    !type(timer_list_type) :: timer

    integer :: num_proj = 0
    integer :: optimisation = 3

    logical :: eig_found = .false.
    !Projections
    logical :: lhasproj = .false.
    ! RS: symmetry-adapted Wannier functions
    logical :: lsitesymmetry = .false.
    logical :: use_bloch_phases = .false.
    logical :: calc_only_A = .false.

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

    !plot
    type(band_plot_type) :: band_plot
    type(wannier_plot_type) :: wann_plot
    type(fermi_surface_plot_type) :: fermi_surface_data

    !transport
    type(transport_type) :: tran
  end type lib_w90_type

  public:: input_reader, create_kmesh, checkpoint, overlaps, wannierise, plot_files, &
           transport, print_times, get_fortran_stdout

  public :: set_option
  interface set_option
    module procedure set_option_bool
    !module procedure set_option_cplx
    module procedure set_option_text
    module procedure set_option_dble
    module procedure set_option_int
  end interface set_option

contains

  subroutine get_fortran_stdout(output)
    use iso_fortran_env, only: output_unit
    implicit none
    integer, intent(out) :: output

    output = output_unit
  end subroutine get_fortran_stdout

  subroutine set_option_bool(string, bool)
    implicit none
    character(*), intent(in) :: string
    logical, intent(in) :: bool
    call update_settings(string, bool, "", 0.d0, 0)
  endsubroutine set_option_bool

  !subroutine set_option_cplx(string,cval)
  !  implicit none
  !  character(*), intent(in) :: string
  !  complex(kind=dp), intent(in) :: cval
  !  call update_settings(string, .false., cval, 0.d0, 0)
  !endsubroutine set_option_cplx

  subroutine set_option_dble(string, rval)
    implicit none
    character(*), intent(in) :: string
    real(kind=dp), intent(in) :: rval
    call update_settings(string, .false., "", rval, 0)
  endsubroutine set_option_dble

  subroutine set_option_int(string, ival)
    implicit none
    character(*), intent(in) :: string
    integer, intent(in) :: ival
    call update_settings(string, .false., "", 0.d0, ival)
  endsubroutine set_option_int

  subroutine set_option_text(string, text)
    implicit none
    character(*), intent(in) :: string
    character(*), intent(in) :: text
    call update_settings(string, .false., text, 0.d0, 0)
  endsubroutine set_option_text

  subroutine input_reader(helper, wan90, seedname, output, status, comm)
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_uppercase, &
      w90_readwrite_clean_infile, w90_readwrite_read_final_alloc
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    implicit none
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: output
    character(len=*), intent(in) :: seedname
    integer, intent(out) :: status
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    logical :: cp_pp, disentanglement

    status = 0
    call w90_readwrite_in_file(seedname, error, comm)
    if (allocated(error)) then
      write (error_unit, *) 'Error in input file access', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    else
      call w90_wannier90_readwrite_read(helper%atom_data, wan90%band_plot, wan90%dis_control, &
                                        wan90%dis_spheres, helper%dis_manifold, &
                                        helper%exclude_bands, helper%fermi_energy_list, &
                                        wan90%fermi_surface_data, helper%kmesh_input, &
                                        helper%kmesh_info, helper%kpt_latt, wan90%output_file, &
                                        wan90%wvfn_read, wan90%wann_control, wan90%proj, &
                                        wan90%proj_input, wan90%real_space_ham, wan90%select_proj, &
                                        helper%kpoint_path, helper%w90_system, wan90%tran, &
                                        helper%print_output, wan90%wann_plot, io_params, &
                                        helper%ws_region, wan90%w90_calculation, helper%eigval, &
                                        helper%real_lattice, physics%bohr, &
                                        wan90%sitesym%symmetrize_eps, helper%mp_grid, &
                                        helper%num_bands, helper%num_kpts, wan90%num_proj, &
                                        helper%num_wann, wan90%optimisation, wan90%eig_found, &
                                        wan90%calc_only_A, cp_pp, helper%gamma_only, &
                                        wan90%lhasproj, .false., .false., wan90%lsitesymmetry, &
                                        wan90%use_bloch_phases, seedname, output, error, comm)
      if (allocated(error)) then
        write (error_unit, *) 'Error in reader', error%code, error%message
        status = sign(1, error%code)
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
          status = sign(1, error%code)
          deallocate (error)
        endif
      endif
      call w90_readwrite_clean_infile(output, seedname, error, comm)
      if (allocated(error)) then
        write (error_unit, *) 'Error in input close', error%code, error%message
        status = sign(1, error%code)
        deallocate (error)
      endif
      helper%seedname = seedname
      if (mpirank(comm) /= 0) helper%print_output%iprint = 0
    endif
  end subroutine input_reader

  subroutine checkpoint(helper, wan90, label, output, comm)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write_chkpt
    !use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    implicit none
    type(lib_global_type), intent(inout) :: helper ! matrices
    type(lib_w90_type), intent(in) :: wan90
    character(len=*), intent(in) :: label
    !complex(kind=dp), intent(in) :: u_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    integer, intent(in) :: output
    !integer, intent(out) :: status
    type(w90comm_type), intent(in) :: comm

    if (mpirank(comm) == 0) then
      ! e.g. label = 'postwann' after wannierisation
      call w90_wannier90_readwrite_write_chkpt(label, helper%exclude_bands, helper%wannier_data, &
                                               helper%kmesh_info, helper%kpt_latt, helper%num_kpts, helper%dis_manifold, &
                                               helper%num_bands, helper%num_wann, helper%u_matrix, &
                                               helper%u_opt, wan90%m_matrix, helper%mp_grid, &
                                               helper%real_lattice, wan90%omega%invariant, helper%have_disentangled, output, &
                                               helper%seedname)
    endif
  end subroutine checkpoint

  subroutine create_kmesh(helper, output, status, comm)
    use w90_kmesh, only: kmesh_get
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    integer, intent(in) :: output
    integer, intent(out) :: status
    type(w90comm_type), intent(in) :: comm
    !
    type(w90_error_type), allocatable :: error

    status = 0
    call kmesh_get(helper%kmesh_input, helper%kmesh_info, helper%print_output, helper%kpt_latt, &
                   helper%real_lattice, helper%num_kpts, helper%gamma_only, output, helper%timer, &
                   error, comm)
    if (allocated(error)) then
      write (error_unit, *) 'Error in kmesh_get', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine create_kmesh

  subroutine overlaps(helper, wan90, output, status, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    use w90_overlap, only: overlap_read
    implicit none
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: output
    integer, intent(out) :: status
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

    if ((.not. associated(wan90%m_matrix)) .or. (.not. associated(wan90%m_orig)) .or. &
        (.not. associated(wan90%a_matrix)) .or. (.not. associated(helper%u_matrix))) then
      write (error_unit, *) 'Matrices not set for overlap call'
      status = 1
      return
    endif
    status = 0
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
    allocate (m_matrix_local(helper%num_wann, helper%num_wann, helper%kmesh_info%nntot, &
                             helper%num_kpts))
    call overlap_read(helper%kmesh_info, wan90%select_proj, wan90%sitesym, wan90%a_matrix, &
                      wan90%m_matrix, m_matrix_local, wan90%m_orig, m_matrix_orig_local, &
                      helper%u_matrix, u_matrix_opt, helper%num_bands, helper%num_kpts, &
                      wan90%num_proj, helper%num_wann, helper%print_output%timing_level, cp_pp, &
                      helper%gamma_only, wan90%lsitesymmetry, wan90%use_bloch_phases, &
                      helper%seedname, output, helper%timer, error, comm)
    deallocate (m_matrix_local)
    if (helper%num_bands > helper%num_wann) then
      deallocate (m_matrix_orig_local)
      !deallocate (m_matrix_orig)
    endif
    if (allocated(u_matrix_opt)) deallocate (u_matrix_opt)
    if (allocated(error)) then
      write (error_unit, *) 'Error in overlaps', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine overlaps

  subroutine disentangle(helper, wan90, output, status, comm)

    ! fixme, probably good to switch from m_matrix_local to m_matrix when mpisize==1

    use w90_disentangle, only: dis_main, splitm
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank

    implicit none

    ! arguments
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: output
    integer, intent(out) :: status
    type(w90comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error
    integer :: optimisation = 3

    status = 0

    if (.not. associated(wan90%m_orig)) then  ! m_matrix_orig_local (nband*nwann for disentangle)
      write (error_unit, *) 'm_orig not set for disentangle call'
      status = 1
      return
    else if (.not. associated(wan90%m_matrix) .and. mpirank(comm) == 0) then ! (nband*nwann*allk root only for wannierise)
      write (error_unit, *) 'm_matrix not set for disentangle call'
      status = 1
      return
    else if (.not. associated(wan90%m_matrix_local)) then ! (nband*nwann*nknode for wannierise)
      write (error_unit, *) 'm_matrix_local not set for disentangle call'
      status = 1
      return
    else if (.not. associated(wan90%a_matrix)) then
      write (error_unit, *) 'a_matrix not set for disentangle call'
      status = 1
      return
    else if (.not. associated(helper%u_matrix)) then
      write (error_unit, *) 'u_matrix not set for disentangle call'
      status = 1
      return
    else if (.not. associated(helper%u_opt)) then
      write (error_unit, *) 'u_opt not set for disentangle call'
      status = 1
      return
    else if (.not. associated(helper%counts)) then
      write (error_unit, *) 'kpt decomp not set for disentangle call'
      status = 1
      return
    else if (.not. associated(helper%displs)) then
      write (error_unit, *) 'kpt decomp not set for disentangle call'
      status = 1
      return
    endif

    call dis_main(wan90%dis_control, wan90%dis_spheres, helper%dis_manifold, helper%kmesh_info, &
                  helper%kpt_latt, wan90%sitesym, helper%print_output, wan90%a_matrix, &
                  wan90%m_orig, helper%u_matrix, helper%u_opt, helper%eigval, &
                  helper%real_lattice, wan90%omega%invariant, helper%num_bands, helper%num_kpts, &
                  helper%num_wann, helper%gamma_only, wan90%lsitesymmetry, output, helper%timer, &
                  helper%counts, helper%displs, error, comm)

    ! copy to m_matrix_local and m_matrix (on root) from m_matrix_orig_local
    call splitm(helper%kmesh_info, helper%print_output, wan90%m_matrix_local, wan90%m_orig, &
                wan90%m_matrix, helper%u_matrix, helper%num_bands, helper%num_kpts, &
                helper%num_wann, optimisation, helper%timer, helper%counts, &
                helper%displs, error, comm)

    helper%have_disentangled = .true.

    if (allocated(error)) then
      write (error_unit, *) 'Error in disentangle', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine disentangle

  subroutine wannierise(helper, wan90, output, status, comm)
    use w90_wannierise, only: wann_main, wann_main_gamma
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: output
    integer, intent(out) :: status
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(in) :: u_opt(:, :, :)
    !complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error

    if ((.not. associated(wan90%m_matrix)) .or. (.not. associated(helper%u_opt)) .or. &
        (.not. associated(helper%u_matrix)) .or. (.not. associated(helper%dist_kpoints))) then
      write (error_unit, *) 'Matrices not set for wannierise call'
      status = 1
      return
    endif
    status = 0
    if (helper%gamma_only) then
      call wann_main_gamma(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                           helper%kmesh_info, helper%kpt_latt, wan90%output_file, wan90%wann_control, &
                           wan90%omega, helper%w90_system, helper%print_output, helper%wannier_data, wan90%m_matrix, &
                           wan90%m_matrix_local, &
                           helper%u_matrix, helper%u_opt, helper%eigval, helper%real_lattice, helper%mp_grid, &
                           helper%num_bands, helper%num_kpts, helper%num_wann, helper%have_disentangled, &
                           wan90%real_space_ham%translate_home_cell, helper%seedname, output, &
                           helper%timer, error, comm)
    else
      call wann_main(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                     wan90%ham_logical, helper%kmesh_info, helper%kpt_latt, wan90%output_file, &
                     wan90%real_space_ham, wan90%wann_control, wan90%omega, wan90%sitesym, &
                     helper%w90_system, helper%print_output, helper%wannier_data, helper%ws_region, &
                     wan90%w90_calculation, wan90%ham_k, wan90%ham_r, wan90%m_matrix, &
                     wan90%m_matrix_local, helper%u_matrix, helper%u_opt, &
                     helper%eigval, helper%real_lattice, wan90%wannier_centres_translated, wan90%irvec, &
                     helper%mp_grid, wan90%ndegen, wan90%shift_vec, wan90%nrpts, helper%num_bands, &
                     helper%num_kpts, wan90%num_proj, helper%num_wann, wan90%optimisation, &
                     wan90%rpt_origin, wan90%band_plot%mode, wan90%tran%mode, &
                     helper%have_disentangled, wan90%lsitesymmetry, helper%seedname, output, &
                     helper%timer, helper%dist_kpoints, error, comm)
    endif
    if (allocated(error)) then
      write (error_unit, *) 'Error in wannierise', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine wannierise

  subroutine plot_files(helper, wan90, output, status, comm)
    use w90_plot, only: plot_main
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type
    implicit none
    type(lib_global_type), intent(inout) :: helper ! inout due to ham_logical
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: output
    integer, intent(out) :: status
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    !complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    !complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    !
    type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error

    status = 0
    call plot_main(helper%atom_data, wan90%band_plot, helper%dis_manifold, helper%fermi_energy_list, &
                   wan90%fermi_surface_data, wan90%ham_logical, helper%kmesh_info, helper%kpt_latt, &
                   wan90%output_file, wan90%wvfn_read, wan90%real_space_ham, helper%kpoint_path, &
                   helper%print_output, helper%wannier_data, wan90%wann_plot, helper%ws_region, &
                   wan90%w90_calculation, wan90%ham_k, wan90%ham_r, wan90%m_matrix, helper%u_matrix, &
                   helper%u_opt, helper%eigval, helper%real_lattice, wan90%wannier_centres_translated, &
                   physics%bohr, wan90%irvec, helper%mp_grid, wan90%ndegen, wan90%shift_vec, wan90%nrpts, &
                   helper%num_bands, helper%num_kpts, helper%num_wann, wan90%rpt_origin, &
                   wan90%tran%mode, helper%have_disentangled, wan90%lsitesymmetry, helper%w90_system%spinors, &
                   helper%seedname, output, helper%timer, error, comm)
    if (allocated(error)) then
      write (error_unit, *) 'Error in plotting', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine plot_files

  subroutine transport(helper, wan90, output, status, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90comm_type, mpirank
    use w90_transport, only: tran_main
    implicit none
    type(lib_global_type), intent(inout) :: helper ! because of ham_logical
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: output
    integer, intent(out) :: status
    type(w90comm_type), intent(in) :: comm
    !complex(kind=dp), intent(in) :: u_opt(:, :, :)
    !complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    !
    !type(w90_physical_constants_type) :: physics
    type(w90_error_type), allocatable :: error

    ! should these tests be done outside?
    status = 0
    if (mpirank(comm) == 0) then
      if (wan90%w90_calculation%transport) then
        call tran_main(helper%atom_data, helper%dis_manifold, helper%fermi_energy_list, &
                       wan90%ham_logical, helper%kpt_latt, wan90%output_file, wan90%real_space_ham, &
                       wan90%tran, helper%print_output, helper%wannier_data, helper%ws_region, &
                       wan90%w90_calculation, wan90%ham_k, wan90%ham_r, helper%u_matrix, helper%u_opt, &
                       helper%eigval, helper%real_lattice, wan90%wannier_centres_translated, wan90%irvec, &
                       helper%mp_grid, wan90%ndegen, wan90%shift_vec, wan90%nrpts, helper%num_bands, &
                       helper%num_kpts, helper%num_wann, wan90%rpt_origin, wan90%band_plot%mode, &
                       helper%have_disentangled, wan90%lsitesymmetry, helper%seedname, output, &
                       helper%timer, error, comm)
        if (allocated(error)) then
          write (error_unit, *) 'Error in transport', error%code, error%message
          status = sign(1, error%code)
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
    type(lib_w90_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: a_matrix(:, :, :)

    helper%a_matrix => a_matrix
  end subroutine set_a_matrix

  subroutine set_m_matrix(helper, m_matrix)
    implicit none
    type(lib_w90_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: m_matrix(:, :, :, :)

    helper%m_matrix => m_matrix
  end subroutine set_m_matrix

  subroutine set_m_matrix_local(helper, m_matrix_local) ! scattered m-matrix
    implicit none
    type(lib_w90_type), intent(inout) :: helper
    complex(kind=dp), intent(inout), target :: m_matrix_local(:, :, :, :)

    helper%m_matrix_local => m_matrix_local
  end subroutine set_m_matrix_local

  subroutine set_m_orig(helper, m_orig) ! m_matrix_local_orig
    implicit none
    type(lib_w90_type), intent(inout) :: helper
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

  subroutine set_kpoint_distribution(helper, dist)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    integer, intent(inout), target :: dist(:)

    helper%dist_kpoints => dist
  end subroutine set_kpoint_distribution

  ! for "old style" blockwise distribution
  subroutine set_kpoint_block(helper, counts, displs)
    implicit none
    type(lib_global_type), intent(inout) :: helper
    integer, intent(inout), target :: counts(:)
    integer, intent(inout), target :: displs(:)

    helper%counts => counts
    helper%displs => displs
  end subroutine set_kpoint_block

end module w90_helper_types
