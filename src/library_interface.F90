
module w90_helper_types

  ! as with fortran routines like allocate, the status variable indicates an error if non-zero
  ! positive is an error, negative is a warning (such as non-convergence) which is recoverable.
  use w90_constants
  use w90_types
  use w90_wannier90_types
  use w90_readwrite, only: update_settings

  implicit none

  ! should we have a lib_wannierise_type?

  type lib_global_type
    ! matrices
    complex(kind=dp), pointer :: u_opt(:, :, :) => null()
    complex(kind=dp), pointer :: u_matrix(:, :, :) => null()

    integer, pointer :: dist_kpoints(:) => null()
    !! distk(i) = node operating on k-point i

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
    type(w90_physical_constants_type) :: physics

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
    type(w90_calculation_type) :: w90_calculation ! separate this? ... maybe yes (JJ)
    ! matrices
    complex(kind=dp), pointer :: a_matrix(:, :, :) => null()
    complex(kind=dp), pointer :: m_matrix(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_matrix_local(:, :, :, :) => null()
    complex(kind=dp), pointer :: m_orig(:, :, :, :) => null()  !m_matrix_orig_local

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
    complex(kind=dp), allocatable :: ham_k(:, :, :)
    complex(kind=dp), allocatable :: ham_r(:, :, :)
    integer, allocatable :: irvec(:, :)
    integer, allocatable :: ndegen(:)
    integer, allocatable :: shift_vec(:, :)
    integer :: nrpts
    integer :: rpt_origin
    real(kind=dp), allocatable :: wannier_centres_translated(:, :)
    type(ham_logical_type) :: ham_logical
    type(sitesym_type) :: sitesym
    type(wann_omega_type) :: omega

    !plot
    type(band_plot_type) :: band_plot
    type(wannier_plot_type) :: wann_plot
    type(fermi_surface_plot_type) :: fermi_surface_data

    !transport
    type(transport_type) :: tran
  end type lib_w90_type

  public:: checkpoint, create_kmesh, get_fortran_stdout, get_fortran_stderr, input_reader, &
           overlaps, plot_files, print_times, transport, wannierise, write_kmesh

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

  subroutine get_fortran_stderr(output)
    use iso_fortran_env, only: error_unit
    implicit none
    integer, intent(out) :: output

    output = error_unit
  end subroutine get_fortran_stderr

  subroutine get_fortran_file(output, name)
    implicit none
    integer, intent(out) :: output
    character(len=*), intent(in) :: name

    open (newunit=output, file=name, form='formatted', status='unknown')
  end subroutine get_fortran_file

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

  subroutine input_reader(helper, wan90, seedname, output, outerr, status, comm)
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_uppercase, &
      w90_readwrite_clean_infile, w90_readwrite_read_final_alloc
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_extra_io_type
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_input
    use w90_comms, only: w90_comm_type, mpirank, comms_sync_err

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params
    logical :: cp_pp, disentanglement
    integer :: eunit

    status = 0

    call w90_readwrite_in_file(seedname, error, comm)
    if (allocated(error)) then
      write (outerr, *) 'Error in input file access', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
      return
    endif

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
                                      helper%real_lattice, helper%physics%bohr, &
                                      wan90%sitesym%symmetrize_eps, helper%mp_grid, &
                                      helper%num_bands, helper%num_kpts, wan90%num_proj, &
                                      helper%num_wann, wan90%optimisation, wan90%eig_found, &
                                      wan90%calc_only_A, cp_pp, helper%gamma_only, &
                                      wan90%lhasproj, .false., .false., wan90%lsitesymmetry, &
                                      wan90%use_bloch_phases, seedname, output, error, comm)

    ! test mpi error handling using "unlucky" input token
    ! this machinery used to sit in w90_wannier90_readwrite_dist
    ! but that routine is obsolete if input file is read on all ranks
    ! fixme, this should be moved to wannier90 main routine (definately doesn't belong here)
    if (helper%print_output%timing_level < 0 .and. mpirank(comm) == abs(helper%print_output%timing_level)) then
      write (*, *) "here...", mpirank(comm)
      call set_error_input(error, 'received unlucky_rank', comm)
    else
      call comms_sync_err(comm, error, 0) ! this is necessary since non-root may never enter an mpi collective if root has exited here
    endif
    if (allocated(error)) then ! applies (is t) for all ranks now
      open (newunit=eunit, file=seedname//'.werr')
      call prterr(error, output, eunit, comm)
      close (eunit)
      status = sign(1, error%code)
      deallocate (error)
    endif
    !!!!! end unlucky code

    if (allocated(error)) then
      write (outerr, *) 'Error in reader', error%code, error%message
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
        write (outerr, *) 'Error in read alloc', error%code, error%message
        status = sign(1, error%code)
        deallocate (error)
      endif
    endif
    call w90_readwrite_clean_infile(output, seedname, error, comm)
    if (allocated(error)) then
      write (outerr, *) 'Error in input close', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
    helper%seedname = seedname ! maybe not keep this separate from "blob"? JJ

    if (mpirank(comm) /= 0) helper%print_output%iprint = 0 ! supress printing non-rank-0
  end subroutine input_reader

  subroutine checkpoint(helper, wan90, label, output, comm)
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write_chkpt
    use w90_comms, only: w90_comm_type, mpirank

    ! write_chkpt never fails?  remarkable.
    ! either before or at start of write_chkpt, a reduction on m_matrix is necessary
    ! fixme JJ check

    ! fixme, seedname might rather prefer to be an argument

    implicit none

    character(len=*), intent(in) :: label
    integer, intent(in) :: output
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(in) :: wan90
    type(w90_comm_type), intent(in) :: comm

    if (mpirank(comm) == 0) then
      ! e.g. label = 'postwann' after wannierisation
      call w90_wannier90_readwrite_write_chkpt(label, helper%exclude_bands, helper%wannier_data, &
                                               helper%kmesh_info, helper%kpt_latt, &
                                               helper%num_kpts, helper%dis_manifold, &
                                               helper%num_bands, helper%num_wann, helper%u_matrix, &
                                               helper%u_opt, wan90%m_matrix, helper%mp_grid, &
                                               helper%real_lattice, wan90%omega%invariant, &
                                               helper%have_disentangled, output, helper%seedname)
    endif
  end subroutine checkpoint

  subroutine create_kmesh(helper, output, outerr, status, comm)
    use w90_kmesh, only: kmesh_get
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(lib_global_type), intent(inout) :: helper
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    status = 0
    call kmesh_get(helper%kmesh_input, helper%kmesh_info, helper%print_output, helper%kpt_latt, &
                   helper%real_lattice, helper%num_kpts, helper%gamma_only, output, helper%timer, &
                   error, comm)
    if (allocated(error)) then
      write (outerr, *) 'Error in kmesh_get', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine create_kmesh

  subroutine write_kmesh(helper, wan90, seedname, output, outerr, status, comm)
    use w90_kmesh, only: kmesh_get, kmesh_write
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname ! needed for nnkp filenames
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    type(w90_error_type), allocatable :: error
    logical :: calc_only_A = .false.! what does this do?

    status = 0
    call kmesh_write(helper%exclude_bands, helper%kmesh_info, wan90%proj, helper%print_output, &
                     helper%kpt_latt, helper%real_lattice, helper%num_kpts, wan90%num_proj, &
                     calc_only_A, helper%w90_system%spinors, seedname, helper%timer)

    if (allocated(error)) then
      write (outerr, *) 'Error in kmesh_write', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine write_kmesh

  subroutine overlaps(helper, wan90, output, outerr, status, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type
    use w90_overlap, only: overlap_read

    implicit none

    ! arguments
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    logical :: cp_pp = .false.
    type(w90_error_type), allocatable :: error

    status = 0

    if (helper%num_bands > helper%num_wann) then ! disentanglement case
      if ((.not. associated(wan90%a_matrix)) .or. (.not. associated(wan90%m_orig))) then
        write (outerr, *) 'Matrices not set for overlap call (disentanglement case)'
        status = 1
        return
      endif
      call overlap_read(helper%kmesh_info, wan90%select_proj, wan90%sitesym, wan90%a_matrix, &
                        wan90%m_orig, helper%num_bands, helper%num_kpts, wan90%num_proj, &
                        helper%num_wann, helper%print_output%timing_level, cp_pp, &
                        helper%gamma_only, wan90%lsitesymmetry, wan90%use_bloch_phases, &
                        helper%seedname, output, helper%timer, helper%dist_kpoints, error, comm)
    else
      if ((.not. associated(helper%u_matrix)) .or. (.not. associated(wan90%m_matrix_local))) then
        write (outerr, *) 'Matrices not set for overlap call'
        status = 1
        return
      endif
      call overlap_read(helper%kmesh_info, wan90%select_proj, wan90%sitesym, helper%u_matrix, &
                        wan90%m_matrix_local, helper%num_bands, helper%num_kpts, wan90%num_proj, &
                        helper%num_wann, helper%print_output%timing_level, cp_pp, &
                        helper%gamma_only, wan90%lsitesymmetry, wan90%use_bloch_phases, &
                        helper%seedname, output, helper%timer, helper%dist_kpoints, error, comm)
    endif
    if (allocated(error)) then
      write (outerr, *) 'Error in overlaps', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine overlaps

  subroutine disentangle(helper, wan90, output, outerr, status, comm)

    use w90_disentangle, only: dis_main, setup_m_loc
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank

    implicit none

    ! arguments
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error
    integer :: optimisation = 3

    status = 0

    if (.not. associated(wan90%m_orig)) then  ! m_matrix_orig_local (nband*nwann for disentangle)
      write (outerr, *) 'm_orig not set for disentangle call'
      status = 1
      return
    else if (.not. associated(wan90%m_matrix_local)) then ! (nband*nwann*nknode for wannierise)
      write (outerr, *) 'm_matrix_local not set for disentangle call'
      status = 1
      return
    else if (.not. associated(wan90%a_matrix)) then
      write (outerr, *) 'a_matrix not set for disentangle call'
      status = 1
      return
    else if (.not. associated(helper%u_matrix)) then
      write (outerr, *) 'u_matrix not set for disentangle call'
      status = 1
      return
    else if (.not. associated(helper%u_opt)) then
      write (outerr, *) 'u_opt not set for disentangle call'
      status = 1
      return
    else if (.not. associated(helper%dist_kpoints)) then
      write (outerr, *) 'kpt decomp not set for disentangle call'
      status = 1
      return
    endif

    call dis_main(wan90%dis_control, wan90%dis_spheres, helper%dis_manifold, helper%kmesh_info, &
                  helper%kpt_latt, wan90%sitesym, helper%print_output, wan90%a_matrix, &
                  wan90%m_orig, helper%u_matrix, helper%u_opt, helper%eigval, &
                  helper%real_lattice, wan90%omega%invariant, helper%num_bands, helper%num_kpts, &
                  helper%num_wann, helper%gamma_only, wan90%lsitesymmetry, output, helper%timer, &
                  helper%dist_kpoints, error, comm)
    if (allocated(error)) then
      write (outerr, *) 'Error in disentangle (dis_main)', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
      return
    endif

    ! copy to m_matrix_local from m_matrix_orig_local (aka m_orig)
    call setup_m_loc(helper%kmesh_info, helper%print_output, wan90%m_matrix_local, wan90%m_orig, &
                     helper%u_matrix, helper%num_bands, helper%num_kpts, helper%num_wann, &
                     optimisation, helper%timer, helper%dist_kpoints, error, comm)
    if (allocated(error)) then
      write (outerr, *) 'Error in disentangle (setup_m_loc)', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
      return
    endif

    helper%have_disentangled = .true.
  end subroutine disentangle

  subroutine wannierise(helper, wan90, output, outerr, status, comm)

    use w90_wannierise, only: wann_main, wann_main_gamma
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(lib_global_type), intent(inout) :: helper
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    if ((.not. associated(wan90%m_matrix_local)) .or. (.not. associated(helper%u_opt)) .or. &
        (.not. associated(helper%u_matrix)) .or. (.not. associated(helper%dist_kpoints))) then
      write (outerr, *) 'Matrices not set for wannierise call'
      status = 1
      return
    endif
    status = 0
    if (helper%gamma_only) then
      call wann_main_gamma(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                           helper%kmesh_info, helper%kpt_latt, wan90%output_file, wan90%wann_control, &
                           wan90%omega, helper%w90_system, helper%print_output, helper%wannier_data, wan90%m_matrix_local, &
                           helper%u_matrix, helper%u_opt, helper%eigval, helper%real_lattice, helper%mp_grid, &
                           helper%num_bands, helper%num_kpts, helper%num_wann, helper%have_disentangled, &
                           wan90%real_space_ham%translate_home_cell, helper%seedname, output, &
                           helper%timer, error, comm)
    else
      call wann_main(helper%atom_data, helper%dis_manifold, helper%exclude_bands, &
                     wan90%ham_logical, helper%kmesh_info, helper%kpt_latt, wan90%output_file, &
                     wan90%real_space_ham, wan90%wann_control, wan90%omega, wan90%sitesym, &
                     helper%w90_system, helper%print_output, helper%wannier_data, helper%ws_region, &
                     wan90%w90_calculation, wan90%ham_k, wan90%ham_r, &
                     wan90%m_matrix_local, helper%u_matrix, helper%u_opt, &
                     helper%eigval, helper%real_lattice, wan90%wannier_centres_translated, wan90%irvec, &
                     helper%mp_grid, wan90%ndegen, wan90%shift_vec, wan90%nrpts, helper%num_bands, &
                     helper%num_kpts, wan90%num_proj, helper%num_wann, wan90%optimisation, &
                     wan90%rpt_origin, wan90%band_plot%mode, wan90%tran%mode, &
                     helper%have_disentangled, wan90%lsitesymmetry, helper%seedname, output, &
                     helper%timer, helper%dist_kpoints, error, comm)
    endif
    if (allocated(error)) then
      write (outerr, *) 'Error in wannierise', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine wannierise

  subroutine plot_files(helper, wan90, output, outerr, status, comm)
    use w90_plot, only: plot_main
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(lib_global_type), intent(inout) :: helper ! inout due to ham_logical -- JJ: eugh, that's nasty, can we change it?
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
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
      write (outerr, *) 'Error in plotting', error%code, error%message
      status = sign(1, error%code)
      deallocate (error)
    endif
  end subroutine plot_files

  subroutine transport(helper, wan90, output, outerr, status, comm)
    use w90_error_base, only: w90_error_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_transport, only: tran_main

    implicit none

    ! arguments
    integer, intent(in) :: output, outerr
    integer, intent(out) :: status
    type(lib_global_type), intent(inout) :: helper ! because of ham_logical
    type(lib_w90_type), intent(inout) :: wan90
    type(w90_comm_type), intent(in) :: comm

    ! local
    type(w90_error_type), allocatable :: error

    ! should these tests be done outside? -- JJ: is tran_main pllel at all?
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
          write (outerr, *) 'Error in transport', error%code, error%message
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

    if (helper%print_output%iprint > 0) call io_print_timings(helper%timer, output)
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

  subroutine prterr(error, stdout, stderr, comm)
    use w90_error_base, only: code_remote, w90_error_type
    use w90_comms, only: comms_no_sync_send, comms_no_sync_recv, w90_comm_type, mpirank, mpisize

    type(w90_error_type), allocatable, intent(in) :: error
    integer, intent(in) :: stderr, stdout
    type(w90_comm_type), intent(in) :: comm

    type(w90_error_type), allocatable :: le ! unchecked error state for calls made in this routine
    integer :: ie ! global error value to be returned
    integer :: je ! error value on remote ranks
    integer :: j ! rank index
    integer :: failrank ! lowest rank reporting an error
    character(len=128) :: mesg ! only print 128 chars of error

    ie = 0
    mesg = 'not set'

    if (mpirank(comm) == 0) then
      ! fixme, report all failing ranks
      do j = mpisize(comm) - 1, 1, -1
        call comms_no_sync_recv(je, 1, j, le, comm)

        if (je /= code_remote .and. je /= 0) then
          failrank = j
          ie = je
          call comms_no_sync_recv(mesg, 128, j, le, comm)
        endif
      enddo
      ! if the error is on rank0
      if (error%code /= code_remote .and. error%code /= 0) then
        failrank = 0
        ie = error%code
        mesg = error%message
      endif

      !if (ie == 0) write (stderr, *) "logic error" ! to arrive here requires this

      write (stderr, *) 'Exiting.......'
      write (stderr, '(1x,a)') trim(mesg)
      write (stderr, '(1x,a,i0,a)') '(rank: ', failrank, ')'
      write (stdout, '(1x,a)') ' error encountered; check error .werr error log'

    else ! non 0 ranks
      je = error%code
      call comms_no_sync_send(je, 1, 0, le, comm)
      if (je /= code_remote .and. je /= 0) then
        mesg = error%message
        call comms_no_sync_send(mesg, 128, 0, le, comm)
      endif
    endif

#ifdef MPI
    call mpi_finalize(je) ! je overwritten here
#endif
    !call exit(ie) ! return true fail code (gnu extension)
    stop
  end subroutine prterr

end module w90_helper_types
