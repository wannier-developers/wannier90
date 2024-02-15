!-*- mode: F90 -*-!
!------------------------------------------------------------!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Please cite                                                !
!                                                            !
!  [ref] "Wannier90 as a community code:                     !
!        new features and applications",                     !
!        G. Pizzi et al.,  J. Phys. Cond. Matt. 32,          !
!        165902 (2020).                                      !
!        http://doi.org/10.1088/1361-648X/ab51ff             !
!                                                            !
! in any publications arising from the use of this code.     !
!                                                            !
! Wannier90 is based on Wannier77, written by N. Marzari,    !
! I. Souza and D. Vanderbilt. For the method please cite     !
!                                                            !
! [ref] N. Marzari and D. Vanderbilt,                        !
!       Phys. Rev. B 56 12847 (1997)                         !
!       http://dx.doi.org/10.1103/PhysRevB.56.12847          !
!                                                            !
! [ref] I. Souza, N. Marzari and D. Vanderbilt,              !
!       Phys. Rev. B 65 035109 (2001)                        !
!       http://dx.doi.org/10.1103/PhysRevB.65.035109         !
!                                                            !
! [ref] N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza,    !
!       D. Vanderbilt, "Maximally localized Wannier          !
!       functions: theory and applications",                 !
!       Rev. Mod. Phys. 84, 1419 (2012)                      !
!       http://dx.doi.org/10.1103/RevModPhys.84.1419         !
!                                                            !
! For a full list of authors and contributors, please        !
! see the README file in the root directory of the           !
! distribution.                                              !
!                                                            !
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

program wannier
  !! The main Wannier90 program

#ifdef MPI08
  use mpi_f08
#endif
#ifdef MPI90
  use mpi
#endif

  use w90_library

  use w90_comms, only: w90_comm_type, comms_sync_err
  use w90_io, only: io_print_timings, io_commandline, io_date, prterr
  use w90_readwrite, only: w90_readwrite_write_header
  use w90_sitesym, only: sitesym_read
  use w90_error, only: w90_error_type, set_error_input

  implicit none

  integer, parameter :: dp = kind(0.d0)

  character(len=:), allocatable :: seedname, progname, cpstatus
  character(len=:), pointer :: restart
  complex(kind=dp), allocatable :: mloc(:, :, :, :)
  complex(kind=dp), allocatable :: u(:, :, :)
  complex(kind=dp), allocatable :: uopt(:, :, :)
  real(kind=dp), allocatable :: eigval(:, :)
  integer, allocatable :: distk(:)
  integer :: i, ctr
  integer :: mpisize, rank, ierr, nkl
  integer, pointer :: nb, nk, nw, nn
  integer :: stdout, stderr
  logical, pointer :: pp
  logical :: ld, lovlp, ldsnt, lwann, lplot, ltran, need_eigvals
  type(lib_common_type), target :: common_data
  type(w90_error_type), allocatable :: error
  character(len=9) :: cdate, ctime

  pp => common_data%w90_calculation%postproc_setup
  restart => common_data%w90_calculation%restart
  nw => common_data%num_wann
  nb => common_data%num_bands
  nk => common_data%num_kpts
  nn => common_data%kmesh_info%nntot

  progname = 'wannier90' ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=91442
  call io_commandline(progname, ld, pp, seedname)

#ifdef MPI
  call mpi_init(ierr)
  ! setup pplel decomp
  call mpi_comm_rank(mpi_comm_world, rank, ierr) ! the type of comm_world depends on interface used
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  call set_parallel_comms(common_data, mpi_comm_world)
#else
  rank = 0
  mpisize = 1
#endif

  ! open main output file
  open (newunit=stdout, file=seedname//'.wout', status="replace")
  ! open main error file
  open (newunit=stderr, file=seedname//'.werr', status="replace")

  call io_date(cdate, ctime)
  write (stderr, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

  call input_reader_special(common_data, seedname, stdout, stderr, ierr)
  if (ierr /= 0) stop

  call w90_input_reader(common_data, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! write useful info (includes jazzy header info)
  call input_print_details(common_data, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! test mpi error handling using "unlucky" input token
  if (rank == -common_data%print_output%timing_level) then
    call set_error_input(error, 'received unlucky_rank', common_data%comm)
  else
    call comms_sync_err(common_data%comm, error, 0) ! this is necessary since non-root may never enter an mpi collective if root has exited here
  endif
  if (allocated(error)) then ! applies (is t) for all ranks now
    call prterr(error, ierr, stdout, stderr, common_data%comm)
#ifdef MPI
    call mpi_finalize(ierr) ! let's be nice
#endif
    stop
  endif
  !!!!! end unlucky code

  ! setup kpoint distribution
  allocate (distk(nk))
  ctr = 0
  do i = 0, mpisize - 1
    nkl = nk/mpisize ! number of kpoints per rank
    if (mod(nk, mpisize) > i) nkl = nkl + 1
    if (nkl > 0) then
      distk(ctr + 1:ctr + nkl) = i
      ctr = ctr + nkl
    endif
  enddo

  ! copy distribution to library
  call set_kpoint_distribution(common_data, distk, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! special branch for writing nnkp file
  if (pp) then
    ! please only invoke on rank 0
    call write_kmesh(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    if (rank == 0) close (unit=stderr, status='delete')
#ifdef MPI
    call mpi_finalize(ierr)
#endif
    stop
  endif

  if (common_data%lsitesymmetry) then
    call sitesym_read(common_data%sitesym, nb, nk, nw, seedname, error, common_data%comm) ! (not a library call)
    if (allocated(error)) then
      write (stderr, *) 'failed to setup symmetry'
      deallocate (error)
      stop
    endif
  endif

  call w90_get_nn(common_data, nn, stdout, stderr, ierr)
  nkl = count(distk == rank) ! number of kpoints this rank
  write (*, *) 'rank, nw, nb, nk, nn, nk(rank): ', rank, nw, nb, nk, nn, nkl

  allocate (mloc(nb, nb, nn, nkl))
  allocate (u(nw, nw, nk))
  allocate (uopt(nb, nw, nk))

  call w90_set_m_local(common_data, mloc)  ! we don't need global m
  call w90_set_u_matrix(common_data, u)
  call w90_set_u_opt(common_data, uopt)

  uopt = 0.d0 ! required when not disentangling

! restart system
  lovlp = .true.
  ldsnt = .true.
  lwann = .true.
  lplot = .true.
  ltran = .false.

  if (restart == '') then
    if (rank == 0) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else
    cpstatus = ''
    call read_chkpt(common_data, cpstatus, stdout, stderr, ierr)
    if (ierr /= 0) stop

    if (restart == 'wannierise' .or. (restart == 'default' .and. cpstatus == 'postdis')) then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .true.
      lplot = .true.
      ltran = .false.
    elseif (restart == 'plot' .or. (restart == 'default' .and. cpstatus == 'postwann')) then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .false.
      lplot = .true.
      ltran = .false.
    elseif (restart == 'transport') then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from transport routines ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .false.
      lplot = .false.
      ltran = .true.
      !else
      ! illegitimate restart choice, should declaim the acceptable choices
    endif
  endif
  ltran = (ltran .or. common_data%w90_calculation%transport)
  ldsnt = (ldsnt .and. (nw < nb)) ! disentanglement only needed if space reduced

  ! circumstances where eigenvalues are needed are a little overcomplicated
  need_eigvals = .false.
  need_eigvals = common_data%w90_calculation%bands_plot
  need_eigvals = (need_eigvals .or. common_data%w90_calculation%fermi_surface_plot)
  need_eigvals = (need_eigvals .or. common_data%output_file%write_hr)
  need_eigvals = (need_eigvals .or. ldsnt) ! disentanglement anyway requires evals

  if (need_eigvals) then
    allocate (eigval(nb, nk))
    call read_eigvals(common_data, eigval, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call w90_set_eigval(common_data, eigval)
  endif

  ! ends setup

  if (lovlp) then
    call overlaps(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop

    if (.not. ldsnt) then
      call w90_project_overlap(common_data, stdout, stderr, ierr)
      if (ierr /= 0) stop
    endif
  endif

  if (ldsnt) then
    call w90_disentangle(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(common_data, 'postdis', stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lwann) then
    call w90_wannierise(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(common_data, 'postwann', stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lplot) then
    call w90_plot(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (ltran) then
    call w90_transport(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  call print_times(common_data, stdout)
  if (rank == 0) close (unit=stderr, status='delete')
#ifdef MPI
  call mpi_finalize(ierr)
#endif

contains
  subroutine input_reader_special(common_data, seedname, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_input, set_error_fatal, set_error_alloc
    use w90_readwrite, only: w90_readwrite_in_file, w90_readwrite_clean_infile
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read_special, w90_extra_io_type

    implicit none

    ! arguments
    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    logical :: disentanglement

    ierr = 0

    ! read data from .win file to internal string array
    call w90_readwrite_in_file(common_data%settings, seedname, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    call w90_wannier90_readwrite_read_special(common_data%settings, common_data%atom_data, &
                                              common_data%kmesh_input, common_data%kmesh_info, &
                                              common_data%kpt_latt, common_data%wann_control, &
                                              common_data%proj, common_data%proj_input, &
                                              common_data%select_proj, common_data%w90_system, &
                                              common_data%w90_calculation, &
                                              common_data%real_lattice, common_data%physics%bohr, &
                                              common_data%mp_grid, common_data%num_bands, &
                                              common_data%num_kpts, common_data%num_proj, &
                                              common_data%num_wann, common_data%gamma_only, &
                                              common_data%lhasproj, &
                                              common_data%use_bloch_phases, &
                                              common_data%dist_kpoints, istdout, error, &
                                              common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    common_data%seedname = seedname

    disentanglement = (common_data%num_bands > common_data%num_wann)

    if (disentanglement) then
      allocate (common_data%dis_manifold%ndimwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating ndimwin in input_reader_special() call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%nfirstwin(common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nfirstwin in input_reader_special() call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      allocate (common_data%dis_manifold%lwindow(common_data%num_bands, common_data%num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating lwindow in input_reader_special() call', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif
    allocate (common_data%wannier_data%centres(3, common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wannier_centres in input_reader_special() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%centres = 0.0_dp
    allocate (common_data%wannier_data%spreads(common_data%num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wannier_spreads in input_reader_special() call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%wannier_data%spreads = 0.0_dp

    ! remove any remaining acceptable keywords; anything that remains is an input error
    call w90_readwrite_clean_infile(common_data%settings, istdout, seedname, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    if (allocated(common_data%settings%in_data)) deallocate (common_data%settings%in_data)
  end subroutine input_reader_special

  subroutine input_print_details(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_readwrite, only: w90_readwrite_write_header
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write, w90_extra_io_type
    use w90_comms, only: mpisize

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error
    type(w90_extra_io_type) :: io_params ! what is this? fixme
    integer :: mpi_size

    ierr = 0

    mpi_size = mpisize(common_data%comm)

    ! write jazzy header info
    call w90_readwrite_write_header(common_data%physics%bohr_version_str, &
                                    common_data%physics%constants_version_str1, &
                                    common_data%physics%constants_version_str2, &
                                    mpi_size, istdout)

    ! write simulation details
    call w90_wannier90_readwrite_write(common_data%atom_data, common_data%band_plot, &
                                       common_data%dis_control, common_data%dis_spheres, &
                                       common_data%fermi_energy_list, &
                                       common_data%fermi_surface_data, common_data%kpt_latt, &
                                       common_data%output_file, common_data%wvfn_read, &
                                       common_data%wann_control, common_data%proj, &
                                       common_data%proj_input, common_data%real_space_ham, &
                                       common_data%select_proj, common_data%kpoint_path, &
                                       common_data%tran, common_data%print_output, &
                                       common_data%wannier_data, &
                                       common_data%wann_plot, io_params, &
                                       common_data%w90_calculation, common_data%real_lattice, &
                                       common_data%sitesym%symmetrize_eps, common_data%mp_grid, &
                                       common_data%num_bands, common_data%num_kpts, &
                                       common_data%num_proj, common_data%num_wann, &
                                       common_data%optimisation, .false., &
                                       common_data%gamma_only, common_data%lsitesymmetry, &
                                       common_data%w90_system%spinors, &
                                       common_data%use_bloch_phases, istdout)

    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine input_print_details

  subroutine write_kmesh(common_data, istdout, istderr, ierr)
    use w90_comms, only: mpirank, comms_sync_err
    use w90_error_base, only: w90_error_type
    use w90_kmesh, only: kmesh_get, kmesh_write

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (mpirank(common_data%comm) == 0) then
      call kmesh_write(common_data%exclude_bands, common_data%kmesh_info, &
                       common_data%select_proj%auto_projections, common_data%proj_input, &
                       common_data%print_output, common_data%kpt_latt, common_data%real_lattice, &
                       common_data%num_kpts, common_data%num_proj, common_data%calc_only_A, &
                       common_data%w90_system%spinors, common_data%seedname, common_data%timer)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif
    call comms_sync_err(common_data%comm, error, 0) ! this is necessary since non-root may never enter an mpi collective if root has exited here
  end subroutine write_kmesh

  subroutine overlaps(common_data, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal
    use w90_overlap, only: overlap_read

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    logical :: cp_pp = .false.
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (.not. common_data%setup_complete) then
      !call w90_create_kmesh(common_data, istdout, istderr, ierr)
      call set_error_fatal(error, 'kmesh is not setup before calling overlaps read routine', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
      if (ierr > 0) return
    endif

    if (common_data%num_bands > common_data%num_wann) then ! disentanglement case
      call overlap_read(common_data%kmesh_info, common_data%select_proj, common_data%sitesym, common_data%u_opt, &
                        common_data%m_matrix_local, common_data%num_bands, common_data%num_kpts, common_data%num_proj, &
                        common_data%num_wann, common_data%print_output, common_data%print_output%timing_level, cp_pp, &
                        common_data%gamma_only, common_data%lsitesymmetry, common_data%use_bloch_phases, &
                        common_data%seedname, istdout, common_data%timer, common_data%dist_kpoints, error, &
                        common_data%comm)
    else
      if (.not. associated(common_data%u_matrix)) then
        call set_error_fatal(error, 'u_matrix not associated at overlaps read routine', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
      if (.not. associated(common_data%m_matrix_local)) then
        call set_error_fatal(error, 'm_matrix_local not associated at overlaps read routine', common_data%comm)
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif

      call overlap_read(common_data%kmesh_info, common_data%select_proj, common_data%sitesym, &
                        common_data%u_matrix, common_data%m_matrix_local, common_data%num_bands, &
                        common_data%num_kpts, common_data%num_proj, common_data%num_wann, &
                        common_data%print_output, common_data%print_output%timing_level, cp_pp, &
                        common_data%gamma_only, common_data%lsitesymmetry, &
                        common_data%use_bloch_phases, common_data%seedname, istdout, &
                        common_data%timer, common_data%dist_kpoints, error, common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine overlaps

  subroutine read_eigvals(common_data, eigval, istdout, istderr, ierr)
    use w90_error, only: w90_error_type, set_error_fatal
    use w90_readwrite, only: w90_readwrite_read_eigvals

    implicit none

    ! arguments
    integer, intent(in) :: istdout, istderr
    real(kind=dp), intent(inout) :: eigval(:, :)
    type(lib_common_type), intent(inout) :: common_data
    integer, intent(out) :: ierr

    ! local vars
    type(w90_error_type), allocatable :: error
    logical :: eig_found

    ierr = 0

    if (size(eigval, 1) /= common_data%num_bands) then
      call set_error_fatal(error, 'eigval not dimensioned correctly (num_bands,num_kpts) in read_eigvals', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    elseif (size(eigval, 2) /= common_data%num_kpts) then
      call set_error_fatal(error, 'eigval not dimensioned correctly (num_bands,num_kpts) in read_eigvals', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    call w90_readwrite_read_eigvals(eig_found, eigval, common_data%num_bands, common_data%num_kpts, &
                                    istdout, common_data%seedname, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    else if (.not. eig_found) then
      call set_error_fatal(error, 'failed to read eigenvalues file in read_eigvals call', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine read_eigvals

  subroutine write_chkpt(common_data, label, istdout, istderr, ierr)
    use w90_comms, only: comms_reduce, mpirank
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_alloc, set_error_dealloc, set_error_fatal
    use w90_wannier90_readwrite, only: w90_wannier90_readwrite_write_chkpt

    implicit none

    ! arguments
    character(len=*), intent(in) :: label ! e.g. 'postdis' or 'postwann' after disentanglement, wannierisation
    integer, intent(in) :: istdout, istderr
    integer, intent(inout) :: ierr
    type(lib_common_type), target, intent(in) :: common_data

    ! local variables
    complex(kind=dp), allocatable :: u(:, :, :), uopt(:, :, :), m(:, :, :, :)
    integer, allocatable :: global_k(:)
    integer, pointer :: nw, nb, nk, nn
    integer :: rank, nkrank, ikg, ikl, istat
    type(w90_error_type), allocatable :: error

    ierr = 0
    rank = mpirank(common_data%comm)

    nb => common_data%num_bands
    nk => common_data%num_kpts
    nn => common_data%kmesh_info%nntot
    nw => common_data%num_wann

    if (.not. associated(common_data%u_opt)) then
      call set_error_fatal(error, 'u_opt not set for write_chkpt call', common_data%comm)
    else if (.not. associated(common_data%u_matrix)) then
      call set_error_fatal(error, 'u_matrix not set for write_chkpt call', common_data%comm)
    else if (.not. associated(common_data%m_matrix_local)) then
      call set_error_fatal(error, 'm_matrix_local not set for write_chkpt call', common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    nkrank = count(common_data%dist_kpoints == rank)
    allocate (global_k(nkrank), stat=istat)
    if (istat /= 0) then
      call set_error_alloc(error, 'Error allocating global_k in write_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    global_k = huge(1); ikl = 1
    do ikg = 1, nk
      if (rank == common_data%dist_kpoints(ikg)) then
        global_k(ikl) = ikg
        ikl = ikl + 1
      endif
    enddo

    ! allocating and partially assigning the full matrix on all ranks and reducing is a terrible idea
    ! alternatively, allocate on root and use point-to-point
    ! or, if required only for checkpoint file writing, then use mpi-io (but needs to be ordered io, alas)
    ! or, even better, use parallel hdf5. JJ Nov 22
    allocate (u(nw, nw, nk), stat=istat) ! all kpts
    if (istat /= 0) call set_error_alloc(error, 'Error allocating u in write_chkpt', common_data%comm)
    allocate (uopt(nb, nw, nk), stat=istat) ! all kpts
    if (istat /= 0) call set_error_alloc(error, 'Error allocating uopt in write_chkpt', common_data%comm)
    allocate (m(nw, nw, nn, nk), stat=istat) ! all kpts
    if (istat /= 0) call set_error_alloc(error, 'Error allocating m in write_chkpt', common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    u(:, :, :) = 0.d0
    uopt(:, :, :) = 0.d0
    m(:, :, :, :) = 0.d0

    do ikl = 1, nkrank
      ikg = global_k(ikl)
      u(:, :, ikg) = common_data%u_matrix(:, :, ikl)
      uopt(:, :, ikg) = common_data%u_opt(:, :, ikl)
      m(:, :, :, ikg) = common_data%m_matrix_local(1:nw, 1:nw, :, ikl)
    enddo

    call comms_reduce(u(1, 1, 1), nw*nw*nk, 'SUM', error, common_data%comm)
    call comms_reduce(uopt(1, 1, 1), nb*nw*nk, 'SUM', error, common_data%comm)
    call comms_reduce(m(1, 1, 1, 1), nw*nw*nn*nk, 'SUM', error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (rank == 0) then
      call w90_wannier90_readwrite_write_chkpt(label, common_data%exclude_bands, &
                                               common_data%wannier_data, common_data%kmesh_info, &
                                               common_data%kpt_latt, nk, common_data%dis_manifold, &
                                               nb, nw, u, uopt, m, common_data%mp_grid, &
                                               common_data%real_lattice, &
                                               common_data%omega%invariant, &
                                               common_data%have_disentangled, &
                                               common_data%print_output%iprint, istdout, &
                                               common_data%seedname)
    endif

    deallocate (u, stat=istat)
    if (istat /= 0) then
      call set_error_dealloc(error, 'Error deallocating u in write_chkpt', common_data%comm)
    endif
    deallocate (uopt, stat=istat)
    if (istat /= 0) then
      call set_error_dealloc(error, 'Error deallocating uopt in write_chkpt', common_data%comm)
    endif
    deallocate (m, stat=istat)
    if (istat /= 0) then
      call set_error_dealloc(error, 'Error deallocating m in write_chkpt', common_data%comm)
    endif
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine write_chkpt

  subroutine read_chkpt(common_data, checkpoint, istdout, istderr, ierr)
    use w90_comms, only: mpirank
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_alloc, set_error_dealloc
    use w90_readwrite, only: w90_readwrite_read_chkpt, w90_readwrite_chkpt_dist

    implicit none

    ! arguments
    character(len=*), intent(out) :: checkpoint
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), target, intent(inout) :: common_data

    ! local variables
    complex(kind=dp), allocatable :: m(:, :, :, :)
    integer, pointer :: nw, nb, nk, nn
    integer :: rank, nexclude = 0, istat
    logical :: ispostw90 = .false. ! ispostw90 is used to print a different error message in case the chk file is missing (did you run w90 first?)
    type(w90_error_type), allocatable :: error

    ierr = 0
    rank = mpirank(common_data%comm)

    nb => common_data%num_bands
    nk => common_data%num_kpts
    nn => common_data%kmesh_info%nntot
    nw => common_data%num_wann

    ! allocating and partially assigning the full matrix on all ranks and reducing is a terrible idea at scale
    ! alternatively, allocate on root and use point-to-point
    ! or, if required only for checkpoint file writing, then use mpi-io (but needs to be ordered io, alas)
    ! or, even better, use parallel hdf5
    ! fixme, check allocation status of u, uopt?
    allocate (m(nw, nw, nn, nk), stat=istat) ! all kpts
    if (istat /= 0) then
      call set_error_alloc(error, 'Error allocating m in read_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    if (rank == 0) then
      if (allocated(common_data%exclude_bands)) nexclude = size(common_data%exclude_bands)

      call w90_readwrite_read_chkpt(common_data%dis_manifold, common_data%exclude_bands, &
                                    common_data%kmesh_info, common_data%kpt_latt, &
                                    common_data%wannier_data, m, common_data%u_matrix, &
                                    common_data%u_opt, common_data%real_lattice, &
                                    common_data%omega%invariant, common_data%mp_grid, nb, &
                                    nexclude, nk, nw, checkpoint, common_data%have_disentangled, &
                                    ispostw90, common_data%seedname, istdout, error, &
                                    common_data%comm)
      if (allocated(error)) then
        call prterr(error, ierr, istdout, istderr, common_data%comm)
        return
      endif
    endif

    ! scatter from m_matrix to m_matrix_local (cf overlap_read)
    call w90_readwrite_chkpt_dist(common_data%dis_manifold, common_data%wannier_data, &
                                  common_data%u_matrix, common_data%u_opt, m, &
                                  common_data%m_matrix_local, common_data%omega%invariant, &
                                  nb, nk, nw, nn, checkpoint, common_data%have_disentangled, &
                                  common_data%dist_kpoints, error, common_data%comm)
    if (allocated(error)) then
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif

    deallocate (m, stat=istat)
    if (istat /= 0) then
      call set_error_alloc(error, 'Error deallocating m in read_chkpt', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
  end subroutine read_chkpt

  subroutine print_times(common_data, istdout)
    use w90_io, only: io_print_timings
    implicit none

    type(lib_common_type), intent(in) :: common_data
    integer, intent(in) :: istdout

    if (common_data%print_output%iprint > 0) call io_print_timings(common_data%timer, istdout)
  end subroutine print_times

  subroutine set_parallel_comms(common_data, comm)
#ifdef MPI08
    use mpi_f08
#endif
    implicit none
    type(lib_common_type), intent(inout) :: common_data
#ifdef MPI08
    type(mpi_comm), intent(in) :: comm
#else
    integer, intent(in) :: comm
#endif
    common_data%comm%comm = comm
  end subroutine set_parallel_comms

  subroutine set_kpoint_distribution(common_data, dist, istdout, istderr, ierr)
    use w90_error_base, only: w90_error_type
    use w90_error, only: set_error_fatal

    implicit none

    ! arguments
    integer, intent(in) :: istderr, istdout, dist(:)
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data

    ! local variables
    type(w90_error_type), allocatable :: error

    ierr = 0

    if (size(dist) < 1) then
      call set_error_fatal(error, 'Error in k-point distribution, mpisize < 1', common_data%comm)
      call prterr(error, ierr, istdout, istderr, common_data%comm)
      return
    endif
    common_data%dist_kpoints = dist
  end subroutine set_kpoint_distribution

end program wannier
