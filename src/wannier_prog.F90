program libv2

#ifdef MPI08
  use mpi_f08
#endif
#ifdef MPI90
  use mpi
#endif

  use w90_library

  use w90_comms, only: w90_comm_type, comms_sync_err
  use w90_io, only: io_print_timings, io_commandline
  use w90_readwrite, only: w90_readwrite_write_header
  use w90_sitesym, only: sitesym_read
  use w90_error, only: w90_error_type

  implicit none

  character(len=:), allocatable :: seedname, progname, cpstatus
  character(len=:), pointer :: restart
  complex(kind=dp), allocatable :: a(:, :, :)
  !complex(kind=dp), allocatable :: m(:,:,:,:)
  complex(kind=dp), allocatable :: mloc(:, :, :, :)
  complex(kind=dp), allocatable :: morig(:, :, :, :)
  complex(kind=dp), allocatable :: u(:, :, :)
  complex(kind=dp), allocatable :: uopt(:, :, :)
  real(kind=dp), allocatable :: eigval(:, :)
  integer, allocatable :: distk(:)
  integer :: i
  integer :: mpisize, rank, ierr, nkl
  integer, pointer :: nb, nk, nw, nn
  integer :: stdout, stderr
  logical, pointer :: pp
  logical :: ld, lovlp, ldsnt, lwann, lplot, ltran, need_eigvals
  type(lib_common_type), target :: common_data
  type(lib_wannier_type), target :: wannier_data
  type(w90_comm_type) :: comm
  type(w90_error_type), allocatable :: error

  pp => wannier_data%w90_calculation%postproc_setup
  restart => wannier_data%w90_calculation%restart
  nw => common_data%num_wann
  nb => common_data%num_bands
  nk => common_data%num_kpts
  nn => common_data%kmesh_info%nntot

  progname = 'wannier90' ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=91442
  call io_commandline(progname, ld, pp, seedname, stderr)

#ifdef MPI
  comm%comm = mpi_comm_world
  call mpi_init(ierr)
  ! setup pplel decomp
  call mpi_comm_rank(comm%comm, rank, ierr)
  call mpi_comm_size(comm%comm, mpisize, ierr)
  call set_parallel_comms(common_data, comm%comm)
#else
  rank = 0
  mpisize = 1
#endif

  ! open main output file
  open (newunit=stdout, file=seedname//'.wout', status="replace")
  ! open main error file
  open (newunit=stderr, file=seedname//'.werr', status="replace")

  call input_reader_special(common_data, wannier_data, seedname, stdout, stderr, ierr)
  if (ierr /= 0) stop

  call input_reader(common_data, wannier_data, seedname, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! test mpi error handling using "unlucky" input token
  ! this machinery used to sit in w90_wannier90_readwrite_dist
  ! but that routine is obsolete if input file is read on all ranks
  ! fixme, this should be moved to wannier90 main routine (definately doesn't belong here)
  if (common_data%print_output%timing_level < 0 &
      .and. rank == abs(common_data%print_output%timing_level)) then
!    call set_error_input(error, 'received unlucky_rank', common_data%comm)
  else
    call comms_sync_err(common_data%comm, error, 0) ! this is necessary since non-root may never enter an mpi collective if root has exited here
  endif
  if (allocated(error)) then ! applies (is t) for all ranks now
!    call prterr(error, ierr, istdout, istderr, common_data%comm)
!    return
  endif
  !!!!! end unlucky code

  ! write jazzy header info
  call w90_readwrite_write_header(common_data%physics%bohr_version_str, &
                                  common_data%physics%constants_version_str1, &
                                  common_data%physics%constants_version_str2, stdout) ! (not a library call)

  ! setup k mesh
  if (.not. common_data%kmesh_info%explicit_nnkpts) then
    call create_kmesh(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  ! setup kpoint distribution
  allocate (distk(nk))
  nkl = nk/mpisize ! number of kpoints per rank
  if (mod(nk, mpisize) > 0) nkl = nkl + 1
  do i = 1, nk
    distk(i) = (i - 1)/nkl ! contiguous blocks with potentially fewer processes on last rank
  enddo
  ! copy distribution to library
  call set_kpoint_distribution(common_data, distk, stdout, stderr, ierr)
  if (ierr /= 0) stop

  nkl = count(distk == rank) ! number of kpoints this rank
  write (*, *) 'rank, nw, nb, nk, nn, nk(rank): ', rank, nw, nb, nk, nn, nkl

  ! special branch for writing nnkp file
  if (pp) then
    ! please only invoke on rank 0
    call write_kmesh(common_data, wannier_data, seedname, stdout, stderr, ierr)
    if (ierr /= 0) stop
    if (rank == 0) close (unit=stderr, status='delete')
#ifdef MPI
    call mpi_finalize(ierr)
#endif
    stop
  endif

  if (wannier_data%lsitesymmetry) then
    call sitesym_read(wannier_data%sitesym, nb, nk, nw, seedname, error, comm) ! (not a library call)
    if (allocated(error)) then
      write (stderr, *) 'failed to setup symmetry'
      deallocate (error)
      stop
    endif
  endif

  allocate (mloc(nw, nw, nn, nkl))
  !allocate (m(nw, nw, nn, nk)) ! we don't need global m
  allocate (u(nw, nw, nk))
  allocate (uopt(nb, nw, nk))

  call set_m_matrix_local(wannier_data, mloc)
  !call set_m_matrix(wannier_data, m)  ! we don't need global m
  call set_u_matrix(common_data, u)
  call set_u_opt(common_data, uopt)

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
    call read_chkpt(common_data, wannier_data, cpstatus, seedname, stdout, stderr, ierr)
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
  ltran = (ltran .or. wannier_data%w90_calculation%transport)
  ldsnt = (ldsnt .and. (nw < nb)) ! disentanglement only needed if space reduced

  if (ldsnt) then
    ! setup arrays before overlap reads
    allocate (a(nb, nw, nk))
    allocate (morig(nb, nb, nn, nkl))
    call set_a_matrix(wannier_data, a)
    call set_m_orig(wannier_data, morig)
  endif

  ! circumstances where eigenvalues are needed are a little overcomplicated
  need_eigvals = .false.
  need_eigvals = wannier_data%w90_calculation%bands_plot
  need_eigvals = (need_eigvals .or. wannier_data%w90_calculation%fermi_surface_plot)
  need_eigvals = (need_eigvals .or. wannier_data%output_file%write_hr)
  need_eigvals = (need_eigvals .or. ldsnt) ! disentanglement anyway requires evals

  if (need_eigvals) then
    allocate (eigval(nb, nk))
    call read_eigvals(common_data, eigval, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call set_eigval(common_data, eigval)
  endif

  ! ends setup

  if (lovlp) then
    call overlaps(common_data, wannier_data, stdout, stderr, ierr)
    if (ierr /= 0) stop

    if (.not. ldsnt) then
      call projovlp(common_data, wannier_data, stdout, stderr, ierr)
      if (ierr /= 0) stop
    endif
  endif

  if (ldsnt) then
    call disentangle(common_data, wannier_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(common_data, wannier_data, 'postdis', seedname, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lwann) then
    call wannierise(common_data, wannier_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(common_data, wannier_data, 'postwann', seedname, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lplot) then
    call plot_files(common_data, wannier_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (ltran) then
    call transport(common_data, wannier_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  call print_times(common_data, stdout)
  if (rank == 0) close (unit=stderr, status='delete')
#ifdef MPI
  call mpi_finalize(ierr)
#endif
end program libv2
