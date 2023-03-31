program libv2

#ifdef MPI08
  use mpi_f08
#endif
#ifdef MPI90
  use mpi
#endif

  use w90_helper_types

  use w90_comms, only: w90_comm_type
  use w90_io, only: io_print_timings
  use w90_readwrite, only: w90_readwrite_write_header
  use w90_sitesym, only: sitesym_read
  use w90_error, only: w90_error_type

  implicit none

  character(len=100) :: seedname
  character(len=:), allocatable :: fn, cpstatus
  character(len=:), pointer :: restart
  complex(kind=dp), allocatable :: a(:, :, :)
  !complex(kind=dp), allocatable :: m(:,:,:,:)
  complex(kind=dp), allocatable :: mloc(:, :, :, :)
  complex(kind=dp), allocatable :: morig(:, :, :, :)
  complex(kind=dp), allocatable :: u(:, :, :)
  complex(kind=dp), allocatable :: uopt(:, :, :)
  real(kind=dp), allocatable :: eigval(:, :)
  integer, allocatable :: distk(:)
  integer :: length, len2, i
  integer :: mpisize, rank, ierr, nkl
  integer, pointer :: nb, nk, nw, nn
  integer :: stdout, stderr
  logical, pointer :: pp
  logical :: lovlp, ldsnt, lwann, lplot, ltran, need_eigvals
  type(lib_global_type), target :: w90main
  type(lib_w90_type), target :: w90dat
  type(w90_comm_type) :: comm
  type(w90_error_type), allocatable :: error

  pp => w90dat%w90_calculation%postproc_setup
  restart => w90dat%w90_calculation%restart
  nw => w90main%num_wann
  nb => w90main%num_bands
  nk => w90main%num_kpts
  nn => w90main%kmesh_info%nntot

  ! get seedname and pp flag, if present
  call get_command_argument(1, seedname, length, ierr)
  if (ierr /= 0) then
    write (*, *) 'failed to parse seedname'
    stop
  endif
  if (seedname == '-pp') then
    pp = .true.
    call get_command_argument(2, seedname, length, ierr)
    if (ierr /= 0) then
      write (*, *) 'failed to parse seedname'
      stop
    endif
  endif
  do i = 1, length
    if (seedname(i:i) == '.') exit
    len2 = i
  enddo
  fn = trim(seedname(1:len2))
  ! end get seedname

#ifdef MPI
  comm%comm = mpi_comm_world
  call mpi_init(ierr)
  ! setup pplel decomp
  call mpi_comm_rank(comm%comm, rank, ierr)
  call mpi_comm_size(comm%comm, mpisize, ierr)
  call set_parallel_comms(w90main, comm%comm)
#else
  rank = 0
  mpisize = 1
#endif

  ! open main output file
  open (newunit=stdout, file=fn//'.wout', status="replace")
  ! open main error file
  open (newunit=stderr, file=fn//'.werr', status="replace")

  call input_reader(w90main, w90dat, fn, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! write jazzy header info
  call w90_readwrite_write_header(w90main%physics%bohr_version_str, &
                                  w90main%physics%constants_version_str1, &
                                  w90main%physics%constants_version_str2, stdout) ! (not a library call)

  ! setup k mesh
  if (.not. w90main%kmesh_info%explicit_nnkpts) then
    call create_kmesh(w90main, stdout, stderr, ierr)
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
  call set_kpoint_distribution(w90main, distk, stdout, stderr, ierr)
  if (ierr /= 0) stop

  nkl = count(distk == rank) ! number of kpoints this rank
  write (*, *) 'rank, nw, nb, nk, nn, nk(rank): ', rank, nw, nb, nk, nn, nkl

  ! special branch for writing nnkp file
  if (pp) then
    ! please only invoke on rank 0
    call write_kmesh(w90main, w90dat, fn, stdout, stderr, ierr)
    if (ierr /= 0) stop
    if (rank == 0) close (unit=stderr, status='delete')
#ifdef MPI
    call mpi_finalize(ierr)
#endif
    stop
  endif

  if (w90dat%lsitesymmetry) then
    call sitesym_read(w90dat%sitesym, nb, nk, nw, fn, error, comm) ! (not a library call)
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

  call set_m_matrix_local(w90dat, mloc)
  !call set_m_matrix(w90dat, m)  ! we don't need global m
  call set_u_matrix(w90main, u)
  call set_u_opt(w90main, uopt)

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
    call read_chkpt(w90main, w90dat, cpstatus, fn, stdout, stderr, ierr)
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
  ltran = (ltran .or. w90dat%w90_calculation%transport)
  ldsnt = (ldsnt .and. (nw < nb)) ! disentanglement only needed if space reduced

  if (ldsnt) then
    ! setup arrays before overlap reads
    allocate (a(nb, nw, nk))
    allocate (morig(nb, nb, nn, nkl))
    call set_a_matrix(w90dat, a)
    call set_m_orig(w90dat, morig)
  endif

  ! circumstances where eigenvalues are needed are a little overcomplicated
  need_eigvals = .false.
  need_eigvals = w90dat%w90_calculation%bands_plot
  need_eigvals = (need_eigvals .or. w90dat%w90_calculation%fermi_surface_plot)
  need_eigvals = (need_eigvals .or. w90dat%output_file%write_hr)
  need_eigvals = (need_eigvals .or. ldsnt) ! disentanglement anyway requires evals

  if (need_eigvals) then
    allocate (eigval(nb, nk))
    call read_eigvals(w90main, eigval, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call set_eigval(w90main, eigval)
  endif

  ! ends setup

  if (lovlp) then
    call overlaps(w90main, w90dat, stdout, stderr, ierr)
    if (ierr /= 0) stop

    if (.not. ldsnt) then
      call projovlp(w90main, w90dat, stdout, stderr, ierr)
      if (ierr /= 0) stop
    endif
  endif

  if (ldsnt) then
    call disentangle(w90main, w90dat, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(w90main, w90dat, 'postdis', fn, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lwann) then
    call wannierise(w90main, w90dat, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(w90main, w90dat, 'postwann', fn, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lplot) then
    call plot_files(w90main, w90dat, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (ltran) then
    call transport(w90main, w90dat, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  call print_times(w90main, stdout)
  if (rank == 0) close (unit=stderr, status='delete')
#ifdef MPI
  call mpi_finalize(ierr)
#endif
end program libv2
