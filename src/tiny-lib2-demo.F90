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
  use w90_readwrite, only: w90_readwrite_write_header, w90_readwrite_read_eigvals
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
  integer, allocatable :: distk(:)
  integer :: length, len2, i
  integer :: mpisize, rank, ierr, nkl
  integer, pointer :: nb, nk, nw, nn
  integer :: stdout, stderr
  logical, pointer :: pp
  logical :: lovlp, ldsnt, lwann, lplot, ltran
  type(lib_global_type), target :: w90main
  type(lib_w90_type), target :: w90dat
  type(w90_comm_type) :: comm
  type(w90_error_type), allocatable :: error

  !jj temporary home for read_eigvals
  logical :: eig_found

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
#endif

  call get_fortran_stdout(stdout)
  call get_fortran_stderr(stderr)
  call input_reader(w90main, w90dat, fn, stdout, stderr, ierr, comm)

  ! special branch for writing nnkp file
  if (pp) then
    call write_kmesh(w90main, w90dat, fn, stdout, stderr, ierr, comm)
#ifdef MPI
    call mpi_finalize(ierr)
#endif
    stop
  endif

  ! open main output file
  open (newunit=stdout, file=fn//'.wout')
  ! open main error file
  open (newunit=stderr, file=fn//'.werr')

  ! write jazzy header info
  call w90_readwrite_write_header(w90main%physics%bohr_version_str, &
                                  w90main%physics%constants_version_str1, &
                                  w90main%physics%constants_version_str2, stdout) ! (not a library call)

#ifdef MPI
  ! setup pplel decomp
  call mpi_comm_rank(comm%comm, rank, ierr)
  call mpi_comm_size(comm%comm, mpisize, ierr)
#else
  rank = 0
  mpisize = 1
#endif

  allocate (distk(nk))
  nkl = nk/mpisize ! number of kpoints per rank
  if (mod(nk, mpisize) > 0) nkl = nkl + 1
  do i = 1, nk
    distk(i) = (i - 1)/nkl ! contiguous blocks with potentially fewer processes on last rank
  enddo
  call set_kpoint_distribution(w90main, distk)
  nkl = count(distk == rank) ! number of kpoints this rank

  ! setup k mesh
  call create_kmesh(w90main, stdout, stderr, ierr, comm)
  write (*, *) 'rank, nw, nb, nk, nn, nk(rank): ', rank, nw, nb, nk, nn, nkl

  if (w90dat%lsitesymmetry) then
    call sitesym_read(w90dat%sitesym, nb, nk, nw, fn, error, comm) ! (not a library call)
    if (allocated(error)) then
      ierr = error%code
      deallocate (error)
      error stop
    endif
  endif

  if (.not. (w90dat%w90_calculation%transport .and. w90dat%tran%read_ht)) then
    call w90_readwrite_read_eigvals(.false., .false., .false., &
                                    w90dat%w90_calculation%bands_plot .or. w90dat%w90_calculation%fermi_surface_plot .or. &
                                    w90dat%output_file%write_hr, nw < nb, eig_found, &
                                    w90main%eigval, w90dat%w90_calculation%postproc_setup, nb, &
                                    nk, stdout, fn, error, comm)
    if (allocated(error)) then
      ierr = error%code
      deallocate (error)
      error stop
    endif
    ! test for equality just a hack for now to avoid overwriting an assigned variable
    if (eig_found .and. w90main%dis_manifold%win_min == -huge(0.0_dp)) w90main%dis_manifold%win_min = minval(w90main%eigval)
    if (eig_found .and. w90main%dis_manifold%win_max == huge(0.0_dp)) w90main%dis_manifold%win_max = maxval(w90main%eigval)
  endif

  if (nw < nb) then ! disentanglement reqired
    allocate (a(nb, nw, nk))
    allocate (morig(nb, nb, nn, nkl))
    call set_a_matrix(w90dat, a)
    call set_m_orig(w90dat, morig)
  endif

  allocate (mloc(nw, nw, nn, nkl))
  !allocate (m(nw, nw, nn, nk)) ! we don't need global m
  allocate (u(nw, nw, nk))
  allocate (uopt(nb, nw, nk))

  call set_m_matrix_local(w90dat, mloc)
  !call set_m_matrix(w90dat, m)  ! we don't need global m
  call set_u_matrix(w90main, u)
  call set_u_opt(w90main, uopt)

! restart system
  lovlp = .true.
  ldsnt = .true.
  lwann = .true.
  lplot = .true.
  ltran = .true.

  if (restart == '') then
    if (rank == 0) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else
    call read_chkpt(w90main, w90dat, cpstatus, fn, stdout, stderr, ierr, comm)
    if (restart == 'wannierise' .or. (restart == 'default' .and. cpstatus == 'postdis')) then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .true.
      lplot = .true.
      ltran = .true.
    elseif (restart == 'plot' .or. (restart == 'default' .and. cpstatus == 'postwann')) then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .false.
      lplot = .true.
      ltran = .true.
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
! end restart system

  if (lovlp) then
    call overlaps(w90main, w90dat, stdout, stderr, ierr, comm)
    if (ierr /= 0) error stop
  endif

  if (ldsnt) then
    if (nw < nb) then ! disentanglement reqired
      call disentangle(w90main, w90dat, stdout, stderr, ierr, comm)
      if (ierr /= 0) error stop
      call write_chkpt(w90main, w90dat, 'postdis', fn, stdout, stderr, ierr, comm)
      if (ierr /= 0) error stop
    endif
  endif

  if (lwann) then
    call wannierise(w90main, w90dat, stdout, stderr, ierr, comm)
    if (ierr /= 0) error stop
    call write_chkpt(w90main, w90dat, 'postwann', fn, stdout, stderr, ierr, comm)
    if (ierr /= 0) error stop
  endif

  if (lplot) then
    call plot_files(w90main, w90dat, stdout, stderr, ierr, comm)
    if (ierr /= 0) error stop
  endif

  call print_times(w90main, stdout)
  if (rank == 0) close (unit=stderr, status='delete')
#ifdef MPI
  call mpi_finalize(ierr)
#endif
end program libv2
