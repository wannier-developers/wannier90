program libv2

  use mpi_f08
  use w90_helper_types

  use w90_comms, only: w90comm_type, comms_array_split
  use w90_io, only: io_print_timings
  use w90_readwrite, only: w90_readwrite_write_header
  use w90_sitesym, only: sitesym_read
  use w90_error, only: w90_error_type

  implicit none

  character(len=100) :: seedname
  character(len=:), allocatable :: fn
  complex(kind=dp), allocatable :: a(:,:,:)
  complex(kind=dp), allocatable :: m(:,:,:,:)
  complex(kind=dp), allocatable :: mloc(:,:,:,:)
  complex(kind=dp), allocatable :: morig(:,:,:,:)
  complex(kind=dp), allocatable :: u(:,:,:)
  complex(kind=dp), allocatable :: uopt(:,:,:)
  integer, allocatable :: counts(:), displs(:), distk(:)
  integer :: length, len2, i, j
  integer :: mpisize, rank, ierr, stat
  integer, pointer :: nb, nk, nw, nn
  integer :: stdout, stderr
  logical, pointer :: pp 
  type(lib_global_type), target :: w90main
  type(lib_w90_type), target :: w90dat
  type(w90comm_type) :: comm 
  type(w90_error_type), allocatable :: error

  pp => w90dat%w90_calculation%postproc_setup
  nw => w90main%num_wann
  nb => w90main%num_bands
  nk => w90main%num_kpts
  nn => w90main%kmesh_info%nntot

  ! get seedname and pp flag, if present
  call get_command_argument(1, seedname, length, stat)
  if (stat /= 0) then
    write(*,*)'failed to parse seedname'
    stop
  endif
  if (seedname == '-pp') then
    pp = .true.
    call get_command_argument(2, seedname, length, stat)
    if (stat /= 0) then
      write(*,*)'failed to parse seedname'
      stop
    endif
  endif
  do i = 1, length
    if (seedname(i:i) == '.') exit
    len2 = i
  enddo
  fn = trim(seedname(1:len2))
  ! end get seedname

  comm%comm= mpi_comm_world
  call mpi_init(ierr)

  call input_reader(w90main, w90dat, fn, 6, stat, comm)

  ! special branch for writing nnkp file
  if (pp) then
    call write_kmesh(w90main, w90dat, fn, 6, stat, comm)
    call mpi_finalize()
    stop
  endif

  ! open main output file
  open(newunit=stdout, file=fn//'.wout')
  ! open main error file 
  !open(newunit=stderr, file=fn//'.werr') ! fixme: need to add stderr to lib calls

  ! write jazzy header info
  call w90_readwrite_write_header(w90main%physics%bohr_version_str, &
                                  w90main%physics%constants_version_str1, &
                                  w90main%physics%constants_version_str2, stdout) ! (not a library call)

  ! setup pplel decomp
  call mpi_comm_rank(comm%comm, rank, ierr)
  call mpi_comm_size(comm%comm, mpisize, ierr)
  allocate(counts(0:mpisize-1))
  allocate(displs(0:mpisize-1))
  allocate(distk(nk))
  call comms_array_split(nk, counts, displs, comm) ! (not a library call) fixme
  distk(:) = -1
  do i = 0, mpisize - 1
    do j = 1, counts(i)
      distk(displs(i) + j) = i
    enddo
  enddo
  call set_kpoint_block(w90main, counts, displs) 
  call set_kpoint_distribution(w90main, distk)
  ! end setup pplel decomp

  call create_kmesh(w90main, stdout, stat, comm)
  write(*,*)'rank, nw, nb, nk, nn, nk(rank): ', rank, nw, nb, nk, nn, counts(rank)

  if (w90dat%lsitesymmetry) then
    call sitesym_read(w90dat%sitesym, nb, nk, nw, fn, error, comm) ! (not a library call)
  endif

  if (nw < nb ) then
    allocate(a(nb, nw, nk))
    allocate(morig(nb, nb, nn, counts(rank)))
    call set_a_matrix(w90dat, a)
    call set_m_orig(w90dat, morig)
  endif

  allocate(mloc(nw, nw, nn, counts(rank)))
!  allocate(m(nw, nw, nn, nk)) ! we don't need global m
  allocate(u(nw, nw, nk))
  allocate(uopt(nb, nw, nk))

  call set_m_matrix_local(w90dat, mloc)
!  call set_m_matrix(w90dat, m)  ! we don't need global m
  call set_u_matrix(w90main, u)
  call set_u_opt(w90main, uopt)

  call overlaps(w90main, w90dat, stdout, stat, comm)

  if (nw < nb) then ! disentanglement reqired
    call disentangle(w90main, w90dat, stdout, stat, comm)
  endif
  call wannierise(w90main, w90dat, stdout, stat, comm)
  call plot_files(w90main, w90dat, stdout, stat, comm)

  call print_times(w90main, stdout)
  call mpi_finalize()
end program libv2

