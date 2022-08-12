program libv2

  use mpi_f08
  use w90_comms, only: w90comm_type
  use w90_helper_types
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
  integer :: mpisize, rank, ierr, stat
  integer :: nb, nk, nw, nn, ik, length, len2, i
  integer :: stdout
  logical :: pp = .false.
  type(lib_global_type) :: w90main
  type(lib_w90_type) :: w90dat
  type(w90comm_type) :: comm 
  type(w90_error_type), allocatable :: error

  ! get seedname and pp flag, if present
  call get_command_argument(1, seedname, length, stat)
  if (stat /= 0) then
    write(*,*)'failed to parse seedname'
    stop
  endif
  if (seedname == '-pp') then
    pp = .true.
    w90dat%w90_calculation%postproc_setup = .true.
    call get_command_argument(2, seedname, length, stat)
    ! check stat
  endif
  do i = 1, length
    if (seedname(i:i) == '.') exit
    len2 = i
  enddo
  fn = trim(seedname(1:len2))
  write(*,*) fn//'.ext'
  ! end get seedname

  comm%comm= mpi_comm_world
  call mpi_init(ierr)
  call mpi_comm_rank(comm%comm, rank, ierr)
  call mpi_comm_size(comm%comm, mpisize, ierr)

  if (mpisize > 1) then
    write(*,*)'not yet parallel'
    stop
  endif

  call input_reader(w90main, w90dat, fn, 6, stat, comm)

  if (pp) then
    call write_kmesh(w90main, w90dat, fn, stat, comm)
    stop
  endif

  ! open main output file
  open(newunit=stdout, file=fn//'.wout')

  call w90_readwrite_write_header(w90main%physics%bohr_version_str, &
                                  w90main%physics%constants_version_str1, &
                                  w90main%physics%constants_version_str2, stdout)

  nw = w90main%num_wann
  nb = w90main%num_bands
  nk = w90main%num_kpts

  ! setup pplel decomp
  allocate(counts(0:mpisize-1))
  allocate(displs(0:mpisize-1))
  allocate(distk(nk))
  counts(rank) = nk
  displs(rank) = 0
  do ik = 1, nk
    distk(ik) = 0
  end do
  call set_kpoint_block(w90main, counts, displs)
  call set_kpoint_distribution(w90main, distk)
  ! end setup pplel decomp

  call create_kmesh(w90main, stdout, stat, comm)
  nn = w90main%kmesh_info%nntot
  write(*,*)'nw, nb, nk, nn: ', nw, nb, nk, nn


  if (w90dat%lsitesymmetry) then
    write(*,*) "w sym"
    call sitesym_read(w90dat%sitesym, nb, nk, nw, fn, error, comm)
  endif

  allocate(a(nb, nw, nk))
  allocate(mloc(nw, nw, nn, nk))
  allocate(m(nw, nw, nn, nk))
  allocate(morig(nb, nb, nn, nk))
  allocate(u(nw, nw, nk))
  allocate(uopt(nb, nw, nk))

  call set_a_matrix(w90dat, a)
  call set_m_matrix_local(w90dat, mloc)
  call set_m_matrix(w90dat, m)
  call set_m_orig(w90dat, morig)
  call set_u_matrix(w90main, u)
  call set_u_opt(w90main, uopt)

  call overlaps(w90main, w90dat, stdout, stat, comm)

  if (nw == nb ) then 
    uopt = u
    mloc = m
    call wannierise(w90main, w90dat, stdout, stat, comm)
  else
    call disentangle(w90main, w90dat, stdout, stat, comm)
    call wannierise(w90main, w90dat, stdout, stat, comm)
  endif

  call plot_files(w90main, w90dat, stdout, stat, comm)
  call print_times(w90main, stdout)
end program libv2
