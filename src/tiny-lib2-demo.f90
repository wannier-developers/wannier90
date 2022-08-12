program libv2

  use mpi_f08
  use w90_comms, only: w90comm_type
  use w90_helper_types

  implicit none

  character(len=100) :: seedname
  complex(kind=dp), allocatable :: a(:,:,:)
  complex(kind=dp), allocatable :: m(:,:,:,:)
  complex(kind=dp), allocatable :: mloc(:,:,:,:)
  complex(kind=dp), allocatable :: morig(:,:,:,:)
  complex(kind=dp), allocatable :: u(:,:,:)
  complex(kind=dp), allocatable :: uopt(:,:,:)
  integer :: mpisize, rank, ierr, stat
  integer :: nb, nk, nw, nn, ik, length, len2, i
  integer, parameter :: stdout = 6
  integer, allocatable :: counts(:), displs(:), distk(:)
  type(lib_global_type) :: w90main
  type(lib_w90_type) :: w90dat
  type(w90comm_type) :: comm 

  call get_command_argument(1, seedname, length, stat)
  if (stat /= 0) then
    write(*,*)'failed to parse seedname'
    stop
  endif
  do i = 1, length
    if (seedname(i:i) == '.') exit
    len2 = i
  enddo
  seedname = seedname(1:len2)

  comm%comm= mpi_comm_world

  call mpi_init(ierr)
  call mpi_comm_rank(comm%comm, rank, ierr)
  call mpi_comm_size(comm%comm, mpisize, ierr)

  if (mpisize > 1) then
    write(*,*)'not yet parallel'
    stop
  endif

  allocate(counts(0:mpisize-1))
  allocate(displs(0:mpisize-1))

  call input_reader(w90main, w90dat, seedname, stdout, stat, comm)

  nw = w90main%num_wann
  nb = w90main%num_bands
  nk = w90main%num_kpts

  ! setup pplel decomp
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

  write(*,*) stat
end program libv2
