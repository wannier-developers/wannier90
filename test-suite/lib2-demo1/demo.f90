module aux
  integer, parameter :: dp = kind(1.0d0)
  private dp
contains

  subroutine read_kp(nk, kpcart)
    integer, intent(in) :: nk
    real(kind=dp), allocatable, intent(out) :: kpcart(:, :)
    allocate (kpcart(3, nk))
    open (newunit=iu, file='KPOINTS')
    do ik = 1, nk
      read (iu, *) (kpcart(ic, ik), ic=1, 3)
    enddo
  end subroutine read_kp

  subroutine read_cell(nbas, uccart, xcart)
    character(len=100) :: junk
    integer, intent(inout) :: nbas
    real(kind=dp), allocatable, intent(inout) :: xcart(:, :)
    real(kind=dp), intent(inout) :: uccart(:, :)
    open (newunit=iu, file='CELL')
    do ib = 1, 3
      read (iu, *) (uccart(ic, ib), ic=1, 3)
    enddo
    open (newunit=iu, file='POSITIONS')
    read (iu, *) nbas
    allocate (xcart(3, nbas))
    do ib = 1, nbas
      read (iu, '(a)') junk
    enddo
    do ib = 1, nbas
      read (iu, *) (xcart(ic, ib), ic=1, 3)
    enddo
  end subroutine read_cell
end module

program ok
  use aux
  use mpi_f08
  use w90_library

  implicit none

  complex(kind=dp), allocatable :: m(:, :, :, :)
  complex(kind=dp), allocatable :: u(:, :, :)
  complex(kind=dp), allocatable :: uopt(:, :, :)
  integer :: ib, ic, nbas, ierr, stdout, stderr
  integer :: nb, nk, nn, nw, nkabc(3)
  integer, allocatable :: distk(:)
  real(kind=dp), allocatable :: kpcart(:, :) ! cartesian k list
  real(kind=dp), allocatable :: xcart(:, :) ! cartesian atom positions
  real(kind=dp) :: uccart(3, 3) ! cartesian unit cell

  type(lib_common_type), target :: w90main
  type(lib_wannier_type), target :: w90dat

  nb = 4
  nw = 4
  nkabc = (/2, 2, 2/)
  nk = nkabc(1)*nkabc(2)*nkabc(3)
  call read_kp(nk, kpcart)
  call read_cell(nbas, uccart, xcart)

  ! MPI
  call mpi_init(ierr)

  ! io and k distribution
  call get_fortran_stdout(stdout)
  !open (newunit=stdout, file='/dev/null')
  call get_fortran_stderr(stderr)
  allocate (distk(nk))
  distk = 0
  call set_kpoint_distribution(w90main, distk, stdout, stderr, ierr)

  ! required settings
  call set_option(w90main, 'num_wann', nw)
  call set_option(w90main, 'num_bands', nb)
  call set_option(w90main, 'num_kpts', nk)  ! why are we providing nk and nkabc and kpcart (aka kpt_latt)?
  call set_option(w90main, 'mp_grid', nkabc)
  call set_option(w90main, 'kpoints', kpcart)
  call set_option(w90main, 'unit_cell_cart', uccart)
  call input_setopt_special(w90main, w90dat, 'wannier', stdout, stderr, ierr)

  ! read some other settings from win file
  call input_reader(w90main, w90dat, 'wannier', stdout, stderr, ierr)

  ! k setup needs attention
  ! must be done before reading overlaps
  call create_kmesh(w90main, stdout, stderr, ierr)
  nn = w90main%kmesh_info%nntot

  ! setup all matrices
  allocate (m(nw, nw, nn, nk)) ! m also known as m_matrix_local
  call set_m_matrix_local(w90dat, m)
  allocate (u(nw, nw, nk))
  call set_u_matrix(w90main, u)

  ! read from ".mmn" and ".amn"
  ! and assign to m and u (or m_orig and a if disentangling)
  call overlaps(w90main, w90dat, stdout, stderr, ierr)
  call projovlp(w90main, w90dat, stdout, stderr, ierr)
  allocate (uopt(nb, nw, nk))
  call set_u_opt(w90main, uopt)

  call wannierise(w90main, w90dat, stdout, stderr, ierr)

  ! printout result
  nbas = size(w90main%wannier_data%centres, 2)
  write (stdout, *) "wannier centres and spreads"
  do ib = 1, nbas
    write (stdout, *) (w90main%wannier_data%centres(ic, ib), ic=1, 3), w90main%wannier_data%spreads(ib)
  enddo
  write (stdout, *) "sum"
  write (stdout, *) (sum(w90main%wannier_data%centres(ic, :)), ic=1, 3), sum(w90main%wannier_data%spreads(:))

  call mpi_finalize()
end program
