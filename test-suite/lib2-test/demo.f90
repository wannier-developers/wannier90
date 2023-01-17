module aux
  integer, parameter :: dp = kind(1.0d0)
  private dp
contains

  subroutine read_a(nb, nk, nw, a)
    character(len=:), allocatable :: junk
    complex(kind=dp), allocatable, intent(inout) :: a(:, :, :)
    integer, intent(inout) :: nb, nk, nw
    real(kind=dp) :: rea, ima
    open (newunit=iu, file='gaas.amn')
    read (iu, '(a)') junk
    read (iu, *) nb, nk, nw
    allocate (a(nb, nw, nk))
    do ik = 1, nk
    do iw = 1, nw
    do ib = 1, nb
      read (iu, *) ijunk, ijunk, ijunk, rea, ima
      a(ib, iw, ik) = cmplx(rea, ima, kind=dp)
    enddo
    enddo
    enddo
  end subroutine read_a

  subroutine read_eval(nb, nk, eval)
    integer, intent(in) :: nb, nk
    real(kind=dp), allocatable, intent(inout) :: eval(:, :)
    allocate (eval(nb, nk))
    open (newunit=iu, file='gaas.eig')
    do ik = 1, nk
    do ib = 1, nb
      read (iu, *) ijunk, ijunk, eval(ib, ik)
    enddo
    enddo
  end subroutine read_eval

  subroutine read_m(nb, nk, nn, m)
    character(len=:), allocatable :: junk
    complex(kind=dp), allocatable, intent(inout) :: m(:, :, :, :)
    integer, intent(in) :: nb, nk
    integer, intent(out) :: nn
    real(kind=dp) :: rea, ima
    open (newunit=iu, file='gaas.mmn')
    read (iu, '(a)') junk
    read (iu, *) ijunk, ijunk, nn
    allocate (m(nb, nb, nn, nk))
    do ik = 1, nk
    do in = 1, nn
      read (iu, *) ijunk, ijunk, ijunk, ijunk, ijunk
      do jb = 1, nb
      do ib = 1, nb
        read (iu, *) rea, ima
        m(ib, jb, in, ik) = cmplx(rea, ima, kind=dp)
      enddo
      enddo
    enddo
    enddo
  end subroutine read_m

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
  use w90_helper_types
  use w90_comms, only: w90_comm_type

  implicit none

  complex(kind=dp), allocatable :: a(:, :, :)
  complex(kind=dp), allocatable :: m(:, :, :, :)
  complex(kind=dp), allocatable :: morig(:, :, :, :)
  complex(kind=dp), allocatable :: u(:, :, :)
  complex(kind=dp), allocatable :: uopt(:, :, :)
  integer :: nbas, ierr, stdout, stderr
  integer :: nb, nk, nn, nw
  integer, allocatable :: distk(:)
  real(kind=dp), allocatable :: eval(:, :)
  real(kind=dp), allocatable :: kpcart(:, :) ! cartesian k list
  real(kind=dp), allocatable :: xcart(:, :) ! cartesian atom positions
  real(kind=dp) :: uccart(3, 3) ! cartesian unit cell

  type(lib_global_type), target :: w90main
  type(lib_w90_type), target :: w90dat
  type(w90_comm_type) :: comm

  ! MPI
  comm%comm = mpi_comm_world
  call mpi_init(ierr)

  ! accumulate data
  call read_a(nb, nk, nw, a) ! also assigns to nb, nk, nw
  call read_eval(nb, nk, eval)
  call read_kp(nk, kpcart)
  call read_m(nb, nk, nn, morig)
  call read_cell(nbas, uccart, xcart)
  write (*, *) nb, nk, nw, nn, nbas

  ! io and k distribution
  call get_fortran_stdout(stdout)
  call get_fortran_stderr(stderr)
  allocate (distk(nk))
  distk = 0
  call set_kpoint_distribution(w90main, distk)

  ! required settings
  call set_option('num_wann', 8)
  call set_option('num_bands', nb)
  call set_option('num_kpts', nk)
  call set_option('unit_cell_cart', uccart)
  call set_option('mp_grid', (/4, 4, 4/))
  call set_option('kpoints', kpcart)

  call set_option('dis_froz_max', 14.0d0) ! not optional for disentanglement? fixme

  ! optional settings
  call set_option('conv_tol', 1.d-10)
  call set_option('conv_window', 3)
  call set_option('dis_mix_ratio', 1.d0)
  call set_option('num_iter', 1000)
  call set_option('num_print_cycles', 40)

  ! apply and forget settings
  call input_setopt(w90main, w90dat, 'jaja', stdout, stderr, ierr, comm)

  ! k setup needs attention
  call create_kmesh(w90main, stdout, stderr, ierr, comm)

  ! setup all matrices
  allocate (m(nw, nw, nn, nk)) !m=mloc, morig=morig
  allocate (u(nw, nw, nk))
  allocate (uopt(nb, nw, nk))

  call set_m_matrix_local(w90dat, m)
  call set_u_matrix(w90main, u)
  call set_u_opt(w90main, uopt)
  call set_a_matrix(w90dat, a)
  call set_m_orig(w90dat, morig)
  call set_eigval(w90main, eval)

  call disentangle(w90main, w90dat, stdout, stderr, ierr, comm)
  call wannierise(w90main, w90dat, stdout, stderr, ierr, comm)
end program
