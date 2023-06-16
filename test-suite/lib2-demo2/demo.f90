module aux
  integer, parameter :: dp = kind(1.0d0)
  private dp
contains
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

  complex(kind=dp), allocatable :: a(:, :, :)
  complex(kind=dp), allocatable :: m(:, :, :, :)
  complex(kind=dp), allocatable :: morig(:, :, :, :)
  complex(kind=dp), allocatable :: u(:, :, :)
  complex(kind=dp), allocatable :: uopt(:, :, :)
  integer :: ib, ic, nbas, ierr, stdout, stderr
  integer :: i, nkl, mpisize, mpirank
  integer :: nb, nk, nn, nw, nkabc(3), exclude(5)
  integer, allocatable :: distk(:)
  real(kind=dp), allocatable :: eval(:, :)
  real(kind=dp), allocatable :: kpcart(:, :) ! cartesian k list
  real(kind=dp), allocatable :: xcart(:, :) ! cartesian atom positions
  real(kind=dp) :: uccart(3, 3) ! cartesian unit cell

  type(lib_common_type), target :: w90main
  type(lib_wannier_type), target :: w90dat

  ! collect data
  nb = 12
  nw = 8
  nkabc = (/4, 4, 4/)
  exclude = (/1, 2, 3, 4, 5/)
  nk = nkabc(1)*nkabc(2)*nkabc(3)
  call read_cell(nbas, uccart, xcart)
  call read_eval(nb, nk, eval)
  call read_kp(nk, kpcart)

  ! setup MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
  call set_parallel_comms(w90main, mpi_comm_world)

  ! io and k distribution
  call get_fortran_stdout(stdout)
  call get_fortran_stderr(stderr)
  allocate (distk(nk))
  nkl = nk/mpisize ! number of kpoints per rank
  if (mod(nk, mpisize) > 0) nkl = nkl + 1
  do i = 1, nk
    distk(i) = (i - 1)/nkl ! contiguous blocks with potentially fewer processes on last rank
  enddo
  call set_kpoint_distribution(w90main, distk, stdout, stderr, ierr)

  ! required settings
  call set_option(w90main, 'num_wann', nw)
  call set_option(w90main, 'num_bands', nb)
  call set_option(w90main, 'num_kpts', nk)
  call set_option(w90main, 'unit_cell_cart', uccart)
  call set_option(w90main, 'mp_grid', nkabc)
  call set_option(w90main, 'kpoints', kpcart)

  ! optional settings
  call set_option(w90main, 'conv_tol', 1.d-10)
  call set_option(w90main, 'conv_window', 3)
  call set_option(w90main, 'dis_froz_max', 14.0d0) ! not optional for disentanglement? fixme
  call set_option(w90main, 'dis_mix_ratio', 1.d0)
  call set_option(w90main, 'dis_num_iter', 1200)
  call set_option(w90main, 'dis_win_max', 24.d0)
  call set_option(w90main, 'num_iter', 1000)
  call set_option(w90main, 'num_print_cycles', 40)
  call set_option(w90main, 'exclude_bands', exclude)
  call set_option(w90main, 'iprint', 0)
  ! apply settings (requires full set of options)
  call input_setopt(w90main, w90dat, 'gaas', stdout, stderr, ierr)

  ! k setup needs attention
  ! must be done before reading overlaps
  call create_kmesh(w90main, stdout, stderr, ierr)
  nn = w90main%kmesh_info%nntot

  ! setup all matrices
  allocate (a(nb, nw, nk))
  call set_a_matrix(w90dat, a)
  allocate (morig(nb, nb, nn, nk))
  call set_m_orig(w90dat, morig) ! m_matrix_local_orig
  allocate (m(nw, nw, nn, nk)) ! m also known as m_matrix_local
  call set_m_matrix_local(w90dat, m)
  allocate (u(nw, nw, nk))
  call set_u_matrix(w90main, u)
  allocate (uopt(nb, nw, nk))
  call set_u_opt(w90main, uopt)

  ! read from ".mmn" and ".amn"
  ! and assign to m_orig and a
  ! (or m and u if not disentangling)
  call overlaps(w90main, w90dat, stdout, stderr, ierr)

  call set_eigval(w90main, eval)

  call disentangle(w90main, w90dat, stdout, stderr, ierr)
  call wannierise(w90main, w90dat, stdout, stderr, ierr)

  if (mpirank == 0) then ! we print the results
    write (stderr, *) "nb, nw, nk: ", nb, nw, nk
    do ib = 1, nw
      write (stdout, *) (w90main%wannier_data%centres(ic, ib), ic=1, 3), w90main%wannier_data%spreads(ib)
    enddo
  endif

  call mpi_finalize()
end program
