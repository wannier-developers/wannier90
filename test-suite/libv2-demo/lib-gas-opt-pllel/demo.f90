program ok
  use mpi_f08
  use w90_library

  implicit none

  complex(8), allocatable :: m_matrix(:, :, :, :)
  complex(8), allocatable :: u_matrix(:, :, :)
  complex(8), allocatable :: u_matrix_opt(:, :, :)
  integer, allocatable :: distk(:), nnkp(:, :)
  integer :: exclude(5)
  integer :: ierr, iu, stdout, stderr
  integer :: i, ib, ic, ik, nkl
  integer :: ika, ikb, ikc, nkabc(3)
  integer :: mpisize, mpirank
  integer :: nb, nk, nn, nw
  real(8), allocatable :: eval(:, :), kpt(:, :)
  real(8) :: uccart(3, 3) ! cartesian unit cell
  type(lib_common_type), target :: w90main

  ! collect data
  nb = 12
  nw = 8
  nkabc = (/4, 4, 4/)
  exclude = (/1, 2, 3, 4, 5/)
  nk = nkabc(1)*nkabc(2)*nkabc(3)
  uccart(1, 1) = -2.8258062938705995d0
  uccart(2, 1) = 0d0
  uccart(3, 1) = 2.8258062938705995d0
  uccart(1, 2) = 0d0
  uccart(2, 2) = 2.8258062938705995d0
  uccart(3, 2) = 2.8258062938705995d0
  uccart(1, 3) = -2.8258062938705995d0
  uccart(2, 3) = 2.8258062938705995d0
  uccart(3, 3) = 0d0

  ! gather eigenvalues
  allocate (eval(nb, nk))
  open (newunit=iu, file='gaas.eig')
  do ik = 1, nk
  do ib = 1, nb
    read (iu, *) i, i, eval(ib, ik)
  enddo
  enddo
  close (iu)

  ! kpoint vectors in w90 order
  i = 0
  allocate (kpt(3, nk))
  do ika = 0, nkabc(1) - 1
    do ikb = 0, nkabc(2) - 1
      do ikc = 0, nkabc(3) - 1
        i = i + 1
        kpt(1, i) = dble(ika)/dble(nkabc(1)); 
        kpt(2, i) = dble(ikb)/dble(nkabc(2)); 
        kpt(3, i) = dble(ikc)/dble(nkabc(3)); 
      enddo
    enddo
  enddo

  ! setup MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  call mpi_comm_rank(mpi_comm_world, mpirank, ierr)

  ! crude k distribution
  allocate (distk(nk))
  nkl = nk/mpisize ! number of kpoints per rank
  if (mod(nk, mpisize) > 0) nkl = nkl + 1
  do i = 1, nk
    distk(i) = (i - 1)/nkl ! contiguous blocks with potentially fewer processes on last rank
  enddo

  ! wannier interface starts
  ! stdout/err
  call get_fortran_stdout(stdout)
  call get_fortran_stderr(stderr)

  ! required settings
  call set_option(w90main, 'kpoints', kpt)
  call set_option(w90main, 'mp_grid', nkabc)
  call set_option(w90main, 'num_bands', nb)
  call set_option(w90main, 'num_kpts', nk)
  call set_option(w90main, 'num_wann', nw)
  call set_option(w90main, 'unit_cell_cart', uccart)

  ! optional settings
  call set_option(w90main, 'conv_tol', 1.d-13)
  call set_option(w90main, 'conv_window', 3)
  call set_option(w90main, 'dis_froz_max', 14.0d0)
  call set_option(w90main, 'dis_mix_ratio', 1.d0)
  call set_option(w90main, 'dis_num_iter', 1200)
  call set_option(w90main, 'distk', distk)
  call set_option(w90main, 'dis_win_max', 24.d0)
  call set_option(w90main, 'exclude_bands', exclude)
  call set_option(w90main, 'fixed_step', 50.d0)
  call set_option(w90main, 'iprint', 0)
  call set_option(w90main, 'num_iter', 1000)

  call input_setopt(w90main, 'gaas', mpi_comm_world, stdout, stderr, ierr) ! apply settings

  call get_nn(w90main, nn, stdout, stderr, ierr); 
  allocate (nnkp(nk, nn))
  call get_nnkp(w90main, nnkp, stdout, stderr, ierr); 
  allocate (m_matrix(nb, nb, nn, nk))
  allocate (u_matrix_opt(nb, nw, nk))
  call set_m_local(w90main, m_matrix) ! m_matrix_local_orig
  call set_u_opt(w90main, u_matrix_opt)

  ! read from ".mmn" and ".amn"
  ! and assign to m and a (now called u)
  ! a dft code would calculate the overlaps here instead
  call overlaps(w90main, stdout, stderr, ierr)

  ! final u matrix
  allocate (u_matrix(nw, nw, nk))
  call set_u_matrix(w90main, u_matrix)

  call set_eigval(w90main, eval)
  call disentangle(w90main, stdout, stderr, ierr)
  call wannierise(w90main, stdout, stderr, ierr)

  if (mpirank == 0) then
    do ib = 1, nw
      write (stdout, '(4f20.10)') (w90main%wannier_data%centres(ic, ib), ic=1, 3), w90main%wannier_data%spreads(ib)
    enddo
  endif
  call mpi_finalize()
end program
