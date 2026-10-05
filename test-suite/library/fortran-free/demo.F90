!-*- mode: F90 -*-!
!------------------------------------------------------------!
! Copyright (C) 2026 Wannier Developer Group                 !
!                                                            !
! This library is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU Lesser General !
! Public License as published by the Free Software           !
! Foundation; either version 2.1 of the License, or (at your !
! option) any later version.                                 !
!                                                            !
! This library is distributed in the hope that it will be    !
! useful,but WITHOUT ANY WARRANTY; without even the implied  !
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR    !
! PURPOSE.  See the GNU Lesser General Public License for    !
! more details.                                              !
!                                                            !
! You should have received a copy of the GNU Lesser General  !
! Public License along with this library; if not, see        !
! <https://www.gnu.org/licenses/>.                           !
!                                                            !
! The webpage of the Wannier90 code is                       !
! <https://www.wannier.org>.                                 !
!                                                            !
! The Wannier90 code is hosted on GitHub                     !
! <https://github.com/wannier-developers/wannier90>          !
!------------------------------------------------------------!

program ok

! fallback to f90 mpi interface for this test's simplified cmake setup
#ifdef W90_MPI08
  use mpi_f08
#else
  use mpi
#endif
  use w90_library
  use w90_library_extra, only: overlaps

  implicit none

  complex(8), allocatable :: m_matrix(:, :, :, :)
  complex(8), allocatable :: u_matrix(:, :, :)
  complex(8), allocatable :: u_matrix_opt(:, :, :)
  integer, allocatable :: distk(:), nnkp(:, :)
  character(len=256) :: exclude
  integer :: ierr, iu, stdout, stderr
  integer :: i, ib, ic, ik, nkl
  integer :: ika, ikb, ikc, nkabc(3)
  integer :: mpisize, mpirank
  integer :: nb, nk, nn, nw, nkloc, ifu, lib_iter, mx_lib_iter, nwarm, size_pages, rss_pages, rss_warm
  integer, parameter :: rss_tol = 4 ! pages; a leaked lib_common_type is several pages per pass
  real(8), allocatable :: centres1(:, :), spreads1(:)
  real(8) :: dev
  logical :: fail
  real(8), allocatable :: eval(:, :), kpt(:, :)
  real(8) :: uccart(3, 3) ! cartesian unit cell
  type(lib_common_type) :: w90main

  ! collect data
  nb = 12
  nw = 8
  nkabc = (/4, 4, 4/)
  exclude = '1-5'
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
  end do
  end do
  close (iu)

  ! kpoint vectors in w90 order
  i = 0
  allocate (kpt(3, nk))
  do ika = 0, nkabc(1) - 1
    do ikb = 0, nkabc(2) - 1
      do ikc = 0, nkabc(3) - 1
        i = i + 1
        kpt(1, i) = dble(ika)/dble(nkabc(1))
        kpt(2, i) = dble(ikb)/dble(nkabc(2))
        kpt(3, i) = dble(ikc)/dble(nkabc(3))
      end do
    end do
  end do

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
  end do
  nkloc = count(distk(:) == mpirank)

  ! stdout/err
  call w90_get_fortran_stdout(stdout)
  call w90_get_fortran_stderr(stderr)

  ! the library instance is set up, run and freed repeatedly; resident memory must not grow
  ! after the allocator warm-up, and every pass must reproduce the first
  mx_lib_iter = 20
  nwarm = 5
  fail = .false.
  allocate (centres1(3, nw), spreads1(nw))
  do lib_iter = 1, mx_lib_iter

    ! required settings
    call w90_set_option(w90main, 'kpoints', kpt)
    call w90_set_option(w90main, 'mp_grid', nkabc)
    call w90_set_option(w90main, 'num_bands', nb)
    !call w90_set_option(w90main, 'num_kpts', nk)
    call w90_set_option(w90main, 'num_wann', nw)
    call w90_set_option(w90main, 'unit_cell_cart', uccart)

    ! optional settings
    call w90_set_option(w90main, 'conv_tol', 1.d-10)
    call w90_set_option(w90main, 'conv_window', 3)
    call w90_set_option(w90main, 'dis_froz_max', 14.0d0)
    call w90_set_option(w90main, 'dis_mix_ratio', 1.d0)
    call w90_set_option(w90main, 'dis_num_iter', 1200)
    call w90_set_option(w90main, 'distk', distk)
    call w90_set_option(w90main, 'dis_win_max', 24.d0)
    call w90_set_option(w90main, 'exclude_bands', exclude)
    call w90_set_option(w90main, 'iprint', 0) ! disable printout
    call w90_set_option(w90main, 'num_iter', 1000)
    call w90_set_option(w90main, 'num_print_cycles', 40)

    call w90_set_comm(w90main, mpi_comm_world)
    call w90_input_setopt(w90main, 'gaas', stdout, stderr, ierr) ! apply settings

    call w90_get_nn(w90main, nn, stdout, stderr, ierr)
    allocate (nnkp(nk, nn))
    call w90_get_nnkp(w90main, nnkp, stdout, stderr, ierr)
    allocate (m_matrix(nb, nb, nn, nkloc))
    allocate (u_matrix_opt(nb, nw, nk))
    call w90_set_m_local(w90main, m_matrix) ! m_matrix_local_orig
    call w90_set_u_opt(w90main, u_matrix_opt)

    ! read from ".mmn" and ".amn"
    ! and assign to m and a (now called u)
    ! a dft code would calculate the overlaps here instead
    call overlaps(w90main, stdout, stderr, ierr) ! from library-extra

    ! pass pointer to eval array
    call w90_set_eigval(w90main, eval)

    ! final u matrix
    allocate (u_matrix(nw, nw, nk))
    call w90_set_u_matrix(w90main, u_matrix)

    call w90_disentangle(w90main, stdout, stderr, ierr)
    call w90_wannierise(w90main, stdout, stderr, ierr)

    ! reuse check: identical input must give identical output on every pass
    if (lib_iter == 1) then
      centres1 = w90main%wannier_data%centres
      spreads1 = w90main%wannier_data%spreads
    else
      dev = max(maxval(abs(w90main%wannier_data%centres - centres1)), &
                maxval(abs(w90main%wannier_data%spreads - spreads1)))
      if (dev > 0d0) then
        fail = .true.
        if (mpirank == 0) write (*, '(a,i0,a,es10.2)') 'pass ', lib_iter, ' differs from pass 1 by ', dev
      end if
    end if

    ! leak check: resident pages (linux only)
    open (newunit=ifu, file='/proc/self/statm', action='read')
    read (ifu, *) size_pages, rss_pages
    close (ifu)
    if (mpirank == 0) write (*, '(a,i4,2(a,i8))') 'pass ', lib_iter, '  vsz pages ', size_pages, '  rss pages ', rss_pages
    if (lib_iter == nwarm) rss_warm = rss_pages

    if (mpirank == 0 .and. lib_iter == mx_lib_iter) then
      open (newunit=ifu, file="results.dat")
      do ib = 1, nw
        write (ifu, '(4f20.10)') (w90main%wannier_data%centres(ic, ib), ic=1, 3), w90main%wannier_data%spreads(ib)
      end do
      close (ifu)
    end if

    deallocate (nnkp)
    deallocate (m_matrix)
    deallocate (u_matrix_opt)
    deallocate (u_matrix)

    call w90_free(w90main)
  end do

  if (rss_pages - rss_warm > rss_tol) then
    fail = .true.
    if (mpirank == 0) write (*, '(a,i0,a,i0,a)') 'rss grew by ', rss_pages - rss_warm, ' pages over ', &
      mx_lib_iter - nwarm, ' passes'
  end if

  if (fail) then
    if (mpirank == 0) write (*, *) 'FAIL'
    call mpi_abort(mpi_comm_world, 1, ierr)
  end if
  if (mpirank == 0) write (*, *) 'free/reuse test OK'

  call mpi_finalize(ierr)
end program
