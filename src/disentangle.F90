!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

module w90_disentangle
  !! This module contains the core routines to extract an optimal
  !! subspace from a set of entangled bands.

  use w90_constants, only: dp, cmplx_0, cmplx_1
  use w90_io, only: io_error, stdout, io_stopwatch
  use w90_parameters, only: num_bands, num_wann, a_matrix, u_matrix_opt, &
    u_matrix, m_matrix_orig, lwindow, dis_conv_window, devel_flag, &
    nntot, timing_level, omega_invariant, u_matrix, lsitesymmetry, &
    lenconfac, iprint, wbtot, dis_num_iter, dis_mix_ratio, dis_win_min, &
    dis_win_max, dis_froz_min, dis_froz_max, dis_spheres_num, &
    dis_spheres_first_wann, num_kpts, nnlist, ndimwin, wb, gamma_only, &
    eigval, length_unit, dis_spheres, m_matrix, dis_conv_tol, frozen_states, &
    optimisation, recip_lattice, kpt_latt, &
    m_matrix_orig_local, m_matrix_local

  use w90_comms, only: on_root, my_node_id, num_nodes, &
    comms_bcast, comms_array_split, &
    comms_gatherv, comms_allreduce

  use w90_sitesym, only: sitesym_slim_d_matrix_band, &
    sitesym_replace_d_matrix_band, sitesym_symmetrize_u_matrix, &
    sitesym_symmetrize_zmatrix, sitesym_dis_extract_symmetry !RS:

  implicit none

  private

  real(kind=dp), allocatable :: eigval_opt(:, :)
  !! At input it contains a large set of eigenvalues. At
  !! it is slimmed down to contain only those inside the energy window.

  logical                :: linner
  !! Is there a frozen window
  logical, allocatable   :: lfrozen(:, :)
  !! true if the i-th band inside outer window is frozen
  integer, allocatable   :: nfirstwin(:)
  !! index of lowest band inside outer window at nkp-th
  integer, allocatable   :: ndimfroz(:)
  !! number of frozen bands at nkp-th k point
  integer, allocatable   :: indxfroz(:, :)
  !! number of bands inside outer window at nkp-th k point
  integer, allocatable   :: indxnfroz(:, :)
  !!   outer-window band index for the i-th non-frozen state
  !! (equals 1 if it is the bottom of outer window)

  public :: dis_main

contains

  !==================================================================!
  subroutine dis_main()
    !==================================================================!
    !! Main disentanglement routine
    !                                                                  !
    !                                                                  !
    !                                                                  !
    !==================================================================!
    use w90_io, only: io_file_unit

    ! internal variables
    integer                       :: nkp, nkp2, nn, j, ierr, page_unit
    integer                       :: nkp_global
    complex(kind=dp), allocatable :: cwb(:, :), cww(:, :)
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs

    if (timing_level > 0) call io_stopwatch('dis: main', 1)

    call comms_array_split(num_kpts, counts, displs)

    if (on_root) write (stdout, '(/1x,a)') &
      '*------------------------------- DISENTANGLE --------------------------------*'

    ! Allocate arrays
    allocate (eigval_opt(num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating eigval_opt in dis_main')
    eigval_opt = eigval

    ! Set up energy windows
    call dis_windows()

    ! Construct the unitarized projection
    call dis_project()

    ! If there is an inner window, need to modify projection procedure
    ! (Sec. III.G SMV)
    if (linner) then
      if (lsitesymmetry) call io_error('in symmetry-adapted mode, frozen window not implemented yet') !YN: RS:
      if (on_root) write (stdout, '(3x,a)') 'Using an inner window (linner = T)'
      call dis_proj_froz()
    else
      if (on_root) write (stdout, '(3x,a)') 'No inner window (linner = F)'
    endif

    ! Debug
    call internal_check_orthonorm()

    ! Slim down the original Mmn(k,b)
    call internal_slim_m()

    lwindow = .false.
    do nkp = 1, num_kpts
      do j = nfirstwin(nkp), nfirstwin(nkp) + ndimwin(nkp) - 1
        lwindow(j, nkp) = .true.
      end do
    end do

    if (lsitesymmetry) call sitesym_slim_d_matrix_band(lwindow)                         !RS: calculate initial U_{opt}(Rk) from U_{opt}(k)
    if (lsitesymmetry) call sitesym_symmetrize_u_matrix(num_bands, u_matrix_opt, lwindow) !RS:
    ! Extract the optimally-connected num_wann-dimensional subspaces
![ysl-b]
    if (.not. gamma_only) then
      call dis_extract()
    else
      call dis_extract_gamma()
    end if
![ysl-e]

    ! Allocate workspace
    allocate (cwb(num_wann, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwb in dis_main')
    allocate (cww(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cww in dis_main')

    ! Find the num_wann x num_wann overlap matrices between
    ! the basis states of the optimal subspaces
    do nkp = 1, counts(my_node_id)
      nkp_global = nkp + displs(my_node_id)
      do nn = 1, nntot
        nkp2 = nnlist(nkp_global, nn)
        call zgemm('C', 'N', num_wann, ndimwin(nkp2), ndimwin(nkp_global), cmplx_1, &
                   u_matrix_opt(:, :, nkp_global), num_bands, m_matrix_orig_local(:, :, nn, nkp), num_bands, &
                   cmplx_0, cwb, num_wann)
        call zgemm('N', 'N', num_wann, num_wann, ndimwin(nkp2), cmplx_1, &
                   cwb, num_wann, u_matrix_opt(:, :, nkp2), num_bands, &
                   cmplx_0, cww, num_wann)
        m_matrix_orig_local(1:num_wann, 1:num_wann, nn, nkp) = cww(:, :)
      enddo
    enddo

    ! Find the initial u_matrix
    if (lsitesymmetry) call sitesym_replace_d_matrix_band() !RS: replace d_matrix_band here
![ysl-b]
    if (.not. gamma_only) then
      call internal_find_u()
    else
      call internal_find_u_gamma()
    end if
![ysl-e]

    if (optimisation <= 0) then
      page_unit = io_file_unit()
      open (unit=page_unit, form='unformatted', status='scratch')
      ! Update the m_matrix accordingly
      do nkp = 1, counts(my_node_id)
        nkp_global = nkp + displs(my_node_id)
        do nn = 1, nntot
          nkp2 = nnlist(nkp_global, nn)
          call zgemm('C', 'N', num_wann, num_wann, num_wann, cmplx_1, &
                     u_matrix(:, :, nkp_global), num_wann, m_matrix_orig_local(:, :, nn, nkp), num_bands, &
                     cmplx_0, cwb, num_wann)
          call zgemm('N', 'N', num_wann, num_wann, num_wann, cmplx_1, &
                     cwb, num_wann, u_matrix(:, :, nkp2), num_wann, &
                     cmplx_0, cww, num_wann)
          write (page_unit) cww(:, :)
        enddo
      enddo
      rewind (page_unit)
      deallocate (m_matrix_orig_local, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_orig_local in dis_main')
      if (on_root) then
        allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating m_matrix in dis_main')
      endif
      allocate (m_matrix_local(num_wann, num_wann, nntot, counts(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating m_matrix_local in dis_main')
      do nkp = 1, counts(my_node_id)
        do nn = 1, nntot
          read (page_unit) m_matrix_local(:, :, nn, nkp)
        end do
      end do
      call comms_gatherv(m_matrix_local, num_wann*num_wann*nntot*counts(my_node_id), &
                         m_matrix, num_wann*num_wann*nntot*counts, num_wann*num_wann*nntot*displs)
      close (page_unit)

    else

      if (on_root) then
        allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating m_matrix in dis_main')
      endif
      allocate (m_matrix_local(num_wann, num_wann, nntot, counts(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating m_matrix_local in dis_main')
      ! Update the m_matrix accordingly
      do nkp = 1, counts(my_node_id)
        nkp_global = nkp + displs(my_node_id)
        do nn = 1, nntot
          nkp2 = nnlist(nkp_global, nn)
          call zgemm('C', 'N', num_wann, num_wann, num_wann, cmplx_1, &
                     u_matrix(:, :, nkp_global), num_wann, m_matrix_orig_local(:, :, nn, nkp), num_bands, &
                     cmplx_0, cwb, num_wann)
          call zgemm('N', 'N', num_wann, num_wann, num_wann, cmplx_1, &
                     cwb, num_wann, u_matrix(:, :, nkp2), num_wann, &
                     cmplx_0, cww, num_wann)
          m_matrix_local(:, :, nn, nkp) = cww(:, :)
        enddo
      enddo
      call comms_gatherv(m_matrix_local, num_wann*num_wann*nntot*counts(my_node_id), &
                         m_matrix, num_wann*num_wann*nntot*counts, num_wann*num_wann*nntot*displs)
      deallocate (m_matrix_orig_local, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_orig_local in dis_main')

    endif

    deallocate (a_matrix, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating a_matrix in dis_main')

    ! Deallocate workspace
    deallocate (cww, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cww in dis_main')
    deallocate (cwb, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cwb in dis_main')

    !zero the unused elements of u_matrix_opt (just in case...)
    do nkp = 1, num_kpts
      do j = 1, num_wann
        if (ndimwin(nkp) < num_bands) &
          u_matrix_opt(ndimwin(nkp) + 1:, j, nkp) = cmplx_0
      enddo
    enddo

!~![ysl-b]
!~!   Apply phase factor ph_g if gamma_only
!~    if (.not. gamma_only) then
!~       do nkp = 1, num_kpts
!~          do j = 1, num_wann
!~             u_matrix_opt(1:ndimwin(nkp),j,nkp)  = u_matrix_opt(1:ndimwin(nkp),j,nkp)
!~          enddo
!~       enddo
!~    else
!~       do nkp = 1, num_kpts
!~          do j = 1, ndimwin(nkp)
!~             u_matrix_opt(j,1:num_wann,nkp)  = conjg(ph_g(j))*u_matrix_opt(j,1:num_wann,nkp)
!~          enddo
!~       enddo
!~    endif
!~![ysl-e]

    ! Deallocate module arrays
    call internal_dealloc()

    if (timing_level > 0 .and. on_root) call io_stopwatch('dis: main', 2)

    return

  contains

    !================================================================!
    subroutine internal_check_orthonorm()
      !================================================================!
      !                                                                !
      !! This subroutine checks that the states in the columns of the
      !! final matrix U_opt are orthonormal at every k-point, i.e.,
      !! that the matrix is unitary in the sense that
      !! conjg(U_opt).U_opt = 1  (but not  U_opt.conjg(U_opt) = 1).
      !!
      !! In particular, this checks whether the projected gaussians
      !! are indeed orthogonal to the frozen states, at those k-points
      !! where both are present in the trial subspace.
      !                                                                !
      !================================================================!

      use w90_constants, only: eps8

      implicit none

      integer          :: nkp, l, m, j
      complex(kind=dp) :: ctmp

      if (timing_level > 1) call io_stopwatch('dis: main: check_orthonorm', 1)

      do nkp = 1, num_kpts
        do l = 1, num_wann
          do m = 1, l
            ctmp = cmplx_0
            do j = 1, ndimwin(nkp)
              ctmp = ctmp + conjg(u_matrix_opt(j, m, nkp))*u_matrix_opt(j, l, nkp)
            enddo
            if (l .eq. m) then
              if (abs(ctmp - cmplx_1) .gt. eps8) then
                if (on_root) write (stdout, '(3i6,2f16.12)') nkp, l, m, ctmp
                if (on_root) write (stdout, '(1x,a)') 'The trial orbitals for disentanglement are not orthonormal'
!                     write(stdout,'(1x,a)') 'Try re-running the calculation with the input keyword'
!                     write(stdout,'(1x,a)') '  devel_flag=orth-fix'
!                     write(stdout,'(1x,a)') 'Please report the sucess or failure of this to the Wannier90 developers'
                call io_error('Error in dis_main: orthonormal error 1')
              endif
            else
              if (abs(ctmp) .gt. eps8) then
                if (on_root) write (stdout, '(3i6,2f16.12)') nkp, l, m, ctmp
                if (on_root) write (stdout, '(1x,a)') 'The trial orbitals for disentanglement are not orthonormal'
!                     write(stdout,'(1x,a)') 'Try re-running the calculation with the input keyword'
!                     write(stdout,'(1x,a)') '  devel_flag=orth-fix'
!                     write(stdout,'(1x,a)') 'Please report the sucess or failure of this to the Wannier90 developers'
                call io_error('Error in dis_main: orthonormal error 2')
              endif
            endif
          enddo
        enddo
      enddo

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: main: check_orthonorm', 2)

      return

    end subroutine internal_check_orthonorm

    !================================================================!
    subroutine internal_slim_m()
      !================================================================!
      !                                                                !
      !! This subroutine slims down the original Mmn(k,b), removing
      !! rows and columns corresponding to u_nks that fall outside
      !! the outer energy window.
      !                                                                !
      !================================================================!

      implicit none

      integer                       :: nkp, nkp2, nn, i, j, m, n, ierr
      integer                       :: nkp_global
      complex(kind=dp), allocatable :: cmtmp(:, :)
      ! Needed to split an array on different nodes
      integer, dimension(0:num_nodes - 1) :: counts
      integer, dimension(0:num_nodes - 1) :: displs

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: main: slim_m', 1)

      call comms_array_split(num_kpts, counts, displs)

      allocate (cmtmp(num_bands, num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating cmtmp in dis_main')

      do nkp = 1, counts(my_node_id)
        nkp_global = nkp + displs(my_node_id)
        do nn = 1, nntot
          nkp2 = nnlist(nkp_global, nn)
          do j = 1, ndimwin(nkp2)
            n = nfirstwin(nkp2) + j - 1
            do i = 1, ndimwin(nkp_global)
              m = nfirstwin(nkp_global) + i - 1
              cmtmp(i, j) = m_matrix_orig_local(m, n, nn, nkp)
            enddo
          enddo
          m_matrix_orig_local(:, :, nn, nkp) = cmplx_0
          do j = 1, ndimwin(nkp2)
            do i = 1, ndimwin(nkp_global)
              m_matrix_orig_local(i, j, nn, nkp) = cmtmp(i, j)
            enddo
          enddo
        enddo
      enddo

      deallocate (cmtmp, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating cmtmp in dis_main')

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: main: slim_m', 2)

      return

    end subroutine internal_slim_m

    !================================================================!
    subroutine internal_find_u()
      !================================================================!
      !                                                                !
      !! This subroutine finds the initial guess for the square unitary
      !! rotation matrix u_matrix. The method is similar to Sec. III.D
      !! of SMV, but with square instead of rectangular matrices:
      !!
      !! First find caa, the square overlap matrix <psitilde_nk|g_m>,
      !! where psitilde is an eigenstate of the optimal subspace.
      !!
      !! Note that, contrary to what is implied in Sec. III.E of SMV,
      !! this does *not* need to be computed by brute: instead we take
      !! advantage of the previous computation of overlaps with the
      !! same projections that are used to initiate the minimization of
      !! Omega.
      !!
      !! Note: |psi> U_opt = |psitilde> and obviously
      !! <psitilde| = (U_opt)^dagger <psi|
      !                                                                !
      !================================================================!

      use w90_sitesym, only: ir2ik, ik2ir !YN: RS:
      implicit none

      integer                       :: nkp, info, ierr
      complex(kind=dp), allocatable :: caa(:, :, :)
      ! For ZGESVD
      real(kind=dp), allocatable :: svals(:)
      real(kind=dp), allocatable :: rwork(:)
      complex(kind=dp), allocatable :: cv(:, :)
      complex(kind=dp), allocatable :: cz(:, :)
      complex(kind=dp), allocatable :: cwork(:)

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: main: find_u', 1)

      ! Currently, this part is not parallelized; thus, we perform the task only on root and then broadcast the result.
      if (on_root) then
        ! Allocate arrays needed for ZGESVD
        allocate (svals(num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating svals in dis_main')
        allocate (rwork(5*num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating rwork in dis_main')
        allocate (cv(num_wann, num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating cv in dis_main')
        allocate (cz(num_wann, num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating cz in dis_main')
        allocate (cwork(4*num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating cwork in dis_main')
        allocate (caa(num_wann, num_wann, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating caa in dis_main')

        do nkp = 1, num_kpts
          if (lsitesymmetry) then                 !YN: RS:
            if (ir2ik(ik2ir(nkp)) .ne. nkp) cycle  !YN: RS:
          endif                                   !YN: RS:
          call zgemm('C', 'N', num_wann, num_wann, ndimwin(nkp), cmplx_1, &
                     u_matrix_opt(:, :, nkp), num_bands, a_matrix(:, :, nkp), num_bands, &
                     cmplx_0, caa(:, :, nkp), num_wann)
          ! Singular-value decomposition
          call ZGESVD('A', 'A', num_wann, num_wann, caa(:, :, nkp), num_wann, &
                      svals, cz, num_wann, cv, num_wann, cwork, 4*num_wann, rwork, info)
          if (info .ne. 0) then
            if (on_root) write (stdout, *) ' ERROR: IN ZGESVD IN dis_main'
            if (on_root) write (stdout, *) 'K-POINT NKP=', nkp, ' INFO=', info
            if (info .lt. 0) then
              if (on_root) write (stdout, *) 'THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
            endif
            call io_error('dis_main: problem in ZGESVD 1')
          endif
          ! u_matrix is the initial guess for the unitary rotation of the
          ! basis states given by the subroutine extract
          call zgemm('N', 'N', num_wann, num_wann, num_wann, cmplx_1, &
                     cz, num_wann, cv, num_wann, cmplx_0, u_matrix(:, :, nkp), num_wann)
        enddo
      endif
      call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts)
!      if (lsitesymmetry) call sitesym_symmetrize_u_matrix(num_wann,u_matrix) !RS:

      if (on_root) then
        ! Deallocate arrays for ZGESVD
        deallocate (caa, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating caa in dis_main')
        deallocate (cwork, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating cwork in dis_main')
        deallocate (cz, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating cz in dis_main')
        deallocate (cv, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating cv in dis_main')
        deallocate (rwork, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating rwork in dis_main')
        deallocate (svals, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating svals in dis_main')
      endif

      if (lsitesymmetry) call sitesym_symmetrize_u_matrix(num_wann, u_matrix) !RS:

      if (timing_level > 1) call io_stopwatch('dis: main: find_u', 2)

      return

    end subroutine internal_find_u

![ysl-b]
    !================================================================!
    subroutine internal_find_u_gamma()
      !================================================================!
      !                                                                !
      !! Make initial u_matrix real
      !! Must be the case when gamma_only = .true.
      !                                                                !
      !================================================================!

      implicit none

      integer                       :: info, ierr
      real(kind=dp), allocatable :: u_opt_r(:, :)
      real(kind=dp), allocatable :: a_matrix_r(:, :)
      real(kind=dp), allocatable :: raa(:, :)
      ! For DGESVD
      real(kind=dp), allocatable :: svals(:)
      real(kind=dp), allocatable :: work(:)
      real(kind=dp), allocatable :: rv(:, :)
      real(kind=dp), allocatable :: rz(:, :)

      if (timing_level > 1) call io_stopwatch('dis: main: find_u_gamma', 1)

      ! Allocate arrays needed for getting a_matrix_r
      allocate (u_opt_r(ndimwin(1), num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating u_opt_r in dis_main')
      allocate (a_matrix_r(ndimwin(1), num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating a_matrix_r in dis_main')

      ! Allocate arrays needed for DGESVD
      allocate (svals(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating svals in dis_main')
      allocate (work(5*num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating rwork in dis_main')
      allocate (rv(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating cv in dis_main')
      allocate (rz(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating cz in dis_main')
      allocate (raa(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating raa in dis_main')

      u_opt_r(:, :) = real(u_matrix_opt(1:ndimwin(1), 1:num_wann, 1), dp)

      a_matrix_r(:, :) = real(a_matrix(1:ndimwin(1), 1:num_wann, 1), kind=dp)

      call dgemm('T', 'N', num_wann, num_wann, ndimwin(1), 1.0_dp, &
                 u_opt_r, ndimwin(1), a_matrix_r, ndimwin(1), &
                 0.0_dp, raa, num_wann)
      ! Singular-value decomposition
      call DGESVD('A', 'A', num_wann, num_wann, raa, num_wann, &
                  svals, rz, num_wann, rv, num_wann, work, 5*num_wann, info)
      if (info .ne. 0) then
        write (stdout, *) ' ERROR: IN DGESVD IN dis_main'
        write (stdout, *) 'K-POINT = Gamma', ' INFO=', info
        if (info .lt. 0) then
          write (stdout, *) 'THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
        endif
        call io_error('dis_main: problem in DGESVD 1')
      endif
      ! u_matrix is the initial guess for the unitary rotation of the
      ! basis states given by the subroutine extract
      call dgemm('N', 'N', num_wann, num_wann, num_wann, 1.0_dp, &
                 rz, num_wann, rv, num_wann, 0.0_dp, raa, num_wann)

      u_matrix(:, :, 1) = cmplx(raa(:, :), 0.0_dp, dp)

      ! Deallocate arrays for DGESVD
      deallocate (raa, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating raa in dis_main')
      deallocate (rz, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating rz in dis_main')
      deallocate (rv, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating rv in dis_main')
      deallocate (work, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating work in dis_main')
      deallocate (svals, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating svals in dis_main')

      ! Deallocate arrays for a_matrix_r
      deallocate (a_matrix_r, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating a_matrix_r in dis_main')
      deallocate (u_opt_r, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating u_opt_r in dis_main')

      if (timing_level > 1) call io_stopwatch('dis: main: find_u_gamma', 2)

      return

    end subroutine internal_find_u_gamma
![ysl-e]

    !==================================!
    subroutine internal_dealloc()
      !==================================!
      !! Deallocate module data
      !                                  !
      !==================================!

      implicit none

      integer :: ierr

      ! Module arrays allocated in dis_windows
      deallocate (lfrozen, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating lfrozen in dis_main')
      deallocate (indxnfroz, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating indxnfroz in dis_main')
      deallocate (indxfroz, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating indxfroz in dis_main')
      deallocate (ndimfroz, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating ndimfroz in dis_main')
      deallocate (nfirstwin, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating nfirstwin in dis_main')

      ! Module arrays allocated in dis_main
      deallocate (eigval_opt, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating eigval_opt in dis_main')

      return

    end subroutine internal_dealloc

  end subroutine dis_main

  !==================================================================!
  subroutine dis_windows()
    !==================================================================!
    !                                                                  !
    !! This subroutine selects the states that are inside the outer
    !! window (ie, the energy window out of which we fish out the
    !! optimally-connected subspace) and those that are inside the
    !! inner window (that make up the frozen manifold, and are
    !! straightfowardly included as they are). This, in practice,
    !! amounts to slimming down the original num_wann x num_wann overlap
    !! matrices, removing rows and columns that belong to u_nks that
    !! have been excluded forever, and marking (indexing) the rows and
    !! columns that correspond to frozen states.
    !!
    !! Note - in windows eigval_opt are shifted, so the lowest ones go
    !! from nfirstwin(nkp) to nfirstwin(nkp)+ndimwin(nkp)-1, and above
    !! they are set to zero.
    !                                                                  !
    !==================================================================!

    implicit none

    ! internal variables
    integer :: i, j, nkp, ierr
    integer :: imin, imax, kifroz_min, kifroz_max
    !~~ GS-start
    real(kind=dp) :: dk(3), kdr2
    logical :: dis_ok
    !~~ GS-end

    ! OUTPUT:
    !     ndimwin(nkp)   number of bands inside outer window at nkp-th k poi
    !     ndimfroz(nkp)  number of frozen bands at nkp-th k point
    !     lfrozen(i,nkp) true if the i-th band inside outer window is frozen
    !     linner         true if there is an inner window
    !     indxfroz(i,nkp) outer-window band index for the i-th frozen state
    !                     (equals 1 if it is the bottom of outer window)
    !     indxnfroz(i,nkp) outer-window band index for the i-th non-frozen s
    !                     (equals 1 if it is the bottom of outer window)
    !     nfirstwin(nkp) index of lowest band inside outer window at nkp-th
    ! MODIFIED:
    !     eigval_opt(nb,nkp) At input it contains a large set of eigenvalues. At
    !                    it is slimmed down to contain only those inside the
    !                    energy window, stored in nb=1,...,ndimwin(nkp)

    if (timing_level > 1 .and. on_root) call io_stopwatch('dis: windows', 1)

    ! Allocate module arrays
    allocate (nfirstwin(num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating nfirstwin in dis_windows')
    allocate (ndimfroz(num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ndimfroz in dis_windows')
    allocate (indxfroz(num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating indxfroz in dis_windows')
    allocate (indxnfroz(num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating indxnfroz in dis_windows')
    allocate (lfrozen(num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating lfrozen in dis_windows')

    linner = .false.

    if (on_root) write (stdout, '(1x,a)') &
      '+----------------------------------------------------------------------------+'
    if (on_root) write (stdout, '(1x,a)') &
      '|                              Energy  Windows                               |'
    if (on_root) write (stdout, '(1x,a)') &
      '|                              ---------------                               |'
    if (on_root) write (stdout, '(1x,a,f10.5,a,f10.5,a)') &
      '|                   Outer: ', dis_win_min, '  to ', dis_win_max, &
      '  (eV)                   |'
    if (frozen_states) then
      if (on_root) write (stdout, '(1x,a,f10.5,a,f10.5,a)') &
        '|                   Inner: ', dis_froz_min, '  to ', dis_froz_max, &
        '  (eV)                   |'
    else
      if (on_root) write (stdout, '(1x,a)') &
        '|                   No frozen states were specified                          |'
    endif
    if (on_root) write (stdout, '(1x,a)') &
      '+----------------------------------------------------------------------------+'

    do nkp = 1, num_kpts
      ! Check which eigenvalues fall within the outer window
      if ((eigval_opt(1, nkp) .gt. dis_win_max) .or. &
          (eigval_opt(num_bands, nkp) .lt. dis_win_min)) then
        if (on_root) write (stdout, *) ' ERROR AT K-POINT: ', nkp
        if (on_root) write (stdout, *) ' ENERGY WINDOW (eV):    [', dis_win_min, ',', dis_win_max, ']'
        if (on_root) write (stdout, *) ' EIGENVALUE RANGE (eV): [', &
          eigval_opt(1, nkp), ',', eigval_opt(num_bands, nkp), ']'
        call io_error('dis_windows: The outer energy window contains no eigenvalues')
      endif

      ! Note: we assume that eigvals are ordered from the bottom up
      imin = 0
      do i = 1, num_bands
        if (imin .eq. 0) then
          if ((eigval_opt(i, nkp) .ge. dis_win_min) .and. &
              (eigval_opt(i, nkp) .le. dis_win_max)) imin = i
          imax = i
        endif
        if (eigval_opt(i, nkp) .le. dis_win_max) imax = i
      enddo

      ndimwin(nkp) = imax - imin + 1

      nfirstwin(nkp) = imin

      !~~ GS-start
      ! disentangle at the current k-point only if it is within one of the
      ! spheres centered at the k-points listed in kpt_dis
      if (dis_spheres_num .gt. 0) then
        dis_ok = .false.
        ! loop on the sphere centers
        do i = 1, dis_spheres_num
          dk = kpt_latt(:, nkp) - dis_spheres(1:3, i)
          dk = matmul(anint(dk) - dk, recip_lattice(:, :))
          ! if the current k-point is included in at least one sphere,
          ! then perform disentanglement as usual
          if (abs(dot_product(dk, dk)) .lt. dis_spheres(4, i)**2) then
            dis_ok = .true.
            exit
          endif
        enddo
        ! this kpoint is not included in any sphere: no disentaglement
        if (.not. dis_ok) then
          ndimwin(nkp) = num_wann
          nfirstwin(nkp) = dis_spheres_first_wann
        endif
      endif
      !~~ GS-end

      if (ndimwin(nkp) .lt. num_wann) then
        if (on_root) write (stdout, 483) 'Error at k-point ', nkp, &
          ' ndimwin=', ndimwin(nkp), ' num_wann=', num_wann
483     format(1x, a17, i4, a8, i3, a9, i3)
        call io_error('dis_windows: Energy window contains fewer states than number of target WFs')
      endif

      do i = 1, ndimwin(nkp)
        lfrozen(i, nkp) = .false.
      enddo

      ! Check which eigenvalues fall within the inner window
      kifroz_min = 0
      kifroz_max = -1

      ! (note that the above obeys kifroz_max-kifroz_min+1=kdimfroz=0, as we w
      if (frozen_states) then
        do i = imin, imax
          if (kifroz_min .eq. 0) then
            if ((eigval_opt(i, nkp) .ge. dis_froz_min) .and. &
                (eigval_opt(i, nkp) .le. dis_froz_max)) then
              ! Relative to bottom of outer window
              kifroz_min = i - imin + 1
              kifroz_max = i - imin + 1
            endif
          elseif (eigval_opt(i, nkp) .le. dis_froz_max) then
            kifroz_max = kifroz_max + 1
            ! DEBUG
            ! if(kifroz_max.ne.i-imin+1) stop 'something wrong...'
            ! ENDDEBUG
          endif
        enddo
      endif

      ndimfroz(nkp) = kifroz_max - kifroz_min + 1

      if (ndimfroz(nkp) .gt. num_wann) then
        if (on_root) write (stdout, 401) nkp, ndimfroz(nkp), num_wann
401     format(' ERROR AT K-POINT ', i4, ' THERE ARE ', i2, &
               ' BANDS INSIDE THE INNER WINDOW AND ONLY', i2, &
               ' TARGET BANDS')
        if (on_root) write (stdout, 402) (eigval_opt(i, nkp), i=imin, imax)
402     format('BANDS: (eV)', 10(F10.5, 1X))
        call io_error('dis_windows: More states in the frozen window than target WFs')
      endif

      if (ndimfroz(nkp) .gt. 0) linner = .true.
      ! DEBUG
      !         write(*,'(a,i4,a,i2,a,i2)') 'k point ',nkp,
      !     &    ' lowest band in outer win is # ',imin,
      !     &    '   # frozen states is ',ndimfroz(nkp)
      ! ENDDEBUG
      ! Generate index array for frozen states (those inside inner window)
      if (ndimfroz(nkp) .gt. 0) then
        do i = 1, ndimfroz(nkp)
          indxfroz(i, nkp) = kifroz_min + i - 1
          lfrozen(indxfroz(i, nkp), nkp) = .true.
        enddo
        if (indxfroz(ndimfroz(nkp), nkp) .ne. kifroz_max) then
          if (on_root) write (stdout, *) ' Error at k-point ', nkp, ' frozen band #', i
          if (on_root) write (stdout, *) ' ndimfroz=', ndimfroz(nkp)
          if (on_root) write (stdout, *) ' kifroz_min=', kifroz_min
          if (on_root) write (stdout, *) ' kifroz_max=', kifroz_max
          if (on_root) write (stdout, *) ' indxfroz(i,nkp)=', indxfroz(i, nkp)
          call io_error('dis_windows: Something fishy...')
        endif
      endif

      ! Generate index array for non-frozen states
      i = 0
      do j = 1, ndimwin(nkp)
        !           if (lfrozen(j,nkp).eqv..false.) then
        if (.not. lfrozen(j, nkp)) then
          i = i + 1
          indxnfroz(i, nkp) = j
        endif
      enddo

      if (i .ne. ndimwin(nkp) - ndimfroz(nkp)) then
        if (on_root) write (stdout, *) ' Error at k-point: ', nkp
        if (on_root) write (stdout, '(3(a,i5))') ' i: ', i, ' ndimwin: ', ndimwin(nkp), &
          ' ndimfroz: ', ndimfroz(nkp)
        call io_error('dis_windows: i .ne. (ndimwin-ndimfroz) at k-point')
      endif

      ! Slim down eigval vector at present k
      do i = 1, ndimwin(nkp)
        j = nfirstwin(nkp) + i - 1
        eigval_opt(i, nkp) = eigval_opt(j, nkp)
      enddo

      do i = ndimwin(nkp) + 1, num_bands
        eigval_opt(i, nkp) = 0.0_dp
      enddo

    enddo
    ! [k-point loop (nkp)]

![ysl-b]
!~    if (gamma_only) then
!~       if (.not. allocated(ph_g)) then
!~          allocate(  ph_g(num_bands),stat=ierr )
!~          if (ierr/=0) call io_error('Error in allocating ph_g in dis_windows')
!~          ph_g = cmplx_1
!~       endif
!~       ! Apply same operation to ph_g
!~       do i = 1, ndimwin(1)
!~          j = nfirstwin(1) + i - 1
!~          ph_g(i) = ph_g(j)
!~       enddo
!~       do i = ndimwin(1) + 1, num_bands
!~          ph_g(i) = cmplx_0
!~       enddo
!~    endif
!~![ysl-e]

    if (iprint > 1) then
      if (on_root) write (stdout, '(1x,a)') &
        '|                        K-points with Frozen States                         |'
      if (on_root) write (stdout, '(1x,a)') &
        '|                        ---------------------------                         |'
      i = 0
      do nkp = 1, num_kpts
        if (ndimfroz(nkp) .gt. 0) then
          i = i + 1
          if (i .eq. 1) then
            if (on_root) write (stdout, '(1x,a,i6)', advance='no') '|', nkp
          else if ((i .gt. 1) .and. (i .lt. 12)) then
            if (on_root) write (stdout, '(i6)', advance='no') nkp
          else if (i .eq. 12) then
            if (on_root) write (stdout, '(i6,a)') nkp, '    |'
            i = 0
          endif
        endif
      enddo
      if (i .ne. 0) then
        do j = 1, 12 - i
          if (on_root) write (stdout, '(6x)', advance='no')
        enddo
        if (on_root) write (stdout, '(a)') '    |'
      endif
      if (on_root) write (stdout, '(1x,a)') &
        '+----------------------------------------------------------------------------+'
    endif

    if (on_root) write (stdout, '(3x,a,i4)') 'Number of target bands to extract: ', num_wann
    if (iprint > 1) then
      if (on_root) write (stdout, '(1x,a)') &
        '+----------------------------------------------------------------------------+'
      if (on_root) write (stdout, '(1x,a)') &
        '|                                  Windows                                   |'
      if (on_root) write (stdout, '(1x,a)') &
        '|                                  -------                                   |'
      if (on_root) write (stdout, '(1x,a)') &
        '|               K-point      Ndimwin     Ndimfroz    Nfirstwin               |'
      if (on_root) write (stdout, '(1x,a)') &
        '|               ----------------------------------------------               |'
      do nkp = 1, num_kpts
        if (on_root) write (stdout, 403) nkp, ndimwin(nkp), ndimfroz(nkp), nfirstwin(nkp)
      enddo
403   format(1x, '|', 14x, i6, 7x, i6, 7x, i6, 6x, i6, 18x, '|')
      if (on_root) write (stdout, '(1x,a)') &
        '+----------------------------------------------------------------------------+'
    endif

    if (timing_level > 1) call io_stopwatch('dis: windows', 2)

    return

  end subroutine dis_windows

  !==================================================================!
  subroutine dis_project()
    !==================================================================!
    !                                                                  !
    !! Construct projections for the start of the disentanglement routine
    !!
    !! Original notes from Nicola (refers only to the square case)
    !!
    !! This subroutine calculates the transformation matrix
    !! CU = CS^(-1/2).CA, where CS = CA.CA^dagger.
    !! CS is diagonalized with a Schur factorization, to be on the safe
    !! side of numerical stability.
    !!
    !! ZGEES computes for an N-by-N complex nonsymmetric matrix Y, the
    !! eigenvalues, the Schur form T, and, optionally, the matrix of
    !! Schur vectors. This gives the Schur factorization Y = Z*T*(Z**H).
    !!
    !! Optionally, it also orders the eigenvalues on the diagonal of the
    !! Schur form so that selected eigenvalues are at the top left.
    !! The leading components of Z then form an orthonormal basis for
    !! the invariant subspace corresponding to the selected eigenvalues.
    !!
    !! A complex matrix is in Schur form if it is upper triangular.
    !!
    !! Notes from Ivo disentangling (e.g. non-square) projection
    !! (See Sec. III.D of SMV paper)
    !! Compute the ndimwin(k) x num_wann matrix cu that yields,
    !! from the ndimwin original Bloch states, the num_wann Bloch-like
    !! states with maximal projection onto the num_wann localised
    !! functions:
    !!
    !! CU = CA.CS^{-1/2}, CS = transpose(CA).CA
    !!
    !! Use the singular-calue decomposition of the matrix CA:
    !!
    !! CA = CZ.CD.CVdagger (note: zgesvd spits out CVdagger)
    !!
    !! which yields
    !!
    !! CU = CZ.CD.CD^{-1}.CVdag
    !!
    !! where CZ is ndimwin(NKP) x ndimwin(NKP) and unitary, CD is
    !! ndimwin(NKP) x num_wann and diagonal, CD^{-1} is
    !! num_wann x num_wann and diagonal, and CVdag is
    !! num_wann x num_wann and unitary.
    !!
    !==================================================================!

    use w90_constants, only: eps5

    implicit none

    ! internal variables
    integer :: i, j, l, m, nkp, info, ierr
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp), allocatable :: rwork(:)
    complex(kind=dp)              :: ctmp2
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cvdag(:, :)
!    complex(kind=dp), allocatable :: catmpmat(:,:,:)

    if (timing_level > 1) call io_stopwatch('dis: project', 1)

    if (on_root) write (stdout, '(/1x,a)') &
      '                  Unitarised projection of Wannier functions                  '
    if (on_root) write (stdout, '(1x,a)') &
      '                  ------------------------------------------                  '
    if (on_root) write (stdout, '(3x,a)') 'A_mn = <psi_m|g_n> --> S = A.A^+ --> U = S^-1/2.A'
    if (on_root) write (stdout, '(3x,a)', advance='no') 'In dis_project...'

!    allocate(catmpmat(num_bands,num_bands,num_kpts),stat=ierr)
!    if (ierr/=0) call io_error('Error in allocating catmpmat in dis_project')
    allocate (svals(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating svals in dis_project')
    allocate (rwork(5*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rwork in dis_project')
    allocate (cvdag(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cvdag in dis_project')
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cz in dis_project')
    allocate (cwork(4*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwork in dis_project')

    ! here we slim down the ca matrix
    ! up to here num_bands(=num_bands) X num_wann(=num_wann)
!    do nkp = 1, num_kpts
!       do j = 1, num_wann
!          do i = 1, ndimwin(nkp)
!             catmpmat(i,j,nkp) = a_matrix(nfirstwin(nkp)+i-1,j,nkp)
!          enddo
!       enddo
!       do j = 1, num_wann
!          a_matrix(1:ndimwin(nkp),j,nkp) = catmpmat(1:ndimwin(nkp),j,nkp)
!       enddo
!       do j = 1, num_wann
!          a_matrix(ndimwin(nkp)+1:num_bands,j,nkp) = cmplx_0
!       enddo
!    enddo
    ! in order to reduce the memory usage, we don't use catmpmat.
    do nkp = 1, num_kpts
      if (ndimwin(nkp) .ne. num_bands) then
        do j = 1, num_wann
          do i = 1, ndimwin(nkp)
            ctmp2 = a_matrix(nfirstwin(nkp) + i - 1, j, nkp)
            a_matrix(i, j, nkp) = ctmp2
          enddo
          a_matrix(ndimwin(nkp) + 1:num_bands, j, nkp) = cmplx_0
        enddo
      endif
    enddo

    do nkp = 1, num_kpts
      ! SINGULAR VALUE DECOMPOSITION
      call ZGESVD('A', 'A', ndimwin(nkp), num_wann, a_matrix(:, :, nkp), &
                  num_bands, svals, cz, num_bands, cvdag, num_bands, cwork, &
                  4*num_bands, rwork, info)
      if (info .ne. 0) then
        if (on_root) write (stdout, *) ' ERROR: IN ZGESVD IN dis_project'
        if (on_root) write (stdout, *) ' K-POINT NKP=', nkp, ' INFO=', info
        if (info .lt. 0) then
          if (on_root) write (stdout, *) ' THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
        endif
        call io_error('dis_project: problem in ZGESVD 1')
      endif

      ! NOTE THAT - AT LEAST FOR LINUX MKL LAPACK - THE OUTPUT OF ZGESVD
      ! GIVES ALREADY Vdagger, SO A=Z.S.Vdagger IS ACTUALLY GIVEN BY cz.s.cvda
      !
      ! also, cu is cz.cd.cd^-1.cvdag, and the asymmetric structure of cd.cd^-
      ! allows us to say cu=cz.cvdag, but where the sum on the inner index
      ! goes only from 1 to num_wann (i.e. DO l=1,num_wann). This is because
      ! cz.cd.cd^-1 is a matrix ndimwin(nkp) X num_wann, identical to the
      ! first num_wann columns of the ndimwin(nkp) X ndimwin(nkp) matrix cz
      !
      ! same for ca: s is a ndimwin(nkp) X num_wann matrix that is
      ! zero everywhere but for its num_wann X num_wann top square part,
      ! that is diagonal. Multiplying cz by the s matrix is equivalent
      ! to moltiplying the first num_wann columns of cz, each by the correspondin
      ! diagonal element of s, that is s(L)
      ! I'm not sure why we reconstruct ca in what follows - in one explicit t
      ! [ aam: it is because a_matrix is overwritten by ZGESVD ]
      ! it seemed to be identical to the input ca (as it should be)

      u_matrix_opt(:, :, nkp) = cmplx_0
      a_matrix(:, :, nkp) = cmplx_0
      do j = 1, num_wann
        do i = 1, ndimwin(nkp)
          do l = 1, num_wann
            u_matrix_opt(i, j, nkp) = u_matrix_opt(i, j, nkp) + cz(i, l)*cvdag(l, j)
            a_matrix(i, j, nkp) = a_matrix(i, j, nkp) + cz(i, l)*svals(l)*cvdag(l, j)
          enddo
        enddo
      enddo

      !
      ! CHECK UNITARITY
      !
      ! note that cu.transpose(cu) is *NOT* an identity ndimwin(nkp) by ndimwi
      ! matrix, but transpose(cu).cu is a num_wann by num_wann identity matrix.
      ! I have once checked the former statement, now I will just leave here t
      ! for the latter (what this means is that the columns of cu are orthonor
      ! vectors).
      do i = 1, num_wann
        do j = 1, num_wann
          ctmp2 = cmplx_0
          do m = 1, ndimwin(nkp)
            ctmp2 = ctmp2 + u_matrix_opt(m, j, nkp)*conjg(u_matrix_opt(m, i, nkp))
          enddo
          if ((i .eq. j) .and. (abs(ctmp2 - cmplx_1) .gt. eps5)) then
            if (on_root) write (stdout, *) ' ERROR: unitarity of initial U'
            if (on_root) write (stdout, '(1x,a,i2)') 'nkp= ', nkp
            if (on_root) write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
            if (on_root) write (stdout, '(1x,a,f12.6,1x,f12.6)') &
              '[u_matrix_opt.transpose(u_matrix_opt)]_ij= ', &
              real(ctmp2, dp), aimag(ctmp2)
            call io_error('dis_project: Error in unitarity of initial U in dis_project')
          endif
          if ((i .ne. j) .and. (abs(ctmp2) .gt. eps5)) then
            if (on_root) write (stdout, *) ' ERROR: unitarity of initial U'
            if (on_root) write (stdout, '(1x,a,i2)') 'nkp= ', nkp
            if (on_root) write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
            if (on_root) write (stdout, '(1x,a,f12.6,1x,f12.6)') &
              '[u_matrix_opt.transpose(u_matrix_opt)]_ij= ', &
              real(ctmp2, dp), aimag(ctmp2)
            call io_error('dis_project: Error in unitarity of initial U in dis_project')
          endif
        enddo
      enddo
    enddo
    ! NKP

    deallocate (cwork, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cwork in dis_project')
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cz in dis_project')
    deallocate (cvdag, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cvdag in dis_project')
    deallocate (rwork, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rwork in dis_project')
    deallocate (svals, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating svals in dis_project')
!    deallocate(catmpmat,stat=ierr)
!    if (ierr/=0) call io_error('Error in deallocating catmpmat in dis_project')

    if (on_root) write (stdout, '(a)') ' done'

    if (timing_level > 1) call io_stopwatch('dis: project', 2)

    return

  end subroutine dis_project

  !==================================================================!
  subroutine dis_proj_froz()
    !==================================================================!
    !                                                                  !
    !! COMPUTES THE LEADING EIGENVECTORS OF Q_froz . P_s . Q_froz,
    !! WHERE P_s PROJECTOR OPERATOR ONTO THE SUBSPACE S OF THE PROJECTED
    !! GAUSSIANS, P_f THE PROJECTOR ONTO THE FROZEN STATES, AND
    !! Q_froz = 1 - P_froz, ALL EXP IN THE BASIS OF THE BLOCH
    !! EIGENSTATES INSIDE THE OUTER ENERGY WINDOW
    !! (See Eq. (27) in Sec. III.G of SMV)
    !                                                                  !
    !==================================================================!

    use w90_constants, only: eps8

    implicit none

    ! INPUT: num_wann,ndimwin,ndimfroz,indxfroz,lfrozen
    ! MODIFIED: u_matrix_opt (At input it contains the gaussians projected onto
    !             the window states in the routine project.f. At output
    !             the entries with the second index from 1 to ndimfroz(nkp)
    !             contain the frozen (inner window) states, while those
    !             from ndimfroz(nkp)+1 to num_wann have been replaced by
    !             the new trial states outside the inner window.)

    ! *********************************************************
    ! VARIABLES USED BY LAPACK'S ZHPEVX DIAGONALIZATION ROUTINE
    ! *********************************************************
    integer, allocatable :: iwork(:)
    integer, allocatable :: ifail(:)
    real(kind=dp), allocatable :: w(:)
    real(kind=dp), allocatable :: rwork(:)
    complex(kind=dp), allocatable :: cap(:)
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)

    ! *********
    ! INTERNAL:
    ! *********
    !
    ! CP_S(M,N)      PROJECTION OPERATOR ONTO THE SUBSPACE OF THE PROJEC
    !                  GAUSSIANS, EXPRESSED IN THE BASIS OF THE ORIGINAL BL
    !                  EIGENSTATES INSIDE THE OUTER WINDOW (FOR THE PRESENT
    !                  K-POINT)
    ! CQ_FROZ(M,N)   PROJECTION OPERATOR ONTO THE SUBSPACE OF THE STATES
    !                  THE SPACE OF FROZEN STATES (BUT INSIDE THE OUTER WIN
    !                  EXPRESSED IN THE BASIS OF THE ORIGINAL BLOCH EIGENST
    !                  INSIDE THE OUTER WINDOW (FOR THE PRESENT K-POINT)
    ! CPQ(M,N)       THE MATRIX cp_s . cq_froz FOR THE PRESENT K-POINT
    ! CQPQ(M,N)      THE MATRIX cq_froz . cp_s . cq_froz FOR THE PRESENT
    !

    integer :: goods, il, iu, nkp, l, j, n, m, info, ierr
    integer :: counter, loop_f, loop_v, vmap(num_bands)
    integer :: nzero
    logical :: take
    character(len=4) :: rep
    complex(kind=dp) :: ctmp
    complex(kind=dp), allocatable :: cp_s(:, :)
    complex(kind=dp), allocatable :: cq_froz(:, :)
    complex(kind=dp), allocatable :: cpq(:, :)
    complex(kind=dp), allocatable :: cqpq(:, :)

    if (timing_level > 1) call io_stopwatch('dis: proj_froz', 1)

    if (on_root) write (stdout, '(3x,a)', advance='no') 'In dis_proj_froz...'

    allocate (iwork(5*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating iwork in dis_proj_froz')
    allocate (ifail(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating ifail in dis_proj_froz')
    allocate (w(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating w in dis_proj_froz')
    allocate (rwork(7*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating rwork in dis_proj_froz')
    allocate (cap((num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cap in dis_proj_froz')
    allocate (cwork(2*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cwork in dis_proj_froz')
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cz in dis_proj_froz')

    allocate (cp_s(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cp_s in dis_proj_froz')
    allocate (cq_froz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cq_froz in dis_proj_froz')
    allocate (cpq(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cpq in dis_proj_froz')
    allocate (cqpq(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cqpq in dis_proj_froz')

    do nkp = 1, num_kpts

      ! aam: this should be done at the end, otherwise important
      !      projection info is lost
!~         ! Put the frozen states in the lowest columns of u_matrix_opt
!~         if (ndimfroz(nkp).gt.0) then
!~            do l = 1, ndimfroz(nkp)
!~               u_matrix_opt(:,l,nkp)=cmplx_0
!~               u_matrix_opt(indxfroz(l,nkp),l,nkp) = cmplx_1
!~            enddo
!~         endif

      ! If there are non-frozen states, compute the num_wann-ndimfroz(nkp) leadin
      ! eigenvectors of cqpq
      if (num_wann .gt. ndimfroz(nkp)) then
        cq_froz = cmplx_0
        cp_s = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, ndimwin(nkp)
            do l = 1, num_wann
              cp_s(m, n) = cp_s(m, n) + u_matrix_opt(m, l, nkp)*conjg(u_matrix_opt(n, l, nkp))
            enddo
          enddo
          if (.not. lfrozen(n, nkp)) cq_froz(n, n) = cmplx_1
        enddo

        cpq = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, ndimwin(nkp)
            do l = 1, ndimwin(nkp)
              cpq(m, n) = cpq(m, n) + cp_s(m, l)*cq_froz(l, n)
            enddo
          enddo
        enddo

        cqpq = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, ndimwin(nkp)
            do l = 1, ndimwin(nkp)
              cqpq(m, n) = cqpq(m, n) + cq_froz(m, l)*cpq(l, n)
            enddo
          enddo
        enddo

        ! DEBUG
        ! check hermiticity of cqpq
        do n = 1, ndimwin(nkp)
          do m = 1, n
            if (abs(cqpq(m, n) - conjg(cqpq(n, m))) .gt. eps8) then
              if (on_root) write (stdout, *) ' matrix CQPQ is not hermitian'
              if (on_root) write (stdout, *) ' k-point ', nkp
              call io_error('dis_proj_froz: error')
            endif
          enddo
        enddo
        ! ENDDEBUG

        cap = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, n
            cap(m + (n - 1)*n/2) = cqpq(m, n)
          enddo
        enddo
        il = ndimwin(nkp) - (num_wann - ndimfroz(nkp)) + 1
        iu = ndimwin(nkp)
        call ZHPEVX('V', 'A', 'U', ndimwin(nkp), cap, 0.0_dp, 0.0_dp, il, &
                    iu, -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)

!~            write(stdout,*) 'w:'
!~            do n=1,ndimwin(nkp)
!~               write(stdout,'(f14.10)') w(n)
!~            enddo
!~            write(stdout,*) 'cz:'
!~            do n=1,ndimwin(nkp)
!~               write(stdout,'(6f12.8)') cz(n,il), cz(n,iu)
!~            enddo

        ! DEBUG
        if (info .lt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING CQPQ MATRIX'
          if (on_root) write (stdout, *) ' THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
          call io_error('dis_proj_frozen: error')
        elseif (info .gt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING CQPQ MATRIX'
          if (on_root) write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
          call io_error('dis_proj_frozen: error')
        endif
        ! ENDDEBUG

        ! DEBUG
        if (m .ne. ndimwin(nkp)) then
          if (on_root) write (stdout, *) ' *** ERROR *** in dis_proj_froz'
          if (on_root) write (stdout, *) ' Number of eigenvalues/vectors obtained is', &
            m, ' not equal to the number asked,', ndimwin(nkp)
          call io_error('dis_proj_frozen: error')
        endif
        ! ENDDEBUG

        ! DEBUG
        ! check that the eigenvalues are between 0 and 1
        if (iprint > 2) then
          if (on_root) write (stdout, '(/a,i3,a,i3,a,i3,a)') ' K-point ', nkp, ' ndimwin: ', &
            ndimwin(nkp), ' we want the', num_wann - ndimfroz(nkp), &
            ' leading eigenvector(s) of QPQ'
        endif
        do j = 1, ndimwin(nkp)
          if (iprint > 2 .and. on_root) write (stdout, '(a,i3,a,f16.12)') '  lambda(', j, ')=', w(j)
!~[aam]        if ( (w(j).lt.eps8).or.(w(j).gt.1.0_dp + eps8) ) then
          if ((w(j) .lt. -eps8) .or. (w(j) .gt. 1.0_dp + eps8)) then
            call io_error('dis_proj_frozen: error - Eigenvalues not between 0 and 1')
          endif
        enddo
        ! ENDDEBUG

        ! [ aam: sometimes the leading eigenvalues form a degenerate set that is
        !        of higher dimensionality than (num_wann - ndimfroz). May need to
        !        fix this at some point. ]

        ! For certain cases we have found that one of the required eigenvectors of cqpq
        ! has a zero eigenvalue (ie it forms a degenerate set with the frozen states).
        ! It depends on floating point math as whether this eigenvalue is positive
        ! or negative (ie +/- 1e-17). If it's positive everything is ok. If negative we
        ! can end up putting one of the frozen states into U_opt (and failing the later
        ! orthogonality check).
        ! This fix detects this situation. If applies we choose the eigenvectors by
        ! checking their orthogonality to the frozen states.
        ! === For version 1.0.1 we make this the default ===

        if (index(devel_flag, 'no-orth-fix') == 0) then
          nzero = 0; goods = 0
          do j = ndimwin(nkp), ndimwin(nkp) - (num_wann - ndimfroz(nkp)) + 1, -1
            if (w(j) < eps8) then
              nzero = nzero + 1
            else
              goods = goods + 1
            end if
          end do
          if (nzero > 0) then
            if (iprint > 2 .and. on_root) then
              write (stdout, *) ' '
              write (stdout, '(1x,a,i0,a)') 'An eigenvalue of QPQ is close to zero at kpoint ' &
                , nkp, '. Using safety check.'
              write (stdout, '(1x,a,i4,a,i4)') 'We must find ', nzero, &
                ' eigenvectors with zero eigenvalues out of a set of ', ndimwin(nkp) - goods
            endif
            !First lets put the 'good' states into vamp
            vmap = 0
            counter = 1
            do j = ndimwin(nkp), ndimwin(nkp) - goods + 1, -1
              vmap(counter) = j
              counter = counter + 1
            end do

            if (iprint > 2 .and. on_root) then
              do loop_f = 1, ndimwin(nkp)
                write (stdout, '(1x,a,i4,a,es13.6)') 'Eigenvector number', loop_f, '    Eigenvalue: ', w(loop_f)
                do loop_v = 1, ndimwin(nkp)
                  write (stdout, '(20x,2f12.8)') cz(loop_v, loop_f)
                end do
                write (stdout, *)
              end do
            end if

            ! We need to find nzero vectors out of the remining ndimwin(nkp)-goods vectors

            do loop_f = 1, nzero
              do loop_v = ndimwin(nkp), 1, -1 !loop backwards for efficiency only
                if (any(vmap == loop_v)) cycle
                !check to see if vector is orthogonal to frozen states in u_matrix_opt
                take = .true.
                do m = 1, ndimfroz(nkp)
                  ctmp = cmplx_0
                  do j = 1, ndimwin(nkp)
                    ctmp = ctmp + conjg(u_matrix_opt(j, m, nkp))*cz(j, loop_v)
                  enddo
                  if (abs(ctmp) .gt. eps8) then
                    take = .false.
                  endif
                enddo
                if (take) then !vector is good and we add it to vmap
                  vmap(goods + loop_f) = loop_v
                  exit
                end if
              end do
            end do

            if (iprint > 2 .and. on_root) then
              write (rep, '(i4)') num_wann - ndimfroz(nkp)
              write (stdout, '(1x,a,'//trim(rep)//'(i0,1x))') 'We use the following eigenvectors: ' &
                , vmap(1:(num_wann - ndimfroz(nkp)))
            end if
            do l = 1, num_wann - ndimfroz(nkp)
              if (vmap(l) == 0) call io_error('dis_proj_froz: Ortho-fix failed to find enough vectors')
            end do

            ! put the correct eigenvectors into u_matrix_opt, and we're all done!
            counter = 1
            do l = ndimfroz(nkp) + 1, num_wann
              u_matrix_opt(1:ndimwin(nkp), l, nkp) = cz(1:ndimwin(nkp), vmap(counter))
              counter = counter + 1
            enddo

          else ! we don't need to use the fix

            do l = ndimfroz(nkp) + 1, num_wann
              u_matrix_opt(1:ndimwin(nkp), l, nkp) = cz(1:ndimwin(nkp), il)
              il = il + 1
            enddo

            if (il - 1 .ne. iu) then
              call io_error('dis_proj_frozen: error -  il-1.ne.iu  (in ortho-fix)')
            endif

          end if

        else ! if .not. using ortho-fix

          ! PICK THE num_wann-nDIMFROZ(NKP) LEADING EIGENVECTORS AS TRIAL STATES
          ! and PUT THEM RIGHT AFTER THE FROZEN STATES IN u_matrix_opt
          do l = ndimfroz(nkp) + 1, num_wann
            if (on_root) write (stdout, *) 'il=', il
            u_matrix_opt(1:ndimwin(nkp), l, nkp) = cz(1:ndimwin(nkp), il)
            il = il + 1
          enddo

          ! DEBUG
          if (il - 1 .ne. iu) then
            call io_error('dis_proj_frozen: error -  il-1.ne.iu')
          endif
          ! ENDDEBUG

        end if

      endif   ! num_wann>nDIMFROZ(NKP)

      ! Put the frozen states in the lowest columns of u_matrix_opt
      if (ndimfroz(nkp) .gt. 0) then
        do l = 1, ndimfroz(nkp)
          u_matrix_opt(:, l, nkp) = cmplx_0
          u_matrix_opt(indxfroz(l, nkp), l, nkp) = cmplx_1
        enddo
      endif

!~         write(stdout,*) 'u_matrix_opt:'
!~         do m=1,ndimwin(nkp)
!~            write(stdout,'(6f12.8)') u_matrix_opt(m,1,nkp), &
!~                 u_matrix_opt(m,ndimfroz(nkp),nkp), u_matrix_opt(m,num_wann,nkp)
!~         enddo

    enddo   ! NKP

    deallocate (cqpq, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cqpq in dis_proj_froz')
    deallocate (cpq, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cpq in dis_proj_froz')
    deallocate (cq_froz, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cq_froz in dis_proj_froz')
    deallocate (cp_s, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cp_s in dis_proj_froz')

    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cz in dis_proj_froz')
    deallocate (cwork, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cwork in dis_proj_froz')
    deallocate (cap, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cap in dis_proj_froz')
    deallocate (rwork, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating rwork in dis_proj_froz')
    deallocate (w, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating w in dis_proj_froz')
    deallocate (ifail, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating ifail in dis_proj_froz')
    deallocate (iwork, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating iwork in dis_proj_froz')

    if (on_root) write (stdout, '(a)') ' done'

    if (timing_level > 1) call io_stopwatch('dis: proj_froz', 2)

    return

  end subroutine dis_proj_froz

  !==================================================================!
  subroutine dis_extract()
    !==================================================================!
    !                                                                  !
    !! Extracts an num_wann-dimensional subspace at each k by
    !! minimizing Omega_I
    !                                                                  !
    !==================================================================!

    use w90_io, only: io_wallclocktime
    use w90_sitesym, only: ir2ik, ik2ir, nkptirr, nsymmetry, kptsym !YN: RS:

    implicit none

    ! MODIFIED:
    !           u_matrix_opt (At input it contains the initial guess for the optima
    ! subspace (expressed in terms of the original states inside the window)
    ! output it contains the  states that diagonalize the hamiltonian inside
    ! optimal subspace (again expressed in terms of the original window stat
    ! Giving out states that diagonalize the hamiltonian inside the optimal
    ! subspace (instead of the eigenstates of the Z matrix) is useful for
    ! performing the Wannier interpolation of the energy bands as described
    ! Sec. III.F of SMV)
    !
    !           eigval (At input: original energy eigenvalues.
    ! At output: eigenvalues of the hamiltonian inside optimal subspace)

    ! ----------------------------------------------------------------------
    ! TO DO: The complement subspace is computed but is not saved anywhere!
    ! (Check what was done with it in original code space.f)
    ! Diagonalize Z matrix only at those k points where ndimwin>num_wann?
    ! ----------------------------------------------------------------------

    ! *******************
    ! SHELLS OF K-VECTORS
    ! *******************
    ! nshells           number of shells of k-points to be used in the
    !                   finite-difference formulas for the k-derivatives
    ! aam: wb is now wb(1:nntot) 09/04/2006
    ! wb(nkp,nnx)       weight of the nnx-th b-vector (ordered along shells
    !                   of increasing length) associated with the nkp-th k-p
    ! wbtot             sum of the weights of all b-vectors associated with
    !                   given k-point (k-point 1 is used in calculation)
    ! nnlist(nkp,nnx)   vkpt(1:3,nnlist(nkp,nnx)) is the nnx-th neighboring
    !                   k-point of the nkp-th k-point vkpt(1:3,nkp) (or its
    !                   periodic image in the "home Brillouin zone")
    ! cm(n,m,nkp,nnx)   Overlap matrix <u_nk|u_{m,k+b}>

    ! Internal variables
    integer :: i, j, l, m, n, nn, nkp, nkp2, info, ierr, ndimk, p
    integer :: icompflag, iter, ndiff
    real(kind=dp) :: womegai, wkomegai, womegai1, rsum, delta_womegai
    real(kind=dp), allocatable :: wkomegai1(:)

    ! for MPI
    real(kind=dp), allocatable :: wkomegai1_loc(:)
    complex(kind=dp), allocatable :: camp_loc(:, :, :)
    complex(kind=dp), allocatable :: u_matrix_opt_loc(:, :, :)

    complex(kind=dp), allocatable :: ceamp(:, :, :)
    complex(kind=dp), allocatable :: camp(:, :, :)
    complex(kind=dp), allocatable :: czmat_in(:, :, :)
    complex(kind=dp), allocatable :: czmat_out(:, :, :)
    ! the z-matrices are now stored in local arrays
    complex(kind=dp), allocatable :: czmat_in_loc(:, :, :)
    complex(kind=dp), allocatable :: czmat_out_loc(:, :, :)
    complex(kind=dp), allocatable :: cham(:, :, :)

    integer, allocatable :: iwork(:)
    integer, allocatable :: ifail(:)
    real(kind=dp), allocatable :: w(:)
    real(kind=dp), allocatable :: rwork(:)
    complex(kind=dp), allocatable :: cap(:)
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)

    complex(kind=dp), allocatable :: cwb(:, :), cww(:, :), cbw(:, :)

    real(kind=dp), allocatable :: history(:)
    logical                       :: dis_converged
    complex(kind=dp) :: lambda(num_wann, num_wann) !RS:

    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer :: nkp_loc

    if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract', 1)

    if (on_root) write (stdout, '(/1x,a)') &
      '                  Extraction of optimally-connected subspace                  '
    if (on_root) write (stdout, '(1x,a)') &
      '                  ------------------------------------------                  '

    allocate (cwb(num_wann, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cwb in dis_extract')
    allocate (cww(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cww in dis_extract')
    allocate (cbw(num_bands, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cbw in dis_extract')
    cwb = cmplx_0; cww = cmplx_0; cbw = cmplx_0

    allocate (iwork(5*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating iwork in dis_extract')
    allocate (ifail(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating ifail in dis_extract')
    allocate (w(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating w in dis_extract')
    allocate (rwork(7*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating rwork in dis_extract')
    allocate (cap((num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cap in dis_extract')
    allocate (cwork(2*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cwork in dis_extract')
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cz in dis_extract')

    ! for MPI
    call comms_array_split(num_kpts, counts, displs)
    allocate (u_matrix_opt_loc(num_bands, num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating u_matrix_opt_loc in dis_extract')
    ! Copy matrix elements from global U matrix to local U matrix
    do nkp_loc = 1, counts(my_node_id)
      nkp = nkp_loc + displs(my_node_id)
      u_matrix_opt_loc(:, :, nkp_loc) = u_matrix_opt(:, :, nkp)
    enddo
    allocate (wkomegai1_loc(max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wkomegai1_loc in dis_extract')
    allocate (czmat_in_loc(num_bands, num_bands, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating czmat_in_loc in dis_extract')
    allocate (czmat_out_loc(num_bands, num_bands, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating czmat_out_loc in dis_extract')

    allocate (wkomegai1(num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wkomegai1 in dis_extract')
    allocate (czmat_in(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating czmat_in in dis_extract')
    allocate (czmat_out(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating czmat_out in dis_extract')

    allocate (history(dis_conv_window), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating history in dis_extract')

    ! ********************************************
    ! ENERGY WINDOWS AND SUBSPACES AT EACH K-POINT
    ! ********************************************
    ! num_wann             dimensionality of the subspace at each k-point
    !                   (number of Wannier functions per unit cell that we w
    ! NDIMWIN(NKP)      number of bands at the nkp-th k-point that fall
    !                   within the outer energy window
    ! NDIMFROZ(NKP)     number of frozen bands at the nkp-th k-point
    ! INDXNFROZ(I,NKP)  INDEX (BETWEEN 1 AND NDIMWIN(NKP)) OF THE I-TH NON-F
    !                   ORIGINAL BAND STATE AT THE NKP-TH K-POINT
    ! U_MATRIX_OPT(J,L,NKP)    AMPLITUDE OF THE J-TH ENERGY EIGENVECTOR INSIDE THE
    !                   ENERGY WINDOW AT THE NKP-TH K-POINT IN THE EXPANSION
    !                   THE L-TH LEADING RLAMBDA EIGENVECTOR AT THE SAME K-P
    !                   If there are M_k frozen states, they occupy the lowe
    !                   entries of the second index of u_matrix_opt, and the leadin
    !                   nabnds-M_k eigenvectors of the Z matrix occupy the
    !                   remaining slots
    ! CAMP(J,L,NKP)     SAME AS U_MATRIX_OPT, BUT FOR THE COMPLEMENT SUBSPACE INSID
    !                   ENERGY WINDOW (I.E., THE NON-LEADING RLAMBDA EIGENVE
    ! CEAMP(J,L,NKPTS)  SAME AS U_MATRIX_OPT, BUT INSTEAD OF RLAMBDA EIGENVECTOR, I
    !                   FOR THE ENERGY EIGENVECTOR OBTAINED BY DIAGONALIZING
    !                   HAMILTONIAN IN THE OPTIMIZED SUBSPACE
    ! CZMAT_IN(M,N,NKP) Z-MATRIX [Eq. (21) SMV]
    ! CZMAT_OUT(M,N,NKP) OUTPUT Z-MATRIX FROM THE PRESENT ITERATION
    ! RLAMBDA           An eigenvalue of the Z matrix
    ! womegai           Gauge-invariant Wannier spread, computed usinf all s
    !                   from current iteration
    ! wkomegai1(NKP)    Eq. (18) of SMV
    ! womegai1          Eq.(11) of SMV (like wowmegai, but neighboring state
    !                   for computing overlaps are from previous iteration.
    !                   become equal at self-consistency)
    ! alphafixe         mixing parameter for the iterative procedure
    ! nitere            total number of iterations

    ! DEBUG
    if (iprint > 2) then
      if (on_root) write (stdout, '(a,/)') '  Original eigenvalues inside outer window:'
      do nkp = 1, num_kpts
        if (on_root) write (stdout, '(a,i3,3x,20(f9.5,1x))') '  K-point ', nkp, &
          (eigval_opt(i, nkp), i=1, ndimwin(nkp))
      enddo
    endif
    ! ENDDEBUG

    ! TO DO: Check if this is the best place to initialize icompflag
    icompflag = 0

    if (on_root) write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'
    if (on_root) write (stdout, '(1x,a)') &
      '|  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS'
    if (on_root) write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'

    dis_converged = .false.

    ! ------------------
    ! BIG ITERATION LOOP
    ! ------------------
    do iter = 1, dis_num_iter

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_1', 1)

      if (iter .eq. 1) then
        ! Initialize Z matrix at k points w/ non-frozen states
        do nkp_loc = 1, counts(my_node_id)
          nkp = nkp_loc + displs(my_node_id)
          if (num_wann .gt. ndimfroz(nkp)) call internal_zmatrix(nkp, nkp_loc, czmat_in_loc(:, :, nkp_loc))
        enddo

        if (lsitesymmetry) then
          call comms_gatherv(czmat_in_loc, num_bands*num_bands*counts(my_node_id), &
                             czmat_in, num_bands*num_bands*counts, num_bands*num_bands*displs)
          call comms_bcast(czmat_in(1, 1, 1), num_bands*num_bands*num_kpts)
          call sitesym_symmetrize_zmatrix(czmat_in, lwindow) !RS:
          do nkp_loc = 1, counts(my_node_id)
            nkp = nkp_loc + displs(my_node_id)
            czmat_in_loc(:, :, nkp_loc) = czmat_in(:, :, nkp)
          end do
        end if

      else
        ! [iter.ne.1]
        ! Update Z matrix at k points with non-frozen states, using a mixing sch
        do nkp_loc = 1, counts(my_node_id)
          nkp = nkp_loc + displs(my_node_id)
          if (lsitesymmetry) then                !YN: RS:
            if (ir2ik(ik2ir(nkp)) .ne. nkp) cycle !YN: RS:
          endif                                  !YN: RS:
          if (num_wann .gt. ndimfroz(nkp)) then
            ndimk = ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, ndimk
              do j = 1, i
                czmat_in_loc(j, i, nkp_loc) = &
                  cmplx(dis_mix_ratio, 0.0_dp, dp)*czmat_out_loc(j, i, nkp_loc) &
                  + cmplx(1.0_dp - dis_mix_ratio, 0.0_dp, dp)*czmat_in_loc(j, i, nkp_loc)
                ! hermiticity
                czmat_in_loc(i, j, nkp_loc) = conjg(czmat_in_loc(j, i, nkp_loc))
              enddo
            enddo
          endif
        enddo
      endif
      ! [if iter=1]
      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_1', 2)

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_2', 1)

      womegai1 = 0.0_dp
      ! wkomegai1 is defined by Eq. (18) of SMV.
      ! Contribution to wkomegai1 from frozen states should be calculated now
      ! every k (before updating any k), so that for iter>1 overlaps are with
      ! non-frozen neighboring states from the previous iteration

      wkomegai1 = real(num_wann, dp)*wbtot
      if (lsitesymmetry) then                                                                        !RS:
        do nkp = 1, nkptirr                                                                            !RS:
          wkomegai1(ir2ik(nkp)) = wkomegai1(ir2ik(nkp))*nsymmetry/count(kptsym(:, nkp) .eq. ir2ik(nkp)) !RS:
        enddo                                                                                       !RS:
      endif                                                                                          !RS:
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        wkomegai1_loc(nkp_loc) = wkomegai1(nkp)
      end do

      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        if (ndimfroz(nkp) .gt. 0) then
          if (lsitesymmetry) call io_error('not implemented in symmetry-adapted mode') !YN: RS:
          do nn = 1, nntot
            nkp2 = nnlist(nkp, nn)
            call zgemm('C', 'N', ndimfroz(nkp), ndimwin(nkp2), ndimwin(nkp), cmplx_1, &
                       u_matrix_opt(:, :, nkp), num_bands, m_matrix_orig_local(:, :, nn, nkp_loc), num_bands, cmplx_0, &
                       cwb, num_wann)
            call zgemm('N', 'N', ndimfroz(nkp), num_wann, ndimwin(nkp2), cmplx_1, &
                       cwb, num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
            rsum = 0.0_dp
            do n = 1, num_wann
              do m = 1, ndimfroz(nkp)
                rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
              enddo
            enddo
            wkomegai1_loc(nkp_loc) = wkomegai1_loc(nkp_loc) - wb(nn)*rsum
          enddo
        endif
      enddo
      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_2', 2)

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_3', 1)

      !! ! send chunks of wkomegai1 to root node
      !! call comms_gatherv(wkomegai1_loc, counts(my_node_id), wkomegai1, counts, displs)
      !! ! send back the whole wkomegai1 array to other nodes
      !! call comms_bcast(wkomegai1(1), num_kpts)

      ! Refine optimal subspace at k points w/ non-frozen states
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        if (lsitesymmetry) then                                                     !RS:
          if (ir2ik(ik2ir(nkp)) .ne. nkp) cycle                                      !RS:
        end if                                                                      !RS:
        if (lsitesymmetry) then                                                     !RS:

          call sitesym_dis_extract_symmetry(nkp, ndimwin(nkp), czmat_in_loc(:, :, nkp_loc), lambda, u_matrix_opt_loc(:, :, nkp_loc)) !RS:

          do j = 1, num_wann                                                          !RS:
            wkomegai1_loc(nkp_loc) = wkomegai1_loc(nkp_loc) - real(lambda(j, j), kind=dp)               !RS:
          enddo                                                                    !RS:
        else                                                                        !RS:
          if (num_wann .gt. ndimfroz(nkp)) then
            ! Diagonalize Z matrix
            do j = 1, ndimwin(nkp) - ndimfroz(nkp)
              do i = 1, j
                cap(i + ((j - 1)*j)/2) = czmat_in_loc(i, j, nkp_loc)
              enddo
            enddo
            ndiff = ndimwin(nkp) - ndimfroz(nkp)
            call ZHPEVX('V', 'A', 'U', ndiff, cap, 0.0_dp, 0.0_dp, 0, 0, &
                        -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)
            if (info .lt. 0) then
              if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING Z MATRIX'
              if (on_root) write (stdout, *) ' THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
              call io_error(' dis_extract: error')
            endif
            if (info .gt. 0) then
              if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING Z MATRIX'
              if (on_root) write (stdout, *) info, ' EIGENVECTORS FAILED TO CONVERGE'
              call io_error(' dis_extract: error')
            endif

            ! Update the optimal subspace by incorporating the num_wann-ndimfroz(nkp) l
            ! eigenvectors of the Z matrix into u_matrix_opt. Also, add contribution from
            ! non-frozen states to wkomegai1(nkp) (minus the corresponding eigenvalu
            m = ndimfroz(nkp)
            do j = ndimwin(nkp) - num_wann + 1, ndimwin(nkp) - ndimfroz(nkp)
              m = m + 1
              wkomegai1_loc(nkp_loc) = wkomegai1_loc(nkp_loc) - w(j)
              u_matrix_opt_loc(1:ndimwin(nkp), m, nkp_loc) = cmplx_0
              ndimk = ndimwin(nkp) - ndimfroz(nkp)
              do i = 1, ndimk
                p = indxnfroz(i, nkp)
                u_matrix_opt_loc(p, m, nkp_loc) = cz(i, j)
              enddo
            enddo
          endif
          ! [if num_wann>ndimfroz(nkp)]
        endif !RS:

        ! Now that we have contribs. from both frozen and non-frozen states to
        ! wkomegai1(nkp), add it to womegai1
        womegai1 = womegai1 + wkomegai1_loc(nkp_loc)

        if (index(devel_flag, 'compspace') > 0) then

          ! AT THE LAST ITERATION FIND A BASIS FOR THE (NDIMWIN(NKP)-num_wann)-DIMENS
          ! COMPLEMENT SPACE

          if (iter .eq. dis_num_iter) then
            allocate (camp(num_bands, num_bands, num_kpts), stat=ierr)
            if (ierr /= 0) call io_error('Error allocating camp in dis_extract')
            allocate (camp_loc(num_bands, num_bands, max(1, counts(my_node_id))), stat=ierr)
            if (ierr /= 0) call io_error('Error allocating ucamp_loc in dis_extract')

            if (ndimwin(nkp) .gt. num_wann) then
              do j = 1, ndimwin(nkp) - num_wann
                if (num_wann .gt. ndimfroz(nkp)) then
                  ! USE THE NON-LEADING EIGENVECTORS OF THE Z-MATRIX
                  camp_loc(1:ndimwin(nkp), j, nkp_loc) = cz(1:ndimwin(nkp), j)
                else
                  ! Then num_wann=NDIMFROZ(NKP)
                  ! USE THE ORIGINAL NON-FROZEN BLOCH EIGENSTATES
                  do i = 1, ndimwin(nkp)
                    camp_loc(i, j, nkp_loc) = cmplx_0
                    if (i .eq. indxnfroz(j, nkp)) camp_loc(i, j, nkp_loc) = cmplx_1
                  enddo
                endif
              enddo
            else
              icompflag = 1
            endif
          endif

        end if ! index(devel_flag,'compspace')>0

      enddo
      ! [Loop over k points (nkp)]

      !! ! send chunks of wkomegai1 to root node
      !! call comms_gatherv(wkomegai1_loc, counts(my_node_id), wkomegai1, counts, displs)
      !! ! send back the whole wkomegai1 array to other nodes
      !! call comms_bcast(wkomegai1(1), num_kpts)

      call comms_allreduce(womegai1, 1, 'SUM')

      call comms_gatherv(u_matrix_opt_loc, num_bands*num_wann*counts(my_node_id), &
                         u_matrix_opt, num_bands*num_wann*counts, num_bands*num_wann*displs)
      call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts)
      if (lsitesymmetry) call sitesym_symmetrize_u_matrix(num_bands, u_matrix_opt, lwindow) !RS:

      if (index(devel_flag, 'compspace') > 0) then
        if (iter .eq. dis_num_iter) then
          call comms_gatherv(camp_loc, num_bands*num_bands*counts(my_node_id), &
                             camp, num_bands*num_bands*counts, num_bands*num_bands*displs)

          call comms_bcast(camp(1, 1, 1), num_bands*num_bands*num_kpts)
        endif
      endif

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_3', 2)

      womegai1 = womegai1/real(num_kpts, dp)

      ! DEBUG
      ! Orthonormality check
      !         do nkp=1,nkpts
      !           write(*,*) ' '
      !           write(*,'(a8,i4)') 'k-point ',nkp
      !           do l=1,num_wann
      !           do m=1,l
      !             ctmp=czero
      !             do j=1,ndimwin(nkp)
      !               ctmp=ctmp+conjg(u_matrix_opt(j,m,nkp))*u_matrix_opt(j,l,nkp)
      !             enddo
      !             write(*,'(i2,2x,i2,f16.12,1x,f16.12)') l,m,ctmp
      !             if(l.eq.m) then
      !               if(abs(ctmp-cmplx(1.0d0,0.0d0)).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             else
      !               if(abs(ctmp).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             endif
      !           enddo
      !           enddo
      !         enddo
      ! ENDDEBUG

      ! Compute womegai  using the updated subspaces at all k, i.e.,
      ! replacing (i-1) by (i) in Eq. (12) SMV
      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_4', 1)

      womegai = 0.0_dp
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        wkomegai = 0.0_dp
        do nn = 1, nntot
          nkp2 = nnlist(nkp, nn)
          call zgemm('C', 'N', num_wann, ndimwin(nkp2), ndimwin(nkp), cmplx_1, &
                     u_matrix_opt(:, :, nkp), num_bands, m_matrix_orig_local(:, :, nn, nkp_loc), num_bands, cmplx_0, &
                     cwb, num_wann)
          call zgemm('N', 'N', num_wann, num_wann, ndimwin(nkp2), cmplx_1, &
                     cwb, num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
          rsum = 0.0_dp
          do n = 1, num_wann
            do m = 1, num_wann
              rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
            enddo
          enddo
          wkomegai = wkomegai + wb(nn)*rsum
        enddo
        wkomegai = real(num_wann, dp)*wbtot - wkomegai
        womegai = womegai + wkomegai
      enddo

      call comms_allreduce(womegai, 1, 'SUM')

      womegai = womegai/real(num_kpts, dp)
      ! [Loop over k (nkp)]
      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract_4', 2)

      delta_womegai = womegai1/womegai - 1.0_dp

      if (on_root) write (stdout, 124) iter, womegai1*lenconfac**2, womegai*lenconfac**2, &
        delta_womegai, io_wallclocktime()

124   format(2x, i6, 3x, f14.8, 3x, f14.8, 6x, es10.3, 2x, f8.2, 4x, '<-- DIS')

      ! Construct the updated Z matrix, CZMAT_OUT, at k points w/ non-frozen s
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        if (num_wann .gt. ndimfroz(nkp)) call internal_zmatrix(nkp, nkp_loc, czmat_out_loc(:, :, nkp_loc))
      enddo

      if (lsitesymmetry) then
        call comms_gatherv(czmat_out_loc, num_bands*num_bands*counts(my_node_id), &
                           czmat_out, num_bands*num_bands*counts, num_bands*num_bands*displs)
        call comms_bcast(czmat_out(1, 1, 1), num_bands*num_bands*num_kpts)
        call sitesym_symmetrize_zmatrix(czmat_out, lwindow) !RS:
        do nkp_loc = 1, counts(my_node_id)
          nkp = nkp_loc + displs(my_node_id)
          czmat_out_loc(:, :, nkp_loc) = czmat_out(:, :, nkp)
        end do
      end if

      call internal_test_convergence()

      if (dis_converged) then
        if (on_root) write (stdout, '(/13x,a,es10.3,a,i2,a)') &
          '<<<      Delta <', dis_conv_tol, &
          '  over ', dis_conv_window, ' iterations     >>>'
        if (on_root) write (stdout, '(13x,a)') '<<< Disentanglement convergence criteria satisfied >>>'
        exit
      endif

    enddo
    ! [BIG ITERATION LOOP (iter)]

    deallocate (czmat_out, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating czmat_out in dis_extract')
    deallocate (czmat_in, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating czmat_in in dis_extract')
    deallocate (czmat_out_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating czmat_out_loc in dis_extract')
    deallocate (czmat_in_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating czmat_in_loc in dis_extract')

    if (on_root) then
      allocate (ceamp(num_bands, num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating ceamp in dis_extract')
      allocate (cham(num_bands, num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating cham in dis_extract')
    endif

    if (.not. dis_converged) then
      if (on_root) write (stdout, '(/5x,a)') &
        '<<< Warning: Maximum number of disentanglement iterations reached >>>'
      if (on_root) write (stdout, '(10x,a)') '<<< Disentanglement convergence criteria not satisfied >>>'
    endif

    if (index(devel_flag, 'compspace') > 0) then

      if (icompflag .eq. 1) then
        if (iprint > 2) then
          if (on_root) write (stdout, ('(/4x,a)')) &
            'WARNING: Complement subspace has zero dimensions at the following k-points:'
          i = 0
          if (on_root) write (stdout, '(4x)', advance='no')
          do nkp = 1, num_kpts
            if (ndimwin(nkp) .eq. num_wann) then
              i = i + 1
              if (i .le. 12) then
                if (on_root) write (stdout, '(i6)', advance='no') nkp
              else
                i = 1
                if (on_root) write (stdout, '(/4x)', advance='no')
                if (on_root) write (stdout, '(i6)', advance='no') nkp
              endif
            endif
          enddo
        endif
      endif

    endif

    ! Write the final womegai. This should remain unchanged during the
    ! subsequent minimization of Omega_tilde in wannierise.f90
    ! We store it in the checkpoint file as a sanity check
    if (on_root) write (stdout, '(/8x,a,f14.8,a/)') 'Final Omega_I ', &
      womegai*lenconfac**2, ' ('//trim(length_unit)//'^2)'

    ! Set public variable omega_invariant
    omega_invariant = womegai

    ! Currently, this part is not parallelized; thus, we perform the task only on root and then broadcast the result.
    if (on_root) then
      ! DIAGONALIZE THE HAMILTONIAN WITHIN THE OPTIMIZED SUBSPACES
      do nkp = 1, num_kpts

        do j = 1, num_wann
          do i = 1, num_wann
            cham(i, j, nkp) = cmplx_0
            do l = 1, ndimwin(nkp)
              cham(i, j, nkp) = cham(i, j, nkp) + conjg(u_matrix_opt(l, i, nkp)) &
                                *u_matrix_opt(l, j, nkp)*eigval_opt(l, nkp)
            enddo
          enddo
        enddo

        do j = 1, num_wann
          do i = 1, j
            cap(i + ((j - 1)*j)/2) = cham(i, j, nkp)
          enddo
        enddo

        call ZHPEVX('V', 'A', 'U', num_wann, cap, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                    m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)

        if (info .lt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
          if (on_root) write (stdout, *) ' THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
          call io_error(' dis_extract: error')
        endif
        if (info .gt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
          if (on_root) write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
          call io_error(' dis_extract: error')
        endif

        ! Store the energy eigenvalues of the optimal subspace (used in wann_ban
        eigval_opt(1:num_wann, nkp) = w(1:num_wann)

        ! CALCULATE AMPLITUDES OF THE CORRESPONDING ENERGY EIGENVECTORS IN TERMS
        ! THE ORIGINAL ("WINDOW SPACE") ENERGY EIGENVECTORS
        do j = 1, num_wann
          do i = 1, ndimwin(nkp)
            ceamp(i, j, nkp) = cmplx_0
            do l = 1, num_wann
              ceamp(i, j, nkp) = ceamp(i, j, nkp) + cz(l, j)*u_matrix_opt(i, l, nkp)
            enddo
          enddo
        enddo
        ! NKP
      enddo

      ! DEBUG
      if (iprint > 2) then
        if (on_root) write (stdout, '(/,a,/)') '  Eigenvalues inside optimal subspace:'
        do nkp = 1, num_kpts
          if (on_root) write (stdout, '(a,i3,2x,20(f9.5,1x))') '  K-point ', &
            nkp, (eigval_opt(i, nkp), i=1, num_wann)
        enddo
      endif
      ! ENDDEBUG

      ! Replace u_matrix_opt by ceamp. Both span the
      ! same space, but the latter is more convenient for the purpose of obtai
      ! an optimal Fourier-interpolated band structure: see Sec. III.E of SMV.
      if (.not. lsitesymmetry) then                                                                         !YN:
        do nkp = 1, num_kpts
          do j = 1, num_wann
            u_matrix_opt(1:ndimwin(nkp), j, nkp) = ceamp(1:ndimwin(nkp), j, nkp)
          enddo
        enddo
        !else                                                                                                        !YN:
        ! Above is skipped as we require Uopt(Rk) to be related to Uopt(k)                                        !YN: RS:
        !write(stdout,"(a)")  &                                                                                   !YN: RS:
        !  'Note(symmetry-adapted mode): u_matrix_opt are no longer the eigenstates of the subspace Hamiltonian.' !RS:
      endif                                                                                                        !YN:
    endif
    call comms_bcast(eigval_opt(1, 1), num_bands*num_kpts)
    call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts)

    if (index(devel_flag, 'compspace') > 0) then

      if (icompflag .eq. 1) then
        if (iprint > 2) then
          if (on_root) write (stdout, *) 'AT SOME K-POINT(S) COMPLEMENT SUBSPACE HAS ZERO DIMENSIONALITY'
          if (on_root) write (stdout, *) '=> DID NOT CREATE FILE COMPSPACE.DAT'
        endif
      else
        ! DIAGONALIZE THE HAMILTONIAN IN THE COMPLEMENT SUBSPACE, WRITE THE
        ! CORRESPONDING EIGENFUNCTIONS AND ENERGY EIGENVALUES
        do nkp = 1, num_kpts
          do j = 1, ndimwin(nkp) - num_wann
            do i = 1, ndimwin(nkp) - num_wann
              cham(i, j, nkp) = cmplx_0
              do l = 1, ndimwin(nkp)
                cham(i, j, nkp) = cham(i, j, nkp) + conjg(camp(l, i, nkp)) &
                                  *camp(l, j, nkp)*eigval_opt(l, nkp)
              enddo
            enddo
          enddo
          do j = 1, ndimwin(nkp) - num_wann
            do i = 1, j
              cap(i + ((j - 1)*j)/2) = cham(i, j, nkp)
            enddo
          enddo
          ndiff = ndimwin(nkp) - num_wann
          call ZHPEVX('V', 'A', 'U', ndiff, cap, 0.0_dp, 0.0_dp, 0, 0, &
                      -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)
          if (info .lt. 0) then
            if (on_root) write (stdout, *) '*** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
            if (on_root) write (stdout, *) 'THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
            call io_error(' dis_extract: error')
          endif
          if (info .gt. 0) then
            if (on_root) write (stdout, *) '*** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
            if (on_root) write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
            call io_error(' dis_extract: error')
          endif
          ! CALCULATE AMPLITUDES OF THE ENERGY EIGENVECTORS IN THE COMPLEMENT SUBS
          ! TERMS OF THE ORIGINAL ENERGY EIGENVECTORS
          do j = 1, ndimwin(nkp) - num_wann
            do i = 1, ndimwin(nkp)
              camp(i, j, nkp) = cmplx_0
              do l = 1, ndimwin(nkp) - num_wann
!write(stdout,*) 'i=',i,'   j=',j,'   l=',l
!write(stdout,*) '           camp(i,j,nkp)=',camp(i,j,nkp)
!write(stdout,*) '           cz(l,j)=',cz(l,j)
!write(stdout,*) '           u_matrix_opt(i,l,nkp)=',u_matrix_opt(i,l,nkp)

! aam: 20/10/2006 -- the second dimension of u_matrix_opt is out of bounds (allocated as num_wann)!
! commenting this line out.
!                     camp(i,j,nkp) = camp(i,j,nkp) + cz(l,j) * u_matrix_opt(i,l,nkp)
              enddo
            enddo
          enddo
        enddo   ! [loop over k points (nkp)]

      endif   ! [if icompflag=1]

    endif     ![if(index(devel_flag,'compspace')>0)]

    deallocate (history, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating history in dis_extract')

    if (on_root) then
      deallocate (cham, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating cham in dis_extract')
    endif
    if (allocated(camp)) then
      deallocate (camp, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating camp in dis_extract')
    end if
    if (allocated(camp_loc)) then
      deallocate (camp_loc, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating camp_loc in dis_extract')
    endif
    if (on_root) then
      deallocate (ceamp, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating ceamp in dis_extract')
    endif
    deallocate (u_matrix_opt_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating u_matrix_opt_loc in dis_extract')
    deallocate (wkomegai1_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating wkomegai1_loc in dis_extract')
    deallocate (wkomegai1, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating wkomegai1 in dis_extract')

    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cz in dis_extract')
    deallocate (cwork, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cwork in dis_extract')
    deallocate (cap, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cap in dis_extract')
    deallocate (rwork, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating rwork in dis_extract')
    deallocate (w, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating w in dis_extract')
    deallocate (ifail, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating ifail in dis_extract')
    deallocate (iwork, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating iwork in dis_extract')

    deallocate (cbw, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cbw in dis_extract')
    deallocate (cww, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cww in dis_extract')
    deallocate (cwb, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cwb in dis_extract')

    if (on_root) write (stdout, '(1x,a/)') &
      '+----------------------------------------------------------------------------+'

    if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract', 2)

    return

  contains

    subroutine internal_test_convergence()
      !! Check if we have converged

      implicit none

      integer :: ierr
      real(kind=dp), allocatable :: temp_hist(:)

      allocate (temp_hist(dis_conv_window), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating temp_hist in dis_extract')

      if (iter .le. dis_conv_window) then
        history(iter) = delta_womegai
      else
        temp_hist = eoshift(history, 1, delta_womegai)
        history = temp_hist
      endif

      dis_converged = .false.
      if (iter .ge. dis_conv_window) then
        dis_converged = all(abs(history) .lt. dis_conv_tol)
      endif

      deallocate (temp_hist, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating temp_hist in dis_extract')

      return

    end subroutine internal_test_convergence

    !==================================================================!
    subroutine internal_zmatrix(nkp, nkp_loc, cmtrx)
      !==================================================================!
      !! Compute the Z-matrix
      !                                                                  !
      !                                                                  !
      !                                                                  !
      !==================================================================!

      implicit none

      integer, intent(in) :: nkp
      integer, intent(in) :: nkp_loc
      !! Which kpoint
      complex(kind=dp), intent(out) :: cmtrx(num_bands, num_bands)
      !! (M,N)-TH ENTRY IN THE (NDIMWIN(NKP)-NDIMFROZ(NKP)) x (NDIMWIN(NKP)-NDIMFRO
      !! HERMITIAN MATRIX AT THE NKP-TH K-POINT

      ! Internal variables
      integer          :: l, m, n, p, q, nn, nkp2, ndimk
      complex(kind=dp) :: csum

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract: zmatrix', 1)

      cmtrx = cmplx_0
      ndimk = ndimwin(nkp) - ndimfroz(nkp)
      do nn = 1, nntot
        nkp2 = nnlist(nkp, nn)
        call zgemm('N', 'N', num_bands, num_wann, ndimwin(nkp2), cmplx_1, &
                   m_matrix_orig_local(:, :, nn, nkp_loc), num_bands, u_matrix_opt(:, :, nkp2), num_bands, &
                   cmplx_0, cbw, num_bands)
        do n = 1, ndimk
          q = indxnfroz(n, nkp)
          do m = 1, n
            p = indxnfroz(m, nkp)
            csum = cmplx_0
            do l = 1, num_wann
              csum = csum + cbw(p, l)*conjg(cbw(q, l))
            enddo
            cmtrx(m, n) = cmtrx(m, n) + cmplx(wb(nn), 0.0_dp, kind=dp)*csum
            cmtrx(n, m) = conjg(cmtrx(m, n))
          enddo
        enddo
      enddo

      if (timing_level > 1 .and. on_root) call io_stopwatch('dis: extract: zmatrix', 2)

      return

    end subroutine internal_zmatrix

!~      !==================================================================!
!~!      function dis_zeig(nkp,m,cmk)
!~      function dis_zeig(nkp,m)
!~      !==================================================================!
!~      !                                                                  !
!~      !                                                                  !
!~      !                                                                  !
!~      !                                                                  !
!~      !==================================================================!
!~
!~        ! Computes <lambda>_mk = sum_{n=1}^N sum_b w_b |<u_{mk}|u_{n,k+b}>|^2
!~        ! [See Eqs. (12) and (17) of SMV]
!~
!~        implicit none
!~
!~        integer, intent(in) :: nkp
!~        integer, intent(in) :: m
!~!        complex(kind=dp), intent(in) :: cmk(num_bands,num_bands,nntot)
!~
!~        ! Internal variables
!~        real(kind=dp) :: dis_zeig
!~        complex(kind=dp) :: cdot_bloch
!~        integer :: n,nn,ndnnx,ndnn,nnsh,nkp2,l,j
!~
!~        dis_zeig=0.0_dp
!~
!~!        do nn=1,nntot
!~!           nkp2=nnlist(nkp,nn)
!~        do n = 1, num_wann
!~           do nn = 1, nntot
!~                 nkp2 = nnlist(nkp,nn)
!~                 ! Dotproduct
!~                 cdot_bloch = cmplx_0
!~                 do l = 1, ndimwin(nkp)
!~                    do j = 1, ndimwin(nkp2)
!~                       cdot_bloch = cdot_bloch + &
!~!                            conjg(u_matrix_opt(l,m,nkp)) * u_matrix_opt(j,n,nkp2) * cmk(l,j,nn)
!~                            conjg(u_matrix_opt(l,m,nkp)) * u_matrix_opt(j,n,nkp2) * m_matrix_orig(l,j,nn,nkp)
!~                    enddo
!~                 enddo
!~                 write(stdout,'(a,4i5,2f15.10)') 'zeig:',nkp,nn,m,n,cdot_bloch
!~!                 call zgemm('C','N',num_wann,ndimwin(nkp2),ndimwin(nkp),cmplx_1,&
!~!                      u_matrix_opt(:,:,nkp),num_bands,m_matrix_orig(:,:,nn,nkp),num_bands,cmplx_0,&
!~!                      cwb,num_wann)
!~!                 call zgemm('N','N',num_wann,num_wann,ndimwin(nkp2),cmplx_1,&
!~!                      cwb,num_wann,u_matrix_opt(:,:,nkp),num_bands,cmplx_0,cww,num_wann)
!~
!~                 dis_zeig = dis_zeig + wb(nn) * abs(cdot_bloch)**2
!~
!~!                 do n=1,num_wann
!~!                    dis_zeig = dis_zeig + wb(nn) * abs(cww(m,n))**2
!~!                 enddo
!~
!~              enddo
!~        enddo
!~
!~        return
!~
!~      end function dis_zeig

  end subroutine dis_extract

![ysl-b]
  !==================================================================!
  subroutine dis_extract_gamma()
    !==================================================================!
    !                                                                  !
    !! Extracts an num_wann-dimensional subspace at each k by
    !! minimizing Omega_I (Gamma point version)
    !                                                                  !
    !==================================================================!

    use w90_io, only: io_time

    implicit none

    ! MODIFIED:
    !           u_matrix_opt (At input it contains the initial guess for the optimal
    ! subspace (expressed in terms of the original states inside the window). At
    ! output it contains the  states that diagonalize the hamiltonian inside the
    ! optimal subspace (again expressed in terms of the original window states).
    ! Giving out states that diagonalize the hamiltonian inside the optimal
    ! subspace (instead of the eigenstates of the Z matrix) is useful for
    ! performing the Wannier interpolation of the energy bands as described in
    ! Sec. III.F of SMV)
    !
    !           eigval (At input: original energy eigenvalues.
    ! At output: eigenvalues of the hamiltonian inside optimal subspace)

    ! ----------------------------------------------------------------------
    ! TO DO: The complement subspace is computed but is not saved anywhere!
    ! (Check what was done with it in original code space.f)
    ! Diagonalize Z matrix only at those k points where ndimwin>num_wann?
    ! ----------------------------------------------------------------------

    ! *******************
    ! SHELLS OF K-VECTORS
    ! *******************
    ! nshells           number of shells of k-points to be used in the
    !                   finite-difference formulas for the k-derivatives
    ! aam: wb is now wb(1:nntot) 09/04/2006
    ! wb(nkp,nnx)       weight of the nnx-th b-vector (ordered along shells
    !                   of increasing length) associated with the nkp-th k-p
    ! wbtot             sum of the weights of all b-vectors associated with
    !                   given k-point (k-point 1 is used in calculation)
    ! nnlist(nkp,nnx)   vkpt(1:3,nnlist(nkp,nnx)) is the nnx-th neighboring
    !                   k-point of the nkp-th k-point vkpt(1:3,nkp) (or its
    !                   periodic image in the "home Brillouin zone")
    ! cm(n,m,nkp,nnx)   Overlap matrix <u_nk|u_{m,k+b}>

    ! Internal variables
    integer :: i, j, l, m, n, nn, nkp, nkp2, info, ierr, ndimk, p
    integer :: icompflag, iter, ndiff
    real(kind=dp) :: womegai, wkomegai, womegai1, rsum, delta_womegai
    real(kind=dp), allocatable :: wkomegai1(:)
    complex(kind=dp), allocatable :: ceamp(:, :, :)
    complex(kind=dp), allocatable :: camp(:, :, :)
    complex(kind=dp), allocatable :: cham(:, :, :)
!@@@
    real(kind=dp), allocatable :: rzmat_in(:, :, :)
    real(kind=dp), allocatable :: rzmat_out(:, :, :)
!@@@
    integer, allocatable :: iwork(:)
    integer, allocatable :: ifail(:)
!@@@
    real(kind=dp), allocatable :: work(:)
    real(kind=dp), allocatable :: cap_r(:)
    real(kind=dp), allocatable :: rz(:, :)
!@@@
    real(kind=dp), allocatable :: w(:)
    complex(kind=dp), allocatable :: cz(:, :)

    complex(kind=dp), allocatable :: cwb(:, :), cww(:, :), cbw(:, :)

    real(kind=dp), allocatable :: history(:)
    logical                       :: dis_converged

    if (timing_level > 1) call io_stopwatch('dis: extract', 1)

    write (stdout, '(/1x,a)') &
      '                  Extraction of optimally-connected subspace                  '
    write (stdout, '(1x,a)') &
      '                  ------------------------------------------                  '

    allocate (cwb(num_wann, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cwb in dis_extract_gamma')
    allocate (cww(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cww in dis_extract_gamma')
    allocate (cbw(num_bands, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cbw in dis_extract_gamma')
    cwb = cmplx_0; cww = cmplx_0; cbw = cmplx_0

    allocate (iwork(5*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating iwork in dis_extract_gamma')
    allocate (ifail(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating ifail in dis_extract_gamma')
    allocate (w(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating w in dis_extract_gamma')
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cz in dis_extract_gamma')
!@@@
    allocate (work(8*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating work in dis_extract_gamma')
    allocate (cap_r((num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cap_r in dis_extract_gamma')
    allocate (rz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating rz in dis_extract_gamma')
!@@@

    allocate (wkomegai1(num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wkomegai1 in dis_extract_gamma')
!@@@
    allocate (rzmat_in(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating rzmat_in in dis_extract')
    allocate (rzmat_out(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating rzmat_out in dis_extract')
!@@@
    allocate (history(dis_conv_window), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating history in dis_extract_gamma')

    ! ********************************************
    ! ENERGY WINDOWS AND SUBSPACES AT EACH K-POINT
    ! ********************************************
    ! num_wann             dimensionality of the subspace at each k-point
    !                   (number of Wannier functions per unit cell that we w
    ! NDIMWIN(NKP)      number of bands at the nkp-th k-point that fall
    !                   within the outer energy window
    ! NDIMFROZ(NKP)     number of frozen bands at the nkp-th k-point
    ! INDXNFROZ(I,NKP)  INDEX (BETWEEN 1 AND NDIMWIN(NKP)) OF THE I-TH NON-F
    !                   ORIGINAL BAND STATE AT THE NKP-TH K-POINT
    ! U_MATRIX_OPT(J,L,NKP) AMPLITUDE OF THE J-TH ENERGY EIGENVECTOR INSIDE THE
    !                   ENERGY WINDOW AT THE NKP-TH K-POINT IN THE EXPANSION OF
    !                   THE L-TH LEADING RLAMBDA EIGENVECTOR AT THE SAME K-POINT
    !                   If there are M_k frozen states, they occupy the lowest
    !                   entries of the second index of u_matrix_opt, and the leading
    !                   nbands-M_k eigenvectors of the Z matrix occupy the
    !                   remaining slots
    ! CAMP(J,L,NKP)     SAME AS U_MATRIX_OPT, BUT FOR THE COMPLEMENT SUBSPACE INSIDE THE
    !                   ENERGY WINDOW (I.E., THE NON-LEADING RLAMBDA EIGENVECTORS)
    ! CEAMP(J,L,NKPTS)  SAME AS U_MATRIX_OPT, BUT INSTEAD OF RLAMBDA EIGENVECTOR, I
    !                   FOR THE ENERGY EIGENVECTOR OBTAINED BY DIAGONALIZING
    !                   HAMILTONIAN IN THE OPTIMIZED SUBSPACE
    ! RZMAT_IN(M,N,NKP) Z-MATRIX [Eq. (21) SMV]
    ! RZMAT_OUT(M,N,NKP) OUTPUT Z-MATRIX FROM THE PRESENT ITERATION
    ! RLAMBDA           An eigenvalue of the Z matrix
    ! womegai           Gauge-invariant Wannier spread, computed usinf all s
    !                   from current iteration
    ! wkomegai1(NKP)    Eq. (18) of SMV
    ! womegai1          Eq.(11) of SMV (like wowmegai, but neighboring state
    !                   for computing overlaps are from previous iteration.
    !                   become equal at self-consistency)
    ! alphafixe         mixing parameter for the iterative procedure
    ! nitere            total number of iterations

    ! DEBUG
    if (iprint > 2) then
      write (stdout, '(a,/)') '  Original eigenvalues inside outer window:'
      do nkp = 1, num_kpts
        write (stdout, '(a,i3,3x,20(f9.5,1x))') '  K-point ', nkp, &
          (eigval_opt(i, nkp), i=1, ndimwin(nkp))
      enddo
    endif
    ! ENDDEBUG

    ! TO DO: Check if this is the best place to initialize icompflag
    icompflag = 0

    write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'
    write (stdout, '(1x,a)') &
      '|  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS'
    write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'

    dis_converged = .false.

    ! ------------------
    ! BIG ITERATION LOOP
    ! ------------------
    do iter = 1, dis_num_iter

      if (iter .eq. 1) then
        ! Initialize Z matrix at k points w/ non-frozen states
        do nkp = 1, num_kpts
          if (num_wann .gt. ndimfroz(nkp)) call internal_zmatrix_gamma(nkp, rzmat_in(:, :, nkp))
        enddo
      else
        ! [iter.ne.1]
        ! Update Z matrix at k points with non-frozen states, using a mixing sch
        do nkp = 1, num_kpts
          if (num_wann .gt. ndimfroz(nkp)) then
            ndimk = ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, ndimk
              do j = 1, i
                rzmat_in(j, i, nkp) = &
                  dis_mix_ratio*rzmat_out(j, i, nkp) &
                  + (1.0_dp - dis_mix_ratio)*rzmat_in(j, i, nkp)
                ! hermiticity
                rzmat_in(i, j, nkp) = rzmat_in(j, i, nkp)
              enddo
            enddo
          endif
        enddo
      endif
      ! [if iter=1]

      womegai1 = 0.0_dp
      ! wkomegai1 is defined by Eq. (18) of SMV.
      ! Contribution to wkomegai1 from frozen states should be calculated now
      ! every k (before updating any k), so that for iter>1 overlaps are with
      ! non-frozen neighboring states from the previous iteration

      wkomegai1 = real(num_wann, dp)*wbtot
      do nkp = 1, num_kpts
        if (ndimfroz(nkp) .gt. 0) then
          do nn = 1, nntot
            nkp2 = nnlist(nkp, nn)
            call zgemm('C', 'N', ndimfroz(nkp), ndimwin(nkp2), ndimwin(nkp), cmplx_1, &
                       u_matrix_opt(:, :, nkp), num_bands, m_matrix_orig(:, :, nn, nkp), num_bands, cmplx_0, &
                       cwb, num_wann)
            call zgemm('N', 'N', ndimfroz(nkp), num_wann, ndimwin(nkp2), cmplx_1, &
                       cwb, num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
            rsum = 0.0_dp
            do n = 1, num_wann
              do m = 1, ndimfroz(nkp)
                rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
              enddo
            enddo
            wkomegai1(nkp) = wkomegai1(nkp) - wb(nn)*rsum
          enddo
        endif
      enddo

      ! Refine optimal subspace at k points w/ non-frozen states
      do nkp = 1, num_kpts
        if (num_wann .gt. ndimfroz(nkp)) then
          ! Diagonalize Z matrix
          do j = 1, ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, j
              cap_r(i + ((j - 1)*j)/2) = rzmat_in(i, j, nkp)
            enddo
          enddo
          ndiff = ndimwin(nkp) - ndimfroz(nkp)
          call DSPEVX('V', 'A', 'U', ndiff, cap_r, 0.0_dp, 0.0_dp, 0, 0, &
                      -1.0_dp, m, w, rz, num_bands, work, iwork, ifail, info)
          if (info .lt. 0) then
            write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING Z MATRIX'
            write (stdout, *) ' THE ', -info, ' ARGUMENT OF DSPEVX HAD AN ILLEGAL VALUE'
            call io_error(' dis_extract_gamma: error')
          endif
          if (info .gt. 0) then
            write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING Z MATRIX'
            write (stdout, *) info, ' EIGENVECTORS FAILED TO CONVERGE'
            call io_error(' dis_extract_gamma: error')
          endif
          cz(:, :) = cmplx(rz(:, :), 0.0_dp, dp)
          !
          ! Update the optimal subspace by incorporating the num_wann-ndimfroz(nkp) l
          ! eigenvectors of the Z matrix into u_matrix_opt. Also, add contribution from
          ! non-frozen states to wkomegai1(nkp) (minus the corresponding eigenvalu
          m = ndimfroz(nkp)
          do j = ndimwin(nkp) - num_wann + 1, ndimwin(nkp) - ndimfroz(nkp)
            m = m + 1
            wkomegai1(nkp) = wkomegai1(nkp) - w(j)
            u_matrix_opt(1:ndimwin(nkp), m, nkp) = cmplx_0
            ndimk = ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, ndimk
              p = indxnfroz(i, nkp)
              u_matrix_opt(p, m, nkp) = cz(i, j)
            enddo
          enddo
        endif
        ! [if num_wann>ndimfroz(nkp)]

        ! Now that we have contribs. from both frozen and non-frozen states to
        ! wkomegai1(nkp), add it to womegai1
        womegai1 = womegai1 + wkomegai1(nkp)

        ! AT THE LAST ITERATION FIND A BASIS FOR THE (NDIMWIN(NKP)-num_wann)-DIMENS
        ! COMPLEMENT SPACE
        if (index(devel_flag, 'compspace') > 0) then

          if (iter .eq. dis_num_iter) then
            allocate (camp(num_bands, num_bands, num_kpts), stat=ierr)
            camp = cmplx_0
            if (ierr /= 0) call io_error('Error allocating camp in dis_extract_gamma')
            if (ndimwin(nkp) .gt. num_wann) then
              do j = 1, ndimwin(nkp) - num_wann
                if (num_wann .gt. ndimfroz(nkp)) then
                  ! USE THE NON-LEADING EIGENVECTORS OF THE Z-MATRIX
                  camp(1:ndimwin(nkp), j, nkp) = cz(1:ndimwin(nkp), j)
                else
                  ! Then num_wann=NDIMFROZ(NKP)
                  ! USE THE ORIGINAL NON-FROZEN BLOCH EIGENSTATES
                  do i = 1, ndimwin(nkp)
                    camp(i, j, nkp) = cmplx_0
                    if (i .eq. indxnfroz(j, nkp)) camp(i, j, nkp) = cmplx_1
                  enddo
                endif
              enddo
            else
              icompflag = 1
            endif
          endif
        end if

      enddo
      ! [Loop over k points (nkp)]

      womegai1 = womegai1/real(num_kpts, dp)

      ! DEBUG
      ! Orthonormality check
      !         do nkp=1,nkpts
      !           write(*,*) ' '
      !           write(*,'(a8,i4)') 'k-point ',nkp
      !           do l=1,num_wann
      !           do m=1,l
      !             ctmp=czero
      !             do j=1,ndimwin(nkp)
      !               ctmp=ctmp+conjg(u_matrix_opt(j,m,nkp))*u_matrix_opt(j,l,nkp)
      !             enddo
      !             write(*,'(i2,2x,i2,f16.12,1x,f16.12)') l,m,ctmp
      !             if(l.eq.m) then
      !               if(abs(ctmp-cmplx(1.0d0,0.0d0)).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             else
      !               if(abs(ctmp).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             endif
      !           enddo
      !           enddo
      !         enddo
      ! ENDDEBUG

      ! Compute womegai  using the updated subspaces at all k, i.e.,
      ! replacing (i-1) by (i) in Eq. (12) SMV

      womegai = 0.0_dp
      do nkp = 1, num_kpts
        wkomegai = 0.0_dp
        do nn = 1, nntot
          nkp2 = nnlist(nkp, nn)
          call zgemm('C', 'N', num_wann, ndimwin(nkp2), ndimwin(nkp), cmplx_1, &
                     u_matrix_opt(:, :, nkp), num_bands, m_matrix_orig(:, :, nn, nkp), num_bands, cmplx_0, &
                     cwb, num_wann)
          call zgemm('N', 'N', num_wann, num_wann, ndimwin(nkp2), cmplx_1, &
                     cwb, num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
          rsum = 0.0_dp
          do n = 1, num_wann
            do m = 1, num_wann
              rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
            enddo
          enddo
          wkomegai = wkomegai + wb(nn)*rsum
        enddo
        wkomegai = real(num_wann, dp)*wbtot - wkomegai
        womegai = womegai + wkomegai
      enddo
      womegai = womegai/real(num_kpts, dp)
      ! [Loop over k (nkp)]

      delta_womegai = womegai1/womegai - 1.0_dp

      write (stdout, 124) iter, womegai1*lenconfac**2, womegai*lenconfac**2, &
        delta_womegai, io_time()

124   format(2x, i6, 3x, f14.8, 3x, f14.8, 6x, es10.3, 2x, f8.2, 4x, '<-- DIS')

      ! Construct the updated Z matrix, CZMAT_OUT, at k points w/ non-frozen s
      do nkp = 1, num_kpts
        if (num_wann .gt. ndimfroz(nkp)) call internal_zmatrix_gamma(nkp, rzmat_out(:, :, nkp))
      enddo

      call internal_test_convergence()

      if (dis_converged) then
        write (stdout, '(/13x,a,es10.3,a,i2,a)') &
          '<<<      Delta <', dis_conv_tol, &
          '  over ', dis_conv_window, ' iterations     >>>'
        write (stdout, '(13x,a)') '<<< Disentanglement convergence criteria satisfied >>>'
        exit
      endif

    enddo
    ! [BIG ITERATION LOOP (iter)]

    deallocate (rzmat_out, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating rzmat_out in dis_extract_gamma')
    deallocate (rzmat_in, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating rzmat_in in dis_extract_gamma')

    allocate (ceamp(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating ceamp in dis_extract_gamma')
    allocate (cham(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating cham in dis_extract_gamma')

    if (.not. dis_converged) then
      write (stdout, '(/5x,a)') &
        '<<< Warning: Maximum number of disentanglement iterations reached >>>'
      write (stdout, '(10x,a)') '<<< Disentanglement convergence criteria not satisfied >>>'
    endif

    if (icompflag .eq. 1) then
      if (iprint > 2) then
        write (stdout, ('(/4x,a)')) &
          'WARNING: Complement subspace has zero dimensions at the following k-points:'
        i = 0
        write (stdout, '(4x)', advance='no')
        do nkp = 1, num_kpts
          if (ndimwin(nkp) .eq. num_wann) then
            i = i + 1
            if (i .le. 12) then
              write (stdout, '(i6)', advance='no') nkp
            else
              i = 1
              write (stdout, '(/4x)', advance='no')
              write (stdout, '(i6)', advance='no') nkp
            endif
          endif
        enddo
      endif
    endif

    ! Write the final womegai. This should remain unchanged during the
    ! subsequent minimization of Omega_tilde in wannierise.f90
    ! We store it in the checkpoint file as a sanity check
    write (stdout, '(/8x,a,f14.8,a/)') 'Final Omega_I ', &
      womegai*lenconfac**2, ' ('//trim(length_unit)//'^2)'

    ! Set public variable omega_invariant
    omega_invariant = womegai

    ! DIAGONALIZE THE HAMILTONIAN WITHIN THE OPTIMIZED SUBSPACES
    do nkp = 1, num_kpts

      do j = 1, num_wann
        do i = 1, num_wann
          cham(i, j, nkp) = cmplx_0
          do l = 1, ndimwin(nkp)
            cham(i, j, nkp) = cham(i, j, nkp) + conjg(u_matrix_opt(l, i, nkp)) &
                              *u_matrix_opt(l, j, nkp)*eigval_opt(l, nkp)
          enddo
        enddo
      enddo
!@@@
      do j = 1, num_wann
        do i = 1, j
          cap_r(i + ((j - 1)*j)/2) = real(cham(i, j, nkp), dp)
        enddo
      enddo

      call DSPEVX('V', 'A', 'U', num_wann, cap_r, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                  m, w, rz, num_bands, work, iwork, ifail, info)

      if (info .lt. 0) then
        write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
        write (stdout, *) ' THE ', -info, ' ARGUMENT OF DSPEVX HAD AN ILLEGAL VALUE'
        call io_error(' dis_extract_gamma: error')
      endif
      if (info .gt. 0) then
        write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
        write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
        call io_error(' dis_extract_gamma: error')
      endif

      cz = cmplx_0
      cz(1:num_wann, 1:num_wann) = cmplx(rz(1:num_wann, 1:num_wann), 0.0_dp, dp)

!@@@
      ! Store the energy eigenvalues of the optimal subspace (used in wann_ban
      eigval_opt(1:num_wann, nkp) = w(1:num_wann)

      ! CALCULATE AMPLITUDES OF THE CORRESPONDING ENERGY EIGENVECTORS IN TERMS
      ! THE ORIGINAL ("WINDOW SPACE") ENERGY EIGENVECTORS
      do j = 1, num_wann
        do i = 1, ndimwin(nkp)
          ceamp(i, j, nkp) = cmplx_0
          do l = 1, num_wann
            ceamp(i, j, nkp) = ceamp(i, j, nkp) + cz(l, j)*u_matrix_opt(i, l, nkp)
          enddo
        enddo
      enddo
      ! NKP
    enddo
    ! DEBUG
    if (iprint > 2) then
      write (stdout, '(/,a,/)') '  Eigenvalues inside optimal subspace:'
      do nkp = 1, num_kpts
        write (stdout, '(a,i3,2x,20(f9.5,1x))') '  K-point ', &
          nkp, (eigval_opt(i, nkp), i=1, num_wann)
      enddo
    endif
    ! ENDDEBUG

    ! Replace u_matrix_opt by ceamp. Both span the
    ! same space, but the latter is more convenient for the purpose of obtai
    ! an optimal Fourier-interpolated band structure: see Sec. III.E of SMV.
    do nkp = 1, num_kpts
      do j = 1, num_wann
        u_matrix_opt(1:ndimwin(nkp), j, nkp) = ceamp(1:ndimwin(nkp), j, nkp)
      enddo
    enddo

    ! aam: 01/05/2009: added devel_flag if statement as the complementary
    !      subspace code was causing catastrophic seg-faults
    if (index(devel_flag, 'compspace') > 0) then

      ! The compliment subspace code needs work: jry
      if (icompflag .eq. 1) then
        if (iprint > 2) then
          write (stdout, *) 'AT SOME K-POINT(S) COMPLEMENT SUBSPACE HAS ZERO DIMENSIONALITY'
          write (stdout, *) '=> DID NOT CREATE FILE COMPSPACE.DAT'
        endif
      else
        ! DIAGONALIZE THE HAMILTONIAN IN THE COMPLEMENT SUBSPACE, WRITE THE
        ! CORRESPONDING EIGENFUNCTIONS AND ENERGY EIGENVALUES
        do nkp = 1, num_kpts
          do j = 1, ndimwin(nkp) - num_wann
            do i = 1, ndimwin(nkp) - num_wann
              cham(i, j, nkp) = cmplx_0
              do l = 1, ndimwin(nkp)
                cham(i, j, nkp) = cham(i, j, nkp) + conjg(camp(l, i, nkp)) &
                                  *camp(l, j, nkp)*eigval_opt(l, nkp)
              enddo
            enddo
          enddo
!@@@
          do j = 1, ndimwin(nkp) - num_wann
            do i = 1, j
              cap_r(i + ((j - 1)*j)/2) = real(cham(i, j, nkp), dp)
            enddo
          enddo
          ndiff = ndimwin(nkp) - num_wann

          call DSPEVX('V', 'A', 'U', ndiff, cap_r, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                      m, w, rz, num_bands, work, iwork, ifail, info)

          if (info .lt. 0) then
            write (stdout, *) '*** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
            write (stdout, *) 'THE ', -info, ' ARGUMENT OF DSPEVX HAD AN ILLEGAL VALUE'
            call io_error(' dis_extract_gamma: error')
          endif
          if (info .gt. 0) then
            write (stdout, *) '*** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
            write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
            call io_error(' dis_extract_gamma: error')
          endif

          cz = cmplx_0
          cz(1:ndiff, 1:ndiff) = cmplx(rz(1:ndiff, 1:ndiff), 0.0_dp, dp)

!@@@
          ! CALCULATE AMPLITUDES OF THE ENERGY EIGENVECTORS IN THE COMPLEMENT SUBS
          ! TERMS OF THE ORIGINAL ENERGY EIGENVECTORS
          do j = 1, ndimwin(nkp) - num_wann
            do i = 1, ndimwin(nkp)
              camp(i, j, nkp) = cmplx_0
              do l = 1, ndimwin(nkp) - num_wann
!write(stdout,*) 'i=',i,'   j=',j,'   l=',l
!write(stdout,*) '           camp(i,j,nkp)=',camp(i,j,nkp)
!write(stdout,*) '           cz(l,j)=',cz(l,j)
!write(stdout,*) '           u_matrix_opt(i,l,nkp)=',u_matrix_opt(i,l,nkp)

! aam: 20/10/2006 -- the second dimension of u_matrix_opt is out of bounds (allocated as num_wann)!
! commenting this line out.
!                     camp(i,j,nkp) = camp(i,j,nkp) + cz(l,j) * u_matrix_opt(i,l,nkp)
              enddo
            enddo
          enddo
        enddo
        ! [loop over k points (nkp)]

      endif
      ! [if icompflag=1]

    endif
    ! [if index(devel_flag,'compspace')>0]

    deallocate (history, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating history in dis_extract_gamma')

    deallocate (cham, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cham in dis_extract_gamma')
    if (allocated(camp)) then
      deallocate (camp, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating camp in dis_extract_gamma')
    end if
    deallocate (ceamp, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating ceamp in dis_extract_gamma')
    deallocate (wkomegai1, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating wkomegai1 in dis_extract_gamma')

!@@@
    deallocate (rz, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating rz in dis_extract_gamma')
    deallocate (cap_r, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cap_r in dis_extract_gamma')
    deallocate (work, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating work in dis_extract_gamma')
!@@@
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cz in dis_extract_gamma')
    deallocate (w, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating w in dis_extract_gamma')
    deallocate (ifail, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating ifail in dis_extract_gamma')
    deallocate (iwork, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating iwork in dis_extract_gamma')

    deallocate (cbw, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cbw in dis_extract_gamma')
    deallocate (cww, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cww in dis_extract_gamma')
    deallocate (cwb, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating cwb in dis_extract_gamma')

    write (stdout, '(1x,a/)') &
      '+----------------------------------------------------------------------------+'

    if (timing_level > 1) call io_stopwatch('dis: extract_gamma', 2)

    return

  contains

    subroutine internal_test_convergence()
      !! Test for convergence (Gamma point routine)

      implicit none

      integer :: ierr
      real(kind=dp), allocatable :: temp_hist(:)

      allocate (temp_hist(dis_conv_window), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating temp_hist in dis_extract_gamma')

      if (iter .le. dis_conv_window) then
        history(iter) = delta_womegai
      else
        temp_hist = eoshift(history, 1, delta_womegai)
        history = temp_hist
      endif

      dis_converged = .false.
      if (iter .ge. dis_conv_window) then
        dis_converged = all(abs(history) .lt. dis_conv_tol)
      endif

      deallocate (temp_hist, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating temp_hist in dis_extract_gamma')

      return

    end subroutine internal_test_convergence

    !==================================================================!
    subroutine internal_zmatrix_gamma(nkp, rmtrx)
      !==================================================================!
      !! Compute Z-matrix (Gamma point routine)
      !                                                                  !
      !                                                                  !
      !                                                                  !
      !==================================================================!

      implicit none

      integer, intent(in) :: nkp
      !! Which k-point
      real(kind=dp), intent(out) :: rmtrx(num_bands, num_bands)
      !!(M,N)-TH ENTRY IN THE (NDIMWIN(NKP)-NDIMFROZ(NKP)) x (NDIMWIN(NKP)-NDIMFRO
      !! HERMITIAN MATRIX AT THE NKP-TH K-POINT

      ! Internal variables
      integer          :: l, m, n, p, q, nn, nkp2, ndimk
      complex(kind=dp) :: csum

      if (timing_level > 1) call io_stopwatch('dis: extract_gamma: zmatrix_gamma', 1)

      rmtrx = 0.0_dp
      ndimk = ndimwin(nkp) - ndimfroz(nkp)
      do nn = 1, nntot
        nkp2 = nnlist(nkp, nn)
        call zgemm('N', 'N', num_bands, num_wann, ndimwin(nkp2), cmplx_1, &
                   m_matrix_orig(:, :, nn, nkp), num_bands, u_matrix_opt(:, :, nkp2), num_bands, &
                   cmplx_0, cbw, num_bands)
        do n = 1, ndimk
          q = indxnfroz(n, nkp)
          do m = 1, n
            p = indxnfroz(m, nkp)
            csum = cmplx_0
            do l = 1, num_wann
              csum = csum + cbw(p, l)*conjg(cbw(q, l))
            enddo
            rmtrx(m, n) = rmtrx(m, n) + wb(nn)*real(csum, dp)
            rmtrx(n, m) = rmtrx(m, n)
          enddo
        enddo
      enddo

      if (timing_level > 1) call io_stopwatch('dis: extract_gamma: zmatrix_gamma', 2)

      return

    end subroutine internal_zmatrix_gamma

  end subroutine dis_extract_gamma

![ysl-e]

end module w90_disentangle
