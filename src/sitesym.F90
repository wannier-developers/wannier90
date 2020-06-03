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
!                                                            !
! Reference:                                                 !
!    R. Sakuma, Symmetry-adapted Wannier functions           !
!    in the maximal localization procedure,                  !
!    Phys Rev B 87, 235109 (2013)                            !
!                                                            !
!------------------------------------------------------------!

module w90_sitesym
  !! Routines to impose the site symmetry during minimisation of spread

  use w90_constants, only: dp, cmplx_1, cmplx_0
  use w90_io, only: io_error, stdout

  implicit none

  private :: symmetrize_ukirr
  private :: orthogonalize_u

  public  :: sitesym_slim_d_matrix_band
  public  :: sitesym_replace_d_matrix_band
  public  :: sitesym_symmetrize_u_matrix
  public  :: sitesym_symmetrize_gradient
  public  :: sitesym_symmetrize_rotation
  public  :: sitesym_symmetrize_zmatrix
  public  :: sitesym_dis_extract_symmetry
  public  :: sitesym_read
  public  :: sitesym_dealloc

  ! Variables and parameters needed by other modules
  integer, public, save :: nkptirr = 9999
  integer, public, save :: nsymmetry = 9999
  integer, allocatable, public, save :: kptsym(:, :), ir2ik(:), ik2ir(:)
  complex(kind=dp), allocatable, public, save :: d_matrix_band(:, :, :, :)
  complex(kind=dp), allocatable, public, save :: d_matrix_wann(:, :, :, :)

contains

  !==================================================================!
  subroutine sitesym_slim_d_matrix_band(lwindow_in)
    !==================================================================!
    use w90_parameters, only: num_bands, num_kpts

    implicit none

    logical, optional, intent(in) :: lwindow_in(num_bands, num_kpts)
    integer :: ik, i, j, nb, ir
    integer :: nindx(num_bands)

    !write(stdout,"(a)") '-- sitesym_slim_d_matrix_band --'

    do ir = 1, nkptirr
      ik = ir2ik(ir)
      j = 0
      do i = 1, num_bands
        if (lwindow_in(i, ik)) then
          j = j + 1
          nindx(j) = i
        endif
      enddo
      nb = j
      do j = 1, nb
        i = nindx(j)
        d_matrix_band(1:nb, j, :, ir) = d_matrix_band(nindx(1:nb), i, :, ir)
        if (nb .lt. num_bands) then
          d_matrix_band(nb + 1:, j, :, ir) = 0
        endif
      enddo
    enddo

    return
  end subroutine sitesym_slim_d_matrix_band

  !==================================================================!
  subroutine sitesym_replace_d_matrix_band()
    !==================================================================!
    use w90_parameters, only: num_wann

    implicit none

    !write(stdout,"(a)") '-- sitesym_replace_d_matrix_band --'
    !write(stdout,"(a)") 'd_matrix_band is replaced by d_matrix_wann'
    deallocate (d_matrix_band)
    allocate (d_matrix_band(num_wann, num_wann, nsymmetry, nkptirr))
    d_matrix_band = d_matrix_wann

    return
  end subroutine sitesym_replace_d_matrix_band

  !==========================================================================!
  subroutine sitesym_symmetrize_u_matrix(ndim, umat, lwindow_in)
    !==========================================================================!
    !                                                                          !
    ! calculate U(Rk)=d(R,k)*U(k)*D^{\dagger}(R,k) in the following two cases: !
    !                                                                          !
    ! 1. "disentanglement" phase (present(lwindow))                            !
    !    ndim=num_bands                                                        !
    !                                                                          !
    ! 2. Minimization of Omega_{D+OD} (.not.present(lwindow))                  !
    !    ndim=num_wann,  d=d_matrix_band                                       !
    !                                                                          !
    !==========================================================================!
    use w90_parameters, only: num_wann, num_bands, num_kpts

    implicit none

    integer, intent(in) :: ndim
    complex(kind=dp), intent(inout) :: umat(ndim, num_wann, num_kpts)
    logical, optional, intent(in) :: lwindow_in(num_bands, num_kpts)

    ! local
    integer :: ik, ir, isym, irk, n
    logical :: ldone(num_kpts)
    complex(kind=dp) :: cmat(ndim, num_wann)

    if (present(lwindow_in) .and. (ndim .ne. num_bands)) call io_error('ndim!=num_bands')
    if (.not. present(lwindow_in)) then
      if (ndim .ne. num_wann) call io_error('ndim!=num_wann')
    endif

    ldone = .false.
    do ir = 1, nkptirr
      ik = ir2ik(ir)
      ldone(ik) = .true.
      if (present(lwindow_in)) then
        n = count(lwindow_in(:, ik))
      else
        n = ndim
      endif
      if (present(lwindow_in)) then
        call symmetrize_ukirr(ir, ndim, umat(:, :, ik), n)
      else
        call symmetrize_ukirr(ir, ndim, umat(:, :, ik))
      endif
      do isym = 2, nsymmetry
        irk = kptsym(isym, ir)
        if (ldone(irk)) cycle
        ldone(irk) = .true.
        ! cmat = d(R,k) * U(k)
        call zgemm('N', 'N', n, num_wann, n, cmplx_1, &
                   d_matrix_band(:, :, isym, ir), ndim, &
                   umat(:, :, ik), ndim, cmplx_0, cmat, ndim)

        ! umat(Rk) = cmat*D^{+}(R,k) = d(R,k) * U(k) * D^{+}(R,k)
        call zgemm('N', 'C', n, num_wann, num_wann, cmplx_1, cmat, ndim, &
                   d_matrix_wann(:, :, isym, ir), num_wann, cmplx_0, umat(:, :, irk), ndim)
      enddo
    enddo
    if (any(.not. ldone)) call io_error('error in sitesym_symmetrize_u_matrix')

    return
  end subroutine sitesym_symmetrize_u_matrix

  !==================================================================!
  subroutine sitesym_symmetrize_gradient(imode, grad)
    !==================================================================!
    use w90_parameters, only: num_wann, num_kpts
    use w90_utility, only: utility_zgemm

    implicit none

    integer, intent(in) :: imode
    complex(kind=dp), intent(inout) :: grad(num_wann, num_wann, num_kpts)
    integer :: ik, ir, isym, irk, ngk

    complex(kind=dp) :: grad_total(num_wann, num_wann)
    complex(kind=dp) :: cmat1(num_wann, num_wann)
    complex(kind=dp) :: cmat2(num_wann, num_wann)

    logical :: lfound(num_kpts)

    if (imode .eq. 1) then
      lfound = .false.
      do ir = 1, nkptirr
        ik = ir2ik(ir)
        grad_total = grad(:, :, ik)
        lfound(ik) = .true.
        do isym = 2, nsymmetry
          irk = kptsym(isym, ir)
          if (lfound(irk)) cycle
          lfound(irk) = .true.
          !
          ! cmat1 = D(R,k)^{+} G(Rk) D(R,k)
          ! cmat2 = D(R,k)^{\dagger} G(Rk)
          !
          call utility_zgemm(cmat2, d_matrix_wann(:, :, isym, ir), 'C', &
                             grad(:, :, irk), 'N', num_wann)
          call utility_zgemm(cmat1, cmat2, 'N', &
                             d_matrix_wann(:, :, isym, ir), 'N', num_wann)
          grad_total = grad_total + cmat1
        enddo
        grad(:, :, ik) = grad_total
      enddo
      do ik = 1, num_kpts
        if (ir2ik(ik2ir(ik)) .ne. ik) grad(:, :, ik) = 0
      enddo
    endif ! if (imode.eq.1)
    !
    ! grad -> 1/N_{R'} \sum_{R'} D^{+}(R',k) grad D(R',k)
    ! where R' k = k
    !
    do ir = 1, nkptirr
      ik = ir2ik(ir)
      ngk = count(kptsym(:, ir) .eq. ik)
      if (ngk .eq. 1) cycle
      grad_total = grad(:, :, ik)
      do isym = 2, nsymmetry
        if (kptsym(isym, ir) .ne. ik) cycle
        !
        ! calculate cmat1 = D^{+}(R,k) G(Rk) D(R,k)
        !
        ! step 1: cmat2 =  G(Rk) D(R,k)
        call utility_zgemm(cmat2, grad(:, :, ik), 'N', &
                           d_matrix_wann(:, :, isym, ir), 'N', num_wann)
        ! step 2: cmat1 = D^{+}(R,k) * cmat2
        call utility_zgemm(cmat1, d_matrix_wann(:, :, isym, ir), 'C', &
                           cmat2, 'N', num_wann)
        grad_total = grad_total + cmat1
      enddo
      grad(:, :, ik) = grad_total/ngk
    enddo

    return
  end subroutine sitesym_symmetrize_gradient

  !==================================================================!
  subroutine sitesym_symmetrize_rotation(urot)
    !==================================================================!
    use w90_parameters, only: num_wann, num_kpts, u_matrix
    use w90_utility, only: utility_zgemm

    implicit none

    integer :: ik, ir, isym, irk
    complex(kind=dp), intent(inout) :: urot(num_wann, num_wann, num_kpts)

    complex(kind=dp) :: cmat1(num_wann, num_wann)
    complex(kind=dp) :: cmat2(num_wann, num_wann)

    logical :: ldone(num_kpts)

    ldone = .false.
    do ir = 1, nkptirr
      ik = ir2ik(ir)
      ldone(ik) = .true.
      do isym = 2, nsymmetry
        irk = kptsym(isym, ir)
        if (irk .eq. ik) cycle
        if (ldone(irk)) cycle
        ldone(irk) = .true.
        ! cmat2 = UROT(k)*D(R,k)^{\dagger}
        call utility_zgemm(cmat2, urot(:, :, ik), 'N', &
                           d_matrix_wann(:, :, isym, ir), 'C', num_wann)
        ! cmat1 = D(R,k)*cmat2
        call utility_zgemm(cmat1, d_matrix_wann(:, :, isym, ir), 'N', &
                           cmat2, 'N', num_wann)
        urot(:, :, irk) = cmat1(:, :)
      enddo
    enddo
    if (any(.not. ldone)) call io_error('error in sitesym_symmetrize_rotation')

    return
  end subroutine sitesym_symmetrize_rotation

  !==================================================================!
  subroutine sitesym_symmetrize_zmatrix(czmat, lwindow_in)
    !==================================================================!
    !                                                                  !
    !    Z(k) <- \sum_{R} d^{+}(R,k) Z(Rk) d(R,k)                      !
    !                                                                  !
    !==================================================================!
    use w90_parameters, only: num_bands, num_kpts

    implicit none

    complex(kind=dp), intent(inout) :: czmat(num_bands, num_bands, num_kpts)
    logical, intent(in) :: lwindow_in(num_bands, num_kpts)
    logical :: lfound(num_kpts)
    integer :: ik, ir, isym, irk, nd
    complex(kind=dp) :: cztmp(num_bands, num_bands)
    complex(kind=dp) :: cmat1(num_bands, num_bands)
    complex(kind=dp) :: cmat2(num_bands, num_bands)

    lfound = .false.
    do ir = 1, nkptirr
      ik = ir2ik(ir)
      nd = count(lwindow_in(:, ik))
      lfound(ik) = .true.
      do isym = 2, nsymmetry
        irk = kptsym(isym, ir)
        if (lfound(irk)) cycle
        lfound(irk) = .true.
        ! cmat1 = Z(R,k)*d(R,k)
        call zgemm('N', 'N', nd, nd, nd, cmplx_1, czmat(:, :, irk), num_bands, &
                   d_matrix_band(:, :, isym, ir), num_bands, cmplx_0, cmat1, num_bands)
        ! cmat2 = d^{+}(R,k) Z(R,k) d(R,k) = d^{+}(R,k) cmat1
        call zgemm('C', 'N', nd, nd, nd, cmplx_1, d_matrix_band(:, :, isym, ir), num_bands, &
                   cmat1, num_bands, cmplx_0, cmat2, num_bands)
        czmat(:, :, ik) = czmat(:, :, ik) + cmat2(:, :)
      enddo

      cztmp(:, :) = czmat(:, :, ik)
      do isym = 2, nsymmetry
        irk = kptsym(isym, ir)
        if (irk .ne. ik) cycle
        call zgemm('N', 'N', nd, nd, nd, cmplx_1, cztmp, num_bands, &
                   d_matrix_band(:, :, isym, ir), num_bands, cmplx_0, cmat1, num_bands)
        ! cmat2 = d^{+}(R,k) Z(R,k) d(R,k) = d^{+}(R,k) cmat1
        call zgemm('C', 'N', nd, nd, nd, cmplx_1, d_matrix_band(:, :, isym, ir), num_bands, &
                   cmat1, num_bands, cmplx_0, cmat2, num_bands)
        czmat(:, :, ik) = czmat(:, :, ik) + cmat2(:, :)
      enddo
      czmat(:, :, ik) = czmat(:, :, ik)/count(kptsym(:, ir) .eq. ik)
    enddo

    return
  end subroutine sitesym_symmetrize_zmatrix

  !==================================================================!
  subroutine symmetrize_ukirr(ir, ndim, umat, n)
    !==================================================================!
    !                                                                  !
    !  calculate u~(k)=1/N_{R'} \sum_{R'} d^{+}(R',k) u(k) D(R',k)     !
    !  where R'k = k                                                   !
    !  and orthonormalize it                                           !
    !                                                                  !
    !==================================================================!
    use w90_parameters, only: num_wann, num_bands, symmetrize_eps

    implicit none

    integer, intent(in) :: ir, ndim
    complex(kind=dp), intent(inout) :: umat(ndim, num_wann)
    integer, optional, intent(in) :: n

    integer :: isym, ngk, i, iter, ntmp
    integer, parameter :: niter = 100

    real(kind=dp)    :: diff
    complex(kind=dp) :: usum(ndim, num_wann)
    complex(kind=dp) :: cmat_sub(ndim, num_wann)
    complex(kind=dp) :: cmat(ndim, num_wann)
    complex(kind=dp) :: cmat2(num_wann, num_wann)

    !write(stdout,"(a)") '-- symmetrize_ukirr --'
    if (present(n)) then
      if (ndim .ne. num_bands) call io_error('ndim!=num_bands')
      ntmp = n
    else
      if (ndim .ne. num_wann) call io_error('ndim!=num_wann')
      ntmp = ndim
    endif

    ngk = count(kptsym(:, ir) .eq. ir2ik(ir))
    if (ngk .eq. 1) then
      call orthogonalize_u(ndim, num_wann, umat, ntmp)
      return
    endif

    do iter = 1, niter
      usum(:, :) = 0
      cmat2(:, :) = 0
      do i = 1, num_wann
        cmat2(i, i) = cmat2(i, i) + ngk
      enddo
      do isym = 1, nsymmetry
        if (kptsym(isym, ir) .ne. ir2ik(ir)) cycle
        !
        ! cmat = d^{+}(R,k) U(k) D(R,k)
        ! size of umat: umat(ndim,num_wann)
        !
        ! cmat_sub = U(k) D(R,k)
        call zgemm('N', 'N', ntmp, num_wann, num_wann, cmplx_1, &
                   umat, ndim, d_matrix_wann(:, :, isym, ir), num_wann, &
                   cmplx_0, cmat_sub, ndim)
        ! cmat = d^{+}(R,k) * cmat_sub
        call zgemm('C', 'N', ntmp, num_wann, ntmp, cmplx_1, &
                   d_matrix_band(:, :, isym, ir), ndim, cmat_sub, ndim, &
                   cmplx_0, cmat, ndim)
        usum(:, :) = usum(:, :) + cmat(:, :)
        ! check
        cmat2(:, :) = cmat2(:, :) - &
                      matmul(conjg(transpose(umat(:ntmp, :))), cmat(:ntmp, :))
      enddo ! isym
      diff = sum(abs(cmat2))
      if (diff .lt. symmetrize_eps) exit
      if (iter .eq. niter) then
        write (stdout, "(a)") 'Error in symmetrize_u: not converged'
        write (stdout, "(a)") 'Either eps is too small or specified irreps is not'
        write (stdout, "(a)") '  compatible with the bands'
        write (stdout, "(a,2e20.10)") 'diff,eps=', diff, symmetrize_eps
        call io_error('symmetrize_ukirr: not converged')
      endif
      usum = usum/ngk
      call orthogonalize_u(ndim, num_wann, usum, ntmp)
      umat(:, :) = usum
    enddo ! iter

    return
  end subroutine symmetrize_ukirr

  !==================================================================!
  subroutine orthogonalize_u(ndim, m, u, n)
    !==================================================================!

    implicit none

    integer, intent(in) :: ndim, m
    complex(kind=dp), intent(inout) :: u(ndim, m)
    integer, intent(in) :: n

    complex(kind=dp), allocatable :: smat(:, :), evecl(:, :), evecr(:, :)
    complex(kind=dp), allocatable :: WORK(:)
    real(kind=dp), allocatable :: eig(:), RWORK(:)
    integer :: INFO, i, j, l
    integer :: LWORK

    if (n .lt. m) call io_error('n<m')
    allocate (smat(n, m)); smat(1:n, 1:m) = u(1:n, 1:m)
    allocate (evecl(n, n), evecr(m, m))
    allocate (eig(min(m, n)))

    allocate (RWORK(5*min(n, m)))
    LWORK = 2*min(m, n) + max(m, n)
    allocate (WORK(LWORK))

    ! Singular-value decomposition
    call ZGESVD('A', 'A', n, m, smat, n, &
                eig, evecl, n, evecr, m, WORK, LWORK, RWORK, INFO)
    if (info .ne. 0) then
      call io_error(' ERROR: IN ZGESVD IN orthogonalize_u')
    endif
    deallocate (smat, eig, WORK, RWORK)
    ! u_matrix is the initial guess for the unitary rotation of the
    ! basis states given by the subroutine extract
    u = 0
    do j = 1, m
    do l = 1, m
      do i = 1, n
        u(i, j) = u(i, j) + evecl(i, l)*evecr(l, j)
      enddo
    enddo
    enddo
    deallocate (evecl, evecr)

    return
  end subroutine orthogonalize_u

  !==================================================================!
  subroutine sitesym_dis_extract_symmetry(ik, n, zmat, lambda, umat)
    !==================================================================!
    !                                                                  !
    !   minimize Omega_I by steepest descendent                        !
    !                                                                  !
    !   delta U_{mu I}(k) = Z_{mu mu'}*U_{mu' I}(k)                    !
    !                        - \sum_{J} lambda_{JI} U_{mu J}(k)        !
    !   lambda_{JI}=U^{*}_{mu J} Z_{mu mu'} U_{mu' I}                  !
    !                                                                  !
    !==================================================================!
    use w90_parameters, only: num_bands, num_wann

    implicit none

    integer, intent(in) :: ik, n
    complex(kind=dp), intent(in) :: zmat(num_bands, num_bands)
    complex(kind=dp), intent(out) :: lambda(num_wann, num_wann)
    complex(kind=dp), intent(inout) :: umat(num_bands, num_wann)

    complex(kind=dp) :: umatnew(num_bands, num_wann)
    complex(kind=dp) :: ZU(num_bands, num_wann)
    complex(kind=dp) :: deltaU(num_bands, num_wann), carr(num_bands)
    integer :: i, m, INFO, IFAIL(2), IWORK(5*2)
    complex(kind=dp) :: HP(3), SP(3), V(2, 2), CWORK(2*2)
    real(kind=dp)    :: W(2), RWORK(7*2), sp3
    integer :: iter
    integer, parameter :: niter = 50

    do iter = 1, niter
      !  Z*U
      call zgemm('N', 'N', n, num_wann, n, cmplx_1, &
                 zmat, num_bands, umat, num_bands, &
                 cmplx_0, ZU, num_bands)
      ! lambda = U^{+}*Z*U
      call zgemm('C', 'N', num_wann, num_wann, n, cmplx_1, &
                 umat, num_bands, ZU, num_bands, &
                 cmplx_0, lambda, num_wann)

      deltaU(:, :) = ZU(:, :) - matmul(umat, lambda)
      if (sum(abs(deltaU(:n, :))) .lt. 1e-10) return
      ! band-by-band minimization
      do i = 1, num_wann
        ! diagonalize 2x2 matrix
        HP(1) = real(dot_product(umat(1:n, i), ZU(1:n, i)), kind=dp)
        HP(2) = dot_product(ZU(1:n, i), deltaU(1:n, i)) ! (1,2) matrix element
        carr(1:n) = matmul(zmat(1:n, 1:n), deltaU(1:n, i))
        HP(3) = real(dot_product(deltaU(1:n, i), carr(1:n)), kind=dp) ! (2,2)

        SP(1) = real(dot_product(umat(1:n, i), umat(1:n, i)), kind=dp)
        SP(2) = dot_product(umat(1:n, i), deltaU(1:n, i))
        SP(3) = real(dot_product(deltaU(1:n, i), deltaU(1:n, i)), kind=dp)

        sp3 = real(SP(3), kind=dp)
        if (abs(sp3) .lt. 1e-10) then
          umatnew(:, i) = umat(:, i)
          cycle
        endif
        call ZHPGVX(1, 'V', 'A', 'U', 2, HP, SP, 0.0_dp, 0.0_dp, 0, 0, &
                    -1.0_dp, m, W, V, 2, CWORK, RWORK, IWORK, IFAIL, INFO)
        if (INFO .ne. 0) then
          write (stdout, *) 'error in sitesym_dis_extract_symmetry: INFO=', INFO
          if (INFO .gt. 0) then
            if (INFO .le. 2) then
              write (stdout, *) INFO, ' eigenvectors failed to converge'
              write (stdout, *) IFAIL(1:INFO)
            else
              write (stdout, *) ' S is not positive definite'
              write (stdout, *) 'sp3=', sp3
            endif
            call io_error('error at sitesym_dis_extract_symmetry')
          endif
        endif
        ! choose the larger eigenstate
        umatnew(:, i) = V(1, 2)*umat(:, i) + V(2, 2)*deltaU(:, i)
      enddo ! i
      call symmetrize_ukirr(ik2ir(ik), num_bands, umatnew, n)
      umat(:, :) = umatnew(:, :)
    enddo ! iter

    return
  end subroutine sitesym_dis_extract_symmetry

  !==================================================================!
  subroutine sitesym_read()
    !==================================================================!
    use w90_parameters, only: num_bands, num_wann, num_kpts
    use w90_io, only: io_file_unit, io_error, seedname

    implicit none

    integer :: iu, ibnum, iknum, ierr

    iu = io_file_unit()
    open (unit=iu, file=trim(seedname)//".dmn", form='formatted', status='old', action='read')
    read (iu, *)
    read (iu, *) ibnum, nsymmetry, nkptirr, iknum
    if (ibnum .ne. num_bands) call io_error("Error: Number of bands is not correct (sitesym_read)")
    if (iknum .ne. num_kpts) call io_error("Error: Number of k-points is not correct (sitesym_read)")

    allocate (ik2ir(num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ik2ir in sitesym_read')
    allocate (ir2ik(nkptirr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ir2ik in sitesym_read')
    allocate (kptsym(nsymmetry, nkptirr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating kptsym in sitesym_read')
    allocate (d_matrix_band(num_bands, num_bands, nsymmetry, nkptirr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating d_matrix_band in sitesym_read')
    allocate (d_matrix_wann(num_wann, num_wann, nsymmetry, nkptirr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating d_matrix_wann in sitesym_read')

    read (iu, *) ik2ir
    read (iu, *) ir2ik
    read (iu, *) kptsym
    read (iu, *) d_matrix_wann
    read (iu, *) d_matrix_band
    close (iu)

    return
  end subroutine sitesym_read

  !==================================================================!
  subroutine sitesym_dealloc
    !==================================================================!
    use w90_io, only: io_error

    implicit none

    integer :: ierr

    deallocate (ik2ir, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ik2ir in sitesym_dealloc')
    deallocate (ir2ik, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ir2ik in sitesym_dealloc')
    deallocate (kptsym, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating kptsym in sitesym_dealloc')
    deallocate (d_matrix_band, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating d_matrix_band in sitesym_dealloc')
    deallocate (d_matrix_wann, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating d_matrix_wann in sitesym_dealloc')

    return
  end subroutine sitesym_dealloc

end module w90_sitesym

