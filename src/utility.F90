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

module w90_utility
  !! Module contains lots of useful general routines

  use w90_constants, only: dp

  implicit none

  private

  public :: utility_inv3
  public :: utility_inv2
  public :: utility_det3
  public :: utility_recip_lattice
  public :: utility_metric
  public :: utility_compar
  public :: utility_cart_to_frac
  public :: utility_frac_to_cart
  public :: utility_string_to_coord
  public :: utility_lowercase
  public :: utility_strip
  public :: utility_zgemm
  public :: utility_zgemm_new
  public :: utility_zgemmm
  public :: utility_translate_home
  public :: utility_rotate
  public :: utility_rotate_new
  public :: utility_matmul_diag
  public :: utility_rotate_diag
  public :: utility_commutator_diag
  public :: utility_re_tr
  public :: utility_re_tr_prod
  public :: utility_im_tr
  public :: utility_im_tr_prod
  public :: utility_w0gauss
  public :: utility_w0gauss_vec
  public :: utility_wgauss
  public :: utility_zdotu
  public :: utility_diagonalize

contains

  !=============================================================!
  subroutine utility_zgemm(c, a, transa, b, transb, n)
    !=============================================================!
    !                                                             !
    !! Return matrix product of complex n x n matrices a and b:
    !!
    !!                       C = Op(A) Op(B)
    !!
    !! transa = 'N'  ==> Op(A) = A
    !! transa = 'T'  ==> Op(A) = transpose(A)
    !! transa = 'C'  ==> Op(A) = congj(transpose(A))
    !!
    !! similarly for B
    !                                                             !
    !=============================================================!

    use w90_constants, only: cmplx_0, cmplx_1

    implicit none

    integer, intent(in)  :: n
    character(len=1), intent(in)  :: transa
    character(len=1), intent(in)  :: transb
    complex(kind=dp), intent(in)  :: a(n, n)
    complex(kind=dp), intent(in)  :: b(n, n)
    complex(kind=dp), intent(out) :: c(n, n)

    call zgemm(transa, transb, n, n, n, cmplx_1, a, n, b, n, cmplx_0, c, n)

    return

  end subroutine utility_zgemm

  !===================================================================
  function utility_det3(A)                   !
    !==================================================================!
    !                                                                  !
    !    Return determinant of a  3x3 matrix A                         !
    !                                                                  !
    !===================================================================

    real(kind=dp), intent(in)  :: a(3, 3)
    real(kind=dp)  :: utility_det3
    utility_det3 = A(1, 1)*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2)) + &
                   A(1, 2)*(A(2, 3)*A(3, 1) - A(2, 1)*A(3, 3)) + &
                   A(1, 3)*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
    return
  end function utility_det3

  !=============================================================!
  subroutine utility_zgemm_new(a, b, c, transa_opt, transb_opt)
    !=============================================================!
    !                                                             !
    ! Return matrix product of complex matrices a and b:          !
    !                                                             !
    !                       C = Op(A) Op(B)                       !
    !                                                             !
    ! transa = 'N'  ==> Op(A) = A                                 !
    ! transa = 'T'  ==> Op(A) = transpose(A)                      !
    ! transa = 'C'  ==> Op(A) = congj(transpose(A))               !
    !                                                             !
    ! similarly for B                                             !
    !                                                             !
    ! Due to the use of assumed shape arrays, this routine is a   !
    ! safer and more general replacement for the above routine    !
    ! utility_zgemm. Consider removing utility_zgemm and using    !
    ! utility_zgemm_new throughout.                               !
    !                                                             !
    !=============================================================!

    use w90_constants, only: cmplx_0, cmplx_1

    implicit none

    complex(kind=dp), intent(in)            :: a(:, :)
    complex(kind=dp), intent(in)            :: b(:, :)
    complex(kind=dp), intent(out)           :: c(:, :)
    character(len=1), intent(in), optional  :: transa_opt
    character(len=1), intent(in), optional  :: transb_opt

    integer          :: m, n, k
    character(len=1) :: transa, transb

    transa = 'N'
    transb = 'N'
    if (present(transa_opt)) transa = transa_opt
    if (present(transb_opt)) transb = transb_opt

    ! m ... number of rows in Op(A) and C
    ! n ... number of columns in Op(B) and C
    ! k ... number of columns in Op(A) resp. rows in Op(B)
    m = size(c, 1)
    n = size(c, 2)

    if (transa /= 'N') then
      k = size(a, 1)
    else
      k = size(a, 2)
    end if

    call zgemm(transa, transb, m, n, k, cmplx_1, a, size(a, 1), b, size(b, 1), cmplx_0, c, m)

  end subroutine utility_zgemm_new
  !=============================================================!
  function utility_zdotu(a, b)
    complex(kind=dp), intent(in), dimension(:)  :: a, b
    complex(kind=dp) :: utility_zdotu
    utility_zdotu = sum(a*b)
    return
  end function utility_zdotu

  !=============================================================!
  subroutine utility_zgemmm(a, transa, b, transb, c, transc, &
                            prod1, eigval, prod2)
    !===============================================================!
    ! Returns the complex matrix-matrix-matrix product              !
    ! --> prod1 = op(a).op(b).op(c),                                !
    ! where op(a/b/c) are defined according to transa/transb/transc !
    ! (see also documentation of utility_zgemm above)               !
    !                                                               !
    ! If eigval and prod2 are present, also                         !
    ! --> prod2 = op(a).diag(eigval).op(b).op(c)                    !
    ! is returned.                                                  !
    !===============================================================!

    complex(kind=dp), dimension(:, :), intent(in)  :: a, b, c
    character(len=1), intent(in)                  :: transa, transb, transc
    real(kind=dp), dimension(:), optional, &
      intent(in)       :: eigval
    complex(kind=dp), dimension(:, :), optional, &
      intent(out) :: prod1, prod2

    complex(kind=dp), dimension(:, :), allocatable :: tmp
    integer                                       :: nb, mc, i, j

    ! query matrix sizes
    ! naming convention:
    ! matrix op(a) [resp. op(b) and op(c)] is of size na x ma [resp. nb x mb and nc x mc]
    ! only nb (=ma) and mc are explicitly needed
    if (transb /= 'N') then
      nb = size(b, 2)
    else
      nb = size(b, 1)
    end if
    if (transc /= 'N') then
      mc = size(c, 1)
    else
      mc = size(c, 2)
    end if

    ! tmp = op(b).op(c)
    allocate (tmp(nb, mc))
    call utility_zgemm_new(b, c, tmp, transb, transc)

    ! prod1 = op(a).tmp
    if (present(prod1)) then
      call utility_zgemm_new(a, tmp, prod1, transa, 'N')
    end if

    if (present(prod2) .and. present(eigval)) then
      ! tmp = diag(eigval).tmp
      forall (i=1:nb, j=1:mc)
      tmp(i, j) = eigval(i)*tmp(i, j)
      end forall
      ! prod2 = op(a).tmp
      call utility_zgemm_new(a, tmp, prod2, transa, 'N')
    end if
  end subroutine

  !===================================================================
  subroutine utility_inv3(a, b, det)                   !
    !==================================================================!
    !                                                                  !
    !! Return in b the adjoint of the 3x3 matrix a, and its
    !! determinant.
    !! The inverse is defined as the adjoind divided by the
    !! determinant, so that inverse(a) = b/det
    !                                                                  !
    !===================================================================

    implicit none
    real(kind=dp), intent(in)  :: a(3, 3)
    real(kind=dp), intent(out) :: b(3, 3)
    real(kind=dp), intent(out) :: det

    real(kind=dp):: work(6, 6)
    integer :: i, j, k, l, ll, kk

    do i = 1, 2
      do j = 1, 2
        do k = 1, 3
          do l = 1, 3
            kk = 3*(i - 1) + k
            ll = 3*(j - 1) + l
            work(kk, ll) = a(k, l)
          end do
        end do
      end do
    end do

    det = 0.0_dp
    do i = 1, 3
      det = det + work(1, i)*work(2, i + 1)*work(3, i + 2)
    end do

    do i = 4, 6
      det = det - work(1, i)*work(2, i - 1)*work(3, i - 2)
    end do

    do j = 1, 3
      do i = 1, 3
        b(j, i) = (work(i + 1, j + 1)*work(i + 2, j + 2) - work(i + 1, j + 2) &
                   *work(i + 2, j + 1))
      end do
    end do

    return

  end subroutine utility_inv3

  !===================================================================
  subroutine utility_inv2(a, b, det)                   !
    !==================================================================!
    !                                                                  !
    !! Return in b the adjoint of the 2x2 matrix
    !! a, together with the determinant of a.
    !! The inverse is defined as the adjoind divided by the
    !! determinant, so that inverse(a) = b/det
    !                                                                  !
    !===================================================================

    implicit none
    real(kind=dp), intent(in)  :: a(2, 2)
    real(kind=dp), intent(out) :: b(2, 2)
    real(kind=dp), intent(out) :: det

    det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)

    b(1, 1) = a(2, 2)
    b(1, 2) = -a(1, 2)
    b(2, 1) = -a(2, 1)
    b(2, 2) = a(1, 1)

    return

  end subroutine utility_inv2

  !===================================================================
  subroutine utility_recip_lattice(real_lat, recip_lat, volume)  !
    !==================================================================!
    !                                                                  !
    !!  Calculates the reciprical lattice vectors and the cell volume
    !                                                                  !
    !===================================================================

    use w90_constants, only: dp, twopi, eps5
    use w90_io, only: io_error

    implicit none
    real(kind=dp), intent(in)  :: real_lat(3, 3)
    real(kind=dp), intent(out) :: recip_lat(3, 3)
    real(kind=dp), intent(out) :: volume

    recip_lat(1, 1) = real_lat(2, 2)*real_lat(3, 3) - real_lat(3, 2)*real_lat(2, 3)
    recip_lat(1, 2) = real_lat(2, 3)*real_lat(3, 1) - real_lat(3, 3)*real_lat(2, 1)
    recip_lat(1, 3) = real_lat(2, 1)*real_lat(3, 2) - real_lat(3, 1)*real_lat(2, 2)
    recip_lat(2, 1) = real_lat(3, 2)*real_lat(1, 3) - real_lat(1, 2)*real_lat(3, 3)
    recip_lat(2, 2) = real_lat(3, 3)*real_lat(1, 1) - real_lat(1, 3)*real_lat(3, 1)
    recip_lat(2, 3) = real_lat(3, 1)*real_lat(1, 2) - real_lat(1, 1)*real_lat(3, 2)
    recip_lat(3, 1) = real_lat(1, 2)*real_lat(2, 3) - real_lat(2, 2)*real_lat(1, 3)
    recip_lat(3, 2) = real_lat(1, 3)*real_lat(2, 1) - real_lat(2, 3)*real_lat(1, 1)
    recip_lat(3, 3) = real_lat(1, 1)*real_lat(2, 2) - real_lat(2, 1)*real_lat(1, 2)

    volume = real_lat(1, 1)*recip_lat(1, 1) + &
             real_lat(1, 2)*recip_lat(1, 2) + &
             real_lat(1, 3)*recip_lat(1, 3)

    if (abs(volume) < eps5) then
      call io_error(' Found almost zero Volume in utility_recip_lattice')
    end if

    recip_lat = twopi*recip_lat/volume
    volume = abs(volume)

    return

  end subroutine utility_recip_lattice

  !===================================================================
  subroutine utility_compar(a, b, ifpos, ifneg)
    !==================================================================!
    !                                                                  !
    !! Compares two vectors
    !                                                                  !
    !===================================================================
    use w90_constants, only: eps8

    implicit none

    real(kind=dp), intent(in) :: a(3)
    real(kind=dp), intent(in) :: b(3)
    integer, intent(out) :: ifpos, ifneg

    real(kind=dp) :: rrp, rrm

    rrp = (a(1) - b(1))**2 + (a(2) - b(2))**2 + (a(3) - b(3))**2
    rrm = (a(1) + b(1))**2 + (a(2) + b(2))**2 + (a(3) + b(3))**2
    ifpos = 0
    if (abs(rrp) .lt. eps8) ifpos = 1
    ifneg = 0
    if (abs(rrm) .lt. eps8) ifneg = 1

    return

  end subroutine utility_compar

  !===================================================================
  subroutine utility_metric(real_lat, recip_lat, &
                            real_metric, recip_metric)
    !==================================================================!
    !                                                                  !
    !!  Calculate the real and reciprical space metrics
    !                                                                  !
    !===================================================================
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3, 3)
    real(kind=dp), intent(in)  :: recip_lat(3, 3)
    real(kind=dp), intent(out) :: real_metric(3, 3)
    real(kind=dp), intent(out) :: recip_metric(3, 3)

    integer :: i, j, l

    real_metric = 0.0_dp; recip_metric = 0.0_dp

    do j = 1, 3
      do i = 1, j
        do l = 1, 3
          real_metric(i, j) = real_metric(i, j) + real_lat(i, l)*real_lat(j, l)
          recip_metric(i, j) = recip_metric(i, j) + recip_lat(i, l)*recip_lat(j, l)
        enddo
        if (i .lt. j) then
          real_metric(j, i) = real_metric(i, j)
          recip_metric(j, i) = recip_metric(i, j)
        endif
      enddo
    enddo

  end subroutine utility_metric

  !===================================================================
  subroutine utility_frac_to_cart(frac, cart, real_lat)
    !==================================================================!
    !                                                                  !
    !!  Convert from fractional to Cartesian coordinates
    !                                                                  !
    !===================================================================
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3, 3)
    real(kind=dp), intent(in)  :: frac(3)
    real(kind=dp), intent(out) :: cart(3)

    integer :: i

    do i = 1, 3
      cart(i) = real_lat(1, i)*frac(1) + real_lat(2, i)*frac(2) + real_lat(3, i)*frac(3)
    end do

    return

  end subroutine utility_frac_to_cart

  !===================================================================
  subroutine utility_cart_to_frac(cart, frac, recip_lat)
    !==================================================================!
    !                                                                  !
    !!  Convert from Cartesian to fractional coordinates
    !                                                                  !
    !===================================================================
    use w90_constants, only: twopi
    implicit none

    real(kind=dp), intent(in)  :: recip_lat(3, 3)
    real(kind=dp), intent(out)  :: frac(3)
    real(kind=dp), intent(in)  :: cart(3)

    integer :: i

    do i = 1, 3
      frac(i) = recip_lat(i, 1)*cart(1) + recip_lat(i, 2)*cart(2) + recip_lat(i, 3)*cart(3)
    end do

    frac = frac/twopi

    return

  end subroutine utility_cart_to_frac

  !=============================!
  function utility_strip(string)!
    !=============================!
    !                             !
    !! Strips string of all blank spaces
    !                             !
    !=============================!

    use w90_io, only: maxlen

    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen) :: utility_strip

    integer :: ispc, ipos, ilett, icount

    ! Initialise
    utility_strip = repeat(' ', maxlen)

    ispc = ichar(' ')
    icount = 0
    do ipos = 1, len(string)
      ilett = ichar(string(ipos:ipos))
      if (ilett .ne. ispc) then
        icount = icount + 1
        utility_strip(icount:icount) = string(ipos:ipos)
      endif
    enddo

    utility_strip = trim(utility_strip)

    return

  end function utility_strip

  !=================================!
  function utility_lowercase(string)!
    !=================================!
    !                                 !
    !! Takes a string and converts to
    !!  lowercase characters
    !                                 !
    !=================================!

    use w90_io, only: maxlen

    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen) :: utility_lowercase

    integer :: iA, iZ, idiff, ipos, ilett

    iA = ichar('A')
    iZ = ichar('Z')
    idiff = iZ - ichar('z')

    utility_lowercase = string

    do ipos = 1, len(string)
      ilett = ichar(string(ipos:ipos))
      if ((ilett .ge. iA) .and. (ilett .le. iZ)) &
        utility_lowercase(ipos:ipos) = char(ilett - idiff)
    enddo

    utility_lowercase = trim(adjustl(utility_lowercase))

    return

  end function utility_lowercase

  !====================================================!
  subroutine utility_string_to_coord(string_tmp, outvec)!
    !====================================================!
    !                                                    !
    !! Takes a string in the form 0.0,1.0,0.5
    !! and returns an array of the real num
    !                                                    !
    !====================================================!
    use w90_io, only: io_error, maxlen

    implicit none

    character(len=maxlen), intent(in)  :: string_tmp
    real(kind=dp), intent(out) :: outvec(3)

    integer :: pos
    character(len=maxlen)  :: ctemp
    character(len=maxlen)  :: ctemp2

    ctemp = string_tmp
    pos = index(ctemp, ',')
    if (pos <= 0) call io_error('utility_string_to_coord: Problem reading string into real number '//trim(string_tmp))
    ctemp2 = ctemp(1:pos - 1)
    read (ctemp2, *, err=100, end=100) outvec(1)
    ctemp = ctemp(pos + 1:)
    pos = index(ctemp, ',')
    ctemp2 = ctemp(1:pos - 1)
    read (ctemp2, *, err=100, end=100) outvec(2)
    ctemp = ctemp(pos + 1:)
    read (ctemp, *, err=100, end=100) outvec(3)

    return

100 call io_error('utility_string_to_coord: Problem reading string into real number '//trim(string_tmp))

  end subroutine utility_string_to_coord

!~  !===========================================!
!~  function utility_string_to_coord(string_tmp)!
!~  !===========================================!
!~  !                                           !
!~  !  Takes a string in the form 0.0,1.0,0.5   !
!~  !   and returns an array of the real num    !
!~  !                                           !
!~  !===========================================!
!~
!~    implicit none
!~
!~    character(len=80), intent(in)  :: string_tmp
!~    real(kind=dp) :: utility_string_to_coord(3)
!~
!~    integer :: pos,pos2
!~    character(len=80)  :: ctemp
!~    character(len=80)  :: ctemp2
!~
!~
!~    ctemp=string_tmp
!~    pos2=index(ctemp,',')
!~    ctemp2=ctemp(1:pos2-1)
!~    read(ctemp2,*) utility_string_to_coord(1)
!~    ctemp=ctemp(pos2+1:)
!~    pos2=index(ctemp,',')
!~    ctemp2=ctemp(1:pos2-1)
!~    read(ctemp2,*) utility_string_to_coord(2)
!~    ctemp=ctemp(pos2+1:)
!~    read(ctemp,*) utility_string_to_coord(3)
!~
!~  end function utility_string_to_coord

  !========================================================!
  subroutine utility_translate_home(vec, real_lat, recip_lat)
    !========================================================!
    !                                                        !
    !! Translate a vector to the home unit cell
    !                                                        !
    !========================================================!

    implicit none

    real(kind=dp), intent(inout) :: vec(3)
    real(kind=dp), intent(in)    :: real_lat(3, 3)
    real(kind=dp), intent(in)    :: recip_lat(3, 3)

    ! <<<local variables>>>
    integer       :: ind
    real(kind=dp) :: r_home(3), r_frac(3)
    real(kind=dp) :: shift

    r_home = 0.0_dp; r_frac = 0.0_dp

    ! Cartesian --> fractional
    call utility_cart_to_frac(vec, r_frac, recip_lat)
    ! Rationalise to interval [0,1]
    do ind = 1, 3
      if (r_frac(ind) .lt. 0.0_dp) then
        shift = real(ceiling(abs(r_frac(ind))), kind=dp)
        r_frac(ind) = r_frac(ind) + shift
      endif
      if (r_frac(ind) .gt. 1.0_dp) then
        shift = -real(int(r_frac(ind)), kind=dp)
        r_frac(ind) = r_frac(ind) + shift
      endif
    enddo
    ! Fractional --> Cartesian
    call utility_frac_to_cart(r_frac, r_home, real_lat)

    vec = r_home

    return
  end subroutine utility_translate_home

  !============================================================!
  subroutine utility_diagonalize(mat, dim, eig, rot)
    !============================================================!
    !                                                            !
    !! Diagonalize the dim x dim  hermitian matrix 'mat' and
    !! return the eigenvalues 'eig' and the unitary rotation 'rot'
    !                                                            !
    !============================================================!

    use w90_constants, only: dp, cmplx_0
    use w90_io, only: io_error, stdout

    integer, intent(in)           :: dim
    complex(kind=dp), intent(in)  :: mat(dim, dim)
    real(kind=dp), intent(out)    :: eig(dim)
    complex(kind=dp), intent(out) :: rot(dim, dim)

    complex(kind=dp)   :: mat_pack((dim*(dim + 1))/2), cwork(2*dim)
    real(kind=dp)      :: rwork(7*dim)
    integer            :: i, j, info, nfound, iwork(5*dim), ifail(dim)

    do j = 1, dim
      do i = 1, j
        mat_pack(i + ((j - 1)*j)/2) = mat(i, j)
      enddo
    enddo
    rot = cmplx_0; eig = 0.0_dp; cwork = cmplx_0; rwork = 0.0_dp; iwork = 0
    call ZHPEVX('V', 'A', 'U', dim, mat_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                nfound, eig(1), rot, dim, cwork, rwork, iwork, ifail, info)
    if (info < 0) then
      write (stdout, '(a,i3,a)') 'THE ', -info, &
        ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
      call io_error('Error in utility_diagonalize')
    endif
    if (info > 0) then
      write (stdout, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
      call io_error('Error in utility_diagonalize')
    endif

  end subroutine utility_diagonalize

  !===========================================================!
  function utility_rotate(mat, rot, dim)
    !==========================================================!
    !                                                           !
    !! Rotates the dim x dim matrix 'mat' according to
    !! (rot)^dagger.mat.rot, where 'rot' is a unitary matrix
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp

    integer          :: dim
    complex(kind=dp) :: utility_rotate(dim, dim)
    complex(kind=dp) :: mat(dim, dim)
    complex(kind=dp) :: rot(dim, dim)

    utility_rotate = matmul(matmul(transpose(conjg(rot)), mat), rot)

  end function utility_rotate

  !===========================================================!
  subroutine utility_rotate_new(mat, rot, N, reverse)
    !==============================================================!
    !                                                              !
    ! Rotates the N x N matrix 'mat' according to                  !
    ! * (rot)^dagger.mat.rot (reverse = .false. or not present) OR !
    ! * rot.mat.(rot)^dagger (reverse = .true.),                   !
    ! where 'rot' is a unitary matrix.                             !
    ! The matrix 'mat' is overwritten.                             !
    !                                                              !
    !==============================================================!

    use w90_constants, only: dp

    integer, intent(in)             :: N
    logical, optional, intent(in)   :: reverse
    complex(kind=dp), intent(inout) :: mat(N, N)
    complex(kind=dp), intent(in)    :: rot(N, N)
    complex(kind=dp)                :: tmp(N, N)
    logical                         :: rev

    if (.not. present(reverse)) then
      rev = .false.
    else
      rev = reverse
    end if

    if (rev) then
      call utility_zgemm_new(rot, mat, tmp, 'N', 'C')
      call utility_zgemm_new(rot, tmp, mat, 'N', 'C')
    else
      call utility_zgemm_new(mat, rot, tmp, 'C', 'N')
      call utility_zgemm_new(tmp, rot, mat, 'C', 'N')
    end if

  end subroutine utility_rotate_new

  !===========================================================!
  function utility_matmul_diag(mat1, mat2, dim)
    !===========================================================!
    !                                                           !
    !! Computes the diagonal elements of the matrix mat1.mat2
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp, cmplx_0

    integer          :: dim
    complex(kind=dp) :: utility_matmul_diag(dim)
    complex(kind=dp) :: mat1(dim, dim)
    complex(kind=dp) :: mat2(dim, dim)

    integer i, j

    utility_matmul_diag = cmplx_0
    do i = 1, dim
      do j = 1, dim
        utility_matmul_diag(i) = utility_matmul_diag(i) + mat1(i, j)*mat2(j, i)
      end do
    end do

  end function utility_matmul_diag

  !===========================================================!
  function utility_rotate_diag(mat, rot, dim)
    !===========================================================!
    !                                                           !
    !! Rotates the dim x dim matrix 'mat' according to
    !! (rot)^dagger.mat.rot, where 'rot' is a unitary matrix.
    !! Computes only the diagonal elements of rotated matrix.
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp

    integer          :: dim
    complex(kind=dp) :: utility_rotate_diag(dim)
    complex(kind=dp) :: mat(dim, dim)
    complex(kind=dp) :: rot(dim, dim)
    complex(kind=dp) :: tmp(dim, dim)

    call utility_zgemm_new(rot, mat, tmp, 'C', 'N')
    utility_rotate_diag = utility_matmul_diag(tmp, rot, dim)

  end function utility_rotate_diag

  !===========================================================!
  function utility_commutator_diag(mat1, mat2, dim)
    !===========================================================!
    !                                                           !
    !! Computes diagonal elements of
    !! [mat1,mat2]=mat1.mat2-mat2.mat1
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp

    integer          :: dim
    complex(kind=dp) :: utility_commutator_diag(dim)
    complex(kind=dp) :: mat1(dim, dim)
    complex(kind=dp) :: mat2(dim, dim)

    utility_commutator_diag = utility_matmul_diag(mat1, mat2, dim) - utility_matmul_diag(mat2, mat1, dim)

  end function utility_commutator_diag

  !===================================================!
  function utility_re_tr_prod(a, b)
    !================================================!
    !                                                !
    ! Return Re(tr(a.b)), i.e. the real part of the  !
    ! trace of the matrix product of a and b.        !
    !                                                !
    !================================================!
    use w90_constants, only: dp, cmplx_0, cmplx_i

    complex(kind=dp), dimension(:, :), intent(in) :: a, b
    real(kind=dp) :: utility_re_tr_prod
    real(kind=dp) :: s
    integer       :: i, j, n, m

    n = min(size(a, 1), size(b, 2))
    m = min(size(a, 2), size(b, 1))

    s = 0
    do i = 1, n
      do j = 1, m
        s = s + dble(a(i, j)*b(j, i))
      end do
    end do
    utility_re_tr_prod = s
  end function

  !===================================================!
  function utility_im_tr_prod(a, b)
    !====================================================!
    !                                                    !
    ! Return Im(tr(a.b)), i.e. the imaginary part of the !
    ! trace of the matrix product of a and b.            !
    !                                                    !
    !====================================================!
    use w90_constants, only: dp, cmplx_0, cmplx_i

    complex(kind=dp), dimension(:, :), intent(in) :: a, b

    real(kind=dp) :: utility_im_tr_prod
    real(kind=dp) :: s
    integer       :: i, j, n, m

    n = min(size(a, 1), size(b, 2))
    m = min(size(a, 2), size(b, 1))

    s = 0
    do i = 1, n
      do j = 1, m
        s = s + aimag(a(i, j)*b(j, i))
      end do
    end do
    utility_im_tr_prod = s
  end function

  !===================================================!
  function utility_re_tr(mat)
    !========================!
    !                        !
    !! Real part of the trace
    !                        !
    !========================!

    use w90_constants, only: dp, cmplx_0, cmplx_i

    real(kind=dp)                    :: utility_re_tr
    complex(kind=dp), dimension(:, :) :: mat

    integer          :: i, mydim
    complex(kind=dp) :: cdum

    mydim = size(mat, 1)

    cdum = cmplx_0
    do i = 1, mydim
      cdum = cdum + mat(i, i)
    enddo
    utility_re_tr = aimag(cmplx_i*cdum)

  end function utility_re_tr

  function utility_im_tr(mat)
    !=============================!
    !                             !
    !! Imaginary part of the trace
    !                             !
    !=============================!

    use w90_constants, only: dp, cmplx_0

    real(kind=dp)                    :: utility_im_tr
    complex(kind=dp), dimension(:, :) :: mat

    integer          :: i, mydim
    complex(kind=dp) :: cdum

    mydim = size(mat, 1)

    cdum = cmplx_0
    do i = 1, mydim
      cdum = cdum + mat(i, i)
    enddo
    utility_im_tr = aimag(cdum)

  end function utility_im_tr

  function utility_wgauss(x, n)
    !-----------------------------------------------------------------------
    !
    !! this function computes the approximate theta function for the
    !! given order n, at the point x.
    !!
    !! (n>=0) : Methfessel-Paxton case. See PRB 40, 3616 (1989).
    !!
    !! (n=-1 ): Cold smearing (Marzari-Vanderbilt). See PRL 82, 3296 (1999)
    !!       1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
    !!
    !! (n=-99): Fermi-Dirac case: 1.0/(1.0+exp(-x)).
    !
    use w90_constants, only: dp, pi

    implicit none
    real(kind=dp) :: utility_wgauss, x
    !! output: the value of the function
    !! input: the argument of the function
    integer :: n
    !! input: the order of the function
    !
    !    the local variables
    !

    real(kind=dp) :: a, hp, arg, hd, xp
    ! the coefficient a_n
    ! the hermitean function
    ! the argument of the exponential
    ! the hermitean function
    ! auxiliary variable (cold smearing)
    integer :: i, ni
    ! counter on the n indices
    ! counter on 2n
!    real(kind=dp), external :: gauss_freq, qe_erf
    real(kind=dp), parameter :: maxarg = 200.0_dp
    ! maximum value for the argument of the exponential

    ! Fermi-Dirac smearing
    if (n .eq. -99) then
      if (x .lt. -maxarg) then
        utility_wgauss = 0.0_dp
      elseif (x .gt. maxarg) then
        utility_wgauss = 1.0_dp
      else
        utility_wgauss = 1.00_dp/(1.00_dp + exp(-x))
      endif
      return

    endif
    ! Cold smearing
    if (n .eq. -1) then
      xp = x - 1.00_dp/sqrt(2.00_dp)
      arg = min(maxarg, xp**2)
      utility_wgauss = 0.50_dp*qe_erf(xp) + 1.00_dp/sqrt(2.00_dp*pi)*exp(- &
                                                                         arg) + 0.50_dp
      return

    endif
    ! Methfessel-Paxton
    utility_wgauss = gauss_freq(x*sqrt(2.00_dp))
    if (n .eq. 0) return
    hd = 0.0_dp
    arg = min(maxarg, x**2)
    hp = exp(-arg)
    ni = 0
    a = 1.0_dp/sqrt(pi)
    do i = 1, n
      hd = 2.00_dp*x*hp - 2.00_dp*DBLE(ni)*hd
      ni = ni + 1
      a = -a/(DBLE(i)*4.00_dp)
      utility_wgauss = utility_wgauss - a*hd
      hp = 2.00_dp*x*hd - 2.00_dp*DBLE(ni)*hp
      ni = ni + 1
    enddo
    return
  end function utility_wgauss

  function utility_w0gauss(x, n)
    !-----------------------------------------------------------------------
    !
    !! the derivative of utility_wgauss:  an approximation to the delta function
    !!
    !! (n>=0) : derivative of the corresponding Methfessel-Paxton utility_wgauss
    !!
    !! (n=-1 ): derivative of cold smearing:
    !!              1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
    !!
    !! (n=-99): derivative of Fermi-Dirac function: 0.5/(1.0+cosh(x))
    !
    use w90_constants, only: dp, pi
    use w90_io, only: io_error
    implicit none
    real(kind=dp) :: utility_w0gauss, x
    !! output: the value of the function
    !! input: the point where to compute the function

    integer :: n
    !! input: the order of the smearing function
    !
    !    here the local variables
    !
    real(kind=dp) :: a, arg, hp, hd, sqrtpm1
    ! the coefficients a_n
    ! the argument of the exponential
    ! the hermite function
    ! the hermite function

    integer :: i, ni
    ! counter on n values
    ! counter on 2n values

    ! Fermi-Dirac smearing

    sqrtpm1 = 1.0_dp/sqrt(pi)

    if (n .eq. -99) then
      if (abs(x) .le. 36.0) then
        utility_w0gauss = 1.00_dp/(2.00_dp + exp(-x) + exp(+x))
        ! in order to avoid problems for large values of x in the e
      else
        utility_w0gauss = 0.0_dp
      endif
      return

    endif
    ! cold smearing  (Marzari-Vanderbilt)
    if (n .eq. -1) then
      arg = min(200.0_dp, (x - 1.00_dp/sqrt(2.00_dp))**2)
      utility_w0gauss = sqrtpm1*exp(-arg)*(2.00_dp - sqrt(2.00_dp)*x)
      return

    endif

    if (n .gt. 10 .or. n .lt. 0) &
      call io_error('utility_w0gauss higher order smearing is untested and unstable')

    ! Methfessel-Paxton
    arg = min(200.0_dp, x**2)
    utility_w0gauss = exp(-arg)*sqrtpm1
    if (n .eq. 0) return
    hd = 0.00_dp
    hp = exp(-arg)
    ni = 0
    a = sqrtpm1
    do i = 1, n
      hd = 2.00_dp*x*hp - 2.00_dp*DBLE(ni)*hd
      ni = ni + 1
      a = -a/(DBLE(i)*4.00_dp)
      hp = 2.00_dp*x*hd - 2.00_dp*DBLE(ni)*hp
      ni = ni + 1
      utility_w0gauss = utility_w0gauss + a*hp
    enddo
    return
  end function utility_w0gauss

  function utility_w0gauss_vec(x, n) result(res)
    !-----------------------------------------------------------------------
    !  Stepan Tsirkin: a vectorized version of the outine, gets x as an array.
    !
    !! the derivative of utility_wgauss:  an approximation to the delta function
    !!
    !! (n>=0) : derivative of the corresponding Methfessel-Paxton utility_wgauss
    !!
    !! (n=-1 ): derivative of cold smearing:
    !!              1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
    !!
    !! (n=-99): derivative of Fermi-Dirac function: 0.5/(1.0+cosh(x))
    !
    use w90_constants, only: dp, pi
    use w90_io, only: io_error
    implicit none
    real(kind=dp), intent(in)   ::  x(:)
    real(kind=dp), allocatable  :: res(:), arg(:)

    !! output: the value of the function
    !! input: the point where to compute the function

    integer :: n
    !! input: the order of the smearing function
    !
    !    here the local variables
    !
    real(kind=dp) :: sqrtpm1

    allocate (res(size(x)))
    allocate (arg(size(x)))
    sqrtpm1 = 1.0_dp/sqrt(pi)

    if (n .eq. -99) then
      call io_error('utility_w0gauss_vec not implemented for n == 99')
    endif

    ! cold smearing  (Marzari-Vanderbilt)
    if (n .eq. -1) then
      call io_error('utility_w0gauss_vec not implemented for n == -1')
    endif

    if (n .gt. 10 .or. n .lt. 0) &
      call io_error('utility_w0gauss higher order smearing is untested and unstable')

    ! Methfessel-Paxton
    arg = min(200.0_dp, x**2)
    res = exp(-arg)*sqrtpm1
    if (n .eq. 0) return
    call io_error('utility_w0gauss_vec not implemented for n >0 ')
    return
  end function utility_w0gauss_vec

  function qe_erf(x)
    !---------------------------------------------------------------------
    !
    !! Error function - computed from the rational approximations of
    !! W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
    !!
    !!     for abs(x) le 0.47 erf is calculated directly
    !!     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
    !
    use w90_constants, only: dp

    implicit none
    real(kind=dp), intent(in) :: x
    real(kind=dp) :: x2, p1(4), q1(4)
    real(kind=dp) :: qe_erf
    data p1/2.426679552305318E2_dp, 2.197926161829415E1_dp, &
      6.996383488619136_dp, -3.560984370181538E-2_dp/
    data q1/2.150588758698612E2_dp, 9.116490540451490E1_dp, &
      1.508279763040779E1_dp, 1.000000000000000_dp/
    !
    if (abs(x) > 6.0_dp) then
      !
      !  erf(6)=1-10^(-17) cannot be distinguished from 1
      !
      qe_erf = sign(1.0_dp, x)
    else
      if (abs(x) <= 0.47_dp) then
        x2 = x**2
        qe_erf = x*(p1(1) + x2*(p1(2) + x2*(p1(3) + x2*p1(4)))) &
                 /(q1(1) + x2*(q1(2) + x2*(q1(3) + x2*q1(4))))
      else
        qe_erf = 1.0_dp - qe_erfc(x)
      endif
    endif
    !
    return
  end function qe_erf
  !
  !---------------------------------------------------------------------
  function qe_erfc(x)
    !---------------------------------------------------------------------
    !
    !! erfc(x) = 1-erf(x)  - See comments in erf
    !
    use w90_constants, only: dp
    implicit none
    real(kind=dp), intent(in) :: x
    real(kind=dp)            :: qe_erfc
    real(kind=dp) :: ax, x2, xm2, p2(8), q2(8), p3(5), q3(5), pim1
    data p2/3.004592610201616E2_dp, 4.519189537118719E2_dp, &
      3.393208167343437E2_dp, 1.529892850469404E2_dp, &
      4.316222722205674E1_dp, 7.211758250883094_dp, &
      5.641955174789740E-1_dp, -1.368648573827167E-7_dp/
    data q2/3.004592609569833E2_dp, 7.909509253278980E2_dp, &
      9.313540948506096E2_dp, 6.389802644656312E2_dp, &
      2.775854447439876E2_dp, 7.700015293522947E1_dp, &
      1.278272731962942E1_dp, 1.000000000000000_dp/
    data p3/-2.996107077035422E-3_dp, -4.947309106232507E-2_dp, &
      -2.269565935396869E-1_dp, -2.786613086096478E-1_dp, &
      -2.231924597341847E-2_dp/
    data q3/1.062092305284679E-2_dp, 1.913089261078298E-1_dp, &
      1.051675107067932_dp, 1.987332018171353_dp, &
      1.000000000000000_dp/

    data pim1/0.56418958354775629_dp/
    !        ( pim1= sqrt(1/pi) )
    ax = abs(x)
    if (ax > 26.0_dp) then
      !
      !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
      !
      qe_erfc = 0.0_dp
    elseif (ax > 4.0_dp) then
      x2 = x**2
      xm2 = (1.0_dp/ax)**2
      qe_erfc = (1.0_dp/ax)*exp(-x2)*(pim1 + xm2*(p3(1) + xm2*(p3(2) + xm2*(p3(3) + xm2*(p3(4) + xm2*p3(5)))))/ &
                                      (q3(1) + xm2*(q3(2) + xm2*(q3(3) + xm2*(q3(4) + xm2*q3(5))))))
    elseif (ax > 0.47_dp) then
      x2 = x**2
      qe_erfc = exp(-x2)* &
                (p2(1) + ax*(p2(2) + ax*(p2(3) + ax*(p2(4) + ax*(p2(5) + ax*(p2(6) + ax*(p2(7) + ax*p2(8))))))))/ &
                (q2(1) + ax*(q2(2) + ax*(q2(3) + ax*(q2(4) + ax*(q2(5) + ax*(q2(6) + ax*(q2(7) + ax*q2(8))))))))
    else
      qe_erfc = 1.0_dp - qe_erf(ax)
    endif
    !
    ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
    !
    if (x < 0.0_dp) qe_erfc = 2.0_dp - qe_erfc
    !
    return
  end function qe_erfc
  !
  !---------------------------------------------------------------------
  function gauss_freq(x)
    !---------------------------------------------------------------------
    !
    !! gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
    !! - See comments in erf
    !
    use w90_constants, only: dp

    implicit none
    real(kind=dp), intent(in) :: x
    real(kind=dp)            :: gauss_freq
    real(kind=dp), parameter :: c = 0.7071067811865475_dp
    !        ( c= sqrt(1/2) )
    !
    gauss_freq = 0.5_dp*qe_erfc(-x*c)
    !
    return
  end function gauss_freq

end module w90_utility
