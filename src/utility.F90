!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,          !
!  Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

module w90_utility

  use w90_constants, only : dp

  implicit none

  private

  public :: utility_inv3
  public :: utility_recip_lattice
  public :: utility_compar
  public :: utility_compute_metric
  public :: utility_cart_to_frac
  public :: utility_frac_to_cart
  public :: utility_string_to_coord
  public :: utility_lowercase
  public :: utility_strip
  public :: utility_zgemm
  public :: utility_translate_home
  public :: utility_rotate
  public :: utility_matmul_diag
  public :: utility_rotate_diag
  public :: utility_commutator_diag
  public :: utility_re_tr
  public :: utility_im_tr
  public :: w0gauss
  public :: wgauss
  public :: utility_diagonalize


contains


  !=============================================================!
  subroutine utility_zgemm(c,a,transa,b,transb,n)
    !=============================================================!
    !                                                             !
    ! Return matrix product of complex n x n matrices a and b:    !
    !                                                             !
    !                       C = Op(A) Op(B)                       !
    !                                                             !
    ! transa = 'N'  ==> Op(A) = A                                 !
    ! transa = 'T'  ==> Op(A) = transpose(A)                      !
    ! transa = 'C'  ==> Op(A) = congj(transpose(A))               !
    !                                                             !
    ! similarly for B                                             !
    !                                                             !
    !=============================================================!

    use w90_constants, only: cmplx_0,cmplx_1

    implicit none

    integer,           intent(in)  :: n
    character(len=1),  intent(in)  :: transa
    character(len=1),  intent(in)  :: transb
    complex(kind=dp),  intent(in)  :: a(n,n)
    complex(kind=dp),  intent(in)  :: b(n,n)
    complex(kind=dp),  intent(out) :: c(n,n)

    call zgemm(transa,transb,n,n,n,cmplx_1,a,n,b,n,cmplx_0,c,n)

    return

  end subroutine utility_zgemm



    !===================================================================
                 subroutine utility_inv3 (a, b, det)                   !
    !==================================================================!
    !                                                                  !
    !    Return the inverse of a 3x3 matrix and its determinant        !
    !                                                                  !
    !===================================================================


    implicit none
    real(kind=dp), intent(in)  :: a (3, 3)
    real(kind=dp), intent(out) :: b (3, 3)  
    real(kind=dp), intent(out) :: det

    real(kind=dp):: work (6, 6)  
    integer :: i,j,k,l,ll,kk


    do  i = 1, 2  
       do  j = 1, 2  
          do  k = 1, 3  
             do  l = 1, 3  
                kk = 3 * (i - 1) + k  
                ll = 3 * (j - 1) + l  
                work (kk, ll) = a (k, l)  
             end do
          end do
       end do
    end do

    det = 0.0_dp
    do i = 1, 3  
       det = det + work (1, i) * work (2, i + 1) * work (3, i + 2)  
    end do

    do i = 4, 6  
       det = det - work (1, i) * work (2, i - 1) * work (3, i - 2)  
    end do

    do j = 1, 3  
       do i = 1, 3  
          b (j, i) = (work (i + 1, j + 1) * work (i + 2, j + 2) - work (i + 1, j + 2) &
               * work (i + 2, j + 1) )
       end do
    end do

    return

  end subroutine utility_inv3


    !===================================================================
         subroutine utility_recip_lattice (real_lat,recip_lat,volume)  !
    !==================================================================!
    !                                                                  !
    !  Calculates the reciprical lattice vectors and the cell volume   !
    !                                                                  !
    !===================================================================

    use w90_constants,  only : dp,twopi,eps5
    use w90_io,         only : io_error

    implicit none
    real(kind=dp), intent(in)  :: real_lat (3, 3)
    real(kind=dp), intent(out) :: recip_lat (3, 3)  
    real(kind=dp), intent(out) :: volume

    recip_lat(1,1)=real_lat(2,2)*real_lat(3,3)-real_lat(3,2)*real_lat(2,3)
    recip_lat(1,2)=real_lat(2,3)*real_lat(3,1)-real_lat(3,3)*real_lat(2,1)
    recip_lat(1,3)=real_lat(2,1)*real_lat(3,2)-real_lat(3,1)*real_lat(2,2)
    recip_lat(2,1)=real_lat(3,2)*real_lat(1,3)-real_lat(1,2)*real_lat(3,3)
    recip_lat(2,2)=real_lat(3,3)*real_lat(1,1)-real_lat(1,3)*real_lat(3,1)
    recip_lat(2,3)=real_lat(3,1)*real_lat(1,2)-real_lat(1,1)*real_lat(3,2)
    recip_lat(3,1)=real_lat(1,2)*real_lat(2,3)-real_lat(2,2)*real_lat(1,3)
    recip_lat(3,2)=real_lat(1,3)*real_lat(2,1)-real_lat(2,3)*real_lat(1,1)
    recip_lat(3,3)=real_lat(1,1)*real_lat(2,2)-real_lat(2,1)*real_lat(1,2)

    volume=real_lat(1,1)*recip_lat(1,1) + &
         real_lat(1,2)*recip_lat(1,2) + &
         real_lat(1,3)*recip_lat(1,3)  


    if( abs(volume) < eps5 ) then
       call io_error(' Found almost zero Volume in utility_recip_lattice')
    end if

    recip_lat=twopi*recip_lat/volume
    volume=abs(volume)

    return

  end subroutine utility_recip_lattice


    !===================================================================
  subroutine utility_compar(a,b,ifpos,ifneg)
    !==================================================================!
    !                                                                  !
    !                                                                  !
    !                                                                  !
    !===================================================================
    use w90_constants, only: eps8

    implicit none

    real(kind=dp), intent(in) :: a(3)
    real(kind=dp), intent(in) :: b(3)
    integer, intent(out) :: ifpos,ifneg

    real(kind=dp) :: rrp,rrm

    rrp=(a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2
    rrm=(a(1)+b(1))**2+(a(2)+b(2))**2+(a(3)+b(3))**2
    ifpos=0
    if (abs(rrp).lt.eps8) ifpos=1
    ifneg=0
    if (abs(rrm).lt.eps8) ifneg=1

    return

  end subroutine utility_compar



    !===================================================================
         subroutine utility_compute_metric(real_lat,recip_lat, &
                                              real_metric,recip_metric)
    !==================================================================!
    !                                                                  !
    !  Calculate the real and reciprical space metrics                 !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out) :: real_metric(3,3)
    real(kind=dp), intent(out) :: recip_metric(3,3)

    integer :: i,j,l

    real_metric=0.0_dp ; recip_metric=0.0_dp

    do j=1,3
       do i=1,j
          do l=1,3
             real_metric(i,j)=real_metric(i,j)+real_lat(i,l)*real_lat(j,l)
             recip_metric(i,j)=recip_metric(i,j)+recip_lat(i,l)*recip_lat(j,l)
          enddo
          if(i.lt.j) then
             real_metric(j,i)=real_metric(i,j)
             recip_metric(j,i)=recip_metric(i,j)
          endif
       enddo
    enddo

  end subroutine utility_compute_metric



    !===================================================================
         subroutine utility_frac_to_cart(frac,cart,real_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: frac(3)
    real(kind=dp), intent(out) :: cart(3)

    integer :: i

    do i=1,3
       cart(i)=real_lat(1,i)*frac(1) + real_lat(2,i)*frac(2) + real_lat(3,i)*frac(3) 
    end do

    return

  end subroutine utility_frac_to_cart


    !===================================================================
         subroutine utility_cart_to_frac(cart,frac,recip_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    use w90_constants, only : twopi
    implicit none

    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out)  :: frac(3)
    real(kind=dp), intent(in)  :: cart(3)

    integer :: i

    do i=1,3
       frac(i)=recip_lat(i,1)*cart(1) + recip_lat(i,2)*cart(2) + recip_lat(i,3)*cart(3) 
    end do

    frac=frac/twopi


    return

  end subroutine utility_cart_to_frac


  !=============================!
  function utility_strip(string)!
  !=============================!
  !                             !
  !    Strips string of all     !
  !        blank spaces         !
  !                             !
  !=============================!

    use w90_io, only : maxlen

    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen) :: utility_strip

    integer :: ispc,ipos,ilett,icount

    ! Initialise
    utility_strip=repeat(' ',maxlen)

    ispc = ichar(' ')
    icount=0
    do ipos=1,len(string)
       ilett=ichar(string(ipos:ipos))
       if (ilett.ne.ispc) then
          icount=icount+1
          utility_strip(icount:icount)=string(ipos:ipos)
       endif
    enddo

    utility_strip=trim(utility_strip)

    return

  end function utility_strip


  !=================================!
  function utility_lowercase(string)!
  !=================================!
  !                                 !
  ! Takes a string and converts to  !
  !      lowercase characters       !
  !                                 !
  !=================================!

    use w90_io, only : maxlen

    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen) :: utility_lowercase

    integer :: iA,iZ,idiff,ipos,ilett

    iA = ichar('A')
    iZ = ichar('Z')
    idiff = iZ-ichar('z')

    utility_lowercase = string

    do ipos=1,len(string)
       ilett = ichar(string(ipos:ipos))
       if ((ilett.ge.iA).and.(ilett.le.iZ)) &
            utility_lowercase(ipos:ipos)=char(ilett-idiff)
    enddo

    utility_lowercase = trim(adjustl(utility_lowercase))

    return

  end function utility_lowercase



  !====================================================!
  subroutine utility_string_to_coord(string_tmp,outvec)!
  !====================================================!
  !                                                    !
  !      Takes a string in the form 0.0,1.0,0.5        !
  !       and returns an array of the real num         !
  !                                                    !
  !====================================================!
    use w90_io, only :io_error,maxlen

    implicit none

    character(len=maxlen), intent(in)  :: string_tmp
    real(kind=dp), intent(out) :: outvec(3)
    
    integer :: pos
    character(len=maxlen)  :: ctemp
    character(len=maxlen)  :: ctemp2


    ctemp=string_tmp
    pos=index(ctemp,',')
    if(pos<=0) call io_error('utility_string_to_coord: Problem reading string into real number '//trim(string_tmp))
    ctemp2=ctemp(1:pos-1)
    read(ctemp2,*,err=100,end=100) outvec(1)
    ctemp=ctemp(pos+1:)
    pos=index(ctemp,',')
    ctemp2=ctemp(1:pos-1)
    read(ctemp2,*,err=100,end=100) outvec(2)
    ctemp=ctemp(pos+1:)
    read(ctemp,*,err=100,end=100) outvec(3)

    return

100 call io_error('utility_string_to_coord: Problem reading string into real number '//trim(string_tmp))

  end subroutine  utility_string_to_coord


!!$  !===========================================!
!!$  function utility_string_to_coord(string_tmp)!
!!$  !===========================================!
!!$  !                                           !
!!$  !  Takes a string in the form 0.0,1.0,0.5   !
!!$  !   and returns an array of the real num    !
!!$  !                                           !
!!$  !===========================================!
!!$
!!$    implicit none
!!$
!!$    character(len=80), intent(in)  :: string_tmp
!!$    real(kind=dp) :: utility_string_to_coord(3)
!!$    
!!$    integer :: pos,pos2
!!$    character(len=80)  :: ctemp
!!$    character(len=80)  :: ctemp2
!!$
!!$
!!$    ctemp=string_tmp
!!$    pos2=index(ctemp,',')
!!$    ctemp2=ctemp(1:pos2-1)
!!$    read(ctemp2,*) utility_string_to_coord(1)
!!$    ctemp=ctemp(pos2+1:)
!!$    pos2=index(ctemp,',')
!!$    ctemp2=ctemp(1:pos2-1)
!!$    read(ctemp2,*) utility_string_to_coord(2)
!!$    ctemp=ctemp(pos2+1:)
!!$    read(ctemp,*) utility_string_to_coord(3)
!!$
!!$  end function utility_string_to_coord

    !========================================================!
    subroutine utility_translate_home(vec,real_lat,recip_lat)
    !========================================================!
    !                                                        !
    !        Translate a vector to the home unit cell        !
    !                                                        !
    !========================================================!

      implicit none

      real(kind=dp), intent(inout) :: vec(3)
      real(kind=dp), intent(in)    :: real_lat(3,3)
      real(kind=dp), intent(in)    :: recip_lat(3,3)
      
      ! <<<local variables>>>
      integer       :: ind
      real(kind=dp) :: r_home(3),r_frac(3)
      real(kind=dp) :: shift

      r_home=0.0_dp;r_frac=0.0_dp

      ! Cartesian --> fractional
      call utility_cart_to_frac(vec,r_frac,recip_lat)
      ! Rationalise to interval [0,1]
      do ind=1,3
         if (r_frac(ind).lt.0.0_dp) then
            shift=real(ceiling(abs(r_frac(ind))),kind=dp)
            r_frac(ind)=r_frac(ind)+shift
         endif
         if (r_frac(ind).gt.1.0_dp) then
            shift=-real(int(r_frac(ind)),kind=dp)
            r_frac(ind)=r_frac(ind)+shift
         endif
      enddo
      ! Fractional --> Cartesian
      call utility_frac_to_cart(r_frac,r_home,real_lat)
      
      vec = r_home

      return
    end subroutine utility_translate_home

  !============================================================!
  subroutine utility_diagonalize(mat,dim,eig,rot)
  !============================================================!
  !                                                            !
  ! Diagonalize the dim x dim  hermitian matrix 'mat' and      !
  ! return the eigenvalues 'eig' and the unitary rotation 'rot'!
  !                                                            !
  !============================================================!

    use w90_constants, only : dp,cmplx_0
    use w90_io, only        : io_error,stdout
   
    integer, intent(in)           :: dim
    complex(kind=dp), intent(in)  :: mat(dim,dim)
    real(kind=dp), intent(out)    :: eig(dim)
    complex(kind=dp), intent(out) :: rot(dim,dim)

    complex(kind=dp)   :: mat_pack((dim*(dim+1))/2),cwork(2*dim)
    real(kind=dp)      :: rwork(7*dim)
    integer            :: i,j,info,nfound,iwork(5*dim),ifail(dim)

    do j=1,dim
       do i=1,j
          mat_pack(i+((j-1)*j)/2)=mat(i,j)
       enddo
    enddo
    rot=cmplx_0;eig=0.0_dp;cwork=cmplx_0;rwork=0.0_dp;iwork=0
    call ZHPEVX('V','A','U',dim,mat_pack,0.0_dp,0.0_dp,0,0,-1.0_dp, &
         nfound,eig(1),rot,dim,cwork,rwork,iwork,ifail,info)
    if(info < 0) then
       write(stdout,'(a,i3,a)') 'THE ',-info,&
            ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
       call io_error('Error in utility_diagonalize')
    endif
    if(info > 0) then
       write(stdout,'(i3,a)') info,' EIGENVECTORS FAILED TO CONVERGE'
       call io_error('Error in utility_diagonalize')
    endif

  end subroutine utility_diagonalize
 
  !===========================================================!
  function utility_rotate(mat,rot,dim)
  !==========================================================!
  !                                                           !
  ! Rotates the dim x dim matrix 'mat' according to           !
  ! (rot)^dagger.mat.rot, where 'rot' is a unitary matrix     !
  !                                                           !
  !===========================================================!

    use w90_constants, only : dp
   
    integer          :: dim
    complex(kind=dp) :: utility_rotate(dim,dim)
    complex(kind=dp) :: mat(dim,dim)
    complex(kind=dp) :: rot(dim,dim)

    utility_rotate=matmul(matmul(transpose(conjg(rot)),mat),rot)

  end function utility_rotate
 
  !===========================================================!
  function utility_matmul_diag(mat1,mat2,dim)
  !===========================================================!
  !                                                           !
  ! Computes the diagonal elements of the matrix mat1.mat2    !
  !                                                           !
  !===========================================================!

    use w90_constants, only : dp,cmplx_0
   
    integer          :: dim
    complex(kind=dp) :: utility_matmul_diag(dim)
    complex(kind=dp) :: mat1(dim,dim)
    complex(kind=dp) :: mat2(dim,dim)

    integer i,j

    utility_matmul_diag=cmplx_0
   do i=1,dim
       do j=1,dim
          utility_matmul_diag(i)=utility_matmul_diag(i)+mat1(i,j)*mat2(j,i)
       end do
    end do

    end function utility_matmul_diag
 
 
  !===========================================================!
  function utility_rotate_diag(mat,rot,dim)
  !===========================================================!
  !                                                           !
  ! Rotates the dim x dim matrix 'mat' according to           !
  ! (rot)^dagger.mat.rot, where 'rot' is a unitary matrix.    !
  ! Computes only the diagonal elements of rotated matrix.    !
  !                                                           !
  !===========================================================!

    use w90_constants, only : dp
   
    integer          :: dim
    complex(kind=dp) :: utility_rotate_diag(dim)
    complex(kind=dp) :: mat(dim,dim)
    complex(kind=dp) :: rot(dim,dim)


    utility_rotate_diag=utility_matmul_diag(matmul(transpose(conjg(rot)),mat),rot,dim)

  end function utility_rotate_diag
 
  !===========================================================!
  function utility_commutator_diag(mat1,mat2,dim)
  !===========================================================!
  !                                                           !
  ! Computes diagonal elements of                             !
  ! [mat1,mat2]=mat1.mat2-mat2.mat1                           ! 
  !                                                           !
  !===========================================================!

    use w90_constants, only : dp
   
    integer          :: dim
    complex(kind=dp) :: utility_commutator_diag(dim)
    complex(kind=dp) :: mat1(dim,dim)
    complex(kind=dp) :: mat2(dim,dim)

    utility_commutator_diag=utility_matmul_diag(mat1,mat2,dim)-utility_matmul_diag(mat2,mat1,dim)

  end function utility_commutator_diag

 
  !===================================================!
  function utility_re_tr(mat)
  !========================!
  !                        !
  ! Real part of the trace !
  !                        !
  !========================!

    use w90_constants, only  : dp,cmplx_0,cmplx_i
    
    real(kind=dp)                    :: utility_re_tr
    complex(kind=dp), dimension(:,:) :: mat
    
    integer          :: i,mydim
    complex(kind=dp) :: cdum

    mydim=size(mat,1)
    
    cdum=cmplx_0
    do i=1,mydim
       cdum=cdum+mat(i,i)
    enddo
    utility_re_tr=aimag(cmplx_i*cdum)
    
  end function utility_re_tr

  function utility_im_tr(mat)
  !=============================!
  !                             !
  ! Imaginary part of the trace !
  !                             !
  !=============================!

    use w90_constants, only  : dp,cmplx_0

    real(kind=dp)                    :: utility_im_tr
    complex(kind=dp), dimension(:,:) :: mat
    
    integer          :: i,mydim
    complex(kind=dp) :: cdum
    
    mydim=size(mat,1)

    cdum=cmplx_0
    do i=1,mydim
       cdum=cdum+mat(i,i)
    enddo
    utility_im_tr=aimag(cdum)

  end function utility_im_tr


  function w0gauss (x)!, n)
!-----------------------------------------------------------------------
!
!     the derivative of wgauss:  an approximation to the delta function
!
! --> (n>=0) : derivative of the corresponding Methfessel-Paxton wgauss
!
! --> (n=-1 ): derivative of cold smearing:
!              1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
!
! --> (n=-99): derivative of Fermi-Dira! function: 0.5/(1.0+cosh(x))
!
                                                                                                                                            
      implicit none
      real(kind=dp) :: w0gauss, x
! output: the value of the function
! input: the point where to compute the function
                                                                                                                                            
      integer :: n
! input: the order of the smearing function
!
!    here the local variables
!
      real(kind=dp) :: a, arg, hp, hd, pi
! the coefficients a_n
! the argument of the exponential
! the hermite function
! the hermite function
! pi
                                                                                                                                            
      integer :: i, ni
! counter on n values
! counter on 2n values
n=1
      pi = 3.14159265358979d0
! Fermi-Dira! smearing
!      if (n.eq. - 99) then
!         if (abs (x) .le.36.0) then
!            w0gauss = 1.0d0 / (2.0 + exp ( - x) + exp ( + x) )
!! in order to avoid problems for large values of x in the e
!         else
!            w0gauss = 0.d0
!         endif
!         return
                                                                                                                                            
!      endif
!     cold smearing  (Marzari-Vanderbilt)
!      if (n.eq. - 1) then
!         arg = min (200.d0, (x - 1.0d0 / sqrt (2.0d0) ) **2)
!         w0gauss = 1.0d0 / sqrt (pi) * exp ( - arg) * (2.0d0 - sqrt (2.0d0) * x)
!         return
                                                                                                                                            
!      endif
!     Methfessel-Paxton
      arg = min (200.d0, x**2)
      w0gauss = exp ( - arg) / sqrt (pi)
      if (n.eq.0) return
      hd = 0.0d0
      hp = exp ( - arg)
      ni = 0
      a = 1.0 / sqrt (pi)
      do i = 1, n
         hd = 2.0d0 * x * hp - 2.0d0 * dble (ni) * hd
         ni = ni + 1
         a = - a / (dble (i) * 4.0d0)
         hp = 2.0d0 * x * hd-2.0d0 * dble (ni) * hp
         ni = ni + 1
         w0gauss = w0gauss + a * hp
      enddo
      return
      end function w0gauss




! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Publi! License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
      function wgauss (x, n)
!-----------------------------------------------------------------------
!
!     this function computes the approximate theta function for the
!     given order n, at the point x.
!
! --> (n>=0) : Methfessel-Paxton case. See PRB 40, 3616 (1989).
!
! --> (n=-1 ): Cold smearing (Marzari-Vanderbilt). See PRL 82, 3296 (199
!       1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
!
! --> (n=-99): Fermi-Dira! case: 1.0/(1.0+exp(-x)).
!
use w90_io
      implicit none

      real(kind=dp) :: wgauss, x
! output: the value of the function
! input: the argument of the function
      integer :: n
! input: the order of the function
!
!    the local variables
!

      real(kind=dp) :: a, hp, arg, maxarg, hd, pi, xp
! the coefficient a_n
! the hermitean function
! the argument of the exponential
! maximum value for the argument of the exponen
! the hermitean function
! pi
! the freq function
! the erf function
! auxiliary variable (cold smearing)
      integer :: i, ni
! counter on the n indices
! counter on 2n
      real(kind=dp), parameter :: invsqrtpi=.39894228040143288356_dp
      real(kind=dp), parameter :: invsqrttwo =0.70710678118654752440_dp
      parameter (maxarg = 200.d0)
      
      pi = 3.14159265358979d0
call io_stopwatch('wgauss',1)
      ! Fermi-Dira! smearing
      if (n.eq. - 99) then
         if (x.lt. - maxarg) then
            wgauss = 0.d0
         elseif (x.gt.maxarg) then
            wgauss = 1.d0
         else
            wgauss = 1.0d0 / (1.0 + exp ( - x) )
         endif
         return
         
      endif
      ! Cold smearing
      if (n.eq. - 1) then
         xp = x - invsqrttwo!1.0d0 / sqrt (2.0d0)
!         arg = min (maxarg, xp**2)
!         write(51,*) x,xp
         wgauss = 0.5d0 * algor_erf (xp) + invsqrtpi *  exp ( - xp*xp) + 0.5d0
        
!         wgauss = 0.5d0 * erf (xp) + invsqrtpi *  exp ( - xp*xp) + 0.5d0
!         wgauss = 0.5d0 * erf (xp) + 1.0d0 / sqrt (2.0d0 * pi) *  exp ( - arg) + 0.5d0
call io_stopwatch('wgauss',2)
         return
      endif
! Methfessel-Paxton
!      wgauss = gauss_freq (x * sqrt (2.0d0) )
      wgauss = 0.5d0*(1.d0-algor_erf(-1.d0*x))
      if (n.eq.0) return
      hd = 0.d0
      arg = min (maxarg, x**2)
      hp = exp ( - arg)
      ni = 0
      a = 1.d0 / sqrt (pi)
      do i = 1, n
         hd = 2.0d0 * x * hp - 2.0d0 * dble (ni) * hd
         ni = ni + 1
         a = - a / (dble (i) * 4.0d0)
         wgauss = wgauss - a * hd
         hp = 2.0d0 * x * hd-2.0d0 * dble (ni) * hp
         ni = ni + 1
      enddo
call io_stopwatch('wgauss',2)
      return
      end function wgauss



 

  function algor_erf(x)
    !=========================================================================!
    ! Calculate an accurage approximation to the error function, erf(x),      !
    !    erf(x)=(2/sqrt(pi))*integral(0->x)[exp(-t^2)]dt.                     !
    !    Based upon parameterization given in NSWC Mathematics Library.       !
    !    NB This is machine portable and not dependent on a system function.  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x, intent=in, argument of erf to evaluate.                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 25/09/2001                               !
    !=========================================================================!
    implicit none
                                                                                                                                            
    real(kind=dp), intent(in) :: x
    real(kind=dp)             :: algor_erf
                                                                                                                                            
    !expansion parameters
    real(kind=dp), dimension(1:5) :: a=(/ 0.771058495001320E-04_dp, -0.133733772997339E-02_dp,   &
         & 0.323076579225834E-01_dp, 0.479137145607681E-01_dp, 0.128379167095513E+00_dp /)
    real(kind=dp), dimension(1:3) :: b=(/ 0.301048631703895E-02_dp,  0.538971687740286E-01_dp,   &
         & 0.375795757275549E+00_dp /)
    real(kind=dp), dimension(1:8) :: p=(/-1.36864857382717E-07_dp,  5.64195517478974E-01_dp, &
         & 7.21175825088309E+00_dp,  4.31622272220567E+01_dp,  1.52989285046940E+02_dp, &
         & 3.39320816734344E+02_dp,  4.51918953711873E+02_dp,  3.00459261020162E+02_dp /)
    real(kind=dp), dimension(1:8) :: q=(/ 1.00000000000000E+00_dp,  1.27827273196294E+01_dp, &
         & 7.70001529352295E+01_dp,  2.77585444743988E+02_dp,  6.38980264465631E+02_dp,  &
         & 9.31354094850610E+02_dp,  7.90950925327898E+02_dp,  3.00459260956983E+02_dp /)
    real(kind=dp), dimension(1:5) :: r=(/ 2.10144126479064E+00_dp,  2.62370141675169E+01_dp, &
         & 2.13688200555087E+01_dp,  4.65807828718470E+00_dp,  2.82094791773523E-01_dp /)
    real(kind=dp), dimension(1:4) :: s=(/ 9.41537750555460E+01_dp,  1.87114811799590E+02_dp, &
         & 9.90191814623914E+01_dp,  1.80124575948747E+01_dp /)
    real(kind=dp)                 :: c=0.564189583547756_dp

    !local vars
    real(kind=dp) :: ax, bot, t, top, x2
                                                                                                                                            
    ax=abs(x)
    x2=x*x
    if (ax<0.5_dp) then
       t=x2
       top=((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0_dp
       bot=((  b(1)*t + b(2))*t + b(3))*t + 1.0_dp
       algor_erf=ax*(top/bot)
                                                                                                                                            
    else if (ax<4.0_dp) then
       top=((((((p(1) *ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
            &  + p(6))*ax + p(7))*ax + p(8)
       bot=((((((q(1) *ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
            &  + q(6))*ax + q(7))*ax + q(8)
       algor_erf=0.5_dp + (0.5_dp - exp(-x2)*top/bot)
                                                                                                                                            
    else if (ax<5.8_dp) then
       t=1.0_dp/x2
       top=(((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
       bot=(((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0_dp
       algor_erf=0.5_dp + (0.5_dp - exp(-x2)*((c - top/(x2*bot))/ax))
                                                                                                                                            
    else
       algor_erf=1.0_dp            !large |x| limit
                                                                                                                                            
    end if
                                                                                                                                            
    !now put the sign back in
    algor_erf=sign(algor_erf,x)


    return
  end function algor_erf
                                                                                                                                            
  function algor_erfc(x)
    !=========================================================================!
    ! Calculate an accurage approximation to the complemenary error function, !
    !    erfc(x)=(2/sqrt(pi))*integral(x->infinity)[exp(-t^2)]dt.             !
    !    Based upon parameterization given in NSWC Mathematics Library.       !
    !    NB This is machine portable and not dependent on a system function.  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x, intent=in, argument of erfc to evaluate.                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 25/09/2001                               !
    !=========================================================================!
    implicit none
                                                                                                                                            
    real(kind=dp), intent(in) :: x
    real(kind=dp)             :: algor_erfc

    !expansion parameters
    real(kind=dp), dimension(1:5) :: a=(/  0.771058495001320E-04_dp, -0.133733772997339E-02_dp,  &
         &  0.323076579225834E-01_dp,  0.479137145607681E-01_dp,  0.128379167095513E+00_dp /)
    real(kind=dp), dimension(1:3) :: b=(/  0.301048631703895E-02_dp,  0.538971687740286E-01_dp,  &
         &  0.375795757275549E+00_dp /)
    real(kind=dp), dimension(1:8) :: p=(/ -1.36864857382717E-07_dp,  5.64195517478974E-01_dp, &
         &  7.21175825088309E+00_dp,  4.31622272220567E+01_dp,  1.52989285046940E+02_dp,  &
         &  3.39320816734344E+02_dp,  4.51918953711873E+02_dp,  3.00459261020162E+02_dp /)
    real(kind=dp), dimension(1:8) :: q=(/  1.00000000000000E+00_dp,  1.27827273196294E+01_dp, &
         &  7.70001529352295E+01_dp, 2.77585444743988E+02_dp,  6.38980264465631E+02_dp, &
         &  9.31354094850610E+02_dp, 7.90950925327898E+02_dp,  3.00459260956983E+02_dp /)
    real(kind=dp), dimension(1:5) :: r=(/  2.10144126479064E+00_dp,  2.62370141675169E+01_dp, &
         &  2.13688200555087E+01_dp, 4.65807828718470E+00_dp,  2.82094791773523E-01_dp /)
    real(kind=dp), dimension(1:4) :: s=(/ 9.41537750555460E+01_dp, 1.87114811799590E+02_dp, &
         &  9.90191814623914E+01_dp, 1.80124575948747E+01_dp /)
    real(kind=dp)                 :: c = .564189583547756_dp
                                                                                                                                            
    ! Local variables
    real(kind=dp) :: ax, x2, t, top, bot
                                                                                                                                            
    ax=abs(x)
    x2=x*x
    if (ax<0.5_dp) then
       t=x2
       top=((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0_dp
       bot=((  b(1)*t + b(2))*t + b(3)) * t + 1.0_dp
       algor_erfc=0.5_dp + (0.5_dp-x*(top/bot))
                                                                                                                                            
    else if (ax<4.0_dp) then
       top=((((((p(1) *ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax  &
            &  + p(7))*ax + p(8)
       bot=((((((q(1) *ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax  &
            &  + q(7))*ax + q(8)
       algor_erfc=exp(-x2)*top/bot
       if (x<0.0_dp) algor_erfc=2.0_dp - algor_erfc
                                                                                                                                            
    else
       if (x<=-5.6_dp) then
          algor_erfc=2.0_dp    !large negative x limit

       else if (x>100.0_dp) then
          algor_erfc=0.0_dp    !large positive x limit
                                                                                                                                            
       else
          t=1.0_dp/x2
          top=(((r(1)*t + r(2))*t + r(3))*t + r(4)) * t + r(5)
          bot=(((s(1)*t + s(2))*t + s(3))*t + s(4)) * t + 1.0_dp
          algor_erfc=exp(-x2)*(c-t*top/bot)/ax
          if (x<0.0_dp) algor_erfc=2.0_dp - algor_erfc
                                                                                                                                            
       end if
    end if
                                                                                                                                            
    return
  end function algor_erfc
                       



end module w90_utility
