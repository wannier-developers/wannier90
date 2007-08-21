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

end module w90_utility
