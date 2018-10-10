!-*- mode: F90 -*-!
!                                                            !
! Copyright (C) 2012 Lampros Andrinopoulos, Nicholas Hine,   !
!                    and Arash A Mostofi                     !
!                                                            !
! In publications arising from use of this code, please cite !
!                                                            !
!     L Andrinopoulos, NDM Hine and AA Mostofi,              !
!                  J. Chem. Phys. 135, 154105 (2011)         !
!                                                            !
! The methodology coded here is based on the original work   !
! of PL Silvestrelli, Phys Rev Lett 100, 053002 (2008)       !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!
! Code to calculate vdW energies from MLWFs                  !
!------------------------------------------------------------!
program vdw

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: k = 1500         ! for integration scheme
  real(dp) :: fnl, fnl2, fnl4, c6nl, evdw, rnl, f, c6_eff, shift, &
              rc_coeff, degeneracy, tol_dist, tol_occ
  integer :: n, l, num_frag, iwann, iwann2, ifrag, ifrag2
  real(dp), allocatable :: r(:, :, :), rc(:, :), s(:, :), rvdw(:, :), fw(:, :), &
                           num_proxim(:, :), distance(:, :, :)
!!  real(dp), allocatable :: s_amalg(:,:),r_amalg(:,:)
  real(dp), parameter :: pi = 3.141592653589793238462643383279d0
  real(dp), parameter :: ang2a0 = 0.52917720859d0 ! NIST CODATA 2006 value
!!  real(dp), parameter :: ang2a0=0.5291772108d0 ! old W90 value
  character(len=30) :: filename, units, dummy
  integer, allocatable :: num_wann(:), num_wann_new(:), num_orb(:)
  integer :: num_wann_old, num_wann_tot
  integer :: max_num_wann
  logical, allocatable :: p(:, :)
  logical :: amalgamate, disentangle
  real(kind=DP) :: rnp
  integer :: iostat

  call get_command_argument(1, filename)

  open (unit=1, file=filename, action='read', iostat=iostat)
  if (iostat /= 0) then
    write (0, '(a)') "Error while reading file '"//trim(filename)//"'."
    write (0, '(a)') "Pass the correct file name on the command line!"
    stop 1
  end if
  open (unit=2, file=adjustl(trim(filename))//".out", action='write')
  write (2, '(a/)') ' Running w90_vdw. Written by L Andrinopoulos et al. (c) 2012.'
  write (2, '(a/)') ' L Andrinopoulos, NDM Hine and AA Mostofi, J Chem Phys 135, 154105 (2011).'
  write (2, '(a,a,a/)') ' Reading input file ', adjustl(trim(filename)), '...'

  ! read input parameters
  read (1, *) dummy, disentangle ! whether MLWFs are disentangled
  write (2, '(a,l2)') ' disentangle: ', disentangle
  read (1, *) dummy, amalgamate  ! whether MLWFs need to be amalgamated
  write (2, '(a,l2)') ' amalgamate : ', amalgamate
  read (1, *) dummy, degeneracy  ! number of electrons per MLWF
  write (2, '(a,f5.3)') ' degeneracy : ', degeneracy
  read (1, '(9x,i3)') num_frag   ! number of distinct molecular fragments
  write (2, '(a,i3)') ' num_frag   : ', num_frag
  allocate (num_wann(num_frag), num_wann_new(num_frag))
  allocate (p(num_frag, 3), num_orb(num_frag))
  read (1, '(12x)')
  read (1, *) num_wann(1:num_frag) ! number of MLWFs per fragment
  write (2, '(a)', advance='no') ' num_wann   : '
  do n = 1, num_frag
    write (2, '(i4,1x)', advance='no') num_wann(n)
  enddo
  read (1, *) dummy, tol_occ ! tolerance for splitting p-orbitals
  write (2, '(/a,f6.3)') ' tol_occ    : ', tol_occ
  if (.not. disentangle) write (2, '(a)') ' [disentanglement not used; ignoring tol_occ]'
  read (1, '(12x)')
  do ifrag = 1, num_frag
    read (1, '(l1,1x,l1,1x,l1)') p(ifrag, 1:3) ! directions in which to split p-orbitals
  enddo
  do ifrag = 1, num_frag
    write (2, '(a,i2,1x,l2,1x,l2,1x,l2)') ' pxyz       : ', ifrag, (p(ifrag, n), n=1, 3)
  enddo
  if (.not. disentangle) write (2, '(a)') ' [disentanglement not used; ignoring pxyz]'
  ! tolerance for amalgamating cocentric MLWFs
  ! not used unless amalgamate = T
  read (1, *) dummy, tol_dist
  write (2, '(a,f6.3)') ' tol_dist   : ', tol_dist
  if (.not. amalgamate) write (2, '(a)') ' [amalgamate not used; ignoring tol_dist]'

  max_num_wann = maxval(num_wann(:))
  allocate (r(num_frag, 3, 2*max_num_wann), &
            rc(num_frag, 2*max_num_wann), &
            s(num_frag, 2*max_num_wann), &
            rvdw(num_frag, 2*max_num_wann), &
            fw(num_frag, 2*max_num_wann), &
            distance(num_frag, max_num_wann, max_num_wann), &
            num_proxim(num_frag, max_num_wann))

  read (1, '(19x)')
  read (1, *) units ! units in which centres and spreads are given
  if (adjustl(trim(units)) == 'ang') then
    write (2, '(a,a)') ' units      :  ', adjustl(trim(units))
  else
    write (2, '(a)') ' units      :   bohr (assumed)'
  endif
  do ifrag = 1, num_frag
    do iwann = 1, num_wann(ifrag)
      read (1, *) r(ifrag, 1:3, iwann), s(ifrag, iwann), fw(ifrag, iwann)
    enddo
  enddo
  close (1)
  ! end read input parameters

  ! set electronic occupancy of each MLWF
  fw(:, :) = degeneracy*fw(:, :)
  ! convert spread^2 to spread
  s(:, :) = sqrt(s(:, :))

  ! aam: convert centres and spreads to bohr if necessary
  !      (code works in atomic units)
  if (adjustl(trim(units)) == 'ang') then
    r(:, :, :) = r(:, :, :)/ang2a0
    s(:, :) = s(:, :)/ang2a0
  endif

  write (2, '(a)') ' centres, sqrt(quadratic spreads) [in Ang], and electronic occupancies:'
  do ifrag = 1, num_frag
    do iwann = 1, num_wann(ifrag)
      write (2, '(i3,1x,i3,1x,3f13.8,f12.8,f12.8)') &
        ifrag, iwann, r(ifrag, 1:3, iwann)*ang2a0, s(ifrag, iwann)*ang2a0, fw(ifrag, iwann)
    enddo
  enddo

  write (2, '(/a)') ' done reading input file'

  if (amalgamate) then

    write (2, '(/a)', advance='no') ' amalgamating MLWF centres...'

    num_proxim(:, :) = 0
    do ifrag = 1, num_frag

      do
        num_wann_old = num_wann(ifrag)

        iwannloop: do iwann = 1, num_wann_old
          do iwann2 = iwann + 1, num_wann_old

            distance(ifrag, iwann, iwann2) = sqrt( &
                                             (r(ifrag, 1, iwann) - r(ifrag, 1, iwann2))**(2.0d0) + &
                                             (r(ifrag, 2, iwann) - r(ifrag, 2, iwann2))**(2.0d0) + &
                                             (r(ifrag, 3, iwann) - r(ifrag, 3, iwann2))**(2.0d0))

            if (distance(ifrag, iwann, iwann2) <= tol_dist) then
              num_proxim(ifrag, iwann) = num_proxim(ifrag, iwann) + 1
              num_wann(ifrag) = num_wann(ifrag) - 1
              rnp = real(num_proxim(ifrag, iwann), DP)
              r(ifrag, 1:3, iwann) = r(ifrag, 1:3, iwann)*rnp/(rnp + 1.0_DP) + r(ifrag, 1:3, iwann2)/(rnp + 1.0_DP)
              r(ifrag, :, iwann2:num_wann(ifrag)) = r(ifrag, :, iwann2 + 1:num_wann(ifrag) + 1)
              s(ifrag, iwann) = s(ifrag, iwann)*rnp/(rnp + 1.0_DP) + s(ifrag, iwann2)/(rnp + 1.0_DP)
              s(ifrag, iwann2:num_wann(ifrag)) = s(ifrag, iwann2 + 1:num_wann(ifrag) + 1)
              fw(ifrag, iwann) = fw(ifrag, iwann) + fw(ifrag, iwann2)
              fw(ifrag, iwann2:num_wann(ifrag)) = fw(ifrag, iwann2 + 1:num_wann(ifrag) + 1)
              exit iwannloop
            endif
          enddo
        enddo iwannloop
        if (num_wann(ifrag) == num_wann_old) exit
      enddo

    enddo

    write (2, '(a)') ' ... done'

  endif

  num_orb(:) = 0
  do ifrag = 1, num_frag
    do l = 1, 3
      if (p(ifrag, l)) then
        do iwann = 1, num_wann(ifrag)
          if (fw(ifrag, iwann) <= 2.0d0*tol_occ) then
            num_orb(ifrag) = num_orb(ifrag) + 1
            s(ifrag, iwann) = 7.0d0*sqrt(2.0d0)/16.0d0*s(ifrag, iwann)
            shift = sqrt(15.0d0)/7.0d0*s(ifrag, iwann)
            r(ifrag, l, iwann) = r(ifrag, l, iwann) + shift
            s(ifrag, num_wann(ifrag) + num_orb(ifrag)) = s(ifrag, iwann)
            r(ifrag, 1:3, num_wann(ifrag) + num_orb(ifrag)) = r(ifrag, 1:3, iwann)
            r(ifrag, l, num_wann(ifrag) + num_orb(ifrag)) = r(ifrag, l, iwann) - 2.0d0*shift
            fw(ifrag, iwann) = 0.5d0*fw(ifrag, iwann)
            fw(ifrag, num_wann(ifrag) + num_orb(ifrag)) = fw(ifrag, iwann)
          endif
        enddo
      endif
    enddo
    num_wann(ifrag) = num_wann(ifrag) + num_orb(ifrag)
  enddo

  num_wann_tot = 0
  do ifrag = 1, num_frag
    num_wann_tot = num_wann_tot + num_wann(ifrag)
  enddo

  !Cutoff radius
  rc_coeff = 1.0_dp
  rc(:, :) = rc_coeff*s(:, :)*sqrt(3.0d0)*(0.769d0 + 0.5d0*log(s(:, :)))

!VdW radii obtained by 0.01 electron density contour;
!if (radcut=='vdw') then
!rvdw(:,:)=rc_coeff*s(:,:)*(1.475d0-0.866d0*log(s(:,:)))
!else
  rvdw(:, :) = rc(:, :)
!endif

  evdw = 0.0_dp
  c6_eff = 0.0_dp

  write (2, '(/a)') ' calculating E_vdW and C6_effective ...'

  write (2, '(/a)') ' centres, sqrt(quadratic spreads) [in Ang], and electronic occupancies'
  write (2, '(a)') ' after p-splitting (if any) and amalgamation of co-centric MLWFs (if any):'

  do ifrag = 1, num_frag
    do iwann = 1, num_wann(ifrag)
      ! write centres, spreads and occupancies used in calculating E_vdW
      write (2, '(i3,1x,i3,1x,3f13.8,f12.8,f12.8)') &
        ifrag, iwann, r(ifrag, 1:3, iwann)*ang2a0, s(ifrag, iwann)*ang2a0, fw(ifrag, iwann)

      do ifrag2 = ifrag + 1, num_frag
        do iwann2 = 1, num_wann(ifrag2) !No double counting, no self-interacting terms

          ! |rn-rl|=rnl
          rnl = sqrt((r(ifrag, 1, iwann) - r(ifrag2, 1, iwann2))**(2.0d0) + &
                     (r(ifrag, 2, iwann) - r(ifrag2, 2, iwann2))**(2.0d0) + &
                     (r(ifrag, 3, iwann) - r(ifrag2, 3, iwann2))**(2.0d0))

          !Integration scheme (Romberg)
          call c(fnl, k, s(ifrag, iwann), s(ifrag2, iwann2), fw(ifrag, iwann), fw(ifrag2, iwann2), rc_coeff)
          call c(fnl2, 2*k, s(ifrag, iwann), s(ifrag2, iwann2), fw(ifrag, iwann), fw(ifrag2, iwann2), rc_coeff)
          call c(fnl4, 4*k, s(ifrag, iwann), s(ifrag2, iwann2), fw(ifrag, iwann), fw(ifrag2, iwann2), rc_coeff)
          f = 1.0d0/45.0d0*(64.0d0*fnl4 - 20.0d0*fnl2 + fnl)

          !C6 Coefficients and VdW Energy calculation
          c6nl = (s(ifrag, iwann)**(1.5d0))*(s(ifrag2, iwann2)**(3.0d0))*0.5d0*(3.0d0**(-1.25d0))*f
          c6_eff = c6_eff + c6nl
          evdw = evdw - 1.0d0/(1.0d0 + exp(-20.0d0*(rnl/(rvdw(ifrag, iwann) + rvdw(ifrag2, iwann2)) - 1))) &
                 *c6nl/(rnl**(6.0d0))

        enddo
      enddo
    enddo
  enddo

!  close(3)
  write (2, '(/a)') repeat('=', 50)
  write (2, '(a,f16.8,a)') ' vdW energy = ', evdw, ' Ha'
  write (2, '(a,f16.8,a)') ' C6_eff     = ', c6_eff, ' Ha*Bohr^6'
  write (2, '(a)') repeat('=', 50)

  write (2, '(/a)') ' Have a nice day.'

  close (2)

contains

  !Integration subroutine
  subroutine c(s, n, sn, sl, fn, fl, rcoeff)

    implicit none
    integer, intent(in) :: n
    real(dp) :: x, y, hx, hy, xc, yc, beta, integral, f
    integer :: i, j
    real(dp), allocatable :: a(:), b(:)
    real(dp), intent(in) :: sn, sl, fn, fl, rcoeff
    real(dp), intent(out) :: s

    beta = (sn/sl)**(1.5d0)

    !Integration boundaries
    xc = 3.0d0*rcoeff*(0.769d0 + 0.5d0*log(sn))
    yc = 3.0d0*rcoeff*(0.769d0 + 0.5d0*log(sl))

    allocate (a(n), b(n))

    !Sum coefficients
    a(:) = 1.0d0
    b(:) = 1.0d0
    a(1) = 13.0d0/12d0
    b(1) = a(1)
    a(n - 1) = a(1)
    b(n - 1) = a(1)
    a(n) = 5.0d0/12d0
    b(n) = a(n)

    hx = xc/n
    hy = yc/n

    integral = 0

    do i = 1, n
      x = i*hx
      do j = 1, n
        y = j*hy
        f = (x*x*y*y*exp(-x)*exp(-y))/(exp(-x)/(beta*sqrt(fl)) + exp(-y)/sqrt(fn))
        integral = integral + a(i)*b(j)*hx*hy*f
      end do
    end do

    s = integral
    deallocate (a, b)
  end subroutine c

end program vdw
