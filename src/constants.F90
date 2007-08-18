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

module w90_constants

  implicit none

  private

  integer, parameter, public          :: dp = selected_real_kind(15,300)
  real(kind=dp), parameter, public    :: pi=3.141592653589793238462643383279_dp
  real(kind=dp), parameter, public    :: twopi = 2*pi
  complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)
  complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)

  real(kind=dp), parameter, public :: bohr = 0.5291772108_dp

  real(kind=dp), parameter, public    :: eps2  = 1.0e-2_dp
  real(kind=dp), parameter, public    :: eps5  = 1.0e-5_dp
  real(kind=dp), parameter, public    :: eps6  = 1.0e-6_dp
  real(kind=dp), parameter, public    :: eps7  = 1.0e-7_dp
  real(kind=dp), parameter, public    :: eps8  = 1.0e-8_dp
  real(kind=dp), parameter, public    :: eps10 = 1.0e-10_dp

end module w90_constants
