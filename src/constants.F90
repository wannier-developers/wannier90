!
! Copyright (C) 2004,2006 Jonathan Yates, Arash Mostofi,
!            Nicola Marzari, Ivo Souza, David Vanderbilt
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------


module constants

  implicit none

  private

  integer, parameter, public          :: dp = selected_real_kind(15,300)
  real(kind=dp), parameter, public    :: pi=3.141592653589793238462643383279_dp
  real(kind=dp), parameter, public    :: twopi = 2*pi
  complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)
  complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)

  real(kind=dp), parameter, public :: bohr = 0.5291772108_dp


end module constants
