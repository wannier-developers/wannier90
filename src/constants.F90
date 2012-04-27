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

  !!! GENERIC CONSTANTS !!!
!aam_2012-04-11; fix to run on MacBook Air 
  integer, parameter, public          :: dp = kind(1.0d0)
!  integer, parameter, public          :: dp = selected_real_kind(14,200)
!  integer, parameter, public          :: dp = selected_real_kind(15,300)
  real(kind=dp), parameter, public    :: pi=3.141592653589793238462643383279_dp
  real(kind=dp), parameter, public    :: twopi = 2*pi
  complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)
  complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)

  !!! NUMERICAL CONVERGENCE CONSTANTS !!!
  real(kind=dp), parameter, public    :: eps2  = 1.0e-2_dp
  real(kind=dp), parameter, public    :: eps5  = 1.0e-5_dp
  real(kind=dp), parameter, public    :: eps6  = 1.0e-6_dp
  real(kind=dp), parameter, public    :: eps7  = 1.0e-7_dp
  real(kind=dp), parameter, public    :: eps8  = 1.0e-8_dp
  real(kind=dp), parameter, public    :: eps10 = 1.0e-10_dp
  ! Cutoff for the smearing functions
  real(kind=dp), parameter, public    :: smearing_cutoff = 10._dp 
  ! [GP, Apr 20, 2012] Don't smear but simply add the contribution to the
  ! relevant bin if the smearing/binwidth ratio is smaller than this value
  real(kind=dp), parameter, public    :: min_smearing_binwidth_ratio = 2._dp 


  !!! PHYSICAL CONSTANTS !!!
  real(kind=dp), parameter, public :: bohr = 0.5291772108_dp

  ! Added by Ivo, to be used in interpolation code
  !
  ! Values of the fundamental constants taken from 
  ! http://physics.nist.gov/cuu/Constants/index.html
  ! 
  real(kind=dp), parameter, public :: elem_charge_SI=1.602176487e-19_dp
  real(kind=dp), parameter, public :: elec_mass_SI=9.10938215e-31_dp
  real(kind=dp), parameter, public :: hbar_SI=1.054571628e-34_dp
  real(kind=dp), parameter, public :: k_B_SI=1.3806504e-23_dp
  real(kind=dp), parameter, public :: bohr_magn_SI=927.400915e-26_dp
  real(kind=dp), parameter, public :: eps0_SI=8.854187817e-12_dp
  real(kind=dp), parameter, public :: speedlight_SI=299792458.0_dp
  real(kind=dp), parameter, public :: eV_au=3.674932540e-2_dp

end module w90_constants
