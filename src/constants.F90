!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2007-13 Jonathan Yates, Arash Mostofi,       !
!                Giovanni Pizzi, Young-Su Lee,               !
!                Nicola Marzari, Ivo Souza, David Vanderbilt !
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
  !
  ! Values of the fundamental constants taken from 
  ! http://physics.nist.gov/cuu/Constants/index.html
  ! 
!! Pick up default value; unfortunately, it doesn't work with indentation...
#ifndef CODATA2006
#ifndef CODATA2010
#define CODATA2006
#endif
#endif

! We don't check for multiply-defined flags: if it is the case, the code will not
! compile because of multiple definitions of the same flags

#ifdef CODATA2006
! ##### CODATA 2006 ##### !
!#warning "WANNIER90 INFO: Using CODATA 2006 constant values"
  real(kind=dp), parameter, public :: elem_charge_SI=1.602176487e-19_dp   ! C
  real(kind=dp), parameter, public :: elec_mass_SI=9.10938215e-31_dp      ! kg
  real(kind=dp), parameter, public :: hbar_SI=1.054571628e-34_dp          ! J * s
  real(kind=dp), parameter, public :: k_B_SI=1.3806504e-23_dp             ! J / K
  real(kind=dp), parameter, public :: bohr_magn_SI=927.400915e-26_dp      ! J / T
  real(kind=dp), parameter, public :: eps0_SI=8.854187817e-12_dp          ! F / m
  real(kind=dp), parameter, public :: speedlight_SI=299792458.0_dp        ! m / s
  real(kind=dp), parameter, public :: eV_au=3.674932540e-2_dp              ! (see table of Conv. Factors)
  real(kind=dp), parameter, public :: bohr_angstrom_internal=0.52917720859_dp 
  ! Leave the length to this value, and don't exceed in length (needed for output formatting)
  character(len=75), parameter, public :: constants_version_str1 = "-> Using CODATA 2006 constant values"
  character(len=75), parameter, public :: constants_version_str2 = "   (http://physics.nist.gov/cuu/Constants/index.html)"
#endif

#ifdef CODATA2010
! ##### CODATA 2010 ##### !
!#warning "WANNIER90 INFO: Using CODATA 2010 constant values"
  real(kind=dp), parameter, public :: elem_charge_SI=1.602176565e-19_dp
  real(kind=dp), parameter, public :: elec_mass_SI=9.10938291e-31_dp
  real(kind=dp), parameter, public :: hbar_SI=1.054571726e-34_dp
  real(kind=dp), parameter, public :: k_B_SI=1.3806488e-23_dp
  real(kind=dp), parameter, public :: bohr_magn_SI=927.400968e-26_dp
  real(kind=dp), parameter, public :: eps0_SI=8.854187817e-12_dp
  real(kind=dp), parameter, public :: speedlight_SI=299792458.0_dp
  real(kind=dp), parameter, public :: eV_au=3.674932379e-2_dp
  real(kind=dp), parameter, public :: bohr_angstrom_internal=0.52917721092_dp
  ! Leave the length to this value, and don't exceed in length (needed for output formatting)
  character(len=75), parameter, public :: constants_version_str1 = "-> Using CODATA 2010 constant values"
  character(len=75), parameter, public :: constants_version_str2 = "   (http://physics.nist.gov/cuu/Constants/index.html)"
#endif

#ifdef USE_WANNIER90_V1_BOHR
!#warning "WANNIER90 INFO: Using WANNIER ver. 1 version of bohr"
  real(kind=dp), parameter, public :: bohr = 0.5291772108_dp
  ! Leave the length to this value, and don't exceed in length (needed for output formatting)
  character(len=75), parameter, public :: bohr_version_str = "-> Using Bohr value from Wannier90 ver. 1.x (DEPRECATED!)"
#else
!#warning "WANNIER90 INFO: Using CODATA version of bohr"
  real(kind=dp), parameter, public :: bohr = bohr_angstrom_internal
  ! Leave the length to this value, and don't exceed in length (needed for output formatting)
  character(len=75), parameter, public :: bohr_version_str = "-> Using Bohr value from CODATA"
#endif


end module w90_constants
