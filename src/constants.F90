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

module w90_constants

  !! This module contains the definitions of constants
  !! used in Wannier90 - both numerical constants such as pi
  !! and numerical convergence tolerances, but also physical
  !! constant such as the speed of light
  !!
  !! Values of the fundamental constants are taken from
  !! http://physics.nist.gov/cuu/Constants/index.html
  !! By default CODATA2010 is used (CODATA2006 can be selected
  !! using an appropriate compile-time flag (see INSTALL guide)

  implicit none

  private

  !~~ GENERIC CONSTANTS ~~!
  integer, parameter, public          :: i64 = selected_int_kind(15)
  !! 64bit integer
!aam_2012-04-11; fix to run on MacBook Air
  integer, parameter, public          :: dp = kind(1.0d0)
  !! double precision
!  integer, parameter, public          :: dp = selected_real_kind(14,200)
!  integer, parameter, public          :: dp = selected_real_kind(15,300)
  real(kind=dp), parameter, public    :: pi = 3.141592653589793238462643383279_dp
  !! $$\pi$$
  real(kind=dp), parameter, public    :: twopi = 2*pi
  !! $$2\pi$$
  complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp, 1.0_dp)
  !! i as a complex variable
  complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp, 0.0_dp)
  !! 0 as a complex variable
  complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp, 0.0_dp)
  !! 1 as a complex variable

  !~~ NUMERICAL CONVERGENCE CONSTANTS ~~!
  real(kind=dp), parameter, public    :: eps2 = 1.0e-2_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps5 = 1.0e-5_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps6 = 1.0e-6_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps7 = 1.0e-7_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps8 = 1.0e-8_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps10 = 1.0e-10_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: smearing_cutoff = 10._dp
  !! Cutoff for the smearing functions
  real(kind=dp), parameter, public    :: min_smearing_binwidth_ratio = 2._dp
  !! Don't smear but simply add the contribution to the
  !! relevant bin if the smearing/binwidth ratio is smaller than this value

  !~~ PHYSICAL CONSTANTS ~~!
  !
  ! Values of the fundamental constants taken from
  ! http://physics.nist.gov/cuu/Constants/index.html
  !
!~ Pick up default value; unfortunately, it doesn't work with indentation...
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
  real(kind=dp), parameter, public :: elem_charge_SI = 1.602176487e-19_dp   ! C
  !! e
  real(kind=dp), parameter, public :: elec_mass_SI = 9.10938215e-31_dp      ! kg
  !! $$m_e$$
  real(kind=dp), parameter, public :: hbar_SI = 1.054571628e-34_dp          ! J * s
  !! $$\hbar$$
  real(kind=dp), parameter, public :: k_B_SI = 1.3806504e-23_dp             ! J / K
  !! $$k_B$$
  real(kind=dp), parameter, public :: bohr_magn_SI = 927.400915e-26_dp      ! J / T
  !! $$\mu_B$$
  real(kind=dp), parameter, public :: eps0_SI = 8.854187817e-12_dp          ! F / m
  !! $$\epsilon_0$$
  real(kind=dp), parameter, public :: speedlight_SI = 299792458.0_dp        ! m / s
  !! $$c$$
  real(kind=dp), parameter, public :: eV_au = 3.674932540e-2_dp              ! (see table of Conv. Factors)
  !! eV in atomic units
  real(kind=dp), parameter, public :: eV_seconds = 6.582119e-16_dp
  !! Electron Volt in seconds
  real(kind=dp), parameter, public :: bohr_angstrom_internal = 0.52917720859_dp
  !! Bohr to $$\AA$$
  ! Leave the length to this value, and don't exceed in length (needed for output formatting)
  character(len=75), parameter, public :: constants_version_str1 = "-> Using CODATA 2006 constant values"
  character(len=75), parameter, public :: constants_version_str2 = "   (http://physics.nist.gov/cuu/Constants/index.html)"
#endif

#ifdef CODATA2010
! ##### CODATA 2010 ##### !
!#warning "WANNIER90 INFO: Using CODATA 2010 constant values"
  real(kind=dp), parameter, public :: elem_charge_SI = 1.602176565e-19_dp
  !! elemental charge
  real(kind=dp), parameter, public :: elec_mass_SI = 9.10938291e-31_dp
  !! electron mass
  real(kind=dp), parameter, public :: hbar_SI = 1.054571726e-34_dp
  !! hbar
  real(kind=dp), parameter, public :: k_B_SI = 1.3806488e-23_dp
  !! Boltzman Constant
  real(kind=dp), parameter, public :: bohr_magn_SI = 927.400968e-26_dp
  !! Bohr magneton
  real(kind=dp), parameter, public :: eps0_SI = 8.854187817e-12_dp
  !! Vacuum Dielectric Constant
  real(kind=dp), parameter, public :: speedlight_SI = 299792458.0_dp
  !! Speed of light
  real(kind=dp), parameter, public :: eV_au = 3.674932379e-2_dp
  !! Electron Volt in atomic units
  real(kind=dp), parameter, public :: eV_seconds = 6.582119e-16_dp
  !! Electron Volt in seconds
  real(kind=dp), parameter, public :: bohr_angstrom_internal = 0.52917721092_dp
  !! Bohr to Anstrom Conversion factor
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
