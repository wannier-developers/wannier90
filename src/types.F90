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

module w90_types
  !! This module contains types common to wannier90/postw90
  !! Instances of these types may be defined in one or more parameters
  use w90_constants, only: dp

  public

  type projection_type
    integer, allocatable :: l(:)
    integer, allocatable :: m(:)
    integer, allocatable :: s(:)
    real(kind=dp), allocatable :: s_qaxis(:, :)
    real(kind=dp), allocatable :: z(:, :)
    real(kind=dp), allocatable :: x(:, :)
    integer, allocatable :: radial(:)
    real(kind=dp), allocatable :: zona(:)
  end type projection_type

end module w90_types
