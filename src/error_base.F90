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

module w90_error_base
  !! Module to handle errors

  implicit none

  public

  !! Codify error state with integer code and human readable string
  type w90_error_type
    integer :: code
    character(len=128) :: message
#ifdef W90DEV
  contains
    final :: untrapped_error
#endif
  end type w90_error_type

  integer, parameter :: code_remote = -99 ! special code for error triggered by other mpi rank
  integer, parameter :: code_deactivated = -888 ! special code for error triggered by other mpi rank

contains

  subroutine untrapped_error(err)
    type(w90_error_type), intent(in) :: err
    ! this routine should never be called, so write to stderr and call "stop" in desparation
    if (err%code == code_deactivated) return
    write (0, *) "UNTRAPPED ERROR: ", err%code
    write (0, *) "UNTRAPPED ERROR: ", err%message
    stop
  end subroutine untrapped_error

  subroutine set_base_error(err, mesg, code)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    integer, intent(in) :: code
    allocate (err)
    err%message = trim(mesg)
    err%code = code
  end subroutine set_base_error

  ! the following is required by subroutines of the comms module
  subroutine recv_error(err)
    type(w90_error_type), allocatable, intent(out) :: err
    allocate (err)
    err%code = code_remote
  end subroutine recv_error

end module w90_error_base
