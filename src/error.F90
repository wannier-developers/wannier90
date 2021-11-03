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

module w90_error
  !! Module to handle errors

  use w90_constants, only: dp
  !use w90_types, only: w90_error_type

  implicit none

  public

  integer, parameter :: code_ok = 0 
  integer, parameter :: code_fatal = 1 
  integer, parameter :: code_alloc = 2 
  integer, parameter :: code_dealloc = 3 
  integer, parameter :: code_unconv = -1

  type w90_error_type
    !! Codify error state with integer code and human readable string
    integer :: code
    character(len=120) :: message
    contains
    final :: untrapped_error
  end type w90_error_type

  contains

  subroutine untrapped_error(err)
    type(w90_error_type), intent(in) :: err
    ! this routine should never be called, writing to stderr in desparation
    write(0,*) "UNTRAPPED ERROR: ", err%code
    write(0,*) "UNTRAPPED ERROR: ", err%message
    ! this routine should never be called, call stop() in utter desparation FIXME?
    stop
  end subroutine untrapped_error

  subroutine set_error_alloc(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate(err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_alloc
  end subroutine set_error_alloc

  subroutine set_error_fatal(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate(err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_fatal
  end subroutine set_error_fatal

  subroutine set_error_dealloc(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate(err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_dealloc
  end subroutine set_error_dealloc

  subroutine set_error_unconv(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate(err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_unconv
  end subroutine set_error_unconv




end module w90_error

  


