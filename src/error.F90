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

  use w90_error_base
  use w90_comms

  implicit none

  public

  integer, parameter :: code_ok = 0

  integer, parameter :: code_fatal = 1
  integer, parameter :: code_alloc = 2
  integer, parameter :: code_dealloc = 3
  integer, parameter :: code_mpi = 4
  integer, parameter :: code_input = 5
  integer, parameter :: code_file = 6

  integer, parameter :: code_unconv = -1
  integer, parameter :: code_warning = -2 ! mostly for failing to plot something

contains

  subroutine set_error_fatal(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_fatal)
    call comms_sync_error(comm, err, code_fatal)
  end subroutine set_error_fatal

  subroutine set_error_alloc(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_alloc)
    call comms_sync_error(comm, err, code_alloc)
  end subroutine set_error_alloc

  subroutine set_error_dealloc(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_dealloc)
    call comms_sync_error(comm, err, code_dealloc)
  end subroutine set_error_dealloc

  ! note, this is not used in comms routines
  ! (otherwise circular dependency)
  subroutine set_error_mpi(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_mpi)
    call comms_sync_error(comm, err, code_mpi)
  end subroutine set_error_mpi

  subroutine set_error_input(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_input)
    call comms_sync_error(comm, err, code_input)
  end subroutine set_error_input

  subroutine set_error_file(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_file)
    call comms_sync_error(comm, err, code_file)
  end subroutine set_error_file

  subroutine set_error_unconv(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_unconv)
    call comms_sync_error(comm, err, code_unconv)
  end subroutine set_error_unconv

  subroutine set_error_warn(err, mesg, comm)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    type(w90_comm_type), intent(in) :: comm
    call set_base_error(err, mesg, code_warning)
    call comms_sync_error(comm, err, code_warning)
  end subroutine set_error_warn

end module w90_error
