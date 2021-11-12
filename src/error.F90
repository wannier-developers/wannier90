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

  !use w90_constants, only: dp
  !use w90_types, only: w90_error_type

  implicit none

  public

  integer, parameter :: code_ok = 0
  integer, parameter :: code_fatal = 1
  integer, parameter :: code_alloc = 2
  integer, parameter :: code_dealloc = 3
  integer, parameter :: code_mpi = 4
  integer, parameter :: code_input_command = 5 !BGS is it a warning?
  integer, parameter :: code_file_read = 6
  integer, parameter :: code_file_open = 7
  integer, parameter :: code_matrix_lib = 8
  integer, parameter :: code_not_unitary = 9
  integer, parameter :: code_sitesym = 10 !BGS should probably just be fatal?
  integer, parameter :: code_unconv = -1
  integer, parameter :: code_plot = -2 ! failing to plot something isn't fatal?
  integer, parameter :: code_warning = -3
  integer, parameter :: code_disentangle = 11

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
    write (0, *) "UNTRAPPED ERROR: ", err%code
    write (0, *) "UNTRAPPED ERROR: ", err%message
    ! this routine should never be called, call stop() in utter desparation FIXME?
    stop
  end subroutine untrapped_error

  subroutine set_error_alloc(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_alloc
  end subroutine set_error_alloc

  subroutine set_error_fatal(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_fatal
  end subroutine set_error_fatal

  subroutine set_error_dealloc(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_dealloc
  end subroutine set_error_dealloc

  subroutine set_error_mpi(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_mpi
  end subroutine set_error_mpi

  subroutine set_error_input(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_input_command
  end subroutine set_error_input

  subroutine set_error_file(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_file_read
  end subroutine set_error_file

  subroutine set_error_open(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_file_open
  end subroutine set_error_open

  ! BGS should this just be unconv if lapack ierr /= 0?
  subroutine set_error_lapack(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_matrix_lib
  end subroutine set_error_lapack

  ! BGS what should this be? fatal?
  subroutine set_error_not_unitary(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_not_unitary
  end subroutine set_error_not_unitary

  ! BGS what should this be? fatal?
  subroutine set_error_sym(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_sitesym
  end subroutine set_error_sym

  subroutine set_error_unconv(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    ! BGS should io_stopwatch call this or something else?
    err%message = mesg !FIXME, trim to 120
    err%code = code_unconv
  end subroutine set_error_unconv

  subroutine set_error_plot(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_plot
  end subroutine set_error_plot

  subroutine set_warning(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_warning ! trivial coding error in io_stopwatch etc...
  end subroutine set_warning

  subroutine set_error_dis(err, mesg)
    type(w90_error_type), allocatable, intent(out) :: err
    character(len=*), intent(in) :: mesg
    allocate (err)
    err%message = mesg !FIXME, trim to 120
    err%code = code_disentangle
  end subroutine set_error_dis

end module w90_error

