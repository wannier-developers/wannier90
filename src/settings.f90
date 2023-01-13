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
!                                                            !
!------------------------------------------------------------!

module w90_settings

  use w90_constants, only: maxlen, dp
  implicit none

  !public

  type settings_data
    !!==================================================
    !! structure to hold a scalar and array setting
    !!==================================================
    ! for simplicity, consider arrays of different rank
    ! as different types; otherwise reshape, etc.
    character(len=maxlen) :: keyword ! token
    character(len=maxlen) :: txtdata ! text data item
    ! integer data
    integer, allocatable :: i1d(:)
    integer, allocatable :: i2d(:, :)
    integer :: idata
    ! logical data
    logical, allocatable :: l1d(:)
    logical :: ldata
    ! fp data
    real(kind=dp), allocatable :: r1d(:)
    real(kind=dp), allocatable :: r2d(:, :)
    real(kind=dp) :: rdata
  end type settings_data
  type settings_type
    integer :: num_entries = 0, num_entries_max = 0 ! number of keywords stored and max
    type(settings_data), allocatable :: entries(:)
    ! data for processing input file
    integer :: num_lines
    character(len=maxlen), allocatable :: in_data(:) ! contents of .win file
  end type settings_type
  type(settings_type) :: settings

  private :: expand_settings
  private :: init_settings
  public :: set_option
  public :: update_settings

  interface set_option
    module procedure set_option_logical
    !module procedure set_option_b1d
    module procedure set_option_text
    module procedure set_option_i1d
    module procedure set_option_i2d
    module procedure set_option_int
    module procedure set_option_r1d
    module procedure set_option_r2d
    module procedure set_option_real
  end interface set_option

contains
  !================================================!
  subroutine init_settings()
    implicit none
    integer :: defsize = 10
    allocate (settings%entries(defsize))
    settings%num_entries = 0
    settings%num_entries_max = defsize
  end subroutine init_settings

  subroutine expand_settings() ! this is a compromise to avoid a fixed size
    type(settings_data), allocatable :: nentries(:); 
    integer :: n, m ! old, new sizes
    integer :: incsize = 10
    n = settings%num_entries_max
    m = n + incsize
    allocate (nentries(m)); nentries(1:n) = settings%entries(1:n); call move_alloc(nentries, settings%entries) !f2003, note that "new" space not initialised
    settings%num_entries_max = m
  end subroutine expand_settings

  subroutine update_settings(keyword, bool, text, rval, ival, i1d, i2d, r1d, r2d)
    implicit none
    character(*), optional, intent(in) :: keyword
    character(*), optional, intent(in) :: text
    integer, optional, intent(in) :: i1d(:)
    integer, optional, intent(in) :: i2d(:, :)
    integer, optional, intent(in) :: ival
    logical, optional, intent(in) :: bool
    real(kind=dp), optional, intent(in) :: r1d(:)
    real(kind=dp), optional, intent(in) :: r2d(:, :)
    real(kind=dp), optional, intent(in) :: rval
    integer :: i
    if (.not. allocated(settings%entries)) call init_settings()
    i = settings%num_entries
    settings%entries(i)%keyword = keyword
    settings%entries(i)%txtdata = text
    settings%entries(i)%i1d = i1d ! this causes an automatic allocation
    settings%entries(i)%i2d = i2d
    settings%entries(i)%idata = ival
    settings%entries(i)%ldata = bool
    settings%entries(i)%r1d = r1d
    settings%entries(i)%r2d = r2d
    settings%entries(i)%rdata = rval
    settings%num_entries = i + 1
    if (settings%num_entries == settings%num_entries_max) call expand_settings()
  end subroutine update_settings

  subroutine set_option_text(string, text)
    implicit none
    character(*), intent(in) :: string
    character(*), intent(in) :: text
    call update_settings(keyword=string, text=text)
  endsubroutine set_option_text

  subroutine set_option_logical(string, bool)
    implicit none
    character(*), intent(in) :: string
    logical, intent(in) :: bool
    call update_settings(keyword=string, bool=bool)
  endsubroutine set_option_logical

  subroutine set_option_i1d(string, arr)
    implicit none
    character(*), intent(in) :: string
    integer, intent(in) :: arr(:)
    call update_settings(keyword=string, i1d=arr)
  endsubroutine set_option_i1d

  subroutine set_option_i2d(string, arr)
    implicit none
    character(*), intent(in) :: string
    integer, intent(in) :: arr(:, :)
    call update_settings(keyword=string, i2d=arr)
  endsubroutine set_option_i2d

  subroutine set_option_int(string, ival)
    implicit none
    character(*), intent(in) :: string
    integer, intent(in) :: ival
    call update_settings(keyword=string, ival=ival)
  endsubroutine set_option_int

  subroutine set_option_r1d(string, arr)
    implicit none
    character(*), intent(in) :: string
    real(kind=dp), intent(in) :: arr(:)
    call update_settings(keyword=string, r1d=arr)
  endsubroutine set_option_r1d

  subroutine set_option_r2d(string, arr)
    implicit none
    character(*), intent(in) :: string
    real(kind=dp), intent(in) :: arr(:, :)
    call update_settings(keyword=string, r2d=arr)
  endsubroutine set_option_r2d

  subroutine set_option_real(string, rval)
    implicit none
    character(*), intent(in) :: string
    real(kind=dp), intent(in) :: rval
    call update_settings(keyword=string, rval=rval)
  endsubroutine set_option_real
end module w90_settings
