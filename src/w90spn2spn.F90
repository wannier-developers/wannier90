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
!  w90spn2spn: reformat spn files                            !
!                                                            !
!------------------------------------------------------------!

module w90spn_parameters

  use w90_types

  implicit none

  public

  integer, save :: num_kpts !BGS put in k_point_type?
  integer, save :: num_bands !! Number of bands

end module w90spn_parameters

module w90_conv_spn
  !! Module to convert spn files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use w90_constants, only: dp

  implicit none

  logical, save :: export_flag
  ! header length is identical to that of pw2wannier90.f90
  character(len=60), save :: header
  complex(kind=dp), allocatable, save :: spn_o(:, :, :, :)

contains

  subroutine io_error(error_msg, stdout)
    !================================================
    !
    !! Abort the code giving an error message
    !
    !================================================

    implicit none

    character(len=*), intent(in) :: error_msg
    integer :: stdout

    close (stdout)
    write (*, '(1x,a)') trim(error_msg)
    write (*, '(A)') "Error: examine the output/error file for details"
    stop
  end subroutine io_error

  !================================================!
  subroutine print_usage(stdout)
    !================================================!
    !
    !! Writes the usage of the program to stdout
    !
    !================================================!

    implicit none

    integer, intent(in) :: stdout

    write (stdout, '(A)') "Usage:"
    write (stdout, '(A)') "  w90spn2spn.x ACTION [SEEDNAME]"
    write (stdout, '(A)') "where ACTION can be one of the following:"
    write (stdout, '(A)') "  -export"
    write (stdout, '(A)') "  -u2f"
    write (stdout, '(A)') "      Convert from unformatted (standard) format to formatted format, to export"
    write (stdout, '(A)') "      the spn file on a different machine."
    write (stdout, '(A)') "      The seedname.spn file is read and the seedname.spn.fmt file is generated."
    write (stdout, '(A)') "  -import"
    write (stdout, '(A)') "  -f2u"
    write (stdout, '(A)') "      Convert from formatted format to unformatted (standard) format, to import"
    write (stdout, '(A)') "      the spn file seedname.spn.fmt from a different machine."
    write (stdout, '(A)') "      The seedname.spn.fmt file is read and the seedname.spn file is generated."
  end subroutine print_usage

  !================================================!
  subroutine conv_get_seedname(stdout, seedname)
    !================================================!
    !
    !! Set the seedname from the command line
    !
    !================================================!
    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(inout)  :: seedname

    integer :: num_arg
    character(len=50) :: ctemp

    num_arg = command_argument_count()
    if (num_arg == 1) then
      seedname = 'wannier'
    elseif (num_arg == 2) then
      call get_command_argument(2, seedname)
    else
      call print_usage(stdout)
      call io_error('Wrong command line arguments, see logfile for usage', stdout)
    end if

    ! If on the command line the whole seedname.win was passed, I strip the last ".win"
    if (len(trim(seedname)) .ge. 5) then
      if (seedname(len(trim(seedname)) - 4 + 1:) .eq. ".win") then
        seedname = seedname(:len(trim(seedname)) - 4)
      end if
    end if

    call get_command_argument(1, ctemp)
    if (index(ctemp, '-import') > 0) then
      export_flag = .false.
    elseif (index(ctemp, '-f2u') > 0) then
      export_flag = .false.
    elseif (index(ctemp, '-export') > 0) then
      export_flag = .true.
    elseif (index(ctemp, '-u2f') > 0) then
      export_flag = .true.
    else
      write (stdout, '(A)') 'Wrong command line action: '//trim(ctemp)
      call print_usage(stdout)
      call io_error('Wrong command line arguments, see logfile for usage', stdout)
    end if

  end subroutine conv_get_seedname

  !================================================!
  subroutine conv_read_spn(stdout, seedname)
    !================================================!
    !
    !! Read unformatted spn file
    !
    !================================================!

    use w90_constants, only: eps6, dp
    use w90spn_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: spn_unit, m, n, ik, ierr, s, counter
    complex(kind=dp), allocatable :: spn_temp(:, :)

    write (stdout, '(3a)') 'Reading information from unformatted file ', trim(seedname), '.spn :'

    open (newunit=spn_unit, file=trim(seedname)//'.spn', status='old', form='unformatted', err=109)

    ! Read comment line
    read (spn_unit, err=110, end=110) header
    header = ADJUSTL(header)
    write (stdout, '(1x,a)') trim(header)

    ! Consistency checks
    read (spn_unit, err=110, end=110) num_bands, num_kpts
    write (stdout, '(1x,a,i0)') "Number of bands: ", num_bands
    write (stdout, '(1x,a,i0)') "Number of k-points: ", num_kpts

    allocate (spn_o(num_bands, num_bands, num_kpts, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating spm_temp in conv_read_spn', stdout)

    allocate (spn_temp(3, (num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating spm_temp in conv_read_spn', stdout)
    do ik = 1, num_kpts
      read (spn_unit) ((spn_temp(s, m), s=1, 3), m=1, (num_bands*(num_bands + 1))/2)
      counter = 0
      do m = 1, num_bands
        do n = 1, m
          counter = counter + 1
          spn_o(n, m, ik, 1) = spn_temp(1, counter)
          spn_o(n, m, ik, 2) = spn_temp(2, counter)
          spn_o(n, m, ik, 3) = spn_temp(3, counter)
          ! Although each diagonal element of spn_o should be a real number,
          ! actually it has a very small imaginary part.
          ! We skip the conjugation on the diagonal elements so that
          ! the file after formatted <==> unformatted conversions is exactly
          ! the same as the original file, otherwise the diagonal elements
          ! are the conjugations of those of the original file.
          if (m == n) cycle
          spn_o(m, n, ik, 1) = conjg(spn_temp(1, counter))
          spn_o(m, n, ik, 2) = conjg(spn_temp(2, counter))
          spn_o(m, n, ik, 3) = conjg(spn_temp(3, counter))
        end do
      end do
    end do

    close (spn_unit)
    write (stdout, '(1x,a)') "spn: read."

    deallocate (spn_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating spm_temp in conv_read_spn', stdout)

    write (stdout, '(1x,a)') 'read done.'

    return

109 call io_error('Error opening '//trim(seedname)//'.spn.fmt in conv_read_spn', stdout)
110 call io_error('Error reading '//trim(seedname)//'.spn.fmt in conv_read_spn', stdout)

  end subroutine conv_read_spn

  !================================================!
  subroutine conv_read_spn_fmt(stdout, seedname)
    !================================================!
    !
    !! Read formatted spn file
    !
    !================================================!

    use w90_constants, only: eps6, dp
    use w90spn_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: spn_unit, m, n, ik, ierr
    real(kind=dp) :: s_real, s_img

    write (stdout, '(3a)') 'Reading information from formatted file ', trim(seedname), '.spn.fmt :'

    open (newunit=spn_unit, file=trim(seedname)//'.spn.fmt', status='old', position='rewind', form='formatted', err=109)

    ! Read comment line
    read (spn_unit, '(a)') header
    header = ADJUSTL(header)
    write (stdout, '(1x,a)') trim(header)

    ! Consistency checks
    read (spn_unit, *, err=110, end=110) num_bands, num_kpts
    write (stdout, '(1x,a,i0)') "Number of bands: ", num_bands
    write (stdout, '(1x,a,i0)') "Number of k-points: ", num_kpts

    allocate (spn_o(num_bands, num_bands, num_kpts, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating spn_o in conv_read_spn_fmt', stdout)

    do ik = 1, num_kpts
      do m = 1, num_bands
        do n = 1, m
          read (spn_unit, *, err=110, end=110) s_real, s_img
          spn_o(n, m, ik, 1) = cmplx(s_real, s_img, dp)
          read (spn_unit, *, err=110, end=110) s_real, s_img
          spn_o(n, m, ik, 2) = cmplx(s_real, s_img, dp)
          read (spn_unit, *, err=110, end=110) s_real, s_img
          spn_o(n, m, ik, 3) = cmplx(s_real, s_img, dp)
          ! Although each diagonal element of spn_o should be a real number,
          ! actually it has a very small imaginary part.
          ! We skip the conjugation on the diagonal elements so that
          ! the file after formatted <==> unformatted conversions is exactly
          ! the same as the original file, otherwise the diagonal elements
          ! are the conjugations of those of the original file.
          if (m == n) cycle
          ! Read upper-triangular part, now build the rest
          spn_o(m, n, ik, 1) = conjg(spn_o(n, m, ik, 1))
          spn_o(m, n, ik, 2) = conjg(spn_o(n, m, ik, 2))
          spn_o(m, n, ik, 3) = conjg(spn_o(n, m, ik, 3))
        end do
      end do
    enddo

    close (spn_unit)

    write (stdout, '(1x,a)') 'read done.'

    return

109 call io_error('Error opening '//trim(seedname)//'.spn.fmt in conv_read_spn_fmt', stdout)
110 call io_error('Error reading '//trim(seedname)//'.spn.fmt in conv_read_spn_fmt', stdout)

  end subroutine conv_read_spn_fmt

  !================================================!
  subroutine conv_write_spn(stdout, seedname)
    !================================================!
    !
    !! Write unformatted spn file
    !
    !================================================!

    use w90_io, only: io_date
    use w90spn_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: spn_unit, m, n, ik, counter, s, ierr
    complex(kind=dp), allocatable :: spn_temp(:, :)

    write (stdout, '(3a)') 'Writing information to unformatted file ', trim(seedname), '.spn :'

    open (newunit=spn_unit, file=trim(seedname)//'.spn', form='unformatted')

    allocate (spn_temp(3, (num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating spm_temp in conv_write_spn', stdout)

    write (spn_unit) header
    write (spn_unit) num_bands, num_kpts

    do ik = 1, num_kpts
      counter = 0
      do m = 1, num_bands
        do n = 1, m
          counter = counter + 1
          do s = 1, 3
            spn_temp(s, counter) = spn_o(n, m, ik, s)
          end do
        end do
      end do
      write (spn_unit) ((spn_temp(s, m), s=1, 3), m=1, ((num_bands*(num_bands + 1))/2))
    end do

    close (spn_unit)

    write (stdout, '(1x,a)') 'write done.'

  end subroutine conv_write_spn

  !================================================!
  subroutine conv_write_spn_fmt(stdout, seedname)
    !================================================!
    !
    !! Write formatted spn file
    !
    !================================================!

    use w90_io, only: io_date
    use w90spn_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: spn_unit, m, n, ik, s

    write (stdout, '(3a)') 'Writing information to formatted file ', trim(seedname), '.spn.fmt :'

    open (newunit=spn_unit, file=trim(seedname)//'.spn.fmt', form='formatted', status='replace', position='rewind')

    write (spn_unit, *) header
    write (spn_unit, *) num_bands, num_kpts

    do ik = 1, num_kpts
      do m = 1, num_bands
        do n = 1, m
          do s = 1, 3
            write (spn_unit, '(2es26.16)') spn_o(n, m, ik, s)
          enddo
        enddo
      enddo
    end do

    close (spn_unit)

    write (stdout, '(1x, a)') 'write done.'

  end subroutine conv_write_spn_fmt

end module w90_conv_spn

program w90spn2spn
  !! Program to convert spn files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use w90_constants, only: dp
  use w90_conv_spn

  implicit none

  ! Export mode:
  !  TRUE:  create formatted .spn.fmt from unformatted .spn ('-export')
  !  FALSE: create unformatted .spn from formatted .spn.fmt ('-import')

  !logical :: file_found
  !integer :: file_unit
  integer :: stdout !, ierr, num_nodes
  character(len=50) :: seedname

  open (newunit=stdout, file='w90spn2spn.log')

  call conv_get_seedname(stdout, seedname)

  if (export_flag .eqv. .true.) then
    call conv_read_spn(stdout, seedname)
    write (stdout, '(a)') ''
    call conv_write_spn_fmt(stdout, seedname)
  else
    call conv_read_spn_fmt(stdout, seedname)
    write (stdout, '(a)') ''
    call conv_write_spn(stdout, seedname)
  end if

  close (unit=stdout)

end program w90spn2spn
