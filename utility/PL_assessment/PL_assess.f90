!---------------------------------------------
!               PL_assess.f90                !
!---------------------------------------------
!        Written by M.Shelley 09/09          !
!         Imperial College, London           !
!                                            !
!    Copyright (c) 2009 Matthew Shelley      !
!                                            !
!---------------------------------------------
!                                            !
! Function: To assess principal layer size   !
!           of a system from a single cell   !
!           many k-point calculation.        !
!                                            !
!--------------------------------------------!

program PL_assessment

  implicit none

  integer, parameter :: dp = kind(1.d0)

  integer :: nk(3), num_wann, tran_dir, ik1, ik2, ik3, ik_counter, i, j, k1, k2, k3, ni, nj, k, ioerror, ik, i_counter
  real(dp) :: imag_h
  real(dp), allocatable :: hr(:, :, :), ave(:), st_dev(:), max_abs_val(:)
  character(40) :: seedname
  character(33) :: null_read

  call get_data(nk, num_wann)

  ioerror = 0

  write (*, *) 'Enter seedname: '
  read (*, *) seedname
  open (unit=10, file=trim(seedname)//'_hr.dat', action='read', iostat=ioerror)
  if (ioerror /= 0) then
    write (*, *) 'Error opening ', trim(seedname)//'_hr.dat : ', ioerror
    stop
  endif

  i_counter = 0
  do i = 1, 3
    if (nk(i) .lt. 1) then
      write (*, *) ' Error: Incorrect k-point input'
      stop
    endif
    !
    if (nk(i) .eq. 1) then
      nk(i) = 0
      i_counter = i_counter + 1
    else
      tran_dir = i
    endif
    !
    if ((i_counter .ne. 2) .and. (i .eq. 3)) then
      write (*, *) ' Error: Only transport direction may have more than one k-point'
      stop
    endif
    if (i_counter .eq. 3) then
      write (*, *) ' Error: Transport direction must have more than one k-point'
      stop
    endif
  enddo

  allocate (hr(num_wann, num_wann, nk(tran_dir)/2 + 1))
  allocate (ave(nk(tran_dir)/2 + 1))
  allocate (st_dev(nk(tran_dir)/2 + 1))
  allocate (max_abs_val(nk(tran_dir)/2 + 1))

  read (10, '(33a)', iostat=ioerror) null_read
  if (ioerror /= 0) then
    write (*, *) 'Error reading ', trim(seedname)//'_hr.dat : ', ioerror
    stop
  endif

  do ik1 = -nk(1)/2, 0
    do ik2 = -nk(2)/2, 0
      do ik3 = -nk(3)/2, 0
        if (tran_dir .eq. 1) then
          ik_counter = -ik1 + 1
        elseif (tran_dir .eq. 2) then
          ik_counter = -ik2 + 1
        elseif (tran_dir .eq. 3) then
          ik_counter = -ik3 + 1
        endif
        do i = 1, num_wann
          do j = 1, num_wann
            read (10, *, iostat=ioerror) k1, k2, k3, ni, nj, hr(i, j, ik_counter), imag_h
            if (ioerror /= 0) then
              write (*, *) 'Error reading file: ', trim(seedname)//'_hr.dat : ', ioerror
              stop
            endif
          enddo
        enddo
      enddo
    enddo
  enddo

  ave = 0.d0
  do ik = 1, nk(tran_dir)/2 + 1
    do i = 1, num_wann
      ave(ik) = ave(ik) + abs(hr(i, i, ik))
    enddo
  enddo
  ave = ave/num_wann

  st_dev = 0
  do ik = 1, nk(tran_dir)/2 + 1
    do i = 1, num_wann
      st_dev(ik) = st_dev(ik) + (abs(hr(i, i, ik)) - ave(ik))**2
    enddo
  enddo
  st_dev = dsqrt(st_dev/num_wann)

  max_abs_val = 0.d0
  do ik = 1, nk(tran_dir)/2 + 1
    max_abs_val(ik) = maxval(abs(hr(:, :, ik)))
  enddo

  open (unit=11, file=trim(seedname)//'_pl.dat', action='write', iostat=ioerror)
  if (ioerror /= 0) then
    write (*, *) 'Error opening ', trim(seedname)//'_pl.dat : ', ioerror
    stop
  endif

  do ik = 1, nk(tran_dir)/2 + 1
    write (11, '(i4,3(2x,e12.6))') ik - 1, ave(ik), st_dev(ik), max_abs_val(ik)
  enddo

  close (10)
  close (11)

  write (*, *) trim(seedname)//'_pl.dat written'

  deallocate (ave)
  deallocate (st_dev)
  deallocate (hr)

contains

  subroutine get_data(nk, num_wann)

    implicit none

    integer, intent(out) :: nk(3), num_wann

    integer :: num_arg
    character(4) :: n1_char, n2_char, n3_char, num_wann_char

    num_arg = command_argument_count()

    if (num_arg == 4) then
      call get_command_argument(1, n1_char)
      call get_command_argument(2, n2_char)
      call get_command_argument(3, n3_char)
      call get_command_argument(4, num_wann_char)
    else
      write (*, *) 'Input incorrect'
      stop
    endif

    nk = 0
    num_wann = 0

    read (n1_char, '(i4)') nk(1)
    read (n2_char, '(i4)') nk(2)
    read (n3_char, '(i4)') nk(3)
    read (num_wann_char, '(i4)') num_wann

    if (num_wann .lt. 1) then
      write (*, *) ' Error: Incorrect wannier function input'
      stop
    endif

  endsubroutine get_data

end program PL_assessment
