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

module w90_conv
  !! Module to convert checkpoint files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use w90_constants, only: dp
  use w90_io, only: stdout, io_error, seedname
  implicit none

  logical, save :: export_flag
  character(len=33), save :: header
contains

  subroutine print_usage()
    !! Writes the usage of the program to stdout
    write (stdout, '(A)') "Usage:"
    write (stdout, '(A)') "  w90chk2chk.x ACTION [SEEDNAME]"
    write (stdout, '(A)') "where ACTION can be one of the following:"
    write (stdout, '(A)') "  -export"
    write (stdout, '(A)') "  -u2f"
    write (stdout, '(A)') "      Convert from unformatted (standard) format to formatted format, to export"
    write (stdout, '(A)') "      the checkpoint file on a different machine."
    write (stdout, '(A)') "      The seedname.chk file is read and the seedname.chk.fmt file is generated."
    write (stdout, '(A)') "  -import"
    write (stdout, '(A)') "  -f2u"
    write (stdout, '(A)') "      Convert from formatted format to unformatted (standard) format, to import"
    write (stdout, '(A)') "      the checkpoint file seedname.chk.fmt from a different machine."
    write (stdout, '(A)') "      The seedname.chk.fmt file is read and the seedname.chk file is generated."
  end subroutine print_usage

  subroutine conv_get_seedname
    !! Set the seedname from the command line
    implicit none

    integer :: num_arg
    character(len=50) :: ctemp

    num_arg = command_argument_count()
    if (num_arg == 1) then
      seedname = 'wannier'
    elseif (num_arg == 2) then
      call get_command_argument(2, seedname)
    else
      call print_usage
      call io_error('Wrong command line arguments, see logfile for usage')
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
      call print_usage
      call io_error('Wrong command line arguments, see logfile for usage')
    end if

  end subroutine conv_get_seedname

  !=======================================!
  subroutine conv_read_chkpt()
    !=======================================!
    !! Read formatted checkpoint file
    !=======================================!

    use w90_constants, only: eps6
    use w90_io, only: io_error, io_file_unit, stdout, seedname
    use w90_parameters

    implicit none

    integer :: chk_unit, i, j, k, l, nkp, ierr

    write (stdout, '(1x,3a)') 'Reading information from file ', trim(seedname), '.chk :'

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk', status='old', form='unformatted', err=121)

    ! Read comment line
    read (chk_unit) header
    write (stdout, '(1x,a)') trim(header)

    ! Consistency checks
    read (chk_unit) num_bands                           ! Number of bands
    write (stdout, '(a,i0)') "Number of bands: ", num_bands
    read (chk_unit) num_exclude_bands                   ! Number of excluded bands
    if (num_exclude_bands < 0) then
      call io_error('Invalid value for num_exclude_bands')
    endif
    allocate (exclude_bands(num_exclude_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating exclude_bands in conv_read_chkpt')
    read (chk_unit) (exclude_bands(i), i=1, num_exclude_bands) ! Excluded bands
    write (stdout, '(a)', advance='no') "Excluded bands: "
    if (num_exclude_bands == 0) then
      write (stdout, '(a)') "none."
    else
      do i = 1, num_exclude_bands - 1
        write (stdout, '(I0,a)', advance='no') exclude_bands(i), ','
      end do
      write (stdout, '(I0,a)') exclude_bands(num_exclude_bands), '.'
    end if
    read (chk_unit) ((real_lattice(i, j), i=1, 3), j=1, 3)  ! Real lattice
    write (stdout, '(a)') "Real lattice: read."
    read (chk_unit) ((recip_lattice(i, j), i=1, 3), j=1, 3)  ! Reciprocal lattice
    write (stdout, '(a)') "Reciprocal lattice: read."
    read (chk_unit) num_kpts                ! K-points
    write (stdout, '(a,I0)') "Num kpts:", num_kpts
    read (chk_unit) (mp_grid(i), i=1, 3)         ! M-P grid
    write (stdout, '(a)') "mp_grid: read."
    if (.not. allocated(kpt_latt)) then
      allocate (kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpt_latt in conv_read_chkpt')
    endif
    read (chk_unit) ((kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts)
    write (stdout, '(a)') "kpt_latt: read."
    read (chk_unit) nntot                ! nntot
    write (stdout, '(a,I0)') "nntot:", nntot
    read (chk_unit) num_wann                ! num_wann
    write (stdout, '(a,I0)') "num_wann:", num_wann

    read (chk_unit) checkpoint             ! checkpoint
    checkpoint = adjustl(trim(checkpoint))
    write (stdout, '(a,I0)') "checkpoint: "//trim(checkpoint)

    read (chk_unit) have_disentangled      ! whether a disentanglement has been performed

    if (have_disentangled) then
      write (stdout, '(a)') "have_disentangled: TRUE"

      read (chk_unit) omega_invariant     ! omega invariant
      write (stdout, '(a)') "omega_invariant: read."

      ! lwindow
      if (.not. allocated(lwindow)) then
        allocate (lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating lwindow in conv_read_chkpt')
      endif
      read (chk_unit, err=122) ((lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
      write (stdout, '(a)') "lwindow: read."

      ! ndimwin
      if (.not. allocated(ndimwin)) then
        allocate (ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ndimwin in conv_read_chkpt')
      endif
      read (chk_unit, err=123) (ndimwin(nkp), nkp=1, num_kpts)
      write (stdout, '(a)') "ndimwin: read."

      ! U_matrix_opt
      if (.not. allocated(u_matrix_opt)) then
        allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating u_matrix_opt in conv_read_chkpt')
      endif
      read (chk_unit, err=124) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
      write (stdout, '(a)') "U_matrix_opt: read."

    else
      write (stdout, '(a)') "have_disentangled: FALSE"
    endif

    ! U_matrix
    if (.not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating u_matrix in conv_read_chkpt')
    endif
    read (chk_unit, err=125) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)
    write (stdout, '(a)') "U_matrix: read."

    ! M_matrix
    if (.not. allocated(m_matrix)) then
      allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating m_matrix in conv_read_chkpt')
    endif
    read (chk_unit, err=126) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, nntot), l=1, num_kpts)
    write (stdout, '(a)') "M_matrix: read."

    ! wannier_centres
    if (.not. allocated(wannier_centres)) then
      allocate (wannier_centres(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in conv_read_chkpt')
    end if
    read (chk_unit, err=127) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)
    write (stdout, '(a)') "wannier_centres: read."

    ! wannier spreads
    if (.not. allocated(wannier_spreads)) then
      allocate (wannier_spreads(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in conv_read_chkpt')
    end if
    read (chk_unit, err=128) (wannier_spreads(i), i=1, num_wann)
    write (stdout, '(a)') "wannier_spreads: read."

    close (chk_unit)

    write (stdout, '(a/)') ' ... done'

    return

121 call io_error('Error opening '//trim(seedname)//'.chk in conv_read_chkpt')
122 call io_error('Error reading lwindow from '//trim(seedname)//'.chk in conv_read_chkpt')
123 call io_error('Error reading ndimwin from '//trim(seedname)//'.chk in conv_read_chkpt')
124 call io_error('Error reading u_matrix_opt from '//trim(seedname)//'.chk in conv_read_chkpt')
125 call io_error('Error reading u_matrix from '//trim(seedname)//'.chk in conv_read_chkpt')
126 call io_error('Error reading m_matrix from '//trim(seedname)//'.chk in conv_read_chkpt')
127 call io_error('Error reading wannier_centres from '//trim(seedname)//'.chk in conv_read_chkpt')
128 call io_error('Error reading wannier_spreads from '//trim(seedname)//'.chk in conv_read_chkpt')

  end subroutine conv_read_chkpt

  subroutine conv_read_chkpt_fmt()
    !=======================================!
    !! Read formatted checkpoint file
    !=======================================!

    use w90_constants, only: eps6
    use w90_io, only: io_error, io_file_unit, stdout, seedname
    use w90_parameters

    implicit none

    integer :: chk_unit, i, j, k, l, nkp, ierr, idum
    character(len=30) :: cdum
    real(kind=dp) :: rreal, rimag

    write (stdout, '(1x,3a)') 'Reading information from formatted file ', trim(seedname), '.chk.fmt :'

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk.fmt', status='old', position='rewind', form='formatted', err=121)

    ! Read comment line
    read (chk_unit, '(A)') header
    write (stdout, '(1x,a)') trim(header)

    ! Consistency checks
    read (chk_unit, *) num_bands                           ! Number of bands
    write (stdout, '(a,i0)') "Number of bands: ", num_bands
    read (chk_unit, *) num_exclude_bands                   ! Number of excluded bands
    if (num_exclude_bands < 0) then
      call io_error('Invalid value for num_exclude_bands')
    endif
    allocate (exclude_bands(num_exclude_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating exclude_bands in conv_read_chkpt_fmt')
    do i = 1, num_exclude_bands
      read (chk_unit, *) exclude_bands(i) ! Excluded bands
    end do
    write (stdout, '(a)', advance='no') "Excluded bands: "
    if (num_exclude_bands == 0) then
      write (stdout, '(a)') "none."
    else
      do i = 1, num_exclude_bands - 1
        write (stdout, '(I0,a)', advance='no') exclude_bands(i), ','
      end do
      write (stdout, '(I0,a)') exclude_bands(num_exclude_bands), '.'
    end if
    read (chk_unit, *) ((real_lattice(i, j), i=1, 3), j=1, 3)  ! Real lattice
    write (stdout, '(a)') "Real lattice: read."
    read (chk_unit, *) ((recip_lattice(i, j), i=1, 3), j=1, 3)  ! Reciprocal lattice
    write (stdout, '(a)') "Reciprocal lattice: read."
    read (chk_unit, *) num_kpts                ! K-points
    write (stdout, '(a,I0)') "Num kpts:", num_kpts
    read (chk_unit, *) (mp_grid(i), i=1, 3)         ! M-P grid
    write (stdout, '(a)') "mp_grid: read."
    if (.not. allocated(kpt_latt)) then
      allocate (kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpt_latt in conv_read_chkpt_fmt')
    endif
    do nkp = 1, num_kpts
      read (chk_unit, *, err=115) (kpt_latt(i, nkp), i=1, 3)
    end do
    write (stdout, '(a)') "kpt_latt: read."
    read (chk_unit, *) nntot                ! nntot
    write (stdout, '(a,I0)') "nntot:", nntot
    read (chk_unit, *) num_wann                ! num_wann
    write (stdout, '(a,I0)') "num_wann:", num_wann

    read (chk_unit, *) checkpoint             ! checkpoint
    checkpoint = adjustl(trim(checkpoint))
    write (stdout, '(a,I0)') "checkpoint: "//trim(checkpoint)

    read (chk_unit, *) idum
    if (idum == 1) then
      have_disentangled = .true.
    elseif (idum == 0) then
      have_disentangled = .false.
    else
      write (cdum, '(I0)') idum
      call io_error('Error reading formatted chk: have_distenangled should be 0 or 1, it is instead '//cdum)
    end if

    if (have_disentangled) then
      write (stdout, '(a)') "have_disentangled: TRUE"

      read (chk_unit, *) omega_invariant     ! omega invariant
      write (stdout, '(a)') "omega_invariant: read."

      ! lwindow
      if (.not. allocated(lwindow)) then
        allocate (lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating lwindow in conv_read_chkpt_fmt')
      endif
      do nkp = 1, num_kpts
        do i = 1, num_bands
          read (chk_unit, *) idum
          if (idum == 1) then
            lwindow(i, nkp) = .true.
          elseif (idum == 0) then
            lwindow(i, nkp) = .false.
          else
            write (cdum, '(I0)') idum
            call io_error('Error reading formatted chk: lwindow(i,nkp) should be 0 or 1, it is instead '//cdum)
          end if
        end do
      end do
      write (stdout, '(a)') "lwindow: read."

      ! ndimwin
      if (.not. allocated(ndimwin)) then
        allocate (ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating ndimwin in conv_read_chkpt_fmt')
      endif
      do nkp = 1, num_kpts
        read (chk_unit, *, err=123) ndimwin(nkp)
      end do
      write (stdout, '(a)') "ndimwin: read."

      ! U_matrix_opt
      if (.not. allocated(u_matrix_opt)) then
        allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating u_matrix_opt in conv_read_chkpt_fmt')
      endif
      do nkp = 1, num_kpts
        do j = 1, num_wann
          do i = 1, num_bands
            read (chk_unit, *, err=124) rreal, rimag
            u_matrix_opt(i, j, nkp) = cmplx(rreal, rimag, kind=dp)
          end do
        end do
      end do
      write (stdout, '(a)') "U_matrix_opt: read."

    else
      write (stdout, '(a)') "have_disentangled: FALSE"
    endif

    ! U_matrix
    if (.not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating u_matrix in conv_read_chkpt_fmt')
    endif
    do k = 1, num_kpts
      do j = 1, num_wann
        do i = 1, num_wann
          read (chk_unit, *, err=124) rreal, rimag
          u_matrix(i, j, k) = cmplx(rreal, rimag, kind=dp)
        end do
      end do
    end do
    write (stdout, '(a)') "U_matrix: read."

    ! M_matrix
    if (.not. allocated(m_matrix)) then
      allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating m_matrix in conv_read_chkpt_fmt')
    endif
    do l = 1, num_kpts
      do k = 1, nntot
        do j = 1, num_wann
          do i = 1, num_wann
            read (chk_unit, *, err=124) rreal, rimag
            m_matrix(i, j, k, l) = cmplx(rreal, rimag, kind=dp)
          end do
        end do
      end do
    end do
    write (stdout, '(a)') "M_matrix: read."

    ! wannier_centres
    if (.not. allocated(wannier_centres)) then
      allocate (wannier_centres(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in conv_read_chkpt_fmt')
    end if
    do j = 1, num_wann
      read (chk_unit, *, err=127) (wannier_centres(i, j), i=1, 3)
    end do
    write (stdout, '(a)') "wannier_centres: read."

    ! wannier spreads
    if (.not. allocated(wannier_spreads)) then
      allocate (wannier_spreads(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in conv_read_chkpt_fmt')
    end if
    do i = 1, num_wann
      read (chk_unit, *, err=128) wannier_spreads(i)
    end do
    write (stdout, '(a)') "wannier_spreads: read."

    close (chk_unit)

    write (stdout, '(a/)') ' ... done'

    return

115 call io_error('Error reading variable from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
121 call io_error('Error opening '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
122 call io_error('Error reading lwindow from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
123 call io_error('Error reading ndimwin from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
124 call io_error('Error reading u_matrix_opt from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
125 call io_error('Error reading u_matrix from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
126 call io_error('Error reading m_matrix from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
127 call io_error('Error reading wannier_centres from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')
128 call io_error('Error reading wannier_spreads from '//trim(seedname)//'.chk.fmt in conv_read_chkpt_fmt')

  end subroutine conv_read_chkpt_fmt

  subroutine conv_write_chkpt()
    !=======================================!
    !! Write formatted checkpoint file
    !=======================================!

    use w90_io, only: io_file_unit, io_date, seedname
    use w90_parameters

    implicit none

    integer :: chk_unit, nkp, i, j, k, l
    character(len=9) :: cdate, ctime

    write (stdout, '(/1x,3a)', advance='no') 'Writing checkpoint file ', trim(seedname), '.chk...'

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk', form='unformatted')

    write (chk_unit) header                                   ! Date and time from the read file
    write (chk_unit) num_bands                                ! Number of bands
    write (chk_unit) num_exclude_bands                        ! Number of excluded bands
    write (chk_unit) (exclude_bands(i), i=1, num_exclude_bands) ! Excluded bands
    write (chk_unit) ((real_lattice(i, j), i=1, 3), j=1, 3)        ! Real lattice
    write (chk_unit) ((recip_lattice(i, j), i=1, 3), j=1, 3)       ! Reciprocal lattice
    write (chk_unit) num_kpts                                 ! Number of k-points
    write (chk_unit) (mp_grid(i), i=1, 3)                       ! M-P grid
    write (chk_unit) ((kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts) ! K-points
    write (chk_unit) nntot                                    ! Number of nearest k-point neighbours
    write (chk_unit) num_wann                                 ! Number of wannier functions
    ! Next is correct: it always print out 20 characters
    write (chk_unit) checkpoint                               ! Position of checkpoint
    write (chk_unit) have_disentangled      ! Whether a disentanglement has been performed
    if (have_disentangled) then
      write (chk_unit) omega_invariant     ! Omega invariant
      ! lwindow, ndimwin and U_matrix_opt
      write (chk_unit) ((lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
      write (chk_unit) (ndimwin(nkp), nkp=1, num_kpts)
      write (chk_unit) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
    endif
    write (chk_unit) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)               ! U_matrix
    write (chk_unit) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, nntot), l=1, num_kpts) ! M_matrix
    write (chk_unit) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)
    write (chk_unit) (wannier_spreads(i), i=1, num_wann)
    close (chk_unit)

    write (stdout, '(a/)') ' done'

  end subroutine conv_write_chkpt

  subroutine conv_write_chkpt_fmt()
    !=======================================!
    !! Write formatted checkpoint file
    !=======================================!

    use w90_io, only: io_file_unit, io_date, seedname
    use w90_parameters

    implicit none

    integer :: chk_unit, nkp, i, j, k, l
    character(len=9) :: cdate, ctime

    write (stdout, '(/1x,3a)', advance='no') 'Writing formatted checkpoint file ', trim(seedname), '.chk.fmt...'

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk.fmt', form='formatted', status='replace', position='rewind')

    write (chk_unit, '(A33)') header                                   ! Date and time from the read file
    write (chk_unit, '(I0)') num_bands                                ! Number of bands
    write (chk_unit, '(I0)') num_exclude_bands                        ! Number of excluded bands
    do i = 1, num_exclude_bands
      write (chk_unit, '(I0)') exclude_bands(i) ! Excluded bands
    end do
    write (chk_unit, '(9G25.17)') ((real_lattice(i, j), i=1, 3), j=1, 3)        ! Real lattice
    write (chk_unit, '(9G25.17)') ((recip_lattice(i, j), i=1, 3), j=1, 3)       ! Reciprocal lattice
    write (chk_unit, '(I0)') num_kpts                                 ! Number of k-points
    write (chk_unit, '(I0," ",I0," ",I0)') (mp_grid(i), i=1, 3)                       ! M-P grid
    do nkp = 1, num_kpts
      write (chk_unit, '(3G25.17)') (kpt_latt(i, nkp), i=1, 3) ! K-points
    end do
    write (chk_unit, '(I0)') nntot                                    ! Number of nearest k-point neighbours
    write (chk_unit, '(I0)') num_wann                                 ! Number of wannier functions
    write (chk_unit, '(A20)') checkpoint                               ! Position of checkpoint
    if (have_disentangled) then
      write (chk_unit, '(I1)') 1      ! Whether a disentanglement has been performed
    else
      write (chk_unit, '(I1)') 0      ! Whether a disentanglement has been performed
    end if
    if (have_disentangled) then
      write (chk_unit, '(G25.17)') omega_invariant     ! Omega invariant
      ! lwindow, ndimwin and U_matrix_opt
      do nkp = 1, num_kpts
        do i = 1, num_bands
          if (lwindow(i, nkp)) then
            write (chk_unit, '(I1)') 1
          else
            write (chk_unit, '(I1)') 0
          end if
        end do
      end do
      do nkp = 1, num_kpts
        write (chk_unit, '(I0)') ndimwin(nkp)
      end do
      do nkp = 1, num_kpts
        do j = 1, num_wann
          do i = 1, num_bands
            write (chk_unit, '(2G25.17)') u_matrix_opt(i, j, nkp)
          end do
        end do
      end do
    endif
    do k = 1, num_kpts
      do j = 1, num_wann
        do i = 1, num_wann
          write (chk_unit, '(2G25.17)') u_matrix(i, j, k)
        end do
      end do
    end do
    do l = 1, num_kpts
      do k = 1, nntot
        do j = 1, num_wann
          do i = 1, num_wann
            write (chk_unit, '(2G25.17)') m_matrix(i, j, k, l)
          end do
        end do
      end do
    end do
    do j = 1, num_wann
      write (chk_unit, '(3G25.17)') (wannier_centres(i, j), i=1, 3)
    end do
    do i = 1, num_wann
      write (chk_unit, '(G25.17)') wannier_spreads(i)
    end do
    close (chk_unit)

    write (stdout, '(a/)') ' done'

  end subroutine conv_write_chkpt_fmt

end module w90_conv

program w90chk2chk
  !! Program to convert checkpoint files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use w90_constants, only: dp
  use w90_io, only: io_file_unit, stdout, io_error, seedname
  use w90_conv
  use w90_comms, only: num_nodes, comms_setup, comms_end
  implicit none

  ! Export mode:
  !  TRUE:  create formatted .chk.fmt from unformatted .chk ('-export')
  !  FALSE: create unformatted .chk from formatted .chk.fmt ('-import')
  logical :: file_found
  integer :: file_unit

  call comms_setup

  stdout = io_file_unit()
  open (unit=stdout, file='w90chk2chk.log')

  if (num_nodes /= 1) then
    call io_error('w90chk2chk can only be used in serial...')
  endif

  call conv_get_seedname

  if (export_flag .eqv. .true.) then
    call conv_read_chkpt()
    call conv_write_chkpt_fmt()
  else
    call conv_read_chkpt_fmt()
    call conv_write_chkpt()
  end if

!  close(unit=stdout,status='delete')
  close (unit=stdout)

  call comms_end

end program w90chk2chk

