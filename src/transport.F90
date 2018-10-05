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
!
!-----------------------------------------------------------------------!
!  Based on                                                             !
!  < dosqc_1.0 >                                                        !
!  Density Of States and Quantum Conductance - Version 1.0              !
!  Marco Buongiorno Nardelli, January 2000.                             !
!                                                                       !
!  Reference:                                                           !
!  - M. Buongiorno Nardelli, "Electronic transport in extended systems: !
!  application to carbon nanotubes", Phys. Rev. B, vol. 60(11), 7828    !
!  (1999)                                                               !
!-----------------------------------------------------------------------!

!=======================================================================!
! Definition of parameters used in w90_transport                        !
!=======================================================================!
!                                                                       !
! transport_mode       = 'bulk' or 'lcr'                                !
! tran_win_min         = minimum E                                      !
! tran_win_max         = maximum E                                      !
! tran_energy_step     = delta E                                        !
! tran_num_bb          = # of WFs in a principal layer of a perfectly   !
!                        periodic bulk system                           !
! tran_num_ll          = # of WFs in a principal layer of a left lead   !
! tran_num_rr          = # of WFs in a principal layer of a right lead  !
! tran_num_cc          = # of WFs in a disordered conducter cell        !
! tran_num_lc          = # of WFs in a disordered conducter cell that   !
!                        are used to calculate interaction with a       !
!                        left lead                                      !
! tran_num_cr          = # of WFs in a disordered conducter cell that   !
!                        are used to calculate interaction with a       !
!                        right lead                                     !
! tran_num_bandc       = width of band-diagonal hC matrix               !
! tran_read_ht         = .true. => read H matrix from h*.dat files      !
! tran_write_ht        = .true. => write H matrix from h*.dat files     !
! tran_use_same_lead   = .true. => in L-C-R construction, left and      !
!                                  right lead are the same kind         !
! tran_num_cell_ll     = # of unit cells in a PL of left lead           !
! tran_num_cell_rr     = # of unit cells in a PL of right lead          !
!                        (equal to tran_num_cell_ll for now)            !
! tran_group_threshold = distance defining the grouping of WFs          !
!=======================================================================!

module w90_transport
  !! Module to handle ballistic transport.
  !!Based on
  !!  < dosqc_1.0 >
  !!  Density Of States and Quantum Conductance - Version 1.0
  !!  Marco Buongiorno Nardelli, January 2000.
  !!  Reference:
  !!  - M. Buongiorno Nardelli, "Electronic transport in extended systems:
  !!  application to carbon nanotubes", Phys. Rev. B, vol. 60(11), 7828
  !!  (1999)

  use w90_constants, only: dp

  implicit none

  private

  complex(kind=dp), parameter :: eta = (0.0_dp, 0.0005_dp)
  !! small complex number

  integer, parameter :: nterx = 50
  !! nterx  = # of maximum iteration to calculate transfer matrix
  integer :: one_dim_vec
  !! cartesian axis to which real_lattice(:,one_dim_vec) is parallel
  integer :: nrpts_one_dim
  integer :: num_pl
  !! number of unit cell in a principal layer
  integer, dimension(3) :: coord
  !! coord : coord(1) defines the conduction direction according to 1=x,2=y,3=z,
  !! coord(2),coord(3) define the other directions during sorting routines
  integer, allocatable :: tran_sorted_idx(:)
  !! index of sorted WF centres to unsorted

  real(kind=dp), allocatable :: hr_one_dim(:, :, :)
  real(kind=dp), allocatable :: hB0(:, :)
  real(kind=dp), allocatable :: hB1(:, :)
  real(kind=dp), allocatable :: hL0(:, :)
  real(kind=dp), allocatable :: hL1(:, :)
  real(kind=dp), allocatable :: hR0(:, :)
  real(kind=dp), allocatable :: hR1(:, :)
  real(kind=dp), allocatable :: hC(:, :)
  real(kind=dp), allocatable :: hLC(:, :)
  real(kind=dp), allocatable :: hCR(:, :)

  public :: tran_main
  public :: tran_dealloc

contains
  !==================================================================!
  subroutine tran_main()
    !! Main transport subroutine
    !==================================================================!

    use w90_io, only: stdout, io_stopwatch
    use w90_parameters, only: transport_mode, tran_read_ht, timing_level, write_hr, &
      write_xyz
    use w90_hamiltonian, only: hamiltonian_get_hr, hamiltonian_write_hr, hamiltonian_setup

    implicit none

    real(kind=dp), allocatable, dimension(:, :)     :: signatures
    integer                                      :: num_G
    logical                                      :: pl_warning

    if (timing_level > 0) call io_stopwatch('tran: main', 1)

    write (stdout, '(/1x,a)') '*---------------------------------------------------------------------------*'
    write (stdout, '(1x,a)') '|                              TRANSPORT                                    |'
    write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
    write (stdout, *)

    if (index(transport_mode, 'bulk') > 0) then
      write (stdout, '(/1x,a/)') 'Calculation of Quantum Conductance and DoS: bulk mode'
      if (.not. tran_read_ht) then
        call hamiltonian_setup()
        call hamiltonian_get_hr()
        if (write_hr) call hamiltonian_write_hr()
        call tran_reduce_hr()
        call tran_cut_hr_one_dim()
        call tran_get_ht()
        if (write_xyz) call tran_write_xyz()
      end if
      call tran_bulk()
    end if

    if (index(transport_mode, 'lcr') > 0) then
      write (stdout, '(/1x,a/)') 'Calculation of Quantum Conductance and DoS: lead-conductor-lead mode'
      if (.not. tran_read_ht) then
        call hamiltonian_setup()
        call hamiltonian_get_hr()
        if (write_hr) call hamiltonian_write_hr()
        call tran_reduce_hr()
        call tran_cut_hr_one_dim()
        write (stdout, *) '------------------------- 2c2 Calculation Type: ------------------------------'
        write (stdout, *) ' '
        call tran_find_integral_signatures(signatures, num_G)
        call tran_lcr_2c2_sort(signatures, num_G, pl_warning)
        if (write_xyz) call tran_write_xyz()
        call tran_parity_enforce(signatures)
        call tran_lcr_2c2_build_ham(pl_warning)
      endif
      call tran_lcr()
    end if

    if (timing_level > 0) call io_stopwatch('tran: main', 2)

  end subroutine tran_main

  !==================================================================!
  subroutine tran_reduce_hr()
    !==================================================================!
    !
    ! reduce ham_r from 3-d to 1-d
    !
    use w90_constants, only: dp, eps8
    use w90_io, only: io_error, io_stopwatch, stdout
    use w90_parameters, only: one_dim_dir, real_lattice, num_wann, &
      mp_grid, timing_level
    use w90_hamiltonian, only: irvec, nrpts, ham_r

    implicit none

    integer :: ierr
    integer :: irvec_max, irvec_tmp(3), two_dim_vec(2)
    integer :: i, j
    integer :: i1, i2, i3, n1, nrpts_tmp, loop_rpt

    if (timing_level > 1) call io_stopwatch('tran: reduce_hr', 1)

    ! Find one_dim_vec which is parallel to one_dim_dir
    ! two_dim_vec - the other two lattice vectors
    j = 0
    do i = 1, 3
      if (abs(abs(real_lattice(one_dim_dir, i)) &
              - sqrt(dot_product(real_lattice(:, i), real_lattice(:, i)))) .lt. eps8) then
        one_dim_vec = i
        j = j + 1
      end if
    end do
    if (j .ne. 1) then
      write (stdout, '(i3,a)') j, ' : 1-D LATTICE VECTOR NOT DEFINED'
      call io_error('Error: 1-d lattice vector not defined in tran_reduce_hr')
    end if

    j = 0
    do i = 1, 3
      if (i .ne. one_dim_vec) then
        j = j + 1
        two_dim_vec(j) = i
      end if
    end do

    ! starting H matrix should include all W-S supercell where
    ! the center of the cell spans the full space of the home cell
    ! adding one more buffer layer when mp_grid(one_dim_vec) is an odd number

    !irvec_max = (mp_grid(one_dim_vec)+1)/2
    irvec_tmp = maxval(irvec, DIM=2) + 1
    irvec_max = irvec_tmp(one_dim_vec)
    nrpts_one_dim = 2*irvec_max + 1

    allocate (hr_one_dim(num_wann, num_wann, -irvec_max:irvec_max), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hr_one_dim in tran_reduce_hr')
    hr_one_dim = 0.0_dp

    ! check imaginary part
    write (stdout, '(1x,a,F12.6)') 'Maximum imaginary part of the real-space Hamiltonian: ', maxval(abs(aimag(ham_r)))

    ! select a subset of ham_r, where irvec is 0 along the two other lattice vectors

    nrpts_tmp = 0
    loop_n1: do n1 = -irvec_max, irvec_max
      do loop_rpt = 1, nrpts
        i1 = mod(n1 - irvec(one_dim_vec, loop_rpt), mp_grid(one_dim_vec))
        i2 = irvec(two_dim_vec(1), loop_rpt)
        i3 = irvec(two_dim_vec(2), loop_rpt)
        if (i1 .eq. 0 .and. i2 .eq. 0 .and. i3 .eq. 0) then
          nrpts_tmp = nrpts_tmp + 1
          hr_one_dim(:, :, n1) = real(ham_r(:, :, loop_rpt), dp)
          cycle loop_n1
        end if
      end do
    end do loop_n1

    if (nrpts_tmp .ne. nrpts_one_dim) then
      write (stdout, '(a)') 'FAILED TO EXTRACT 1-D HAMILTONIAN'
      call io_error('Error: cannot extract 1d hamiltonian in tran_reduce_hr')
    end if

    if (timing_level > 1) call io_stopwatch('tran: reduce_hr', 2)

    return

  end subroutine tran_reduce_hr

  !==================================================================!
  subroutine tran_cut_hr_one_dim()
    !==================================================================!
    !
    use w90_constants, only: dp
    use w90_io, only: io_stopwatch, stdout
    use w90_parameters, only: num_wann, mp_grid, timing_level, real_lattice, &
      hr_cutoff, dist_cutoff, dist_cutoff_mode, &
      one_dim_dir, length_unit, transport_mode, &
      tran_num_cell_ll, tran_num_ll, dist_cutoff_hc
    use w90_hamiltonian, only: wannier_centres_translated

    implicit none
    !
    integer :: irvec_max
    integer :: i, j, n1
    real(kind=dp) :: hr_max
    real(kind=dp) :: dist
    real(kind=dp) :: dist_vec(3)
    real(kind=dp) :: dist_ij_vec(3)
    real(kind=dp) :: shift_vec(3, -nrpts_one_dim/2:nrpts_one_dim/2)
    real(kind=dp) :: hr_tmp(num_wann, num_wann)

    !
    if (timing_level > 1) call io_stopwatch('tran: cut_hr_one_dim', 1)
    !
    irvec_max = nrpts_one_dim/2
    ! maximum possible dist_cutoff
    dist = real(mp_grid(one_dim_vec), dp)*abs(real_lattice(one_dim_dir, one_dim_vec))/2.0_dp

    if (dist_cutoff .gt. dist) then
      write (stdout, '(1x,a,1x,F10.5,1x,a)') 'dist_cutoff', dist_cutoff, trim(length_unit), 'is too large'
      dist_cutoff = dist
      ! aam_2012-04-13
      dist_cutoff_hc = dist
      write (stdout, '(4x,a,1x,F10.5,1x,a)') 'reset to', dist_cutoff, trim(length_unit)
    end if

    do n1 = -irvec_max, irvec_max
      shift_vec(:, n1) = real(n1, dp)*(real_lattice(:, one_dim_vec))
!           write(stdout,'(a,3f10.6)') 'shift_vec', shift_vec(:,n1)
    end do

    ! apply dist_cutoff first
    if (index(dist_cutoff_mode, 'one_dim') > 0) then
      do i = 1, num_wann
        do j = 1, num_wann
          dist_ij_vec(one_dim_dir) = wannier_centres_translated(one_dim_dir, i) - wannier_centres_translated(one_dim_dir, j)
          do n1 = -irvec_max, irvec_max
            dist_vec(one_dim_dir) = dist_ij_vec(one_dim_dir) + shift_vec(one_dim_dir, n1)
            !
            !MS: Add special case for lcr: We must not cut the elements that are within
            !    dist_cutoff under PBC's (single kpt assumed) in order to build
            !    hamiltonians correctly in tran_2c2_build_hams
            !
            if ((index(transport_mode, 'lcr') > 0) .and. &
                !~                    (tran_num_cell_ll .eq. 1)        .and. &
                (abs(dist_vec(one_dim_dir)) .gt. dist_cutoff)) then
              ! Move to right
              dist_vec(one_dim_dir) = dist_ij_vec(one_dim_dir) + real_lattice(one_dim_dir, one_dim_vec)
              ! Move to left
              if (abs(dist_vec(one_dim_dir)) .gt. dist_cutoff) &
                dist_vec(one_dim_dir) = dist_ij_vec(one_dim_dir) - real_lattice(one_dim_dir, one_dim_vec)
            endif
            !
            !end MS
            !
            dist = abs(dist_vec(one_dim_dir))
            if (dist .gt. dist_cutoff) hr_one_dim(j, i, n1) = 0.0_dp
          end do
        end do
      end do
    else
      do i = 1, num_wann
        do j = 1, num_wann
          dist_ij_vec(:) = wannier_centres_translated(:, i) - wannier_centres_translated(:, j)
          do n1 = -irvec_max, irvec_max
            dist_vec(:) = dist_ij_vec(:) + shift_vec(:, n1)
            dist = sqrt(dot_product(dist_vec, dist_vec))
            !
            ! MS: Special case (as above) equivalent for alternate definition of cut off
            !
            if ((index(transport_mode, 'lcr') > 0) .and. &
                !~                   (tran_num_cell_ll .eq. 1)         .and. &
                (dist .gt. dist_cutoff)) then
              ! Move to right
              dist_vec(:) = dist_ij_vec(:) + real_lattice(:, one_dim_vec)
              dist = sqrt(dot_product(dist_vec, dist_vec))
              ! Move to left
              if (dist .gt. dist_cutoff) then
                dist_vec(:) = dist_ij_vec(:) - real_lattice(:, one_dim_vec)
                dist = sqrt(dot_product(dist_vec, dist_vec))
              endif
            endif
            !
            ! End MS
            !
            if (dist .gt. dist_cutoff) hr_one_dim(j, i, n1) = 0.0_dp
          end do
        end do
      end do
    end if

    ! output maximum to check a decay of H as a function of lattice vector R
    write (stdout, '(/1x,a78)') repeat('-', 78)
    write (stdout, '(1x,4x,a)') &
      'Maximum real part of the real-space Hamiltonian at each lattice point'
    write (stdout, '(1x,8x,a62)') repeat('-', 62)
    write (stdout, '(1x,11x,a,11x,a)') 'Lattice point R', 'Max |H_ij(R)|'
    ! calculate number of units inside a principal layer
    num_pl = 0
    do n1 = -irvec_max, irvec_max
      hr_tmp(:, :) = abs(hr_one_dim(:, :, n1))
      hr_max = maxval(hr_tmp)
      if (hr_max .gt. hr_cutoff) then
        if (abs(n1) .gt. num_pl) num_pl = abs(n1)
      else
        hr_one_dim(:, :, n1) = 0.0_dp
      end if
      write (stdout, '(1x,9x,5x,I5,5x,12x,F12.6)') n1, hr_max
    end do
    write (stdout, '(1x,8x,a62)') repeat('-', 62)
    if (index(transport_mode, 'lcr') > 0) then
      write (stdout, '(/1x,a,I6)') 'Number of unit cells inside the principal layer:', tran_num_cell_ll
      write (stdout, '(1x,a,I6)') 'Number of Wannier Functions inside the principal layer:', tran_num_ll
    elseif (index(transport_mode, 'bulk') > 0) then
      write (stdout, '(/1x,a,I6)') 'Number of unit cells inside the principal layer:', num_pl
      write (stdout, '(1x,a,I6)') 'Number of Wannier Functions inside the principal layer:', num_pl*num_wann
    endif
    ! apply hr_cutoff to each element inside the principal layer
    do n1 = -num_pl, num_pl
      do i = 1, num_wann
        do j = 1, num_wann
          if (abs(hr_one_dim(j, i, n1)) .lt. hr_cutoff) hr_one_dim(j, i, n1) = 0.0_dp
        end do
      end do
    end do

    if (timing_level > 1) call io_stopwatch('tran: cut_hr_one_dim', 2)

    return

  end subroutine tran_cut_hr_one_dim

  !==================================================================!
  subroutine tran_get_ht()
    !==================================================================!
    !  construct h00 and h01
    !==================================================================!
    !
    use w90_constants, only: dp
    use w90_io, only: io_error, io_stopwatch, seedname, io_date, &
      io_file_unit
    use w90_parameters, only: num_wann, tran_num_bb, tran_write_ht, &
      nfermi, fermi_energy_list, timing_level
    !
    implicit none
    !
    integer :: ierr, file_unit
    integer :: i, j, n1, im, jm
    character(len=9)   :: cdate, ctime
    !
    if (timing_level > 1) call io_stopwatch('tran: get_ht', 1)
    !
    if (nfermi > 1) call io_error("Error in tran_get_ht: nfermi>1. " &
                                  //"Set the fermi level using the input parameter 'fermi_evel'")
    !
    !
    tran_num_bb = num_pl*num_wann
    !
    allocate (hB0(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hB0 in tran_get_ht')
    allocate (hB1(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hB1 in tran_get_ht')
    !
    hB0 = 0.0_dp
    hB1 = 0.0_dp
    !
    ! h00
    do j = 0, num_pl - 1
      do i = 0, num_pl - 1
        n1 = i - j
        im = i*num_wann
        jm = j*num_wann
        hB0(jm + 1:jm + num_wann, im + 1:im + num_wann) = hr_one_dim(:, :, n1)
      end do
    end do

    ! h01
    do j = 1, num_pl
      do i = 0, j - 1
        n1 = i - (j - 1) + num_pl
        im = i*num_wann
        jm = (j - 1)*num_wann
        hB1(jm + 1:jm + num_wann, im + 1:im + num_wann) = hr_one_dim(:, :, n1)
      end do
    end do

    ! shift by fermi_energy
    do i = 1, tran_num_bb
      hB0(i, i) = hB0(i, i) - fermi_energy_list(1)
    end do

    if (tran_write_ht) then

      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_htB.dat', status='unknown', form='formatted', action='write')

      call io_date(cdate, ctime)
      write (file_unit, *) 'written on '//cdate//' at '//ctime ! Date and time
      write (file_unit, '(I6)') tran_num_bb
      write (file_unit, '(6F12.6)') ((hB0(j, i), j=1, tran_num_bb), i=1, tran_num_bb)
      write (file_unit, '(I6)') tran_num_bb
      write (file_unit, '(6F12.6)') ((hB1(j, i), j=1, tran_num_bb), i=1, tran_num_bb)

      close (file_unit)

    end if

    if (timing_level > 1) call io_stopwatch('tran: get_ht', 2)

    return

  end subroutine tran_get_ht

  !==================================================================!
  subroutine tran_bulk()
    !==================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_1, cmplx_i, pi
    use w90_io, only: io_error, io_stopwatch, seedname, io_date, &
      io_file_unit, stdout
    use w90_parameters, only: tran_num_bb, tran_read_ht, &
      tran_win_min, tran_win_max, tran_energy_step, &
      timing_level

    implicit none

    integer :: qc_unit, dos_unit
    integer :: ierr
    integer :: n_e, n, i
    real(kind=dp) ::  qc, dos
    real(kind=dp) ::  e_scan
    complex(kind=dp) :: e_scan_cmp
    complex(kind=dp), allocatable, dimension(:, :) :: tot, tott
    complex(kind=dp), allocatable, dimension(:, :) :: g_B, gR, gL
    complex(kind=dp), allocatable, dimension(:, :) :: sLr, sRr
    complex(kind=dp), allocatable, dimension(:, :) :: s1, s2, c1
    character(len=50) :: filename
    character(len=9)  :: cdate, ctime

    if (timing_level > 1) call io_stopwatch('tran: bulk', 1)

    allocate (tot(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tot in tran_bulk')
    allocate (tott(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tott in tran_bulk')
    allocate (g_B(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating g_B in tran_bulk')
    allocate (gL(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating gL in tran_bulk')
    allocate (gR(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating gR in tran_bulk')
    allocate (sLr(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating sLr in tran_bulk')
    allocate (sRr(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating sRr in tran_bulk')
    allocate (s1(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating s1 in tran_bulk')
    allocate (s2(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating s2 in tran_bulk')
    allocate (c1(tran_num_bb, tran_num_bb), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating c1 in tran_bulk')

    call io_date(cdate, ctime)

    qc_unit = io_file_unit()
    open (qc_unit, file=trim(seedname)//'_qc.dat', status='unknown', &
          form='formatted', action='write')
    write (qc_unit, *) '## written on '//cdate//' at '//ctime ! Date and time

    dos_unit = io_file_unit()
    open (dos_unit, file=trim(seedname)//'_dos.dat', status='unknown', &
          form='formatted', action='write')
    write (dos_unit, *) '## written on '//cdate//' at '//ctime ! Date and time

    !   set up the layer hamiltonians

    if (tran_read_ht) then
      allocate (hB0(tran_num_bb, tran_num_bb), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating hB0 in tran_bulk')
      allocate (hB1(tran_num_bb, tran_num_bb), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating hB1 in tran_bulk')
      filename = trim(seedname)//'_htB.dat'
      call tran_read_htX(tran_num_bb, hB0, hB1, filename)
    end if

    !   loop over the energies

    n_e = floor((tran_win_max - tran_win_min)/tran_energy_step) + 1

    write (stdout, '(/1x,a)', advance='no') 'Calculating quantum&
        & conductance and density of states...'

    do n = 1, n_e
      e_scan = tran_win_min + real(n - 1, dp)*tran_energy_step

!       if (mod(n,nint(0.1*n_e)).eq.0) write(stdout,'(a)',advance='no') '.'

      ! compute conductance according to Fisher and Lee
      ! retarded Green

      e_scan_cmp = e_scan + eta
      call tran_transfer(tot, tott, hB0, hB1, e_scan_cmp, tran_num_bb)
      call tran_green(tot, tott, hB0, hB1, e_scan, g_B, 0, 1, tran_num_bb)

      ! compute S_Lr and S_Rr

      c1(:, :) = cmplx(hB1(:, :), kind=dp)

      ! Self-energy (Sigma_L^r) : sLr = (hB1)^+ * tott
      ! Self-energy (Sigma_R^r) : sRr = (hB1)   * tot
      sLr = cmplx_0
      sRr = cmplx_0
      call ZGEMM('C', 'N', tran_num_bb, tran_num_bb, tran_num_bb, cmplx_1, c1, &
                 tran_num_bb, tott, tran_num_bb, cmplx_0, sLr, tran_num_bb)
      call ZGEMM('N', 'N', tran_num_bb, tran_num_bb, tran_num_bb, cmplx_1, c1, &
                 tran_num_bb, tot, tran_num_bb, cmplx_0, sRr, tran_num_bb)

      ! Gamma_L = i(Sigma_L^r-Sigma_L^a)
      gL = cmplx_i*(sLr - conjg(transpose(sLr)))
      ! Gamma_R = i(Sigma_R^r-Sigma_R^a)
      gR = cmplx_i*(sRr - conjg(transpose(sRr)))

      s1 = cmplx_0
      s2 = cmplx_0
      c1 = cmplx_0
      ! s1 = Gamma_L * g_B^r
      call ZGEMM('N', 'N', tran_num_bb, tran_num_bb, tran_num_bb, cmplx_1, gL, &
                 tran_num_bb, g_B, tran_num_bb, cmplx_0, s1, tran_num_bb)
      ! s2 = Gamma_L * g_B^r * Gamma_R
      call ZGEMM('N', 'N', tran_num_bb, tran_num_bb, tran_num_bb, cmplx_1, s1, &
                 tran_num_bb, gR, tran_num_bb, cmplx_0, s2, tran_num_bb)
      ! c1 = Gamma_L * g_B^r * Gamma_R * g_B^a
      call ZGEMM('N', 'C', tran_num_bb, tran_num_bb, tran_num_bb, cmplx_1, s2, &
                 tran_num_bb, g_B, tran_num_bb, cmplx_0, c1, tran_num_bb)

      qc = 0.0_dp
      do i = 1, tran_num_bb
        qc = qc + real(c1(i, i), dp)
      end do
!       write(qc_unit,'(f12.6,f15.6)') e_scan, qc
      write (qc_unit, '(f15.9,f18.9)') e_scan, qc

      dos = 0.0_dp
      do i = 1, tran_num_bb
        dos = dos - aimag(g_B(i, i))
      end do
      dos = dos/pi
      write (dos_unit, '(f15.9,f18.9)') e_scan, dos

    end do

    write (stdout, '(a/)') ' done'

    close (qc_unit)
    close (dos_unit)

    deallocate (c1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating c1 in tran_bulk')
    deallocate (s2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating s2 in tran_bulk')
    deallocate (s1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating s1 in tran_bulk')
    deallocate (sRr, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating sRr in tran_bulk')
    deallocate (sLr, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating sLr in tran_bulk')
    deallocate (gR, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating gR in tran_bulk')
    deallocate (gL, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating gL in tran_bulk')
    deallocate (g_B, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating g_B in tran_bulk')
    deallocate (tott, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tott in tran_bulk')
    deallocate (tot, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tot in tran_bulk')

    if (timing_level > 1) call io_stopwatch('tran: bulk', 2)

    return

  end subroutine tran_bulk

  !==================================================================!
  subroutine tran_lcr()
    !==================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_1, cmplx_i, pi
    use w90_io, only: io_error, io_stopwatch, seedname, io_date, &
      stdout, io_file_unit
    use w90_parameters, only: tran_num_ll, tran_num_rr, tran_num_cc, tran_num_lc, &
      tran_num_cr, tran_num_bandc, &
      tran_win_min, tran_win_max, tran_energy_step, &
      tran_use_same_lead, timing_level, tran_read_ht

    implicit none

    integer :: qc_unit, dos_unit
    integer :: ierr
    integer :: KL, KU, KC
    integer :: n_e, n, i, j, k, info
    integer, allocatable :: ipiv(:)
    real(kind=dp) ::  qc, dos
    real(kind=dp) ::  e_scan
    real(kind=dp), allocatable, dimension(:, :) :: hCband
    complex(kind=dp) :: e_scan_cmp
    complex(kind=dp), allocatable, dimension(:, :) :: hLC_cmp, hCR_cmp, &
                                                      totL, tottL, totR, tottR, &
                                                      g_surf_L, g_surf_R, g_C, g_C_inv, &
                                                      gR, gL, sLr, sRr, s1, s2, c1, c2
    character(len=50) :: filename
    character(len=9)  :: cdate, ctime

    if (timing_level > 1) call io_stopwatch('tran: lcr', 1)

    call io_date(cdate, ctime)

    qc_unit = io_file_unit()
    open (qc_unit, file=trim(seedname)//'_qc.dat', status='unknown', &
          form='formatted', action='write')
    write (qc_unit, *) '## written on '//cdate//' at '//ctime ! Date and time

    dos_unit = io_file_unit()
    open (dos_unit, file=trim(seedname)//'_dos.dat', status='unknown', &
          form='formatted', action='write')
    write (dos_unit, *) '## written on '//cdate//' at '//ctime ! Date and time

    KL = max(tran_num_lc, tran_num_cr, tran_num_bandc) - 1
    KU = KL
    KC = max(tran_num_lc, tran_num_cr)

    allocate (hCband(2*KL + KU + 1, tran_num_cc), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hCband in tran_lcr')
    allocate (hLC_cmp(tran_num_ll, tran_num_lc), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hLC_cmp in tran_lcr')
    allocate (hCR_cmp(tran_num_cr, tran_num_rr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hCR_cmp in tran_lcr')

    !If construct used only when reading matrices from file
    if (tran_read_ht) then
      allocate (hL0(tran_num_ll, tran_num_ll), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating hL0 in tran_lcr')
      allocate (hL1(tran_num_ll, tran_num_ll), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating hL1 in tran_lcr')
      allocate (hC(tran_num_cc, tran_num_cc), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating hC in tran_lcr')
      allocate (hLC(tran_num_ll, tran_num_lc), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating hLC in tran_lcr')
      allocate (hCR(tran_num_cr, tran_num_rr), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating hCR in tran_lcr')

      filename = trim(seedname)//'_htL.dat'
      call tran_read_htX(tran_num_ll, hL0, hL1, filename)

      if (.not. tran_use_same_lead) then
        allocate (hR0(tran_num_rr, tran_num_rr), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating hR0 in tran_lcr')
        allocate (hR1(tran_num_rr, tran_num_rr), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating hR1 in tran_lcr')
        filename = trim(seedname)//'_htR.dat'
        call tran_read_htX(tran_num_rr, hR0, hR1, filename)
      end if

      filename = trim(seedname)//'_htC.dat'
      call tran_read_htC(tran_num_cc, hC, filename)
      filename = trim(seedname)//'_htLC.dat'
      call tran_read_htXY(tran_num_ll, tran_num_lc, hLC, filename)
      filename = trim(seedname)//'_htCR.dat'
      call tran_read_htXY(tran_num_cr, tran_num_rr, hCR, filename)
    endif

    !  Banded matrix H_C  :  save memory !
    do j = 1, tran_num_cc
      do i = max(1, j - KU), min(tran_num_cc, j + KL)
        hCband(KL + KU + 1 + i - j, j) = hC(i, j)
      end do
    end do
    deallocate (hC, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating hC in tran_lcr')

    !  H_LC : to a complex matrix
    hLC_cmp(:, :) = cmplx(hLC(:, :), kind=dp)
    deallocate (hLC, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating hLC in tran_lcr')

    !  H_CR : to a complex matrix
    hCR_cmp(:, :) = cmplx(hCR(:, :), kind=dp)
    deallocate (hCR, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating hCR in tran_lcr')

    allocate (totL(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating totL in tran_lcr')
    allocate (tottL(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tottL in tran_lcr')
    if (.not. tran_use_same_lead) then
      allocate (totR(tran_num_rr, tran_num_rr), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating totR in tran_lcr')
      allocate (tottR(tran_num_rr, tran_num_rr), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating tottR in tran_lcr')
    end if
    allocate (g_surf_L(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating g_surf_L in tran_lcr')
    allocate (g_surf_R(tran_num_rr, tran_num_rr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating g_surf_R in tran_lcr')
    allocate (g_C_inv(2*KL + KU + 1, tran_num_cc), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating g_C_inv in tran_lcr')
    allocate (g_C(tran_num_cc, tran_num_cc), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating g_C in tran_lcr')
    allocate (sLr(tran_num_lc, tran_num_lc), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating sLr in tran_lcr')
    allocate (sRr(tran_num_cr, tran_num_cr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating sRr in tran_lcr')
    allocate (gL(tran_num_lc, tran_num_lc), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating gL in tran_lcr')
    allocate (gR(tran_num_cr, tran_num_cr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating gR in tran_lcr')
    allocate (c1(tran_num_lc, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating c1 in tran_lcr')
    allocate (c2(tran_num_cr, tran_num_rr), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating c2 in tran_lcr')
    allocate (s1(KC, KC), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating s1 in tran_lcr')
    allocate (s2(KC, KC), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating s2 in tran_lcr')
    allocate (ipiv(tran_num_cc), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ipiv in tran_lcr')

    !  Loop over the energies
    n_e = floor((tran_win_max - tran_win_min)/tran_energy_step) + 1

    write (stdout, '(/1x,a)', advance='no') 'Calculating quantum conductance and density of states...'

    do n = 1, n_e

      e_scan = tran_win_min + real(n - 1, dp)*tran_energy_step

      !    compute conductance according to Fisher and Lee
      !    compute self-energies following Datta

      e_scan_cmp = e_scan + eta

      ! Surface green function for the left lead : g_surf_L
      call tran_transfer(totL, tottL, hL0, hL1, e_scan_cmp, tran_num_ll)
      call tran_green(totL, tottL, hL0, hL1, e_scan, g_surf_L, -1, 1, tran_num_ll)

      ! Self-energy (Sigma_L) : sLr = (hLC_cmp)^+ * g_surf_L * hLC_cmp
      c1 = cmplx_0
      sLr = cmplx_0
      call ZGEMM('C', 'N', tran_num_lc, tran_num_ll, tran_num_ll, cmplx_1, &
                 hLC_cmp, tran_num_ll, g_surf_L, tran_num_ll, cmplx_0, c1, tran_num_lc)
      call ZGEMM('N', 'N', tran_num_lc, tran_num_lc, tran_num_ll, cmplx_1, &
                 c1, tran_num_lc, hLC_cmp, tran_num_ll, cmplx_0, sLr, tran_num_lc)

      ! Surface green function for the right lead : g_surf_R
      if (tran_use_same_lead) then
        call tran_green(totL, tottL, hL0, hL1, e_scan, g_surf_R, 1, 1, tran_num_rr)
      else
        call tran_transfer(totR, tottR, hR0, hR1, e_scan_cmp, tran_num_rr)
        call tran_green(totR, tottR, hR0, hR1, e_scan, g_surf_R, 1, 1, tran_num_rr)
      end if

      ! Self-energy (Sigma_R) : sRr = hCR_cmp * g_surf_R * (hCR_cmp)^+
      c2 = cmplx_0
      sRr = cmplx_0
      call ZGEMM('N', 'N', tran_num_cr, tran_num_rr, tran_num_rr, cmplx_1, &
                 hCR_cmp, tran_num_cr, g_surf_R, tran_num_rr, cmplx_0, c2, tran_num_cr)
      call ZGEMM('N', 'C', tran_num_cr, tran_num_cr, tran_num_rr, cmplx_1, &
                 c2, tran_num_cr, hCR_cmp, tran_num_cr, cmplx_0, sRr, tran_num_cr)

      ! g_C^-1 = -H
      g_C_inv(:, :) = cmplx(-hCband(:, :), kind=dp)

      ! g_C^-1 = -H - Sigma_L^r
      do j = 1, tran_num_lc
        do i = max(1, j - KU), min(tran_num_lc, j + KL)
          g_C_inv(KL + KU + 1 + i - j, j) = g_C_inv(KL + KU + 1 + i - j, j) - sLr(i, j)
        end do
      end do

      ! g_C^-1 = -H - Sigma_L^r - Sigma_R^r
      do j = (tran_num_cc - tran_num_cr) + 1, tran_num_cc
        do i = max((tran_num_cc - tran_num_cr) + 1, j - (tran_num_cr - 1)), min(tran_num_cc, j + (tran_num_cr - 1))
          g_C_inv(KL + KU + 1 + i - j, j) = &
            g_C_inv(KL + KU + 1 + i - j, j) - &
            sRr(i - (tran_num_cc - tran_num_cr), j - (tran_num_cc - tran_num_cr))
        end do
      end do

      ! g_C^-1 = eI - H - Sigma_L^r - Sigma_R^r
      do i = 1, tran_num_cc
        g_C_inv(KL + KU + 1, i) = e_scan + g_C_inv(KL + KU + 1, i)
      end do

      ! invert g_C^-1 => g_C
      g_C = cmplx_0
      do i = 1, tran_num_cc
        g_C(i, i) = cmplx_1
      end do

      call ZGBSV(tran_num_cc, KL, KU, tran_num_cc, g_C_inv, 2*KL + KU + 1, ipiv, g_C, tran_num_cc, info)
      if (info .ne. 0) then
        write (stdout, *) 'ERROR: IN ZGBSV IN tran_lcr, INFO=', info
        call io_error('tran_lcr: problem in ZGBSV')
      end if

      ! Gamma_L = i(Sigma_L^r-Sigma_L^a)
      gL = cmplx_i*(sLr - conjg(transpose(sLr)))

      ! s1 = Gamma_L * g_C^r
      s1 = cmplx_0
      do j = 1, KC
        do i = 1, tran_num_lc
          do k = 1, tran_num_lc
            s1(i, j) = s1(i, j) + gL(i, k)*g_C(k, j + (tran_num_cc - KC))
          end do
        end do
      end do

      ! Gamma_R = i(Sigma_R^r-Sigma_R^a)
      gR = cmplx_i*(sRr - conjg(transpose(sRr)))

      ! s2 = Gamma_R * g_C^a
      s2 = cmplx_0
      do j = 1, KC
        do i = 1, tran_num_cr
          do k = 1, tran_num_cr
            s2(i + (KC - tran_num_cr), j) = s2(i + (KC - tran_num_cr), j) &
                                            + gR(i, k)*conjg(g_C(j, k + (tran_num_cc - tran_num_cr)))
          end do
        end do
      end do

      qc = 0.0_dp
      do i = 1, KC
        do j = 1, KC
          qc = qc + real(s1(i, j)*s2(j, i), dp)
        end do
      end do
      write (qc_unit, '(f15.9,f18.9)') e_scan, qc

      ! compute density of states for the conductor layer

      dos = 0.0_dp
      do i = 1, tran_num_cc
        dos = dos - aimag(g_C(i, i))
      end do
      dos = dos/pi
      write (dos_unit, '(f15.9,f18.9)') e_scan, dos

    end do

    write (stdout, '(a)') ' done'

    close (qc_unit)
    close (dos_unit)

    deallocate (ipiv, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ipiv in tran_lcr')
    deallocate (s2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating s2 in tran_lcr')
    deallocate (s1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating s1 in tran_lcr')
    deallocate (c2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating c2 in tran_lcr')
    deallocate (c1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating c1 in tran_lcr')
    deallocate (gR, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating gR in tran_lcr')
    deallocate (gL, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating gL in tran_lcr')
    deallocate (sRr, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating sRr in tran_lcr')
    deallocate (sLr, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating sLr in tran_lcr')
    deallocate (g_C, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating g_C in tran_lcr')
    deallocate (g_C_inv, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating g_C_inv in tran_lcr')
    deallocate (g_surf_R, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating g_surf_R in tran_lcr')
    deallocate (g_surf_L, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating g_surf_L in tran_lcr')
    if (allocated(tottR)) deallocate (tottR, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tottR in tran_lcr')
    if (allocated(totR)) deallocate (totR, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating totR in tran_lcr')
    deallocate (tottL, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tottL in tran_lcr')
    deallocate (totL, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating totL in tran_lcr')
    deallocate (hCR_cmp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating hCR_cmp in tran_lcr')
    deallocate (hLC_cmp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating hLC_cmp in tran_lcr')
    deallocate (hCband, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating hCband in tran_lcr')

    if (timing_level > 1) call io_stopwatch('tran: lcr', 2)

    return

  end subroutine tran_lcr

  !==================================================================!
  subroutine tran_transfer(tot, tott, h_00, h_01, e_scan_cmp, nxx)
    !==================================================================!
    !                                                                  !
    ! iterative construction of the transfer matrix                    !
    ! as Lopez-Sancho^2&Rubio, J.Phys.F:Met.Phys., v.14, 1205 (1984)   !
    ! and ibid. v.15, 851 (1985)                                       !
    !                                                                  !
    !===================================================================

    use w90_constants, only: dp, cmplx_0, cmplx_1, eps7
    use w90_io, only: stdout, io_error

    implicit none

    integer, intent(in) :: nxx
    complex(kind=dp), intent(in) ::  e_scan_cmp
    complex(kind=dp), intent(out) ::  tot(nxx, nxx)
    complex(kind=dp), intent(out) ::  tott(nxx, nxx)
    real(kind=dp), intent(in) :: h_00(nxx, nxx)
    real(kind=dp), intent(in) :: h_01(nxx, nxx)
    !
    integer  :: ierr, info
    integer  :: i, j, n, nxx2
    integer, allocatable :: ipiv(:)
    real(kind=dp) :: conver, conver2
    complex(kind=dp), allocatable, dimension(:, :) :: tsum, tsumt
    complex(kind=dp), allocatable, dimension(:, :) :: t11, t12
    complex(kind=dp), allocatable, dimension(:, :) :: s1, s2
    complex(kind=dp), allocatable, dimension(:, :, :) :: tau, taut

    allocate (ipiv(nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ipiv in tran_transfer')
    allocate (tsum(nxx, nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tsum in tran_transfer')
    allocate (tsumt(nxx, nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tsumt in tran_transfer')
    allocate (t11(nxx, nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating t11 in tran_transfer')
    allocate (t12(nxx, nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating t12 in tran_transfer')
    allocate (s1(nxx, nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating s1 in tran_transfer')
    allocate (s2(nxx, nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating s2 in tran_transfer')
    allocate (tau(nxx, nxx, 2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tau in tran_transfer')
    allocate (taut(nxx, nxx, 2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating taut in tran_transfer')

    nxx2 = nxx*nxx

    tot = cmplx_0
    tott = cmplx_0

    ! construction of the transfer matrix
    ! t12 = e - h_00
    t12(:, :) = cmplx(-h_00(:, :), kind=dp)
    do i = 1, nxx
      t12(i, i) = e_scan_cmp + t12(i, i)
    end do

    ! compute (e - h_00)^-1 and store it in t11
    t11 = cmplx_0
    do i = 1, nxx
      t11(i, i) = cmplx_1
    end do

    ! inverse of t12 -> t11
    call ZGESV(nxx, nxx, t12, nxx, ipiv, t11, nxx, info)
    if (info .ne. 0) then
      write (stdout, *) 'ERROR:  IN ZGESV IN tran_transfer, INFO=', info
      call io_error('tran_transfer: problem in ZGESV 1')
    end if

    ! compute intermediate t-matrices (defined as tau(nxx,nxx,niter)
    ! and taut(...)):
    tau = cmplx_0
    taut = cmplx_0

    ! t_0:
    t12(:, :) = cmplx(h_01(:, :), kind=dp)

    !  tau  = ( e - H_00 )^-1 * H_01^+
    call ZGEMM('N', 'C', nxx, nxx, nxx, cmplx_1, t11, nxx, t12, nxx, cmplx_0, tau(1, 1, 1), nxx)
    !  taut = ( e - H_00 )^-1 * H_01
    call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, t11, nxx, t12, nxx, cmplx_0, taut(1, 1, 1), nxx)

    !   initialize T:
    tot(:, :) = tau(:, :, 1)
    tsum(:, :) = taut(:, :, 1)

    !   initialize T^bar:
    tott(:, :) = taut(:, :, 1)
    tsumt(:, :) = tau(:, :, 1)

    !   main loop:
    do n = 1, nterx

      t11 = cmplx_0
      t12 = cmplx_0

      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, tau(1, 1, 1), nxx, taut(1, 1, 1), nxx, cmplx_0, t11, nxx)
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, taut(1, 1, 1), nxx, tau(1, 1, 1), nxx, cmplx_0, t12, nxx)

      s1(:, :) = -t11(:, :) - t12(:, :)
      do i = 1, nxx
        s1(i, i) = cmplx_1 + s1(i, i)
      end do

      s2 = cmplx_0
      do i = 1, nxx
        s2(i, i) = cmplx_1
      end do

      call ZGESV(nxx, nxx, s1, nxx, ipiv, s2, nxx, info)
      if (info .ne. 0) then
        write (stdout, *) 'ERROR:  IN ZGESV IN tran_transfer, INFO=', info
        call io_error('tran_transfer: problem in ZGESV 2')
      end if

      t11 = cmplx_0
      t12 = cmplx_0

      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, tau(1, 1, 1), nxx, tau(1, 1, 1), nxx, cmplx_0, t11, nxx)
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, taut(1, 1, 1), nxx, taut(1, 1, 1), nxx, cmplx_0, t12, nxx)
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, s2, nxx, t11, nxx, cmplx_0, tau(1, 1, 2), nxx)
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, s2, nxx, t12, nxx, cmplx_0, taut(1, 1, 2), nxx)

      !   put the transfer matrices together

      t11 = cmplx_0
      s1 = cmplx_0

      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, tsum, nxx, tau(1, 1, 2), nxx, cmplx_0, t11, nxx)
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, tsum, nxx, taut(1, 1, 2), nxx, cmplx_0, s1, nxx)
      call ZCOPY(nxx2, t11, 1, s2, 1)
      call ZAXPY(nxx2, cmplx_1, tot, 1, s2, 1)

      tot(:, :) = s2(:, :)
      tsum(:, :) = s1(:, :)

      t11 = cmplx_0
      s1 = cmplx_0

      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, tsumt, nxx, taut(1, 1, 2), nxx, cmplx_0, t11, nxx)
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, tsumt, nxx, tau(1, 1, 2), nxx, cmplx_0, s1, nxx)
      call ZCOPY(nxx2, t11, 1, s2, 1)
      call ZAXPY(nxx2, cmplx_1, tott, 1, s2, 1)

      tott(:, :) = s2(:, :)
      tsumt(:, :) = s1(:, :)

      tau(:, :, 1) = tau(:, :, 2)
      taut(:, :, 1) = taut(:, :, 2)

      ! convergence check on the t-matrices

      conver = 0.0_dp
      conver2 = 0.0_dp

      do j = 1, nxx
        do i = 1, nxx
          conver = conver + sqrt(real(tau(i, j, 2), dp)**2 + aimag(tau(i, j, 2))**2)
          conver2 = conver2 + sqrt(real(taut(i, j, 2), dp)**2 + aimag(taut(i, j, 2))**2)
        end do
      end do

      if (conver .lt. eps7 .and. conver2 .lt. eps7) return
    end do

    if (conver .gt. eps7 .or. conver2 .gt. eps7) &
      call io_error('Error in converging transfer matrix in tran_transfer')

    deallocate (taut, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating taut in tran_transfer')
    deallocate (tau, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tau in tran_transfer')
    deallocate (s2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating s2 in tran_transfer')
    deallocate (s1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating s1 in tran_transfer')
    deallocate (t12, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating t12 in tran_transfer')
    deallocate (t11, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating t11 in tran_transfer')
    deallocate (tsumt, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tsumt in tran_transfer')
    deallocate (tsum, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tsum in tran_transfer')
    deallocate (ipiv, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ipiv in tran_transfer')

    return

  end subroutine tran_transfer

  !==================================================================!
  subroutine tran_green(tot, tott, h_00, h_01, e_scan, g, igreen, invert, nxx)
    !==================================================================!
    !   construct green's functions
    !
    !   igreen = -1  left surface
    !   igreen =  1  right surface
    !   igreen =  0  bulk

    !   invert = 0 computes g^-1
    !   invert = 1 computes g^-1 and g
    !==================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_io, only: stdout, io_error

    implicit none

    integer, intent(in) :: nxx
    integer, intent(in) :: igreen
    integer, intent(in) :: invert
    real(kind=dp), intent(in) :: e_scan
    real(kind=dp), intent(in) :: h_00(nxx, nxx), h_01(nxx, nxx)
    complex(kind=dp), intent(in) :: tot(nxx, nxx), tott(nxx, nxx)
    complex(kind=dp), intent(out) :: g(nxx, nxx)

    integer :: ierr, info
    integer :: i
    integer, allocatable :: ipiv(:)
    complex(kind=dp), allocatable, dimension(:, :) :: g_inv, eh_00
    complex(kind=dp), allocatable, dimension(:, :) :: s1, s2, c1

    allocate (ipiv(nxx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ipiv in tran_green')
    allocate (g_inv(nxx, nxx))
    if (ierr /= 0) call io_error('Error in allocating g_inv in tran_green')
    allocate (eh_00(nxx, nxx))
    if (ierr /= 0) call io_error('Error in allocating eh_00 in tran_green')
    allocate (c1(nxx, nxx))
    if (ierr /= 0) call io_error('Error in allocating c1 in tran_green')
    allocate (s1(nxx, nxx))
    if (ierr /= 0) call io_error('Error in allocating s1 in tran_green')
    allocate (s2(nxx, nxx))
    if (ierr /= 0) call io_error('Error in allocating s2 in tran_green')

    c1(:, :) = cmplx(h_01(:, :), kind=dp)

    select case (igreen)

    case (1)

      ! construct the surface green's function g00

      s1 = cmplx_0
      ! s1 = H_01 * T
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, c1, nxx, tot, nxx, cmplx_0, s1, nxx)

      ! eh_00 =  -H_00 - H_01*T
      eh_00(:, :) = cmplx(-h_00(:, :), kind=dp) - s1(:, :)
      ! eh_00 = e_scan -H_00 - H_01*T
      do i = 1, nxx
        eh_00(i, i) = cmplx(e_scan, kind=dp) + eh_00(i, i)
      end do

      g_inv(:, :) = eh_00(:, :)

      ! identity
      g = cmplx_0
      do i = 1, nxx
        g(i, i) = cmplx_1
      end do

      if (invert .eq. 1) then
        call ZGESV(nxx, nxx, eh_00, nxx, ipiv, g, nxx, info)
        if (info .ne. 0) then
          write (stdout, *) 'ERROR:  IN ZGESV IN tran_green, INFO=', info
          call io_error('tran_green: problem in ZGESV 1')
        end if
      end if

    case (-1)

      !  construct the dual surface green's function gbar00

      s1 = cmplx_0
      ! s1 = H_01^+ * T^bar
      call ZGEMM('C', 'N', nxx, nxx, nxx, cmplx_1, c1, nxx, tott, nxx, cmplx_0, s1, nxx)

      ! s1 = -H_00 - H_01^+ * T^bar
      eh_00(:, :) = cmplx(-h_00(:, :), kind=dp) - s1(:, :)
      ! s1 = e_scan - H_00 - H_01^+ * T^bar
      do i = 1, nxx
        eh_00(i, i) = cmplx(e_scan, kind=dp) + eh_00(i, i)
      end do

      g_inv(:, :) = eh_00(:, :)

      ! identity
      g = cmplx_0
      do i = 1, nxx
        g(i, i) = cmplx_1
      end do

      if (invert .eq. 1) then
        call ZGESV(nxx, nxx, eh_00, nxx, ipiv, g, nxx, info)
        if (info .ne. 0) then
          write (stdout, *) 'ERROR:  IN ZGESV IN tran_green, INFO=', info
          call io_error('tran_green: problem in ZGESV 2')
        end if
      end if

    case (0)

      !  construct the bulk green's function gnn or (if surface=.true.) the
      !  sub-surface green's function

      s1 = cmplx_0
      s2 = cmplx_0
      ! s1 = H_01 * T
      call ZGEMM('N', 'N', nxx, nxx, nxx, cmplx_1, c1, nxx, tot, nxx, cmplx_0, s1, nxx)
      ! s2 = H_01^+ * T^bar
      call ZGEMM('C', 'N', nxx, nxx, nxx, cmplx_1, c1, nxx, tott, nxx, cmplx_0, s2, nxx)

      eh_00(:, :) = cmplx(-h_00(:, :), kind=dp) - s1(:, :) - s2(:, :)
      do i = 1, nxx
        eh_00(i, i) = cmplx(e_scan, kind=dp) + eh_00(i, i)
      end do

      g_inv(:, :) = eh_00(:, :)

      ! identity
      g = cmplx_0
      do i = 1, nxx
        g(i, i) = cmplx_1
      end do

      if (invert .eq. 1) then
        call ZGESV(nxx, nxx, eh_00, nxx, ipiv, g, nxx, info)
        if (info .ne. 0) then
          write (stdout, *) 'ERROR:  IN ZGESV IN tran_green, INFO=', info
          call io_error('tran_green: problem in ZGESV 3')
        end if
      end if

    end select

    deallocate (s2)
    if (ierr /= 0) call io_error('Error in deallocating s2 in tran_green')
    deallocate (s1)
    if (ierr /= 0) call io_error('Error in deallocating s1 in tran_green')
    deallocate (c1)
    if (ierr /= 0) call io_error('Error in deallocating c1 in tran_green')
    deallocate (eh_00)
    if (ierr /= 0) call io_error('Error in deallocating eh_00 in tran_green')
    deallocate (g_inv)
    if (ierr /= 0) call io_error('Error in deallocating g_inv in tran_green')
    deallocate (ipiv)
    if (ierr /= 0) call io_error('Error in deallocating ipiv in tran_green')

    return

  end subroutine tran_green

  !============================================!
  subroutine tran_read_htX(nxx, h_00, h_01, h_file)
    !============================================!

    use w90_constants, only: dp
    use w90_io, only: stdout, io_file_unit, io_error, maxlen

    implicit none

    integer, intent(in) ::  nxx
    real(kind=dp), intent(out) :: h_00(nxx, nxx), h_01(nxx, nxx)
    character(len=50), intent(in) :: h_file
    !
    integer :: i, j, nw, file_unit
    character(len=maxlen) :: dummy

    file_unit = io_file_unit()

    open (unit=file_unit, file=h_file, form='formatted', &
          status='old', action='read', err=101)

    write (stdout, '(/a)', advance='no') ' Reading H matrix from '//h_file//'  : '

    read (file_unit, '(a)', err=102, end=102) dummy
    write (stdout, '(a)') trim(dummy)

    read (file_unit, *, err=102, end=102) nw
    if (nw .ne. nxx) call io_error('wrong matrix size in transport: read_htX')
    read (file_unit, *) ((h_00(i, j), i=1, nxx), j=1, nxx)
    read (file_unit, *, err=102, end=102) nw
    if (nw .ne. nxx) call io_error('wrong matrix size in transport: read_htX')
    read (file_unit, *, err=102, end=102) ((h_01(i, j), i=1, nxx), j=1, nxx)

    close (unit=file_unit)

    return

101 call io_error('Error: Problem opening input file '//h_file)
102 call io_error('Error: Problem reading input file '//h_file)

  end subroutine tran_read_htX

  !============================================!
  subroutine tran_read_htC(nxx, h_00, h_file)
    !============================================!

    use w90_constants, only: dp
    use w90_io, only: stdout, io_file_unit, io_error, maxlen

    implicit none

    integer, intent(in) ::  nxx
    real(kind=dp), intent(out) :: h_00(nxx, nxx)
    character(len=50), intent(in) :: h_file
    !
    integer :: i, j, nw, file_unit
    character(len=maxlen) :: dummy

    file_unit = io_file_unit()

    open (unit=file_unit, file=h_file, form='formatted', &
          status='old', action='read', err=101)

    write (stdout, '(/a)', advance='no') ' Reading H matrix from '//h_file//'  : '

    read (file_unit, '(a)', err=102, end=102) dummy
    write (stdout, '(a)') trim(dummy)

    read (file_unit, *, err=102, end=102) nw
    if (nw .ne. nxx) call io_error('wrong matrix size in transport: read_htC')
    read (file_unit, *, err=102, end=102) ((h_00(i, j), i=1, nxx), j=1, nxx)

    close (unit=file_unit)

    return

101 call io_error('Error: Problem opening input file '//h_file)
102 call io_error('Error: Problem reading input file '//h_file)

  end subroutine tran_read_htC

  !============================================!
  subroutine tran_read_htXY(nxx1, nxx2, h_01, h_file)
    !============================================!

    use w90_constants, only: dp
    use w90_io, only: stdout, io_file_unit, io_error, maxlen

    implicit none

    integer, intent(in) ::  nxx1, nxx2
    real(kind=dp), intent(out) :: h_01(nxx1, nxx2)
    character(len=50), intent(in) :: h_file
    !
    integer :: i, j, nw1, nw2, file_unit
    character(len=maxlen) :: dummy

    file_unit = io_file_unit()

    open (unit=file_unit, file=h_file, form='formatted', &
          status='old', action='read', err=101)

    write (stdout, '(/a)', advance='no') ' Reading H matrix from '//h_file//'  : '

    read (file_unit, '(a)', err=102, end=102) dummy
    write (stdout, '(a)') trim(dummy)

    read (file_unit, *, err=102, end=102) nw1, nw2

    if (nw1 .ne. nxx1 .or. nw2 .ne. nxx2) call io_error('wrong matrix size in transport: read_htXY')

    read (file_unit, *, err=102, end=102) ((h_01(i, j), i=1, nxx1), j=1, nxx2)

    close (unit=file_unit)

    return

101 call io_error('Error: Problem opening input file '//h_file)
102 call io_error('Error: Problem reading input file '//h_file)

  end subroutine tran_read_htXY

!========================================
  subroutine tran_find_integral_signatures(signatures, num_G)
    !=========================================================================!
    ! Reads <seedname>.unkg file that contains the u_nk(G) and calculate      !
    ! Fourier components of each wannier function. Linear combinations of     !
    ! these provide integral of different spatial dependence.                 !
    ! The set of these integrals provide a signature for distinguishing the   !
    ! type and 'parity' of each wannier function.                             !
    !=========================================================================!
    use w90_constants, only: dp, cmplx_0, twopi, cmplx_i
    use w90_io, only: io_error, stdout, seedname, io_file_unit, io_date, &
      io_stopwatch

    use w90_parameters, only: num_wann, have_disentangled, num_bands, u_matrix, u_matrix_opt, &
      real_lattice, iprint, timing_level

    use w90_hamiltonian, only: wannier_centres_translated

    implicit none
    integer, intent(out)                                    :: num_G
    real(kind=dp), allocatable, dimension(:, :), intent(out)   :: signatures

    complex(kind=dp), allocatable                           :: unkg(:, :), tran_u_matrix(:, :)
    complex(kind=dp)                                       :: phase_factor, signature_basis(32)

    real(kind=dp)                                          :: i_unkg, r_unkg, wf_frac(3), det_rl, inv_t_rl(3, 3), &
                                                              mag_signature_sq

!~     character(len=11)                                      :: unkg_file

    logical                                                :: have_file

    integer, allocatable, dimension(:, :)                     :: g_abc
    integer                                                :: i, ibnd, file_unit, ierr, p, p_max, n, m, ig, a, b, c, ig_idx(32)

    if (timing_level > 1) call io_stopwatch('tran: find_sigs_unkg_int', 1)
    !
    file_unit = io_file_unit()
    inquire (file=trim(seedname)//'.unkg', exist=have_file)
    if (.not. have_file) call io_error('tran_hr_parity_unkg: file '//trim(seedname)//'.unkg not found')
    open (file_unit, file=trim(seedname)//'.unkg', form='formatted', action='read')
    !
    !Read unkg file
    !
    write (stdout, '(3a)') ' Reading '//trim(seedname)//'.unkg  file'
    read (file_unit, *) num_G

    allocate (signatures(20, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating signatures in tran_find_sigs_unkg_int')
    allocate (unkg(num_G, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating unkg in tran_find_sigs_unkg_int')
    allocate (g_abc(num_G, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating g_abc in tran_find_sigs_unkg_int')
    if (have_disentangled) then
      allocate (tran_u_matrix(num_bands, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating tran_u_matrix in tran_find_sigs_unkg_int')
    else
      allocate (tran_u_matrix(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating tran_u_matrix in tran_find_sigs_unkg_int')
    endif
    !
    unkg = cmplx_0
    do m = 1, num_bands
      do i = 1, num_G
        read (file_unit, *) ibnd, ig, a, b, c, r_unkg, i_unkg
        if ((ig .ne. i) .OR. (ibnd .ne. m)) then
          call io_error('tran_find_sigs_unkg_int: Incorrect bands or g vectors')
        endif
        unkg(i, m) = cmplx(r_unkg, i_unkg, kind=dp)
        g_abc(i, :) = (/a, b, c/)
      enddo
    enddo
    !
    close (file_unit)
    !
    ! Computing inverse of transpose of real_lattice
    !
    det_rl = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - real_lattice(2, 3)*real_lattice(3, 2)) &
             - real_lattice(2, 2)*(real_lattice(2, 1)*real_lattice(3, 3) - real_lattice(2, 3)*real_lattice(3, 1)) &
             + real_lattice(3, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - real_lattice(2, 2)*real_lattice(3, 1))

    inv_t_rl(1, 1) = (real_lattice(2, 2)*real_lattice(3, 3) - real_lattice(3, 2)*real_lattice(2, 3))
    inv_t_rl(1, 2) = (real_lattice(2, 1)*real_lattice(3, 3) - real_lattice(3, 1)*real_lattice(2, 3))
    inv_t_rl(1, 3) = (real_lattice(2, 1)*real_lattice(3, 2) - real_lattice(3, 1)*real_lattice(2, 2))

    inv_t_rl(2, 1) = (real_lattice(1, 2)*real_lattice(3, 3) - real_lattice(3, 2)*real_lattice(1, 3))
    inv_t_rl(2, 2) = (real_lattice(1, 1)*real_lattice(3, 3) - real_lattice(3, 1)*real_lattice(1, 3))
    inv_t_rl(2, 3) = (real_lattice(1, 1)*real_lattice(3, 2) - real_lattice(3, 1)*real_lattice(1, 2))

    inv_t_rl(3, 1) = (real_lattice(1, 2)*real_lattice(2, 3) - real_lattice(2, 2)*real_lattice(1, 3))
    inv_t_rl(3, 2) = (real_lattice(1, 1)*real_lattice(2, 3) - real_lattice(2, 1)*real_lattice(1, 3))
    inv_t_rl(3, 3) = (real_lattice(1, 1)*real_lattice(2, 2) - real_lattice(2, 1)*real_lattice(1, 2))

    inv_t_rl = inv_t_rl/det_rl
    !
    !Loop over wannier functions to compute generalised U matrix
    !
    signatures = 0.0_dp
    tran_u_matrix = cmplx_0
    do n = 1, num_wann
      !
      !Disentanglement step
      !
      if (have_disentangled) then
        do p = 1, num_bands
          do m = 1, num_wann
            tran_u_matrix(p, n) = tran_u_matrix(p, n) + u_matrix(m, n, 1)*u_matrix_opt(p, m, 1)
          enddo
        enddo
        p_max = num_bands
      else
        tran_u_matrix = u_matrix(:, :, 1)
        p_max = num_wann
      endif
    enddo

    if (iprint .ge. 5) write (stdout, *) 'Printing integral signatures for each wannier function:'
    !
    ! Loop over all wannier functions
    !
    do n = 1, num_wann
      !
      ! Find fraction coordinate of wannier function in lattice vector basis
      ! wf_frac(:)=(transpose(real_lattice))^(-1) * wannier_centres_translated(:,n)
      !
      wf_frac = 0.0_dp
      wf_frac = matmul(inv_t_rl, wannier_centres_translated(:, n))
      !
      ! Loop over all g vectors, find a,b,c's required
      !
      do ig = 1, num_G
        ! 0th Order
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(1) = ig   ! 1
        ! 1st Order
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(2) = ig   ! x
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(3) = ig   ! y
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 1)) ig_idx(4) = ig   ! z
        ! 2nd Order
        if ((g_abc(ig, 1) .eq. 2) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(5) = ig   ! x^2
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(6) = ig   ! xy
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. -1) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(7) = ig   ! xy
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 1)) ig_idx(8) = ig   ! xz
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. -1)) ig_idx(9) = ig   ! xz
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 2) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(10) = ig  ! y^2
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. 1)) ig_idx(11) = ig  ! yz
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. -1)) ig_idx(12) = ig  ! yz
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 2)) ig_idx(13) = ig  ! z^2
        ! 3rd Order
        if ((g_abc(ig, 1) .eq. 3) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(14) = ig  ! x^3
        if ((g_abc(ig, 1) .eq. 2) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(15) = ig  ! x^2y
        if ((g_abc(ig, 1) .eq. 2) .and. (g_abc(ig, 2) .eq. -1) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(16) = ig  ! x^2y
        if ((g_abc(ig, 1) .eq. 2) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 1)) ig_idx(17) = ig  ! x^2z
        if ((g_abc(ig, 1) .eq. 2) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. -1)) ig_idx(18) = ig  ! x^2z
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 2) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(19) = ig  ! xy^2
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. -2) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(20) = ig  ! xy^2
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. 1)) ig_idx(21) = ig  ! xyz
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. -1)) ig_idx(22) = ig  ! xyz
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. -1) .and. (g_abc(ig, 3) .eq. 1)) ig_idx(23) = ig  ! xyz
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. -1) .and. (g_abc(ig, 3) .eq. -1)) ig_idx(24) = ig  ! xyz
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 2)) ig_idx(25) = ig  ! xz^2
        if ((g_abc(ig, 1) .eq. 1) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. -2)) ig_idx(26) = ig  ! xz^2
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 3) .and. (g_abc(ig, 3) .eq. 0)) ig_idx(27) = ig  ! y^3
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 2) .and. (g_abc(ig, 3) .eq. 1)) ig_idx(28) = ig  ! y^2z
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 2) .and. (g_abc(ig, 3) .eq. -1)) ig_idx(29) = ig  ! y^2z
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. 2)) ig_idx(30) = ig  ! yz^2
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 1) .and. (g_abc(ig, 3) .eq. -2)) ig_idx(31) = ig  ! yz^2
        if ((g_abc(ig, 1) .eq. 0) .and. (g_abc(ig, 2) .eq. 0) .and. (g_abc(ig, 3) .eq. 3)) ig_idx(32) = ig  ! z^3
      enddo
      !
      ! Loop over the 32 required g-vectors
      !
      signature_basis = cmplx_0
      do ig = 1, 32
        phase_factor = cmplx_0
        !
        ! Compute the phase factor exp(-i*G*x_c)
        !
        phase_factor = exp(-twopi*cmplx_i*(g_abc(ig_idx(ig), 1)*wf_frac(1) &
                                           + g_abc(ig_idx(ig), 2)*wf_frac(2) &
                                           + g_abc(ig_idx(ig), 3)*wf_frac(3)))
        !
        ! Compute integrals that form the basis of the spatial integrals that form the signature
        do p = 1, p_max
          signature_basis(ig) = signature_basis(ig) + tran_u_matrix(p, n)*conjg(unkg(ig_idx(ig), p))
        enddo
        signature_basis(ig) = signature_basis(ig)*phase_factor
      enddo
      !
      ! Definitions of the signature integrals
      !
      ! 0th Order
      signatures(1, n) = real(signature_basis(1))                                                                  ! 1
      ! 1st Order
      signatures(2, n) = aimag(signature_basis(2))                                                                 ! x
      signatures(3, n) = aimag(signature_basis(3))                                                                 ! y
      signatures(4, n) = aimag(signature_basis(4))                                                                 ! z
      ! 2nd Orde r
      signatures(5, n) = real(signature_basis(1) - signature_basis(5))/2                                           ! x^2
      signatures(6, n) = real(signature_basis(7) - signature_basis(6))/2                                           ! xy
      signatures(7, n) = real(signature_basis(9) - signature_basis(8))/2                                           ! xz
      signatures(8, n) = real(signature_basis(1) - signature_basis(10))/2                                          ! y^2
      signatures(9, n) = real(signature_basis(12) - signature_basis(11))/2                                          ! yz
      signatures(10, n) = real(signature_basis(1) - signature_basis(13))/2                                          ! z^2
      ! 3rd Order
      signatures(11, n) = aimag(3*signature_basis(2) - signature_basis(14))/4                                        ! x^3
      signatures(12, n) = aimag(2*signature_basis(3) + signature_basis(16) - signature_basis(15))/4                    ! x^2y
      signatures(13, n) = aimag(2*signature_basis(4) + signature_basis(18) - signature_basis(17))/4                    ! x^2z
      signatures(14, n) = aimag(2*signature_basis(2) - signature_basis(20) - signature_basis(19))/4                    ! xy^2
      signatures(15, n) = aimag(signature_basis(23) + signature_basis(22) - signature_basis(21) - signature_basis(24))/4 ! xyz
      signatures(16, n) = aimag(2*signature_basis(2) - signature_basis(26) - signature_basis(25))/4                    ! xz^2
      signatures(17, n) = aimag(3*signature_basis(3) - signature_basis(27))/4                                        ! y^3
      signatures(18, n) = aimag(2*signature_basis(4) + signature_basis(29) - signature_basis(28))/4                    ! y^2z
      signatures(19, n) = aimag(2*signature_basis(3) - signature_basis(31) - signature_basis(30))/4                    ! yz^2
      signatures(20, n) = aimag(3*signature_basis(4) - signature_basis(32))/4                                        ! z^3

      if (iprint .ge. 5) then
        write (stdout, *) ' '
        write (stdout, *) ' Wannier function: ', n
        do ig = 1, 20
          write (stdout, *) ig - 1, signatures(ig, n)
        enddo
      endif
      !
      !Normalise signature of each wannier function to a unit vector
      !
      mag_signature_sq = 0.0_dp
      do ig = 1, 20
        mag_signature_sq = mag_signature_sq + abs(signatures(ig, n))**2
      enddo
      signatures(:, n) = signatures(:, n)/sqrt(mag_signature_sq)
      !
    enddo ! Wannier Function loop

    !
    ! Set num_G = 20 to ensure later subroutines work correctly
    !
    num_G = 20

    deallocate (tran_u_matrix, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating tran_u_matrix in tran_find_signatures')
    deallocate (g_abc, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating g_abc in tran_find_signatures')
    deallocate (unkg, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating unkg in tran_find_signatures')

    if (timing_level > 1) call io_stopwatch('tran: find_sigs_unkg_int', 2)

    return

  end subroutine tran_find_integral_signatures

  !========================================!
  subroutine tran_lcr_2c2_sort(signatures, num_G, pl_warning)
    !=======================================================!
    ! This is the main subroutine controling the sorting    !
    ! for the 2c2 geometry. We first sort in the conduction !
    ! direction, group, sort in 2nd direction, group and    !
    ! sort in 3rd direction. Rigourous checks are performed !
    ! to ensure group and subgroup structure is consistent  !
    ! between principal layers (PLs), and unit cells. Once  !
    ! checks are passed we consider the possibility of      !
    ! multiple wannier functions are of similar centre, and !
    ! sort those                                            !
    !=======================================================!

    use w90_constants, only: dp
    use w90_io, only: io_error, stdout, io_stopwatch
    use w90_parameters, only: one_dim_dir, tran_num_ll, num_wann, tran_num_cell_ll, &
      real_lattice, tran_group_threshold, iprint, timing_level, lenconfac, &
      wannier_spreads, write_xyz, dist_cutoff
    use w90_hamiltonian, only: wannier_centres_translated

    implicit none

    integer, intent(in)                                :: num_G
    real(dp), intent(in), dimension(:, :)                :: signatures
    logical, intent(out)                               :: pl_warning

    real(dp), dimension(2, num_wann)                    :: centres_non_sorted, centres_initial_sorted
    real(dp), dimension(2, tran_num_ll)                 :: PL1, PL2, PL3, PL4, PL
    real(dp), dimension(2, num_wann - (4*tran_num_ll))    :: central_region
    real(dp)                                          :: reference_position, &
                                                         cell_length, distance, PL_max_val, PL_min_val

!~    integer                                           :: l,max_i,iterator !aam: unused variables
    integer                                           :: i, j, k, PL_selector, &
                                                         sort_iterator, sort_iterator2, ierr, temp_coord_2, temp_coord_3, n, &
                                                         num_wann_cell_ll, num_wf_group1, num_wf_last_group
    integer, allocatable, dimension(:)                  :: PL_groups, &
                                                           PL1_groups, PL2_groups, PL3_groups, PL4_groups, central_region_groups
    integer, allocatable, dimension(:, :)                :: PL_subgroup_info, &
                                                            PL1_subgroup_info, PL2_subgroup_info, PL3_subgroup_info, &
                                                            PL4_subgroup_info, central_subgroup_info, temp_subgroup

    character(30)                                     :: fmt_1

    allocate (tran_sorted_idx(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tran_sorted_idx in tran_lcr_2c2_sort')

    num_wann_cell_ll = tran_num_ll/tran_num_cell_ll

    if (timing_level > 1) call io_stopwatch('tran: lcr_2c2_sort', 1)

    sort_iterator = 0
    !
    !Check translated centres have been found
    !
    if (size(wannier_centres_translated) .eq. 0) then
      call io_error('Translated centres not known : required perform lcr transport, try restart=plot')
    endif

    !read one_dim_dir and creates an array (coord) that correspond to the
    !conduction direction (coord(1)) and the two perpendicular directions
    !(coord(2),coord(3)), such that a right-handed set is formed
    !
    if (one_dim_dir .eq. 1) then
      coord(1) = 1
      coord(2) = 2
      coord(3) = 3
    elseif (one_dim_dir .eq. 2) then
      coord(1) = 2
      coord(2) = 3
      coord(3) = 1
    elseif (one_dim_dir .eq. 3) then
      coord(1) = 3
      coord(2) = 1
      coord(3) = 2
    endif
    !
    !Check
    !
    if (((real_lattice(coord(1), coord(2)) .ne. 0) .or. (real_lattice(coord(1), coord(3)) .ne. 0)) .or. &
        ((real_lattice(coord(2), coord(1)) .ne. 0) .or. (real_lattice(coord(3), coord(1)) .ne. 0))) then
      call io_error( &
      'Lattice vector in conduction direction must point along x,y or z &
      & direction and be orthogonal to the remaining lattice vectors.')
    endif
    !
    !Check
    !
    if (num_wann .le. 4*tran_num_ll) then
      call io_error('Principle layers are too big.')
    endif

100 continue
    !
    !Extract a 2d array of the wannier_indices and their coord(1) from wannier_centers_translated
    !
    do i = 1, num_wann
      centres_non_sorted(1, i) = i
      centres_non_sorted(2, i) = wannier_centres_translated(coord(1), i)
    enddo
    write (stdout, '(/a)') ' Sorting WFs into principal layers'
    !
    !Initial sorting according to coord(1).
    !
    call sort(centres_non_sorted, centres_initial_sorted)
    !
    !Extract principal layers. WARNING: This extraction implies the structure of the supercell is
    !2 principal layers of lead on the left and on the right of a central conductor.
    !
    PL1 = centres_initial_sorted(:, 1:tran_num_ll)
    PL2 = centres_initial_sorted(:, tran_num_ll + 1:2*tran_num_ll)
    PL3 = centres_initial_sorted(:, num_wann - (2*tran_num_ll - 1):num_wann - (tran_num_ll))
    PL4 = centres_initial_sorted(:, num_wann - (tran_num_ll - 1):)
    !
    if (sort_iterator .eq. 1) then
      temp_coord_2 = coord(2)
      temp_coord_3 = coord(3)
      coord(2) = temp_coord_3
      coord(3) = temp_coord_2
    endif
    !
    if (iprint .ge. 4) then
      write (stdout, *) ' Group Breakdown of each principal layer'
    endif
    !
    !Loop over principal layers
    !
    do i = 1, 4
      !
      !Creating a variable PL_selector which choose the appropriate PL
      !
      PL_selector = i
      select case (PL_selector)
      case (1)
        PL = PL1
      case (2)
        PL = PL2
      case (3)
        PL = PL3
      case (4)
        PL = PL4
      endselect
      !
      !Grouping wannier functions with similar coord(1)
      !
      call group(PL, PL_groups)

      if (iprint .ge. 4) then
        !
        !Print group breakdown
        !
        write (fmt_1, '(i5)') size(PL_groups)
        fmt_1 = adjustl(fmt_1)
        fmt_1 = '(a3,i1,a1,i5,a2,'//trim(fmt_1)//'i4,a1)'
        write (stdout, fmt_1) ' PL', i, ' ', size(PL_groups), ' (', (PL_groups(j), j=1, size(PL_groups)), ')'
      endif
      !
      !Returns the sorted PL and informations on this PL
      !
      allocate (PL_subgroup_info(size(PL_groups), maxval(PL_groups)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating PL_subgroup_info in tran_lcr_2c2_sort')
      call master_sort_and_group(PL, PL_groups, tran_num_ll, PL_subgroup_info)

      select case (PL_selector)
      case (1)
        allocate (PL1_groups(size(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL1_groups in tran_lcr_2c2_sort')
        allocate (PL1_subgroup_info(size(PL_groups), maxval(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL1_subgroup_info in tran_lcr_2c2_sort')
        !
        PL1 = PL
        PL1_groups = PL_groups
        PL1_subgroup_info = PL_subgroup_info
        !
        deallocate (PL_subgroup_info, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL1_subgroup_info in tran_lcr_2c2_sort')
      case (2)
        allocate (PL2_groups(size(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL2_groups in tran_lcr_2c2_sort')
        allocate (PL2_subgroup_info(size(PL_groups), maxval(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL2_subgroup_info in tran_lcr_2c2_sort')
        !
        PL2 = PL
        PL2_groups = PL_groups
        PL2_subgroup_info = PL_subgroup_info
        !
        deallocate (PL_subgroup_info, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL2_subgroup_info in tran_lcr_2c2_sort')
      case (3)
        allocate (PL3_groups(size(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL3_groups in tran_lcr_2c2_sort')
        allocate (PL3_subgroup_info(size(PL_groups), maxval(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL3_subgroup_info in tran_lcr_2c2_sort')
        !
        PL3 = PL
        PL3_groups = PL_groups
        PL3_subgroup_info = PL_subgroup_info
        !
        deallocate (PL_subgroup_info, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL3_subgroup_info in tran_lcr_2c2_sort')
      case (4)
        allocate (PL4_groups(size(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL4_groups in tran_lcr_2c2_sort')
        allocate (PL4_subgroup_info(size(PL_groups), maxval(PL_groups)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating PL4_subgroup_info in tran_lcr_2c2_sort')
        !
        PL4 = PL
        PL4_groups = PL_groups
        PL4_subgroup_info = PL_subgroup_info
        !
        deallocate (PL_subgroup_info, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL4_subgroup_info in tran_lcr_2c2_sort')
      endselect

      deallocate (PL_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL_groups in tran_lcr_2c2_sort')
    enddo ! Principal layer loop

    !
    !Grouping and sorting of central conductor region
    !
    !Define central region
    !
    central_region = centres_initial_sorted(:, 2*tran_num_ll + 1:num_wann - (2*tran_num_ll))
    !
    !Group central region
    !
    call group(central_region, central_region_groups)
    !
    !Print central region group breakdown
    !
    if (iprint .ge. 4) then
      write (stdout, *) ' Group Breakdown of central region'
      write (fmt_1, '(i5)') size(central_region_groups)
      fmt_1 = adjustl(fmt_1)
      fmt_1 = '(a5,i5,a2,'//trim(fmt_1)//'i4,a1)'
      write (stdout, fmt_1) '     ', size(central_region_groups), ' (', &
        (central_region_groups(j), j=1, size(central_region_groups)), ')'
    endif
    !
    !Returns sorted central group region
    !
    allocate (central_subgroup_info(size(central_region_groups), maxval(central_region_groups)), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating central_group_info in tran_lcr_2c2_sort')
    call master_sort_and_group(central_region, central_region_groups, num_wann - (4*tran_num_ll), central_subgroup_info)
    deallocate (central_subgroup_info, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating central_group_info in tran_lcr_2c2_sort')
    write (stdout, *) ' '
    !
    !Build the sorted index array
    !
    tran_sorted_idx = nint(centres_initial_sorted(1, :))
    tran_sorted_idx(1:tran_num_ll) = nint(PL1(1, :))
    tran_sorted_idx(tran_num_ll + 1:2*tran_num_ll) = nint(PL2(1, :))
    tran_sorted_idx(2*tran_num_ll + 1:num_wann - (2*tran_num_ll)) = nint(central_region(1, :))
    tran_sorted_idx(num_wann - (2*tran_num_ll - 1):num_wann - (tran_num_ll)) = nint(PL3(1, :))
    tran_sorted_idx(num_wann - (tran_num_ll - 1):) = nint(PL4(1, :))

    sort_iterator = sort_iterator + 1
    !
    !Checks:
    !
    if ((size(PL1_groups) .ne. size(PL2_groups)) .or. &
        (size(PL2_groups) .ne. size(PL3_groups)) .or. &
        (size(PL3_groups) .ne. size(PL4_groups))) then
      if (sort_iterator .ge. 2) then
        if (write_xyz) call tran_write_xyz()
        call io_error('Sorting techniques exhausted:&
          & Inconsistent number of groups among principal layers')
      endif
      write (stdout, *) 'Inconsistent number of groups among principal layers: restarting sorting...'
      deallocate (PL1_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL1_groups in tran_lcr_2c2_sort')
      deallocate (PL2_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL2_groups in tran_lcr_2c2_sort')
      deallocate (PL3_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL3_groups in tran_lcr_2c2_sort')
      deallocate (PL4_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL4_groups in tran_lcr_2c2_sort')
      deallocate (PL1_subgroup_info)
      if (ierr /= 0) call io_error('Error deallocating PL1_subgroup_info in tran_lcr_2c2_sort')
      deallocate (PL2_subgroup_info, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL2_subgroup_info in tran_lcr_2c2_sort')
      deallocate (PL3_subgroup_info, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL3_subgroup_info in tran_lcr_2c2_sort')
      deallocate (PL4_subgroup_info, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL4_subgroup_info in tran_lcr_2c2_sort')
      goto 100
    endif
    !
    do i = 1, size(PL1_groups)
      if ((PL1_groups(i) .ne. PL2_groups(i)) .or. &
          (PL2_groups(i) .ne. PL3_groups(i)) .or. &
          (PL3_groups(i) .ne. PL4_groups(i))) then
        if (sort_iterator .ge. 2) then
          if (write_xyz) call tran_write_xyz()
          call io_error &
           ('Sorting techniques exhausted: Inconsitent number of wannier function among &
             & similar groups within principal layers')
        endif
        write (stdout, *) 'Inconsitent number of wannier function among &
          &similar groups within& principal layers: restarting sorting...'

        deallocate (PL1_groups, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL1_groups in tran_lcr_2c2_sort')
        deallocate (PL2_groups, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL2_groups in tran_lcr_2c2_sort')
        deallocate (PL3_groups, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL3_groups in tran_lcr_2c2_sort')
        deallocate (PL4_groups, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL4_groups in tran_lcr_2c2_sort')
        deallocate (PL1_subgroup_info)
        if (ierr /= 0) call io_error('Error deallocating PL1_subgroup_info in tran_lcr_2c2_sort')
        deallocate (PL2_subgroup_info, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL2_subgroup_info in tran_lcr_2c2_sort')
        deallocate (PL3_subgroup_info, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL3_subgroup_info in tran_lcr_2c2_sort')
        deallocate (PL4_subgroup_info, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating PL4_subgroup_info in tran_lcr_2c2_sort')
        goto 100
      endif
    enddo
    !
    !Now we check that the leftmost group and the rightmost group aren't
    !supposed to be the same group
    !
    reference_position = wannier_centres_translated(coord(1), tran_sorted_idx(1))
    cell_length = real_lattice(coord(1), coord(1))
    sort_iterator2 = 1
    do i = 1, tran_num_ll
      distance = abs(abs(reference_position - wannier_centres_translated(coord(1), tran_sorted_idx(num_wann - i + 1))) &
                     - cell_length)
      if (distance .lt. tran_group_threshold) then
        wannier_centres_translated(coord(1), tran_sorted_idx(num_wann - i + 1)) = &
          wannier_centres_translated(coord(1), tran_sorted_idx(num_wann - i + 1)) - cell_length
        sort_iterator2 = sort_iterator2 + 1
      endif
    enddo

    if (sort_iterator2 .gt. 1) then
      write (stdout, *) ' Grouping inconsistency found between first and last unit cells: '
      write (stdout, *) ' suspect Wannier functions have been translated. '
      write (stdout, *) ' '
      write (stdout, *) ' Rebuilding Hamiltonian...'
      write (stdout, *) ' '
      deallocate (hr_one_dim, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating hr_one_dim in tran_lcr_2c2_sort')
      call tran_reduce_hr()
      call tran_cut_hr_one_dim()
      write (stdout, *) ' '
      write (stdout, *) ' Restarting sorting...'
      write (stdout, *) ' '
      sort_iterator = sort_iterator - 1
      deallocate (PL1_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL1_groups in tran_lcr_2c2_sort')
      deallocate (PL2_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL2_groups in tran_lcr_2c2_sort')
      deallocate (PL3_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL3_groups in tran_lcr_2c2_sort')
      deallocate (PL4_groups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL4_groups in tran_lcr_2c2_sort')
      deallocate (PL1_subgroup_info)
      if (ierr /= 0) call io_error('Error deallocating PL1_subgroup_info in tran_lcr_2c2_sort')
      deallocate (PL2_subgroup_info, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL2_subgroup_info in tran_lcr_2c2_sort')
      deallocate (PL3_subgroup_info, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL3_subgroup_info in tran_lcr_2c2_sort')
      deallocate (PL4_subgroup_info, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating PL4_subgroup_info in tran_lcr_2c2_sort')
      goto 100
    endif
    !
    ! if we reach this point, we don't have any left/right problems anymore. So we now
    ! check for inconsistencies in subgroups
    !
    allocate (temp_subgroup(size(PL1_subgroup_info, 1), size(PL1_subgroup_info, 2)), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tmp_subgroup in tran_lcr_2c2_sort')
    do i = 2, 4
      select case (i)
      case (2)
        temp_subgroup = PL1_subgroup_info - PL2_subgroup_info
      case (3)
        temp_subgroup = PL1_subgroup_info - PL3_subgroup_info
      case (4)
        temp_subgroup = PL1_subgroup_info - PL4_subgroup_info
      end select
      do j = 1, size(temp_subgroup, 1)
        do k = 1, size(temp_subgroup, 2)
          if (temp_subgroup(j, k) .ne. 0) then
            if (sort_iterator .ge. 2) then
              if (write_xyz) call tran_write_xyz()
              call io_error &
                ('Sorting techniques exhausted: Inconsitent subgroup structures among principal layers')
            endif
            write (stdout, *) 'Inconsitent subgroup structure among principal layers: restarting sorting...'
            deallocate (temp_subgroup, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating tmp_subgroup in tran_lcr_2c2_sort')
            deallocate (PL1_groups, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating PL1_groups in tran_lcr_2c2_sort')
            deallocate (PL2_groups, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating PL2_groups in tran_lcr_2c2_sort')
            deallocate (PL3_groups, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating PL3_groups in tran_lcr_2c2_sort')
            deallocate (PL4_groups, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating PL4_groups in tran_lcr_2c2_sort')
            deallocate (PL1_subgroup_info)
            if (ierr /= 0) call io_error('Error deallocating PL1_subgroup_info in tran_lcr_2c2_sort')
            deallocate (PL2_subgroup_info, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating PL2_subgroup_info in tran_lcr_2c2_sort')
            deallocate (PL3_subgroup_info, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating PL3_subgroup_info in tran_lcr_2c2_sort')
            deallocate (PL4_subgroup_info, stat=ierr)
            if (ierr /= 0) call io_error('Error deallocating PL4_subgroup_info in tran_lcr_2c2_sort')
            goto 100
          endif
        enddo
      enddo
    enddo
    !
    ! At this point, every check has been cleared, and we need to use
    ! the parity signatures of the WFs for the possibility of equal centres
    !
    call check_and_sort_similar_centres(signatures, num_G)

    write (stdout, *) ' '
    write (stdout, *) '------------------------- Sorted Wannier Centres -----------------------------'
    write (stdout, *) 'Sorted index   Unsorted index           x           y           z     Spread  '
    write (stdout, *) '==================================PL1========================================='
    n = 0
    do k = 1, 4
      if (k .eq. 2) write (stdout, *) '==================================PL2========================================='
      if (k .eq. 3) then
        write (stdout, *) '===========================Central Region==================================='
        do i = 1, num_wann - 4*tran_num_ll
          n = n + 1
          write (stdout, FMT='(2x,i6,10x,i6,6x,4F12.6)') n, tran_sorted_idx(n), &
            wannier_centres_translated(1, tran_sorted_idx(n)), &
            wannier_centres_translated(2, tran_sorted_idx(n)), &
            wannier_centres_translated(3, tran_sorted_idx(n)), &
            wannier_spreads(tran_sorted_idx(n))*lenconfac**2
        enddo
        write (stdout, *) '==================================PL3========================================='
      endif
      if (k .eq. 4) write (stdout, *) '==================================PL4========================================='
      do i = 1, tran_num_cell_ll
        do j = 1, num_wann_cell_ll
          n = n + 1
          write (stdout, FMT='(2x,i6,10x,i6,6x,4F12.6)') n, tran_sorted_idx(n), &
            wannier_centres_translated(1, tran_sorted_idx(n)), &
            wannier_centres_translated(2, tran_sorted_idx(n)), &
            wannier_centres_translated(3, tran_sorted_idx(n)), &
            wannier_spreads(tran_sorted_idx(n))*lenconfac**2
        enddo
        if (i .ne. tran_num_cell_ll) write (stdout, *) '---------------------&
          &---------------------------------------------------------'

      enddo
    enddo

    write (stdout, *) '=============================================================================='
    write (stdout, *) ' '

    !
    ! MS: Use sorting to assess whether dist_cutoff is suitable for correct PL cut
    !     by using limits of coord(1) values in 1st and last groups of PL1, & 1st group of PL2
    !
    pl_warning = .false.
    num_wf_group1 = size(PL1_subgroup_info(1, :))
    if (size(PL1_groups) .ge. 1) then
      num_wf_last_group = size(PL1_subgroup_info(size(PL1_groups), :))
    else
      !
      !Exception for 1 group in unit cell.
      !
      num_wf_last_group = num_wann_cell_ll
    endif
    PL_min_val = maxval(wannier_centres_translated(coord(1), tran_sorted_idx(tran_num_ll - num_wf_last_group + 1:tran_num_ll))) &
                 - minval(wannier_centres_translated(coord(1), tran_sorted_idx(1:num_wf_group1)))
    PL_max_val = minval(wannier_centres_translated(coord(1), tran_sorted_idx(tran_num_ll + 1:tran_num_ll + num_wf_group1))) &
                 - minval(wannier_centres_translated(coord(1), tran_sorted_idx(1:num_wf_group1)))
    if ((dist_cutoff .lt. PL_min_val) .or. (dist_cutoff .gt. PL_max_val)) then
      write (stdout, '(a)') ' WARNING: Expected dist_cutoff to be a PL length, I think this'
      write (stdout, '(2(a,f10.6),a)') ' WARNING: is somewhere between ', PL_min_val, ' and ', PL_max_val, ' Ang'
      pl_warning = .true.
    endif
    !
    ! End MS.
    !
    if (timing_level > 1) call io_stopwatch('tran: lcr_2c2_sort', 2)

    return

  end subroutine tran_lcr_2c2_sort

  !========================================!
  subroutine master_sort_and_group(Array, Array_groups, Array_size, subgroup_info)
    !=============================================================!
    ! General sorting and grouping subroutine which takes Array,  !
    ! an ordered in conduction direction array of wannier function!
    ! indexes and positions, and returns the ordered (and grouped)!
    ! indexes and positions after considering the other two       !
    ! directions. Sub group info is also return for later checks. !
    !=============================================================!

    use w90_constants, only: dp
    use w90_io, only: io_error, stdout, io_stopwatch
    use w90_parameters, only: iprint, timing_level
    use w90_hamiltonian, only: wannier_centres_translated

    implicit none

    integer, intent(in), dimension(:)                 :: Array_groups
    integer, intent(in)                              :: Array_size

    integer, intent(out), allocatable, dimension(:, :)  :: subgroup_info

    real(dp), intent(inout), dimension(2, Array_size)  :: Array

    integer                                         :: i, j, k, Array_num_groups, increment, ierr, &
                                                       subgroup_increment, group_num_subgroups
    integer, allocatable, dimension(:)                :: group_subgroups

    real(dp), allocatable, dimension(:, :)             :: group_array, sorted_group_array, &
                                                          subgroup_array, sorted_subgroup_array
    character(30)                                   :: fmt_2

    if (timing_level > 2) call io_stopwatch('tran: lcr_2c2_sort: master_sort', 1)

    allocate (subgroup_info(size(Array_groups), maxval(Array_groups)), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating subgroup_info in master_sort_and_group')
    subgroup_info = 0
    !
    !Number of groups inside the principal layer
    !
    Array_num_groups = size(Array_groups)

    !
    !Convenient variable which will be amended later. Used to appropriately extract the group array from the Array
    !
    increment = 1
    !
    !Loop over groups inside Array
    !
    do j = 1, Array_num_groups
      allocate (group_array(2, Array_groups(j)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating group_array in master_sort_and_group')
      allocate (sorted_group_array(2, Array_groups(j)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating sorted_group_array in master_sort_and_group')
      !
      !Extract the group from the Array
      !
      group_array = Array(:, increment:increment + Array_groups(j) - 1)
      !
      !Updating group_array to contain coord(2)
      !
      do k = 1, Array_groups(j)
        group_array(2, k) = wannier_centres_translated(coord(2), int(group_array(1, k)))
      enddo

      call sort(group_array, sorted_group_array)
      call group(sorted_group_array, group_subgroups)

      group_num_subgroups = size(group_subgroups)

      if (iprint .ge. 4) then
        !
        !Printing subgroup breakdown
        !
        write (fmt_2, '(i5)') group_num_subgroups
        fmt_2 = adjustl(fmt_2)
        fmt_2 = '(a7,i3,a1,i5,a2,'//trim(fmt_2)//'i4,a1)'
        write (stdout, fmt_2) ' Group ', j, ' ', group_num_subgroups, ' (', (group_subgroups(i), i=1, group_num_subgroups), ')'
      endif
      !
      ! filling up subgroup_info
      !
      do k = 1, group_num_subgroups
        subgroup_info(j, k) = group_subgroups(k)
      enddo
      !
      !Convenient variable which will be amended later. Used to appropriately extract the subgroup array from the group_array
      !
      subgroup_increment = 1
      !
      !Loop over subgroups inside group
      !
      do k = 1, group_num_subgroups
        allocate (subgroup_array(2, group_subgroups(k)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating subgroup_array in master_sort_and_group')
        allocate (sorted_subgroup_array(2, group_subgroups(k)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating sorted_subgroup_array in master_sort_and_group')
        !
        !Extract the subgroup from the group
        !
        subgroup_array = sorted_group_array(:, subgroup_increment:subgroup_increment + group_subgroups(k) - 1)
        !
        !Updating subgroup_array to contain coord(3)
        !
        do i = 1, group_subgroups(k)
          subgroup_array(2, i) = wannier_centres_translated(coord(3), int(subgroup_array(1, i)))
        enddo

        call sort(subgroup_array, sorted_subgroup_array)
        !
        !Update sorted_group array with the sorted subgroup array
        !
        sorted_group_array(:, subgroup_increment:subgroup_increment + group_subgroups(k) - 1) = sorted_subgroup_array
        !
        !Update the subgroup_increment
        !
        subgroup_increment = subgroup_increment + group_subgroups(k)
        deallocate (sorted_subgroup_array, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating sorted_subgroup_array in master_sort_and_group')
        deallocate (subgroup_array, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating subgroup_array in master_sort_and_group')
      enddo
      !
      !Update Array with the sorted group array
      !
      Array(:, increment:increment + Array_groups(j) - 1) = sorted_group_array
      !
      !Update the group increment
      !
      increment = increment + Array_groups(j)
      deallocate (group_array, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating group_array in master_sort_and_group')
      deallocate (sorted_group_array, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating sorted_group_array in master_sort_and_group')
      deallocate (group_subgroups, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating group_subgroups in master_sort_and_group')
    enddo

    if (timing_level > 2) call io_stopwatch('tran: lcr_2c2_sort: master_sort', 2)

    return

  end subroutine master_sort_and_group

  !========================================!
  subroutine sort(non_sorted, sorted)
    !========================================!

    use w90_constants, only: dp

    implicit none

    real(dp), intent(inout), dimension(:, :)       :: non_sorted
    real(dp), intent(out), dimension(:, :)         :: sorted

    integer, dimension(1)                        :: min_loc
    integer                                     :: num_col, i

    num_col = size(non_sorted, 2)

    do i = 1, num_col
      !
      !look for the location of the minimum value of the coordinates in non_sorted
      !
      min_loc = minloc(non_sorted(2, :))
      !
      !now the index in the first row of sorted is the index non_sorted(1,min_loc)
      !
      sorted(1, i) = non_sorted(1, min_loc(1))
      !
      !here is the corresponding coordinate
      !
      sorted(2, i) = non_sorted(2, min_loc(1))
      !
      !here one replaces the minimum coordinate with 10**10 such that this value
      !will not be picked-up again by minloc
      !
      non_sorted(2, min_loc(1)) = 10.0**10
    enddo

    return

  endsubroutine sort

  !========================================!
  subroutine group(array, array_groups)
    !========================================!

    use w90_constants, only: dp
    use w90_io, only: io_error

    use w90_parameters, only: tran_group_threshold

    implicit none

    real(dp), intent(in), dimension(:, :)           :: array
    integer, intent(out), allocatable, dimension(:) :: array_groups

    integer, allocatable, dimension(:)             :: dummy_array
    logical, allocatable, dimension(:)             :: logic
    integer                                      :: array_idx, i, j, group_number, array_size, ierr

    array_size = size(array, 2)

    allocate (dummy_array(array_size), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating dummy_array in group')
    allocate (logic(array_size), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating logic in group')
    !
    !Initialise dummy array
    !
    dummy_array = 0
    !
    !Initialise logic to false
    !
    logic = .false.
    !
    !Define counter of number of groups
    !
    array_idx = 1
    !
    !Loop over columns of array (ie array_size)
    !
    do i = 1, array_size
      !
      !If an element of logic is true then it means the wannier function has already been grouped
      !
      if (logic(i) .eqv. .false.) then
        !
        !Create a group for the wannier function
        !
        logic(i) = .true.
        !
        !Initialise the number of wannier functions in this group to be 1
        !
        group_number = 1
        !
        !Loop over the rest of wannier functions in array
        !
        do j = min(i + 1, array_size), array_size
          !
          !Special termination cases
          !
          if ((j .eq. 1) .or. (i .eq. array_size)) then
            dummy_array(array_idx) = group_number
            exit
          endif
          if (j .eq. array_size .and. (abs(array(2, j) - array(2, i)) .le. tran_group_threshold)) then
            group_number = group_number + 1
            dummy_array(array_idx) = group_number
            logic(j) = .true.
            exit
          endif
          !
          !Check distance between wannier function_i and wannier function_j
          !
          if (abs(array(2, j) - array(2, i)) .le. tran_group_threshold) then
            !
            !Increment number of wannier functions in group
            !
            group_number = group_number + 1
            !
            !Assigns wannier function to the group
            !
            logic(j) = .true.
          else
            !
            !Group is finished and store number of wanniers in the group to dummy_array
            !
            dummy_array(array_idx) = group_number
            !
            !Increment number of groups
            !
            array_idx = array_idx + 1
            exit
          endif
        enddo
      endif
    enddo
    !
    !Copy elements of dummy_array to array_groups
    !
    allocate (array_groups(array_idx), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating array_groups in group')
    array_groups = dummy_array(:array_idx)

    deallocate (dummy_array, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating dummy_array in group')
    deallocate (logic, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating logic in group')

    return

  end subroutine group

  !=========================================================
  subroutine check_and_sort_similar_centres(signatures, num_G)
    !=======================================================!
    ! Here, we consider the possiblity of wannier functions !
    ! with similar centres, such as a set of d-orbitals     !
    ! centred on an atom. We use tran_group_threshold to    !
    ! identify them, if they exist, then use the signatures !
    ! to dishinguish and sort then consistently from unit   !
    ! cell to unit cell.                                    !
    !                                                       !
    ! MS: For 2two-shot and beyond, some parameters,        !
    ! eg, first_group_element will need changing to consider!
    ! the geometry of the new systems.                      !
    !=======================================================!

    use w90_constants, only: dp
    use w90_io, only: stdout, io_stopwatch, io_error
    use w90_parameters, only: tran_num_ll, num_wann, tran_num_cell_ll, iprint, timing_level, &
      tran_group_threshold, write_xyz
    use w90_hamiltonian, only: wannier_centres_translated

    implicit none

    integer, intent(in)                                :: num_G
    real(dp), intent(in), dimension(:, :)                :: signatures

    integer                                           :: i, j, k, l, ierr, group_iterator, coord_iterator, num_wf_iterator, &
                                                         num_wann_cell_ll, iterator, max_position(1), p, num_wf_cell_iter

    integer, allocatable, dimension(:)                  :: idx_similar_wf, group_verifier, sorted_idx, centre_id
    real(dp), allocatable, dimension(:)                 :: dot_p
    integer, allocatable, dimension(:, :)                :: tmp_wf_verifier, wf_verifier, first_group_element, &
                                                            ref_similar_centres, unsorted_similar_centres
    integer, allocatable, dimension(:, :, :)              :: wf_similar_centres

    logical, allocatable, dimension(:)                  :: has_similar_centres

    if (timing_level > 2) call io_stopwatch('tran: lcr_2c2_sort: similar_centres', 1)

    num_wann_cell_ll = tran_num_ll/tran_num_cell_ll

    allocate (wf_similar_centres(tran_num_cell_ll*4, num_wann_cell_ll, num_wann_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating wf_similar_centre in check_and_sort_similar_centres')
    allocate (idx_similar_wf(num_wann_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating idx_similar_wf in check_and_sort_similar_centres')
    allocate (has_similar_centres(num_wann_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating has_similar_centres in check_and_sort_similar_centres')
    allocate (tmp_wf_verifier(4*tran_num_cell_ll, num_wann_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tmp_wf_verifier in check_and_sort_similar_centres')
    allocate (group_verifier(4*tran_num_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating group_verifier in check_and_sort_similar_centres')
    allocate (first_group_element(4*tran_num_cell_ll, num_wann_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating first_group_element in check_and_sort_similar_centres')
    allocate (centre_id(num_wann_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating centre_id in check_and_sort_similar_centres')

    !
    ! First find WFs with similar centres: store in wf_similar_centres(cell#,group#,WF#)
    !
    group_verifier = 0
    tmp_wf_verifier = 0
    first_group_element = 0
    centre_id = 0
    !
    ! Loop over unit cells in PL1,PL2,PL3 and PL4
    !
    do i = 1, 4*tran_num_cell_ll
      group_iterator = 0
      has_similar_centres = .false.
      !
      ! Loops over wannier functions in present unit cell
      !
      num_wf_cell_iter = 0
      do j = 1, num_wann_cell_ll
        num_wf_iterator = 0
        !
        ! 2nd Loop over wannier functions in the present unit cell
        !
        do k = 1, num_wann_cell_ll
          if ((.not. has_similar_centres(k)) .and. (j .ne. k)) then
            coord_iterator = 0
            !
            ! Loop over x,y,z to find similar centres
            !
            do l = 1, 3
              if (i .le. 2*tran_num_cell_ll) then
                if (abs(wannier_centres_translated(coord(l), tran_sorted_idx(j + (i - 1)*num_wann_cell_ll)) - &
                        wannier_centres_translated(coord(l), tran_sorted_idx(k + (i - 1)*num_wann_cell_ll))) &
                    .le. tran_group_threshold) then
                  coord_iterator = coord_iterator + 1
                else
                  exit
                endif
              else
                if (abs(wannier_centres_translated(coord(l), tran_sorted_idx(num_wann - 2*tran_num_ll + &
                                                                             j + (i - 2*tran_num_cell_ll - 1)*num_wann_cell_ll)) - &
                        wannier_centres_translated(coord(l), tran_sorted_idx(num_wann - 2*tran_num_ll + &
                                                                             k + (i - 2*tran_num_cell_ll - 1)*num_wann_cell_ll))) &
                    .le. tran_group_threshold) then
                  coord_iterator = coord_iterator + 1
                else
                  exit
                endif
              endif
            enddo
            if (coord_iterator .eq. 3) then
              if (.not. has_similar_centres(j)) then
                num_wf_iterator = num_wf_iterator + 1
                if (i .le. 2*tran_num_cell_ll) then
                  idx_similar_wf(num_wf_iterator) = tran_sorted_idx(j + (i - 1)*num_wann_cell_ll)
                else
                  idx_similar_wf(num_wf_iterator) = tran_sorted_idx(j + num_wann - 2*tran_num_ll + &
                                                                    (i - 2*tran_num_cell_ll - 1)*num_wann_cell_ll)
                endif
                if (i .le. 2*tran_num_cell_ll) then
                  first_group_element(i, j) = j + (i - 1)*num_wann_cell_ll
                else
                  first_group_element(i, j) = num_wann - 2*tran_num_ll + &
                                              j + (i - 2*tran_num_cell_ll - 1)*num_wann_cell_ll
                endif
                num_wf_cell_iter = num_wf_cell_iter + 1
                centre_id(num_wf_cell_iter) = j
              endif
              has_similar_centres(k) = .true.
              has_similar_centres(j) = .true.
              num_wf_iterator = num_wf_iterator + 1
              if (i .le. 2*tran_num_cell_ll) then
                idx_similar_wf(num_wf_iterator) = tran_sorted_idx(k + (i - 1)*num_wann_cell_ll)
              else
                idx_similar_wf(num_wf_iterator) = tran_sorted_idx(k + num_wann - 2*tran_num_ll + &
                                                                  (i - 2*tran_num_cell_ll - 1)*num_wann_cell_ll)
              endif
            endif
          endif
        enddo ! loop over k
        if (num_wf_iterator .gt. 0) then
          group_iterator = group_iterator + 1
          wf_similar_centres(i, group_iterator, :) = idx_similar_wf(:)
          !
          !Save number of WFs in each group
          !
          tmp_wf_verifier(i, group_iterator) = num_wf_iterator
        endif
      enddo
      if ((count(has_similar_centres) .eq. 0) .and. (i .eq. 1)) then
        write (stdout, '(a)') ' No wannier functions found with similar centres: sorting completed'
        exit
      elseif (i .eq. 1) then
        write (stdout, *) ' Wannier functions found with similar centres: '
        write (stdout, *) '  -> using signatures to complete sorting '
      endif
      !
      !Save number of group of WFs in each unit cell and compare to previous unit cell
      !
      group_verifier(i) = group_iterator
      if (iprint .ge. 4) write (stdout, '(a11,i4,a13,i4)') ' Unit cell:', i, '  Num groups:', group_verifier(i)
      if (i .ne. 1) then
        if (group_verifier(i) .ne. group_verifier(i - 1)) then
          if (write_xyz) call tran_write_xyz()
          call io_error('Inconsitent number of groups of similar centred wannier functions between unit cells')
        elseif (i .eq. 4*tran_num_cell_ll) then
          write (stdout, *) ' Consistent groups of similar centred wannier functions between '
          write (stdout, *) ' unit cells found'
          write (stdout, *) ' '
        endif
      endif
    enddo  !Loop over all unit cells in PL1,PL2,PL3,PL4
    !
    ! Perform check to ensure consistent number of WFs between equivalent groups in different unit cells
    !
    if (any(has_similar_centres)) then
      !
      !
      allocate (wf_verifier(4*tran_num_cell_ll, group_verifier(1)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating wf_verifier in check_and_sort_similar_centres')
      !
      !
      if (iprint .ge. 4) write (stdout, *) 'Unit cell   Group number   Num WFs'
      wf_verifier = 0
      wf_verifier = tmp_wf_verifier(:, 1:group_verifier(1))
      do i = 1, 4*tran_num_cell_ll
        do j = 1, group_verifier(1)
          if (iprint .ge. 4) write (stdout, '(a3,i4,a9,i4,a7,i4)') '   ', i, '         ', j, '       ', wf_verifier(i, j)
          if (i .ne. 1) then
            if (wf_verifier(i, j) .ne. wf_verifier(i - 1, j)) &
                call io_error('Inconsitent number of wannier &
                  &functions between equivalent groups of similar &
                &centred wannier functions')
          endif
        enddo
      enddo
      write (stdout, *) ' Consistent number of wannier functions between equivalent groups of similar'
      write (stdout, *) ' centred wannier functions'
      write (stdout, *) ' '
      !
      write (stdout, *) ' Fixing order of similar centred wannier functions using parity signatures'
      !
      do i = 2, 4*tran_num_cell_ll
        do j = 1, group_verifier(1)
          !
          ! Make array of WF numbers which act as a reference to sort against
          ! and an array which need sorting
          !
          allocate (ref_similar_centres(group_verifier(1), wf_verifier(1, j)), stat=ierr)
          if (ierr /= 0) call io_error('Error in allocating ref_similar_centres in check_and_sort_similar_centres')
          allocate (unsorted_similar_centres(group_verifier(1), wf_verifier(1, j)), stat=ierr)
          if (ierr /= 0) call io_error('Error in allocating unsorted_similar_centres in check_and_sort_similar_centres')
          allocate (sorted_idx(wf_verifier(1, j)), stat=ierr)
          if (ierr /= 0) call io_error('Error in allocating sorted_idx in check_and_sort_similar_centres')
          allocate (dot_p(wf_verifier(1, j)), stat=ierr)
          if (ierr /= 0) call io_error('Error in allocating dot_p in check_and_sort_similar_centres')
          !
          do k = 1, wf_verifier(1, j)
            ref_similar_centres(j, k) = wf_similar_centres(1, j, k)
            unsorted_similar_centres(j, k) = wf_similar_centres(i, j, k)
          enddo
          !
          sorted_idx = 0
          do k = 1, wf_verifier(1, j)
            dot_p = 0.0_dp
            !
            ! building the array of positive dot products of signatures between unsorted_similar_centres(j,k)
            ! and all the ref_similar_centres(j,:)
            !
            do l = 1, wf_verifier(1, j)
              do p = 1, num_G
                dot_p(l) = dot_p(l) + abs(signatures(p, unsorted_similar_centres(j, k)))* &
                           abs(signatures(p, ref_similar_centres(j, l)))
              enddo
            enddo
            !
            max_position = maxloc(dot_p)
            !
            sorted_idx(max_position(1)) = unsorted_similar_centres(j, k)
          enddo
          !
          ! we have the properly ordered indexes for group j in unit cell i, now we need
          ! to overwrite the tran_sorted_idx array at the proper position
          !
          tran_sorted_idx(first_group_element(i, centre_id(j)):first_group_element(i, centre_id(j)) +&
               &wf_verifier(i, j) - 1) = sorted_idx(:)
          !
          deallocate (dot_p, stat=ierr)
          if (ierr /= 0) call io_error('Error in deallocating dot_p in check_and_sort_similar_centres')
          deallocate (sorted_idx, stat=ierr)
          if (ierr /= 0) call io_error('Error in deallocating sorted_idx in check_and_sort_similar_centres')
          deallocate (unsorted_similar_centres, stat=ierr)
          if (ierr /= 0) call io_error('Error in deallocating unsorted_similar_centres in check_and_sort_similar_centres')
          deallocate (ref_similar_centres, stat=ierr)
          if (ierr /= 0) call io_error('Error in deallocating ref_similar_centres in check_and_sort_similar_centres')
        enddo
      enddo
      !
      ! checking that all the indices of WFs in the new tran_sorted_idx are distinct
      ! Remark: physically, no two WFs with similar centres can have the same type so we should expect
      ! this check to always pass unless the signatures/wf are very weird !!
      !
      do k = 1, num_wann
        iterator = 0
        do l = 1, num_wann
          if (tran_sorted_idx(l) .eq. k) then
            iterator = iterator + 1
          endif
        enddo
        !
        if ((iterator .ge. 2) .or. (iterator .eq. 0)) call io_error( &
        'A Wannier Function appears either zero times or twice after sorting, this may be due to a &
        &poor wannierisation and/or disentanglement')
        !write(stdout,*) ' WF : ',k,' appears ',iterator,' time(s)'
      enddo
      deallocate (wf_verifier, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating wf_verifier in check_and_sort_similar_centres')
    endif

    deallocate (centre_id, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating centre_id in check_and_sort_similar_centres')
    deallocate (first_group_element, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating first_group_element in check_and_sort_similar_centres')
    deallocate (group_verifier, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating group_verifier in check_and_sort_similar_centres')
    deallocate (tmp_wf_verifier, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating tmp_wf_verifier in check_and_sort_similar_centres')
    deallocate (has_similar_centres, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating has_similar_centres in check_and_sort_similar_centres')
    deallocate (idx_similar_wf, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating idx_similar_wf in check_and_sort_similar_centres')
    deallocate (wf_similar_centres, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating wf_similar_centre in check_and_sort_similar_centres')

    if (timing_level > 2) call io_stopwatch('tran: lcr_2c2_sort: similar_centres', 2)

    return

  end subroutine check_and_sort_similar_centres

  !=====================================!
  subroutine tran_write_xyz()
    !=====================================!
    !                                     !
    ! Write xyz file with Wannier centres !
    ! and atomic positions                !
    !                                     !
    !=====================================!

    use w90_io, only: seedname, io_file_unit, io_date, stdout
    use w90_parameters, only: num_wann, &
      atoms_pos_cart, atoms_symbol, num_species, &
      atoms_species_num, num_atoms, transport_mode
    use w90_hamiltonian, only: wannier_centres_translated

    implicit none

    integer          :: iw, ind, xyz_unit, nat, nsp
    character(len=9) :: cdate, ctime
    real(kind=dp)    :: wc(3, num_wann)

    if (index(transport_mode, 'bulk') > 0) wc = wannier_centres_translated
    if (index(transport_mode, 'lcr') > 0) then
      do iw = 1, num_wann
        wc(:, iw) = wannier_centres_translated(:, tran_sorted_idx(iw))
      enddo
    endif

    xyz_unit = io_file_unit()
    open (xyz_unit, file=trim(seedname)//'_centres.xyz', form='formatted')
    !
    write (xyz_unit, '(i6)') num_wann + num_atoms
    !
    call io_date(cdate, ctime)
    write (xyz_unit, '(a84)') 'Wannier centres and atomic positions, written by Wannier90 on '//cdate//' at '//ctime
    !
    do iw = 1, num_wann
      write (xyz_unit, '("X",6x,3(f14.8,3x))') (wc(ind, iw), ind=1, 3)
    end do
    do nsp = 1, num_species
      do nat = 1, atoms_species_num(nsp)
        write (xyz_unit, '(a2,5x,3(f14.8,3x))') atoms_symbol(nsp), atoms_pos_cart(:, nat, nsp)
      end do
    end do

    write (stdout, *) ' Wannier centres written to file '//trim(seedname)//'_centres.xyz'

    return

  end subroutine tran_write_xyz

  !==============================================================!
  subroutine tran_parity_enforce(signatures)
    !==============================================================!
    ! Here, the signatures of the each wannier fucntion (stored in !
    ! signatures) is used to determine its relavite parity         !
    ! with respect to the first unit cell. The parity pattern of   !
    ! first unit cell is then enforced.                            !
    !==============================================================!

    use w90_constants, only: dp
    use w90_io, only: stdout, io_stopwatch
    use w90_parameters, only: tran_num_cell_ll, num_wann, tran_num_ll, &
      timing_level, iprint, tran_easy_fix

    implicit none

    real(dp), intent(inout), dimension(:, :)               :: signatures

    integer                                             :: i, j, k, wf_idx, num_wann_cell_ll
    real(dp)                                            :: signature_dot_p

    if (timing_level > 1) call io_stopwatch('tran: parity_enforce', 1)

    !
    ! NP: special "easy" fix of the parities by switching the sign
    ! of the Wannier Functions if the first element of the signature
    ! is found negative. Then updating the signature and the Hamiltonian
    ! matrix element for the corresponding line and column
    !
    if (tran_easy_fix) then
      do i = 1, num_wann
        if (real(signatures(1, i)) .lt. 0.0_dp) then
          signatures(:, i) = -signatures(:, i)
          do k = 1, num_wann
            hr_one_dim(k, i, 0) = -hr_one_dim(k, i, 0)
            hr_one_dim(i, k, 0) = -hr_one_dim(i, k, 0)
          enddo
        endif
      enddo
    endif

    num_wann_cell_ll = tran_num_ll/tran_num_cell_ll
    if (iprint .eq. 5) write (stdout, '(a101)') 'Unit cell    Sorted WF index    Unsort WF index  &
         &Unsorted WF Equiv       Signature Dot Product'
    !
    ! Loop over unit cell in principal layers
    !
    do i = 2, 4*tran_num_cell_ll
      !
      ! Loop over wannier functions in unit cell
      !
      do j = 1, num_wann_cell_ll
        if (i .le. 2*tran_num_cell_ll) then
          wf_idx = j + (i - 1)*num_wann_cell_ll
        else
          wf_idx = num_wann - 2*tran_num_ll + j + (i - 1 - 2*tran_num_cell_ll)*num_wann_cell_ll
        endif
        signature_dot_p = dot_product(signatures(:, tran_sorted_idx(j)), signatures(:, tran_sorted_idx(wf_idx)))
        if (iprint .eq. 5) then
          write (stdout, '(2x,i4,3(13x,i5),12x,f20.17)') &
            i, wf_idx, tran_sorted_idx(wf_idx), tran_sorted_idx(j), signature_dot_p
        endif
        if (abs(signature_dot_p) .le. 0.8_dp) then
          write (stdout, '(a28,i4,a64,i4,a20)') ' WARNING: Wannier function (', tran_sorted_idx(wf_idx), &
            ') seems to has poor resemblance to equivalent wannier function (', tran_sorted_idx(j), ') in first unit cell'
          if (iprint .lt. 5) write (stdout, *) 'Dot product of signatures: ', signature_dot_p
        endif
        if (signature_dot_p .lt. 0.0_dp) then
          do k = 1, num_wann
            hr_one_dim(k, tran_sorted_idx(wf_idx), 0) = -hr_one_dim(k, tran_sorted_idx(wf_idx), 0)
            hr_one_dim(tran_sorted_idx(wf_idx), k, 0) = -hr_one_dim(tran_sorted_idx(wf_idx), k, 0)
          enddo
        endif
      enddo
    enddo

    if (timing_level > 1) call io_stopwatch('tran: parity_enforce', 2)

    return

  end subroutine tran_parity_enforce

  !========================================!
  subroutine tran_lcr_2c2_build_ham(pl_warning)
    !==============================================!
    ! Builds hamiltonians blocks required for the  !
    ! Greens function caclulations of the quantum  !
    ! conductance according to the 2c2 geometry.   !
    ! Leads are also symmetrised, in that unit cell!
    ! sub-blocks are copied to create truely ideal !
    ! leads.                                       !
    !==============================================!

    use w90_constants, only: dp, eps5
    use w90_io, only: io_error, stdout, seedname, io_file_unit, io_date, io_stopwatch
    use w90_parameters, only: tran_num_cell_ll, num_wann, tran_num_ll, kpt_cart, nfermi, fermi_energy_list, &
      tran_write_ht, tran_num_rr, tran_num_lc, tran_num_cr, tran_num_cc, &
      tran_num_bandc, timing_level, dist_cutoff_mode, dist_cutoff, &
      dist_cutoff_hc
    use w90_hamiltonian, only: wannier_centres_translated

    implicit none

    logical, intent(in)                     :: pl_warning

    integer                                :: i, j, k, num_wann_cell_ll, file_unit, ierr, band_size

    real(dp), allocatable, dimension(:, :)    :: sub_block
    real(dp)                               :: PL_length, dist, dist_vec(3)

    character(len=9)                       :: cdate, ctime

    if (timing_level > 1) call io_stopwatch('tran: lcr_2c2_build_ham', 1)

    if (nfermi > 1) call io_error("Error in tran_lcr_2c2_build_ham: nfermi>1. " &
                                  //"Set the fermi level using the input parameter 'fermi_evel'")

    allocate (hL0(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hL0 in tran_lcr_2c2_build_ham')
    allocate (hL1(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hL1 in tran_lcr_2c2_build_ham')
    allocate (hR0(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hR0 in tran_lcr_2c2_build_ham')
    allocate (hR1(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hR1 in tran_lcr_2c2_build_ham')
    allocate (hLC(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hLC in tran_lcr_2c2_build_ham')
    allocate (hCR(tran_num_ll, tran_num_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hCR in tran_lcr_2c2_build_ham')
    allocate (hC(num_wann - (2*tran_num_ll), num_wann - (2*tran_num_ll)), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating hC in tran_lcr_2c2_build_ham')
    !
    !This checks that only the gamma point is used in wannierisation
    !This is necessary since this calculation only makes sense if we
    !have periodicity over the supercell.
    !
    if ((size(kpt_cart, 2) .ne. 1) .and. (kpt_cart(1, 1) .eq. 0.0_dp) &
        .and. (kpt_cart(2, 1) .eq. 0.0_dp) &
        .and. (kpt_cart(3, 1) .eq. 0.0_dp)) then
      call io_error('Calculation must be performed at gamma only')
    endif

    num_wann_cell_ll = tran_num_ll/tran_num_cell_ll

    allocate (sub_block(num_wann_cell_ll, num_wann_cell_ll), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating sub_block in tran_lcr_2c2_build_ham')
    !
    !Build hL0 & hL1
    !
    hL0 = 0.0_dp
    hL1 = 0.0_dp
    !
    !Loop over the sub_blocks corresponding to distinct unit cells inside the principal layer
    !
    do i = 1, tran_num_cell_ll
      !
      !Each sub_block will be duplicated along the corresponding diagonal. This ensures the correct symmetry for the leads.
      !
      sub_block = 0.0_dp
      !
      !Extract matrix elements from hr_one_dim needed for hL0 (and lower triangular sub_blocks of hL1)
      !
      do j = 1, num_wann_cell_ll
        do k = 1, num_wann_cell_ll
          sub_block(j, k) = hr_one_dim(tran_sorted_idx(j), tran_sorted_idx((i - 1)*num_wann_cell_ll + k), 0)
        enddo
      enddo
      !
      !Filling up hL0 sub_block by sub_block
      !
      do j = 1, tran_num_cell_ll - i + 1
        !
        !Fill diagonal and upper diagonal sub_blocks
        !
        hL0((j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll, &
            (j - 1)*num_wann_cell_ll + 1 + (i - 1)*num_wann_cell_ll:j*num_wann_cell_ll + (i - 1)*num_wann_cell_ll) = sub_block
        !
        !Fill lower diagonal sub_blocks
        !
        if (i .gt. 1) then
          hL0((j - 1)*num_wann_cell_ll + 1 + (i - 1)*num_wann_cell_ll:j*num_wann_cell_ll + (i - 1)*num_wann_cell_ll, &
              (j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll) = transpose(sub_block)
        endif
      enddo
      !
      !Filling up non-diagonal hL1 sub_blocks (nothing need be done for i=1)
      !
      if (i .gt. 1) then
        do j = 1, i - 1
          hL1((tran_num_cell_ll - (i - j))*num_wann_cell_ll + 1:(tran_num_cell_ll - (i - 1 - j))*num_wann_cell_ll, &
              (j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll) = sub_block
        enddo
      endif
      !
      ! MS: Get diagonal and upper triangular sublocks for hL1 - use periodic image of PL4
      !
      sub_block = 0.0_dp
      !
      if (i == 1) then !Do diagonal only
        do j = 1, num_wann_cell_ll
          do k = 1, num_wann_cell_ll
            sub_block(j, k) = hr_one_dim( &
                              tran_sorted_idx(num_wann - tran_num_ll + j), &
                              tran_sorted_idx((i - 1)*num_wann_cell_ll + k), 0)
          enddo
        enddo
        !
        ! MS: Now fill subblocks of hL1
        !
        do j = 1, tran_num_cell_ll - i + 1
          hL1((j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll, &
              (j - 1)*num_wann_cell_ll + 1 + (i - 1)*num_wann_cell_ll:j*num_wann_cell_ll + (i - 1)* &
              num_wann_cell_ll) = sub_block
        enddo
      endif
    enddo
    !
    !Special case tran_num_cell_ll=1, the diagonal sub-block of hL1 is hL1, so cannot be left as zero
    !
    if (tran_num_cell_ll .eq. 1) then
      do j = num_wann - num_wann_cell_ll + 1, num_wann
        do k = 1, num_wann_cell_ll
          hL1(j - num_wann + num_wann_cell_ll, k) = hr_one_dim(tran_sorted_idx(j), tran_sorted_idx(k), 0)
        enddo
      enddo
    endif

    !
    !Build hR0 & hR1
    !
    hR0 = 0.0_dp
    hR1 = 0.0_dp
    !
    !Loop over the sub_blocks corresponding to distinct unit cells inside the principal layer
    !
    do i = 1, tran_num_cell_ll
      !
      !Each sub_block will be duplicated along the corresponding diagonal. This ensures the correct symmetry for the leads.
      !
      sub_block = 0.0_dp
      !
      !Extract matrix elements from hr_one_dim needed for hR0 (and lower triangular sub_blocks of hR1)
      !
      do j = 1, num_wann_cell_ll
        do k = 1, num_wann_cell_ll
          sub_block(j, k) = hr_one_dim(tran_sorted_idx(num_wann - i*(num_wann_cell_ll) + j), &
                                       tran_sorted_idx(num_wann - num_wann_cell_ll + k), 0)
        enddo
      enddo
      !
      !Filling up hR0 sub_block by sub_block
      !
      do j = 1, tran_num_cell_ll - i + 1
        !
        !Fill diagonal and upper diagonal sub_blocks
        !
        hR0((j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll, &
            (j - 1)*num_wann_cell_ll + 1 + (i - 1)*num_wann_cell_ll:j*num_wann_cell_ll + (i - 1)*num_wann_cell_ll) = sub_block
        !
        !Fill lower diagonal sub_blocks
        !
        if (i .gt. 1) then
          hR0((j - 1)*num_wann_cell_ll + 1 + (i - 1)*num_wann_cell_ll:j*num_wann_cell_ll + (i - 1)*num_wann_cell_ll, &
              (j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll) = transpose(sub_block)
        endif
      enddo
      !
      !Filling up non-diagonal hR1 sub_blocks (nothing need be done for i=1)
      !
      if (i .gt. 1) then
        do j = 1, i - 1
          hR1((tran_num_cell_ll - (i - j))*num_wann_cell_ll + 1:(tran_num_cell_ll - (i - 1 - j))*num_wann_cell_ll, &
              (j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll) = sub_block
        enddo
      endif
      !
      ! MS: Get diagonal and upper triangular sublocks for hR1 - use periodic image of PL1
      !
      sub_block = 0.0_dp
      !
      if (i == 1) then  !Do diagonal only
        do j = 1, num_wann_cell_ll
          do k = 1, num_wann_cell_ll
            sub_block(j, k) = hr_one_dim(tran_sorted_idx((i - 1)*num_wann_cell_ll + k), &
                                         tran_sorted_idx(num_wann - tran_num_ll + j), 0)
          enddo
        enddo
        !
        ! MS: Now fill subblocks of hR1
        !
        do j = 1, tran_num_cell_ll - i + 1
          hR1((j - 1)*num_wann_cell_ll + 1:j*num_wann_cell_ll, &
              (j - 1)*num_wann_cell_ll + 1 + (i - 1)*num_wann_cell_ll:j*num_wann_cell_ll + (i - 1)*num_wann_cell_ll) = sub_block
        enddo
      endif
    enddo
    !
    !Special case tran_num_cell_ll=1, the diagonal sub-block of hR1 is hR1, so cannot be left as zero
    !
    if (tran_num_cell_ll .eq. 1) then
      do j = 1, num_wann_cell_ll
        do k = num_wann - num_wann_cell_ll + 1, num_wann
          hR1(k - num_wann + num_wann_cell_ll, j) = hr_one_dim(tran_sorted_idx(j), tran_sorted_idx(k), 0)
        enddo
      enddo
    endif

    !
    !Building hLC
    !
    hLC = 0.0_dp
    do i = 1, tran_num_ll
      do j = tran_num_ll + 1, 2*tran_num_ll
        hLC(i, j - tran_num_ll) = hr_one_dim(tran_sorted_idx(i), tran_sorted_idx(j), 0)
      enddo
    enddo
!----!
! MS ! Rely on dist_cutoff doing the work here, as it cuts element-wise, not block wise (incorrect)
!----!
!    if (tran_num_cell_ll .gt. 1) then
!        do j=1,tran_num_cell_ll
!            do k=1,tran_num_cell_ll
!                if (k .ge. j) then
!                    hLC((j-1)*num_wann_cell_ll+1:j*num_wann_cell_ll,(k-1)*num_wann_cell_ll+1:k*num_wann_cell_ll)=0.0_dp
!                endif
!            enddo
!        enddo
!    endif
!---!
!end!
!---!

    !
    !Building hC
    !
    hC = 0.0_dp
    !
    band_size = 0
    if (dist_cutoff_hc .ne. dist_cutoff) then
      dist_cutoff = dist_cutoff_hc
      write (stdout, *) 'Applying dist_cutoff_hc to Hamiltonian for construction of hC'
      deallocate (hr_one_dim, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating hr_one_dim in tran_lcr_2c2_sort')
      call tran_reduce_hr()
      call tran_cut_hr_one_dim()
    endif

    do i = tran_num_ll + 1, num_wann - tran_num_ll
      do j = tran_num_ll + 1, num_wann - tran_num_ll
        hC(i - tran_num_ll, j - tran_num_ll) = hr_one_dim(tran_sorted_idx(i), tran_sorted_idx(j), 0)
        !
        ! Impose a ham_cutoff of 1e-4 eV to reduce tran_num_bandc (and in turn hCband, and speed up transport)
        !
        if (abs(hC(i - tran_num_ll, j - tran_num_ll)) .lt. 10.0_dp*eps5) then
          hC(i - tran_num_ll, j - tran_num_ll) = 0.0_dp
          band_size = max(band_size, abs(i - j))
        endif
      enddo
    enddo
    !
    !Building hCR
    !
    hCR = 0.0_dp
    do i = num_wann - 2*tran_num_ll + 1, num_wann - tran_num_ll
      do j = num_wann - tran_num_ll + 1, num_wann
        hCR(i - (num_wann - 2*tran_num_ll), j - (num_wann - tran_num_ll)) = hr_one_dim(tran_sorted_idx(i), tran_sorted_idx(j), 0)
      enddo
    enddo
!----!
! MS ! Rely on dist_cutoff doing the work here, as it cuts element-wise, not block wise (incorrect)
!----!
!    if (tran_num_cell_ll .gt. 1) then
!        do j=1,tran_num_cell_ll
!            do k=1,tran_num_cell_ll
!                if (k .ge. j) then
!                    hCR((j-1)*num_wann_cell_ll+1:j*num_wann_cell_ll,(k-1)*num_wann_cell_ll+1:k*num_wann_cell_ll)=0.0_dp
!                endif
!            enddo
!        enddo
!    endif
!---!
!end!
!---!

    !
    !Subtract the Fermi energy from the diagonal elements of hC,hL0,hR0
    !
    do i = 1, tran_num_ll
      hL0(i, i) = hL0(i, i) - fermi_energy_list(1)
      hR0(i, i) = hR0(i, i) - fermi_energy_list(1)
    enddo
    do i = 1, num_wann - (2*tran_num_ll)
      hC(i, i) = hC(i, i) - fermi_energy_list(1)
    enddo
    !
    !Define tran_num_** parameters that are used later in tran_lcr
    !
    tran_num_rr = tran_num_ll
    tran_num_lc = tran_num_ll
    tran_num_cr = tran_num_ll
    tran_num_cc = num_wann - (2*tran_num_ll)

    !
    ! Set appropriate tran_num_bandc if has not been set (0.0_dp is default value)
    !
    if (tran_num_bandc .eq. 0.0_dp) then
      tran_num_bandc = min(band_size + 1, (tran_num_cc + 1)/2 + 1)
    endif

    !
    ! MS: Find and print effective PL length
    !
    if (.not. pl_warning) then
      PL_length = 0.0_dp
      do i = 1, tran_num_ll
        do j = 1, tran_num_ll
          if (abs(hL1(i, j)) .gt. 0.0_dp) then
            if (index(dist_cutoff_mode, 'one_dim') .gt. 0) then
              dist = abs(wannier_centres_translated(coord(1), tran_sorted_idx(i)) &
                         - wannier_centres_translated(coord(1), tran_sorted_idx(j + tran_num_ll)))
            else
              dist_vec(:) = wannier_centres_translated(:, tran_sorted_idx(i)) &
                            - wannier_centres_translated(:, tran_sorted_idx(j + tran_num_ll))
              dist = sqrt(dot_product(dist_vec, dist_vec))
            endif
            PL_length = max(PL_length, dist)
          endif
          if (abs(hR1(i, j)) .gt. 0.0_dp) then
            if (index(dist_cutoff_mode, 'one_dim') .gt. 0) then
              dist = abs(wannier_centres_translated(coord(1), tran_sorted_idx(num_wann - 2*tran_num_ll + i)) &
                         - wannier_centres_translated(coord(1), tran_sorted_idx(num_wann - tran_num_ll + j)))
            else
              dist_vec(:) = wannier_centres_translated(:, tran_sorted_idx(num_wann - 2*tran_num_ll + i)) &
                            - wannier_centres_translated(:, tran_sorted_idx(num_wann - tran_num_ll + j))
              dist = sqrt(dot_product(dist_vec, dist_vec))
            endif
            PL_length = max(PL_length, dist)
          endif
        enddo
      enddo
      write (stdout, '(1x,a,f12.6,a)') 'Approximate effective principal layer length is: ', PL_length, ' Ang.'
    endif

    !
    !Writing to file:
    !
    if (tran_write_ht) then
      write (stdout, *) '------------------------------- Writing ht files  ----------------------------'
      !
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_htL.dat', status='unknown', form='formatted', action='write')

      call io_date(cdate, ctime)
      write (file_unit, *) 'written on '//cdate//' at '//ctime ! Date and time
      write (file_unit, '(I6)') tran_num_ll
      write (file_unit, '(6F12.6)') ((hL0(j, i), j=1, tran_num_ll), i=1, tran_num_ll)
      write (file_unit, '(I6)') tran_num_ll
      write (file_unit, '(6F12.6)') ((hL1(j, i), j=1, tran_num_ll), i=1, tran_num_ll)

      close (file_unit)
      write (stdout, *) ' '//trim(seedname)//'_htL.dat  written'
      !
      !hR
      !
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_htR.dat', status='unknown', form='formatted', action='write')

      call io_date(cdate, ctime)
      write (file_unit, *) 'written on '//cdate//' at '//ctime ! Date and time
      write (file_unit, '(I6)') tran_num_rr
      write (file_unit, '(6F12.6)') ((hR0(j, i), j=1, tran_num_rr), i=1, tran_num_rr)
      write (file_unit, '(I6)') tran_num_rr
      write (file_unit, '(6F12.6)') ((hR1(j, i), j=1, tran_num_rr), i=1, tran_num_rr)

      close (file_unit)
      write (stdout, *) ' '//trim(seedname)//'_htR.dat  written'
      !
      !hLC
      !
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_htLC.dat', status='unknown', form='formatted', action='write')

      call io_date(cdate, ctime)
      write (file_unit, *) 'written on '//cdate//' at '//ctime ! Date and time
      write (file_unit, '(2I6)') tran_num_ll, tran_num_lc
      write (file_unit, '(6F12.6)') ((hLC(j, i), j=1, tran_num_lc), i=1, tran_num_lc)

      close (file_unit)
      write (stdout, *) ' '//trim(seedname)//'_htLC.dat written'
      !
      !hCR
      !
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_htCR.dat', status='unknown', form='formatted', action='write')

      call io_date(cdate, ctime)
      write (file_unit, *) 'written on '//cdate//' at '//ctime ! Date and time
      write (file_unit, '(2I6)') tran_num_cr, tran_num_rr
      write (file_unit, '(6F12.6)') ((hCR(j, i), j=1, tran_num_cr), i=1, tran_num_cr)

      close (file_unit)
      write (stdout, *) ' '//trim(seedname)//'_htCR.dat written'
      !
      !hC
      !
      file_unit = io_file_unit()
      open (file_unit, file=trim(seedname)//'_htC.dat', status='unknown', form='formatted', action='write')

      call io_date(cdate, ctime)
      write (file_unit, *) 'written on '//cdate//' at '//ctime ! Date and time
      write (file_unit, '(I6)') tran_num_cc
      write (file_unit, '(6F12.6)') ((hC(j, i), j=1, tran_num_cc), i=1, tran_num_cc)

      close (file_unit)
      write (stdout, *) ' '//trim(seedname)//'_htC.dat  written'

      write (stdout, *) '------------------------------------------------------------------------------'
    end if

    deallocate (sub_block, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating sub_block in tran_lcr_2c2_build_ham')

    if (timing_level > 1) call io_stopwatch('tran: lcr_2c2_build_ham', 2)

    return

  end subroutine tran_lcr_2c2_build_ham

  !======================================!
  subroutine tran_dealloc()
    !! Dellocate module data
    !====================================!

    use w90_io, only: io_error

    implicit none

    integer :: ierr

    if (allocated(hR1)) then
      deallocate (hR1, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating hR1 in tran_dealloc')
    end if
    if (allocated(hR0)) then
      deallocate (hR0, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating hR0 in tran_dealloc')
    end if
    if (allocated(hL1)) then
      deallocate (hL1, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating hL1 in tran_dealloc')
    end if
    if (allocated(hB1)) then
      deallocate (hB1, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating hB1 in tran_dealloc')
    end if
    if (allocated(hB0)) then
      deallocate (hB0, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating hB0 in tran_dealloc')
    end if
    if (allocated(hr_one_dim)) then
      deallocate (hr_one_dim, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating hr_one_dim in tran_dealloc')
    end if

    return

  end subroutine tran_dealloc

end module w90_transport

