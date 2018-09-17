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
!                                                            !
!------------------------------------------------------------!

module w90_geninterp
  !! Generic Interpolation Routine
  !!
  !! written by Giovanni Pizzi
  !! THEOS, EPFL, Station 12, 1015 Lausanne (Switzerland)
  !! June, 2012

  use w90_constants
  use w90_parameters, only: geninterp_alsofirstder, num_wann, recip_lattice, real_lattice, &
    timing_level, geninterp_single_file
  use w90_io, only: io_error, stdout, io_stopwatch, io_file_unit, seedname, io_stopwatch
  use w90_get_oper, only: get_HH_R, HH_R
  use w90_comms
  use w90_utility, only: utility_diagonalize
  use w90_postw90_common, only: pw90common_fourier_R_to_k
  use w90_wan_ham, only: wham_get_eig_deleig
  use w90_io, only: io_date
  implicit none

  private
  public :: geninterp_main

contains

  subroutine internal_write_header(outdat_unit, commentline)
    !! Writes a header for the output file(s).

    integer, intent(in) :: outdat_unit
    !! Integer with the output file unit. The file must be already open.
    character(len=*)    :: commentline !! no intent?
    !! String with the comment taken from the output, to be written on the output

    character(len=9)   :: cdate, ctime

    call io_date(cdate, ctime)
    write (outdat_unit, '(A)') "# Written on "//cdate//" at "//ctime ! Date and time
    ! I rewrite the comment line on the output
    write (outdat_unit, '(A)') "# Input file comment: "//trim(commentline)

    if (geninterp_alsofirstder) then
      write (outdat_unit, '(A)') "#  Kpt_idx  K_x (1/ang)       K_y (1/ang)        K_z (1/ang)       Energy (eV)"// &
        "      EnergyDer_x       EnergyDer_y       EnergyDer_z"
    else
      write (outdat_unit, '(A)') "#  Kpt_idx  K_x (1/ang)       K_y (1/ang)        K_z (1/ang)       Energy (eV)"
    end if
  end subroutine internal_write_header

  subroutine geninterp_main()
    !! This routine prints the band energies (and possibly the band derivatives)
    !!
    !! This routine is parallel, even if ***the scaling is very bad*** since at the moment
    !! everything must be written by the root node (we need that the output is sorted).
    !! But at least if works independently of the number of processors.
    !! I think that a way to write in parallel to the output would help a lot,
    !! so that we don't have to send all eigenvalues to the root node.
    integer            :: kpt_unit, outdat_unit, num_kpts, ierr, i, j, k, enidx
    character(len=500) :: commentline
    character(len=50)  :: cdum
    integer, dimension(:), allocatable              :: kpointidx, localkpointidx
    real(kind=dp), dimension(:, :), allocatable      :: kpoints, localkpoints
    complex(kind=dp), dimension(:, :), allocatable   :: HH
    complex(kind=dp), dimension(:, :), allocatable   :: UU
    complex(kind=dp), dimension(:, :, :), allocatable :: delHH
    real(kind=dp), dimension(3)                     :: kpt, frac
    real(kind=dp), dimension(:, :, :), allocatable    :: localdeleig
    real(kind=dp), dimension(:, :, :), allocatable    :: globaldeleig
    real(kind=dp), dimension(:, :), allocatable      :: localeig
    real(kind=dp), dimension(:, :), allocatable      :: globaleig
    logical                                         :: absoluteCoords
    character(len=200)                              :: outdat_filename

    integer, dimension(0:num_nodes - 1)               :: counts
    integer, dimension(0:num_nodes - 1)               :: displs

    if (on_root .and. (timing_level > 0)) call io_stopwatch('geninterp_main', 1)

    if (on_root) then
      write (stdout, *)
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '|                      Generic Band Interpolation routines                  |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'

      kpt_unit = io_file_unit()
      open (unit=kpt_unit, file=trim(seedname)//'_geninterp.kpt', form='formatted', status='old', err=105)

      ! First line: comment (e.g. creation date, author, ...)
      read (kpt_unit, '(A500)', err=106, end=106) commentline
      read (kpt_unit, *, err=106, end=106) cdum

      if (index(cdum, 'crystal') > 0) then
        absoluteCoords = .false.
      elseif (index(cdum, 'frac') > 0) then
        absoluteCoords = .false.
      elseif (index(cdum, 'cart') > 0) then
        absoluteCoords = .true.
      elseif (index(cdum, 'abs') > 0) then
        absoluteCoords = .true.
      else
        call io_error('Error on second line of file '//trim(seedname)//'_geninterp.kpt: '// &
                      'unable to recognize keyword')
      end if

      ! Third line: number of following kpoints
      read (kpt_unit, *, err=106, end=106) num_kpts
    end if

    call comms_bcast(num_kpts, 1)

    allocate (HH(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating HH in calcTDF')
    allocate (UU(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating UU in calcTDF')
    if (geninterp_alsofirstder) then
      allocate (delHH(num_wann, num_wann, 3), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating delHH in calcTDF')
    end if

    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>
    call get_HH_R

    if (on_root) then
      allocate (kpointidx(num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpointidx in geinterp_main.')
      allocate (kpoints(3, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpoints in geinterp_main.')
      if (geninterp_single_file) then
        allocate (globaleig(num_wann, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaleig in geinterp_main.')
        allocate (globaldeleig(num_wann, 3, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaldeleig in geinterp_main.')
      end if
    else
      ! On the other nodes, I still allocate them with size 1 to avoid
      ! that some compilers still try to access the memory
      allocate (kpointidx(1), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpointidx in geinterp_main.')
      allocate (kpoints(1, 1), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpoints in geinterp_main.')
      if (geninterp_single_file) then
        allocate (globaleig(num_wann, 1), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaleig in geinterp_main.')
        allocate (globaldeleig(num_wann, 3, 1), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaldeleig in geinterp_main.')
      end if
    end if

    ! I precalculate how to split on different nodes
    call comms_array_split(num_kpts, counts, displs)

    allocate (localkpoints(3, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating localkpoints in geinterp_main.')

    allocate (localeig(num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating localeig in geinterp_main.')
    allocate (localdeleig(num_wann, 3, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating localdeleig in geinterp_main.')

    ! On root, I read numpoints_thischunk points
    if (on_root) then
      ! Lines with integer identifier and three coordinates
      ! (in crystallographic coordinates relative to the reciprocal lattice vectors)
      do i = 1, num_kpts
        read (kpt_unit, *, err=106, end=106) kpointidx(i), kpt
        ! Internally, I need the relative (fractional) coordinates in units of the reciprocal-lattice vectors
        if (absoluteCoords .eqv. .false.) then
          kpoints(:, i) = kpt
        else
          kpoints(:, i) = 0._dp
          ! I use the real_lattice (transposed) and a factor of 2pi instead of inverting again recip_lattice
          do j = 1, 3
            kpoints(j, i) = real_lattice(j, 1)*kpt(1) + real_lattice(j, 2)*kpt(2) + &
                            real_lattice(j, 3)*kpt(3)
          end do
          kpoints(:, i) = kpoints(:, i)/(2._dp*pi)
        end if
      end do
      close (kpt_unit)
    end if

    ! Now, I distribute the kpoints; 3* because I send kx, ky, kz
    call comms_scatterv(localkpoints, 3*counts(my_node_id), kpoints, 3*counts, 3*displs)
    if (.not. geninterp_single_file) then
      ! Allocate at least one entry, even if we don't use it
      allocate (localkpointidx(max(1, counts(my_node_id))), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating localkpointidx in geinterp_main.')
      call comms_scatterv(localkpointidx, counts(my_node_id), kpointidx, counts, displs)
    end if

    ! I open the output file(s)
    if (geninterp_single_file) then
      if (on_root) then
        outdat_filename = trim(seedname)//'_geninterp.dat'
        outdat_unit = io_file_unit()
        open (unit=outdat_unit, file=trim(outdat_filename), form='formatted', err=107)

        call internal_write_header(outdat_unit, commentline)
      end if
    else
      if (num_nodes > 99999) then
        write (outdat_filename, '(a,a,I0,a)') trim(seedname), '_geninterp_', my_node_id, '.dat'
      else
        write (outdat_filename, '(a,a,I5.5,a)') trim(seedname), '_geninterp_', my_node_id, '.dat'
      endif
      outdat_unit = io_file_unit()
      open (unit=outdat_unit, file=trim(outdat_filename), form='formatted', err=107)

      call comms_bcast(commentline, len(commentline))

      call internal_write_header(outdat_unit, commentline)
    end if

    ! And now, each node calculates its own k points
    do i = 1, counts(my_node_id)
      kpt = localkpoints(:, i)
      ! Here I get the band energies and the velocities (if required)
      if (geninterp_alsofirstder) then
        call wham_get_eig_deleig(kpt, localeig(:, i), localdeleig(:, :, i), HH, delHH, UU)
      else
        call pw90common_fourier_R_to_k(kpt, HH_R, HH, 0)
        call utility_diagonalize(HH, num_wann, localeig(:, i), UU)
      end if
    end do

    if (geninterp_single_file) then
      ! Now, I get the results from the different nodes
      call comms_gatherv(localeig, num_wann*counts(my_node_id), globaleig, &
                         num_wann*counts, num_wann*displs)

      if (geninterp_alsofirstder) then
        call comms_gatherv(localdeleig, 3*num_wann*counts(my_node_id), globaldeleig, &
                           3*num_wann*counts, 3*num_wann*displs)
      end if

      ! Now the printing, only on root node
      if (on_root) then
        do i = 1, num_kpts
          kpt = kpoints(:, i)
          ! First calculate the absolute coordinates for printing
          frac = 0._dp
          do j = 1, 3
            frac(j) = recip_lattice(1, j)*kpt(1) + recip_lattice(2, j)*kpt(2) + recip_lattice(3, j)*kpt(3)
          end do

          ! I print each line
          if (geninterp_alsofirstder) then
            do enidx = 1, num_wann
              write (outdat_unit, '(I10,7G18.10)') kpointidx(i), frac, &
                globaleig(enidx, i), globaldeleig(enidx, :, i)
            end do
          else
            do enidx = 1, num_wann
              write (outdat_unit, '(I10,4G18.10)') kpointidx(i), frac, globaleig(enidx, i)
            end do
          end if
        end do
        close (outdat_unit)
      end if
    else
      ! Each node simply writes to its own file
      do i = 1, counts(my_node_id)
        kpt = localkpoints(:, i)
        ! First calculate the absolute coordinates for printing
        frac = 0._dp
        do j = 1, 3
          frac(j) = recip_lattice(1, j)*kpt(1) + recip_lattice(2, j)*kpt(2) + recip_lattice(3, j)*kpt(3)
        end do

        ! I print each line
        if (geninterp_alsofirstder) then
          do enidx = 1, num_wann
            write (outdat_unit, '(I10,7G18.10)') localkpointidx(i), frac, &
              localeig(enidx, i), localdeleig(enidx, :, i)
          end do
        else
          do enidx = 1, num_wann
            write (outdat_unit, '(I10,4G18.10)') localkpointidx(i), frac, localeig(enidx, i)
          end do
        end if
      end do
      close (outdat_unit)
    end if

    ! All k points processed: Final processing/deallocations

    if (on_root) then
      write (stdout, '(1x,a)') '|                                All done.                                  |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'

    end if

    if (allocated(kpointidx)) deallocate (kpointidx)
    if (allocated(kpoints)) deallocate (kpoints)
    if (allocated(localkpoints)) deallocate (localkpoints)
    if (allocated(HH)) deallocate (HH)
    if (allocated(UU)) deallocate (UU)
    if (allocated(delHH)) deallocate (delHH)
    if (allocated(localeig)) deallocate (localeig)
    if (allocated(localdeleig)) deallocate (localdeleig)
    if (allocated(globaleig)) deallocate (globaleig)
    if (allocated(globaldeleig)) deallocate (globaldeleig)

    if (on_root .and. (timing_level > 0)) call io_stopwatch('geninterp_main', 2)

    return

105 call io_error('Error: Problem opening k-point file '//trim(seedname)//'_geninterp.kpt')
106 call io_error('Error: Problem reading k-point file '//trim(seedname)//'_geninterp.kpt')
107 call io_error('Error: Problem opening output file '//trim(outdat_filename))
  end subroutine geninterp_main

end module w90_geninterp
