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
!  w90_geninterp: interpolation functions                    !
!                                                            !
!------------------------------------------------------------!

module w90_geninterp
  !! Generic Interpolation Routine
  !!
  !! written by Giovanni Pizzi
  !! THEOS, EPFL, Station 12, 1015 Lausanne (Switzerland)
  !! June, 2012

  implicit none

  private

  public :: geninterp_main

contains

  subroutine internal_write_header(outdat_unit, commentline, pw90_geninterp)
    !! Writes a header for the output file(s).

    use w90_postw90_types, only: pw90_geninterp_mod_type
    use w90_io, only: io_date

    ! arguments
    type(pw90_geninterp_mod_type), intent(in) :: pw90_geninterp
    integer, intent(in) :: outdat_unit
    !! Integer with the output file unit. The file must be already open.
    character(len=*) :: commentline !! no intent?
    !! String with the comment taken from the output, to be written on the output

    ! local variables
    character(len=9) :: cdate, ctime

    call io_date(cdate, ctime)
    write (outdat_unit, '(A)') "# Written on "//cdate//" at "//ctime ! Date and time
    ! I rewrite the comment line on the output
    write (outdat_unit, '(A)') "# Input file comment: "//trim(commentline)

    if (pw90_geninterp%alsofirstder) then
      write (outdat_unit, '(A)') "#  Kpt_idx  K_x (1/ang)       K_y (1/ang)        K_z (1/ang)       Energy (eV)"// &
        "      EnergyDer_x       EnergyDer_y       EnergyDer_z"
    else
      write (outdat_unit, '(A)') "#  Kpt_idx  K_x (1/ang)       K_y (1/ang)        K_z (1/ang)       Energy (eV)"
    end if
  end subroutine internal_write_header

  subroutine geninterp_main(dis_manifold, pw90_geninterp, kpt_latt, pw90_band_deriv_degen, ws_region, print_output, &
                            wannier_data, ws_distance, wigner_seitz, HH_R, v_matrix, u_matrix, eigval, &
                            real_lattice, scissors_shift, mp_grid, num_bands, &
                            num_kpts, num_wann, num_valence_bands, effective_model, &
                            have_disentangled, seedname, stdout, comm)

    !! This routine prints the band energies (and possibly the band derivatives)
    !!
    !! This routine is parallel, even if ***the scaling is very bad*** since at the moment
    !! everything must be written by the root node (we need that the output is sorted).
    !! But at least if works independently of the number of processors.
    !! I think that a way to write in parallel to the output would help a lot,
    !! so that we don't have to send all eigenvalues to the root node.
    use w90_constants, only: dp, pi
    use w90_postw90_types, only: pw90_geninterp_mod_type, &
      pw90_band_deriv_degen_type, wigner_seitz_type
    use w90_types, only: dis_manifold_type, print_output_type, &
      wannier_data_type, ws_region_type, ws_distance_type
    use w90_io, only: io_error, io_stopwatch, io_file_unit, io_stopwatch
    use w90_postw90_common, only: pw90common_fourier_R_to_k
    use w90_utility, only: utility_diagonalize, utility_recip_lattice_base
    use w90_wan_ham, only: wham_get_eig_deleig
    use w90_get_oper, only: get_HH_R
    use w90_comms, only: mpirank, mpisize, comms_bcast, comms_array_split, comms_scatterv, &
      comms_gatherv, w90comm_type

    ! arguments
    type(dis_manifold_type), intent(in)          :: dis_manifold
    type(pw90_geninterp_mod_type), intent(in)    :: pw90_geninterp
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(ws_region_type), intent(in)             :: ws_region
    type(print_output_type), intent(in)          :: print_output
    type(wannier_data_type), intent(in)          :: wannier_data
    type(ws_distance_type), intent(inout)        :: ws_distance
    type(wigner_seitz_type), intent(inout)       :: wigner_seitz
    type(w90comm_type), intent(in)               :: comm

    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :)
    complex(kind=dp), intent(in) :: v_matrix(:, :, :), u_matrix(:, :, :)

    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift
    real(kind=dp), intent(in) :: kpt_latt(:, :)

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_bands, num_kpts, num_wann, num_valence_bands, stdout

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    ! local variables
    real(kind=dp) :: recip_lattice(3, 3), volume
    integer :: kpt_unit, outdat_unit, ierr, i, j, enidx
    integer :: nkinterp ! number of kpoints for which we perform the interpolation
    integer, allocatable :: counts(:), displs(:)
    integer :: my_node_id, num_nodes
    integer, allocatable :: kpointidx(:), localkpointidx(:)
    real(kind=dp), allocatable :: kpoints(:, :), localkpoints(:, :)
    real(kind=dp)  :: kpt(3), frac(3)
    real(kind=dp), allocatable :: localdeleig(:, :, :)
    real(kind=dp), allocatable :: globaldeleig(:, :, :)
    real(kind=dp), allocatable :: localeig(:, :)
    real(kind=dp), allocatable :: globaleig(:, :)
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    character(len=500) :: commentline
    character(len=50) :: cdum
    character(len=200) :: outdat_filename
    logical :: absoluteCoords
    logical :: on_root = .false.

    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)
    if (my_node_id == 0) on_root = .true.
    allocate (counts(0:num_nodes - 1))
    allocate (displs(0:num_nodes - 1))

    if (print_output%iprint > 0 .and. (print_output%timing_level > 0)) &
      call io_stopwatch('geninterp_main', 1, stdout, seedname)

    if (on_root) then
      write (stdout, *)
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '|                      Generic Band Interpolation routines                  |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'

      kpt_unit = io_file_unit()
      open (unit=kpt_unit, file=trim(seedname)//'_geninterp.kpt', form='formatted', status='old', &
            err=105)

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
                      'unable to recognize keyword', stdout, seedname)
      end if

      ! Third line: number of following kpoints
      read (kpt_unit, *, err=106, end=106) nkinterp
    end if

    call comms_bcast(nkinterp, 1, stdout, seedname, comm)

    allocate (HH(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating HH in calcTDF', stdout, seedname)
    allocate (UU(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating UU in calcTDF', stdout, seedname)
    if (pw90_geninterp%alsofirstder) then
      allocate (delHH(num_wann, num_wann, 3), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating delHH in calcTDF', stdout, seedname)
    end if

    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>
    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, eigval, &
                  real_lattice, scissors_shift, num_bands, num_kpts, num_wann, num_valence_bands, &
                  effective_model, have_disentangled, seedname, stdout, comm)

    if (on_root) then
      allocate (kpointidx(nkinterp), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpointidx in geinterp_main.', stdout, seedname)
      allocate (kpoints(3, nkinterp), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpoints in geinterp_main.', stdout, seedname)
      if (pw90_geninterp%single_file) then
        allocate (globaleig(num_wann, nkinterp), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaleig in geinterp_main.', stdout, &
                                     seedname)
        allocate (globaldeleig(num_wann, 3, nkinterp), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaldeleig in geinterp_main.', stdout, &
                                     seedname)
      end if
    else
      ! On the other nodes, I still allocate them with size 1 to avoid
      ! that some compilers still try to access the memory
      allocate (kpointidx(1), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpointidx in geinterp_main.', stdout, seedname)
      allocate (kpoints(1, 1), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpoints in geinterp_main.', stdout, seedname)
      if (pw90_geninterp%single_file) then
        allocate (globaleig(num_wann, 1), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaleig in geinterp_main.', stdout, &
                                     seedname)
        allocate (globaldeleig(num_wann, 3, 1), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating globaldeleig in geinterp_main.', stdout, &
                                     seedname)
      end if
    end if

    ! I precalculate how to split on different nodes
    call comms_array_split(nkinterp, counts, displs, comm)

    allocate (localkpoints(3, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating localkpoints in geinterp_main.', stdout, &
                                 seedname)

    allocate (localeig(num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating localeig in geinterp_main.', stdout, seedname)
    allocate (localdeleig(num_wann, 3, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating localdeleig in geinterp_main.', stdout, seedname)

    ! On root, I read numpoints_thischunk points
    if (on_root) then
      ! Lines with integer identifier and three coordinates
      ! (in crystallographic coordinates relative to the reciprocal lattice vectors)
      do i = 1, nkinterp
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
    call comms_scatterv(localkpoints, 3*counts(my_node_id), kpoints, 3*counts, 3*displs, stdout, &
                        seedname, comm)
    if (.not. pw90_geninterp%single_file) then
      ! Allocate at least one entry, even if we don't use it
      allocate (localkpointidx(max(1, counts(my_node_id))), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating localkpointidx in geinterp_main.', stdout, &
                                   seedname)
      call comms_scatterv(localkpointidx, counts(my_node_id), kpointidx, counts, displs, stdout, &
                          seedname, comm)
    end if

    ! I open the output file(s)
    if (pw90_geninterp%single_file) then
      if (on_root) then
        outdat_filename = trim(seedname)//'_geninterp.dat'
        outdat_unit = io_file_unit()
        open (unit=outdat_unit, file=trim(outdat_filename), form='formatted', err=107)

        call internal_write_header(outdat_unit, commentline, pw90_geninterp)
      end if
    else
      if (num_nodes > 99999) then
        write (outdat_filename, '(a,a,I0,a)') trim(seedname), '_geninterp_', my_node_id, '.dat'
      else
        write (outdat_filename, '(a,a,I5.5,a)') trim(seedname), '_geninterp_', my_node_id, '.dat'
      endif
      outdat_unit = io_file_unit()
      open (unit=outdat_unit, file=trim(outdat_filename), form='formatted', err=107)

      call comms_bcast(commentline, len(commentline), stdout, seedname, comm)

      call internal_write_header(outdat_unit, commentline, pw90_geninterp)
    end if

    ! And now, each node calculates its own k points
    do i = 1, counts(my_node_id)
      kpt = localkpoints(:, i)
      ! Here I get the band energies and the velocities (if required)
      if (pw90_geninterp%alsofirstder) then
        call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, print_output, wannier_data, &
                                 ws_distance, wigner_seitz, delHH, HH, HH_R, u_matrix, UU, v_matrix, &
                                 localdeleig(:, :, i), localeig(:, i), eigval, kpt, real_lattice, &
                                 scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                                 num_valence_bands, effective_model, have_disentangled, &
                                 seedname, stdout, comm)
      else
        call pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, HH, HH_R, kpt, &
                                       real_lattice, mp_grid, 0, num_wann, seedname, stdout)
        call utility_diagonalize(HH, num_wann, localeig(:, i), UU, stdout, seedname)
      end if
    end do

    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
    if (pw90_geninterp%single_file) then
      ! Now, I get the results from the different nodes
      call comms_gatherv(localeig, num_wann*counts(my_node_id), globaleig, &
                         num_wann*counts, num_wann*displs, stdout, seedname, comm)

      if (pw90_geninterp%alsofirstder) then
        call comms_gatherv(localdeleig, 3*num_wann*counts(my_node_id), globaldeleig, &
                           3*num_wann*counts, 3*num_wann*displs, stdout, seedname, comm)
      end if

      ! Now the printing, only on root node
      if (on_root) then
        do i = 1, nkinterp
          kpt = kpoints(:, i)
          ! First calculate the absolute coordinates for printing
          frac = 0._dp
          do j = 1, 3
            frac(j) = recip_lattice(1, j)*kpt(1) + recip_lattice(2, j)*kpt(2) &
                      + recip_lattice(3, j)*kpt(3)
          end do

          ! I print each line
          if (pw90_geninterp%alsofirstder) then
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
          frac(j) = recip_lattice(1, j)*kpt(1) + recip_lattice(2, j)*kpt(2) &
                    + recip_lattice(3, j)*kpt(3)
        end do

        ! I print each line
        if (pw90_geninterp%alsofirstder) then
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

    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch('geninterp_main', 2, &
                                                                         stdout, seedname)

    return

105 call io_error('Error: Problem opening k-point file '//trim(seedname)//'_geninterp.kpt', &
                  stdout, seedname)
106 call io_error('Error: Problem reading k-point file '//trim(seedname)//'_geninterp.kpt', &
                  stdout, seedname)
107 call io_error('Error: Problem opening output file '//trim(outdat_filename), stdout, seedname)
  end subroutine geninterp_main

end module w90_geninterp
