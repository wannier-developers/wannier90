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

module w90_hamiltonian
  !! Module to obtain the Hamiltonian in a wannier basis
  !! This is a simplified routine, more sophisticated properties
  !! are found in postw90 (e.g. w90_get_oper)
  use w90_constants, only: dp
  use w90_comms, only: on_root

  implicit none

! private
  !
! complex(kind=dp), public, save, allocatable :: ham_r(:, :, :)
  !! Hamiltonian matrix in WF representation
  !
! integer, public, save, allocatable :: irvec(:, :)
  !!  The irpt-th Wigner-Seitz grid point has components
  !! irvec(1:3,irpt) in the basis of the lattice vectors
  !
! integer, public, save, allocatable :: shift_vec(:, :) ? not needed allocation here
  !
! integer, public, save, allocatable :: ndegen(:)
  !! Weight of the irpt-th point is 1/ndegen(irpt)
  !
! integer, public, save              :: nrpts
  !! number of Wigner-Seitz grid points
  !
! integer, public, save              :: rpt_origin
  !! index of R=0
  !
! real(kind=dp), public, save, allocatable :: wannier_centres_translated(:, :)
  !! translated Wannier centres

  public :: hamiltonian_get_hr
  public :: hamiltonian_write_hr
  public :: hamiltonian_setup
  public :: hamiltonian_dealloc
  public :: hamiltonian_write_rmn
  public :: hamiltonian_write_tb

  type ham_logical
    logical :: ham_have_setup = .false.
    logical :: have_translated = .false.
    logical :: use_translation = .false.
    logical :: have_ham_r = .false.
    logical :: have_ham_k = .false.
    logical :: hr_written = .false.
    logical :: tb_written = .false.
  end type ham_logical

contains

  !============================================!
  subroutine hamiltonian_setup(ws_distance_tol, ws_search_size, real_metric, &
                               mp_grid, transport_mode, bands_plot_mode, transport, bands_plot, num_kpts, &
                               num_wann, timing_level, iprint, ham_r, irvec, ndegen, nrpts, rpt_origin, &
                               wannier_centres_translated, hmlg, ham_k)
    !! Allocate arrays and setup data
    !============================================!

    use w90_constants, only: cmplx_0
    use w90_io, only: io_error

    implicit none

    ! passed variables
    integer, intent(in) :: iprint
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: timing_level
    integer, intent(in) :: ws_search_size(3)
    integer, intent(inout) :: nrpts
    integer, intent(inout) :: rpt_origin
    integer, intent(inout), allocatable :: irvec(:, :)
    integer, intent(inout), allocatable :: ndegen(:)

    logical, intent(in) :: bands_plot
    logical, intent(in) :: transport

    real(kind=dp), intent(in) :: real_metric(3, 3)
    real(kind=dp), intent(in) :: ws_distance_tol
    real(kind=dp), intent(inout), allocatable :: wannier_centres_translated(:, :)

    character(len=20), intent(in) :: bands_plot_mode
    character(len=20), intent(in) :: transport_mode

    complex(kind=dp), intent(inout), allocatable :: ham_k(:, :, :)
    complex(kind=dp), intent(inout), allocatable :: ham_r(:, :, :)

    type(ham_logical), intent(inout) :: hmlg

    ! local variables
    integer :: ierr

    if (hmlg%ham_have_setup) return

    !
    ! Determine whether to use translation
    !
    if (bands_plot .and. (index(bands_plot_mode, 'cut') .ne. 0)) hmlg%use_translation = .true.
    if (transport .and. (index(transport_mode, 'bulk') .ne. 0)) hmlg%use_translation = .true.
    if (transport .and. (index(transport_mode, 'lcr') .ne. 0)) hmlg%use_translation = .true.
    !
    ! Set up Wigner-Seitz vectors
    !
    call hamiltonian_wigner_seitz(rpt_origin, nrpts, ndegen, irvec, mp_grid, real_metric, ws_search_size, ws_distance_tol, &
                                  timing_level, iprint, count_pts=.true.)
    !
    allocate (irvec(3, nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating irvec in hamiltonian_setup')
    irvec = 0
    !
    allocate (ndegen(nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ndegen in hamiltonian_setup')
    ndegen = 0
    !
    allocate (ham_r(num_wann, num_wann, nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_r in hamiltonian_setup')
    ham_r = cmplx_0
    !
    allocate (ham_k(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_k in hamiltonian_setup')
    ham_k = cmplx_0
    !
    ! Set up the wigner_seitz vectors
    !
    call hamiltonian_wigner_seitz(rpt_origin, nrpts, ndegen, irvec, mp_grid, real_metric, ws_search_size, ws_distance_tol, &
                                  timing_level, iprint, count_pts=.false.)
    !
    allocate (wannier_centres_translated(3, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating wannier_centres_translated in hamiltonian_setup')
    wannier_centres_translated = 0.0_dp

    hmlg%ham_have_setup = .true.

    return
  end subroutine hamiltonian_setup

  !============================================!
  subroutine hamiltonian_dealloc(ham_r, irvec, ndegen, wannier_centres_translated, &
                                 hmlg, ham_k)
    !! Deallocate module data
    !============================================!

    use w90_io, only: io_error

    implicit none

    complex(kind=dp), allocatable, intent(inout) :: ham_k(:, :, :)
!   logical, intent(inout) :: ham_have_setup
!   logical, intent(inout) :: have_translated
!   logical, intent(inout) :: use_translation
    type(ham_logical), intent(inout) :: hmlg

    ! passed variables
    integer, intent(inout), allocatable :: ndegen(:)
    integer, intent(inout), allocatable :: irvec(:, :)

    real(kind=dp), intent(inout), allocatable :: wannier_centres_translated(:, :)

    complex(kind=dp), intent(inout), allocatable :: ham_r(:, :, :)

    ! local variables
    integer :: ierr

    if (allocated(ham_r)) then
      deallocate (ham_r, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ham_r in hamiltonian_dealloc')
    end if
    if (allocated(ham_k)) then
      deallocate (ham_k, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ham_k in hamiltonian_dealloc')
    end if
    if (allocated(irvec)) then
      deallocate (irvec, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating irvec in hamiltonian_dealloc')
    end if
    if (allocated(ndegen)) then
      deallocate (ndegen, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating ndegen in hamiltonian_dealloc')
    end if
    if (allocated(wannier_centres_translated)) then
      deallocate (wannier_centres_translated, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating wannier_centres_translated in param_dealloc')
    end if

    hmlg%ham_have_setup = .false.
    hmlg%have_translated = .false.
    hmlg%use_translation = .false.
    hmlg%have_ham_r = .false.
    hmlg%have_ham_k = .false.
    hmlg%hr_written = .false.
    hmlg%tb_written = .false.

    return
  end subroutine hamiltonian_dealloc

  !============================================!
  subroutine hamiltonian_get_hr(real_lattice, recip_lattice, wannier_centres, &
                                num_atoms, atoms_pos_cart, translation_centre_frac, automatic_translation, &
                                num_species, atoms_species_num, lenconfac, have_disentangled, ndimwin, &
                                lwindow, u_matrix_opt, kpt_latt, eigval, u_matrix, lsitesymmetry, num_bands, &
                                num_kpts, num_wann, timing_level, ham_r, irvec, shift_vec, nrpts, &
                                wannier_centres_translated, hmlg, ham_k)
    !============================================!
    !                                            !
    !!  Calculate the Hamiltonian in the WF basis
    !                                            !
    !============================================!

    use w90_constants, only: cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_stopwatch

    implicit none

    complex(kind=dp), allocatable, intent(inout) :: ham_k(:, :, :)
!   logical, intent(inout) :: have_translated
!   logical, intent(inout) :: use_translation
    type(ham_logical), intent(inout) :: hmlg

    ! passed variables
    integer, intent(inout), allocatable :: shift_vec(:, :)
    integer, intent(inout) :: nrpts
    integer, intent(inout) :: irvec(:, :)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: timing_level
    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_species
    integer, intent(in) :: atoms_species_num(:)
    integer, intent(in) :: ndimwin(:)

    real(kind=dp), intent(inout) :: wannier_centres_translated(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) :: wannier_centres(:, :)
    real(kind=dp), intent(in) :: atoms_pos_cart(:, :, :)
    real(kind=dp), intent(out) :: translation_centre_frac(3)
    real(kind=dp), intent(in) :: lenconfac
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: eigval(:, :)

    complex(kind=dp), intent(inout) :: ham_r(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)

    logical, intent(in) :: automatic_translation
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: lwindow(:, :)
    logical, intent(in) :: lsitesymmetry  !YN:

    ! local variables
    complex(kind=dp)     :: fac
    real(kind=dp)        :: rdotk
    real(kind=dp)        :: eigval_opt(num_bands, num_kpts)
    real(kind=dp)        :: eigval2(num_wann, num_kpts)
    real(kind=dp)        :: irvec_tmp(3)
    integer              :: loop_kpt, i, j, m, irpt, ideg, ierr, counter
    complex(kind=dp)     :: utmp(num_bands, num_wann) !RS:

    if (timing_level > 1) call io_stopwatch('hamiltonian: get_hr', 1)

    if (hmlg%have_ham_r) then
      if (hmlg%have_translated .eqv. hmlg%use_translation) then
        goto 200
      else
        goto 100
      endif
    end if

    if (hmlg%have_ham_k) go to 100

!~    if (.not. allocated(ham_k)) then
!~       allocate(ham_k(num_wann,num_wann,num_kpts),stat=ierr)
!~       if (ierr/=0) call io_error('Error in allocating ham_k in hamiltonian_get_hr')
!~    end if

    ham_k = cmplx_0
    eigval_opt = 0.0_dp
    eigval2 = 0.0_dp

    if (have_disentangled) then

      ! slim down eigval to contain states within the outer window

      do loop_kpt = 1, num_kpts
        counter = 0
        do j = 1, num_bands
          if (lwindow(j, loop_kpt)) then
            counter = counter + 1
            eigval_opt(counter, loop_kpt) = eigval(j, loop_kpt)
          end if
        end do
      end do

      ! rotate eigval into the optimal subspace
      ! in general eigval would be a matrix at each kpoints
      ! but we choose u_matrix_opt such that the Hamiltonian is
      ! diagonal at each kpoint. (I guess we should check it here)

      if (.not. lsitesymmetry) then                                                                             !YN:
        do loop_kpt = 1, num_kpts
          do j = 1, num_wann
            do m = 1, ndimwin(loop_kpt)
              eigval2(j, loop_kpt) = eigval2(j, loop_kpt) + eigval_opt(m, loop_kpt)* &
                                     real(conjg(u_matrix_opt(m, j, loop_kpt))*u_matrix_opt(m, j, loop_kpt), dp)
            enddo
          enddo
        enddo
      else                                                                                                     !YN:
        ! u_matrix_opt are not the eigenvectors of the Hamiltonian any more                                   !RS:
        ! so we have to calculate ham_k in the following way                                                  !RS:
        do loop_kpt = 1, num_kpts                                                                                !RS:
          utmp(1:ndimwin(loop_kpt), :) = &                                                                     !RS:
            matmul(u_matrix_opt(1:ndimwin(loop_kpt), :, loop_kpt), u_matrix(:, :, loop_kpt))                    !RS:
          do j = 1, num_wann                                                                                    !RS:
            do i = 1, j                                                                                        !RS:
              do m = 1, ndimwin(loop_kpt)                                                                     !RS:
                ham_k(i, j, loop_kpt) = ham_k(i, j, loop_kpt) + eigval_opt(m, loop_kpt)*conjg(utmp(m, i))*utmp(m, j) !RS:
              enddo                                                                                        !RS:
              if (i .lt. j) ham_k(j, i, loop_kpt) = conjg(ham_k(i, j, loop_kpt))                                   !RS:
            enddo                                                                                           !RS:
          enddo                                                                                              !RS:
        enddo                                                                                                 !RS:
      endif                                                                                                    !YN:

    else
      eigval2(1:num_wann, :) = eigval(1:num_wann, :)
    end if

    ! At this point eigval2 contains num_wann values which belong to the wannier subspace.

    ! Rotate Hamiltonian into the basis of smooth bloch states
    !          H(k)=U^{dagger}(k).H_0(k).U(k)
    ! Note: we enforce hermiticity here

    if (.not. lsitesymmetry .or. .not. have_disentangled) then !YN:
      do loop_kpt = 1, num_kpts
        do j = 1, num_wann
          do i = 1, j
            do m = 1, num_wann
              ham_k(i, j, loop_kpt) = ham_k(i, j, loop_kpt) + eigval2(m, loop_kpt)* &
                                      conjg(u_matrix(m, i, loop_kpt))*u_matrix(m, j, loop_kpt)
            enddo
            if (i .lt. j) ham_k(j, i, loop_kpt) = conjg(ham_k(i, j, loop_kpt))
          enddo
        enddo
      enddo
    endif                                                  !YN:

    hmlg%have_ham_k = .true.

100 continue

    ! Fourier transform rotated hamiltonian into WF basis
    ! H_ij(k) --> H_ij(R) = (1/N_kpts) sum_k e^{-ikR} H_ij(k)
!~    if (.not.allocated(ham_r)) then
!~      allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
!~      if (ierr/=0) call io_error('Error in allocating ham_r in hamiltonian_get_hr')
!~    end if

    ham_r = cmplx_0

    if (.not. hmlg%use_translation) then

      do irpt = 1, nrpts
        do loop_kpt = 1, num_kpts
          rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
          fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
          ham_r(:, :, irpt) = ham_r(:, :, irpt) + fac*ham_k(:, :, loop_kpt)
        enddo
      enddo

      hmlg%have_translated = .false.

    else

      allocate (shift_vec(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating shift_vec in hamiltonian_get_hr')
      call internal_translate_centres(num_wann, lenconfac, atoms_species_num, num_species, &
                                      automatic_translation, translation_centre_frac, atoms_pos_cart, num_atoms, &
                                      wannier_centres, recip_lattice, real_lattice, shift_vec, wannier_centres_translated)

      do irpt = 1, nrpts
        do loop_kpt = 1, num_kpts
          do i = 1, num_wann
            do j = 1, num_wann
              ! ham_r(j,i,irpt)
              ! interaction btw j at 0 and i at irvec(:,irpt)
              irvec_tmp(:) = irvec(:, irpt) + shift_vec(:, i) - shift_vec(:, j)
              rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec_tmp(:), dp))
              fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
              ham_r(j, i, irpt) = ham_r(j, i, irpt) + fac*ham_k(j, i, loop_kpt)
            end do
          end do
        enddo
      enddo

      hmlg%have_translated = .true.

    end if

    ! [lp] if required, compute the minimum diistances
!     if (use_ws_distance) then
!         allocate(irdist_ws(3,ndegenx,num_wann,num_wann,nrpts),stat=ierr)
!         if (ierr/=0) call io_error('Error in allocating irdist_ws in hamiltonian_get_hr')
!         allocate(wdist_ndeg(num_wann,num_wann,nrpts),stat=ierr)
!         if (ierr/=0) call io_error('Error in allocating wcenter_ndeg in hamiltonian_get_hr')
    !
!         call ws_translate_dist(nrpts, irvec)
!     endif

    hmlg%have_ham_r = .true.

200 continue

    if (allocated(shift_vec)) then
      deallocate (shift_vec, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating shift_vec in hamiltonian_get_hr')
    end if

    if (timing_level > 1) call io_stopwatch('hamiltonian: get_hr', 2)

    return

  contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !====================================================!
    subroutine internal_translate_centres(num_wann, lenconfac, atoms_species_num, &
                                          num_species, automatic_translation, translation_centre_frac, &
                                          atoms_pos_cart, num_atoms, wannier_centres, recip_lattice, real_lattice, &
                                          shift_vec, wannier_centres_translated)
      !! Translate the centres of the WF into the home cell
      !====================================================!

      use w90_io, only: stdout, io_error
      use w90_utility, only: utility_cart_to_frac, utility_frac_to_cart

      implicit none

      ! passed variables
      integer, intent(inout) :: shift_vec(:, :)
      integer, intent(in) :: num_wann
      integer, intent(in) :: atoms_species_num(:)
      integer, intent(in) :: num_species
      integer, intent(in) :: num_atoms

      real(kind=dp), intent(inout) :: wannier_centres_translated(:, :)
      real(kind=dp), intent(in) :: real_lattice(3, 3)
      real(kind=dp), intent(in) :: lenconfac
      real(kind=dp), intent(out) :: translation_centre_frac(3)
      real(kind=dp), intent(in) :: atoms_pos_cart(:, :, :)
      real(kind=dp), intent(in) :: wannier_centres(:, :)
      real(kind=dp), intent(in) :: recip_lattice(3, 3)

      logical, intent(in) :: automatic_translation

      ! local variables
      integer :: iw, ierr, nat, nsp, ind

      real(kind=dp), allocatable :: r_home(:, :), r_frac(:, :)
      real(kind=dp) :: c_pos_cart(3), c_pos_frac(3)
      real(kind=dp) :: r_frac_min(3)

!~      if (.not.allocated(wannier_centres_translated)) then
!~         allocate(wannier_centres_translated(3,num_wann),stat=ierr)
!~         if (ierr/=0) call io_error('Error in allocating wannier_centres_translated &
!~              &in internal_translate_wannier_centres')
!~      end if

      allocate (r_home(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating r_home in internal_translate_centres')
      allocate (r_frac(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating r_frac in internal_translate_centres')
      r_home = 0.0_dp; r_frac = 0.0_dp

      if (automatic_translation) then
        ! Calculate centre of atomic positions
        c_pos_cart = 0.0_dp; c_pos_frac = 0.0_dp
        do nsp = 1, num_species
          do nat = 1, atoms_species_num(nsp)
            c_pos_cart(:) = c_pos_cart(:) + atoms_pos_cart(:, nat, nsp)
          enddo
        enddo
        c_pos_cart = c_pos_cart/num_atoms
        ! Cartesian --> fractional
        call utility_cart_to_frac(c_pos_cart, translation_centre_frac, recip_lattice)
      end if
      ! Wannier function centres will be in [c_pos_frac-0.5,c_pos_frac+0.5]
      r_frac_min(:) = translation_centre_frac(:) - 0.5_dp

      ! Cartesian --> fractional
      do iw = 1, num_wann
        call utility_cart_to_frac(wannier_centres(:, iw), r_frac(:, iw), recip_lattice)
        ! Rationalise r_frac - r_frac_min to interval [0,1]
        !  by applying shift of -floor(r_frac - r_frac_min)
        shift_vec(:, iw) = -floor(r_frac(:, iw) - r_frac_min(:))
        r_frac(:, iw) = r_frac(:, iw) + real(shift_vec(:, iw), dp)
        ! Fractional --> Cartesian
        call utility_frac_to_cart(r_frac(:, iw), r_home(:, iw), real_lattice)
      end do

      ! NEVER overwrite wannier_centres
      !wannier_centres = r_home

      if (on_root) then
        write (stdout, '(1x,a)') 'Translated centres'
        write (stdout, '(4x,a,3f10.6)') 'translation centre in fractional coordinate:', translation_centre_frac(:)
        do iw = 1, num_wann
          write (stdout, 888) iw, (r_home(ind, iw)*lenconfac, ind=1, 3)
        end do
        write (stdout, '(1x,a78)') repeat('-', 78)
        write (stdout, *)
      endif
      wannier_centres_translated = r_home

      deallocate (r_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating r_frac in internal_translate_centres')
      deallocate (r_home, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating r_home in internal_translate_centres')

      return

888   format(2x, 'WF centre ', i5, 2x, '(', f10.6, ',', f10.6, ',', f10.6, ' )')

    end subroutine internal_translate_centres

  end subroutine hamiltonian_get_hr

  !============================================!
  subroutine hamiltonian_write_hr(num_wann, timing_level, ham_r, irvec, ndegen, nrpts, hmlg)
    !============================================!
    !!  Write the Hamiltonian in the WF basis
    !============================================!

    use w90_io, only: io_error, io_stopwatch, io_file_unit, &
      seedname, io_date

    integer            :: i, j, irpt, file_unit
    character(len=33) :: header
    character(len=9)  :: cdate, ctime

    integer, intent(inout) :: nrpts
    integer, intent(in) :: ndegen(:)
    integer, intent(inout) :: irvec(:, :)
    complex(kind=dp), intent(in) :: ham_r(:, :, :)
    type(ham_logical), intent(inout) :: hmlg

!   from w90_parameters
    integer, intent(in) :: num_wann
    integer, intent(in) :: timing_level
!   end w90_parameters

    if (hmlg%hr_written) return

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_hr', 1)

    ! write the  whole matrix with all the indices

    file_unit = io_file_unit()
    open (file_unit, file=trim(seedname)//'_hr.dat', form='formatted', &
          status='unknown', err=101)

    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime

    write (file_unit, *) header ! Date and time
    write (file_unit, *) num_wann
    write (file_unit, *) nrpts
    write (file_unit, '(15I5)') (ndegen(i), i=1, nrpts)
    do irpt = 1, nrpts
      do i = 1, num_wann
        do j = 1, num_wann
          write (file_unit, '(5I5,2F12.6)') irvec(:, irpt), j, i, &
            ham_r(j, i, irpt)
        end do
      end do
    end do

    close (file_unit)

    hmlg%hr_written = .true.

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_hr', 2)

    return

101 call io_error('Error: hamiltonian_write_hr: problem opening file '//trim(seedname)//'_hr.dat')

  end subroutine hamiltonian_write_hr

  !================================================================================!
  subroutine hamiltonian_wigner_seitz(rpt_origin, nrpts, ndegen, irvec, mp_grid, &
                                      real_metric, ws_search_size, ws_distance_tol, timing_level, iprint, &
                                      count_pts)
    !================================================================================!
    !! Calculates a grid of points that fall inside of (and eventually on the
    !! surface of) the Wigner-Seitz supercell centered on the origin of the B
    !! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3
    !================================================================================!

    use w90_constants, only: eps7, eps8
    use w90_io, only: io_error, io_stopwatch, stdout

    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points

    implicit none

    ! passed variables
    integer, intent(inout) :: nrpts
    integer, intent(inout), allocatable :: ndegen(:)
    integer, intent(inout), allocatable :: irvec(:, :)
    integer, intent(inout) :: rpt_origin
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: ws_search_size(3)
    integer, intent(in) :: timing_level
    integer, intent(in) :: iprint

    logical, intent(in) :: count_pts

    real(kind=dp), intent(in) :: real_metric(3, 3)
    real(kind=dp), intent(in) :: ws_distance_tol

    ! local variables
    integer       :: ndiff(3)
    real(kind=dp) :: tot, dist_min
    real(kind=dp), allocatable :: dist(:)
    integer       :: n1, n2, n3, i1, i2, i3, icnt, i, j, ierr, dist_dim

    if (timing_level > 1) call io_stopwatch('hamiltonian: wigner_seitz', 1)

    dist_dim = 1
    do i = 1, 3
      dist_dim = dist_dim*((ws_search_size(i) + 1)*2 + 1)
    end do
    allocate (dist(dist_dim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating dist in hamiltonian_wigner_seitz')

    ! The Wannier functions live in a supercell of the real space unit cell
    ! this supercell is mp_grid unit cells long in each direction
    !
    ! We loop over grid points r on a unit cell that is (2*ws_search_size+1)**3 times
    ! larger than this primitive supercell.
    !
    ! One of these points is in the W-S cell if it is closer to R=0 than any of the
    ! other points, R (where R are the translation vectors of the supercell)

    ! In the end nrpts contains the total number of grid
    ! points that have been found in the Wigner-Seitz cell

    nrpts = 0
    ! Loop over the lattice vectors of the primitive cell
    ! that live in a supercell which is (2*ws_search_size+1)**2
    ! larger than the Born-von Karman supercell.
    ! We need to find which among these live in the Wigner-Seitz cell
    do n1 = -ws_search_size(1)*mp_grid(1), ws_search_size(1)*mp_grid(1)
      do n2 = -ws_search_size(2)*mp_grid(2), ws_search_size(2)*mp_grid(2)
        do n3 = -ws_search_size(3)*mp_grid(3), ws_search_size(3)*mp_grid(3)
          ! Loop over the lattice vectors R of the Born-von Karman supercell
          ! that contains all the points of the previous loop.
          ! There are (2*(ws_search_size+1)+1)**3 points R. R=0 corresponds to
          ! i1=i2=i3=0, or icnt=((2*(ws_search_size+1)+1)**3 + 1)/2
          icnt = 0
          do i1 = -ws_search_size(1) - 1, ws_search_size(1) + 1
            do i2 = -ws_search_size(2) - 1, ws_search_size(2) + 1
              do i3 = -ws_search_size(3) - 1, ws_search_size(3) + 1
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1*mp_grid(1)
                ndiff(2) = n2 - i2*mp_grid(2)
                ndiff(3) = n3 - i3*mp_grid(3)
                dist(icnt) = 0.0_dp
                do i = 1, 3
                  do j = 1, 3
                    dist(icnt) = dist(icnt) + real(ndiff(i), dp)*real_metric(i, j) &
                                 *real(ndiff(j), dp)
                  enddo
                enddo
              enddo
            enddo
          enddo
          ! AAM: On first pass, we reference unallocated variables (ndegen,irvec)
          dist_min = minval(dist)
          if (abs(dist((dist_dim + 1)/2) - dist_min) .lt. ws_distance_tol**2) then
            nrpts = nrpts + 1
            if (.not. count_pts) then
              ndegen(nrpts) = 0
              do i = 1, dist_dim
                if (abs(dist(i) - dist_min) .lt. ws_distance_tol**2) ndegen(nrpts) = ndegen(nrpts) + 1
              end do
              irvec(1, nrpts) = n1
              irvec(2, nrpts) = n2
              irvec(3, nrpts) = n3
              !
              ! Record index of r=0
              if (n1 == 0 .and. n2 == 0 .and. n3 == 0) rpt_origin = nrpts
            endif
          end if

          !n3
        enddo
        !n2
      enddo
      !n1
    enddo
    !
    deallocate (dist, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating dist hamiltonian_wigner_seitz')
    if (count_pts) then
      if (timing_level > 1) call io_stopwatch('hamiltonian: wigner_seitz', 2)
      return
    end if

    ! Check the "sum rule"
    tot = 0.0_dp
    do i = 1, nrpts
      tot = tot + 1.0_dp/real(ndegen(i), dp)
    enddo

    if (iprint >= 3 .and. on_root) then
      write (stdout, '(1x,i4,a,/)') nrpts, ' lattice points in Wigner-Seitz supercell:'
      do i = 1, nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', irvec(1, i), irvec(2, i), &
          irvec(3, i), '  degeneracy: ', ndegen(i)
      enddo
      write (stdout, '(1x,a,f12.3)') ' tot = ', tot
      write (stdout, '(1x,a,i12)') ' mp_grid product = ', mp_grid(1)*mp_grid(2)*mp_grid(3)
    endif
    if (abs(tot - real(mp_grid(1)*mp_grid(2)*mp_grid(3), dp)) > eps8) then
      call io_error('ERROR in hamiltonian_wigner_seitz: error in finding Wigner-Seitz points')
    endif

    if (timing_level > 1) call io_stopwatch('hamiltonian: wigner_seitz', 2)

    return

  end subroutine hamiltonian_wigner_seitz

  !============================================!
  subroutine hamiltonian_write_rmn(m_matrix, wb, bk, num_wann, num_kpts, kpt_latt, &
                                   nntot, irvec, nrpts)
    !! Write out the matrix elements of r
    !============================================!

    use w90_constants, only: twopi, cmplx_i
    use w90_io, only: io_error, io_file_unit, seedname, io_date

    implicit none

    ! passed variables
    integer, intent(inout) :: nrpts
    integer, intent(inout) :: irvec(:, :)
    integer, intent(in) :: num_wann
    integer, intent(in) :: num_kpts
    integer, intent(in) :: nntot
    real(kind=dp), intent(in) :: wb(:)
    real(kind=dp), intent(in) :: bk(:, :, :)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)

    ! local variables
    complex(kind=dp) :: fac
    real(kind=dp) :: rdotk
    complex(kind=dp) :: position(3)
    character(len=33) :: header
    character(len=9)  :: cdate, ctime
    integer :: loop_rpt, m, n, nkp, ind, nn, file_unit

    file_unit = io_file_unit()
    open (file_unit, file=trim(seedname)//'_r.dat', form='formatted', status='unknown', err=101)
    call io_date(cdate, ctime)

    header = 'written on '//cdate//' at '//ctime
    write (file_unit, *) header ! Date and time
    write (file_unit, *) num_wann
    write (file_unit, *) nrpts

    do loop_rpt = 1, nrpts
      do m = 1, num_wann
        do n = 1, num_wann
          position(:) = 0._dp
          do nkp = 1, num_kpts
            rdotk = twopi*dot_product(kpt_latt(:, nkp), real(irvec(:, loop_rpt), dp))
            fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
            do ind = 1, 3
              do nn = 1, nntot
                if (m .eq. n) then
                  ! For loop_rpt==rpt_origin, this reduces to
                  ! Eq.(32) of Marzari and Vanderbilt PRB 56,
                  ! 12847 (1997). Otherwise, is is Eq.(44)
                  ! Wang, Yates, Souza and Vanderbilt PRB 74,
                  ! 195118 (2006), modified according to
                  ! Eqs.(27,29) of Marzari and Vanderbilt
                  position(ind) = position(ind) - &
                                  wb(nn)*bk(ind, nn, nkp)*aimag(log(m_matrix(n, m, nn, nkp)))*fac
                else
                  ! Eq.(44) Wang, Yates, Souza and Vanderbilt PRB 74, 195118 (2006)
                  position(ind) = position(ind) + &
                                  cmplx_i*wb(nn)*bk(ind, nn, nkp)*m_matrix(n, m, nn, nkp)*fac
                endif
              end do
            end do
          end do
          write (file_unit, '(5I5,6F12.6)') irvec(:, loop_rpt), n, m, position(:)
        end do
      end do
    end do

    close (file_unit)

    return

101 call io_error('Error: hamiltonian_write_rmn: problem opening file '//trim(seedname)//'_r')

  end subroutine hamiltonian_write_rmn

  !============================================!
  subroutine hamiltonian_write_tb(real_lattice, num_wann, wb, bk, m_matrix, &
                                  num_kpts, kpt_latt, nntot, timing_level, ham_r, irvec, ndegen, nrpts, &
                                  hmlg)
    !============================================!
    !! Write in a single file all the information
    !! that is needed to set up a Wannier-based
    !! tight-binding model:
    !! * lattice vectors
    !! * <0n|H|Rn>
    !! * <0n|r|Rn>
    !============================================!

    use w90_io, only: io_error, io_stopwatch, io_file_unit, &
      seedname, io_date

    use w90_constants, only: twopi, cmplx_i

    integer            :: i, j, irpt, ik, nn, idir, file_unit
    character(len=33)  :: header
    character(len=9)   :: cdate, ctime
    complex(kind=dp)   :: fac, pos_r(3)
    real(kind=dp)      :: rdotk

    integer, intent(inout) :: nrpts
    integer, intent(in) :: ndegen(:)
    integer, intent(inout) :: irvec(:, :)
    complex(kind=dp), intent(in) :: ham_r(:, :, :)
    type(ham_logical), intent(inout) :: hmlg

!   from w90_parameters
    integer, intent(in) :: num_wann
    integer, intent(in) :: num_kpts
    integer, intent(in) :: nntot
    integer, intent(in) :: timing_level
    real(kind=dp), intent(in) :: wb(:)
    real(kind=dp), intent(in) :: bk(:, :, :)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
!   end w90_parameters

    if (hmlg%tb_written) return

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_tb', 1)

    file_unit = io_file_unit()
    open (file_unit, file=trim(seedname)//'_tb.dat', form='formatted', &
          status='unknown', err=101)

    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime

    write (file_unit, *) header ! Date and time
    !
    ! lattice vectors
    !
    write (file_unit, *) real_lattice(1, :) !a_1
    write (file_unit, *) real_lattice(2, :) !a_2
    write (file_unit, *) real_lattice(3, :) !a_3
    !
    write (file_unit, *) num_wann
    write (file_unit, *) nrpts
    write (file_unit, '(15I5)') (ndegen(i), i=1, nrpts)
    !
    ! <0n|H|Rm>
    !
    do irpt = 1, nrpts
      write (file_unit, '(/,3I5)') irvec(:, irpt)
      do i = 1, num_wann
        do j = 1, num_wann
          write (file_unit, '(2I5,3x,2(E15.8,1x))') j, i, ham_r(j, i, irpt)
        end do
      end do
    end do
    !
    ! <0n|r|Rm>
    !
    do irpt = 1, nrpts
      write (file_unit, '(/,3I5)') irvec(:, irpt)
      do i = 1, num_wann
        do j = 1, num_wann
          pos_r(:) = 0._dp
          do ik = 1, num_kpts
            rdotk = twopi*dot_product(kpt_latt(:, ik), real(irvec(:, irpt), dp))
            fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
            do idir = 1, 3
              do nn = 1, nntot
                if (i == j) then
                  ! For irpt==rpt_origin, this reduces to
                  ! Eq.(32) of Marzari and Vanderbilt PRB 56,
                  ! 12847 (1997). Otherwise, is is Eq.(44)
                  ! Wang, Yates, Souza and Vanderbilt PRB 74,
                  ! 195118 (2006), modified according to
                  ! Eqs.(27,29) of Marzari and Vanderbilt
                  pos_r(idir) = pos_r(idir) - &
                                wb(nn)*bk(idir, nn, ik)*aimag(log(m_matrix(i, i, nn, ik)))*fac
                else
                  ! Eq.(44) Wang, Yates, Souza and Vanderbilt PRB 74, 195118 (2006)
                  pos_r(idir) = pos_r(idir) + &
                                cmplx_i*wb(nn)*bk(idir, nn, ik)*m_matrix(j, i, nn, ik)*fac
                endif
              end do
            end do
          end do
          write (file_unit, '(2I5,3x,6(E15.8,1x))') j, i, pos_r(:)
        end do
      end do
    end do
    close (file_unit)

    hmlg%tb_written = .true.

    if (timing_level > 1) call io_stopwatch('hamiltonian: write_tb', 2)

    return

101 call io_error('Error: hamiltonian_write_tb: problem opening file ' &
                  //trim(seedname)//'_tb.dat')

  end subroutine hamiltonian_write_tb

end module w90_hamiltonian
