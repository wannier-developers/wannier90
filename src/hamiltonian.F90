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
!  w90_hamiltonian: Hamiltonian in Wannier basis             !
!                                                            !
!------------------------------------------------------------!

module w90_hamiltonian

  !! Module to obtain the Hamiltonian in a Wannier basis
  !! This is a simplified routine, more sophisticated properties
  !! are found in postw90 (e.g. w90_get_oper)

  use w90_constants, only: dp
  use w90_types
  use w90_error

  implicit none

  public :: hamiltonian_dealloc
  public :: hamiltonian_get_hr
  public :: hamiltonian_setup
  public :: hamiltonian_write_hr
  public :: hamiltonian_write_tb

contains

  !================================================!

  subroutine hamiltonian_setup(ham_logical, print_output, ws_region, w90_calculation, ham_k, &
                               ham_r, real_lattice, wannier_centres_translated, irvec, mp_grid, &
                               ndegen, num_kpts, num_wann, nrpts, rpt_origin, bands_plot_mode, &
                               stdout, timer, error, transport_mode, comm)
    !================================================!
    !
    !! Allocate arrays and setup data
    !
    !================================================!

    use w90_constants, only: cmplx_0
    use w90_types, only: print_output_type, ws_region_type, timer_list_type
    use w90_wannier90_types, only: w90_calculation_type, ham_logical_type

    implicit none

    ! arguments
    type(ham_logical_type), intent(inout) :: ham_logical
    type(print_output_type), intent(in) :: print_output
    type(w90_calculation_type), intent(in) :: w90_calculation
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(ws_region_type), intent(in) :: ws_region
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: mp_grid(3)
    integer, intent(inout), allocatable :: irvec(:, :)
    integer, intent(inout), allocatable :: ndegen(:)
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(inout) :: nrpts
    integer, intent(inout) :: rpt_origin
    integer, intent(in) :: stdout

    real(kind=dp), intent(in)                 :: real_lattice(3, 3)
    real(kind=dp), intent(inout), allocatable :: wannier_centres_translated(:, :)

    complex(kind=dp), intent(inout), allocatable :: ham_k(:, :, :)
    complex(kind=dp), intent(inout), allocatable :: ham_r(:, :, :)

    character(len=*), intent(in) :: bands_plot_mode
    character(len=20), intent(in)  :: transport_mode

    ! local variables
    integer :: ierr

    if (ham_logical%ham_have_setup) return
    !
    ! Determine whether to use translation
    !
    if (w90_calculation%bands_plot .and. (index(bands_plot_mode, 'cut') .ne. 0)) &
      ham_logical%use_translation = .true.
    if (w90_calculation%transport .and. (index(transport_mode, 'bulk') .ne. 0)) &
      ham_logical%use_translation = .true.
    if (w90_calculation%transport .and. (index(transport_mode, 'lcr') .ne. 0)) &
      ham_logical%use_translation = .true.
    !
    ! Set up Wigner-Seitz vectors
    !
    call hamiltonian_wigner_seitz(ws_region, print_output, real_lattice, irvec, mp_grid, ndegen, &
                                  nrpts, rpt_origin, stdout, timer, error, .true., comm)
    if (allocated(error)) return

    allocate (irvec(3, nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating irvec in hamiltonian_setup', comm)
      return
    endif
    irvec = 0

    allocate (ndegen(nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ndegen in hamiltonian_setup', comm)
      return
    endif
    ndegen = 0

    allocate (ham_r(num_wann, num_wann, nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ham_r in hamiltonian_setup', comm)
      return
    endif
    ham_r = cmplx_0

    allocate (ham_k(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ham_k in hamiltonian_setup', comm)
      return
    endif
    ham_k = cmplx_0
    !
    ! Set up the wigner_seitz vectors
    !
    call hamiltonian_wigner_seitz(ws_region, print_output, real_lattice, irvec, mp_grid, ndegen, &
                                  nrpts, rpt_origin, stdout, timer, error, .false., comm)
    if (allocated(error)) return

    allocate (wannier_centres_translated(3, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wannier_centres_translated in hamiltonian_setup', comm)
      return
    endif

    wannier_centres_translated = 0.0_dp
    ham_logical%ham_have_setup = .true.

    return
  end subroutine hamiltonian_setup

  !================================================!
  subroutine hamiltonian_dealloc(ham_logical, ham_k, ham_r, wannier_centres_translated, irvec, &
                                 ndegen, error, comm)
    !================================================!
    !
    !! Deallocate module data
    !
    !================================================!

    use w90_wannier90_types, only: ham_logical_type

    implicit none

    ! arguments
    type(ham_logical_type), intent(inout) :: ham_logical
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(inout), allocatable :: ndegen(:)
    integer, intent(inout), allocatable :: irvec(:, :)

    real(kind=dp), intent(inout), allocatable :: wannier_centres_translated(:, :)

    complex(kind=dp), intent(inout), allocatable :: ham_r(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: ham_k(:, :, :)

    ! local variables
    integer :: ierr

    if (allocated(ham_r)) then
      deallocate (ham_r, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ham_r in hamiltonian_dealloc', comm)
        return
      endif
    end if
    if (allocated(ham_k)) then
      deallocate (ham_k, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ham_k in hamiltonian_dealloc', comm)
        return
      endif
    end if
    if (allocated(irvec)) then
      deallocate (irvec, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating irvec in hamiltonian_dealloc', comm)
        return
      endif
    end if
    if (allocated(ndegen)) then
      deallocate (ndegen, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ndegen in hamiltonian_dealloc', comm)
        return
      endif
    end if
    if (allocated(wannier_centres_translated)) then
      deallocate (wannier_centres_translated, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wannier_centres_translated in w90_readwrite_dealloc', comm)
        return
      endif
    end if

    ham_logical%ham_have_setup = .false.
    ham_logical%have_translated = .false.
    ham_logical%use_translation = .false.
    ham_logical%have_ham_r = .false.
    ham_logical%have_ham_k = .false.
    !================================================!
  end subroutine hamiltonian_dealloc

  !================================================!
  subroutine hamiltonian_get_hr(atom_data, dis_manifold, ham_logical, real_space_ham, &
                                print_output, ham_k, ham_r, u_matrix, u_matrix_opt, eigval, &
                                kpt_latt, real_lattice, wannier_centres, &
                                wannier_centres_translated, irvec, shift_vec, nrpts, num_bands, &
                                num_kpts, num_wann, have_disentangled, stdout, timer, error, &
                                lsitesymmetry, comm)
    !================================================!
    !
    !!  Calculate the Hamiltonian in the WF basis
    !
    !================================================!

    use w90_constants, only: cmplx_0, cmplx_i, twopi
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: atom_data_type, dis_manifold_type, print_output_type, timer_list_type
    use w90_wannier90_types, only: real_space_ham_type, ham_logical_type

    implicit none

    ! arguments
    type(ham_logical_type), intent(inout)    :: ham_logical
    type(atom_data_type), intent(in)         :: atom_data
    type(real_space_ham_type), intent(inout) :: real_space_ham
    type(print_output_type), intent(in)      :: print_output
    type(dis_manifold_type), intent(in)      :: dis_manifold
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in)           :: comm
    type(timer_list_type), intent(inout)     :: timer

    integer, intent(inout), allocatable :: shift_vec(:, :)
    integer, intent(inout)              :: irvec(:, :)
    integer, intent(inout)              :: nrpts
    integer, intent(in)                 :: num_bands
    integer, intent(in)                 :: num_kpts
    integer, intent(in)                 :: num_wann
    integer, intent(in)                 :: stdout

    real(kind=dp), intent(inout) :: wannier_centres_translated(:, :)
    real(kind=dp), intent(in)    :: real_lattice(3, 3)
    real(kind=dp), intent(in)    :: wannier_centres(:, :)
    real(kind=dp), intent(in)    :: kpt_latt(:, :)
    real(kind=dp), intent(in)    :: eigval(:, :)

    complex(kind=dp), intent(inout)              :: ham_r(:, :, :)
    complex(kind=dp), intent(in)                 :: u_matrix(:, :, :)
    complex(kind=dp), intent(in)                 :: u_matrix_opt(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: ham_k(:, :, :)

    logical, intent(in) :: lsitesymmetry  !YN:
    logical, intent(in) :: have_disentangled

    ! local variables
    integer          :: loop_kpt, i, j, m, irpt, ierr, counter
    real(kind=dp)    :: rdotk
    real(kind=dp), allocatable    :: eigval_opt(:, :) !(num_bands, num_kpts)
    real(kind=dp), allocatable    :: eigval2(:, :)    !(num_wann, num_kpts)
    real(kind=dp)    :: irvec_tmp(3)
    complex(kind=dp), allocatable :: utmp(:, :)       !(num_bands, num_wann)
    complex(kind=dp) :: fac

    if (print_output%timing_level > 1) call io_stopwatch_start('hamiltonian: get_hr', timer)

    if (ham_logical%have_ham_r) then
      if (ham_logical%have_translated .eqv. ham_logical%use_translation) then
        goto 200
      else
        goto 100
      endif
    end if

    if (ham_logical%have_ham_k) go to 100

    ham_k = cmplx_0

    allocate (eigval2(num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating eigval2 in hamiltonian_get_hr', comm)
      return
    endif

    eigval2 = 0.0_dp

    if (have_disentangled) then

      ! start allocation of eigval_opt, utmp; used only if have_disentangled.
      allocate (eigval_opt(num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating eigval_opt in hamiltonian_get_hr', comm)
        return
      endif

      allocate (utmp(num_bands, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating utmp in hamiltonian_get_hr', comm)
        return
      endif

      eigval_opt = 0.0_dp
      ! end allocation of eigval_opt, utmp

      ! slim down eigval to contain states within the outer window

      do loop_kpt = 1, num_kpts
        counter = 0
        do j = 1, num_bands
          if (dis_manifold%lwindow(j, loop_kpt)) then
            counter = counter + 1
            eigval_opt(counter, loop_kpt) = eigval(j, loop_kpt)
          end if
        end do
      end do

      ! rotate eigval into the optimal subspace
      ! in general eigval would be a matrix at each kpoints
      ! but we choose u_matrix_opt such that the Hamiltonian is
      ! diagonal at each kpoint. (I guess we should check it here)

      if (.not. lsitesymmetry) then
        do loop_kpt = 1, num_kpts
          do j = 1, num_wann
            do m = 1, dis_manifold%ndimwin(loop_kpt)
              eigval2(j, loop_kpt) = eigval2(j, loop_kpt) + eigval_opt(m, loop_kpt)* &
                                     real(conjg(u_matrix_opt(m, j, loop_kpt))* &
                                          u_matrix_opt(m, j, loop_kpt), dp)
            enddo
          enddo
        enddo
      else
        ! u_matrix_opt are not the eigenvectors of the Hamiltonian any more
        ! so we have to calculate ham_k in the following way
        do loop_kpt = 1, num_kpts
          utmp(1:dis_manifold%ndimwin(loop_kpt), :) = &
            matmul(u_matrix_opt(1:dis_manifold%ndimwin(loop_kpt), :, loop_kpt), &
                   u_matrix(:, :, loop_kpt))
          do j = 1, num_wann
            do i = 1, j
              do m = 1, dis_manifold%ndimwin(loop_kpt)
                ham_k(i, j, loop_kpt) = ham_k(i, j, loop_kpt) + eigval_opt(m, loop_kpt)* &
                                        conjg(utmp(m, i))*utmp(m, j)
              enddo
              if (i .lt. j) ham_k(j, i, loop_kpt) = conjg(ham_k(i, j, loop_kpt))
            enddo
          enddo
        enddo
      endif

    else
      eigval2(1:num_wann, :) = eigval(1:num_wann, :)
    end if

    ! At this point eigval2 contains num_wann values which belong to the wannier subspace.

    ! Rotate Hamiltonian into the basis of smooth bloch states
    !          H(k)=U^{dagger}(k).H_0(k).U(k)
    ! Note: we enforce hermiticity here

    if (.not. lsitesymmetry .or. .not. have_disentangled) then
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
    endif

    ham_logical%have_ham_k = .true.

100 continue

    ! Fourier transform rotated hamiltonian into WF basis
    ! H_ij(k) --> H_ij(R) = (1/N_kpts) sum_k e^{-ikR} H_ij(k)
!~    if (.not.allocated(ham_r)) then
!~      allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
!~      if (ierr/=0) call io_error('Error in allocating ham_r in hamiltonian_get_hr')
!~    end if

    ham_r = cmplx_0

    if (.not. ham_logical%use_translation) then

      do irpt = 1, nrpts
        do loop_kpt = 1, num_kpts
          rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
          fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
          ham_r(:, :, irpt) = ham_r(:, :, irpt) + fac*ham_k(:, :, loop_kpt)
        enddo
      enddo

      ham_logical%have_translated = .false.

    else

      allocate (shift_vec(3, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating shift_vec in hamiltonian_get_hr', comm)
        return
      endif
      call internal_translate_centres(atom_data, real_space_ham, real_lattice, wannier_centres, &
                                      wannier_centres_translated, shift_vec, print_output%iprint, &
                                      num_wann, error)
      if (allocated(error)) return

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

      ham_logical%have_translated = .true.

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

    ham_logical%have_ham_r = .true.

200 continue

    if (allocated(shift_vec)) then
      deallocate (shift_vec, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating shift_vec in hamiltonian_get_hr', comm)
        return
      endif
    end if

    if (allocated(eigval2)) then
      deallocate (eigval2, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating eigval2 in hamiltonian_get_hr', comm)
        return
      endif
    end if

    if (allocated(eigval_opt)) then
      deallocate (eigval_opt, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating eigval_opt in hamiltonian_get_hr', comm)
        return
      endif
    end if

    if (allocated(utmp)) then
      deallocate (utmp, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating utmp in hamiltonian_get_hr', comm)
        return
      endif
    end if

    if (print_output%timing_level > 1) call io_stopwatch_stop('hamiltonian: get_hr', timer)

    return

  contains

    !================================================!
    subroutine internal_translate_centres(atom_data, real_space_ham, real_lattice, &
                                          wannier_centres, wannier_centres_translated, shift_vec, &
                                          iprint, num_wann, error)
      !================================================!
      !
      !! Translate the centres of the WF into the home cell
      !
      !================================================!

      use w90_utility, only: utility_cart_to_frac, utility_frac_to_cart, utility_inverse_mat
      use w90_types, only: atom_data_type
      use w90_wannier90_types, only: real_space_ham_type

      implicit none

      ! arguments
      type(atom_data_type), intent(in) :: atom_data
      type(real_space_ham_type), intent(inout) :: real_space_ham
      type(w90_error_type), allocatable, intent(out) :: error

      integer, intent(inout) :: shift_vec(:, :)
      integer, intent(in)    :: iprint
      integer, intent(in)    :: num_wann

      real(kind=dp), intent(inout) :: wannier_centres_translated(:, :)
      real(kind=dp), intent(in)    :: real_lattice(3, 3)
      real(kind=dp), intent(in)    :: wannier_centres(:, :)

      ! local variables
      integer :: iw, ierr, nat, nsp, ind
      real(kind=dp)              :: inv_lattice(3, 3)
      real(kind=dp), allocatable :: r_home(:, :), r_frac(:, :)
      real(kind=dp)              :: c_pos_cart(3), c_pos_frac(3)
      real(kind=dp)              :: r_frac_min(3)

!~      if (.not.allocated(wannier_centres_translated)) then
!~         allocate(wannier_centres_translated(3,num_wann),stat=ierr)
!~         if (ierr/=0) call io_error('Error in allocating wannier_centres_translated &
!~              &in internal_translate_wannier_centres')
!~      end if

      allocate (r_home(3, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating r_home in internal_translate_centres', comm)
        return
      endif
      allocate (r_frac(3, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating r_frac in internal_translate_centres', comm)
        return
      endif
      r_home = 0.0_dp; r_frac = 0.0_dp

      call utility_inverse_mat(real_lattice, inv_lattice)
      if (real_space_ham%automatic_translation) then
        ! Calculate centre of atomic positions
        c_pos_cart = 0.0_dp; c_pos_frac = 0.0_dp
        do nsp = 1, atom_data%num_species
          do nat = 1, atom_data%species_num(nsp)
            c_pos_cart(:) = c_pos_cart(:) + atom_data%pos_cart(:, nat, nsp)
          enddo
        enddo
        c_pos_cart = c_pos_cart/atom_data%num_atoms
        ! Cartesian --> fractional
        call utility_cart_to_frac(c_pos_cart, real_space_ham%translation_centre_frac, inv_lattice)
      end if
      ! Wannier function centres will be in [c_pos_frac-0.5,c_pos_frac+0.5]
      r_frac_min(:) = real_space_ham%translation_centre_frac(:) - 0.5_dp

      ! Cartesian --> fractional
      do iw = 1, num_wann
        call utility_cart_to_frac(wannier_centres(:, iw), r_frac(:, iw), inv_lattice)
        ! Rationalise r_frac - r_frac_min to interval [0,1]
        !  by applying shift of -floor(r_frac - r_frac_min)
        shift_vec(:, iw) = -floor(r_frac(:, iw) - r_frac_min(:))
        r_frac(:, iw) = r_frac(:, iw) + real(shift_vec(:, iw), dp)
        ! Fractional --> Cartesian
        call utility_frac_to_cart(r_frac(:, iw), r_home(:, iw), real_lattice)
      end do

      ! NEVER overwrite wannier_centres
      !wannier_centres = r_home

      if (iprint > 0) then
        write (stdout, '(1x,a)') 'Translated centres'
        write (stdout, '(4x,a,3f10.6)') 'translation centre in fractional coordinate:', &
          real_space_ham%translation_centre_frac(:)
        do iw = 1, num_wann
          write (stdout, 888) iw, (r_home(ind, iw)*print_output%lenconfac, ind=1, 3)
        end do
        write (stdout, '(1x,a78)') repeat('-', 78)
        write (stdout, *)
      endif
      wannier_centres_translated = r_home

      deallocate (r_frac, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating r_frac in internal_translate_centres', comm)
        return
      endif
      deallocate (r_home, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating r_home in internal_translate_centres', comm)
        return
      endif

      return

888   format(2x, 'WF centre ', i5, 2x, '(', f10.6, ',', f10.6, ',', f10.6, ' )')

    end subroutine internal_translate_centres

  end subroutine hamiltonian_get_hr

  !================================================!
  subroutine hamiltonian_write_hr(ham_r, irvec, ndegen, nrpts, num_wann, timing_level, seedname, &
                                  timer, error, comm)
    !================================================!
    !
    !!  Write the Hamiltonian in the WF basis
    !
    !================================================!

    use w90_io, only: io_stopwatch_start, io_stopwatch_stop, io_date
    use w90_types, only: timer_list_type
    use w90_comms, only: w90_comm_type

    ! arguments
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: irvec(:, :)
    integer, intent(in) :: ndegen(:)
    integer, intent(in) :: nrpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: timing_level

    complex(kind=dp), intent(in) :: ham_r(:, :, :)

    character(len=50), intent(in) :: seedname

    ! local variables
    integer :: i, j, irpt, file_unit, ierr
    character(len=33) :: header
    character(len=9) :: cdate, ctime

    if (timing_level > 1) call io_stopwatch_start('hamiltonian: write_hr', timer)

    ! write the  whole matrix with all the indices

    open (newunit=file_unit, file=trim(seedname)//'_hr.dat', form='formatted', status='unknown', &
          iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error: hamiltonian_write_hr: problem opening file '//trim(seedname)//'_hr.dat', comm)
      return
    endif

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
    if (timing_level > 1) call io_stopwatch_stop('hamiltonian: write_hr', timer)
  end subroutine hamiltonian_write_hr

  !================================================!
  subroutine hamiltonian_wigner_seitz(ws_region, print_output, real_lattice, irvec, mp_grid, &
                                      ndegen, nrpts, rpt_origin, stdout, timer, error, count_pts, &
                                      comm)
    !================================================!
    !! Calculates a grid of points that fall inside of (and eventually on the
    !! surface of) the Wigner-Seitz supercell centered on the origin of the B
    !! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3
    !================================================!

    use w90_constants, only: eps8
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_utility, only: utility_metric
    use w90_types, only: print_output_type, ws_region_type, timer_list_type

    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points

    implicit none

    ! arguments
    type(ws_region_type), intent(in)    :: ws_region
    type(print_output_type), intent(in) :: print_output
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(inout)              :: nrpts
    integer, intent(inout), allocatable :: ndegen(:)
    integer, intent(inout), allocatable :: irvec(:, :)
    integer, intent(inout)              :: rpt_origin
    integer, intent(in)                 :: mp_grid(3)
    integer, intent(in)                 :: stdout

    real(kind=dp), intent(in)           :: real_lattice(3, 3)

    logical, intent(in)                 :: count_pts

    ! local variables
    integer       :: ndiff(3)
    integer       :: n1, n2, n3, i1, i2, i3, icnt, i, j, ierr, dist_dim
    real(kind=dp)              :: tot, dist_min
    real(kind=dp), allocatable :: dist(:)
    real(kind=dp)              :: real_metric(3, 3)

    if (print_output%timing_level > 1) &
      call io_stopwatch_start('hamiltonian: wigner_seitz', timer)

    call utility_metric(real_lattice, real_metric)
    dist_dim = 1
    do i = 1, 3
      dist_dim = dist_dim*((ws_region%ws_search_size(i) + 1)*2 + 1)
    end do
    allocate (dist(dist_dim), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating dist in hamiltonian_wigner_seitz', comm)
      return
    endif

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
    do n1 = -ws_region%ws_search_size(1)*mp_grid(1), ws_region%ws_search_size(1)*mp_grid(1)
      do n2 = -ws_region%ws_search_size(2)*mp_grid(2), ws_region%ws_search_size(2)*mp_grid(2)
        do n3 = -ws_region%ws_search_size(3)*mp_grid(3), ws_region%ws_search_size(3)*mp_grid(3)
          ! Loop over the lattice vectors R of the Born-von Karman supercell
          ! that contains all the points of the previous loop.
          ! There are (2*(ws_search_size+1)+1)**3 points R. R=0 corresponds to
          ! i1=i2=i3=0, or icnt=((2*(ws_search_size+1)+1)**3 + 1)/2
          icnt = 0
          do i1 = -ws_region%ws_search_size(1) - 1, ws_region%ws_search_size(1) + 1
            do i2 = -ws_region%ws_search_size(2) - 1, ws_region%ws_search_size(2) + 1
              do i3 = -ws_region%ws_search_size(3) - 1, ws_region%ws_search_size(3) + 1
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
          if (abs(dist((dist_dim + 1)/2) - dist_min) .lt. ws_region%ws_distance_tol**2) then
            nrpts = nrpts + 1
            if (.not. count_pts) then
              ndegen(nrpts) = 0
              do i = 1, dist_dim
                if (abs(dist(i) - dist_min) .lt. ws_region%ws_distance_tol**2) &
                  ndegen(nrpts) = ndegen(nrpts) + 1
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
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating dist hamiltonian_wigner_seitz', comm)
      return
    endif
    if (count_pts) then
      if (print_output%timing_level > 1) &
        call io_stopwatch_stop('hamiltonian: wigner_seitz', timer)
      return
    end if

    ! Check the "sum rule"
    tot = 0.0_dp
    do i = 1, nrpts
      tot = tot + 1.0_dp/real(ndegen(i), dp)
    enddo

    if (print_output%iprint >= 3) then
      write (stdout, '(1x,i4,a,/)') nrpts, ' lattice points in Wigner-Seitz supercell:'
      do i = 1, nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', irvec(1, i), irvec(2, i), &
          irvec(3, i), '  degeneracy: ', ndegen(i)
      enddo
      write (stdout, '(1x,a,f12.3)') ' tot = ', tot
      write (stdout, '(1x,a,i12)') ' mp_grid product = ', mp_grid(1)*mp_grid(2)*mp_grid(3)
    endif
    if (abs(tot - real(mp_grid(1)*mp_grid(2)*mp_grid(3), dp)) > eps8) then
      call set_error_fatal(error, 'ERROR in hamiltonian_wigner_seitz: error in finding Wigner-Seitz points', comm)
      return
    endif

    if (print_output%timing_level > 1) call io_stopwatch_stop('hamiltonian: wigner_seitz', timer)

    return

  end subroutine hamiltonian_wigner_seitz

  !================================================!
  subroutine hamiltonian_write_tb(kmesh_info, ham_r, m_matrix, kpt_latt, real_lattice, irvec, &
                                  ndegen, nrpts, num_kpts, num_wann, timing_level, seedname, &
                                  timer, dist_k, error, comm)
    !================================================!
    !! Write in a single file all the information
    !! that is needed to set up a Wannier-based
    !! tight-binding model:
    !! * lattice vectors
    !! * <0n|H|Rn>
    !! * <0n|r|Rn>
    !================================================!

    use w90_io, only: io_stopwatch_start, io_stopwatch_stop, io_date
    use w90_constants, only: twopi, cmplx_i
    use w90_types, only: kmesh_info_type

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: dist_k(:)
    integer, intent(in) :: ndegen(:)
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: irvec(:, :)
    integer, intent(in) :: nrpts
    integer, intent(in) :: timing_level

    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    complex(kind=dp), intent(in) :: ham_r(:, :, :)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)

    character(len=50), intent(in)  :: seedname

    ! local variables
    integer :: ierr
    integer :: i, j, irpt, ik, nn, idir, file_unit
    integer :: rank, ik_rank
    real(kind=dp) :: rdotk
    complex(kind=dp) :: fac, pos_r(3)
    character(len=33) :: header
    character(len=9) :: cdate, ctime
    logical :: on_root = .false.

    rank = mpirank(comm)

    if (rank == 0) on_root = .true.

    if (on_root) then
      if (timing_level > 1) call io_stopwatch_start('hamiltonian: write_tb', timer)

      open (newunit=file_unit, file=trim(seedname)//'_tb.dat', form='formatted', status='unknown', &
            iostat=ierr)
      if (ierr /= 0) then
        call set_error_file(error, 'Error: hamiltonian_write_tb: problem opening file '//trim(seedname)//'_tb.dat', comm)
        return
      endif

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
    endif ! on_root
    !
    ! <0n|r|Rm>
    !
    do irpt = 1, nrpts
      if (on_root) write (file_unit, '(/,3I5)') irvec(:, irpt)
      do i = 1, num_wann
        do j = 1, num_wann
          pos_r(:) = 0._dp
          ik_rank = 1
          do ik = 1, num_kpts
            if (dist_k(ik) /= rank) cycle

            rdotk = twopi*dot_product(kpt_latt(:, ik), real(irvec(:, irpt), dp))
            fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
            do idir = 1, 3
              do nn = 1, kmesh_info%nntot
                if (i == j) then
                  ! For irpt==rpt_origin, this reduces to
                  ! Eq.(32) of Marzari and Vanderbilt PRB 56,
                  ! 12847 (1997). Otherwise, is is Eq.(44)
                  ! Wang, Yates, Souza and Vanderbilt PRB 74,
                  ! 195118 (2006), modified according to
                  ! Eqs.(27,29) of Marzari and Vanderbilt
                  pos_r(idir) = pos_r(idir) - kmesh_info%wb(nn)*kmesh_info%bk(idir, nn, ik) &
                                *aimag(log(m_matrix(i, i, nn, ik_rank)))*fac
                else
                  ! Eq.(44) Wang, Yates, Souza and Vanderbilt PRB 74, 195118 (2006)
                  pos_r(idir) = pos_r(idir) + cmplx_i*kmesh_info%wb(nn) &
                                *kmesh_info%bk(idir, nn, ik)*m_matrix(j, i, nn, ik_rank)*fac
                endif
              end do
            end do
            ik_rank = ik_rank + 1
          end do
          call comms_reduce(pos_r(1), 3, 'SUM', error, comm)
          if (on_root) write (file_unit, '(2I5,3x,6(E15.8,1x))') j, i, pos_r(:)
        end do
      end do
    end do

    if (on_root) then
      close (file_unit)
      if (timing_level > 1) call io_stopwatch_stop('hamiltonian: write_tb', timer)
    endif
  end subroutine hamiltonian_write_tb
end module w90_hamiltonian
