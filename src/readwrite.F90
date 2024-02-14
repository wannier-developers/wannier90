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
!  w90_readwrite: input parsing and information printout     !
!     routines for input/output for data used/needed by      !
!     *both* wannier90.x and postw90.x                       !
!                                                            !
!------------------------------------------------------------!

module w90_readwrite

  !! Common read/write routines for data needed by both
  !! wannier90.x and postw90.x executables

  use w90_constants, only: dp, maxlen
  use w90_types
  use w90_comms, only: w90_comm_type, mpisize

  implicit none

  private

  public :: w90_readwrite_read_distk
  public :: w90_readwrite_chkpt_dist
  public :: w90_readwrite_dealloc
  public :: w90_readwrite_get_convention_type
  public :: w90_readwrite_get_smearing_type
  public :: w90_readwrite_lib_set_atoms
  public :: w90_readwrite_read_chkpt
  public :: w90_readwrite_read_chkpt_header
  public :: w90_readwrite_read_chkpt_matrices
  public :: w90_readwrite_write_header
  public :: w90_readwrite_get_block_length
  public :: w90_readwrite_get_centre_constraints
  public :: w90_readwrite_get_projections
  public :: w90_readwrite_get_range_vector
  public :: w90_readwrite_get_smearing_index
  public :: w90_readwrite_get_vector_length
  public :: w90_readwrite_in_file
  public :: w90_readwrite_set_kmesh
  public :: w90_readwrite_clean_infile
  public :: w90_readwrite_clear_keywords
  public :: w90_readwrite_read_algorithm_control
  public :: w90_readwrite_read_atoms
  public :: w90_readwrite_read_dis_manifold
  public :: w90_readwrite_read_eigvals
  public :: w90_readwrite_read_exclude_bands
  public :: w90_readwrite_read_fermi_energy
  public :: w90_readwrite_read_final_alloc
  public :: w90_readwrite_read_gamma_only
  public :: w90_readwrite_read_kmesh_data
  public :: w90_readwrite_read_kpath
  public :: w90_readwrite_read_kpoints
  public :: w90_readwrite_read_lattice
  public :: w90_readwrite_read_mp_grid
  public :: w90_readwrite_read_num_bands
  public :: w90_readwrite_read_num_wann
  public :: w90_readwrite_read_system
  public :: w90_readwrite_read_units
  public :: w90_readwrite_read_verbosity
  public :: w90_readwrite_read_ws_data

  public :: w90_readwrite_get_keyword
  public :: w90_readwrite_get_keyword_block
  public :: w90_readwrite_get_keyword_vector

  public :: expand_settings
  public :: init_settings

contains
  !================================================!
  subroutine w90_readwrite_read_verbosity(settings, print_output, svd_omega, error, comm)
    use w90_error, only: w90_error_type
    use w90_comms, only: mpirank
    implicit none
    type(print_output_type), intent(inout) :: print_output
    logical, intent(inout) :: svd_omega
    logical :: found
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: unlucky_rank

    call w90_readwrite_get_keyword(settings, 'timing_level', found, error, comm, &
                                   i_value=print_output%timing_level)
    if (allocated(error)) return

    ! special test case; use the timing_level variable to
    ! communicate a kill to remote process, testing error handling
    call w90_readwrite_get_keyword(settings, 'unlucky', found, error, comm, i_value=unlucky_rank)
    if (found) then
      if (unlucky_rank > 0) then
        print_output%timing_level = -unlucky_rank
      endif
    endif

    call w90_readwrite_get_keyword(settings, 'iprint', found, error, comm, &
                                   i_value=print_output%iprint)
    if (allocated(error)) return

    if (print_output%iprint >= 2) svd_omega = .true. ! a printout that does not have its own option flag

    if (mpirank(comm) /= 0) print_output%iprint = 0 ! supress printing non-rank-0
  end subroutine w90_readwrite_read_verbosity

  subroutine w90_readwrite_read_algorithm_control(settings, optimisation, error, comm)
    use w90_error, only: w90_error_type
    implicit none
    integer, intent(inout) :: optimisation
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(settings_type), intent(inout) :: settings

    logical :: found

    call w90_readwrite_get_keyword(settings, 'optimisation', found, error, comm, &
                                   i_value=optimisation)
    if (allocated(error)) return
  end subroutine w90_readwrite_read_algorithm_control

  subroutine w90_readwrite_read_units(settings, lenconfac, length_unit, energy_unit, bohr, error, &
                                      comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    real(kind=dp), intent(inout) :: lenconfac
    character(len=*), intent(inout) :: length_unit
    character(len=*), intent(inout) :: energy_unit
    real(kind=dp), intent(in) :: bohr
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: ic
    logical :: found

    call w90_readwrite_get_keyword(settings, 'energy_unit', found, error, comm, c_value=energy_unit)
    if (allocated(error)) return

    call w90_readwrite_get_keyword(settings, 'length_unit', found, error, comm, c_value=length_unit)
    if (allocated(error)) return
    if (found) then ! if not specified, it may take already have capitalised values, Ang or Bohr
      if (length_unit .ne. 'ang' .and. length_unit .ne. 'bohr') then
        call set_error_input(error, 'Error: value of length_unit not recognised in w90_readwrite_read_units', comm)
        return
      else if (length_unit .eq. 'bohr') then
        lenconfac = 1.0_dp/bohr
      endif
    endif

    ! Length unit (ang --> Ang, bohr --> Bohr)
    ! this is set to uppercase only for printout...
    ic = ichar(length_unit(1:1))
    if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
      length_unit(1:1) = char(ic + ichar('Z') - ichar('z'))

  end subroutine w90_readwrite_read_units

  subroutine w90_readwrite_read_num_wann(settings, num_wann, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    integer, intent(inout) :: num_wann
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    logical :: found

    call w90_readwrite_get_keyword(settings, 'num_wann', found, error, comm, i_value=num_wann)
    if (allocated(error)) return

    if (.not. found) then
      call set_error_input(error, 'Error: You must specify num_wann', comm)
      return
    endif
    if (num_wann <= 0) then
      call set_error_input(error, 'Error: num_wann must be greater than zero', comm)
      return
    endif
  end subroutine w90_readwrite_read_num_wann

  subroutine w90_readwrite_read_distk(settings, distk, nkin, error, comm)
    ! read distribution of kpoints
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none

    integer, allocatable, intent(inout) :: distk(:)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings
    integer, intent(in) :: nkin

    integer :: nk
    logical :: found

    found = .false.

    call w90_readwrite_get_range_vector(settings, 'distk', found, nk, .true., error, comm)
    if (allocated(error)) return

    if (found) then
      if (nk /= nkin) then
        call set_error_input(error, 'incorrect length of k-distribution', comm)
        return
      endif

      allocate (distk(nkin))

      call w90_readwrite_get_range_vector(settings, 'distk', found, nk, .false., error, comm, distk)
      if (allocated(error)) return
    else
      allocate (distk(nkin))

      distk = 0 ! default to no distribution if not specified
    end if
  endsubroutine w90_readwrite_read_distk

  subroutine w90_readwrite_read_exclude_bands(settings, exclude_bands, num_exclude_bands, error, &
                                              comm)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none

    integer, allocatable, intent(inout) :: exclude_bands(:)
    integer, intent(out) :: num_exclude_bands
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: ierr
    logical :: found = .false.

    num_exclude_bands = 0
    call w90_readwrite_get_range_vector(settings, 'exclude_bands', found, num_exclude_bands, &
                                        .true., error, comm)
    if (allocated(error)) return

    if (found) then
      if (num_exclude_bands < 1) then
        call set_error_input(error, 'Error: problem reading exclude_bands', comm)
        return
      endif
      if (allocated(exclude_bands)) deallocate (exclude_bands)
      allocate (exclude_bands(num_exclude_bands), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating exclude_bands in w90_readwrite_read_exclude_bands', comm)
        return
      endif
      call w90_readwrite_get_range_vector(settings, 'exclude_bands', found, num_exclude_bands, &
                                          .false., error, comm, exclude_bands)
      if (allocated(error)) return
      if (any(exclude_bands < 1)) then
        call set_error_input(error, 'Error: exclude_bands must contain positive numbers', comm)
        return
      endif
    end if
  end subroutine w90_readwrite_read_exclude_bands

  subroutine w90_readwrite_read_num_bands(settings, pw90_effective_model, num_bands, num_wann, &
                                          stdout, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    logical, intent(in) :: pw90_effective_model
    integer, intent(inout) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: i_temp
    logical :: found

    call w90_readwrite_get_keyword(settings, 'num_bands', found, error, comm, i_value=i_temp)
    if (allocated(error)) return

    ! in the pw90_effective_model case, this variable is ignored
    if (.not. pw90_effective_model) then
      if (found) then
        num_bands = i_temp
        if (num_bands < num_wann) then
          write (stdout, *) 'num_bands', num_bands
          write (stdout, *) 'num_wann', num_wann
          call set_error_input(error, 'Error: num_bands must be greater than or equal to num_wann', comm)
          return
        endif
      else
        num_bands = num_wann
      endif
    endif
  end subroutine w90_readwrite_read_num_bands

  subroutine w90_readwrite_read_gamma_only(settings, gamma_only, num_kpts, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    logical, intent(inout) :: gamma_only
    integer, intent(in) :: num_kpts
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    logical :: found, ltmp

    ltmp = .false.
    call w90_readwrite_get_keyword(settings, 'gamma_only', found, error, comm, l_value=ltmp)
    if (allocated(error)) return
    gamma_only = ltmp
    if (gamma_only .and. (num_kpts .ne. 1)) then
      call set_error_input(error, 'Error: gamma_only is true, but num_kpts > 1', comm)
      return
    endif
  end subroutine w90_readwrite_read_gamma_only

  subroutine w90_readwrite_read_mp_grid(settings, pw90_effective_model, mp_grid, num_kpts, error, &
                                        comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    logical, intent(in) :: pw90_effective_model
    integer, intent(inout) :: mp_grid(3), num_kpts
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: iv_temp(3)
    logical :: found

    call w90_readwrite_get_keyword_vector(settings, 'mp_grid', found, 3, error, comm, &
                                          i_value=iv_temp)
    if (allocated(error)) return

    ! ignored in pw90_effective_model case
    if (.not. pw90_effective_model) then
      if (found) mp_grid = iv_temp
      if (.not. found) then
        call set_error_input(error, 'Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid', comm)
        return
      elseif (any(mp_grid < 1)) then
        call set_error_input(error, 'Error: mp_grid must be greater than zero', comm)
        return
      end if
      num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
    end if
  end subroutine w90_readwrite_read_mp_grid

  subroutine w90_readwrite_read_system(settings, w90_system, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    type(w90_system_type), intent(inout) :: w90_system
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    logical :: found, ltmp
    integer :: itmp

    ltmp = .false.  ! by default our WF are not spinors
    call w90_readwrite_get_keyword(settings, 'spinors', found, error, comm, l_value=ltmp)
    if (allocated(error)) return
    if (found) then
      w90_system%spinors = ltmp
    else
      w90_system%spinors = .false. ! specify a default behaviour
    endif
    ! We need to know if the bands are double degenerate due to spin, e.g. when
    ! calculating the DOS
    if (w90_system%spinors) then
      w90_system%num_elec_per_state = 1
    else
      w90_system%num_elec_per_state = 2 ! the default
    endif

    call w90_readwrite_get_keyword(settings, 'num_elec_per_state', found, error, comm, &
                                   i_value=itmp)
    if (allocated(error)) return
    if (found) then
      if (itmp /= 1 .and. itmp /= 2) then
        call set_error_input(error, 'Error: num_elec_per_state can be only 1 or 2', comm)
        return
      else
        if (w90_system%spinors .and. itmp /= 1) then
          call set_error_input(error, 'Error: when spinors = T num_elec_per_state must be 1', comm)
          return
        else
          w90_system%num_elec_per_state = itmp
        endif
      endif
    endif

    call w90_readwrite_get_keyword(settings, 'num_valence_bands', found, error, comm, &
                                   i_value=w90_system%num_valence_bands)
    if (allocated(error)) return
    if (found .and. (w90_system%num_valence_bands .le. 0)) then
      call set_error_input(error, 'Error: num_valence_bands should be greater than zero', comm)
      return
    endif
  end subroutine w90_readwrite_read_system

  subroutine w90_readwrite_read_kpath(settings, kpoint_path, ok, bands_plot, error, comm)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none
    logical, intent(in) :: bands_plot
    type(kpoint_path_type), intent(inout) :: kpoint_path
    logical, intent(out) :: ok
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: i_temp, ierr, bands_num_spec_points
    logical :: found

    bands_num_spec_points = 0
    call w90_readwrite_get_block_length(settings, 'kpoint_path', found, i_temp, error, comm)
    if (allocated(error)) return
    if (found) then
      ok = .true.
      bands_num_spec_points = i_temp*2
      if (allocated(kpoint_path%labels)) deallocate (kpoint_path%labels)
      allocate (kpoint_path%labels(bands_num_spec_points), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating labels in w90_wannier90_readwrite_read', comm)
        return
      endif
      if (allocated(kpoint_path%points)) deallocate (kpoint_path%points)
      allocate (kpoint_path%points(3, bands_num_spec_points), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating points in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_keyword_kpath(settings, kpoint_path, error, comm)
      if (allocated(error)) return
    else
      ok = .false.
    end if
    call w90_readwrite_get_keyword(settings, 'bands_num_points', found, error, comm, &
                                   i_value=kpoint_path%num_points_first_segment)
    if (allocated(error)) return
    if (bands_plot) then
      if (kpoint_path%num_points_first_segment < 0) then
        call set_error_input(error, 'Error: bands_num_points must be positive', comm)
        return
      endif
    endif
  end subroutine w90_readwrite_read_kpath

  subroutine w90_readwrite_read_fermi_energy(settings, found_fermi_energy, fermi_energy_list, &
                                             error, comm)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none
    logical, intent(out) :: found_fermi_energy
    real(kind=dp), allocatable, intent(out) :: fermi_energy_list(:)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    real(kind=dp) :: fermi_energy
    logical :: fermi_energy_scan
    real(kind=dp) :: fermi_energy_min
    real(kind=dp) :: fermi_energy_max
    real(kind=dp) :: fermi_energy_step
    integer :: i, ierr, n
    logical :: found

    n = 0
    found_fermi_energy = .false.
    call w90_readwrite_get_keyword(settings, 'fermi_energy', found, error, comm, &
                                   r_value=fermi_energy)
    if (allocated(error)) return
    if (found) then
      found_fermi_energy = .true.
      n = 1
    endif

    fermi_energy_scan = .false.
    call w90_readwrite_get_keyword(settings, 'fermi_energy_min', found, error, comm, &
                                   r_value=fermi_energy_min)
    if (allocated(error)) return
    if (found) then
      if (found_fermi_energy) then
        call set_error_input(error, 'Error: Cannot specify both fermi_energy and fermi_energy_min', comm)
        return
      endif
      fermi_energy_scan = .true.
      fermi_energy_max = fermi_energy_min + 1.0_dp
      call w90_readwrite_get_keyword(settings, 'fermi_energy_max', found, error, comm, &
                                     r_value=fermi_energy_max)
      if (allocated(error)) return
      if (found .and. fermi_energy_max <= fermi_energy_min) then
        call set_error_input(error, 'Error: fermi_energy_max must be larger than fermi_energy_min', comm)
        return
      endif
      fermi_energy_step = 0.01_dp
      call w90_readwrite_get_keyword(settings, 'fermi_energy_step', found, error, comm, &
                                     r_value=fermi_energy_step)
      if (allocated(error)) return
      if (found .and. fermi_energy_step <= 0.0_dp) then
        call set_error_input(error, 'Error: fermi_energy_step must be positive', comm)
        return
      endif
      n = nint(abs((fermi_energy_max - fermi_energy_min)/fermi_energy_step)) + 1
    endif

    if (found_fermi_energy) then
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(1), stat=ierr)
      fermi_energy_list(1) = fermi_energy
    elseif (fermi_energy_scan) then
      if (n .eq. 1) then
        fermi_energy_step = 0.0_dp
      else
        fermi_energy_step = (fermi_energy_max - fermi_energy_min)/real(n - 1, dp)
      endif
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(n), stat=ierr)
      do i = 1, n
        fermi_energy_list(i) = fermi_energy_min + (i - 1)*fermi_energy_step
      enddo
!! AAM_2017-03-27: if fermi_energy* parameters are not set in input file
!! then allocate fermi_energy_list with length 1 and set to zero as default.
    else
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(1), stat=ierr)
      fermi_energy_list(1) = 0.0_dp
    endif
    if (ierr /= 0) then
      call set_error_alloc(error, &
                           'Error allocating fermi_energy_list in w90_readwrite_read_fermi_energy', comm)
      return
    endif
  end subroutine w90_readwrite_read_fermi_energy

  subroutine w90_readwrite_read_ws_data(settings, ws_region, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    type(ws_region_type), intent(inout) :: ws_region
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: i
    logical :: found

    call w90_readwrite_get_keyword(settings, 'use_ws_distance', found, error, comm, &
                                   l_value=ws_region%use_ws_distance)
    if (allocated(error)) return

    call w90_readwrite_get_keyword(settings, 'ws_distance_tol', found, error, comm, &
                                   r_value=ws_region%ws_distance_tol)
    if (allocated(error)) return

    call w90_readwrite_get_vector_length(settings, 'ws_search_size', found, i, error, comm)
    if (allocated(error)) return
    if (found) then
      if (i .eq. 1) then
        call w90_readwrite_get_keyword_vector(settings, 'ws_search_size', found, 1, error, comm, &
                                              i_value=ws_region%ws_search_size)
        if (allocated(error)) return
        ws_region%ws_search_size(2) = ws_region%ws_search_size(1)
        ws_region%ws_search_size(3) = ws_region%ws_search_size(1)
      elseif (i .eq. 3) then
        call w90_readwrite_get_keyword_vector(settings, 'ws_search_size', found, 3, error, comm, &
                                              i_value=ws_region%ws_search_size)
        if (allocated(error)) return
      else
        call set_error_input(error, &
                             'Error: ws_search_size must be provided as either one integer or a vector of three integers', comm)
        return
      end if
      if (any(ws_region%ws_search_size <= 0)) then
        call set_error_input(error, 'Error: ws_search_size elements must be greater than zero', comm)
        return
      endif
    end if
  end subroutine w90_readwrite_read_ws_data

  subroutine w90_readwrite_read_eigvals(eig_found, eigval, num_bands, num_kpts, stdout, &
                                        seedname, error, comm)
    !! Read the eigenvalues from wannier.eig

    use w90_error, only: w90_error_type, set_error_file, set_error_file, set_error_alloc

    implicit none

    character(len=*), intent(in)  :: seedname
    integer, intent(in) :: num_bands, num_kpts
    integer, intent(in) :: stdout
    logical, intent(inout) :: eig_found
    real(kind=dp), intent(inout) :: eigval(:, :)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    ! local
    integer :: i, j, k, n, eig_unit

    inquire (file=trim(seedname)//'.eig', exist=eig_found)
    if (.not. eig_found) then
      call set_error_file(error, 'No '//trim(seedname)//'.eig file found. Needed for disentanglement', comm)
      return
    else
      open (newunit=eig_unit, file=trim(seedname)//'.eig', form='formatted', status='old', err=105)
      do k = 1, num_kpts
        do n = 1, num_bands
          read (eig_unit, *, err=106, end=106) i, j, eigval(n, k)
          if ((i .ne. n) .or. (j .ne. k)) then
            write (stdout, '(a)') 'Found a mismatch in '//trim(seedname)//'.eig'
            write (stdout, '(a,i0,a,i0)') 'Wanted band  : ', n, ' found band  : ', i
            write (stdout, '(a,i0,a,i0)') 'Wanted kpoint: ', k, ' found kpoint: ', j
            write (stdout, '(a)') ' '
            write (stdout, '(a)') 'A common cause of this error is using the wrong'
            write (stdout, '(a)') 'number of bands. Check your input files.'
            write (stdout, '(a)') 'If your pseudopotentials have shallow core states remember'
            write (stdout, '(a)') 'to account for these electrons.'
            write (stdout, '(a)') ' '
            call set_error_file(error, 'w90_wannier90_readwrite_read: mismatch in '//trim(seedname)//'.eig', comm)
            return
          end if
        enddo
      end do
      close (eig_unit)
    end if

    return

105 call set_error_file(error, 'Error: Problem opening eigenvalue file '//trim(seedname)//'.eig', comm)
    return
106 call set_error_file(error, 'Error: Problem reading eigenvalue file '//trim(seedname)//'.eig', comm)
    return
  end subroutine w90_readwrite_read_eigvals

  subroutine w90_readwrite_read_dis_manifold(settings, eig_found, dis_manifold, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none

    ! arguments
    logical, intent(in) :: eig_found
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    ! local
    logical :: found, found2

    call w90_readwrite_get_keyword(settings, 'dis_win_min', found, error, comm, &
                                   r_value=dis_manifold%win_min)
    if (allocated(error)) return

    call w90_readwrite_get_keyword(settings, 'dis_win_max', found, error, comm, &
                                   r_value=dis_manifold%win_max)
    if (allocated(error)) return
    ! eig_found check because eigenvalue reading may reset win_min,max limits
    if (eig_found .and. (dis_manifold%win_max .lt. dis_manifold%win_min)) then
      call set_error_input(error, &
                           'Error: w90_readwrite_read_dis_manifold: check disentanglement windows (win_max < win_min !)', comm)
      return
    endif

    call w90_readwrite_get_keyword(settings, 'dis_froz_max', found, error, comm, &
                                   r_value=dis_manifold%froz_max)
    if (allocated(error)) return
    if (found) then
      dis_manifold%frozen_states = .true.
    end if
    call w90_readwrite_get_keyword(settings, 'dis_froz_min', found2, error, comm, &
                                   r_value=dis_manifold%froz_min)
    if (allocated(error)) return
    if (eig_found) then
      if (dis_manifold%froz_max .lt. dis_manifold%froz_min) then
        call set_error_input(error, 'Error: w90_readwrite_read_dis_manifold: check disentanglement frozen windows', comm)
        return
      endif
      if (found2 .and. .not. found) then
        call set_error_input(error, 'Error: w90_readwrite_read_dis_manifold: found dis_froz_min but not dis_froz_max', comm)
        return
      endif
    endif
  end subroutine w90_readwrite_read_dis_manifold

  subroutine w90_readwrite_read_kmesh_data(settings, kmesh_input, error, comm)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: itmp, ierr
    logical :: found

    call w90_readwrite_get_keyword(settings, 'search_shells', found, error, comm, &
                                   i_value=kmesh_input%search_shells)
    if (allocated(error)) return
    if (kmesh_input%search_shells < 0) then
      call set_error_input(error, 'Error: search_shells must be positive', comm)
      return
    endif

    call w90_readwrite_get_keyword(settings, 'kmesh_tol', found, error, comm, &
                                   r_value=kmesh_input%tol)
    if (allocated(error)) return
    if (kmesh_input%tol < 0.0_dp) then
      call set_error_input(error, 'Error: kmesh_tol must be positive', comm)
      return
    endif

    call w90_readwrite_get_range_vector(settings, 'shell_list', found, kmesh_input%num_shells, &
                                        .true., error, comm)
    if (allocated(error)) return
    if (found) then
      if (kmesh_input%num_shells < 0 .or. kmesh_input%num_shells > max_shells) then
        call set_error_input(error, 'Error: number of shell in shell_list must be between zero and six', comm)
        return
      endif
      !if (allocated(kmesh_input%shell_list)) deallocate (kmesh_input%shell_list)
      allocate (kmesh_input%shell_list(kmesh_input%num_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating shell_list in w90_wannier90_readwrite_read', comm)
        return
      endif
      call w90_readwrite_get_range_vector(settings, 'shell_list', found, kmesh_input%num_shells, &
                                          .false., error, comm, kmesh_input%shell_list)
      if (allocated(error)) return
      if (any(kmesh_input%shell_list < 1)) then
        call set_error_input(error, 'Error: shell_list must contain positive numbers', comm)
        return
      endif
    else
      !if (allocated(kmesh_input%shell_list)) deallocate (kmesh_input%shell_list)
      ! this is the default allocation of the shell_list--used by kmesh_shell_automatic()
      allocate (kmesh_input%shell_list(max_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating shell_list in w90_readwrite_read_kmesh_data', comm)
        return
      endif
    end if

    call w90_readwrite_get_keyword(settings, 'num_shells', found, error, comm, i_value=itmp)
    if (allocated(error)) return
    if (found .and. (itmp /= kmesh_input%num_shells)) then
      call set_error_input(error, 'Error: Found obsolete keyword num_shells. Its value does not agree with shell_list', comm)
      return
    endif

    ! If .true., does not perform the check of B1 of
    ! Marzari, Vanderbild, PRB 56, 12847 (1997)
    ! in kmesh.F90
    ! mainly needed for the interaction with Z2PACK
    ! By default: .false. (perform the tests)
    call w90_readwrite_get_keyword(settings, 'skip_b1_tests', found, error, comm, &
                                   l_value=kmesh_input%skip_B1_tests)
    if (allocated(error)) return
  end subroutine w90_readwrite_read_kmesh_data

  subroutine w90_readwrite_read_kpoints(settings, pw90_effective_model, kpt_latt, num_kpts, bohr, &
                                        error, comm)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc, set_error_dealloc
    implicit none

    ! arguments
    integer, intent(in) :: num_kpts
    logical, intent(in) :: pw90_effective_model
    real(kind=dp), allocatable, intent(out) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: bohr
    type(settings_type), intent(inout) :: settings
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    ! local variables
    real(kind=dp), allocatable :: kpt_cart(:, :)
    integer :: ierr
    logical :: found

    ! pw90_effective_model ignores kpt_cart
    ! this routine allocates the intent(out) kpt_latt

    ierr = 0

    allocate (kpt_latt(3, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating kpt_latt in w90_readwrite_read_kpoints', comm)
      return
    endif

    if (.not. pw90_effective_model) then
      allocate (kpt_cart(3, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating kpt_cart in w90_readwrite_read_kpoints', comm)
        return
      endif

      call w90_readwrite_get_keyword_block(settings, 'kpoints', found, num_kpts, 3, bohr, error, &
                                           comm, r_value=kpt_cart)
      if (allocated(error)) return
      if (.not. found) then
        call set_error_input(error, 'Error: Did not find the kpoint information in the input file', comm)
        return
      endif
      kpt_latt = kpt_cart

      deallocate (kpt_cart, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating kpt_cart in w90_readwrite_read_kpoints', comm)
        return
      endif
    endif
  end subroutine w90_readwrite_read_kpoints

  subroutine w90_readwrite_read_lattice(settings, real_lattice, bohr, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    real(kind=dp), intent(out) :: real_lattice(3, 3)
    real(kind=dp) :: real_lattice_tmp(3, 3)
    real(kind=dp), intent(in) :: bohr
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    logical :: found

    call w90_readwrite_get_keyword_block(settings, 'unit_cell_cart', found, 3, 3, bohr, error, &
                                         comm, r_value=real_lattice_tmp)
    if (allocated(error)) return
    real_lattice = transpose(real_lattice_tmp)
    if (.not. found) then
      call set_error_input(error, 'Error: Did not find the cell information in the input file', comm)
      return
    endif
  end subroutine w90_readwrite_read_lattice

  subroutine w90_readwrite_read_atoms(settings, atom_data, real_lattice, bohr, error, comm)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    type(atom_data_type), intent(inout) :: atom_data
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: bohr
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: i_temp, i_temp2
    logical :: found, found2, lunits

    if (allocated(settings%entries)) return ! don't attempt this read in library mode

    ! Atoms
    call w90_readwrite_get_block_length(settings, 'atoms_frac', found, i_temp, error, comm)
    if (allocated(error)) return
    call w90_readwrite_get_block_length(settings, 'atoms_cart', found2, i_temp2, error, comm, lunits)
    if (allocated(error)) return
    if (found .and. found2) then
      call set_error_input(error, 'Error: Cannot specify both atoms_frac and atoms_cart', comm)
      return
    elseif (found .and. i_temp > 0) then
      lunits = .false.
      atom_data%num_atoms = i_temp
    elseif (found2 .and. i_temp2 > 0) then
      atom_data%num_atoms = i_temp2
      ! when units are specified, one fewer line than the total contains atom information
      if (lunits) atom_data%num_atoms = atom_data%num_atoms - 1
    end if
    if (atom_data%num_atoms > 0) then
      call readwrite_get_atoms(settings, atom_data, lunits, real_lattice, bohr, error, comm)
      if (allocated(error)) return
    end if
  end subroutine w90_readwrite_read_atoms

  subroutine w90_readwrite_clear_keywords(settings, comm)
    ! wannier90.x and postw90.x now only read their own subset of the valid tokens in the ctrl file
    ! checking of the ctrl file is by testing for the presence of any remaining strings in the file
    ! after removing all valid keys.
    !
    ! this routine hoovers up any remaining keys by scanning the ctrl file for (the union of) all
    ! wannier90.x and postw90.x keywords.  The w90_readwrite_get_keyword* functions only assign to optional
    ! arguments: here we call without any, which has the side effect of clearing the input stream.
    !
    ! these lists have been populated using a grep command on the source; it needs to be updated by
    ! hand when the code changes.  There are a lot of keywords; it's not an ideal solution.
    !
    ! (for _vector: just specify zero length)
    ! (for _block: small modification to skip checking/failure when rows=0 )
    !use w90_io, only: io_error
    use w90_error, only: w90_error_type

    implicit none

    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    ! this error is never returned (i.e. errors here are discarded)
    type(w90_error_type), allocatable :: error

    logical :: found
    integer :: lx
    integer, allocatable :: lxa(:)

    ! keywords for wannier.x

    call clear_block(settings, 'atoms_cart', error, comm)
    call clear_block(settings, 'atoms_frac', error, comm)
    call clear_block(settings, 'dis_spheres', error, comm)
    call clear_block(settings, 'kpoint_path', error, comm)
    call clear_block(settings, 'kpoints', error, comm)
    call clear_block(settings, 'nnkpts', error, comm)
    call clear_block(settings, 'projections', error, comm)
    call clear_block(settings, 'slwf_centres', error, comm)
    call clear_block(settings, 'unit_cell_cart', error, comm)
    call w90_readwrite_get_keyword(settings, 'auto_projections', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'bands_num_points', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'bands_plot_dim', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'bands_plot_format', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'bands_plot', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'bands_plot_mode', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'calc_only_A', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'conv_noise_amp', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'conv_noise_num', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'conv_tol', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'conv_window', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'cp_pp', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'devel_flag', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_conv_tol', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_conv_window', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_froz_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_froz_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_mix_ratio', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_num_iter', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_spheres_first_wann', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_spheres_num', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dist_cutoff', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dist_cutoff_hc', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dist_cutoff_mode', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_win_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dis_win_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'energy_unit', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fermi_energy', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fermi_energy_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fermi_energy_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fermi_energy_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fermi_surface_num_points', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fermi_surface_plot_format', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fermi_surface_plot', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'fixed_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gamma_only', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'guiding_centres', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'hr_cutoff', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'hr_plot', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'iprint', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kmesh_spacing', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kmesh_tol', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'length_unit', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_bands', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_cg_steps', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_dump_cycles', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_elec_per_state', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_guide_cycles', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_iter', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_no_guide_iter', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_print_cycles', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_shells', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_valence_bands', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'num_wann', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'one_dim_axis', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'optimisation', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'postproc_setup', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'precond', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'restart', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'search_shells', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'site_symmetry', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'skip_b1_tests', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'slwf_constrain', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'slwf_lambda', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'slwf_num', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spin', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spinors', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'symmetrize_eps', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'timing_level', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_easy_fix', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_energy_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_group_threshold', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_bandc', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_bb', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_cc', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_cell_ll', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_cell_rr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_cr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_lc', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_ll', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_num_rr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_read_ht', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'translate_home_cell', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'transport', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'transport_mode', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_use_same_lead', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_win_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_win_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'tran_write_ht', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'trial_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'unlucky', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'use_bloch_phases', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'use_ws_distance', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wannier_plot_format', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wannier_plot', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wannier_plot_mode', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wannier_plot_radius', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wannier_plot_scale', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wannier_plot_spinor_mode', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wannier_plot_spinor_phase', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_bvec', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_hr_diag', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_hr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_proj', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_r2mn', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_rmn', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_tb', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_u_matrices', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_vdw_data', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'write_xyz', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'ws_distance_tol', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wvfn_formatted', found, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'kmesh', found, 0, error, comm) ! the absent arrays have zero length ;-)
    call w90_readwrite_get_keyword_vector(settings, 'mp_grid', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'translation_centre_frac', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'wannier_plot_supercell', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'ws_search_size', found, 0, error, comm)
    call w90_readwrite_get_range_vector(settings, 'bands_plot_project', found, lx, .true., error, comm)
    if (allocated(lxa)) deallocate (lxa); allocate (lxa(lx))
    call w90_readwrite_get_range_vector(settings, 'bands_plot_project', found, lx, .false., error, comm, lxa)
    call w90_readwrite_get_range_vector(settings, 'wannier_plot_list', found, lx, .true., error, comm)
    if (allocated(lxa)) deallocate (lxa); allocate (lxa(lx))
    call w90_readwrite_get_range_vector(settings, 'wannier_plot_list', found, lx, .false., error, comm, lxa)
    call w90_readwrite_get_range_vector(settings, 'select_projections', found, lx, .true., error, comm)
    if (allocated(lxa)) deallocate (lxa); allocate (lxa(lx))
    call w90_readwrite_get_range_vector(settings, 'select_projections', found, lx, .false., error, comm, lxa)
    call w90_readwrite_get_range_vector(settings, 'shell_list', found, lx, .true., error, comm)
    if (allocated(lxa)) deallocate (lxa); allocate (lxa(lx))
    call w90_readwrite_get_range_vector(settings, 'shell_list', found, lx, .false., error, comm, lxa)
    call w90_readwrite_read_exclude_bands(settings, lxa, lx, error, comm)
    ! ends list of wannier.x keywords

    ! keywords for postw90.x
    call w90_readwrite_get_keyword(settings, 'adpt_smr_fac', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'adpt_smr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'adpt_smr_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'berry_curv_adpt_kmesh', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'berry_curv_adpt_kmesh_thresh', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'berry_curv_unit', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'berry', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'berry_kmesh_spacing', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'berry_task', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_2d_dir', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_bandshift_energyshift', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_bandshift_firstband', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_bandshift', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_calc_also_dos', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_adpt_smr_fac', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_adpt_smr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_adpt_smr_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_energy_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_energy_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_energy_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_smr_fixed_en_width', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_dos_smr_type', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_kmesh_spacing', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_mu_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_mu_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_mu_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_relax_time', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_tdf_energy_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_tdf_smr_fixed_en_width', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_tdf_smr_type', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_temp_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_temp_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltz_temp_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'boltzwann', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'degen_thr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_adpt_smr_fac', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_adpt_smr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_adpt_smr_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_energy_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_energy_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_energy_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_kmesh_spacing', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_smr_fixed_en_width', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_smr_type', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'dos_task', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'effective_model', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'geninterp_alsofirstder', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'geninterp', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'geninterp_single_file', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_degen_thresh', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_eigval_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_freq_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_freq_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_freq_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_kmesh_spacing', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_smr_fixed_en_width', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_smr_max_arg', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_smr_type', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'gyrotropic_task', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kpath_bands_colour', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kpath', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kpath_num_points', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kpath_task', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kslice_fermi_lines_colour', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kslice', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kslice_task', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kdotp_num_bands', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_adpt_smr_fac', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_adpt_smr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_adpt_smr_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_eigval_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_freq_max', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_freq_min', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_freq_step', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_smr_fixed_en_width', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'kubo_smr_type', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'sc_eta', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'scissors_shift', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'sc_phase_conv', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'sc_use_eta_corr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'sc_w_thr', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_alpha', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_bandshift_energyshift', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_bandshift_firstband', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_bandshift', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_beta', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_freq_scan', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_gamma', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'shc_method', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'smr_fixed_en_width', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'smr_max_arg', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'smr_type', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spin_axis_azimuth', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spin_axis_polar', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spin_decomp', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spin_kmesh_spacing', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spin_moment', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'spn_formatted', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'transl_inv', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'uhu_formatted', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'use_degen_pert', found, error, comm)
    call w90_readwrite_get_keyword(settings, 'wanint_kpoint_file', found, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'berry_kmesh', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'boltz_kmesh', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'dos_kmesh', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'gyrotropic_box_b1', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'gyrotropic_box_b2', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'gyrotropic_box_b3', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'gyrotropic_box_center', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'gyrotropic_kmesh', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'kdotp_kpoint', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'kslice_2dkmesh', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'kslice_b1', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'kslice_b2', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'kslice_corner', found, 0, error, comm)
    call w90_readwrite_get_keyword_vector(settings, 'spin_kmesh', found, 0, error, comm)
    call w90_readwrite_get_range_vector(settings, 'gyrotropic_band_list', found, lx, .true., error, comm)
    if (allocated(lxa)) deallocate (lxa); allocate (lxa(lx))
    call w90_readwrite_get_range_vector(settings, 'gyrotropic_band_list', found, lx, .false., error, comm, lxa)
    call w90_readwrite_get_range_vector(settings, 'kdotp_bands', found, lx, .true., error, comm)
    if (allocated(lxa)) deallocate (lxa); allocate (lxa(lx))
    call w90_readwrite_get_range_vector(settings, 'kdotp_bands', found, lx, .false., error, comm, lxa)
    call w90_readwrite_get_range_vector(settings, 'dos_project', found, lx, .true., error, comm)
    if (allocated(lxa)) deallocate (lxa); allocate (lxa(lx))
    call w90_readwrite_get_range_vector(settings, 'dos_project', found, lx, .false., error, comm, lxa)
    deallocate (lxa)
    ! ends list of postw90 keywords
    if (allocated(error)) deallocate (error)
  end subroutine w90_readwrite_clear_keywords

  subroutine w90_readwrite_clean_infile(settings, stdout, seedname, error, comm)
    use w90_error, only: w90_error_type, set_error_input, set_error_dealloc
    implicit none
    integer, intent(in) :: stdout
    character(len=*), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: loop, ierr

    ! filter out any remaining accepted keywords from both wannier90.x and postw90.x sets
    ! assumes settings%in_data is allocated
    call w90_readwrite_clear_keywords(settings, comm)

    if (any(len_trim(settings%in_data(:)) > 0)) then
      write (stdout, '(1x,a)') 'The following section of file '//trim(seedname)//'.win contained unrecognised keywords'
      write (stdout, *)
      do loop = 1, settings%num_lines
        if (len_trim(settings%in_data(loop)) > 0) then
          write (stdout, '(1x,a)') trim(settings%in_data(loop))
        end if
      end do
      write (stdout, *)
      call set_error_input(error, 'Unrecognised keyword(s) in input file, see also output file', comm)
      return
    end if

    deallocate (settings%in_data, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating in_data in w90_readwrite_clean_infile', comm)
      return
    endif
    settings%num_lines = 0
  end subroutine w90_readwrite_clean_infile

  subroutine w90_readwrite_read_final_alloc(disentanglement, dis_manifold, wannier_data, num_wann, &
                                            num_bands, num_kpts, error, comm)
    !================================================== !
    ! Some checks and initialisations !
    ! conditionally allocates:
    !   dis_manifold%lwindow(num_bands, num_kpts)
    !   dis_manifold%ndimwin(num_kpts)
    !   dis_manifold%nfirstwin(num_kpts)
    !   wannier_data%centres(3, num_wann)
    !   wannier_data%spreads(num_wann)
    !     small arrays... maybe overkill here?
    !
    !   this is currenty only called by the legacy library (Jun 23)
    !================================================== !
    use w90_error, only: w90_error_type, set_error_alloc
    implicit none
    logical, intent(in) :: disentanglement
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(wannier_data_type), intent(inout) :: wannier_data
    integer, intent(in) :: num_wann, num_bands, num_kpts
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer :: ierr

    if (disentanglement) then
      if (allocated(dis_manifold%ndimwin)) deallocate (dis_manifold%ndimwin)
      allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating ndimwin in w90_wannier90_read_final_alloc()', comm)
        return
      endif
      if (allocated(dis_manifold%nfirstwin)) deallocate (dis_manifold%nfirstwin)
      allocate (dis_manifold%nfirstwin(num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating nfirstwin in w90_wannier90_read_final_alloc()', comm)
        return
      endif
      if (allocated(dis_manifold%lwindow)) deallocate (dis_manifold%lwindow)
      allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating lwindow in w90_wannier90_read_final_alloc()', comm)
        return
      endif
    endif

!    if ( wannier_plot .and. (index(wannier_plot_format,'cub').ne.0) ) then
!       cosa(1)=dot_product(real_lattice(1,:),real_lattice(2,:))
!       cosa(2)=dot_product(real_lattice(1,:),real_lattice(3,:))
!       cosa(3)=dot_product(real_lattice(2,:),real_lattice(3,:))
!       cosa = abs(cosa)
!       if (any(cosa.gt.eps6)) &
!            call io_error('Error: plotting in cube format requires orthogonal lattice vectors')
!    endif

    if (allocated(wannier_data%centres)) deallocate (wannier_data%centres)
    allocate (wannier_data%centres(3, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wannier_centres in w90_wannier90_readwrite_read', comm)
      return
    endif
    wannier_data%centres = 0.0_dp
    if (allocated(wannier_data%spreads)) deallocate (wannier_data%spreads)
    allocate (wannier_data%spreads(num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wannier_spreads in w90_wannier90_readwrite_read', comm)
      return
    endif
    wannier_data%spreads = 0.0_dp
  end subroutine w90_readwrite_read_final_alloc

  subroutine w90_readwrite_set_kmesh(spacing, reclat, mesh)
    !! This routines returns the three integers that define the interpolation k-mesh, satisfying
    !! the condition that the spacing between two neighboring points along each of the three
    !! k_x, k_y and k_z directions is at smaller than a given spacing.
    !!
    !! The reclat is defined as:
    !!   * 'b_1' = (recip_lattice(1,I), i=1,3)
    !!   * 'b_2' = (recip_lattice(2,I), i=1,3)
    !!   * 'b_3' = (recip_lattice(3,I), i=1,3)
    !!
    !!  spacing must be > 0 (and in particular different from zero). We don't check this here.
    !!
    implicit none
    real(kind=dp), intent(in) :: spacing
    !! Minimum spacing between neighboring points, in angstrom^(-1)
    real(kind=dp), intent(in) :: reclat(3, 3)
    !! Matrix of the reciprocal lattice vectors in cartesian coordinates, in angstrom^(-1)
    integer, intent(out) :: mesh(3)
    !! Will contain the three integers defining the interpolation k-mesh

    real(kind=dp) :: blen(3)
    integer :: i

    do i = 1, 3
      blen(i) = sqrt(sum(reclat(i, :)**2))
    end do

    do i = 1, 3
      mesh(i) = int(floor(blen(i)/spacing)) + 1
    end do
  end subroutine w90_readwrite_set_kmesh

  function w90_readwrite_get_smearing_type(smearing_index)
    !! This function returns a string describing the type of smearing
    !! associated to a given smr_index integer value.
    integer, intent(in) :: smearing_index
    !! The integer index for which we want to get the string
    character(len=80) :: w90_readwrite_get_smearing_type
    character(len=4) :: orderstr

    if (smearing_index > 0) then
      write (orderstr, '(I0)') smearing_index
      w90_readwrite_get_smearing_type = "Methfessel-Paxton of order "//trim(orderstr)
    else if (smearing_index .eq. 0) then
      w90_readwrite_get_smearing_type = "Gaussian"
    else if (smearing_index .eq. -1) then
      w90_readwrite_get_smearing_type = "Marzari-Vanderbilt cold smearing"
    else if (smearing_index .eq. -99) then
      w90_readwrite_get_smearing_type = "Fermi-Dirac smearing"
    else
      w90_readwrite_get_smearing_type = "Unknown type of smearing"
    end if

  end function w90_readwrite_get_smearing_type

  function w90_readwrite_get_convention_type(sc_phase_conv)
    !! This function returns a string describing the convention
    !! associated to a sc_phase_conv integer value.
    integer, intent(in) :: sc_phase_conv
    !! The integer index for which we want to get the string
    character(len=80)   :: w90_readwrite_get_convention_type

    !character(len=4)   :: orderstr

    if (sc_phase_conv .eq. 1) then
      w90_readwrite_get_convention_type = "Tight-binding convention"
    else if (sc_phase_conv .eq. 2) then
      w90_readwrite_get_convention_type = "Wannier90 convention"
    else
      w90_readwrite_get_convention_type = "Unknown type of convention"
    end if

  end function w90_readwrite_get_convention_type

  function w90_readwrite_get_smearing_index(string, keyword, error, comm)
    !! This function parses a string containing the type of
    !! smearing and returns the correct index for the smearing_index variable
    !
    !! If the string is not valid, an io_error is issued
    use w90_error, only: w90_error_type, set_error_input
    character(len=*), intent(in) :: string
    !! The string read from input
    character(len=*), intent(in) :: keyword
    !! The keyword that was read (e.g., smr_type), so that we can print a more useful error message
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: w90_readwrite_get_smearing_index
    integer :: pos

    w90_readwrite_get_smearing_index = 0 ! To avoid warnings of unset variables

    if (index(string, 'm-v') > 0) then
      w90_readwrite_get_smearing_index = -1
    elseif (index(string, 'm-p') > 0) then
      pos = index(string, 'm-p')
      if (len(trim(string(pos + 3:))) .eq. 0) then
        ! If the string is only 'm-p', we assume that 'm-p1' was intended
        w90_readwrite_get_smearing_index = 1
      else
        read (string(pos + 3:), *, err=337) w90_readwrite_get_smearing_index
        if (w90_readwrite_get_smearing_index < 0) then
          call set_error_input(error, 'Wrong m-p smearing order in keyword '//trim(keyword), comm)
          return
        endif
      end if
    elseif (index(string, 'f-d') > 0) then
      w90_readwrite_get_smearing_index = -99
      ! Some aliases
    elseif (index(string, 'cold') > 0) then
      w90_readwrite_get_smearing_index = -1
    elseif (index(string, 'gauss') > 0) then
      w90_readwrite_get_smearing_index = 0
      ! Unrecognised keyword
    else
      call set_error_input(error, 'Unrecognised value for keyword '//trim(keyword), comm)
      return
    end if

    return

337 call set_error_input(error, 'Wrong m-p smearing order in keyword '//trim(keyword), comm)
    return

  end function w90_readwrite_get_smearing_index

!================================================
  subroutine w90_readwrite_write_header(bohr_version_str, constants_version_str1, &
                                        constants_version_str2, mpi_size, stdout)
    !! Write a suitable header for the calculation - version authors etc
    use w90_io, only: io_date, w90_version

    implicit none

    integer, intent(in) :: stdout, mpi_size
    character(len=*), intent(in) :: bohr_version_str, constants_version_str1, constants_version_str2
    character(len=9) :: cdate, ctime

    call io_date(cdate, ctime)

    write (stdout, *)
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                   WANNIER90                       |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |        Welcome to the Maximally-Localized         |'
    write (stdout, *) '            |        Generalized Wannier Functions code         |'
    write (stdout, *) '            |            http://www.wannier.org                 |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  Wannier90 Developer Group:                       |'
    write (stdout, *) '            |    Giovanni Pizzi    (EPFL)                       |'
    write (stdout, *) '            |    Valerio Vitale    (Cambridge)                  |'
    write (stdout, *) '            |    David Vanderbilt  (Rutgers University)         |'
    write (stdout, *) '            |    Nicola Marzari    (EPFL)                       |'
    write (stdout, *) '            |    Ivo Souza         (Universidad del Pais Vasco) |'
    write (stdout, *) '            |    Arash A. Mostofi  (Imperial College London)    |'
    write (stdout, *) '            |    Jonathan R. Yates (University of Oxford)       |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  For the full list of Wannier90 3.x authors,      |'
    write (stdout, *) '            |  please check the code documentation and the      |'
    write (stdout, *) '            |  README on the GitHub page of the code            |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  Please cite                                      |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  [ref] "Wannier90 as a community code:            |'
    write (stdout, *) '            |        new features and applications",            |'
    write (stdout, *) '            |        G. Pizzi et al., J. Phys. Cond. Matt. 32,  |'
    write (stdout, *) '            |        165902 (2020).                             |'
    write (stdout, *) '            |        http://doi.org/10.1088/1361-648X/ab51ff    |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  in any publications arising from the use of      |'
    write (stdout, *) '            |  this code. For the method please cite            |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  [ref] "Maximally Localized Generalised Wannier   |'
    write (stdout, *) '            |         Functions for Composite Energy Bands"     |'
    write (stdout, *) '            |         N. Marzari and D. Vanderbilt              |'
    write (stdout, *) '            |         Phys. Rev. B 56 12847 (1997)              |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |  [ref] "Maximally Localized Wannier Functions     |'
    write (stdout, *) '            |         for Entangled Energy Bands"               |'
    write (stdout, *) '            |         I. Souza, N. Marzari and D. Vanderbilt    |'
    write (stdout, *) '            |         Phys. Rev. B 65 035109 (2001)             |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | Copyright (c) 1996-2020                           |'
    write (stdout, *) '            |        The Wannier90 Developer Group and          |'
    write (stdout, *) '            |        individual contributors                    |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            |      Release: ', adjustl(w90_version), '   5th March    2020      |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | This program is free software; you can            |'
    write (stdout, *) '            | redistribute it and/or modify it under the terms  |'
    write (stdout, *) '            | of the GNU General Public License as published by |'
    write (stdout, *) '            | the Free Software Foundation; either version 2 of |'
    write (stdout, *) '            | the License, or (at your option) any later version|'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | This program is distributed in the hope that it   |'
    write (stdout, *) '            | will be useful, but WITHOUT ANY WARRANTY; without |'
    write (stdout, *) '            | even the implied warranty of MERCHANTABILITY or   |'
    write (stdout, *) '            | FITNESS FOR A PARTICULAR PURPOSE. See the GNU     |'
    write (stdout, *) '            | General Public License for more details.          |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            | You should have received a copy of the GNU General|'
    write (stdout, *) '            | Public License along with this program; if not,   |'
    write (stdout, *) '            | write to the Free Software Foundation, Inc.,      |'
    write (stdout, *) '            | 675 Mass Ave, Cambridge, MA 02139, USA.           |'
    write (stdout, *) '            |                                                   |'
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) '            |    Execution started on ', cdate, ' at ', ctime, '    |'
    write (stdout, *) '            +---------------------------------------------------+'
    write (stdout, *) ''
    write (stdout, '(1X,A)') '******************************************************************************'
    write (stdout, '(1X,A)') '* '//constants_version_str1//'*'
    write (stdout, '(1X,A)') '* '//constants_version_str2//'*'
    write (stdout, '(1X,A)') '* '//bohr_version_str//'*'
    write (stdout, '(1X,A)') '******************************************************************************'
    write (stdout, *) ''

    ! show parallel/serial execution
    if (mpi_size == 1) then
#ifdef MPI
      write (stdout, '(/,1x,a)') 'Running in serial (with parallel executable)'
#else
      write (stdout, '(/,1x,a)') 'Running in serial (with serial executable)'
#endif
    else
      write (stdout, '(/,1x,a,i3,a/)') 'Running in parallel on ', mpi_size, ' CPUs'
    endif
  end subroutine w90_readwrite_write_header

!================================================!
  subroutine w90_readwrite_dealloc(exclude_bands, wannier_data, input_proj, kmesh_input, kpt_latt, &
                                   dis_manifold, atom_data, eigval, kpoint_path, error, comm)
    !================================================!
    !! release memory from allocated parameters
    !
    !================================================
    use w90_error, only: w90_error_type, set_error_dealloc

    implicit none

    integer, allocatable, intent(inout) :: exclude_bands(:)
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    type(atom_data_type), intent(inout) :: atom_data
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(kpoint_path_type), intent(inout) :: kpoint_path
    type(proj_type), allocatable, intent(inout) :: input_proj(:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(wannier_data_type), intent(inout) :: wannier_data

    integer :: ierr

    if (allocated(dis_manifold%ndimwin)) then
      deallocate (dis_manifold%ndimwin, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ndimwin in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(dis_manifold%lwindow)) then
      deallocate (dis_manifold%lwindow, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating lwindow in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(eigval)) then
      deallocate (eigval, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating eigval in w90_readwrite_dealloc', comm)
        return
      endif
    endif
    if (allocated(kmesh_input%shell_list)) then
      deallocate (kmesh_input%shell_list, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating shell_list in w90_readwrite_dealloc', comm)
        return
      endif
    endif
    if (allocated(kpt_latt)) then
      deallocate (kpt_latt, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating kpt_latt in w90_readwrite_dealloc', comm)
        return
      endif
    endif
    if (allocated(kpoint_path%labels)) then
      deallocate (kpoint_path%labels, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating labels in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(kpoint_path%points)) then
      deallocate (kpoint_path%points, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating points in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(atom_data%label)) then
      deallocate (atom_data%label, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_label in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(atom_data%symbol)) then
      deallocate (atom_data%symbol, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_symbol in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(atom_data%pos_cart)) then
      deallocate (atom_data%pos_cart, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_pos_cart in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(atom_data%species_num)) then
      deallocate (atom_data%species_num, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_species_num in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(input_proj)) then
      deallocate (input_proj, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(exclude_bands)) then
      deallocate (exclude_bands, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating exclude_bands in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(wannier_data%centres)) then
      deallocate (wannier_data%centres, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wannier_centres in w90_readwrite_dealloc', comm)
        return
      endif
    end if
    if (allocated(wannier_data%spreads)) then
      deallocate (wannier_data%spreads, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wannier_spreads in w90_readwrite_dealloc', comm)
        return
      endif
    endif
  end subroutine w90_readwrite_dealloc

!~  !================================================!
!~  subroutine w90_wannier90_readwrite_write_um
!~    !================================================!
!~    !
!~    ! Dump the U and M to *_um.dat   !
!~    !
!~    !================================================!
!~
!~
!~    use w90_io,        only : io_error,seedname,io_date
!~    implicit none
!~
!~    integer :: i,j,k,l,um_unit
!~    character (len=9) :: cdate, ctime
!~    character(len=33) :: header
!~
!~    call io_date(cdate, ctime)
!~    header='written on '//cdate//' at '//ctime
!~
!~    open(newunit=um_unit,file=trim(seedname)//'_um.dat',form='unformatted')
!~    write(um_unit) header
!~    write(um_unit) omega_invariant
!~    write(um_unit) num_wann,num_kpts,num_nnmax
!~    write(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
!~    write(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
!~    close(um_unit)
!~
!~    return
!~
!~  end subroutine w90_wannier90_readwrite_write_um

!~  !================================================!
!~  subroutine w90_wannier90_readwrite_read_um
!~    !================================================!
!~    !                                !
!~    ! Restore U and M from file      !
!~    !                                !
!~    !================================================!
!~
!~    use w90_io,        only : io_error,seedname
!~    implicit none
!~
!~    integer       :: tmp_num_wann,tmp_num_kpts,tmp_num_nnmax
!~    integer       :: i,j,k,l,um_unit,ierr
!~    character(len=33) :: header
!~    real(kind=dp) :: tmp_omi
!~
!~    open(newunit=um_unit,file=trim(seedname)//'_um.dat',status="old",form='unformatted',err=105)
!~    read(um_unit) header
!~    write(stdout,'(1x,4(a))') 'Reading U and M from file ',trim(seedname),'_um.dat ', header
!~    read(um_unit) tmp_omi
!~    if ( have_disentangled ) then
!~       if ( abs(tmp_omi-omega_invariant).gt.1.0e-10_dp )  &
!~            call io_error('Error in restart: omega_invariant in .chk and um.dat files do not match')
!~    endif
!~    read(um_unit) tmp_num_wann,tmp_num_kpts,tmp_num_nnmax
!~    if(tmp_num_wann/=num_wann) call io_error('Error in w90_wannier90_readwrite_read_um: num_wann mismatch')
!~    if(tmp_num_kpts/=num_kpts) call io_error('Error in w90_wannier90_readwrite_read_um: num_kpts mismatch')
!~    if(tmp_num_nnmax/=num_nnmax) call io_error('Error in w90_wannier90_readwrite_read_um: num_nnmax mismatch')
!~    if (.not.allocated(u_matrix)) then
!~       allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
!~       if (ierr/=0) call io_error('Error allocating u_matrix in w90_wannier90_readwrite_read_um')
!~    endif
!~    read(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
!~    if (.not.allocated(m_matrix)) then
!~       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!~       if (ierr/=0) call io_error('Error allocating m_matrix in w90_wannier90_readwrite_read_um')
!~    endif
!~    read(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
!~    close(um_unit)
!~
!~    return
!~
!~105 call io_error('Error: Problem opening file '//trim(seedname)//'_um.dat in w90_wannier90_readwrite_read_um')
!~
! $  end subroutine w90_wannier90_readwrite_read_um

!================================================!
  subroutine w90_readwrite_read_chkpt(dis_manifold, exclude_bands, kmesh_info, kpt_latt, &
                                      wannier_data, m_matrix, u_matrix, u_matrix_opt, &
                                      real_lattice, omega_invariant, mp_grid, num_bands, &
                                      num_exclude_bands, num_kpts, num_wann, checkpoint, &
                                      have_disentangled, ispostw90, seedname, stdout, error, comm)
    !================================================!
    !! Read checkpoint file
    !! This is used to allocate the matrices.
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !================================================!

    use w90_error, only: w90_error_type, set_error_file, set_error_file, set_error_alloc
    use w90_utility, only: utility_recip_lattice

    implicit none

    ! arguments
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(wannier_data_type), intent(inout) :: wannier_data

    integer, allocatable, intent(inout) :: exclude_bands(:)
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_exclude_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)

    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(inout) :: omega_invariant
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    character(len=*), intent(inout) :: checkpoint
    character(len=*), intent(in)  :: seedname

    logical, intent(in) :: ispostw90 ! Are we running postw90?
    logical, intent(out) :: have_disentangled

    ! local variables
    integer :: chk_unit

    call w90_readwrite_read_chkpt_header(exclude_bands, kmesh_info, kpt_latt, real_lattice, &
                                         mp_grid, num_bands, num_exclude_bands, num_kpts, &
                                         num_wann, checkpoint, have_disentangled, ispostw90, &
                                         seedname, chk_unit, stdout, error, comm)
    if (allocated(error)) return

    call w90_readwrite_read_chkpt_matrices(dis_manifold, kmesh_info, wannier_data, m_matrix, &
                                           u_matrix, u_matrix_opt, omega_invariant, num_bands, &
                                           num_kpts, num_wann, have_disentangled, seedname, &
                                           chk_unit, stdout, error, comm)
  end subroutine w90_readwrite_read_chkpt

!================================================!
  subroutine w90_readwrite_read_chkpt_header(exclude_bands, kmesh_info, kpt_latt, real_lattice, &
                                             mp_grid, num_bands, num_exclude_bands, num_kpts, &
                                             num_wann, checkpoint, have_disentangled, ispostw90, &
                                             seedname, io_unit, stdout, error, comm)
    !================================================!
    !! Read checkpoint file
    !! IMPORTANT! If you change the chkpt format, adapt
    !! accordingly also the w90chk2chk.x utility!
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !================================================!

    use w90_constants, only: eps6
    use w90_error, only: w90_error_type, set_error_file, set_error_file, set_error_alloc
    use w90_utility, only: utility_recip_lattice

    implicit none

    integer, allocatable, intent(inout) :: exclude_bands(:)
    type(kmesh_info_type), intent(in) :: kmesh_info
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_exclude_bands
    integer, intent(out) :: io_unit

    real(kind=dp), intent(in) :: real_lattice(3, 3)

    character(len=*), intent(in)  :: seedname
    character(len=*), intent(inout) :: checkpoint

    logical, intent(in) :: ispostw90 ! Are we running postw90?
    logical, intent(out) :: have_disentangled

    ! local variables
    real(kind=dp) :: recip_lattice(3, 3), volume
    integer :: chk_unit, nkp, i, j, ntmp
    character(len=33) :: header
    real(kind=dp) :: tmp_latt(3, 3), tmp_kpt_latt(3, num_kpts)
    integer :: tmp_excl_bands(1:num_exclude_bands), tmp_mp_grid(1:3)

    write (stdout, '(1x,3a)') 'Reading restart information from file ', trim(seedname), '.chk :'

    open (newunit=chk_unit, file=trim(seedname)//'.chk', status='old', form='unformatted', err=121)
    io_unit = chk_unit

    ! Read comment line
    read (chk_unit) header
    write (stdout, '(1x,a)', advance='no') trim(header)

    ! Consistency checks
    read (chk_unit) ntmp                           ! Number of bands
    if (ntmp .ne. num_bands) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_bands', comm)
      return
    endif
    read (chk_unit) ntmp                           ! Number of excluded bands
    if (ntmp .ne. num_exclude_bands) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_exclude_bands', comm)
      return
    endif
    read (chk_unit) (tmp_excl_bands(i), i=1, num_exclude_bands) ! Excluded bands
    do i = 1, num_exclude_bands
      if (tmp_excl_bands(i) .ne. exclude_bands(i)) then
        call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in exclude_bands', comm)
        return
      endif
    enddo
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Real lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - real_lattice(i, j)) .gt. eps6) then
          call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in real_lattice', comm)
          return
        endif
      enddo
    enddo
    call utility_recip_lattice(real_lattice, recip_lattice, volume, error, comm)
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Reciprocal lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - recip_lattice(i, j)) .gt. eps6) then
          call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in recip_lattice', comm)
          return
        endif
      enddo
    enddo
    read (chk_unit) ntmp                ! K-points
    if (ntmp .ne. num_kpts) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_kpts', comm)
      return
    endif
    read (chk_unit) (tmp_mp_grid(i), i=1, 3)         ! M-P grid
    do i = 1, 3
      if (tmp_mp_grid(i) .ne. mp_grid(i)) then
        call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in mp_grid', comm)
        return
      endif
    enddo
    read (chk_unit) ((tmp_kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts)
    do nkp = 1, num_kpts
      do i = 1, 3
        if (abs(tmp_kpt_latt(i, nkp) - kpt_latt(i, nkp)) .gt. eps6) then
          call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in kpt_latt', comm)
          return
        endif
      enddo
    enddo
    read (chk_unit) ntmp                ! nntot
    if (ntmp .ne. kmesh_info%nntot) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in nntot', comm)
      return
    endif
    read (chk_unit) ntmp                ! num_wann
    if (ntmp .ne. num_wann) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_wann', comm)
      return
    endif
    ! End of consistency checks

    read (chk_unit) checkpoint             ! checkpoint
    checkpoint = adjustl(trim(checkpoint))

    read (chk_unit) have_disentangled      ! whether a disentanglement has been performed

    return

121 if (ispostw90) then
      call set_error_file(error, 'Error opening '//trim(seedname) &
                          //'.chk in w90_readwrite_read_chkpt: did you run wannier90.x first?', comm)
    else
      call set_error_file(error, 'Error opening '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    end if
  end subroutine w90_readwrite_read_chkpt_header

!================================================!
  subroutine w90_readwrite_read_chkpt_matrices(dis_manifold, kmesh_info, wannier_data, m_matrix, &
                                               u_matrix, u_matrix_opt, omega_invariant, num_bands, &
                                               num_kpts, num_wann, have_disentangled, seedname, &
                                               chk_unit, stdout, error, comm)
    !================================================!
    !! Read checkpoint file
    !! IMPORTANT! If you change the chkpt format, adapt
    !! accordingly also the w90chk2chk.x utility!
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !================================================!

    use w90_error, only: w90_error_type, set_error_file, set_error_file, set_error_alloc
    use w90_utility, only: utility_recip_lattice

    implicit none

    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(wannier_data_type), intent(inout) :: wannier_data

    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout, chk_unit

    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)

    real(kind=dp), intent(inout) :: omega_invariant

    character(len=*), intent(in)  :: seedname

    logical, intent(in) :: have_disentangled

    ! local variables
    integer :: nkp, i, j, k, l, ierr

    if (have_disentangled) then

      read (chk_unit) omega_invariant     ! omega invariant

      ! lwindow
      if (.not. allocated(dis_manifold%lwindow)) then
        allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating lwindow in w90_readwrite_read_chkpt', comm)
          return
        endif
      endif
      read (chk_unit, err=122) ((dis_manifold%lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)

      ! ndimwin
      if (.not. allocated(dis_manifold%ndimwin)) then
        allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating ndimwin in w90_readwrite_read_chkpt', comm)
          return
        endif
      endif
      read (chk_unit, err=123) (dis_manifold%ndimwin(nkp), nkp=1, num_kpts)

      ! U_matrix_opt
      read (chk_unit, err=124) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)

    endif

    ! U_matrix
    read (chk_unit, err=125) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)

    ! M_matrix
    read (chk_unit, err=126) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, kmesh_info%nntot), l=1, num_kpts)

    ! wannier_centres
    read (chk_unit, err=127) ((wannier_data%centres(i, j), i=1, 3), j=1, num_wann)

    ! wannier spreads
    read (chk_unit, err=128) (wannier_data%spreads(i), i=1, num_wann)

    close (chk_unit)

    write (stdout, '(a/)') ' ... done'

    return

122 call set_error_file(error, 'Error reading lwindow from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    return
123 call set_error_file(error, 'Error reading ndimwin from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    return
124 call set_error_file(error, 'Error reading u_matrix_opt from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    return
125 call set_error_file(error, 'Error reading u_matrix from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    return
126 call set_error_file(error, 'Error reading m_matrix from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    return
127 call set_error_file(error, 'Error reading wannier_centres from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    return
128 call set_error_file(error, 'Error reading wannier_spreads from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt', comm)
    return
  end subroutine w90_readwrite_read_chkpt_matrices

!================================================!
  subroutine w90_readwrite_chkpt_dist(dis_manifold, wannier_data, u_matrix, u_matrix_opt, &
                                      m_matrix, m_matrix_local, omega_invariant, num_bands, &
                                      num_kpts, num_wann, nntot, checkpoint, have_disentangled, &
                                      distk, error, comm)
    !================================================!
    !
    !! Distribute the chk files
    !
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_date, io_time
    use w90_comms, only: comms_bcast, w90_comm_type, mpirank
    use w90_error, only: w90_error_type, set_error_alloc

    implicit none

    ! arguments
    type(wannier_data_type), intent(inout) :: wannier_data
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: num_kpts
    integer, intent(in) :: nntot
    integer, intent(in) :: distk(:)

    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :) !only alloc/assigned on root
    real(kind=dp), intent(inout) :: omega_invariant

    character(len=*), intent(inout) :: checkpoint
    logical, intent(inout) :: have_disentangled

    ! local variables
    integer :: ierr, ikl, nkl, ikg, rank
    logical :: on_root = .false.

    rank = mpirank(comm)
    if (rank == 0) on_root = .true.

    call comms_bcast(checkpoint, len(checkpoint), error, comm)
    if (allocated(error)) return

    ! assumes u is alloc'd on all nodes
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, error, comm)
    if (allocated(error)) return

    ! assumes m is alloc'd on all nodes
    call comms_bcast(m_matrix(1, 1, 1, 1), num_wann*num_wann*nntot*num_kpts, error, comm)
    if (allocated(error)) return

    ! copy global m into local m
    nkl = count(distk(:) == rank)
    ikl = 1
    do ikg = 1, num_kpts
      if (distk(ikg) == rank) then
        m_matrix_local(:, :, :, ikl) = m_matrix(:, :, :, ikg)
        ikl = ikl + 1
      endif
    enddo

    call comms_bcast(have_disentangled, 1, error, comm)
    if (allocated(error)) return

    if (have_disentangled) then
      if (.not. on_root) then

        if (.not. allocated(dis_manifold%lwindow)) then
          allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating lwindow in w90_readwrite_chkpt_dist', comm)
            return
          endif
        endif

        if (.not. allocated(dis_manifold%ndimwin)) then
          allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating ndimwin in w90_readwrite_chkpt_dist', comm)
            return
          endif
        endif

      end if

      call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts, error, comm)
      if (allocated(error)) return
      call comms_bcast(dis_manifold%lwindow(1, 1), num_bands*num_kpts, error, comm)
      if (allocated(error)) return
      call comms_bcast(dis_manifold%ndimwin(1), num_kpts, error, comm)
      if (allocated(error)) return
      call comms_bcast(omega_invariant, 1, error, comm)
      if (allocated(error)) return
    end if
    call comms_bcast(wannier_data%centres(1, 1), 3*num_wann, error, comm)
    if (allocated(error)) return
    call comms_bcast(wannier_data%spreads(1), num_wann, error, comm)
    if (allocated(error)) return
  end subroutine w90_readwrite_chkpt_dist

!================================================!
  subroutine w90_readwrite_in_file(settings, seedname, error, comm)
    !================================================!
    !! Load the *.win file into a character
    !! array in_file, ignoring comments and
    !! blank lines and converting everything
    !! to lowercase characters
    !================================================!

    use w90_utility, only: utility_lowercase
    use w90_error, only: w90_error_type, set_error_alloc, set_error_file, set_error_file

    implicit none

    character(len=*), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: in_unit, tot_num_lines, ierr, line_counter, loop, in1, in2
    character(len=maxlen) :: dummy
    integer :: pos
    character, parameter :: TABCHAR = char(9)

    open (newunit=in_unit, file=trim(seedname)//'.win', form='formatted', status='old', err=101)

    settings%num_lines = 0; tot_num_lines = 0
    do
      read (in_unit, '(a)', iostat=ierr, err=200, end=210) dummy
      ! [GP-begin, Apr13, 2012]: I convert all tabulation characters to spaces
      pos = index(dummy, TABCHAR)
      do while (pos .ne. 0)
        dummy(pos:pos) = ' '
        pos = index(dummy, TABCHAR)
      end do
      ! [GP-end]
      dummy = adjustl(dummy)
      tot_num_lines = tot_num_lines + 1
      if (.not. dummy(1:1) == '!' .and. .not. dummy(1:1) == '#') then
        if (len(trim(dummy)) > 0) settings%num_lines = settings%num_lines + 1
      endif

    end do

101 call set_error_file(error, 'Error: Problem opening input file '//trim(seedname)//'.win', comm)
    return
200 call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.win', comm)
    return
210 continue
    rewind (in_unit)

    allocate (settings%in_data(settings%num_lines), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating settings%in_data in w90_readwrite_in_file', comm)
      return
    endif

    line_counter = 0
    do loop = 1, tot_num_lines
      read (in_unit, '(a)', iostat=ierr, err=200) dummy
      ! [GP-begin, Apr13, 2012]: I convert all tabulation characters to spaces
      pos = index(dummy, TABCHAR)
      do while (pos .ne. 0)
        dummy(pos:pos) = ' '
        pos = index(dummy, TABCHAR)
      end do
      ! [GP-end]
      dummy = utility_lowercase(dummy)
      dummy = adjustl(dummy)
      if (dummy(1:1) == '!' .or. dummy(1:1) == '#') cycle
      if (len(trim(dummy)) == 0) cycle
      line_counter = line_counter + 1
      in1 = index(dummy, '!')
      in2 = index(dummy, '#')
      if (in1 == 0 .and. in2 == 0) settings%in_data(line_counter) = dummy
      if (in1 == 0 .and. in2 > 0) settings%in_data(line_counter) = dummy(:in2 - 1)
      if (in2 == 0 .and. in1 > 0) settings%in_data(line_counter) = dummy(:in1 - 1)
      if (in2 > 0 .and. in1 > 0) settings%in_data(line_counter) = dummy(:min(in1, in2) - 1)
    end do

    close (in_unit)
  end subroutine w90_readwrite_in_file

  !================================================!
  subroutine w90_readwrite_get_keyword(settings, keyword, found, error, comm, c_value, l_value, i_value, r_value)
    !================================================!
    !
    !! Finds the value of the required keyword.
    !
    !================================================!

    use w90_error, only: w90_error_type, set_error_input, set_error_fatal

    implicit none

    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    character(*), optional, intent(inout) :: c_value
    !! Keyword value
    logical, optional, intent(inout) :: l_value
    !! Keyword value
    integer, optional, intent(inout) :: i_value
    !! Keyword value
    real(kind=dp), optional, intent(inout) :: r_value
    !! Keyword value
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer           :: kl, in, loop, itmp
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    if (allocated(settings%entries) .and. allocated(settings%in_data)) then
      call set_error_fatal(error, 'Error: (library use) options interface and .win parsing clash.'// &
                           '  See library documentation "setting options." (readwrite.F90)', comm)
      return
    elseif (allocated(settings%entries)) then
      do loop = 1, settings%num_entries  ! this means the first occurance of the variable in settings is used
        ! memory beyond num_entries is not initialised
        if (settings%entries(loop)%keyword == keyword) then
          if (present(i_value)) then
            i_value = settings%entries(loop)%idata
          else if (present(r_value)) then
            r_value = settings%entries(loop)%rdata
          else if (present(l_value)) then
            l_value = settings%entries(loop)%ldata
          else if (present(c_value)) then
            c_value = settings%entries(loop)%txtdata
          else
            call set_error_fatal(error, 'Error: keyword sought, but no variable provided to assign to. (readwrite.F90)', comm)
            return
          endif
          found = .true.
          return
        endif
      enddo
    else if (allocated(settings%in_data)) then
      ! by default, scan the input file

      do loop = 1, settings%num_lines
        in = index(settings%in_data(loop), trim(keyword))
        if (in == 0 .or. in > 1) cycle
        itmp = in + len(trim(keyword))
        if (settings%in_data(loop) (itmp:itmp) /= '=' &
            .and. settings%in_data(loop) (itmp:itmp) /= ':' &
            .and. settings%in_data(loop) (itmp:itmp) /= ' ') cycle
        if (found) then
          call set_error_input(error, 'Error: Found keyword '//trim(keyword)//' more than once in input file', comm)
          return
        endif
        found = .true.
        dummy = settings%in_data(loop) (kl + 1:)
        settings%in_data(loop) (1:maxlen) = ' '
        dummy = adjustl(dummy)
        if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
          dummy = dummy(2:)
          dummy = adjustl(dummy)
        end if
      end do

      if (found) then
        if (present(c_value)) c_value = dummy
        if (present(l_value)) then
          if (index(dummy, 't') > 0) then
            l_value = .true.
          elseif (index(dummy, 'f') > 0) then
            l_value = .false.
          else
            call set_error_input(error, 'Error: Problem reading logical keyword '//trim(keyword), comm)
            return
          endif
        endif
        if (present(i_value)) read (dummy, *, err=220, end=220) i_value
        if (present(r_value)) read (dummy, *, err=220, end=220) r_value
      end if
    else
      ! error condition
    endif

    return
220 call set_error_input(error, 'Error: Problem reading keyword '//trim(keyword), comm)
    return
  end subroutine w90_readwrite_get_keyword

  !================================================!
  subroutine w90_readwrite_get_keyword_vector(settings, keyword, found, length, error, comm, &
                                              c_value, l_value, i_value, r_value)
    !================================================!
    !
    !! Finds the values of the required keyword vector
    !
    !================================================!

    use w90_error, only: w90_error_type, set_error_input, set_error_fatal

    implicit none

    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(in)  :: length
    !! Length of vecotr to read
    character(*), optional, intent(inout) :: c_value(length)
    !! Keyword data
    logical, optional, intent(inout) :: l_value(length)
    !! Keyword data
    integer, optional, intent(inout) :: i_value(length)
    !! Keyword data
    real(kind=dp), optional, intent(inout) :: r_value(length)
    !! Keyword data
    type(settings_type), intent(inout) :: settings

    integer           :: kl, in, loop, i
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    if (allocated(settings%entries) .and. allocated(settings%in_data)) then
      call set_error_fatal(error, 'Error: (library use) options interface and .win parsing clash.'// &
                           '  See library documentation "setting options." (readwrite.F90)', comm)
      return
    elseif (allocated(settings%entries)) then

      do loop = 1, settings%num_entries  ! this means the first occurance of the variable in settings is used
        ! memory beyond num_entries is not initialised
        if (settings%entries(loop)%keyword == keyword) then
          if (present(i_value)) then
            i_value = settings%entries(loop)%i1d
          else if (present(r_value)) then
            r_value = settings%entries(loop)%r1d
          else
            call set_error_fatal(error, 'Error: vector sought, but no variable provided to assign to. (readwrite.F90)', comm)
            return
          end if
          found = .true.
        end if
      enddo

    else if (allocated(settings%in_data)) then

      do loop = 1, settings%num_lines
        in = index(settings%in_data(loop), trim(keyword))
        if (in == 0 .or. in > 1) cycle
        if (found) then
          call set_error_input(error, 'Error: Found keyword '//trim(keyword)//' more than once in input file', comm)
          return
        endif
        found = .true.
        dummy = settings%in_data(loop) (kl + 1:)
        settings%in_data(loop) (1:maxlen) = ' '
        dummy = adjustl(dummy)
        if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
          dummy = dummy(2:)
          dummy = adjustl(dummy)
        end if
      end do

      if (found) then
        if (present(c_value)) read (dummy, *, err=230, end=230) (c_value(i), i=1, length)
        if (present(l_value)) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
          call set_error_input(error, 'w90_readwrite_get_keyword_vector unimplemented for logicals', comm)
          return
        endif
        if (present(i_value)) read (dummy, *, err=230, end=230) (i_value(i), i=1, length)
        if (present(r_value)) read (dummy, *, err=230, end=230) (r_value(i), i=1, length)
      end if
    endif

    return

230 call set_error_input(error, 'Error: Problem reading keyword '//trim(keyword)//' in w90_readwrite_get_keyword_vector', comm)
    return
  end subroutine w90_readwrite_get_keyword_vector

!================================================!
  subroutine w90_readwrite_get_vector_length(settings, keyword, found, length, error, comm)
    !================================================!
    !
    !! Returns the length of a keyword vector
    !
    !================================================!

    use w90_error, only: w90_error_type, set_error_input, set_error_fatal

    implicit none

    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(out)  :: length
    !! length of vector
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer           :: kl, in, loop, pos
    character(len=maxlen) :: dummy

    ! get_vector_length only is meaningful for human text in input file
    ! not suitable for data passed via library interface (data in settings%entries)
    !if (.not.allocated(settings%in_data)) then
    !  call set_error_fatal(error, 'w90_readwrite_get_vector_length called with no input file (seeking '//trim(keyword)//')', comm)
    !  return
    !elseif (allocated(settings%entries)) then
    !  call set_error_fatal(error, 'w90_readwrite_get_vector_length called with unspent option arrays', comm)
    !  return
    !endif

    if (allocated(settings%entries) .and. allocated(settings%in_data)) then
      call set_error_fatal(error, 'Error: (library use) options interface and .win parsing clash.'// &
                           '  See library documentation "setting options." (readwrite.F90)', comm)
      return
    elseif (allocated(settings%entries)) then

      do loop = 1, settings%num_entries  ! this means the first occurance of the variable in settings is used
        ! memory beyond num_entries is not initialised
        if (settings%entries(loop)%keyword == keyword) then
          if (allocated(settings%entries(loop)%i1d)) then
            length = size(settings%entries(loop)%i1d)
          else if (allocated(settings%entries(loop)%r1d)) then
            length = size(settings%entries(loop)%r1d)
          else
            call set_error_input(error, 'lib array not i or r', comm)
          endif
          found = .true.
        end if
      enddo

    else if (allocated(settings%in_data)) then

      kl = len_trim(keyword)
      found = .false.
      do loop = 1, settings%num_lines
        in = index(settings%in_data(loop), trim(keyword))
        if (in == 0 .or. in > 1) cycle
        if (found) then
          call set_error_input(error, 'Error: Found keyword '//trim(keyword)//' more than once in input file', comm)
          return
        endif
        found = .true.
        dummy = settings%in_data(loop) (kl + 1:)
        dummy = adjustl(dummy)
        if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
          dummy = dummy(2:)
          dummy = adjustl(dummy)
        end if
      end do

      length = 0
      if (found) then
        if (len_trim(dummy) == 0) then
          call set_error_input(error, 'Error: keyword '//trim(keyword)//' is blank', comm)
          return
        endif
        length = 1
        dummy = adjustl(dummy)
        do
          pos = index(dummy, ' ')
          dummy = dummy(pos + 1:)
          dummy = adjustl(dummy)
          if (len_trim(dummy) > 0) then
            length = length + 1
          else
            exit
          endif
        end do
      end if
    endif ! in_data
  end subroutine w90_readwrite_get_vector_length

  !================================================!
  subroutine w90_readwrite_get_keyword_block(settings, keyword, found, rows, columns, bohr, error, comm, &
                                             c_value, l_value, i_value, r_value)
    !================================================!
    !
    !!   Finds the values of the required data block
    ! i.e. matrix data
    ! applies to: dis_spheres, kpoints, nnkpts, unit_cell_cart
    !
    !================================================!

    use w90_error, only: w90_error_type, set_error_input, set_error_fatal

    implicit none

    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(in)  :: rows
    !! Number of rows
    integer, intent(in)  :: columns
    !! Number of columns
    character(*), optional, intent(inout) :: c_value(columns, rows)
    !! keyword block data
    logical, optional, intent(inout) :: l_value(columns, rows)
    !! keyword block data
    integer, optional, intent(inout) :: i_value(columns, rows)
    !! keyword block data
    real(kind=dp), optional, intent(inout) :: r_value(columns, rows)
    !! keyword block data
    real(kind=dp), intent(in) :: bohr
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: in, ins, ine, loop, i, line_e, line_s, counter, blen
    logical :: found_e, found_s, lconvert
    character(len=maxlen) :: dummy, end_st, start_st

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    if (allocated(settings%entries) .and. allocated(settings%in_data)) then
      call set_error_fatal(error, 'Error: (library use) options interface and .win parsing clash.'// &
                           '  See library documentation "setting options." (readwrite.F90)', comm)
      return
    elseif (allocated(settings%entries)) then

      do loop = 1, settings%num_entries  ! this means the first occurance of the variable in settings is used
        if (settings%entries(loop)%keyword == keyword) then
          if (present(i_value)) then
            i_value = settings%entries(loop)%i2d
          else if (present(r_value)) then
            r_value = settings%entries(loop)%r2d
          else
            call set_error_fatal(error, 'Error: block sought, but no variable provided to assign to. (readwrite.F90)', comm)
            return
          end if
          found = .true.
        end if
      enddo

    else if (allocated(settings%in_data)) then
      do loop = 1, settings%num_lines
        ins = index(settings%in_data(loop), trim(keyword))
        if (ins == 0) cycle
        in = index(settings%in_data(loop), 'begin')
        if (in == 0 .or. in > 1) cycle
        line_s = loop
        if (found_s) then
          call set_error_input(error, 'Error: Found '//trim(start_st)//' more than once in input file', comm)
          return
        endif
        found_s = .true.
      end do

      if (.not. found_s) then
        found = .false.
        return
      end if

      do loop = 1, settings%num_lines
        ine = index(settings%in_data(loop), trim(keyword))
        if (ine == 0) cycle
        in = index(settings%in_data(loop), 'end')
        if (in == 0 .or. in > 1) cycle
        line_e = loop
        if (found_e) then
          call set_error_input(error, 'Error: Found '//trim(end_st)//' more than once in input file', comm)
          return
        endif
        found_e = .true.
      end do

      if (.not. found_e) then
        call set_error_input(error, 'Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', comm)
        return
      end if

      if (line_e <= line_s) then
        call set_error_input(error, 'Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', comm)
        return
      end if

      ! number of lines of data in block
      blen = line_e - line_s - 1

      !    if( blen /= rows) then
      !       if ( index(trim(keyword),'unit_cell_cart').ne.0 ) then
      !          if ( blen /= rows+1 ) call io_error('Error: Wrong number of lines in block '//trim(keyword))
      !       else
      !          call io_error('Error: Wrong number of lines in block '//trim(keyword))
      !       endif
      !    endif

      if ((blen .ne. rows) .and. (blen .ne. rows + 1) .and. (rows .gt. 0)) then
        call set_error_input(error, 'Error: Wrong number of lines in block '//trim(keyword), comm)
        return
      endif

      if ((blen .eq. rows + 1) .and. (rows .gt. 0) .and. &
          (index(trim(keyword), 'unit_cell_cart') .eq. 0)) then
        call set_error_input(error, 'Error: Wrong number of lines in block '//trim(keyword), comm)
        return
      endif

      found = .true.

      lconvert = .false.
      if (blen == rows + 1) then
        dummy = settings%in_data(line_s + 1)
        if (index(dummy, 'ang') .ne. 0) then
          lconvert = .false.
        elseif (index(dummy, 'bohr') .ne. 0) then
          lconvert = .true.
        else
          call set_error_input(error, 'Error: Units in block '//trim(keyword)//' not recognised', comm)
          return
        endif
        settings%in_data(line_s) (1:maxlen) = ' '
        line_s = line_s + 1
      endif

      !    r_value=1.0_dp
      counter = 0
      do loop = line_s + 1, line_e - 1
        dummy = settings%in_data(loop)
        counter = counter + 1
        if (present(c_value)) read (dummy, *, err=240, end=240) (c_value(i, counter), i=1, columns)
        if (present(l_value)) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
          call set_error_input(error, 'w90_readwrite_get_keyword_block unimplemented for logicals', comm)
          return
        endif
        if (present(i_value)) read (dummy, *, err=240, end=240) (i_value(i, counter), i=1, columns)
        if (present(r_value)) read (dummy, *, err=240, end=240) (r_value(i, counter), i=1, columns)
      end do

      if (lconvert) then
        if (present(r_value)) then
          r_value = r_value*bohr
        endif
      endif

      settings%in_data(line_s:line_e) (1:maxlen) = ' '
    endif
    return

240 call set_error_input(error, 'Error: Problem reading block keyword '//trim(keyword), comm)
    return
  end subroutine w90_readwrite_get_keyword_block

  !================================================!
  subroutine w90_readwrite_get_block_length(settings, keyword, found, rows, error, comm, lunits)
    !================================================!
    !
    !! Finds the length of the data block
    !
    !================================================!

    use w90_error, only: w90_error_type, set_error_input, set_error_fatal

    implicit none

    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(out) :: rows
    !! Number of rows
    logical, optional, intent(out) :: lunits
    !! Have we found a unit specification
    type(settings_type), intent(inout) :: settings

    integer           :: i, in, ins, ine, loop, line_e, line_s
    logical           :: found_e, found_s
    character(len=maxlen) :: end_st, start_st, dummy
    character(len=2)  :: atsym
    real(kind=dp)     :: atpos(3)

    ! get_block_length only is meaningful for human text in input file
    ! not suitable for data passed via library interface (data in settings%entries)
    !if (.not. allocated(settings%in_data)) then
    !  call set_error_fatal(error, 'w90_readwrite_get_block_length called with no input file (seeking '//trim(keyword)//')', comm)
    !  return
    !elseif (allocated(settings%entries)) then
    !  call set_error_fatal(error, 'w90_readwrite_get_block_length called with unspent option arrays', comm)
    !  return
    !endif

    found = .false.
    if (allocated(settings%entries)) return ! don't try to do this in library mode

    rows = 0
    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, settings%num_lines
      ins = index(settings%in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(settings%in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call set_error_input(error, 'Error: Found '//trim(start_st)//' more than once in input file', comm)
        return
      endif
      found_s = .true.
    end do

    if (.not. found_s) then
      found = .false.
      return
    end if

    do loop = 1, settings%num_lines
      ine = index(settings%in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(settings%in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call set_error_input(error, 'Error: Found '//trim(end_st)//' more than once in input file', comm)
        return
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call set_error_input(error, 'Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', comm)
      return
    end if

    if (line_e <= line_s) then
      call set_error_input(error, 'Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', comm)
      return
    end if

    rows = line_e - line_s - 1

    found = .true.

    if (present(lunits)) then
      dummy = settings%in_data(line_s + 1)
      read (dummy, *, end=555) atsym, (atpos(i), i=1, 3)
      lunits = .false.
    endif

    if (rows <= 0) then !cope with empty blocks
      found = .false.
      settings%in_data(line_s:line_e) (1:maxlen) = ' '
    end if

    return

555 lunits = .true.

    if (rows <= 1) then !cope with empty blocks
      found = .false.
      settings%in_data(line_s:line_e) (1:maxlen) = ' '
    end if
  end subroutine w90_readwrite_get_block_length

  !================================================!
  subroutine readwrite_get_atoms(settings, atom_data, lunits, real_lattice, bohr, error, comm)
    !================================================!
    !
    !!   Fills the atom data block
    !
    !================================================!

    use w90_utility, only: utility_frac_to_cart, utility_cart_to_frac, utility_inverse_mat
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none

    type(atom_data_type), intent(inout) :: atom_data
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings
    logical, intent(in) :: lunits
    !! Do we expect a first line with the units
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: bohr

    real(kind=dp)     :: inv_lattice(3, 3)
    real(kind=dp)     :: atoms_pos_frac_tmp(3, atom_data%num_atoms)
    real(kind=dp)     :: atoms_pos_cart_tmp(3, atom_data%num_atoms)
    character(len=20) :: keyword
    integer           :: in, ins, ine, loop, i, line_e, line_s, counter
    integer           :: i_temp, loop2, max_sites, ierr, ic
    logical           :: found_e, found_s, found, frac
    character(len=maxlen) :: dummy, end_st, start_st
    character(len=maxlen) :: ctemp(atom_data%num_atoms)
    character(len=maxlen) :: atoms_label_tmp(atom_data%num_atoms)
    logical           :: lconvert

    call utility_inverse_mat(real_lattice, inv_lattice)

    keyword = "atoms_cart"
    frac = .false.
    call w90_readwrite_get_block_length(settings, "atoms_frac", found, i_temp, error, comm)
    if (allocated(error)) return
    if (found) then
      keyword = "atoms_frac"
      frac = .true.
    end if

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, settings%num_lines
      ins = index(settings%in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(settings%in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call set_error_input(error, 'Error: Found '//trim(start_st)//' more than once in input file', comm)
        return
      endif
      found_s = .true.
    end do

    do loop = 1, settings%num_lines
      ine = index(settings%in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(settings%in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call set_error_input(error, 'Error: Found '//trim(end_st)//' more than once in input file', comm)
        return
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call set_error_input(error, 'Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', comm)
      return
    end if

    if (line_e <= line_s) then
      call set_error_input(error, 'Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', comm)
      return
    end if

    lconvert = .false.
    if (lunits) then
      dummy = settings%in_data(line_s + 1)
      if (index(dummy, 'ang') .ne. 0) then
        lconvert = .false.
      elseif (index(dummy, 'bohr') .ne. 0) then
        lconvert = .true.
      else
        call set_error_input(error, 'Error: Units in block atoms_cart not recognised in readwrite_get_atoms', comm)
        return
      endif
      settings%in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    endif

    counter = 0
    do loop = line_s + 1, line_e - 1
      dummy = settings%in_data(loop)
      counter = counter + 1
      if (frac) then
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_frac_tmp(i, counter), i=1, 3)
      else
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_cart_tmp(i, counter), i=1, 3)
      end if
    end do

    if (lconvert) atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr

    settings%in_data(line_s:line_e) (1:maxlen) = ' '

    if (frac) then
      do loop = 1, atom_data%num_atoms
        call utility_frac_to_cart(atoms_pos_frac_tmp(:, loop), atoms_pos_cart_tmp(:, loop), real_lattice)
      end do
    else
      do loop = 1, atom_data%num_atoms
        call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), atoms_pos_frac_tmp(:, loop), inv_lattice)
      end do
    end if

    ! Now we sort the data into the proper structures
    atom_data%num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, atom_data%num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          atom_data%num_species = atom_data%num_species + 1
          ctemp(atom_data%num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atom_data%species_num(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_species_num in readwrite_get_atoms', comm)
      return
    endif
    allocate (atom_data%label(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_label in readwrite_get_atoms', comm)
      return
    endif
    allocate (atom_data%symbol(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_symbol in readwrite_get_atoms', comm)
      return
    endif
    atom_data%species_num(:) = 0

    do loop = 1, atom_data%num_species
      atom_data%label(loop) = ctemp(loop)
      do loop2 = 1, atom_data%num_atoms
        if (trim(atom_data%label(loop)) == trim(atoms_label_tmp(loop2))) then
          atom_data%species_num(loop) = atom_data%species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atom_data%species_num)
    allocate (atom_data%pos_cart(3, max_sites, atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_pos_cart in readwrite_get_atoms', comm)
      return
    endif

    do loop = 1, atom_data%num_species
      counter = 0
      do loop2 = 1, atom_data%num_atoms
        if (trim(atom_data%label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          !atom_data%pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atom_data%pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, atom_data%num_species
      atom_data%symbol(loop) (1:2) = atom_data%label(loop) (1:2)
      ic = ichar(atom_data%symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atom_data%symbol(loop) (2:2) = ' '
    end do

    return

240 call set_error_alloc(error, 'Error: Problem reading block keyword '//trim(keyword), comm)
    return
  end subroutine readwrite_get_atoms

  !================================================!
  subroutine w90_readwrite_lib_set_atoms(atom_data, atoms_label_tmp, atoms_pos_cart_tmp, &
                                         real_lattice, error, comm)
    !================================================!
    !
    !!   Fills the atom data block during a library call
    !
    !================================================!

    use w90_utility, only: utility_cart_to_frac, utility_inverse_mat, utility_lowercase
    use w90_error, only: w90_error_type, set_error_alloc

    implicit none

    type(atom_data_type), intent(inout) :: atom_data
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    character(len=*), intent(in) :: atoms_label_tmp(atom_data%num_atoms)
    !! Atom labels
    real(kind=dp), intent(in) :: atoms_pos_cart_tmp(3, atom_data%num_atoms)
    !! Atom positions
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    real(kind=dp) :: inv_lattice(3, 3)
    real(kind=dp) :: atoms_pos_frac_tmp(3, atom_data%num_atoms)
    integer :: loop2, max_sites, ierr, ic, loop, counter
    character(len=maxlen) :: ctemp(atom_data%num_atoms)
    character(len=maxlen) :: tmp_string

    call utility_inverse_mat(real_lattice, inv_lattice)
    do loop = 1, atom_data%num_atoms
      call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), &
                                atoms_pos_frac_tmp(:, loop), inv_lattice)
    enddo

    ! Now we sort the data into the proper structures
    atom_data%num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, atom_data%num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          atom_data%num_species = atom_data%num_species + 1
          ctemp(atom_data%num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atom_data%species_num(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_species_num in w90_readwrite_lib_set_atoms', comm)
      return
    endif
    allocate (atom_data%label(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_label in w90_readwrite_lib_set_atoms', comm)
      return
    endif
    allocate (atom_data%symbol(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_symbol in w90_readwrite_lib_set_atoms', comm)
      return
    endif
    atom_data%species_num(:) = 0

    do loop = 1, atom_data%num_species
      atom_data%label(loop) = ctemp(loop)
      do loop2 = 1, atom_data%num_atoms
        if (trim(atom_data%label(loop)) == trim(atoms_label_tmp(loop2))) then
          atom_data%species_num(loop) = atom_data%species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atom_data%species_num)
    allocate (atom_data%pos_cart(3, max_sites, atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_pos_cart in w90_readwrite_lib_set_atoms', comm)
      return
    endif

    do loop = 1, atom_data%num_species
      counter = 0
      do loop2 = 1, atom_data%num_atoms
        if (trim(atom_data%label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          !atom_data%pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atom_data%pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, atom_data%num_species
      atom_data%symbol(loop) (1:2) = atom_data%label(loop) (1:2)
      ic = ichar(atom_data%symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atom_data%symbol(loop) (2:2) = ' '
      tmp_string = trim(adjustl(utility_lowercase(atom_data%symbol(loop))))
      atom_data%symbol(loop) (1:2) = tmp_string(1:2)
      tmp_string = trim(adjustl(utility_lowercase(atom_data%label(loop))))
      atom_data%label(loop) (1:2) = tmp_string(1:2)

      ! Upper case the atom labels (eg, si --> Si)
      ic = ichar(atom_data%label(loop) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atom_data%label(loop) (1:1) = char(ic + ichar('Z') - ichar('z'))
    end do
  end subroutine w90_readwrite_lib_set_atoms

  !================================================!
  subroutine w90_readwrite_get_range_vector(settings, keyword, found, length, lcount, error, comm, i_value)
    !================================================!
    !!   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100
    !!   if(lcount) we return the number of states in length
    !================================================!
    use w90_error, only: w90_error_type, set_error_input, set_error_fatal

    implicit none

    character(*), intent(in)    :: keyword
    !! Keyword to examine
    logical, intent(out)   :: found
    !! Is keyword found
    integer, intent(inout) :: length
    !! Number of states
    logical, intent(in)    :: lcount
    !! If T only count states
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    integer, optional, intent(out)   :: i_value(length)
    !! States specified in range vector
    type(settings_type), intent(inout) :: settings

    integer   :: kl, in, loop, num1, num2, i_punc
    integer   :: counter, i_digit, loop_r, range_size
    character(len=maxlen) :: dummy
    character(len=10), parameter :: c_digit = "0123456789"
    character(len=2), parameter :: c_range = "-:"
    character(len=3), parameter :: c_sep = " ,;"
    character(len=5), parameter :: c_punc = " ,;-:"
    character(len=5)  :: c_num1, c_num2

    if (lcount .and. present(i_value)) then
      call set_error_input(error, 'w90_readwrite_get_range_vector: incorrect call', comm)
      return
    endif

    if (allocated(settings%entries)) then ! shortcut for library case
      if (lcount) then
        call w90_readwrite_get_vector_length(settings, keyword, found, length, error, comm)
        return
      else
        call w90_readwrite_get_keyword_vector(settings, keyword, found, length, error, comm, &
                                              i_value=i_value)
        return
      endif
    endif ! end library branch

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, settings%num_lines
      in = index(settings%in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call set_error_input(error, 'Error: Found keyword '//trim(keyword) &
                             //' more than once in input file', comm)
        return
      endif
      found = .true.
      dummy = settings%in_data(loop) (kl + 1:)
      dummy = adjustl(dummy)
      if (.not. lcount) settings%in_data(loop) (1:maxlen) = ' '
      if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
        dummy = dummy(2:)
        dummy = adjustl(dummy)
      end if
    end do

    if (.not. found) return

    counter = 0
    if (len_trim(dummy) == 0) then
      call set_error_input(error, 'Error: keyword '//trim(keyword)//' is blank', comm)
      return
    endif
    dummy = adjustl(dummy)
    do
      i_punc = scan(dummy, c_punc)
      if (i_punc == 0) then
        call set_error_input(error, 'Error parsing keyword '//trim(keyword), comm)
        return
      endif
      c_num1 = dummy(1:i_punc - 1)
      read (c_num1, *, err=101, end=101) num1
      dummy = adjustl(dummy(i_punc:))
      !look for range
      if (scan(dummy, c_range) == 1) then
        i_digit = scan(dummy, c_digit)
        dummy = adjustl(dummy(i_digit:))
        i_punc = scan(dummy, c_punc)
        c_num2 = dummy(1:i_punc - 1)
        read (c_num2, *, err=101, end=101) num2
        dummy = adjustl(dummy(i_punc:))
        range_size = abs(num2 - num1) + 1
        do loop_r = 1, range_size
          counter = counter + 1
          if (.not. lcount) i_value(counter) = min(num1, num2) + loop_r - 1
        end do
      else
        counter = counter + 1
        if (.not. lcount) i_value(counter) = num1
      end if

      if (scan(dummy, c_sep) == 1) dummy = adjustl(dummy(2:))
      if (scan(dummy, c_range) == 1) then
        call set_error_input(error, 'Error parsing keyword '//trim(keyword)//' incorrect range', comm)
        return
      endif
      if (index(dummy, ' ') == 1) exit
    end do

    if (lcount) length = counter
    if (.not. lcount) then
      do loop = 1, counter - 1
        do loop_r = loop + 1, counter
          if (i_value(loop) == i_value(loop_r)) then
            call set_error_input(error, 'Error parsing keyword '//trim(keyword)//' duplicate values', comm)
            return
          endif
        end do
      end do
    end if

    return

101 call set_error_input(error, 'Error parsing keyword '//trim(keyword), comm)
    return
  end subroutine w90_readwrite_get_range_vector

  subroutine w90_readwrite_get_centre_constraints(settings, ccentres_frac, ccentres_cart, &
                                                  proj_site, num_wann, real_lattice, error, comm)
    !================================================!
    !!  assigns projection centres as default centre constraints and global
    !!  Lagrange multiplier as individual Lagrange multipliers then reads
    !!  the centre_constraints block for individual centre constraint parameters
    !
    !================================================!
    use w90_error, only: w90_error_type, set_error_input
    use w90_utility, only: utility_frac_to_cart
    implicit none
    real(kind=dp), intent(inout) :: ccentres_frac(:, :), ccentres_cart(:, :)
    real(kind=dp), intent(in) :: proj_site(:, :)
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    integer :: loop1, index1, constraint_num, loop2
    integer :: column, start, finish, wann
    character(len=maxlen) :: dummy

    do loop1 = 1, num_wann
      do loop2 = 1, 3
        ccentres_frac(loop1, loop2) = proj_site(loop2, loop1)
      end do
    end do

    constraint_num = 0
    do loop1 = 1, settings%num_lines
      dummy = settings%in_data(loop1)
      if (constraint_num > 0) then
        if (trim(dummy) == '') cycle
        index1 = index(dummy, 'begin')
        if (index1 > 0) then
          call set_error_input(error, "slwf_centres block hasn't ended yet", comm)
          return
        endif
        index1 = index(dummy, 'end')
        if (index1 > 0) then
          index1 = index(dummy, 'slwf_centres')
          if (index1 == 0) then
            call set_error_input(error, 'Wrong ending of block (need to end slwf_centres)', comm)
            return
          endif
          settings%in_data(loop1) (1:maxlen) = ' '
          exit
        end if
        column = 0
        start = 1
        finish = 1
        do loop2 = 1, len_trim(dummy)
          if (start == loop2 .and. dummy(loop2:loop2) == ' ') then
            start = loop2 + 1
          end if
          if (start < loop2) then
            if (dummy(loop2:loop2) == ' ') then
              finish = loop2 - 1
              call get_centre_constraint_from_column(column, start, finish, &
                                                     wann, dummy, ccentres_frac, error, comm)
              if (allocated(error)) return
              start = loop2 + 1
              finish = start
            end if
          end if
          if (loop2 == len_trim(dummy) .and. dummy(loop2:loop2) /= ' ') then
            finish = loop2
            call get_centre_constraint_from_column(column, start, finish, &
                                                   wann, dummy, ccentres_frac, error, comm)
            if (allocated(error)) return
            start = loop2 + 1
            finish = start
          end if
        end do
        settings%in_data(loop1) (1:maxlen) = ' '
        constraint_num = constraint_num + 1
      end if
      index1 = index(dummy, 'slwf_centres')
      if (index1 > 0) then
        index1 = index(dummy, 'begin')
        if (index1 > 0) then
          constraint_num = 1
          settings%in_data(loop1) (1:maxlen) = ' '
        end if
      end if
    end do
    do loop1 = 1, num_wann
      call utility_frac_to_cart(ccentres_frac(loop1, :), &
                                ccentres_cart(loop1, :), real_lattice)
    end do
  end subroutine w90_readwrite_get_centre_constraints

  !================================================!
  subroutine get_centre_constraint_from_column(column, start, finish, wann, dummy, ccentres_frac, &
                                               error, comm)
    !================================================!
    !
    !!  assigns value read to constraint
    !!  parameters based on column
    !
    !================================================!
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    integer, intent(inout):: column, start, finish, wann
    character(len=maxlen), intent(inout):: dummy
    real(kind=dp), intent(inout) :: ccentres_frac(:, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    if (column == 0) then
      read (dummy(start:finish), '(i3)') wann
    end if
    if (column > 0) then
      if (column > 4) then
        call set_error_input(error, "Didn't expect anything else after Lagrange multiplier", comm)
        return
      endif
      if (column < 4) read (dummy(start:finish), '(f10.10)') ccentres_frac(wann, column)
    end if
    column = column + 1
  end subroutine get_centre_constraint_from_column

  !================================================!
  subroutine w90_readwrite_get_projections(settings, num_proj, atom_data, num_wann, input_proj, &
                                           inv_lattice, lcount, spinors, bohr, stdout, error, comm)
    !================================================!
    !
    !!  Fills the projection data block
    !
    !================================================!

    use w90_constants, only: eps6, eps2
    use w90_utility, only: utility_cart_to_frac, utility_string_to_coord, utility_strip
    use w90_error, only: w90_error_type, set_error_alloc, set_error_input

    implicit none

    ! arguments
    type(atom_data_type), intent(in) :: atom_data
    type(proj_type), allocatable, intent(inout) :: input_proj(:)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings
    integer, intent(in) :: num_wann
    integer, intent(inout) :: num_proj
    integer, intent(in) :: stdout
    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(in) :: inv_lattice(3, 3)
    logical, intent(in)    :: lcount
    logical, intent(in) :: spinors

    ! local variables
    real(kind=dp)     :: pos_frac(3)
    real(kind=dp)     :: pos_cart(3)
    character(len=20) :: keyword
    integer           :: in, ins, ine, loop, line_e, line_s, counter
    integer           :: sites, species, line, pos1, pos2, pos3, m_tmp, l_tmp, mstate
    integer           :: loop_l, loop_m, loop_sites, ierr, loop_s, spn_counter
    logical           :: found_e, found_s
    character(len=maxlen) :: dummy, end_st, start_st
    character(len=maxlen) :: ctemp, ctemp2, ctemp3, ctemp4, ctemp5, m_string

    integer, parameter :: min_l = -5
    integer, parameter :: max_l = 3
    integer, parameter :: min_m = 1
    integer, parameter :: max_m = 7
    integer            :: ang_states(min_m:max_m, min_l:max_l)
    ! default values for the optional part of the projection definitions
    real(kind=dp), parameter :: proj_z_def(3) = (/0.0_dp, 0.0_dp, 1.0_dp/)
    real(kind=dp), parameter :: proj_x_def(3) = (/1.0_dp, 0.0_dp, 0.0_dp/)
    real(kind=dp), parameter :: proj_s_qaxis_def(3) = (/0.0_dp, 0.0_dp, 1.0_dp/)
    real(kind=dp), parameter :: proj_zona_def = 1.0_dp
    integer, parameter       :: proj_radial_def = 1

    real(kind=dp) :: proj_z_tmp(3)
    real(kind=dp) :: proj_x_tmp(3)
    real(kind=dp) :: proj_s_qaxis_tmp(3)
    real(kind=dp) :: proj_zona_tmp
    integer       :: proj_radial_tmp
    logical       :: lconvert, lrandom, proj_u_tmp, proj_d_tmp
    logical       :: lpartrandom

    real(kind=dp) :: xnorm, znorm, cosphi, sinphi, xnorm_new, cosphi_new

    keyword = "projections"

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    if (.not. lcount) then
      allocate (input_proj(num_proj), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating input_proj in w90_readwrite_get_projections', comm)
        return
      endif
    endif

    do loop = 1, settings%num_lines
      ins = index(settings%in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(settings%in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call set_error_input(error, 'Error: Found '//trim(start_st)//' more than once in input file', comm)
        return
      endif
      found_s = .true.
    end do

    do loop = 1, settings%num_lines
      ine = index(settings%in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(settings%in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call set_error_input(error, &
                             'w90_readwrite_get_projections: Found '//trim(end_st)//' more than once in input file', comm)
        return
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call set_error_input(error, 'w90_readwrite_get_projections: Found '//trim(start_st) &
                           //' but no '//trim(end_st)//' in input file', comm)
      return
    end if

    if (line_e <= line_s) then
      call set_error_input(error, &
                           'w90_readwrite_get_projections: '//trim(end_st)//' comes before '//trim(start_st) &
                           //' in input file', comm)
      return
    end if

    dummy = settings%in_data(line_s + 1)
    lconvert = .false.
    lrandom = .false.
    lpartrandom = .false.
    if (index(dummy, 'ang') .ne. 0) then
      if (.not. lcount) settings%in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    elseif (index(dummy, 'bohr') .ne. 0) then
      if (.not. lcount) settings%in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
      lconvert = .true.
    elseif (index(dummy, 'random') .ne. 0) then
      if (.not. lcount) settings%in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
      if (index(settings%in_data(line_s + 1), end_st) .ne. 0) then
        lrandom = .true.     ! all projections random
      else
        lpartrandom = .true. ! only some projections random
        if (index(settings%in_data(line_s + 1), 'ang') .ne. 0) then
          if (.not. lcount) settings%in_data(line_s) (1:maxlen) = ' '
          line_s = line_s + 1
        elseif (index(settings%in_data(line_s + 1), 'bohr') .ne. 0) then
          if (.not. lcount) settings%in_data(line_s) (1:maxlen) = ' '
          line_s = line_s + 1
          lconvert = .true.
        endif
      endif
    endif

    counter = 0
    if (.not. lrandom) then
      do line = line_s + 1, line_e - 1
        ang_states = 0
        !Assume the default values
        proj_z_tmp = proj_z_def
        proj_x_tmp = proj_x_def
        proj_zona_tmp = proj_zona_def
        proj_radial_tmp = proj_radial_def
        if (spinors) then
          proj_s_qaxis_tmp = proj_s_qaxis_def
          spn_counter = 2
          proj_u_tmp = .true.
          proj_d_tmp = .true.
        else
          spn_counter = 1
        endif
        ! Strip input line of all spaces
        dummy = utility_strip(settings%in_data(line))
        dummy = adjustl(dummy)
        pos1 = index(dummy, ':')
        if (pos1 == 0) then
          call set_error_input(error, &
                               'w90_wannier90_readwrite_read_projection: malformed projection definition: '//trim(dummy), comm)
          return
        endif
        sites = 0
        ctemp = dummy(:pos1 - 1)
        ! Read the atomic site
        if (index(ctemp, 'c=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_cart, error, comm)
          if (allocated(error)) return
          if (lconvert) pos_cart = pos_cart*bohr
          call utility_cart_to_frac(pos_cart(:), pos_frac(:), inv_lattice)
        elseif (index(ctemp, 'f=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_frac, error, comm)
          if (allocated(error)) return
        else
          if (atom_data%num_species == 0) then
            call set_error_input(error, 'w90_wannier90_readwrite_read_projection: ' &
                                 //'Atom centred projection requested but no atoms defined', comm)
            return
          endif
          do loop = 1, atom_data%num_species
            if (trim(ctemp) == atom_data%label(loop)) then
              species = loop
              sites = atom_data%species_num(loop)
              exit
            end if
            if (loop == atom_data%num_species) then
              call set_error_input(error, 'w90_wannier90_readwrite_read_projection: Atom site not recognised '//trim(ctemp), comm)
              return
            endif
          end do
        end if

        dummy = dummy(pos1 + 1:)

        ! scan for quantisation direction
        pos1 = index(dummy, '[')
        if (spinors) then
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 1:))
            pos2 = index(ctemp, ']')
            if (pos2 == 0) then
              call set_error_input(error, &
                                   'w90_readwrite_get_projections: no closing square bracket for spin quantisation dir', comm)
              return
            endif
            ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_s_qaxis_tmp, error, comm)
            dummy = dummy(:pos1 - 1) ! remove [ ] section
          endif
        else
          if (pos1 > 0) then
            call set_error_input(error, 'w90_readwrite_get_projections: spin qdir is defined but spinors=.false.', comm)
            return
          endif
        endif

        ! scan for up or down
        pos1 = index(dummy, '(')
        if (spinors) then
          if (pos1 > 0) then
            proj_u_tmp = .false.; proj_d_tmp = .false.
            ctemp = (dummy(pos1 + 1:))
            pos2 = index(ctemp, ')')
            if (pos2 == 0) then
              call set_error_input(error, 'w90_readwrite_get_projections: no closing bracket for spin', comm)
              return
            endif
            ctemp = ctemp(:pos2 - 1)
            if (index(ctemp, 'u') > 0) proj_u_tmp = .true.
            if (index(ctemp, 'd') > 0) proj_d_tmp = .true.
            if (proj_u_tmp .and. proj_d_tmp) then
              spn_counter = 2
            elseif (.not. proj_u_tmp .and. .not. proj_d_tmp) then
              call set_error_input(error, 'w90_readwrite_get_projections: found brackets but neither u or d', comm)
              return
            else
              spn_counter = 1
            endif
            dummy = dummy(:pos1 - 1) ! remove ( ) section
          endif
        else
          if (pos1 > 0) then
            call set_error_input(error, 'w90_readwrite_get_projections: spin is defined but spinors=.false.', comm)
            return
          endif
        endif

        !Now we know the sites for this line. Get the angular momentum states
        pos1 = index(dummy, ':')
        if (pos1 > 0) then
          ctemp = dummy(:pos1 - 1)
        else
          ctemp = dummy
        end if
        ctemp2 = ctemp
        do
          pos2 = index(ctemp2, ';')
          if (pos2 == 0) then
            ctemp3 = ctemp2
          else
            ctemp3 = ctemp2(:pos2 - 1)
          endif
          if (index(ctemp3, 'l=') == 1) then
            mstate = index(ctemp3, ',')
            if (mstate > 0) then
              read (ctemp3(3:mstate - 1), *, err=101, end=101) l_tmp
            else
              read (ctemp3(3:), *, err=101, end=101) l_tmp
            end if
            if (l_tmp < -5 .or. l_tmp > 3) then
              call set_error_input(error, 'w90_readwrite_get_projections: Incorrect l state requested', comm)
              return
            endif
            if (mstate == 0) then
              if (l_tmp >= 0) then
                do loop_m = 1, 2*l_tmp + 1
                  ang_states(loop_m, l_tmp) = 1
                end do
              elseif (l_tmp == -1) then !sp
                ang_states(1:2, l_tmp) = 1
              elseif (l_tmp == -2) then !sp2
                ang_states(1:3, l_tmp) = 1
              elseif (l_tmp == -3) then !sp3
                ang_states(1:4, l_tmp) = 1
              elseif (l_tmp == -4) then !sp3d
                ang_states(1:5, l_tmp) = 1
              elseif (l_tmp == -5) then !sp3d2
                ang_states(1:6, l_tmp) = 1
              endif
            else
              if (index(ctemp3, 'mr=') /= mstate + 1) then
                call set_error_input(error, 'w90_readwrite_get_projections: Problem reading m state', comm)
                return
              endif
              ctemp4 = ctemp3(mstate + 4:)
              do
                pos3 = index(ctemp4, ',')
                if (pos3 == 0) then
                  ctemp5 = ctemp4
                else
                  ctemp5 = ctemp4(:pos3 - 1)
                endif
                read (ctemp5(1:), *, err=102, end=102) m_tmp
                if (l_tmp >= 0) then
                  if ((m_tmp > 2*l_tmp + 1) .or. (m_tmp <= 0)) then
                    call set_error_input(error, 'w90_readwrite_get_projections: m is > l !', comm)
                    return
                  endif
                elseif (l_tmp == -1 .and. (m_tmp > 2 .or. m_tmp <= 0)) then
                  call set_error_input(error, 'w90_readwrite_get_projections: m has incorrect value (1)', comm)
                  return
                elseif (l_tmp == -2 .and. (m_tmp > 3 .or. m_tmp <= 0)) then
                  call set_error_input(error, 'w90_readwrite_get_projections: m has incorrect value (2)', comm)
                  return
                elseif (l_tmp == -3 .and. (m_tmp > 4 .or. m_tmp <= 0)) then
                  call set_error_input(error, 'w90_readwrite_get_projections: m has incorrect value (3)', comm)
                  return
                elseif (l_tmp == -4 .and. (m_tmp > 5 .or. m_tmp <= 0)) then
                  call set_error_input(error, 'w90_readwrite_get_projections: m has incorrect value (4)', comm)
                  return
                elseif (l_tmp == -5 .and. (m_tmp > 6 .or. m_tmp <= 0)) then
                  call set_error_input(error, 'w90_readwrite_get_projections: m has incorrect value (5)', comm)
                  return
                endif
                ang_states(m_tmp, l_tmp) = 1
                if (pos3 == 0) exit
                ctemp4 = ctemp4(pos3 + 1:)
              enddo
            end if
          else
            do
              pos3 = index(ctemp3, ',')
              if (pos3 == 0) then
                ctemp4 = ctemp3
              else
                ctemp4 = ctemp3(:pos3 - 1)
              endif
              read (ctemp4(1:), *, err=106, end=106) m_string
              select case (trim(adjustl(m_string)))
              case ('s')
                ang_states(1, 0) = 1
              case ('p')
                ang_states(1:3, 1) = 1
              case ('pz')
                ang_states(1, 1) = 1
              case ('px')
                ang_states(2, 1) = 1
              case ('py')
                ang_states(3, 1) = 1
              case ('d')
                ang_states(1:5, 2) = 1
              case ('dz2')
                ang_states(1, 2) = 1
              case ('dxz')
                ang_states(2, 2) = 1
              case ('dyz')
                ang_states(3, 2) = 1
              case ('dx2-y2')
                ang_states(4, 2) = 1
              case ('dxy')
                ang_states(5, 2) = 1
              case ('f')
                ang_states(1:7, 3) = 1
              case ('fz3')
                ang_states(1, 3) = 1
              case ('fxz2')
                ang_states(2, 3) = 1
              case ('fyz2')
                ang_states(3, 3) = 1
              case ('fxyz')
                ang_states(4, 3) = 1
              case ('fz(x2-y2)')
                ang_states(5, 3) = 1
              case ('fx(x2-3y2)')
                ang_states(6, 3) = 1
              case ('fy(3x2-y2)')
                ang_states(7, 3) = 1
              case ('sp')
                ang_states(1:2, -1) = 1
              case ('sp-1')
                ang_states(1, -1) = 1
              case ('sp-2')
                ang_states(2, -1) = 1
              case ('sp2')
                ang_states(1:3, -2) = 1
              case ('sp2-1')
                ang_states(1, -2) = 1
              case ('sp2-2')
                ang_states(2, -2) = 1
              case ('sp2-3')
                ang_states(3, -2) = 1
              case ('sp3')
                ang_states(1:4, -3) = 1
              case ('sp3-1')
                ang_states(1, -3) = 1
              case ('sp3-2')
                ang_states(2, -3) = 1
              case ('sp3-3')
                ang_states(3, -3) = 1
              case ('sp3-4')
                ang_states(4, -3) = 1
              case ('sp3d')
                ang_states(1:5, -4) = 1
              case ('sp3d-1')
                ang_states(1, -4) = 1
              case ('sp3d-2')
                ang_states(2, -4) = 1
              case ('sp3d-3')
                ang_states(3, -4) = 1
              case ('sp3d-4')
                ang_states(4, -4) = 1
              case ('sp3d-5')
                ang_states(5, -4) = 1
              case ('sp3d2')
                ang_states(1:6, -5) = 1
              case ('sp3d2-1')
                ang_states(1, -5) = 1
              case ('sp3d2-2')
                ang_states(2, -5) = 1
              case ('sp3d2-3')
                ang_states(3, -5) = 1
              case ('sp3d2-4')
                ang_states(4, -5) = 1
              case ('sp3d2-5')
                ang_states(5, -5) = 1
              case ('sp3d2-6')
                ang_states(6, -5) = 1
              case default
                call set_error_input(error, 'w90_readwrite_get_projections: Problem reading l state '//trim(ctemp3), comm)
                return
              end select
              if (pos3 == 0) exit
              ctemp3 = ctemp3(pos3 + 1:)
            enddo
          endif
          if (pos2 == 0) exit
          ctemp2 = ctemp2(pos2 + 1:)
        enddo
        ! check for non-default values
        if (pos1 > 0) then
          dummy = dummy(pos1 + 1:)
          ! z axis
          pos1 = index(dummy, 'z=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_z_tmp, error, comm)
            if (allocated(error)) return
          endif
          ! x axis
          pos1 = index(dummy, 'x=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_x_tmp, error, comm)
            if (allocated(error)) return
          endif
          ! diffusivity of orbital
          pos1 = index(dummy, 'zona=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 5:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            read (ctemp, *, err=104, end=104) proj_zona_tmp
          endif
          ! nodes for the radial part
          pos1 = index(dummy, 'r=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            read (ctemp, *, err=105, end=105) proj_radial_tmp
          endif
        end if
        ! if (sites == -1) then
        !   if (counter + spn_counter*sum(ang_states) > num_proj) &
        !     call io_error('w90_readwrite_get_projections: too many projections defined')
        ! else
        !   if (counter + spn_counter*sites*sum(ang_states) > num_proj) &
        !     call io_error('w90_readwrite_get_projections: too many projections defined')
        ! end if

        if (sites == -1) then
          do loop_l = min_l, max_l
            do loop_m = min_m, max_m
              if (ang_states(loop_m, loop_l) == 1) then
                do loop_s = 1, spn_counter
                  counter = counter + 1
                  if (lcount) cycle
                  input_proj(counter)%site(:) = pos_frac
                  input_proj(counter)%l = loop_l
                  input_proj(counter)%m = loop_m
                  input_proj(counter)%z(:) = proj_z_tmp
                  input_proj(counter)%x(:) = proj_x_tmp
                  input_proj(counter)%radial = proj_radial_tmp
                  input_proj(counter)%zona = proj_zona_tmp
                  if (spinors) then
                    if (spn_counter == 1) then
                      if (proj_u_tmp) input_proj(counter)%s = 1
                      if (proj_d_tmp) input_proj(counter)%s = -1
                    else
                      if (loop_s == 1) input_proj(counter)%s = 1
                      if (loop_s == 2) input_proj(counter)%s = -1
                    endif
                    input_proj(counter)%s_qaxis(:) = proj_s_qaxis_tmp
                  endif
                end do
              endif
            end do
          end do
        else
          do loop_sites = 1, sites
            do loop_l = min_l, max_l
              do loop_m = min_m, max_m
                if (ang_states(loop_m, loop_l) == 1) then
                  do loop_s = 1, spn_counter
                    counter = counter + 1
                    if (lcount) cycle
                    call utility_cart_to_frac(atom_data%pos_cart(:, loop_sites, species), &
                                              pos_frac, inv_lattice)
                    input_proj(counter)%site(:) = pos_frac(:)
                    input_proj(counter)%l = loop_l
                    input_proj(counter)%m = loop_m
                    input_proj(counter)%z(:) = proj_z_tmp
                    input_proj(counter)%x(:) = proj_x_tmp
                    input_proj(counter)%radial = proj_radial_tmp
                    input_proj(counter)%zona = proj_zona_tmp
                    if (spinors) then
                      if (spn_counter == 1) then
                        if (proj_u_tmp) input_proj(counter)%s = 1
                        if (proj_d_tmp) input_proj(counter)%s = -1
                      else
                        if (loop_s == 1) input_proj(counter)%s = 1
                        if (loop_s == 2) input_proj(counter)%s = -1
                      endif
                      input_proj(counter)%s_qaxis(:) = proj_s_qaxis_tmp
                    endif
                  end do
                end if
              end do
            end do
          end do
        end if

      end do !end loop over projection block

      ! check there are enough projections and add random projections if required
      if (.not. lpartrandom) then
        if (counter .lt. num_wann) then
          call set_error_input(error, 'w90_readwrite_get_projections: too few projection functions defined', comm)
          return
        endif
      end if
    end if ! .not. lrandom

    if (lcount) then
      if (counter .lt. num_wann) then
        num_proj = num_wann
      else
        num_proj = counter
      endif
      return
    endif

    if (lpartrandom .or. lrandom) then
      call random_seed()  ! comment out this line for reproducible random positions!
      num_proj = num_wann
      ! input_proj is allocated 1:num_wann when 'random'
      do loop = counter + 1, num_wann
        call random_number(input_proj(loop)%site(:))
        input_proj(loop)%l = 0
        input_proj(loop)%m = 1
        input_proj(loop)%z(:) = proj_z_def
        input_proj(loop)%x(:) = proj_x_def
        input_proj(loop)%zona = proj_zona_def
        input_proj(loop)%radial = proj_radial_def
        if (spinors) then
          if (modulo(loop, 2) == 1) then
            input_proj(loop)%s = 1
          else
            input_proj(loop)%s = -1
          end if
          input_proj(loop)%s_qaxis(1) = 0.0_dp
          input_proj(loop)%s_qaxis(2) = 0.0_dp
          input_proj(loop)%s_qaxis(3) = 1.0_dp
        end if
      enddo
    endif

    ! I shouldn't get here, but just in case
    if (.not. lcount) settings%in_data(line_s:line_e) (1:maxlen) = ' '

!~     ! Check
!~     do loop=1,num_wann
!~        if ( abs(sum(proj_z(:,loop)*proj_x(:,loop))).gt.1.0e-6_dp ) then
!~           write(stdout,*) ' Projection:',loop
!~           call io_error(' Error in projections: z and x axes are not orthogonal')
!~        endif
!~     enddo

    ! Normalise z-axis and x-axis and check/fix orthogonality
    do loop = 1, num_proj

      znorm = sqrt(sum(input_proj(loop)%z(:)*input_proj(loop)%z(:)))
      xnorm = sqrt(sum(input_proj(loop)%x(:)*input_proj(loop)%x(:)))
      input_proj(loop)%z(:) = input_proj(loop)%z(:)/znorm             ! normalise z
      input_proj(loop)%x(:) = input_proj(loop)%x(:)/xnorm             ! normalise x
      cosphi = sum(input_proj(loop)%z(:)*input_proj(loop)%x(:))

      ! Check whether z-axis and z-axis are orthogonal
      if (abs(cosphi) .gt. eps6) then

        ! Special case of circularly symmetric projections (pz, dz2, fz3)
        ! just choose an x-axis that is perpendicular to the given z-axis
        if ((input_proj(loop)%l .ge. 0) .and. (input_proj(loop)%m .eq. 1)) then
          proj_x_tmp(:) = input_proj(loop)%x(:)            ! copy of original x-axis
          call random_seed()
          call random_number(proj_z_tmp(:))         ! random vector
          ! calculate new x-axis as the cross (vector) product of random vector with z-axis
          input_proj(loop)%x(1) = proj_z_tmp(2)*input_proj(loop)%z(3) &
                                  - proj_z_tmp(3)*input_proj(loop)%z(2)
          input_proj(loop)%x(2) = proj_z_tmp(3)*input_proj(loop)%z(1) &
                                  - proj_z_tmp(1)*input_proj(loop)%z(3)
          input_proj(loop)%x(3) = proj_z_tmp(1)*input_proj(loop)%z(2) &
                                  - proj_z_tmp(2)*input_proj(loop)%z(1)
          xnorm_new = sqrt(sum(input_proj(loop)%x(:)*input_proj(loop)%x(:)))
          input_proj(loop)%x(:) = input_proj(loop)%x(:)/xnorm_new   ! normalise
          goto 555
        endif

        ! If projection axes non-orthogonal enough, then
        ! user may have made a mistake and should check
        if (abs(cosphi) .gt. eps2) then
          write (stdout, *) ' Projection:', loop
          call set_error_input(error, ' Error in projections: z and x axes are not orthogonal', comm)
          return
        endif

        ! If projection axes are "reasonably orthogonal", project x-axis
        ! onto plane perpendicular to z-axis to make them more so
        sinphi = sqrt(1 - cosphi*cosphi)
        proj_x_tmp(:) = input_proj(loop)%x(:)               ! copy of original x-axis
        ! calculate new x-axis:
        ! x = z \cross (x_tmp \cross z) / sinphi = ( x_tmp - z(z.x_tmp) ) / sinphi
        input_proj(loop)%x(:) = (proj_x_tmp(:) - cosphi*input_proj(loop)%z(:))/sinphi

        ! Final check
555     cosphi_new = sum(input_proj(loop)%z(:)*input_proj(loop)%x(:))
        if (abs(cosphi_new) .gt. eps6) then
          write (stdout, *) ' Projection:'
          call set_error_input(error, ' Error: z and x axes are still not orthogonal after projection', comm)
          return
        endif

      endif

    enddo

    return

101 call set_error_input(error, 'w90_readwrite_get_projections: Problem reading l state into integer '//trim(ctemp3), comm)
    return
102 call set_error_input(error, 'w90_readwrite_get_projections: Problem reading m state into integer '//trim(ctemp3), comm)
    return
104 call set_error_input(error, 'w90_readwrite_get_projections: Problem reading zona into real '//trim(ctemp), comm)
    return
105 call set_error_input(error, 'w90_readwrite_get_projections: Problem reading radial state into integer '//trim(ctemp), comm)
    return
106 call set_error_input(error, 'w90_readwrite_get_projections: Problem reading m state into string '//trim(ctemp3), comm)
    return
  end subroutine w90_readwrite_get_projections

  !================================================!
  subroutine w90_readwrite_get_keyword_kpath(settings, kpoint_path, error, comm)
    !================================================!
    !
    !!  Fills the kpath data block
    !
    !================================================!
    use w90_error, only: w90_error_type, set_error_input

    implicit none

    type(kpoint_path_type), intent(inout) :: kpoint_path
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    type(settings_type), intent(inout) :: settings

    character(len=20) :: keyword
    integer           :: ic, in, ins, ine, loop, inner_loop, i, line_e, line_s, counter
    logical           :: found_e, found_s
    character(len=maxlen) :: dummy, end_st, start_st

    keyword = "kpoint_path"

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, settings%num_lines
      ins = index(settings%in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(settings%in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call set_error_input(error, 'Error: Found '//trim(start_st)//' more than once in input file', comm)
        return
      endif
      found_s = .true.
    end do

    do loop = 1, settings%num_lines
      ine = index(settings%in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(settings%in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call set_error_input(error, 'Error: Found '//trim(end_st)//' more than once in input file', comm)
        return
      endif
      found_e = .true.
    end do

    if (found_s .and. .not. found_e) then
      call set_error_input(error, 'Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', comm)
      return
    end if

    if (found_s .and. found_e) then
      if (line_e <= line_s) then
        call set_error_input(error, 'Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', comm)
        return
      endif
    else
      return !just not found
    end if

    counter = 0
    do loop = line_s + 1, line_e - 1

      counter = counter + 2
      dummy = settings%in_data(loop)
      read (dummy, *, err=240, end=240) kpoint_path%labels(counter - 1), &
        (kpoint_path%points(i, counter - 1), i=1, 3), &
        kpoint_path%labels(counter), (kpoint_path%points(i, counter), i=1, 3)
    end do

    ! Upper case bands labels (eg, x --> X)
    if (allocated(kpoint_path%labels)) then
      do loop = 1, size(kpoint_path%labels)
        do inner_loop = 1, len(kpoint_path%labels(loop))
          ic = ichar(kpoint_path%labels(loop) (inner_loop:inner_loop))
          if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
            kpoint_path%labels(loop) (inner_loop:inner_loop) = char(ic + ichar('Z') - ichar('z'))
        enddo
      enddo
    endif

    settings%in_data(line_s:line_e) (1:maxlen) = ' '

    return

240 call set_error_input(error, 'w90_readwrite_get_keyword_kpath: Problem reading kpath '//trim(dummy), comm)
    return
  end subroutine w90_readwrite_get_keyword_kpath

  !================================================!
  subroutine clear_block(settings, keyword, error, comm)
    !================================================!
    ! a dummy read routine to remove unused but legitimate input block from input stream
    ! needed to preserve input file error checking (i.e. input stream should be empty after all
    ! legitimate keywords/blocks are read)
    !================================================!
    use w90_error, only: w90_error_type, set_error_input

    implicit none

    ! arguments
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    character(len=*), intent(in) :: keyword
    type(settings_type), intent(inout) :: settings

    ! local variables
    integer :: in, ins, ine, loop, line_e, line_s
    logical :: found_e, found_s
    character(len=maxlen) :: end_st, start_st

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, settings%num_lines
      ins = index(settings%in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(settings%in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call set_error_input(error, 'Error: Found '//trim(start_st)//' more than once in input file', comm)
        return
      endif
      found_s = .true.
    end do

    do loop = 1, settings%num_lines
      ine = index(settings%in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(settings%in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call set_error_input(error, 'Error: Found '//trim(end_st)//' more than once in input file', comm)
        return
      endif
      found_e = .true.
    end do

    if (found_s .and. (.not. found_e)) then
      call set_error_input(error, 'Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', comm)
      return
    end if

    if (found_e .and. (.not. found_s)) then
      call set_error_input(error, 'Error: Found '//trim(end_st)//' but no '//trim(start_st)//' in input file', comm)
      return
    end if

    if (found_s .and. found_e) then
      if (line_e <= line_s) then
        call set_error_input(error, 'Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', comm)
        return
      end if

      settings%in_data(line_s:line_e) (1:maxlen) = ' '  ! clear the block from the input stream
    end if ! found tags
  end subroutine clear_block

  subroutine init_settings(settings)
    implicit none
    type(settings_type), intent(inout) :: settings
    integer, parameter :: defsize = 20 ! default size of settings array
    allocate (settings%entries(defsize))
    settings%num_entries = 0
    settings%num_entries_max = defsize
  end subroutine init_settings

  subroutine expand_settings(settings) ! this is a compromise to avoid a fixed size
    type(settings_data), allocatable :: nentries(:); 
    type(settings_type), intent(inout) :: settings
    integer :: n, m ! old, new sizes
    integer, parameter :: incsize = 20 ! default increment when settings array grows
    n = settings%num_entries_max
    m = n + incsize
    allocate (nentries(m)); nentries(1:n) = settings%entries(1:n); call move_alloc(nentries, settings%entries) !f2003, note that "new" space not initialised
    settings%num_entries_max = m
  end subroutine expand_settings

end module w90_readwrite
