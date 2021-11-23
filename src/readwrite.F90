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

  use w90_constants, only: dp
  use w90_io, only: maxlen
  use w90_types

  implicit none

  private

  ! Private data for processing input file
  integer :: num_lines
  character(len=maxlen), allocatable :: in_data(:)

  public :: w90_readwrite_chkpt_dist
  public :: w90_readwrite_dealloc
  public :: w90_readwrite_get_convention_type
  public :: w90_readwrite_get_smearing_type
  public :: w90_readwrite_lib_set_atoms
  public :: w90_readwrite_read_chkpt
  public :: w90_readwrite_write_header
  ! for postw90 parameters
  public :: w90_readwrite_get_block_length
  public :: w90_readwrite_get_centre_constraints
  public :: w90_readwrite_get_keyword
  public :: w90_readwrite_get_keyword_block
  public :: w90_readwrite_get_keyword_vector
  public :: w90_readwrite_get_projections
  public :: w90_readwrite_get_range_vector
  public :: w90_readwrite_get_smearing_index
  public :: w90_readwrite_get_vector_length
  public :: w90_readwrite_in_file
  public :: w90_readwrite_set_kmesh
  public :: w90_readwrite_uppercase
  ! common read routines
  public :: w90_readwrite_clean_infile
  public :: w90_readwrite_clear_keywords
  public :: w90_readwrite_read_algorithm_control
  public :: w90_readwrite_read_atoms
  public :: w90_readwrite_read_devel
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

  private :: clear_block

contains

  !================================================!
  subroutine w90_readwrite_read_verbosity(print_output, stdout, seedname)
    implicit none
    type(print_output_type), intent(inout) :: print_output
    logical :: found
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    print_output%timing_level = 1 ! Verbosity of timing output info
    call w90_readwrite_get_keyword(stdout, seedname, 'timing_level', found, i_value=print_output%timing_level)

    print_output%iprint = 1 ! Verbosity
    call w90_readwrite_get_keyword(stdout, seedname, 'iprint', found, i_value=print_output%iprint)

  end subroutine w90_readwrite_read_verbosity

  subroutine w90_readwrite_read_algorithm_control(optimisation, stdout, seedname)
    implicit none
    integer, intent(inout) :: optimisation
    logical :: found
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    optimisation = 3 ! Verbosity
    call w90_readwrite_get_keyword(stdout, seedname, 'optimisation', found, i_value=optimisation)

  end subroutine w90_readwrite_read_algorithm_control

  subroutine w90_readwrite_read_units(lenconfac, length_unit, energy_unit, bohr, stdout, seedname, &
                                      error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    real(kind=dp), intent(out) :: lenconfac
    integer, intent(in) :: stdout
    character(len=*), intent(out) :: length_unit
    character(len=*), intent(out) :: energy_unit
    character(len=50), intent(in)  :: seedname
    real(kind=dp), intent(in) :: bohr
    type(w90_error_type), allocatable, intent(out) :: error
    logical :: found

    energy_unit = 'ev'
    call w90_readwrite_get_keyword(stdout, seedname, 'energy_unit', found, c_value=energy_unit)

    length_unit = 'ang'
    lenconfac = 1.0_dp
    call w90_readwrite_get_keyword(stdout, seedname, 'length_unit', found, c_value=length_unit)
    if (length_unit .ne. 'ang' .and. length_unit .ne. 'bohr') then
      call set_error_input(error, 'Error: value of length_unit not recognised in w90_readwrite_read_units')
      return
    endif
    if (length_unit .eq. 'bohr') lenconfac = 1.0_dp/bohr
  end subroutine w90_readwrite_read_units

  subroutine w90_readwrite_read_num_wann(num_wann, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    integer, intent(in) :: stdout
    integer, intent(out) :: num_wann
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    logical :: found

    num_wann = -99
    call w90_readwrite_get_keyword(stdout, seedname, 'num_wann', found, i_value=num_wann)
    if (.not. found) then
      call set_error_input(error, 'Error: You must specify num_wann')
      return
    endif
    if (num_wann <= 0) then
      call set_error_input(error, 'Error: num_wann must be greater than zero')
      return
    endif
  end subroutine w90_readwrite_read_num_wann

  subroutine w90_readwrite_read_exclude_bands(exclude_bands, num_exclude_bands, stdout, seedname, &
                                              error)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none

    integer, allocatable, intent(inout) :: exclude_bands(:)
    integer, intent(out) :: num_exclude_bands
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: ierr
    logical :: found

    num_exclude_bands = 0
    call w90_readwrite_get_range_vector(stdout, seedname, 'exclude_bands', found, &
                                        num_exclude_bands, lcount=.true.)
    if (found) then
      if (num_exclude_bands < 1) then
        call set_error_input(error, 'Error: problem reading exclude_bands')
        return
      endif
      if (allocated(exclude_bands)) deallocate (exclude_bands)
      allocate (exclude_bands(num_exclude_bands), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating exclude_bands in w90_readwrite_read_exclude_bands')
        return
      endif
      call w90_readwrite_get_range_vector(stdout, seedname, 'exclude_bands', found, &
                                          num_exclude_bands, .false., exclude_bands)
      if (any(exclude_bands < 1)) then
        call set_error_input(error, 'Error: exclude_bands must contain positive numbers')
        return
      endif
    end if
  end subroutine w90_readwrite_read_exclude_bands

  subroutine w90_readwrite_read_num_bands(pw90_effective_model, library, num_exclude_bands, &
                                          num_bands, num_wann, library_param_read_first_pass, &
                                          stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    logical, intent(in) :: pw90_effective_model, library
    integer, intent(in) :: num_exclude_bands
    integer, intent(inout) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    logical, intent(in) :: library_param_read_first_pass
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: i_temp
    logical :: found

    call w90_readwrite_get_keyword(stdout, seedname, 'num_bands', found, i_value=i_temp)
    if (found .and. library) write (stdout, '(/a)') ' Ignoring <num_bands> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      if (found) num_bands = i_temp
      if (.not. found) num_bands = num_wann
    end if
    ! GP: I subtract it here, but only the first time when I pass the total number of bands
    ! In later calls, I need to pass instead num_bands already subtracted.
    if (library .and. library_param_read_first_pass) num_bands = num_bands - num_exclude_bands
    if (.not. pw90_effective_model) then
      if (found .and. num_bands < num_wann) then
        write (stdout, *) 'num_bands', num_bands
        write (stdout, *) 'num_wann', num_wann
        call set_error_input(error, 'Error: num_bands must be greater than or equal to num_wann')
        return
      endif
    endif
  end subroutine w90_readwrite_read_num_bands

  subroutine w90_readwrite_read_devel(devel_flag, stdout, seedname)
    implicit none
    integer, intent(in) :: stdout
    character(len=*), intent(out) :: devel_flag
    character(len=50), intent(in)  :: seedname

    logical :: found

    devel_flag = ' '
    call w90_readwrite_get_keyword(stdout, seedname, 'devel_flag', found, c_value=devel_flag)
  end subroutine w90_readwrite_read_devel

  subroutine w90_readwrite_read_gamma_only(gamma_only, num_kpts, library, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    integer, intent(in) :: stdout
    logical, intent(inout) :: gamma_only
    integer, intent(in) :: num_kpts
    logical, intent(in) :: library
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    logical :: found, ltmp

    ltmp = .false.
    call w90_readwrite_get_keyword(stdout, seedname, 'gamma_only', found, l_value=ltmp)
    if (.not. library) then
      gamma_only = ltmp
      if (gamma_only .and. (num_kpts .ne. 1)) then
        call set_error_input(error, 'Error: gamma_only is true, but num_kpts > 1')
        return
      endif
    else
      if (found) write (stdout, '(a)') ' Ignoring <gamma_only> in input file'
    endif
  end subroutine w90_readwrite_read_gamma_only

  subroutine w90_readwrite_read_mp_grid(pw90_effective_model, library, mp_grid, num_kpts, stdout, &
                                        seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: pw90_effective_model, library
    integer, intent(inout) :: mp_grid(3), num_kpts
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: iv_temp(3)
    logical :: found

    call w90_readwrite_get_keyword_vector(stdout, seedname, 'mp_grid', found, 3, i_value=iv_temp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <mp_grid> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      if (found) mp_grid = iv_temp
      if (.not. found) then
        call set_error_input(error, 'Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid')
        return
      elseif (any(mp_grid < 1)) then
        call set_error_input(error, 'Error: mp_grid must be greater than zero')
      end if
      num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
    end if
  end subroutine w90_readwrite_read_mp_grid

  subroutine w90_readwrite_read_system(library, w90_system, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    integer, intent(in) :: stdout
    logical, intent(in) :: library
    type(w90_system_type), intent(inout) :: w90_system
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    logical :: found, ltmp

    ltmp = .false.  ! by default our WF are not spinors
    call w90_readwrite_get_keyword(stdout, seedname, 'spinors', found, l_value=ltmp)
    if (.not. library) then
      w90_system%spinors = ltmp
    else
      if (found) write (stdout, '(a)') ' Ignoring <spinors> in input file'
    endif
!    if(spinors .and. (2*(num_wann/2))/=num_wann) &
!       call io_error('Error: For spinor WF num_wann must be even')

    ! We need to know if the bands are double degenerate due to spin, e.g. when
    ! calculating the DOS
    if (w90_system%spinors) then
      w90_system%num_elec_per_state = 1
    else
      w90_system%num_elec_per_state = 2
    endif
    call w90_readwrite_get_keyword(stdout, seedname, 'num_elec_per_state', found, &
                                   i_value=w90_system%num_elec_per_state)
    if ((w90_system%num_elec_per_state /= 1) .and. (w90_system%num_elec_per_state /= 2)) then
      call set_error_input(error, 'Error: num_elec_per_state can be only 1 or 2')
      return
    endif
    if (w90_system%spinors .and. w90_system%num_elec_per_state /= 1) then
      call set_error_input(error, 'Error: when spinors = T num_elec_per_state must be 1')
      return
    endif

    ! set to a negative default value
    w90_system%num_valence_bands = -99
    call w90_readwrite_get_keyword(stdout, seedname, 'num_valence_bands', found, &
                                   i_value=w90_system%num_valence_bands)
    if (found .and. (w90_system%num_valence_bands .le. 0)) then
      call set_error_input(error, 'Error: num_valence_bands should be greater than zero')
      return
    endif
    ! there is a check on this parameter later

  end subroutine w90_readwrite_read_system

  subroutine w90_readwrite_read_kpath(library, kpoint_path, ok, bands_plot, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none
    logical, intent(in) :: library, bands_plot
    type(kpoint_path_type), intent(out) :: kpoint_path
    integer, intent(in) :: stdout
    logical, intent(out) :: ok
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: i_temp, ierr, bands_num_spec_points
    logical :: found

    bands_num_spec_points = 0
    call w90_readwrite_get_block_length(stdout, seedname, 'kpoint_path', found, i_temp, library)
    if (found) then
      ok = .true.
      bands_num_spec_points = i_temp*2
      if (allocated(kpoint_path%labels)) deallocate (kpoint_path%labels)
      allocate (kpoint_path%labels(bands_num_spec_points), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating labels in w90_wannier90_readwrite_read')
        return
      endif
      if (allocated(kpoint_path%points)) deallocate (kpoint_path%points)
      allocate (kpoint_path%points(3, bands_num_spec_points), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating points in w90_wannier90_readwrite_read')
        return
      endif
      call w90_readwrite_get_keyword_kpath(kpoint_path, stdout, seedname)
    else
      ok = .false.
    end if
    kpoint_path%num_points_first_segment = 100
    call w90_readwrite_get_keyword(stdout, seedname, 'bands_num_points', found, &
                                   i_value=kpoint_path%num_points_first_segment)
    ! checks
    if (bands_plot) then
      if (kpoint_path%num_points_first_segment < 0) then
        call set_error_input(error, 'Error: bands_num_points must be positive')
        return
      endif
    endif
  end subroutine w90_readwrite_read_kpath

  subroutine w90_readwrite_read_fermi_energy(found_fermi_energy, fermi_energy_list, stdout, &
                                             seedname, error)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none
    logical, intent(out) :: found_fermi_energy
    real(kind=dp), allocatable, intent(out) :: fermi_energy_list(:)
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    real(kind=dp) :: fermi_energy
    logical :: fermi_energy_scan
    real(kind=dp) :: fermi_energy_min
    real(kind=dp) :: fermi_energy_max
    real(kind=dp) :: fermi_energy_step
    integer :: i, ierr, n
    logical :: found

    n = 0
    found_fermi_energy = .false.
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy', found, r_value=fermi_energy)
    if (found) then
      found_fermi_energy = .true.
      n = 1
    endif

    fermi_energy_scan = .false.
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy_min', found, r_value=fermi_energy_min)
    if (found) then
      if (found_fermi_energy) then
        call set_error_input(error, 'Error: Cannot specify both fermi_energy and fermi_energy_min')
        return
      endif
      fermi_energy_scan = .true.
      fermi_energy_max = fermi_energy_min + 1.0_dp
      call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy_max', found, &
                                     r_value=fermi_energy_max)
      if (found .and. fermi_energy_max <= fermi_energy_min) then
        call set_error_input(error, 'Error: fermi_energy_max must be larger than fermi_energy_min')
        return
      endif
      fermi_energy_step = 0.01_dp
      call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy_step', found, &
                                     r_value=fermi_energy_step)
      if (found .and. fermi_energy_step <= 0.0_dp) then
        call set_error_input(error, 'Error: fermi_energy_step must be positive')
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
!!    elseif(nfermi==0) then
!!        ! This happens when both found_fermi_energy=.false. and
!!        ! fermi_energy_scan=.false. Functionalities that require
!!        ! specifying a Fermi level should give an error message
!!        allocate(fermi_energy_list(1),stat=ierr) ! helps streamline things
!!
!! AAM_2017-03-27: if nfermi is zero (ie, fermi_energy* parameters are not set in input file)
!! then allocate fermi_energy_list with length 1 and set to zero as default.
    else
      if (allocated(fermi_energy_list)) deallocate (fermi_energy_list)
      allocate (fermi_energy_list(1), stat=ierr)
      fermi_energy_list(1) = 0.0_dp
    endif
    if (ierr /= 0) then
      call set_error_alloc(error, &
                           'Error allocating fermi_energy_list in w90_readwrite_read_fermi_energy')
      return
    endif
  end subroutine w90_readwrite_read_fermi_energy

  subroutine w90_readwrite_read_ws_data(ws_region, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    type(ws_region_type), intent(inout) :: ws_region
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: i
    logical :: found

    ws_region%use_ws_distance = .true.
    call w90_readwrite_get_keyword(stdout, seedname, 'use_ws_distance', found, &
                                   l_value=ws_region%use_ws_distance)

    ws_region%ws_distance_tol = 1.e-5_dp
    call w90_readwrite_get_keyword(stdout, seedname, 'ws_distance_tol', found, &
                                   r_value=ws_region%ws_distance_tol)

    ws_region%ws_search_size = 2

    call w90_readwrite_get_vector_length(stdout, seedname, 'ws_search_size', found, length=i)
    if (found) then
      if (i .eq. 1) then
        call w90_readwrite_get_keyword_vector(stdout, seedname, 'ws_search_size', found, 1, &
                                              i_value=ws_region%ws_search_size)
        ws_region%ws_search_size(2) = ws_region%ws_search_size(1)
        ws_region%ws_search_size(3) = ws_region%ws_search_size(1)
      elseif (i .eq. 3) then
        call w90_readwrite_get_keyword_vector(stdout, seedname, 'ws_search_size', found, 3, &
                                              i_value=ws_region%ws_search_size)
      else
        call set_error_input(error, &
                             'Error: ws_search_size must be provided as either one integer or a vector of three integers')
        return
      end if
      if (any(ws_region%ws_search_size <= 0)) then
        call set_error_input(error, 'Error: ws_search_size elements must be greater than zero')
        return
      endif
    end if
  end subroutine w90_readwrite_read_ws_data

  subroutine w90_readwrite_read_eigvals(pw90_effective_model, pw90_boltzwann, pw90_geninterp, &
                                        w90_plot, disentanglement, eig_found, eigval, library, &
                                        postproc_setup, num_bands, num_kpts, stdout, seedname, &
                                        error)

    use w90_io, only: io_file_unit
    use w90_error, only: w90_error_type, set_error_file, set_error_open, set_error_alloc

    implicit none
    integer, intent(in) :: num_bands, num_kpts
    integer, intent(in) :: stdout
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    character(len=50), intent(in)  :: seedname
    logical, intent(in) :: disentanglement, library, postproc_setup
    logical, intent(in) :: pw90_effective_model, pw90_boltzwann, pw90_geninterp, w90_plot
    logical, intent(out) :: eig_found
    type(w90_error_type), allocatable, intent(out) :: error
    ! local
    integer :: i, j, k, n, eig_unit, ierr

    ! Read the eigenvalues from wannier.eig
    eig_found = .false.
    if (.not. library .and. .not. pw90_effective_model) then

      if (.not. postproc_setup) then
        inquire (file=trim(seedname)//'.eig', exist=eig_found)
        if (.not. eig_found) then
          if (disentanglement) then
            call set_error_open(error, 'No '//trim(seedname)//'.eig file found. Needed for disentanglement')
            return
          else if ((w90_plot .or. pw90_boltzwann .or. pw90_geninterp)) then
            call set_error_open(error, 'No '//trim(seedname)//'.eig file found. Needed for interpolation')
            return
          end if
        else
          ! Allocate only here
          allocate (eigval(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating eigval in w90_wannier90_readwrite_read')
            return
          endif

          eig_unit = io_file_unit()
          open (unit=eig_unit, file=trim(seedname)//'.eig', form='formatted', status='old', err=105)
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
                call set_error_file(error, 'w90_wannier90_readwrite_read: mismatch in '//trim(seedname)//'.eig')
                return
              end if
            enddo
          end do
          close (eig_unit)
        end if
      end if
    end if

    if (library .and. allocated(eigval)) eig_found = .true.

    return

105 call set_error_open(error, 'Error: Problem opening eigenvalue file '//trim(seedname)//'.eig')
    return
106 call set_error_file(error, 'Error: Problem reading eigenvalue file '//trim(seedname)//'.eig')
    return

  end subroutine w90_readwrite_read_eigvals

  subroutine w90_readwrite_read_dis_manifold(eig_found, dis_manifold, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    logical, intent(in) :: eig_found
    type(dis_manifold_type), intent(inout) :: dis_manifold
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error
    ! local
    logical :: found, found2

    call w90_readwrite_get_keyword(stdout, seedname, 'dis_win_min', found, &
                                   r_value=dis_manifold%win_min)

    call w90_readwrite_get_keyword(stdout, seedname, 'dis_win_max', found, &
                                   r_value=dis_manifold%win_max)
    if (eig_found .and. (dis_manifold%win_max .lt. dis_manifold%win_min)) then
      call set_error_input(error, 'Error: w90_readwrite_read_dis_manifold: check disentanglement windows')
      return
    endif

    dis_manifold%froz_min = -1.0_dp; dis_manifold%froz_max = 0.0_dp
    ! no default for dis_froz_max
    dis_manifold%frozen_states = .false.
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_froz_max', found, &
                                   r_value=dis_manifold%froz_max)
    if (found) then
      dis_manifold%frozen_states = .true.
      dis_manifold%froz_min = dis_manifold%win_min ! default value for the bottom of frozen window
    end if
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_froz_min', found2, &
                                   r_value=dis_manifold%froz_min)
    if (eig_found) then
      if (dis_manifold%froz_max .lt. dis_manifold%froz_min) then
        call set_error_input(error, 'Error: w90_readwrite_read_dis_manifold: check disentanglement frozen windows')
        return
      endif
      if (found2 .and. .not. found) then
        call set_error_input(error, 'Error: w90_readwrite_read_dis_manifold: found dis_froz_min but not dis_froz_max')
        return
      endif
    endif
    ! ndimwin/lwindow are not read
  end subroutine w90_readwrite_read_dis_manifold

  subroutine w90_readwrite_read_kmesh_data(kmesh_input, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc
    implicit none
    type(kmesh_input_type), intent(out) :: kmesh_input
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: itmp, ierr
    logical :: found

    kmesh_input%search_shells = 36
    call w90_readwrite_get_keyword(stdout, seedname, 'search_shells', found, &
                                   i_value=kmesh_input%search_shells)
    if (kmesh_input%search_shells < 0) then
      call set_error_input(error, 'Error: search_shells must be positive')
      return
    endif

    kmesh_input%tol = 0.000001_dp
    call w90_readwrite_get_keyword(stdout, seedname, 'kmesh_tol', found, r_value=kmesh_input%tol)
    if (kmesh_input%tol < 0.0_dp) then
      call set_error_input(error, 'Error: kmesh_tol must be positive')
      return
    endif

    kmesh_input%num_shells = 0
    call w90_readwrite_get_range_vector(stdout, seedname, 'shell_list', found, &
                                        kmesh_input%num_shells, lcount=.true.)
    if (found) then
      if (kmesh_input%num_shells < 0 .or. kmesh_input%num_shells > max_shells) then
        call set_error_input(error, 'Error: number of shell in shell_list must be between zero and six')
        return
      endif
      if (allocated(kmesh_input%shell_list)) deallocate (kmesh_input%shell_list)
      allocate (kmesh_input%shell_list(kmesh_input%num_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating shell_list in w90_wannier90_readwrite_read')
        return
      endif
      call w90_readwrite_get_range_vector(stdout, seedname, 'shell_list', found, &
                                          kmesh_input%num_shells, .false., kmesh_input%shell_list)
      if (any(kmesh_input%shell_list < 1)) then
        call set_error_input(error, 'Error: shell_list must contain positive numbers')
        return
      endif
    else
      if (allocated(kmesh_input%shell_list)) deallocate (kmesh_input%shell_list)
      allocate (kmesh_input%shell_list(max_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating shell_list in w90_readwrite_read_kmesh_data')
        return
      endif
    end if

    call w90_readwrite_get_keyword(stdout, seedname, 'num_shells', found, i_value=itmp)
    if (found .and. (itmp /= kmesh_input%num_shells)) then
      call set_error_input(error, 'Error: Found obsolete keyword num_shells. Its value does not agree with shell_list')
      return
    endif

    ! If .true., does not perform the check of B1 of
    ! Marzari, Vanderbild, PRB 56, 12847 (1997)
    ! in kmesh.F90
    ! mainly needed for the interaction with Z2PACK
    ! By default: .false. (perform the tests)
    kmesh_input%skip_B1_tests = .false.
    call w90_readwrite_get_keyword(stdout, seedname, 'skip_b1_tests', found, &
                                   l_value=kmesh_input%skip_B1_tests)

  end subroutine w90_readwrite_read_kmesh_data

  subroutine w90_readwrite_read_kpoints(pw90_effective_model, library, kpt_latt, num_kpts, &
                                        bohr, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input, set_error_alloc, set_error_dealloc
    implicit none

    character(len=50), intent(in)  :: seedname
    integer, intent(in) :: num_kpts
    integer, intent(in) :: stdout
    logical, intent(in) :: pw90_effective_model, library
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: bohr
    type(w90_error_type), allocatable, intent(out) :: error

    real(kind=dp), allocatable :: kpt_cart(:, :)
    integer :: ierr
    logical :: found

    if (.not. pw90_effective_model) allocate (kpt_cart(3, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating kpt_cart in w90_readwrite_read_kpoints')
      return
    endif
    if (.not. library) then
      allocate (kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating kpt_latt in w90_readwrite_read_kpoints')
        return
      endif
    end if

    call w90_readwrite_get_keyword_block(stdout, seedname, 'kpoints', found, num_kpts, 3, bohr, &
                                         r_value=kpt_cart)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <kpoints> in input file'
    if (.not. library .and. .not. pw90_effective_model) then
      kpt_latt = kpt_cart
      if (.not. found) then
        call set_error_input(error, 'Error: Did not find the kpoint information in the input file')
        return
      endif
    end if

    ! Calculate the kpoints in cartesian coordinates
    !if (.not. pw90_effective_model) then
    !  do nkp = 1, num_kpts
    !    k_points%kpt_cart(:, nkp) = matmul(k_points%kpt_latt(:, nkp), recip_lattice(:, :))
    !  end do
    !endif
    deallocate (kpt_cart, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating kpt_cart in w90_readwrite_read_kpoints')
      return
    endif

  end subroutine w90_readwrite_read_kpoints

  subroutine w90_readwrite_read_lattice(library, real_lattice, bohr, stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    logical, intent(in) :: library
    integer, intent(in) :: stdout
    real(kind=dp), intent(out) :: real_lattice(3, 3)
    real(kind=dp) :: real_lattice_tmp(3, 3)
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    logical :: found

    call w90_readwrite_get_keyword_block(stdout, seedname, 'unit_cell_cart', found, 3, 3, bohr, &
                                         r_value=real_lattice_tmp)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <unit_cell_cart> in input file'
    if (.not. library) then
      real_lattice = transpose(real_lattice_tmp)
      if (.not. found) then
        call set_error_input(error, 'Error: Did not find the cell information in the input file')
        return
      endif
    end if
  end subroutine w90_readwrite_read_lattice

  subroutine w90_readwrite_read_atoms(library, atom_data, real_lattice, bohr, stdout, seedname, &
                                      error)
    use w90_error, only: w90_error_type, set_error_input
    implicit none
    logical, intent(in) :: library
    integer, intent(in) :: stdout
    type(atom_data_type), intent(inout) :: atom_data
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: i_temp, i_temp2
    logical :: found, found2, lunits

    ! Atoms
    if (.not. library) atom_data%num_atoms = 0
    call w90_readwrite_get_block_length(stdout, seedname, 'atoms_frac', found, i_temp, library)
    if (found .and. library) write (stdout, '(a)') ' Ignoring <atoms_frac> in input file'
    call w90_readwrite_get_block_length(stdout, seedname, 'atoms_cart', found2, i_temp2, library, &
                                        lunits)
    if (found2 .and. library) write (stdout, '(a)') ' Ignoring <atoms_cart> in input file'
    if (.not. library) then
      if (found .and. found2) then
        call set_error_input(error, 'Error: Cannot specify both atoms_frac and atoms_cart')
        return
      endif
      if (found .and. i_temp > 0) then
        lunits = .false.
        atom_data%num_atoms = i_temp
      elseif (found2 .and. i_temp2 > 0) then
        atom_data%num_atoms = i_temp2
        if (lunits) atom_data%num_atoms = atom_data%num_atoms - 1
      end if
      if (atom_data%num_atoms > 0) then
        call readwrite_get_atoms(atom_data, library, lunits, real_lattice, bohr, stdout, seedname)
      end if
    endif
  end subroutine w90_readwrite_read_atoms

  subroutine w90_readwrite_clear_keywords(stdout, seedname)
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

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in) :: seedname

    logical :: found

    ! keywords for wannier.x
    call w90_readwrite_get_keyword_block(stdout, seedname, 'dis_spheres', found, 0, 0, 0.0_dp)
    call w90_readwrite_get_keyword_block(stdout, seedname, 'kpoints', found, 0, 0, 0.0_dp)
    call w90_readwrite_get_keyword_block(stdout, seedname, 'nnkpts', found, 0, 0, 0.0_dp)
    call w90_readwrite_get_keyword_block(stdout, seedname, 'unit_cell_cart', found, 0, 0, 0.0_dp)
    call clear_block(stdout, seedname, 'projections')
    call clear_block(stdout, seedname, 'kpoint_path')
    call w90_readwrite_get_keyword(stdout, seedname, 'auto_projections', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'bands_num_points', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'bands_plot_dim', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'bands_plot_format', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'bands_plot', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'bands_plot_mode', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'calc_only_A', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'conv_noise_amp', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'conv_noise_num', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'conv_tol', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'conv_window', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'cp_pp', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'devel_flag', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_conv_tol', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_conv_window', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_froz_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_froz_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_mix_ratio', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_num_iter', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_spheres_first_wann', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_spheres_num', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dist_cutoff', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dist_cutoff_hc', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dist_cutoff_mode', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_win_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dis_win_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'energy_unit', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_energy_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_surface_num_points', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_surface_plot_format', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fermi_surface_plot', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'fixed_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gamma_only', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'guiding_centres', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'hr_cutoff', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'hr_plot', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'iprint', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kmesh_spacing', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kmesh_tol', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'length_unit', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_bands', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_cg_steps', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_dump_cycles', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_elec_per_state', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_guide_cycles', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_iter', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_no_guide_iter', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_print_cycles', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_shells', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_valence_bands', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'num_wann', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'one_dim_axis', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'optimisation', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'postproc_setup', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'precond', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'restart', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'search_shells', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'site_symmetry', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'skip_b1_tests', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'slwf_constrain', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'slwf_lambda', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'slwf_num', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spin', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spinors', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'symmetrize_eps', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'timing_level', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_easy_fix', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_energy_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_group_threshold', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_bandc', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_bb', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_cc', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_cell_ll', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_cell_rr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_cr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_lc', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_ll', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_num_rr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_read_ht', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'translate_home_cell', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'transport', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'transport_mode', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_use_same_lead', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_win_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_win_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'tran_write_ht', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'trial_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'use_bloch_phases', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'use_ws_distance', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wannier_plot_format', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wannier_plot', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wannier_plot_mode', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wannier_plot_radius', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wannier_plot_scale', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wannier_plot_spinor_mode', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wannier_plot_spinor_phase', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_bvec', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_hr_diag', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_hr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_proj', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_r2mn', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_rmn', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_tb', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_u_matrices', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_vdw_data', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'write_xyz', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'ws_distance_tol', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wvfn_formatted', found)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'kmesh', found, 0) ! the absent arrays have zero length ;-)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'mp_grid', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'translation_centre_frac', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'wannier_plot_supercell', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'ws_search_size', found, 0)
    ! ends list of wannier.x keywords

    ! keywords for postw90.x
    call w90_readwrite_get_keyword(stdout, seedname, 'adpt_smr_fac', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'adpt_smr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'adpt_smr_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'berry_curv_adpt_kmesh', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'berry_curv_adpt_kmesh_thresh', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'berry_curv_unit', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'berry', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'berry_kmesh_spacing', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'berry_task', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_2d_dir', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_bandshift_energyshift', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_bandshift_firstband', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_bandshift', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_calc_also_dos', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr_fac', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_adpt_smr_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_energy_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_energy_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_energy_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_smr_fixed_en_width', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_dos_smr_type', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_kmesh_spacing', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_mu_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_mu_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_mu_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_relax_time', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_tdf_energy_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_tdf_smr_fixed_en_width', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_tdf_smr_type', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_temp_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_temp_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltz_temp_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'boltzwann', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'degen_thr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_adpt_smr_fac', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_adpt_smr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_adpt_smr_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_energy_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_energy_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_energy_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_kmesh_spacing', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_smr_fixed_en_width', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_smr_type', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'dos_task', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'effective_model', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'geninterp_alsofirstder', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'geninterp', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'geninterp_single_file', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_degen_thresh', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_eigval_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_freq_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_freq_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_freq_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_kmesh_spacing', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_smr_fixed_en_width', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_smr_max_arg', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_smr_type', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'gyrotropic_task', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kpath_bands_colour', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kpath', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kpath_num_points', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kpath_task', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kslice_fermi_lines_colour', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kslice', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kslice_task', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kdotp_num_bands', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_adpt_smr_fac', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_adpt_smr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_adpt_smr_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_eigval_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_freq_max', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_freq_min', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_freq_step', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_smr_fixed_en_width', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'kubo_smr_type', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'sc_eta', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'scissors_shift', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'sc_phase_conv', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'sc_use_eta_corr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'sc_w_thr', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_alpha', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_bandshift_energyshift', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_bandshift_firstband', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_bandshift', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_beta', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_freq_scan', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_gamma', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'shc_method', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'smr_fixed_en_width', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'smr_max_arg', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'smr_type', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spin_axis_azimuth', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spin_axis_polar', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spin_decomp', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spin_kmesh_spacing', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spin_moment', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'spn_formatted', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'transl_inv', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'uhu_formatted', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'use_degen_pert', found)
    call w90_readwrite_get_keyword(stdout, seedname, 'wanint_kpoint_file', found)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'berry_kmesh', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'boltz_kmesh', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'dos_kmesh', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'gyrotropic_box_b1', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'gyrotropic_box_b2', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'gyrotropic_box_b3', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'gyrotropic_box_center', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'gyrotropic_kmesh', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'kdotp_kpoint', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'kslice_2dkmesh', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'kslice_b1', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'kslice_b2', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'kslice_corner', found, 0)
    call w90_readwrite_get_keyword_vector(stdout, seedname, 'spin_kmesh', found, 0)
    ! BGS what about get_range_vectors and gyrotropic_band_list, kdotp_bands etc?
    ! ends list of postw90 keywords

  end subroutine w90_readwrite_clear_keywords

  subroutine w90_readwrite_clean_infile(stdout, seedname, error)
    use w90_error, only: w90_error_type, set_error_input, set_error_dealloc
    implicit none
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: loop, ierr

    ! filter out any remaining accepted keywords from both wannier90.x and postw90.x sets
    call w90_readwrite_clear_keywords(stdout, seedname)

    if (any(len_trim(in_data(:)) > 0)) then
      write (stdout, '(1x,a)') 'The following section of file '//trim(seedname)//'.win contained unrecognised keywords'
      write (stdout, *)
      do loop = 1, num_lines
        if (len_trim(in_data(loop)) > 0) then
          write (stdout, '(1x,a)') trim(in_data(loop))
        end if
      end do
      write (stdout, *)
      call set_error_input(error, 'Unrecognised keyword(s) in input file, see also output file')
      return
    end if

    deallocate (in_data, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating in_data in w90_readwrite_clean_infile')
      return
    endif

  end subroutine w90_readwrite_clean_infile

  subroutine w90_readwrite_read_final_alloc(disentanglement, dis_manifold, wannier_data, &
                                            num_wann, num_bands, num_kpts, error)
    !================================================== !
    ! Some checks and initialisations !
    !================================================== !
    use w90_error, only: w90_error_type, set_error_alloc
    implicit none
    logical, intent(in) :: disentanglement
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(wannier_data_type), intent(inout) :: wannier_data
    integer, intent(in) :: num_wann, num_bands, num_kpts
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: ierr

    if (disentanglement) then
      if (allocated(dis_manifold%ndimwin)) deallocate (dis_manifold%ndimwin)
      allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating ndimwin in w90_wannier90_readwrite_read')
        return
      endif
      if (allocated(dis_manifold%lwindow)) deallocate (dis_manifold%lwindow)
      allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating lwindow in w90_wannier90_readwrite_read')
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
      call set_error_alloc(error, 'Error allocating wannier_centres in w90_wannier90_readwrite_read')
      return
    endif
    wannier_data%centres = 0.0_dp
    if (allocated(wannier_data%spreads)) deallocate (wannier_data%spreads)
    allocate (wannier_data%spreads(num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wannier_spreads in w90_wannier90_readwrite_read')
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

  function w90_readwrite_get_smearing_index(string, keyword, stdout, seedname)
    !! This function parses a string containing the type of
    !! smearing and returns the correct index for the smearing_index variable
    !
    !! If the string is not valid, an io_error is issued
    use w90_io, only: io_error
    integer, intent(in) :: stdout
    character(len=*), intent(in) :: string
    !! The string read from input
    character(len=*), intent(in) :: keyword
    !! The keyword that was read (e.g., smr_type), so that we can print a more useful error message
    character(len=50), intent(in)  :: seedname
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
        if (w90_readwrite_get_smearing_index < 0) &
          call io_error('Wrong m-p smearing order in keyword '//trim(keyword), stdout, seedname)
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
      call io_error('Unrecognised value for keyword '//trim(keyword), stdout, seedname)
    end if

    return

337 call io_error('Wrong m-p smearing order in keyword '//trim(keyword), stdout, seedname)

  end function w90_readwrite_get_smearing_index

!================================================
  subroutine w90_readwrite_uppercase(atom_data, kpoint_path, length_unit)
    !================================================
    !! Convert a few things to uppercase to look nice in the output
    !
    !================================================

    implicit none

    type(atom_data_type), intent(inout) :: atom_data
    type(kpoint_path_type), intent(inout) :: kpoint_path
    character(len=*), intent(inout) :: length_unit
    integer :: nsp, ic, loop, inner_loop

    ! Atom labels (eg, si --> Si)
    do nsp = 1, atom_data%num_species
      ic = ichar(atom_data%label(nsp) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atom_data%label(nsp) (1:1) = char(ic + ichar('Z') - ichar('z'))
    enddo

    do nsp = 1, atom_data%num_species
      ic = ichar(atom_data%symbol(nsp) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atom_data%symbol(nsp) (1:1) = char(ic + ichar('Z') - ichar('z'))
    enddo

    ! Bands labels (eg, x --> X)
    if (allocated(kpoint_path%labels)) then
      do loop = 1, size(kpoint_path%labels)
        do inner_loop = 1, len(kpoint_path%labels(loop))
          ic = ichar(kpoint_path%labels(loop) (inner_loop:inner_loop))
          if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
            kpoint_path%labels(loop) (inner_loop:inner_loop) = char(ic + ichar('Z') - ichar('z'))
        enddo
      enddo
    endif

    ! Length unit (ang --> Ang, bohr --> Bohr)
    ic = ichar(length_unit(1:1))
    if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
      length_unit(1:1) = char(ic + ichar('Z') - ichar('z'))

    return

  end subroutine w90_readwrite_uppercase

  subroutine w90_readwrite_write_header(bohr_version_str, constants_version_str1, &
                                        constants_version_str2, stdout)
    !! Write a suitable header for the calculation - version authors etc
    use w90_io, only: io_date, w90_version

    implicit none

    integer, intent(in) :: stdout
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

  end subroutine w90_readwrite_write_header

!================================================!
  subroutine w90_readwrite_dealloc(exclude_bands, wannier_data, input_proj, kmesh_input, kpt_latt, &
                                   dis_manifold, atom_data, eigval, kpoint_path, error)
    !================================================!
    !! release memory from allocated parameters
    !
    !================================================
    use w90_error, only: w90_error_type, set_error_dealloc

    implicit none

    type(atom_data_type), intent(inout) :: atom_data
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(kpoint_path_type), intent(inout) :: kpoint_path
    type(proj_input_type), intent(inout) :: input_proj
    type(wannier_data_type), intent(inout) :: wannier_data
    integer, allocatable, intent(inout) :: exclude_bands(:)
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    type(w90_error_type), allocatable, intent(out) :: error

    integer :: ierr

    if (allocated(dis_manifold%ndimwin)) then
      deallocate (dis_manifold%ndimwin, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ndimwin in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(dis_manifold%lwindow)) then
      deallocate (dis_manifold%lwindow, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating lwindow in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(eigval)) then
      deallocate (eigval, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating eigval in w90_readwrite_dealloc')
        return
      endif
    endif
    if (allocated(kmesh_input%shell_list)) then
      deallocate (kmesh_input%shell_list, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating shell_list in w90_readwrite_dealloc')
        return
      endif
    endif
    if (allocated(kpt_latt)) then
      deallocate (kpt_latt, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating kpt_latt in w90_readwrite_dealloc')
        return
      endif
    endif
    if (allocated(kpoint_path%labels)) then
      deallocate (kpoint_path%labels, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating labels in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(kpoint_path%points)) then
      deallocate (kpoint_path%points, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating points in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(atom_data%label)) then
      deallocate (atom_data%label, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_label in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(atom_data%symbol)) then
      deallocate (atom_data%symbol, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_symbol in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(atom_data%pos_cart)) then
      deallocate (atom_data%pos_cart, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_pos_cart in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(atom_data%species_num)) then
      deallocate (atom_data%species_num, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating atoms_species_num in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%site)) then
      deallocate (input_proj%site, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_site in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%l)) then
      deallocate (input_proj%l, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_l in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%m)) then
      deallocate (input_proj%m, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_m in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%s)) then
      deallocate (input_proj%s, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_s in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%s_qaxis)) then
      deallocate (input_proj%s_qaxis, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_s_qaxis in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%z)) then
      deallocate (input_proj%z, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_z in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%x)) then
      deallocate (input_proj%x, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_x in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%radial)) then
      deallocate (input_proj%radial, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_radial in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(input_proj%zona)) then
      deallocate (input_proj%zona, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating input_proj_zona in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(exclude_bands)) then
      deallocate (exclude_bands, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating exclude_bands in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(wannier_data%centres)) then
      deallocate (wannier_data%centres, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wannier_centres in w90_readwrite_dealloc')
        return
      endif
    end if
    if (allocated(wannier_data%spreads)) then
      deallocate (wannier_data%spreads, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wannier_spreads in w90_readwrite_dealloc')
        return
      endif
    endif
    return

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
!~    use w90_io,        only : io_file_unit,io_error,seedname,io_date
!~    implicit none
!~
!~    integer :: i,j,k,l,um_unit
!~    character (len=9) :: cdate, ctime
!~    character(len=33) :: header
!~
!~    call io_date(cdate, ctime)
!~    header='written on '//cdate//' at '//ctime
!~
!~    um_unit=io_file_unit()
!~    open(unit=um_unit,file=trim(seedname)//'_um.dat',form='unformatted')
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
!~    use w90_io,        only : io_file_unit,io_error,seedname
!~    implicit none
!~
!~    integer       :: tmp_num_wann,tmp_num_kpts,tmp_num_nnmax
!~    integer       :: i,j,k,l,um_unit,ierr
!~    character(len=33) :: header
!~    real(kind=dp) :: tmp_omi
!~
!~    um_unit=io_file_unit()
!~    open(unit=um_unit,file=trim(seedname)//'_um.dat',status="old",form='unformatted',err=105)
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
                                      have_disentangled, ispostw90, seedname, stdout, error)
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
    use w90_io, only: io_file_unit
    use w90_error, only: w90_error_type, set_error_file, set_error_open, set_error_alloc
    use w90_utility, only: utility_recip_lattice

    implicit none

    integer, allocatable, intent(inout) :: exclude_bands(:)
    type(wannier_data_type), intent(inout) :: wannier_data
    type(kmesh_info_type), intent(in) :: kmesh_info
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_exclude_bands

    complex(kind=dp), allocatable, intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix(:, :, :, :)

    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: omega_invariant

    character(len=50), intent(in)  :: seedname
    character(len=*), intent(inout) :: checkpoint

    logical, intent(in) :: ispostw90 ! Are we running postw90?
    logical, intent(out) :: have_disentangled

    ! local variables
    real(kind=dp) :: recip_lattice(3, 3), volume
    integer :: chk_unit, nkp, i, j, k, l, ntmp, ierr
    character(len=33) :: header
    real(kind=dp) :: tmp_latt(3, 3), tmp_kpt_latt(3, num_kpts)
    integer :: tmp_excl_bands(1:num_exclude_bands), tmp_mp_grid(1:3)

    write (stdout, '(1x,3a)') 'Reading restart information from file ', trim(seedname), '.chk :'

    chk_unit = io_file_unit()
    open (unit=chk_unit, file=trim(seedname)//'.chk', status='old', form='unformatted', err=121)

    ! Read comment line
    read (chk_unit) header
    write (stdout, '(1x,a)', advance='no') trim(header)

    ! Consistency checks
    read (chk_unit) ntmp                           ! Number of bands
    if (ntmp .ne. num_bands) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_bands')
      return
    endif
    read (chk_unit) ntmp                           ! Number of excluded bands
    if (ntmp .ne. num_exclude_bands) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_exclude_bands')
      return
    endif
    read (chk_unit) (tmp_excl_bands(i), i=1, num_exclude_bands) ! Excluded bands
    do i = 1, num_exclude_bands
      if (tmp_excl_bands(i) .ne. exclude_bands(i)) then
        call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in exclude_bands')
        return
      endif
    enddo
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Real lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - real_lattice(i, j)) .gt. eps6) then
          call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in real_lattice')
          return
        endif
      enddo
    enddo
    call utility_recip_lattice(real_lattice, recip_lattice, volume, stdout, seedname)
    read (chk_unit) ((tmp_latt(i, j), i=1, 3), j=1, 3)  ! Reciprocal lattice
    do j = 1, 3
      do i = 1, 3
        if (abs(tmp_latt(i, j) - recip_lattice(i, j)) .gt. eps6) then
          call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in recip_lattice')
          return
        endif
      enddo
    enddo
    read (chk_unit) ntmp                ! K-points
    if (ntmp .ne. num_kpts) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_kpts')
      return
    endif
    read (chk_unit) (tmp_mp_grid(i), i=1, 3)         ! M-P grid
    do i = 1, 3
      if (tmp_mp_grid(i) .ne. mp_grid(i)) then
        call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in mp_grid')
        return
      endif
    enddo
    read (chk_unit) ((tmp_kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts)
    do nkp = 1, num_kpts
      do i = 1, 3
        if (abs(tmp_kpt_latt(i, nkp) - kpt_latt(i, nkp)) .gt. eps6) then
          call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in kpt_latt')
          return
        endif
      enddo
    enddo
    read (chk_unit) ntmp                ! nntot
    if (ntmp .ne. kmesh_info%nntot) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in nntot')
      return
    endif
    read (chk_unit) ntmp                ! num_wann
    if (ntmp .ne. num_wann) then
      call set_error_file(error, 'w90_readwrite_read_chk: Mismatch in num_wann')
      return
    endif
    ! End of consistency checks

    read (chk_unit) checkpoint             ! checkpoint
    checkpoint = adjustl(trim(checkpoint))

    read (chk_unit) have_disentangled      ! whether a disentanglement has been performed

    if (have_disentangled) then

      read (chk_unit) omega_invariant     ! omega invariant

      ! lwindow
      if (.not. allocated(dis_manifold%lwindow)) then
        allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating lwindow in w90_readwrite_read_chkpt')
          return
        endif
      endif
      read (chk_unit, err=122) ((dis_manifold%lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)

      ! ndimwin
      if (.not. allocated(dis_manifold%ndimwin)) then
        allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating ndimwin in w90_readwrite_read_chkpt')
          return
        endif
      endif
      read (chk_unit, err=123) (dis_manifold%ndimwin(nkp), nkp=1, num_kpts)

      ! U_matrix_opt
      if (.not. allocated(u_matrix_opt)) then
        allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error allocating u_matrix_opt in w90_readwrite_read_chkpt')
          return
        endif
      endif
      read (chk_unit, err=124) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)

    endif

    ! U_matrix
    if (.not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating u_matrix in w90_readwrite_read_chkpt')
        return
      endif
    endif
    read (chk_unit, err=125) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)

    ! M_matrix
    if (.not. allocated(m_matrix)) then
      allocate (m_matrix(num_wann, num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating m_matrix in w90_readwrite_read_chkpt')
        return
      endif
    endif
    read (chk_unit, err=126) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, kmesh_info%nntot), l=1, num_kpts)

    ! wannier_centres
    read (chk_unit, err=127) ((wannier_data%centres(i, j), i=1, 3), j=1, num_wann)

    ! wannier spreads
    read (chk_unit, err=128) (wannier_data%spreads(i), i=1, num_wann)

    close (chk_unit)

    write (stdout, '(a/)') ' ... done'

    return

121 if (ispostw90) then
      call set_error_open(error, &
                          'Error opening '//trim(seedname)//'.chk in w90_readwrite_read_chkpt: did you run wannier90.x first?')
    else
      call set_error_open(error, 'Error opening '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    end if
    return
122 call set_error_file(error, 'Error reading lwindow from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    return
123 call set_error_file(error, 'Error reading ndimwin from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    return
124 call set_error_file(error, 'Error reading u_matrix_opt from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    return
125 call set_error_file(error, 'Error reading u_matrix from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    return
126 call set_error_file(error, 'Error reading m_matrix from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    return
127 call set_error_file(error, 'Error reading wannier_centres from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    return
128 call set_error_file(error, 'Error reading wannier_spreads from '//trim(seedname)//'.chk in w90_readwrite_read_chkpt')
    return

  end subroutine w90_readwrite_read_chkpt

!================================================!
  subroutine w90_readwrite_chkpt_dist(dis_manifold, wannier_data, u_matrix, u_matrix_opt, &
                                      omega_invariant, num_bands, num_kpts, num_wann, checkpoint, &
                                      have_disentangled, seedname, stdout, error, comm)
    !================================================!
    !
    !! Distribute the chk files
    !
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_file_unit, io_date, io_time
    use w90_comms, only: comms_bcast, w90comm_type, mpirank
    use w90_error, only: w90_error_type, set_error_alloc

    implicit none

    ! arguments
    type(wannier_data_type), intent(inout) :: wannier_data
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90comm_type), intent(in) :: comm

    integer, intent(in) :: stdout
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_wann
    integer, intent(inout) :: num_kpts

    complex(kind=dp), allocatable, intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :)
    real(kind=dp), intent(inout) :: omega_invariant

    character(len=50), intent(in)  :: seedname
    character(len=*), intent(inout) :: checkpoint
    logical, intent(inout) :: have_disentangled

    ! local variables
    integer :: ierr

    logical :: on_root = .false.

    if (mpirank(comm) == 0) on_root = .true.

    call comms_bcast(checkpoint, len(checkpoint), stdout, seedname, comm)

    if (.not. on_root .and. .not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating u_matrix in w90_readwrite_chkpt_dist')
        return
      endif
    endif
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, stdout, seedname, comm)

!    if (.not.on_root .and. .not.allocated(m_matrix)) then
!       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!       if (ierr/=0)&
!            call io_error('Error allocating m_matrix in w90_readwrite_chkpt_dist')
!    endif
!    call comms_bcast(m_matrix(1,1,1,1),num_wann*num_wann*nntot*num_kpts)

    call comms_bcast(have_disentangled, 1, stdout, seedname, comm)

    if (have_disentangled) then
      if (.not. on_root) then

        if (.not. allocated(u_matrix_opt)) then
          allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating u_matrix_opt in w90_readwrite_chkpt_dist')
            return
          endif
        endif

        if (.not. allocated(dis_manifold%lwindow)) then
          allocate (dis_manifold%lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating lwindow in w90_readwrite_chkpt_dist')
            return
          endif
        endif

        if (.not. allocated(dis_manifold%ndimwin)) then
          allocate (dis_manifold%ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error allocating ndimwin in w90_readwrite_chkpt_dist')
            return
          endif
        endif

      end if

      call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts, stdout, seedname, comm)
      call comms_bcast(dis_manifold%lwindow(1, 1), num_bands*num_kpts, stdout, seedname, comm)
      call comms_bcast(dis_manifold%ndimwin(1), num_kpts, stdout, seedname, comm)
      call comms_bcast(omega_invariant, 1, stdout, seedname, comm)
    end if
    call comms_bcast(wannier_data%centres(1, 1), 3*num_wann, stdout, seedname, comm)
    call comms_bcast(wannier_data%spreads(1), num_wann, stdout, seedname, comm)

  end subroutine w90_readwrite_chkpt_dist

!================================================!
  subroutine w90_readwrite_in_file(seedname, error)
    !================================================!
    !! Load the *.win file into a character
    !! array in_file, ignoring comments and
    !! blank lines and converting everything
    !! to lowercase characters
    !================================================!

    use w90_utility, only: utility_lowercase
    use w90_io, only: io_file_unit
    use w90_error, only: w90_error_type, set_error_alloc, set_error_open, set_error_file

    implicit none

    character(len=50), intent(in)  :: seedname
    type(w90_error_type), allocatable, intent(out) :: error

    integer           :: in_unit, tot_num_lines, ierr, line_counter, loop, in1, in2
    character(len=maxlen) :: dummy
    integer           :: pos
    character, parameter :: TABCHAR = char(9)

    in_unit = io_file_unit()
    open (in_unit, file=trim(seedname)//'.win', form='formatted', status='old', err=101)

    num_lines = 0; tot_num_lines = 0
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
        if (len(trim(dummy)) > 0) num_lines = num_lines + 1
      endif

    end do

101 call set_error_open(error, 'Error: Problem opening input file '//trim(seedname)//'.win')
    return
200 call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.win')
    return
210 continue
    rewind (in_unit)

    allocate (in_data(num_lines), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating in_data in w90_readwrite_in_file')
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
      if (in1 == 0 .and. in2 == 0) in_data(line_counter) = dummy
      if (in1 == 0 .and. in2 > 0) in_data(line_counter) = dummy(:in2 - 1)
      if (in2 == 0 .and. in1 > 0) in_data(line_counter) = dummy(:in1 - 1)
      if (in2 > 0 .and. in1 > 0) in_data(line_counter) = dummy(:min(in1, in2) - 1)
    end do

    close (in_unit)

  end subroutine w90_readwrite_in_file

  !================================================!
  subroutine w90_readwrite_get_keyword(stdout, seedname, keyword, found, c_value, l_value, i_value, r_value)
    !================================================!
    !
    !! Finds the value of the required keyword.
    !
    !================================================!

    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
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

    integer           :: kl, in, loop, itmp
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      itmp = in + len(trim(keyword))
      if (in_data(loop) (itmp:itmp) /= '=' &
          .and. in_data(loop) (itmp:itmp) /= ':' &
          .and. in_data(loop) (itmp:itmp) /= ' ') cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      in_data(loop) (1:maxlen) = ' '
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
          call io_error('Error: Problem reading logical keyword '//trim(keyword), stdout, seedname)
        endif
      endif
      if (present(i_value)) read (dummy, *, err=220, end=220) i_value
      if (present(r_value)) read (dummy, *, err=220, end=220) r_value
    end if

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword), stdout, seedname)

  end subroutine w90_readwrite_get_keyword

  !================================================!
  subroutine w90_readwrite_get_keyword_vector(stdout, seedname, keyword, found, length, c_value, &
                                              l_value, i_value, r_value)
    !================================================!
    !
    !! Finds the values of the required keyword vector
    !
    !================================================!

    use w90_io, only: io_error

    implicit none

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
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer           :: kl, in, loop, i
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      in_data(loop) (1:maxlen) = ' '
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
        call io_error('w90_readwrite_get_keyword_vector unimplemented for logicals', stdout, seedname)
      endif
      if (present(i_value)) read (dummy, *, err=230, end=230) (i_value(i), i=1, length)
      if (present(r_value)) read (dummy, *, err=230, end=230) (r_value(i), i=1, length)
    end if

    return

230 call io_error('Error: Problem reading keyword '//trim(keyword)//' in w90_readwrite_get_keyword_vector', stdout, seedname)

  end subroutine w90_readwrite_get_keyword_vector

!================================================!
  subroutine w90_readwrite_get_vector_length(stdout, seedname, keyword, found, length)
    !================================================!
    !
    !! Returns the length of a keyword vector
    !
    !================================================!

    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(out)  :: length
    !! length of vector

    integer           :: kl, in, loop, pos
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      dummy = adjustl(dummy)
      if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
        dummy = dummy(2:)
        dummy = adjustl(dummy)
      end if
    end do

    length = 0
    if (found) then
      if (len_trim(dummy) == 0) call io_error('Error: keyword '//trim(keyword)//' is blank', stdout, seedname)
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

    return

  end subroutine w90_readwrite_get_vector_length

  !================================================!
  subroutine w90_readwrite_get_keyword_block(stdout, seedname, keyword, found, rows, columns, &
                                             bohr, c_value, l_value, i_value, r_value)
    !================================================!
    !
    !!   Finds the values of the required data block
    !
    !================================================!

    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
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

    integer           :: in, ins, ine, loop, i, line_e, line_s, counter, blen
    logical           :: found_e, found_s, lconvert
    character(len=maxlen) :: dummy, end_st, start_st

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    if (.not. found_s) then
      found = .false.
      return
    end if

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
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

    if ((blen .ne. rows) .and. (blen .ne. rows + 1) .and. (rows .gt. 0)) &
      call io_error('Error: Wrong number of lines in block '//trim(keyword), stdout, seedname)

    if ((blen .eq. rows + 1) .and. (rows .gt. 0) .and. (index(trim(keyword), 'unit_cell_cart') .eq. 0)) &
      call io_error('Error: Wrong number of lines in block '//trim(keyword), stdout, seedname)

    found = .true.

    lconvert = .false.
    if (blen == rows + 1) then
      dummy = in_data(line_s + 1)
      if (index(dummy, 'ang') .ne. 0) then
        lconvert = .false.
      elseif (index(dummy, 'bohr') .ne. 0) then
        lconvert = .true.
      else
        call io_error('Error: Units in block '//trim(keyword)//' not recognised', stdout, seedname)
      endif
      in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    endif

!    r_value=1.0_dp
    counter = 0
    do loop = line_s + 1, line_e - 1
      dummy = in_data(loop)
      counter = counter + 1
      if (present(c_value)) read (dummy, *, err=240, end=240) (c_value(i, counter), i=1, columns)
      if (present(l_value)) then
        ! I don't think we need this. Maybe read into a dummy charater
        ! array and convert each element to logical
        call io_error('w90_readwrite_get_keyword_block unimplemented for logicals', stdout, seedname)
      endif
      if (present(i_value)) read (dummy, *, err=240, end=240) (i_value(i, counter), i=1, columns)
      if (present(r_value)) read (dummy, *, err=240, end=240) (r_value(i, counter), i=1, columns)
    end do

    if (lconvert) then
      if (present(r_value)) then
        r_value = r_value*bohr
      endif
    endif

    in_data(line_s:line_e) (1:maxlen) = ' '

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword), stdout, seedname)

  end subroutine w90_readwrite_get_keyword_block

  !================================================!
  subroutine w90_readwrite_get_block_length(stdout, seedname, keyword, found, rows, library, lunits)
    !================================================!
    !
    !! Finds the length of the data block
    !
    !================================================!

    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)  :: keyword
    !! Keyword to examine
    logical, intent(out) :: found
    !! Is keyword present
    integer, intent(out) :: rows
    !! Number of rows
    logical, intent(in) :: library
    logical, optional, intent(out) :: lunits
    !! Have we found a unit specification

    integer           :: i, in, ins, ine, loop, line_e, line_s
    logical           :: found_e, found_s
    character(len=maxlen) :: end_st, start_st, dummy
    character(len=2)  :: atsym
    real(kind=dp)     :: atpos(3)

    rows = 0
    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    if (.not. found_s) then
      found = .false.
      return
    end if

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    rows = line_e - line_s - 1

    found = .true.

    ! Ignore atoms_cart and atoms_frac blocks if running in library mode
    if (library) then
      if (trim(keyword) .eq. 'atoms_cart' .or. trim(keyword) .eq. 'atoms_frac') then
        in_data(line_s:line_e) (1:maxlen) = ' '
      endif
    endif

    if (present(lunits)) then
      dummy = in_data(line_s + 1)
      read (dummy, *, end=555) atsym, (atpos(i), i=1, 3)
      lunits = .false.
    endif

    if (rows <= 0) then !cope with empty blocks
      found = .false.
      in_data(line_s:line_e) (1:maxlen) = ' '
    end if

    return

555 lunits = .true.

    if (rows <= 1) then !cope with empty blocks
      found = .false.
      in_data(line_s:line_e) (1:maxlen) = ' '
    end if

    return

  end subroutine w90_readwrite_get_block_length

  !================================================!
  subroutine readwrite_get_atoms(atom_data, library, lunits, real_lattice, bohr, stdout, seedname)
    !================================================!
    !
    !!   Fills the atom data block
    !
    !================================================!

    use w90_utility, only: utility_frac_to_cart, utility_cart_to_frac, utility_inverse_mat
    use w90_io, only: io_error
    implicit none

    type(atom_data_type), intent(inout) :: atom_data
    integer, intent(in) :: stdout
    logical, intent(in) :: library
    logical, intent(in) :: lunits
    !! Do we expect a first line with the units
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: bohr
    character(len=50), intent(in)  :: seedname

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

    keyword = "atoms_cart"
    frac = .false.
    call w90_readwrite_get_block_length(stdout, seedname, "atoms_frac", found, i_temp, library)
    if (found) then
      keyword = "atoms_frac"
      frac = .true.
    end if

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    lconvert = .false.
    if (lunits) then
      dummy = in_data(line_s + 1)
      if (index(dummy, 'ang') .ne. 0) then
        lconvert = .false.
      elseif (index(dummy, 'bohr') .ne. 0) then
        lconvert = .true.
      else
        call io_error('Error: Units in block atoms_cart not recognised in readwrite_get_atoms', stdout, seedname)
      endif
      in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    endif

    counter = 0
    do loop = line_s + 1, line_e - 1
      dummy = in_data(loop)
      counter = counter + 1
      if (frac) then
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_frac_tmp(i, counter), i=1, 3)
      else
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_cart_tmp(i, counter), i=1, 3)
      end if
    end do

    if (lconvert) atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr

    in_data(line_s:line_e) (1:maxlen) = ' '

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
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in readwrite_get_atoms', stdout, seedname)
    allocate (atom_data%label(atom_data%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in readwrite_get_atoms', stdout, seedname)
    allocate (atom_data%symbol(atom_data%num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in readwrite_get_atoms', stdout, seedname)
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
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in readwrite_get_atoms', stdout, seedname)

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

240 call io_error('Error: Problem reading block keyword '//trim(keyword), stdout, seedname)

  end subroutine readwrite_get_atoms

  !================================================!
  subroutine w90_readwrite_lib_set_atoms(atom_data, atoms_label_tmp, atoms_pos_cart_tmp, &
                                         real_lattice, error)
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
    character(len=*), intent(in) :: atoms_label_tmp(atom_data%num_atoms)
    !! Atom labels
    real(kind=dp), intent(in)      :: atoms_pos_cart_tmp(3, atom_data%num_atoms)
    !! Atom positions
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    real(kind=dp)     :: inv_lattice(3, 3)
    real(kind=dp)     :: atoms_pos_frac_tmp(3, atom_data%num_atoms)
    integer           :: loop2, max_sites, ierr, ic, loop, counter
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
      call set_error_alloc(error, 'Error allocating atoms_species_num in w90_readwrite_lib_set_atoms')
      return
    endif
    allocate (atom_data%label(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_label in w90_readwrite_lib_set_atoms')
      return
    endif
    allocate (atom_data%symbol(atom_data%num_species), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating atoms_symbol in w90_readwrite_lib_set_atoms')
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
      call set_error_alloc(error, 'Error allocating atoms_pos_cart in w90_readwrite_lib_set_atoms')
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
    end do

    return

  end subroutine w90_readwrite_lib_set_atoms

  !================================================!
  subroutine w90_readwrite_get_range_vector(stdout, seedname, keyword, found, length, lcount, i_value)
    !================================================!
    !!   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100
    !!   if(lcount) we return the number of states in length
    !================================================!
    use w90_io, only: io_error

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    character(*), intent(in)    :: keyword
    !! Keyword to examine
    logical, intent(out)   :: found
    !! Is keyword found
    integer, intent(inout) :: length
    !! Number of states
    logical, intent(in)    :: lcount
    !! If T only count states
    integer, optional, intent(out)   :: i_value(length)
    !! States specified in range vector

    integer   :: kl, in, loop, num1, num2, i_punc
    integer   :: counter, i_digit, loop_r, range_size
    character(len=maxlen) :: dummy
    character(len=10), parameter :: c_digit = "0123456789"
    character(len=2), parameter :: c_range = "-:"
    character(len=3), parameter :: c_sep = " ,;"
    character(len=5), parameter :: c_punc = " ,;-:"
    character(len=5)  :: c_num1, c_num2

    if (lcount .and. present(i_value)) call io_error('w90_readwrite_get_range_vector: incorrect call', stdout, seedname)

    kl = len_trim(keyword)

    found = .false.

    do loop = 1, num_lines
      in = index(in_data(loop), trim(keyword))
      if (in == 0 .or. in > 1) cycle
      if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file', stdout, seedname)
      endif
      found = .true.
      dummy = in_data(loop) (kl + 1:)
      dummy = adjustl(dummy)
      if (.not. lcount) in_data(loop) (1:maxlen) = ' '
      if (dummy(1:1) == '=' .or. dummy(1:1) == ':') then
        dummy = dummy(2:)
        dummy = adjustl(dummy)
      end if
    end do

    if (.not. found) return

    counter = 0
    if (len_trim(dummy) == 0) call io_error('Error: keyword '//trim(keyword)//' is blank', stdout, seedname)
    dummy = adjustl(dummy)
    do
      i_punc = scan(dummy, c_punc)
      if (i_punc == 0) call io_error('Error parsing keyword '//trim(keyword), stdout, seedname)
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
      if (scan(dummy, c_range) == 1) call io_error('Error parsing keyword '//trim(keyword)//' incorrect range', stdout, seedname)
      if (index(dummy, ' ') == 1) exit
    end do

    if (lcount) length = counter
    if (.not. lcount) then
      do loop = 1, counter - 1
        do loop_r = loop + 1, counter
          if (i_value(loop) == i_value(loop_r)) &
            call io_error('Error parsing keyword '//trim(keyword)//' duplicate values', stdout, seedname)
        end do
      end do
    end if

    return

101 call io_error('Error parsing keyword '//trim(keyword), stdout, seedname)

  end subroutine w90_readwrite_get_range_vector

  subroutine w90_readwrite_get_centre_constraints(ccentres_frac, ccentres_cart, &
                                                  proj_site, num_wann, real_lattice, stdout, seedname)
    !================================================!
    !!  assigns projection centres as default centre constraints and global
    !!  Lagrange multiplier as individual Lagrange multipliers then reads
    !!  the centre_constraints block for individual centre constraint parameters
    !
    !================================================!
    use w90_io, only: io_error
    use w90_utility, only: utility_frac_to_cart
    integer, intent(in) :: stdout
    real(kind=dp), intent(inout) :: ccentres_frac(:, :), ccentres_cart(:, :)
    real(kind=dp), intent(in) :: proj_site(:, :)
    integer, intent(in) :: num_wann
    character(len=50), intent(in)  :: seedname
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    integer           :: loop1, index1, constraint_num, loop2
    integer           :: column, start, finish, wann
    !logical           :: found
    character(len=maxlen) :: dummy

    do loop1 = 1, num_wann
      do loop2 = 1, 3
        ccentres_frac(loop1, loop2) = proj_site(loop2, loop1)
      end do
    end do

    constraint_num = 0
    do loop1 = 1, num_lines
      dummy = in_data(loop1)
      if (constraint_num > 0) then
        if (trim(dummy) == '') cycle
        index1 = index(dummy, 'begin')
        if (index1 > 0) call io_error("slwf_centres block hasn't ended yet", stdout, seedname)
        index1 = index(dummy, 'end')
        if (index1 > 0) then
          index1 = index(dummy, 'slwf_centres')
          if (index1 == 0) call io_error('Wrong ending of block (need to end slwf_centres)', stdout, seedname)
          in_data(loop1) (1:maxlen) = ' '
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
                                                     wann, dummy, ccentres_frac, stdout, seedname)
              start = loop2 + 1
              finish = start
            end if
          end if
          if (loop2 == len_trim(dummy) .and. dummy(loop2:loop2) /= ' ') then
            finish = loop2
            call get_centre_constraint_from_column(column, start, finish, &
                                                   wann, dummy, ccentres_frac, stdout, seedname)
            start = loop2 + 1
            finish = start
          end if
        end do
        in_data(loop1) (1:maxlen) = ' '
        constraint_num = constraint_num + 1
      end if
      index1 = index(dummy, 'slwf_centres')
      if (index1 > 0) then
        index1 = index(dummy, 'begin')
        if (index1 > 0) then
          constraint_num = 1
          in_data(loop1) (1:maxlen) = ' '
        end if
      end if
    end do
    do loop1 = 1, num_wann
      call utility_frac_to_cart(ccentres_frac(loop1, :), &
                                ccentres_cart(loop1, :), real_lattice)
    end do
  end subroutine w90_readwrite_get_centre_constraints

  !================================================!
  subroutine get_centre_constraint_from_column(column, start, finish, &
                                               wann, dummy, ccentres_frac, stdout, seedname)
    !================================================!
    !
    !!  assigns value read to constraint
    !!  parameters based on column
    !
    !================================================!
    use w90_io, only: io_error
    integer, intent(in) :: stdout
    integer, intent(inout):: column, start, finish, wann
    character(len=maxlen), intent(inout):: dummy
    character(len=50), intent(in)  :: seedname
    real(kind=dp), intent(inout) :: ccentres_frac(:, :)

    if (column == 0) then
      read (dummy(start:finish), '(i3)') wann
    end if
    if (column > 0) then
      if (column > 4) call io_error("Didn't expect anything else after Lagrange multiplier", stdout, seedname)
      if (column < 4) read (dummy(start:finish), '(f10.10)') ccentres_frac(wann, column)
    end if
    column = column + 1
  end subroutine get_centre_constraint_from_column

  !================================================!
  subroutine w90_readwrite_get_projections(num_proj, atom_data, num_wann, input_proj, proj, &
                                           inv_lattice, lcount, spinors, bohr, stdout, seedname)
    !================================================!
    !
    !!  Fills the projection data block
    !
    !================================================!

    use w90_constants, only: eps6, eps2
    use w90_utility, only: utility_cart_to_frac, utility_string_to_coord, utility_strip
    use w90_io, only: io_error

    implicit none

    ! arguments
    type(atom_data_type), intent(in) :: atom_data
    type(proj_input_type), intent(inout) :: input_proj
    type(proj_input_type), intent(inout) :: proj ! intent(out)?
    integer, intent(in) :: num_wann
    integer, intent(inout) :: num_proj
    integer, intent(in) :: stdout
    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(in) :: inv_lattice(3, 3)
    character(len=50), intent(in)  :: seedname
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

!     if(spinors) num_proj=num_wann/2

    if (.not. lcount) then
      allocate (input_proj%site(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_site in w90_readwrite_get_projections', stdout, seedname)
      allocate (input_proj%l(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_l in w90_readwrite_get_projections', stdout, seedname)
      allocate (input_proj%m(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_m in w90_readwrite_get_projections', stdout, seedname)
      allocate (input_proj%z(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_z in w90_readwrite_get_projections', stdout, seedname)
      allocate (input_proj%x(3, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_x in w90_readwrite_get_projections', stdout, seedname)
      allocate (input_proj%radial(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_radial in w90_readwrite_get_projections', stdout, seedname)
      allocate (input_proj%zona(num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating input_proj_zona in w90_readwrite_get_projections', stdout, seedname)
      if (spinors) then
        allocate (input_proj%s(num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s in w90_readwrite_get_projections', stdout, seedname)
        allocate (input_proj%s_qaxis(3, num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating input_proj_s_qaxis in w90_readwrite_get_projections', stdout, seedname)
      endif

      allocate (proj%site(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_site in w90_readwrite_get_projections', stdout, seedname)
      allocate (proj%l(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_l in w90_readwrite_get_projections', stdout, seedname)
      allocate (proj%m(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_m in w90_readwrite_get_projections', stdout, seedname)
      allocate (proj%z(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_z in w90_readwrite_get_projections', stdout, seedname)
      allocate (proj%x(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_x in w90_readwrite_get_projections', stdout, seedname)
      allocate (proj%radial(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_radial in w90_readwrite_get_projections', stdout, seedname)
      allocate (proj%zona(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating proj_zona in w90_readwrite_get_projections', stdout, seedname)
      if (spinors) then
        allocate (proj%s(num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_s in w90_readwrite_get_projections', stdout, seedname)
        allocate (proj%s_qaxis(3, num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating proj_s_qaxis in w90_readwrite_get_projections', stdout, seedname)
      endif
    endif

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('w90_readwrite_get_projections: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('w90_readwrite_get_projections: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, &
                    seedname)
    end if

    if (line_e <= line_s) then
      call io_error('w90_readwrite_get_projections: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, &
                    seedname)
    end if

    dummy = in_data(line_s + 1)
    lconvert = .false.
    lrandom = .false.
    lpartrandom = .false.
    if (index(dummy, 'ang') .ne. 0) then
      if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
    elseif (index(dummy, 'bohr') .ne. 0) then
      if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
      lconvert = .true.
    elseif (index(dummy, 'random') .ne. 0) then
      if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
      line_s = line_s + 1
      if (index(in_data(line_s + 1), end_st) .ne. 0) then
        lrandom = .true.     ! all projections random
      else
        lpartrandom = .true. ! only some projections random
        if (index(in_data(line_s + 1), 'ang') .ne. 0) then
          if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
          line_s = line_s + 1
        elseif (index(in_data(line_s + 1), 'bohr') .ne. 0) then
          if (.not. lcount) in_data(line_s) (1:maxlen) = ' '
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
        dummy = utility_strip(in_data(line))
        dummy = adjustl(dummy)
        pos1 = index(dummy, ':')
        if (pos1 == 0) &
          call io_error('w90_wannier90_readwrite_read_projection: malformed projection definition: '//trim(dummy), stdout, &
                        seedname)
        sites = 0
        ctemp = dummy(:pos1 - 1)
        ! Read the atomic site
        if (index(ctemp, 'c=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_cart, stdout, seedname)
          if (lconvert) pos_cart = pos_cart*bohr
          call utility_cart_to_frac(pos_cart(:), pos_frac(:), inv_lattice)
        elseif (index(ctemp, 'f=') > 0) then
          sites = -1
          ctemp = ctemp(3:)
          call utility_string_to_coord(ctemp, pos_frac, stdout, seedname)
        else
          if (atom_data%num_species == 0) &
            call io_error('w90_wannier90_readwrite_read_projection: Atom centred projection requested but no atoms defined', &
                          stdout, seedname)
          do loop = 1, atom_data%num_species
            if (trim(ctemp) == atom_data%label(loop)) then
              species = loop
              sites = atom_data%species_num(loop)
              exit
            end if
            if (loop == atom_data%num_species) then
              call io_error('w90_wannier90_readwrite_read_projection: Atom site not recognised '//trim(ctemp), &
                            stdout, seedname)
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
            if (pos2 == 0) call io_error &
              ('w90_readwrite_get_projections: no closing square bracket for spin quantisation dir', stdout, seedname)
            ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_s_qaxis_tmp, stdout, seedname)
            dummy = dummy(:pos1 - 1) ! remove [ ] section
          endif
        else
          if (pos1 > 0) call io_error('w90_readwrite_get_projections: spin qdir is defined but spinors=.false.', stdout, seedname)
        endif

        ! scan for up or down
        pos1 = index(dummy, '(')
        if (spinors) then
          if (pos1 > 0) then
            proj_u_tmp = .false.; proj_d_tmp = .false.
            ctemp = (dummy(pos1 + 1:))
            pos2 = index(ctemp, ')')
            if (pos2 == 0) call io_error('w90_readwrite_get_projections: no closing bracket for spin', stdout, seedname)
            ctemp = ctemp(:pos2 - 1)
            if (index(ctemp, 'u') > 0) proj_u_tmp = .true.
            if (index(ctemp, 'd') > 0) proj_d_tmp = .true.
            if (proj_u_tmp .and. proj_d_tmp) then
              spn_counter = 2
            elseif (.not. proj_u_tmp .and. .not. proj_d_tmp) then
              call io_error('w90_readwrite_get_projections: found brackets but neither u or d', stdout, seedname)
            else
              spn_counter = 1
            endif
            dummy = dummy(:pos1 - 1) ! remove ( ) section
          endif
        else
          if (pos1 > 0) call io_error('w90_readwrite_get_projections: spin is defined but spinors=.false.', stdout, seedname)
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
            if (l_tmp < -5 .or. l_tmp > 3) call io_error('w90_readwrite_get_projections: Incorrect l state requested', stdout, &
                                                         seedname)
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
              if (index(ctemp3, 'mr=') /= mstate + 1) &
                call io_error('w90_readwrite_get_projections: Problem reading m state', stdout, seedname)
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
                  if ((m_tmp > 2*l_tmp + 1) .or. (m_tmp <= 0)) call io_error('w90_readwrite_get_projections: m is > l !', &
                                                                             stdout, seedname)
                elseif (l_tmp == -1 .and. (m_tmp > 2 .or. m_tmp <= 0)) then
                  call io_error('w90_readwrite_get_projections: m has incorrect value (1)', stdout, seedname)
                elseif (l_tmp == -2 .and. (m_tmp > 3 .or. m_tmp <= 0)) then
                  call io_error('w90_readwrite_get_projections: m has incorrect value (2)', stdout, seedname)
                elseif (l_tmp == -3 .and. (m_tmp > 4 .or. m_tmp <= 0)) then
                  call io_error('w90_readwrite_get_projections: m has incorrect value (3)', stdout, seedname)
                elseif (l_tmp == -4 .and. (m_tmp > 5 .or. m_tmp <= 0)) then
                  call io_error('w90_readwrite_get_projections: m has incorrect value (4)', stdout, seedname)
                elseif (l_tmp == -5 .and. (m_tmp > 6 .or. m_tmp <= 0)) then
                  call io_error('w90_readwrite_get_projections: m has incorrect value (5)', stdout, seedname)
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
                call io_error('w90_readwrite_get_projections: Problem reading l state '//trim(ctemp3), stdout, seedname)
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
            call utility_string_to_coord(ctemp, proj_z_tmp, stdout, seedname)
          endif
          ! x axis
          pos1 = index(dummy, 'x=')
          if (pos1 > 0) then
            ctemp = (dummy(pos1 + 2:))
            pos2 = index(ctemp, ':')
            if (pos2 > 0) ctemp = ctemp(:pos2 - 1)
            call utility_string_to_coord(ctemp, proj_x_tmp, stdout, seedname)
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
                  input_proj%site(:, counter) = pos_frac
                  input_proj%l(counter) = loop_l
                  input_proj%m(counter) = loop_m
                  input_proj%z(:, counter) = proj_z_tmp
                  input_proj%x(:, counter) = proj_x_tmp
                  input_proj%radial(counter) = proj_radial_tmp
                  input_proj%zona(counter) = proj_zona_tmp
                  if (spinors) then
                    if (spn_counter == 1) then
                      if (proj_u_tmp) input_proj%s(counter) = 1
                      if (proj_d_tmp) input_proj%s(counter) = -1
                    else
                      if (loop_s == 1) input_proj%s(counter) = 1
                      if (loop_s == 2) input_proj%s(counter) = -1
                    endif
                    input_proj%s_qaxis(:, counter) = proj_s_qaxis_tmp
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
                    call utility_cart_to_frac(atom_data%pos_cart(:, loop_sites, species), pos_frac, &
                                              inv_lattice)
                    input_proj%site(:, counter) = pos_frac(:)
                    input_proj%l(counter) = loop_l
                    input_proj%m(counter) = loop_m
                    input_proj%z(:, counter) = proj_z_tmp
                    input_proj%x(:, counter) = proj_x_tmp
                    input_proj%radial(counter) = proj_radial_tmp
                    input_proj%zona(counter) = proj_zona_tmp
                    if (spinors) then
                      if (spn_counter == 1) then
                        if (proj_u_tmp) input_proj%s(counter) = 1
                        if (proj_d_tmp) input_proj%s(counter) = -1
                      else
                        if (loop_s == 1) input_proj%s(counter) = 1
                        if (loop_s == 2) input_proj%s(counter) = -1
                      endif
                      input_proj%s_qaxis(:, counter) = proj_s_qaxis_tmp
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
        if (counter .lt. num_wann) call io_error( &
          'w90_readwrite_get_projections: too few projection functions defined', stdout, seedname)
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
      do loop = counter + 1, num_wann
        call random_number(input_proj%site(:, loop))
        input_proj%l(loop) = 0
        input_proj%m(loop) = 1
        input_proj%z(:, loop) = proj_z_def
        input_proj%x(:, loop) = proj_x_def
        input_proj%zona(loop) = proj_zona_def
        input_proj%radial(loop) = proj_radial_def
        if (spinors) then
          if (modulo(loop, 2) == 1) then
            input_proj%s(loop) = 1
          else
            input_proj%s(loop) = -1
          end if
          input_proj%s_qaxis(1, loop) = 0.
          input_proj%s_qaxis(2, loop) = 0.
          input_proj%s_qaxis(3, loop) = 1.
        end if
      enddo
    endif

    ! I shouldn't get here, but just in case
    if (.not. lcount) in_data(line_s:line_e) (1:maxlen) = ' '

!~     ! Check
!~     do loop=1,num_wann
!~        if ( abs(sum(proj_z(:,loop)*proj_x(:,loop))).gt.1.0e-6_dp ) then
!~           write(stdout,*) ' Projection:',loop
!~           call io_error(' Error in projections: z and x axes are not orthogonal')
!~        endif
!~     enddo

    ! Normalise z-axis and x-axis and check/fix orthogonality
    do loop = 1, num_proj

      znorm = sqrt(sum(input_proj%z(:, loop)*input_proj%z(:, loop)))
      xnorm = sqrt(sum(input_proj%x(:, loop)*input_proj%x(:, loop)))
      input_proj%z(:, loop) = input_proj%z(:, loop)/znorm             ! normalise z
      input_proj%x(:, loop) = input_proj%x(:, loop)/xnorm             ! normalise x
      cosphi = sum(input_proj%z(:, loop)*input_proj%x(:, loop))

      ! Check whether z-axis and z-axis are orthogonal
      if (abs(cosphi) .gt. eps6) then

        ! Special case of circularly symmetric projections (pz, dz2, fz3)
        ! just choose an x-axis that is perpendicular to the given z-axis
        if ((input_proj%l(loop) .ge. 0) .and. (input_proj%m(loop) .eq. 1)) then
          proj_x_tmp(:) = input_proj%x(:, loop)            ! copy of original x-axis
          call random_seed()
          call random_number(proj_z_tmp(:))         ! random vector
          ! calculate new x-axis as the cross (vector) product of random vector with z-axis
          input_proj%x(1, loop) = proj_z_tmp(2)*input_proj%z(3, loop) &
                                  - proj_z_tmp(3)*input_proj%z(2, loop)
          input_proj%x(2, loop) = proj_z_tmp(3)*input_proj%z(1, loop) &
                                  - proj_z_tmp(1)*input_proj%z(3, loop)
          input_proj%x(3, loop) = proj_z_tmp(1)*input_proj%z(2, loop) &
                                  - proj_z_tmp(2)*input_proj%z(1, loop)
          xnorm_new = sqrt(sum(input_proj%x(:, loop)*input_proj%x(:, loop)))
          input_proj%x(:, loop) = input_proj%x(:, loop)/xnorm_new   ! normalise
          goto 555
        endif

        ! If projection axes non-orthogonal enough, then
        ! user may have made a mistake and should check
        if (abs(cosphi) .gt. eps2) then
          write (stdout, *) ' Projection:', loop
          call io_error(' Error in projections: z and x axes are not orthogonal', stdout, seedname)
        endif

        ! If projection axes are "reasonably orthogonal", project x-axis
        ! onto plane perpendicular to z-axis to make them more so
        sinphi = sqrt(1 - cosphi*cosphi)
        proj_x_tmp(:) = input_proj%x(:, loop)               ! copy of original x-axis
        ! calculate new x-axis:
        ! x = z \cross (x_tmp \cross z) / sinphi = ( x_tmp - z(z.x_tmp) ) / sinphi
        input_proj%x(:, loop) = (proj_x_tmp(:) - cosphi*input_proj%z(:, loop))/sinphi

        ! Final check
555     cosphi_new = sum(input_proj%z(:, loop)*input_proj%x(:, loop))
        if (abs(cosphi_new) .gt. eps6) then
          write (stdout, *) ' Projection:', loop
          call io_error(' Error: z and x axes are still not orthogonal after projection', stdout, seedname)
        endif

      endif

    enddo

    return

101 call io_error('w90_readwrite_get_projections: Problem reading l state into integer '//trim(ctemp3), stdout, seedname)
102 call io_error('w90_readwrite_get_projections: Problem reading m state into integer '//trim(ctemp3), stdout, seedname)
104 call io_error('w90_readwrite_get_projections: Problem reading zona into real '//trim(ctemp), stdout, seedname)
105 call io_error('w90_readwrite_get_projections: Problem reading radial state into integer '//trim(ctemp), stdout, seedname)
106 call io_error('w90_readwrite_get_projections: Problem reading m state into string '//trim(ctemp3), stdout, seedname)

  end subroutine w90_readwrite_get_projections

  !================================================!
  subroutine w90_readwrite_get_keyword_kpath(kpoint_path, stdout, seedname)
    !================================================!
    !
    !!  Fills the kpath data block
    !
    !================================================!
    use w90_io, only: io_error

    implicit none

    type(kpoint_path_type), intent(inout) :: kpoint_path
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    character(len=20) :: keyword
    integer           :: in, ins, ine, loop, i, line_e, line_s, counter
    logical           :: found_e, found_s
    character(len=maxlen) :: dummy, end_st, start_st

    keyword = "kpoint_path"

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (line_e <= line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
    end if

    counter = 0
    do loop = line_s + 1, line_e - 1

      counter = counter + 2
      dummy = in_data(loop)
      read (dummy, *, err=240, end=240) kpoint_path%labels(counter - 1), &
        (kpoint_path%points(i, counter - 1), i=1, 3), &
        kpoint_path%labels(counter), (kpoint_path%points(i, counter), i=1, 3)
    end do

    in_data(line_s:line_e) (1:maxlen) = ' '

    return

240 call io_error('w90_readwrite_get_keyword_kpath: Problem reading kpath '//trim(dummy), stdout, seedname)

  end subroutine w90_readwrite_get_keyword_kpath

  !================================================!
  subroutine clear_block(stdout, seedname, keyword)
    !================================================!
    ! a dummy read routine to remove unused but legitimate input block from input stream
    ! needed to preserve input file error checking (i.e. input stream should be empty after all
    ! legitimate keywords/blocks are read)
    !================================================!
    use w90_io, only: io_error

    implicit none

    ! arguments
    integer, intent(in) :: stdout
    character(len=50), intent(in) :: seedname
    character(len=*), intent(in) :: keyword

    ! local variables
    integer :: in, ins, ine, loop, line_e, line_s
    logical :: found_e, found_s
    character(len=maxlen) :: end_st, start_st

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), 'begin')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file', stdout, seedname)
      endif
      found_s = .true.
    end do

    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), 'end')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file', stdout, seedname)
      endif
      found_e = .true.
    end do

    if (found_s .and. (.not. found_e)) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file', stdout, seedname)
    end if

    if (found_e .and. (.not. found_s)) then
      call io_error('Error: Found '//trim(end_st)//' but no '//trim(start_st)//' in input file', stdout, seedname)
    end if

    if (found_s .and. found_e) then
      if (line_e <= line_s) then
        call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file', stdout, seedname)
      end if

      in_data(line_s:line_e) (1:maxlen) = ' '  ! clear the block from the input stream
    end if ! found tags
  end subroutine clear_block

end module w90_readwrite
