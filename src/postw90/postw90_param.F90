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

module pw90_parameters

  use w90_constants, only: dp

  implicit none

  public

  type pw90_calculation_type !postw90.F90
    logical :: kpath
    logical :: kslice
    logical :: dos
    logical :: berry
    logical :: gyrotropic
    logical :: geninterp
    logical :: boltzwann
    !BGS spin_moment?
  end type pw90_calculation_type
  type(pw90_calculation_type), save :: pw90_calcs

  logical, save :: eig_found ! used to control broadcast of eigval

  type postw90_oper_type ! only in postw90/get_oper.F90
    logical :: spn_formatted
    !! Read the spin from fortran formatted file
    logical :: uHu_formatted
    !! Read the uHu from fortran formatted file
  end type postw90_oper_type
  type(postw90_oper_type), save :: postw90_oper

  type postw90_common_type
    logical :: spin_moment !postw90_common and postw90
    logical :: spin_decomp !postw90_common, berry, dos and boltzwann
    real(kind=dp) :: scissors_shift ! get_oper and berry
    ! IVO
    ! Are we running postw90 starting from an effective model?
    logical :: effective_model = .false.
  end type postw90_common_type
  type(postw90_common_type), save :: pw90_common

! Module  s p i n
  type postw90_spin_type
    !Todo - remove spin_ ?
    real(kind=dp) :: spin_axis_polar
    real(kind=dp) :: spin_axis_azimuth
    real(kind=dp) :: spin_kmesh_spacing
    integer :: spin_kmesh(3)
  end type postw90_spin_type
  type(postw90_spin_type), save :: pw90_spin

  type postw90_ham_type
    logical :: use_degen_pert
    real(kind=dp) :: degen_thr
  end type postw90_ham_type
  type(postw90_ham_type), save :: pw90_ham

  ! module  k p a t h (used by postw90/kpath)
  type kpath_type
    character(len=20) :: task
    integer :: num_points
    character(len=20) :: bands_colour
  end type kpath_type
  type(kpath_type), save :: kpath

  ! module  k s l i c e (postw90/kslice)
  type kslice_type
    character(len=20) :: task
    real(kind=dp) :: corner(3)
    real(kind=dp) :: b1(3)
    real(kind=dp) :: b2(3)
    integer :: kmesh2d(2)
    character(len=20) :: fermi_lines_colour
  end type kslice_type
  type(kslice_type), save :: kslice

  ! module  d o s
  ! No need to save 'dos_plot', only used here (introduced 'dos_task')
  logical          :: dos_plot

  !BGS a generic smr_index/fixed_en_width etc type for here, boltzwann etc?
  type dos_plot_type
    character(len=20)    :: task
    logical    :: adpt_smr
    real(kind=dp)    :: adpt_smr_fac
    integer    :: smr_index
    real(kind=dp)    :: smr_fixed_en_width
    real(kind=dp)    :: adpt_smr_max
    real(kind=dp)    :: energy_max
    real(kind=dp)    :: energy_min
    real(kind=dp)    :: energy_step
    integer    :: num_project
    integer, allocatable :: project(:)
    !  character(len=20)    :: plot_format
    real(kind=dp)    :: kmesh_spacing
    integer    :: kmesh(3)
    !  real(kind=dp) :: gaussian_width
  end type dos_plot_type
  type(dos_plot_type), save :: dos_data

  ! Module  b e r r y (mainly postw90/berry)
  type berry_type
    character(len=120) :: task
    real(kind=dp) :: kmesh_spacing
    integer :: kmesh(3)
    integer :: curv_adpt_kmesh
    real(kind=dp) :: curv_adpt_kmesh_thresh
    character(len=20) :: curv_unit ! postw90/kpath, kslice as well
    logical :: kubo_adpt_smr ! also postw90/kpath
    real(kind=dp) :: kubo_adpt_smr_fac
    integer :: kubo_smr_index
    real(kind=dp) :: kubo_smr_fixed_en_width
    real(kind=dp) :: kubo_adpt_smr_max
    integer :: sc_phase_conv
    real(kind=dp) :: sc_eta ! also postw90/wan_ham
    real(kind=dp) :: sc_w_thr
    logical :: wanint_kpoint_file ! also postw90/spin, postw90/dos, postw90.F90
    logical :: transl_inv !also used in postw90/get_oper, postw90/gyrotropic
    integer :: kubo_nfreq
    complex(kind=dp), allocatable :: kubo_freq_list(:)
    real(kind=dp) :: kubo_eigval_max
  end type berry_type
  type(berry_type), save :: berry

  ! spin Hall conductivity (postw90 - common, get_oper, berry, kpath)
  type spin_hall_type
    logical :: freq_scan
    integer :: alpha
    integer :: beta
    integer :: gamma
    logical :: bandshift
    integer :: bandshift_firstband
    real(kind=dp) :: bandshift_energyshift
  end type spin_hall_type
  type(spin_hall_type), save :: spin_hall

  type gyrotropic_type ! postw90 - common, gyrotropic
    character(len=120) :: task
    integer :: kmesh(3)
    real(kind=dp) :: kmesh_spacing
    integer :: smr_index
    real(kind=dp) :: smr_fixed_en_width
    integer :: nfreq
    complex(kind=dp), allocatable :: freq_list(:)
    real(kind=dp) :: box_corner(3), box(3, 3)
    real(kind=dp) :: degen_thresh
    integer, allocatable :: band_list(:)
    integer :: num_bands
    real(kind=dp) :: smr_max_arg
    real(kind=dp) :: eigval_max
  end type gyrotropic_type
  type(gyrotropic_type), save :: gyrotropic

  ! [gp-begin, Jun 1, 2012]
  ! GeneralInterpolator variables - postw90/geninterp
  type geninterp_type
    logical :: alsofirstder
    logical :: single_file
  end type geninterp_type
  type(geninterp_type), save :: geninterp
  ! [gp-end, Jun 1, 2012]

  ! [gp-begin, Apr 12, 2012]
  ! BoltzWann variables (postw90/boltzwann.F90)
  type boltzwann_type
    logical :: calc_also_dos
    integer :: dir_num_2d
    real(kind=dp) :: dos_energy_step
    real(kind=dp) :: dos_energy_min
    real(kind=dp) :: dos_energy_max
    logical :: dos_adpt_smr
    real(kind=dp) :: dos_smr_fixed_en_width
    real(kind=dp) :: dos_adpt_smr_fac
    real(kind=dp) :: dos_adpt_smr_max
    real(kind=dp) :: mu_min
    real(kind=dp) :: mu_max
    real(kind=dp) :: mu_step
    real(kind=dp) :: temp_min
    real(kind=dp) :: temp_max
    real(kind=dp) :: temp_step
    real(kind=dp) :: kmesh_spacing
    integer :: kmesh(3)
    real(kind=dp) :: tdf_energy_step
    integer :: TDF_smr_index
    integer :: dos_smr_index
    real(kind=dp) :: relax_time
    real(kind=dp) :: TDF_smr_fixed_en_width
    logical :: bandshift
    integer :: bandshift_firstband
    real(kind=dp) :: bandshift_energyshift
  end type boltzwann_type
  type(boltzwann_type), save :: boltz
  ! [gp-end, Apr 12, 2012]

end module pw90_parameters

module pw90_param_methods

  use w90_constants, only: dp
  use w90_io, only: stdout, maxlen
  use w90_parameters, only: param_input, dis_data, num_wann, spec_points, &
    eigval, fermi, atoms, real_lattice, recip_lattice, num_bands
  use w90_param_methods !BGS Todo only:
  use pw90_parameters

  implicit none

  private

  ! [gp-begin, Apr 13, 2012]
  ! Global interpolation k mesh variables
  ! These don't need to be public, since their values are copied in the variables of the
  ! local interpolation meshes. JRY: added save attribute
  real(kind=dp), save             :: kmesh_spacing
  integer, save                   :: kmesh(3)
  logical, save                   :: global_kmesh_set
  ! [gp-end]

  ! Adaptive vs. fixed smearing stuff [GP, Jul 12, 2012]
  ! Only internal, always use the local variables defined by each module
  ! that take this value as default
  logical                         :: adpt_smr
  real(kind=dp)                   :: adpt_smr_fac
  real(kind=dp)                   :: adpt_smr_max
  real(kind=dp)                   :: smr_fixed_en_width

  integer :: smr_index

  ! from gyrotropic section
  real(kind=dp)                               :: gyrotropic_freq_min
  real(kind=dp)                               :: gyrotropic_freq_max
  real(kind=dp)                               :: gyrotropic_freq_step
  real(kind=dp)                   :: gyrotropic_box_tmp(3)
  real(kind=dp)                   :: smr_max_arg

  real(kind=dp)                               :: kubo_freq_min
  real(kind=dp)                               :: kubo_freq_max
  real(kind=dp)                               :: kubo_freq_step

  character(len=4), save :: boltz_2d_dir

  ! private data
  character(len=maxlen)              :: ctmp

  public :: param_postw90_read
  public :: param_postw90_write

contains

  subroutine param_postw90_read
    !==================================================================!
    !                                                                  !
    !! Read parameters and calculate derived values
    !!
    !! Note on parallelization: this function should be called
    !! from the root node only!
    !!
    !                                                                  !
    !===================================================================
    !use w90_constants, only: bohr, eps6, cmplx_i
    !use w90_utility, only: utility_recip_lattice
    !use w90_io, only: io_error, io_file_unit, seedname, post_proc_flag
    implicit none

    !local variables
    !real(kind=dp)  :: real_lattice_tmp(3, 3)
    !integer :: nkp, i, j, n, k, itmp, i_temp, i_temp2, eig_unit, loop, ierr, iv_temp(3), rows
    !logical :: found, found2, lunits, chk_found
    !character(len=6) :: spin_str
    !real(kind=dp) :: rv_temp(3)
    !integer, allocatable, dimension(:, :) :: nnkpts_block
    !integer, allocatable, dimension(:) :: nnkpts_idx
    !real(kind=dp) :: cell_volume
    logical                                  :: found_fermi_energy
    logical :: disentanglement, library
    character(len=20) :: energy_unit

    library = .false.
    call param_in_file
    !call param_w90_read_01
    !call param_w90_read_02
    call param_read_03
    !if (.not.(w90_calcs%transport .and. tran%read_ht)) then
    call param_pw90_read_04
    call param_read_05(energy_unit)
    !call param_w90_read_06
    call param_pw90_read_07
    !call param_w90_read_08
    call param_read_09
    !call param_w90_read_10
    call param_read_11(pw90_common%effective_model, library)
    !call param_w90_read_12
    call param_read_13(pw90_common%effective_model, library)
    !call param_w90_read_14
    !call param_w90_read_15
    call param_read_16(library)
    !call param_w90_read_17
    !call param_w90_read_18
    !call param_w90_read_19
    !call param_w90_read_20
    call param_read_21(.false., library)
    !call param_w90_read_22
    call param_read_23(found_fermi_energy)
    call param_pw90_read_24
    call param_read_25(smr_index)
    call param_pw90_read_26(found_fermi_energy)
    !call param_w90_read_27
    !endif
    call param_read_28
    !call param_w90_read_29
    !if (.not.(w90_calcs%transport .and. tran%read_ht)) then
    !call param_w90_read_30
    call param_pw90_read_31
    disentanglement = (num_bands > num_wann)
    call param_read_32(pw90_common%effective_model, pw90_calcs%boltzwann, &
                       pw90_calcs%geninterp, dos_plot, disentanglement, &
                       eig_found, library, .false.)
    call param_w90_read_33(eig_found) !(dis_data?)
    ! Need to make sure use w90_params are read
    call param_pw90_read_34(smr_index)
    !call param_w90_read_35
    call param_pw90_read_36
    !call param_w90_read_37(... one_dim_axis)
    call param_pw90_read_38
    !call param_w90_read_39
    call param_read_40a(pw90_common%effective_model, library)
    call param_read_40c(global_kmesh_set, kmesh_spacing, kmesh)
    call param_pw90_read_40
    !call param_w90_read_41
    !call param_w90_read_42(.false.)
    !call param_w90_read_43
    !endif
    call param_read_44(.false.)
    !if (.not.(w90_calcs%transport .and. tran%read_ht)) then
    !  call param_w90_read_45(disentanglement)
    !endif
  end subroutine param_postw90_read

  subroutine param_pw90_read_04
    implicit none
    logical :: found

    !ivo
    call param_get_keyword('effective_model', found, l_value=pw90_common%effective_model)
  end subroutine param_pw90_read_04

  subroutine param_pw90_read_07
    implicit none
    logical :: found

    postw90_oper%spn_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword('spn_formatted', found, l_value=postw90_oper%spn_formatted)

    postw90_oper%uHu_formatted = .false.       ! formatted or "binary" file
    call param_get_keyword('uhu_formatted', found, l_value=postw90_oper%uHu_formatted)
  end subroutine param_pw90_read_07

  subroutine param_pw90_read_24
    use w90_io, only: io_error
    implicit none
    integer :: i
    logical :: found

    pw90_calcs%kslice = .false.
    call param_get_keyword('kslice', found, l_value=pw90_calcs%kslice)

    kslice%task = 'fermi_lines'
    call param_get_keyword('kslice_task', found, c_value=kslice%task)
    if (pw90_calcs%kslice .and. index(kslice%task, 'fermi_lines') == 0 .and. &
        index(kslice%task, 'curv') == 0 .and. &
        index(kslice%task, 'morb') == 0 .and. &
        index(kslice%task, 'shc') == 0) call io_error &
      ('Error: value of kslice_task not recognised in param_read')
    if (pw90_calcs%kslice .and. index(kslice%task, 'curv') > 0 .and. &
        index(kslice%task, 'morb') > 0) call io_error &
      ("Error: kslice_task cannot include both 'curv' and 'morb'")
    if (pw90_calcs%kslice .and. index(kslice%task, 'shc') > 0 .and. &
        index(kslice%task, 'morb') > 0) call io_error &
      ("Error: kslice_task cannot include both 'shc' and 'morb'")
    if (pw90_calcs%kslice .and. index(kslice%task, 'shc') > 0 .and. &
        index(kslice%task, 'curv') > 0) call io_error &
      ("Error: kslice_task cannot include both 'shc' and 'curv'")

    kslice%kmesh2d(1:2) = 50
    call param_get_vector_length('kslice_2dkmesh', found, length=i)
    if (found) then
      if (i == 1) then
        call param_get_keyword_vector('kslice_2dkmesh', found, 1, &
                                      i_value=kslice%kmesh2d)
        kslice%kmesh2d(2) = kslice%kmesh2d(1)
      elseif (i == 2) then
        call param_get_keyword_vector('kslice_2dkmesh', found, 2, &
                                      i_value=kslice%kmesh2d)
      else
        call io_error('Error: kslice_2dkmesh must be provided as either' &
                      //' one integer or a vector of two integers')
      endif
      if (any(kslice%kmesh2d <= 0)) &
        call io_error('Error: kslice_2dkmesh elements must be' &
                      //' greater than zero')
    endif

    kslice%corner = 0.0_dp
    call param_get_keyword_vector('kslice_corner', found, 3, r_value=kslice%corner)

    kslice%b1(1) = 1.0_dp
    kslice%b1(2) = 0.0_dp
    kslice%b1(3) = 0.0_dp
    call param_get_keyword_vector('kslice_b1', found, 3, r_value=kslice%b1)

    kslice%b2(1) = 0.0_dp
    kslice%b2(2) = 1.0_dp
    kslice%b2(3) = 0.0_dp
    call param_get_keyword_vector('kslice_b2', found, 3, r_value=kslice%b2)

    kslice%fermi_lines_colour = 'none'
    call param_get_keyword('kslice_fermi_lines_colour', found, &
                           c_value=kslice%fermi_lines_colour)
    if (pw90_calcs%kslice .and. index(kslice%fermi_lines_colour, 'none') == 0 .and. &
        index(kslice%fermi_lines_colour, 'spin') == 0) call io_error &
      ('Error: value of kslice_fermi_lines_colour not recognised ' &
       //'in param_read')

!    slice_plot_format         = 'plotmv'
!    call param_get_keyword('slice_plot_format',found,c_value=slice_plot_format)
  end subroutine param_pw90_read_24

  subroutine param_pw90_read_26(found_fermi_energy)
    !IVO
    use w90_io, only: io_error
    implicit none
    logical, intent(in) :: found_fermi_energy
    integer :: i, ierr, loop
    logical :: found

    pw90_calcs%dos = .false.
    call param_get_keyword('dos', found, l_value=pw90_calcs%dos)

    pw90_calcs%berry = .false.
    call param_get_keyword('berry', found, l_value=pw90_calcs%berry)

    berry%transl_inv = .false.
    call param_get_keyword('transl_inv', found, l_value=berry%transl_inv)

    berry%task = ' '
    call param_get_keyword('berry_task', found, c_value=berry%task)
    if (pw90_calcs%berry .and. .not. found) call io_error &
      ('Error: berry=T and berry_task is not set')
    if (pw90_calcs%berry .and. index(berry%task, 'ahc') == 0 .and. index(berry%task, 'morb') == 0 &
        .and. index(berry%task, 'kubo') == 0 .and. index(berry%task, 'sc') == 0 &
        .and. index(berry%task, 'shc') == 0) call io_error &
      ('Error: value of berry_task not recognised in param_read')

    ! Stepan
    pw90_calcs%gyrotropic = .false.
    call param_get_keyword('gyrotropic', found, l_value=pw90_calcs%gyrotropic)
    gyrotropic%task = 'all'
    call param_get_keyword('gyrotropic_task', found, c_value=gyrotropic%task)
    gyrotropic%box(:, :) = 0.0
    gyrotropic%degen_thresh = 0.0_dp
    call param_get_keyword('gyrotropic_degen_thresh', found, r_value=gyrotropic%degen_thresh)

    do i = 1, 3
      gyrotropic%box(i, i) = 1.0_dp
      gyrotropic_box_tmp(:) = 0.0_dp
      call param_get_keyword_vector('gyrotropic_box_b'//achar(48 + i), found, 3, r_value=gyrotropic_box_tmp)
      if (found) gyrotropic%box(i, :) = gyrotropic_box_tmp(:)
    enddo
    gyrotropic%box_corner(:) = 0.0_dp
    call param_get_keyword_vector('gyrotropic_box_center', found, 3, r_value=gyrotropic_box_tmp)
    if (found) gyrotropic%box_corner(:) = &
      gyrotropic_box_tmp(:) - 0.5*(gyrotropic%box(1, :) + gyrotropic%box(2, :) + gyrotropic%box(3, :))

    call param_get_range_vector('gyrotropic_band_list', found, gyrotropic%num_bands, lcount=.true.)
    if (found) then
      if (gyrotropic%num_bands < 1) call io_error('Error: problem reading gyrotropic_band_list')
      if (allocated(gyrotropic%band_list)) deallocate (gyrotropic%band_list)
      allocate (gyrotropic%band_list(gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating gyrotropic_band_list in param_read')
      call param_get_range_vector('gyrotropic_band_list', found, gyrotropic%num_bands, .false., gyrotropic%band_list)
      if (any(gyrotropic%band_list < 1) .or. any(gyrotropic%band_list > num_wann)) &
        call io_error('Error: gyrotropic_band_list asks for a non-valid bands')
    else
      ! include all bands in the calculation
      gyrotropic%num_bands = num_wann
      if (allocated(gyrotropic%band_list)) deallocate (gyrotropic%band_list)
      allocate (gyrotropic%band_list(gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating gyrotropic_band_list in param_read')
      do loop = 1, num_wann
        gyrotropic%band_list(loop) = loop
      end do
    end if

    smr_max_arg = 5.0
    call param_get_keyword('smr_max_arg', found, r_value=smr_max_arg)
    if (found .and. (smr_max_arg <= 0._dp)) &
      call io_error('Error: smr_max_arg must be greater than zero')

    gyrotropic%smr_max_arg = smr_max_arg
    call param_get_keyword('gyrotropic_smr_max_arg', found, &
                           r_value=gyrotropic%smr_max_arg)
    if (found .and. (gyrotropic%smr_max_arg <= 0._dp)) call io_error &
      ('Error: gyrotropic_smr_max_arg must be greater than zero')

!-------------------------------------------------------
!    alpha=0
!    call param_get_keyword('alpha',found,i_value=alpha)

!    beta=0
!    call param_get_keyword('beta',found,i_value=beta)

!    gamma=0
!    call param_get_keyword('gamma',found,i_value=gamma)
!-------------------------------------------------------

    berry%curv_adpt_kmesh = 1
    call param_get_keyword('berry_curv_adpt_kmesh', found, &
                           i_value=berry%curv_adpt_kmesh)
    if (berry%curv_adpt_kmesh < 1) &
      call io_error( &
      'Error:  berry_curv_adpt_kmesh must be a positive integer')

    berry%curv_adpt_kmesh_thresh = 100.0_dp
    call param_get_keyword('berry_curv_adpt_kmesh_thresh', found, &
                           r_value=berry%curv_adpt_kmesh_thresh)

    berry%curv_unit = 'ang2'
    call param_get_keyword('berry_curv_unit', found, c_value=berry%curv_unit)
    if (berry%curv_unit .ne. 'ang2' .and. berry%curv_unit .ne. 'bohr2') &
      call io_error &
      ('Error: value of berry_curv_unit not recognised in param_read')

    berry%wanint_kpoint_file = .false.
    call param_get_keyword('wanint_kpoint_file', found, &
                           l_value=berry%wanint_kpoint_file)

!    smear_temp = -1.0_dp
!    call param_get_keyword('smear_temp',found,r_value=smear_temp)

    berry%kubo_adpt_smr = adpt_smr
    call param_get_keyword('kubo_adpt_smr', found, l_value=berry%kubo_adpt_smr)

    berry%kubo_adpt_smr_fac = adpt_smr_fac
    call param_get_keyword('kubo_adpt_smr_fac', found, &
                           r_value=berry%kubo_adpt_smr_fac)
    if (found .and. (berry%kubo_adpt_smr_fac <= 0._dp)) call io_error &
      ('Error: kubo_adpt_smr_fac must be greater than zero')

    berry%kubo_adpt_smr_max = adpt_smr_max
    call param_get_keyword('kubo_adpt_smr_max', found, &
                           r_value=berry%kubo_adpt_smr_max)
    if (berry%kubo_adpt_smr_max <= 0._dp) call io_error &
      ('Error: kubo_adpt_smr_max must be greater than zero')

    berry%kubo_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('kubo_smr_fixed_en_width', found, &
                           r_value=berry%kubo_smr_fixed_en_width)
    if (found .and. (berry%kubo_smr_fixed_en_width < 0._dp)) call io_error &
      ('Error: kubo_smr_fixed_en_width must be greater than or equal to zero')

    gyrotropic%smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('gyrotropic_smr_fixed_en_width', found, &
                           r_value=gyrotropic%smr_fixed_en_width)
    if (found .and. (gyrotropic%smr_fixed_en_width < 0._dp)) call io_error &
      ('Error: gyrotropic_smr_fixed_en_width must be greater than or equal to zero')

    berry%sc_phase_conv = 1
    call param_get_keyword('sc_phase_conv', found, i_value=berry%sc_phase_conv)
    if ((berry%sc_phase_conv .ne. 1) .and. ((berry%sc_phase_conv .ne. 2))) &
      call io_error('Error: sc_phase_conv must be either 1 or 2')

    pw90_common%scissors_shift = 0.0_dp
    call param_get_keyword('scissors_shift', found, &
                           r_value=pw90_common%scissors_shift)

    spin_hall%freq_scan = .false.
    call param_get_keyword('shc_freq_scan', found, l_value=spin_hall%freq_scan)

    spin_hall%alpha = 1
    call param_get_keyword('shc_alpha', found, i_value=spin_hall%alpha)
    if (found .and. (spin_hall%alpha < 1 .or. spin_hall%alpha > 3)) call io_error &
      ('Error:  shc_alpha must be 1, 2 or 3')

    spin_hall%beta = 2
    call param_get_keyword('shc_beta', found, i_value=spin_hall%beta)
    if (found .and. (spin_hall%beta < 1 .or. spin_hall%beta > 3)) call io_error &
      ('Error:  shc_beta must be 1, 2 or 3')

    spin_hall%gamma = 3
    call param_get_keyword('shc_gamma', found, i_value=spin_hall%gamma)
    if (found .and. (spin_hall%gamma < 1 .or. spin_hall%gamma > 3)) call io_error &
      ('Error:  shc_gamma must be 1, 2 or 3')

    spin_hall%bandshift = .false.
    call param_get_keyword('shc_bandshift', found, l_value=spin_hall%bandshift)
    spin_hall%bandshift = spin_hall%bandshift .and. pw90_calcs%berry .and. .not. (index(berry%task, 'shc') == 0)
    if ((abs(pw90_common%scissors_shift) > 1.0e-7_dp) .and. spin_hall%bandshift) &
      call io_error('Error: shc_bandshift and scissors_shift cannot be used simultaneously')

    spin_hall%bandshift_firstband = 0
    call param_get_keyword('shc_bandshift_firstband', found, i_value=spin_hall%bandshift_firstband)
    if (spin_hall%bandshift .and. (.not. found)) &
      call io_error('Error: shc_bandshift required but no shc_bandshift_firstband provided')
    if ((spin_hall%bandshift_firstband < 1) .and. found) &
      call io_error('Error: shc_bandshift_firstband must >= 1')

    spin_hall%bandshift_energyshift = 0._dp
    call param_get_keyword('shc_bandshift_energyshift', found, r_value=spin_hall%bandshift_energyshift)
    if (spin_hall%bandshift .and. (.not. found)) &
      call io_error('Error: shc_bandshift required but no shc_bandshift_energyshift provided')

    pw90_common%spin_moment = .false.
    call param_get_keyword('spin_moment', found, &
                           l_value=pw90_common%spin_moment)

    pw90_spin%spin_axis_polar = 0.0_dp
    call param_get_keyword('spin_axis_polar', found, &
                           r_value=pw90_spin%spin_axis_polar)

    pw90_spin%spin_axis_azimuth = 0.0_dp
    call param_get_keyword('spin_axis_azimuth', found, &
                           r_value=pw90_spin%spin_axis_azimuth)

    pw90_common%spin_decomp = .false.
    call param_get_keyword('spin_decomp', found, l_value=pw90_common%spin_decomp)

    if (pw90_common%spin_decomp .and. (param_input%num_elec_per_state .ne. 1)) then
      call io_error('spin_decomp can be true only if num_elec_per_state is 1')
    end if

    pw90_ham%use_degen_pert = .false.
    call param_get_keyword('use_degen_pert', found, &
                           l_value=pw90_ham%use_degen_pert)

    pw90_ham%degen_thr = 1.0d-4
    call param_get_keyword('degen_thr', found, r_value=pw90_ham%degen_thr)

    pw90_calcs%kpath = .false.
    call param_get_keyword('kpath', found, l_value=pw90_calcs%kpath)

    kpath%task = 'bands'
    call param_get_keyword('kpath_task', found, c_value=kpath%task)
    if (pw90_calcs%kpath .and. index(kpath%task, 'bands') == 0 .and. &
        index(kpath%task, 'curv') == 0 .and. &
        index(kpath%task, 'morb') == 0 .and. &
        index(kpath%task, 'shc') == 0) call io_error &
      ('Error: value of kpath_task not recognised in param_read')
    if (spec_points%bands_num_spec_points == 0 .and. pw90_calcs%kpath) &
      call io_error('Error: a kpath plot has been requested but there is no kpoint_path block')

    kpath%num_points = 100
    call param_get_keyword('kpath_num_points', found, &
                           i_value=kpath%num_points)
    if (kpath%num_points < 0) &
      call io_error('Error: kpath_num_points must be positive')

    kpath%bands_colour = 'none'
    call param_get_keyword('kpath_bands_colour', found, &
                           c_value=kpath%bands_colour)
    if (pw90_calcs%kpath .and. index(kpath%bands_colour, 'none') == 0 .and. &
        index(kpath%bands_colour, 'spin') == 0 .and. &
        index(kpath%bands_colour, 'shc') == 0) call io_error &
      ('Error: value of kpath_bands_colour not recognised in param_read')
    if (pw90_calcs%kpath .and. index(kpath%task, 'shc') > 0 .and. &
        index(kpath%task, 'spin') > 0) call io_error &
      ("Error: kpath_task cannot include both 'shc' and 'spin'")

    ! set to a negative default value
    param_input%num_valence_bands = -99
    call param_get_keyword('num_valence_bands', found, &
                           i_value=param_input%num_valence_bands)
    if (found .and. (param_input%num_valence_bands .le. 0)) &
      call io_error('Error: num_valence_bands should be greater than zero')
    ! there is a check on this parameter later

    dos_data%task = 'dos_plot'
    if (pw90_calcs%dos) then
      dos_plot = .true.
    else
      dos_plot = .false.
    endif
    call param_get_keyword('dos_task', found, c_value=dos_data%task)
    if (pw90_calcs%dos) then
      if (index(dos_data%task, 'dos_plot') == 0 .and. &
          index(dos_data%task, 'find_fermi_energy') == 0) call io_error &
        ('Error: value of dos_task not recognised in param_read')
      if (index(dos_data%task, 'dos_plot') > 0) dos_plot = .true.
      if (index(dos_data%task, 'find_fermi_energy') > 0 .and. found_fermi_energy) &
        call io_error &
        ('Error: Cannot set "dos_task = find_fermi_energy" and give a value to "fermi_energy"')
    end if

!    sigma_abc_onlyorb=.false.
!    call param_get_keyword('sigma_abc_onlyorb',found,l_value=sigma_abc_onlyorb)

! -------------------------------------------------------------------

    !IVO_END

    dos_data%energy_step = 0.01_dp
    call param_get_keyword('dos_energy_step', found, &
                           r_value=dos_data%energy_step)

    dos_data%adpt_smr = adpt_smr
    call param_get_keyword('dos_adpt_smr', found, &
                           l_value=dos_data%adpt_smr)

    dos_data%adpt_smr_fac = adpt_smr_fac
    call param_get_keyword('dos_adpt_smr_fac', found, &
                           r_value=dos_data%adpt_smr_fac)
    if (found .and. (dos_data%adpt_smr_fac <= 0._dp)) &
      call io_error('Error: dos_adpt_smr_fac must be greater than zero')

    dos_data%adpt_smr_max = adpt_smr_max
    call param_get_keyword('dos_adpt_smr_max', found, &
                           r_value=dos_data%adpt_smr_max)
    if (dos_data%adpt_smr_max <= 0._dp) call io_error &
      ('Error: dos_adpt_smr_max must be greater than zero')

    dos_data%smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('dos_smr_fixed_en_width', found, &
                           r_value=dos_data%smr_fixed_en_width)
    if (found .and. (dos_data%smr_fixed_en_width < 0._dp)) &
      call io_error('Error: dos_smr_fixed_en_width must be greater than or equal to zero')

!    dos_gaussian_width        = 0.1_dp
!    call param_get_keyword('dos_gaussian_width',found,r_value=dos_gaussian_width)

!    dos_plot_format           = 'gnuplot'
!    call param_get_keyword('dos_plot_format',found,c_value=dos_plot_format)

    call param_get_range_vector('dos_project', found, &
                                dos_data%num_project, lcount=.true.)
    if (found) then
      if (dos_data%num_project < 1) call io_error('Error: problem reading dos_project')
      if (allocated(dos_data%project)) deallocate (dos_data%project)
      allocate (dos_data%project(dos_data%num_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dos_project in param_read')
      call param_get_range_vector('dos_project', found, &
                                  dos_data%num_project, .false., &
                                  dos_data%project)
      if (any(dos_data%project < 1) .or. &
          any(dos_data%project > num_wann)) &
        call io_error('Error: dos_project asks for out-of-range Wannier functions')
    else
      ! by default plot all
      dos_data%num_project = num_wann
      if (allocated(dos_data%project)) deallocate (dos_data%project)
      allocate (dos_data%project(dos_data%num_project), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dos_project in param_read')
      do i = 1, dos_data%num_project
        dos_data%project(i) = i
      end do
    endif
  end subroutine param_pw90_read_26

  subroutine param_pw90_read_31
    ! These must be read here, before the check on the existence of the .eig file!
    implicit none
    logical :: found

    pw90_calcs%geninterp = .false.
    call param_get_keyword('geninterp', found, l_value=pw90_calcs%geninterp)
    pw90_calcs%boltzwann = .false.
    call param_get_keyword('boltzwann', found, l_value=pw90_calcs%boltzwann)
  end subroutine param_pw90_read_31

  subroutine param_pw90_read_34(smr_index)
    ! [gp-begin, Jun 1, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! General band interpolator (geninterp)
    !%%%%%%%%%%%%%%%%%%%%
    use w90_io, only: io_error
    implicit none
    integer, intent(in) :: smr_index
    logical :: found, found2

    geninterp%alsofirstder = .false.
    call param_get_keyword('geninterp_alsofirstder', found, l_value=geninterp%alsofirstder)
    geninterp%single_file = .true.
    call param_get_keyword('geninterp_single_file', found, l_value=geninterp%single_file)
    ! [gp-end, Jun 1, 2012]

    ! [gp-begin, Apr 12, 2012]
    !%%%%%%%%%%%%%%%%%%%%
    ! Boltzmann transport
    !%%%%%%%%%%%%%%%%%%%%
    ! Note: to be put AFTER the disentanglement routines!

    boltz%calc_also_dos = .false.
    call param_get_keyword('boltz_calc_also_dos', found, l_value=boltz%calc_also_dos)

    boltz%calc_also_dos = boltz%calc_also_dos .and. pw90_calcs%boltzwann

    ! 0 means the normal 3d case for the calculation of the Seebeck coefficient
    ! The other valid possibilities are 1,2,3 for x,y,z respectively
    boltz%dir_num_2d = 0
    call param_get_keyword('boltz_2d_dir', found, c_value=boltz_2d_dir)
    if (found) then
      if (trim(boltz_2d_dir) == 'no') then
        boltz%dir_num_2d = 0
      elseif (trim(boltz_2d_dir) == 'x') then
        boltz%dir_num_2d = 1
      elseif (trim(boltz_2d_dir) == 'y') then
        boltz%dir_num_2d = 2
      elseif (trim(boltz_2d_dir) == 'z') then
        boltz%dir_num_2d = 3
      else
        call io_error('Error: boltz_2d_dir can only be "no", "x", "y" or "z".')
      end if
    end if

    boltz%dos_energy_step = 0.001_dp
    call param_get_keyword('boltz_dos_energy_step', found, r_value=boltz%dos_energy_step)
    if (found .and. (boltz%dos_energy_step <= 0._dp)) &
      call io_error('Error: boltz_dos_energy_step must be positive')

    if (allocated(eigval)) then
      boltz%dos_energy_min = minval(eigval) - 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      boltz%dos_energy_min = -1.0_dp
    end if
    call param_get_keyword('boltz_dos_energy_min', found, r_value=boltz%dos_energy_min)
    if (allocated(eigval)) then
      boltz%dos_energy_max = maxval(eigval) + 0.6667_dp
    else
      ! Boltz_dos cannot run if eigval is not allocated.
      ! We just set here a default numerical value.
      boltz%dos_energy_max = 0.0_dp
    end if
    call param_get_keyword('boltz_dos_energy_max', found, r_value=boltz%dos_energy_max)
    if (boltz%dos_energy_max <= boltz%dos_energy_min) &
      call io_error('Error: boltz_dos_energy_max must be greater than boltz_dos_energy_min')

    boltz%dos_adpt_smr = adpt_smr
    call param_get_keyword('boltz_dos_adpt_smr', found, l_value=boltz%dos_adpt_smr)

    boltz%dos_adpt_smr_fac = adpt_smr_fac
    call param_get_keyword('boltz_dos_adpt_smr_fac', found, r_value=boltz%dos_adpt_smr_fac)
    if (found .and. (boltz%dos_adpt_smr_fac <= 0._dp)) &
      call io_error('Error: boltz_dos_adpt_smr_fac must be greater than zero')

    boltz%dos_adpt_smr_max = adpt_smr_max
    call param_get_keyword('boltz_dos_adpt_smr_max', found, r_value=boltz%dos_adpt_smr_max)
    if (boltz%dos_adpt_smr_max <= 0._dp) call io_error &
      ('Error: boltz_dos_adpt_smr_max must be greater than zero')

    boltz%dos_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('boltz_dos_smr_fixed_en_width', found, r_value=boltz%dos_smr_fixed_en_width)
    if (found .and. (boltz%dos_smr_fixed_en_width < 0._dp)) &
      call io_error('Error: boltz_dos_smr_fixed_en_width must be greater than or equal to zero')

    boltz%mu_min = -999._dp
    call param_get_keyword('boltz_mu_min', found, r_value=boltz%mu_min)
    if ((.not. found) .and. pw90_calcs%boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_min provided')
    boltz%mu_max = -999._dp
    call param_get_keyword('boltz_mu_max', found2, r_value=boltz%mu_max)
    if ((.not. found2) .and. pw90_calcs%boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_max provided')
    if (found .and. found2 .and. (boltz%mu_max < boltz%mu_min)) &
      call io_error('Error: boltz_mu_max must be greater than boltz_mu_min')
    boltz%mu_step = 0._dp
    call param_get_keyword('boltz_mu_step', found, r_value=boltz%mu_step)
    if ((.not. found) .and. pw90_calcs%boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_mu_step provided')
    if (found .and. (boltz%mu_step <= 0._dp)) &
      call io_error('Error: boltz_mu_step must be greater than zero')

    boltz%temp_min = -999._dp
    call param_get_keyword('boltz_temp_min', found, r_value=boltz%temp_min)
    if ((.not. found) .and. pw90_calcs%boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_min provided')
    boltz%temp_max = -999._dp
    call param_get_keyword('boltz_temp_max', found2, r_value=boltz%temp_max)
    if ((.not. found2) .and. pw90_calcs%boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_max provided')
    if (found .and. found2 .and. (boltz%temp_max < boltz%temp_min)) &
      call io_error('Error: boltz_temp_max must be greater than boltz_temp_min')
    if (found .and. (boltz%temp_min <= 0._dp)) &
      call io_error('Error: boltz_temp_min must be greater than zero')
    boltz%temp_step = 0._dp
    call param_get_keyword('boltz_temp_step', found, r_value=boltz%temp_step)
    if ((.not. found) .and. pw90_calcs%boltzwann) &
      call io_error('Error: BoltzWann required but no boltz_temp_step provided')
    if (found .and. (boltz%temp_step <= 0._dp)) &
      call io_error('Error: boltz_temp_step must be greater than zero')

    ! The interpolation mesh is read later on

    ! By default, the energy step for the TDF is 1 meV
    boltz%tdf_energy_step = 0.001_dp
    call param_get_keyword('boltz_tdf_energy_step', found, r_value=boltz%tdf_energy_step)
    if (boltz%tdf_energy_step <= 0._dp) &
      call io_error('Error: boltz_tdf_energy_step must be greater than zero')

    ! For TDF: TDF smeared in a NON-adaptive way; value in eV, default = 0._dp
    ! (i.e., no smearing)
    boltz%TDF_smr_fixed_en_width = smr_fixed_en_width
    call param_get_keyword('boltz_tdf_smr_fixed_en_width', found, r_value=boltz%TDF_smr_fixed_en_width)
    if (found .and. (boltz%TDF_smr_fixed_en_width < 0._dp)) &
      call io_error('Error: boltz_TDF_smr_fixed_en_width must be greater than or equal to zero')

    ! By default: use the "global" smearing index
    boltz%TDF_smr_index = smr_index
    call param_get_keyword('boltz_tdf_smr_type', found, c_value=ctmp)
    if (found) boltz%TDF_smr_index = get_smearing_index(ctmp, 'boltz_tdf_smr_type')

    ! By default: use the "global" smearing index
    boltz%dos_smr_index = smr_index
    call param_get_keyword('boltz_dos_smr_type', found, c_value=ctmp)
    if (found) boltz%dos_smr_index = get_smearing_index(ctmp, 'boltz_dos_smr_type')

    ! By default: use the "global" smearing index
    dos_data%smr_index = smr_index
    call param_get_keyword('dos_smr_type', found, c_value=ctmp)
    if (found) dos_data%smr_index = get_smearing_index(ctmp, 'dos_smr_type')

    ! By default: use the "global" smearing index
    berry%kubo_smr_index = smr_index
    call param_get_keyword('kubo_smr_type', found, c_value=ctmp)
    if (found) berry%kubo_smr_index = get_smearing_index(ctmp, 'kubo_smr_type')

    ! By default: use the "global" smearing index
    gyrotropic%smr_index = smr_index
    call param_get_keyword('gyrotropic_smr_type', found, c_value=ctmp)
    if (found) gyrotropic%smr_index = get_smearing_index(ctmp, 'gyrotropic_smr_type')

    ! By default: 10 fs relaxation time
    boltz%relax_time = 10._dp
    call param_get_keyword('boltz_relax_time', found, r_value=boltz%relax_time)

    boltz%bandshift = .false.
    call param_get_keyword('boltz_bandshift', found, l_value=boltz%bandshift)
    boltz%bandshift = boltz%bandshift .and. pw90_calcs%boltzwann

    boltz%bandshift_firstband = 0
    call param_get_keyword('boltz_bandshift_firstband', found, i_value=boltz%bandshift_firstband)
    if (boltz%bandshift .and. (.not. found)) &
      call io_error('Error: boltz_bandshift required but no boltz_bandshift_firstband provided')
    boltz%bandshift_energyshift = 0._dp
    call param_get_keyword('boltz_bandshift_energyshift', found, r_value=boltz%bandshift_energyshift)
    if (boltz%bandshift .and. (.not. found)) &
      call io_error('Error: boltz_bandshift required but no boltz_bandshift_energyshift provided')
    ! [gp-end, Apr 12, 2012]
  end subroutine param_pw90_read_34

  subroutine param_pw90_read_36
    use w90_constants, only: cmplx_i
    use w90_io, only: io_error
    implicit none
    integer :: i, ierr
    logical :: found

    if (dis_data%frozen_states) then
      dos_data%energy_max = dis_data%froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      dos_data%energy_max = maxval(eigval) + 0.6667_dp
    else
      dos_data%energy_max = dis_data%win_max + 0.6667_dp
    end if
    call param_get_keyword('dos_energy_max', found, &
                           r_value=dos_data%energy_max)

    if (allocated(eigval)) then
      dos_data%energy_min = minval(eigval) - 0.6667_dp
    else
      dos_data%energy_min = dis_data%win_min - 0.6667_dp
    end if
    call param_get_keyword('dos_energy_min', found, &
                           r_value=dos_data%energy_min)

    kubo_freq_min = 0.0_dp
    gyrotropic_freq_min = kubo_freq_min
    call param_get_keyword('kubo_freq_min', found, r_value=kubo_freq_min)
    !
    if (dis_data%frozen_states) then
      kubo_freq_max = dis_data%froz_max - fermi%energy_list(1) + 0.6667_dp
    elseif (allocated(eigval)) then
      kubo_freq_max = maxval(eigval) - minval(eigval) + 0.6667_dp
    else
      kubo_freq_max = dis_data%win_max - dis_data%win_min + 0.6667_dp
    end if
    gyrotropic_freq_max = kubo_freq_max
    call param_get_keyword('kubo_freq_max', found, r_value=kubo_freq_max)

    !
    kubo_freq_step = 0.01_dp
    call param_get_keyword('kubo_freq_step', found, r_value=kubo_freq_step)
    if (found .and. kubo_freq_step < 0.0_dp) call io_error( &
      'Error: kubo_freq_step must be positive')
    !
    berry%kubo_nfreq = nint((kubo_freq_max - kubo_freq_min)/kubo_freq_step) + 1
    if (berry%kubo_nfreq <= 1) berry%kubo_nfreq = 2
    kubo_freq_step = (kubo_freq_max - kubo_freq_min)/(berry%kubo_nfreq - 1)
    !
    if (allocated(berry%kubo_freq_list)) deallocate (berry%kubo_freq_list)
    allocate (berry%kubo_freq_list(berry%kubo_nfreq), stat=ierr)
    if (ierr /= 0) &
      call io_error('Error allocating kubo_freq_list in param_read')
    do i = 1, berry%kubo_nfreq
      berry%kubo_freq_list(i) = kubo_freq_min &
                                + (i - 1)*(kubo_freq_max - kubo_freq_min)/(berry%kubo_nfreq - 1)
    enddo
    !
    ! TODO: Alternatively, read list of (complex) frequencies; kubo_nfreq is
    !       the length of the list

    gyrotropic_freq_step = 0.01_dp
    call param_get_keyword('gyrotropic_freq_min', found, r_value=gyrotropic_freq_min)
    call param_get_keyword('gyrotropic_freq_max', found, r_value=gyrotropic_freq_max)
    call param_get_keyword('gyrotropic_freq_step', found, r_value=gyrotropic_freq_step)
    gyrotropic%nfreq = nint((gyrotropic_freq_max - gyrotropic_freq_min)/gyrotropic_freq_step) + 1
    if (gyrotropic%nfreq <= 1) gyrotropic%nfreq = 2
    gyrotropic_freq_step = (gyrotropic_freq_max - gyrotropic_freq_min)/(gyrotropic%nfreq - 1)
    if (allocated(gyrotropic%freq_list)) deallocate (gyrotropic%freq_list)
    allocate (gyrotropic%freq_list(gyrotropic%nfreq), stat=ierr)
    if (ierr /= 0) &
      call io_error('Error allocating gyrotropic_freq_list in param_read')
    do i = 1, gyrotropic%nfreq
      gyrotropic%freq_list(i) = gyrotropic_freq_min &
                                + (i - 1)*(gyrotropic_freq_max - gyrotropic_freq_min)/(gyrotropic%nfreq - 1) &
                                + cmplx_i*gyrotropic%smr_fixed_en_width
    enddo

    if (dis_data%frozen_states) then
      berry%kubo_eigval_max = dis_data%froz_max + 0.6667_dp
    elseif (allocated(eigval)) then
      berry%kubo_eigval_max = maxval(eigval) + 0.6667_dp
    else
      berry%kubo_eigval_max = dis_data%win_max + 0.6667_dp
    end if
    gyrotropic%eigval_max = berry%kubo_eigval_max

    call param_get_keyword('kubo_eigval_max', found, r_value=berry%kubo_eigval_max)
    call param_get_keyword('gyrotropic_eigval_max', found, r_value=gyrotropic%eigval_max)
  end subroutine param_pw90_read_36

  subroutine param_pw90_read_38
    implicit none
    logical :: found

    berry%sc_eta = 0.04
    call param_get_keyword('sc_eta', found, r_value=berry%sc_eta)

    berry%sc_w_thr = 5.0d0
    call param_get_keyword('sc_w_thr', found, r_value=berry%sc_w_thr)
  end subroutine param_pw90_read_38

  subroutine param_pw90_read_40
    implicit none

    ! To be called after having read the global flag
    call get_module_kmesh(moduleprefix='boltz', &
                          should_be_defined=pw90_calcs%boltzwann, &
                          module_kmesh=boltz%kmesh, &
                          module_kmesh_spacing=boltz%kmesh_spacing)

    call get_module_kmesh(moduleprefix='berry', &
                          should_be_defined=pw90_calcs%berry, &
                          module_kmesh=berry%kmesh, &
                          module_kmesh_spacing=berry%kmesh_spacing)

    call get_module_kmesh(moduleprefix='gyrotropic', &
                          should_be_defined=pw90_calcs%gyrotropic, &
                          module_kmesh=gyrotropic%kmesh, &
                          module_kmesh_spacing=gyrotropic%kmesh_spacing)

    call get_module_kmesh(moduleprefix='spin', &
                          should_be_defined=pw90_common%spin_moment, &
                          module_kmesh=pw90_spin%spin_kmesh, &
                          module_kmesh_spacing=pw90_spin%spin_kmesh_spacing)

    call get_module_kmesh(moduleprefix='dos', &
                          should_be_defined=pw90_calcs%dos, &
                          module_kmesh=dos_data%kmesh, &
                          module_kmesh_spacing=dos_data%kmesh_spacing)
  end subroutine param_pw90_read_40

  subroutine get_module_kmesh(moduleprefix, should_be_defined, module_kmesh, module_kmesh_spacing)
    !! This function reads and sets the interpolation mesh variables needed by a given module
    !>
    !!  This function MUST be called after having read the global kmesh and kmesh_spacing!!
    !!  if the user didn't provide an interpolation_mesh_spacing, it is set to -1, so that
    !!       one can check in the code what the user asked for
    !!  The function takes care also of setting the default value to the global one if no local
    !!       keyword is defined
    use w90_io, only: io_error
    character(len=*), intent(in)       :: moduleprefix
    !!The prefix that is appended before the name of the variables. In particular,
    !!if the prefix is for instance XXX, the two variables that are read from the
    !!input file are XXX_kmesh and XXX_kmesh_spacing.
    logical, intent(in)                :: should_be_defined
    !! A logical flag: if it is true, at the end the code stops if no value is specified.
    !! Define it to .false. if no check should be performed.
    !! Often, you can simply pass the flag that activates the module itself.
    integer, dimension(3), intent(out) :: module_kmesh
    !! the integer array (length 3) where the interpolation mesh will be saved
    real(kind=dp), intent(out)         :: module_kmesh_spacing
    !! the real number on which the min mesh spacing is saved. -1 if it the
    !!user specifies in input the mesh and not the mesh_spacing

    logical :: found, found2
    integer :: i

    ! Default values
    module_kmesh_spacing = -1._dp
    module_kmesh = 0
    call param_get_keyword(trim(moduleprefix)//'_kmesh_spacing', found, r_value=module_kmesh_spacing)
    if (found) then
      if (module_kmesh_spacing .le. 0._dp) &
        call io_error('Error: '//trim(moduleprefix)//'_kmesh_spacing must be greater than zero')

      call internal_set_kmesh(module_kmesh_spacing, recip_lattice, module_kmesh)
    end if
    call param_get_vector_length(trim(moduleprefix)//'_kmesh', found2, length=i)
    if (found2) then
      if (found) &
        call io_error('Error: cannot set both '//trim(moduleprefix)//'_kmesh and ' &
                      //trim(moduleprefix)//'_kmesh_spacing')
      if (i .eq. 1) then
        call param_get_keyword_vector(trim(moduleprefix)//'_kmesh', found2, 1, i_value=module_kmesh)
        module_kmesh(2) = module_kmesh(1)
        module_kmesh(3) = module_kmesh(1)
      elseif (i .eq. 3) then
        call param_get_keyword_vector(trim(moduleprefix)//'_kmesh', found2, 3, i_value=module_kmesh)
      else
        call io_error('Error: '//trim(moduleprefix)// &
                      '_kmesh must be provided as either one integer or a vector of 3 integers')
      end if
      if (any(module_kmesh <= 0)) &
        call io_error('Error: '//trim(moduleprefix)//'_kmesh elements must be greater than zero')
    end if

    if ((found .eqv. .false.) .and. (found2 .eqv. .false.)) then
      ! This is the case where no  "local" interpolation k-mesh is provided in the input
      if (global_kmesh_set) then
        module_kmesh = kmesh
        ! I set also boltz_kmesh_spacing so that I can check if it is < 0 or not, and if it is
        ! > 0 I can print on output the mesh spacing that was chosen
        module_kmesh_spacing = kmesh_spacing
      else
        if (should_be_defined) &
          call io_error('Error: '//trim(moduleprefix)//' module required, but no interpolation mesh given.')
      end if
    end if
  end subroutine get_module_kmesh

!===================================================================
  subroutine param_postw90_write
    !==================================================================!
    !                                                                  !
    !! write postw90 parameters to stdout
    !                                                                  !
    !===================================================================

    implicit none

    integer :: i, loop, nat, nsp
    real(kind=dp) :: cell_volume

    ! System
    write (stdout, *)
    write (stdout, '(36x,a6)') '------'
    write (stdout, '(36x,a6)') 'SYSTEM'
    write (stdout, '(36x,a6)') '------'
    write (stdout, *)
    if (param_input%lenconfac .eq. 1.0_dp) then
      write (stdout, '(30x,a21)') 'Lattice Vectors (Ang)'
    else
      write (stdout, '(28x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write (stdout, 101) 'a_1', (real_lattice(1, I)*param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'a_2', (real_lattice(2, I)*param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'a_3', (real_lattice(3, I)*param_input%lenconfac, i=1, 3)
    write (stdout, *)
    cell_volume = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - real_lattice(3, 2)*real_lattice(2, 3)) + &
                  real_lattice(1, 2)*(real_lattice(2, 3)*real_lattice(3, 1) - real_lattice(3, 3)*real_lattice(2, 1)) + &
                  real_lattice(1, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - real_lattice(3, 1)*real_lattice(2, 2))
    write (stdout, '(19x,a17,3x,f11.5)', advance='no') &
      'Unit Cell Volume:', cell_volume*param_input%lenconfac**3
    if (param_input%lenconfac .eq. 1.0_dp) then
      write (stdout, '(2x,a7)') '(Ang^3)'
    else
      write (stdout, '(2x,a8)') '(Bohr^3)'
    endif
    write (stdout, *)
    if (param_input%lenconfac .eq. 1.0_dp) then
      write (stdout, '(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
      write (stdout, '(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write (stdout, 101) 'b_1', (recip_lattice(1, I)/param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'b_2', (recip_lattice(2, I)/param_input%lenconfac, i=1, 3)
    write (stdout, 101) 'b_3', (recip_lattice(3, I)/param_input%lenconfac, i=1, 3)
    write (stdout, *) ' '
    ! Atoms
    if (atoms%num_atoms > 0) then
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
      if (param_input%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
      else
        write (stdout, '(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
      endif
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      do nsp = 1, atoms%num_species
        do nat = 1, atoms%species_num(nsp)
          write (stdout, '(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') &
  &                 '|', atoms%symbol(nsp), nat, atoms%pos_frac(:, nat, nsp),&
  &                 '|', atoms%pos_cart(:, nat, nsp)*param_input%lenconfac, '|'
        end do
      end do
      write (stdout, '(1x,a)') '*----------------------------------------------------------------------------*'
    else
      write (stdout, '(25x,a)') 'No atom positions specified'
    end if
    write (stdout, *) ' '
    ! Main
    write (stdout, *) ' '

    write (stdout, '(1x,a78)') '*-------------------------------- POSTW90 -----------------------------------*'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of Wannier Functions               :', num_wann, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Number of electrons per state             :', &
      param_input%num_elec_per_state, '|'
    if (abs(pw90_common%scissors_shift) > 1.0e-7_dp .or. param_input%iprint > 0) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Scissor shift applied to conduction bands :', pw90_common%scissors_shift, '|'
      if (param_input%num_valence_bands > 0) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Number of valence bands                   :', &
          param_input%num_valence_bands, '|'
      else
        write (stdout, '(1x,a78)') '|  Number of valence bands                   :       not defined             |'
      endif
    endif
    if (pw90_common%spin_decomp .or. param_input%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Spin decomposition                        :', pw90_common%spin_decomp, '|'
    if (pw90_common%spin_moment .or. param_input%iprint > 2) &
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Spin moment                       :', pw90_common%spin_moment, '|'
    if (pw90_common%spin_decomp .or. pw90_common%spin_moment .or. param_input%iprint > 2) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Polar angle of spin quantisation axis     :', pw90_spin%spin_axis_polar, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Azimuthal angle of spin quantisation axis :', pw90_spin%spin_axis_azimuth, '|'
      if (postw90_oper%spn_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  Spn file-type                   :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  Spn file-type                   :', 'unformatted', '|'
      endif
      if (postw90_oper%uHu_formatted) then
        write (stdout, '(1x,a46,9x,a9,13x,a1)') '|  uHu file-type                   :', 'formatted', '|'
      else
        write (stdout, '(1x,a46,7x,a11,13x,a1)') '|  uHu file-type                   :', 'unformatted', '|'
      endif
    end if

    if (size(fermi%energy_list) == 1) then
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi energy (eV)                         :', fermi%energy_list(1), '|'
    else
      write (stdout, '(1x,a21,I8,a12,f8.3,a4,f8.3,a3,13x,a1)') '|  Fermi energy     :', size(fermi%energy_list), &
        ' steps from ', fermi%energy_list(1), ' to ', &
        fermi%energy_list(size(fermi%energy_list)), ' eV', '|'
    end if

    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :', param_input%iprint, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :', param_input%timing_level, '|'
    write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Optimisation (0=memory, 3=speed)          :', param_input%optimisation, '|'
    write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :', trim(param_input%length_unit), '|'
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    write (stdout, '(1x,a78)') '*------------------------ Global Smearing Parameters ------------------------*'
    if (adpt_smr) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', adpt_smr_fac, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width (eV)       :', adpt_smr_max, '|'

    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', smr_fixed_en_width, '|'
    endif
    write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', trim(param_get_smearing_type(smr_index)), '|'
    if (global_kmesh_set) then
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       T', '|'
      if (kmesh_spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          kmesh(1), 'x', kmesh(2), 'x', kmesh(3), ' Spacing = ', kmesh_spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , kmesh(1), 'x', kmesh(2), 'x', kmesh(3), '|'
      endif
    else
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Global interpolation k-points defined     :', '       F', '|'
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    ! DOS
    if (pw90_calcs%dos .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*---------------------------------- DOS -------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plotting Density of States                :', pw90_calcs%dos, '|'
      if (dos_data%num_project > 1) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Wannier Projected DOS             :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', dos_data%energy_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', dos_data%energy_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', dos_data%energy_step, '|'
      if (dos_data%adpt_smr .eqv. adpt_smr .and. &
          dos_data%adpt_smr_fac == adpt_smr_fac .and. &
          dos_data%adpt_smr_max == adpt_smr_max .and. &
          dos_data%smr_fixed_en_width == smr_fixed_en_width .and. &
          smr_index == dos_data%smr_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (dos_data%adpt_smr) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', &
            dos_data%adpt_smr_fac, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', &
            dos_data%adpt_smr_max, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
            dos_data%smr_fixed_en_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', &
          trim(param_get_smearing_type(dos_data%smr_index)), '|'
      endif
      if (kmesh(1) == dos_data%kmesh(1) .and. &
          kmesh(2) == dos_data%kmesh(2) .and. &
          kmesh(3) == dos_data%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (dos_data%kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            dos_data%kmesh(1), 'x', dos_data%kmesh(2), 'x', &
            dos_data%kmesh(3), ' Spacing = ', dos_data%kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :', &
            dos_data%kmesh(1), 'x', dos_data%kmesh(2), 'x', &
            dos_data%kmesh(3), '|'
        endif
      endif
    endif
    write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'

    if (pw90_calcs%kpath .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KPATH ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a path in k-space   :', pw90_calcs%kpath, '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', kpath%num_points, '|'
      if (index(kpath%task, 'bands') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy bands                         :', '       F', '|'
      endif
      if (index(kpath%task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature                      :', '       F', '|'
      endif
      if (index(kpath%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(kpath%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the bands    :', trim(kpath%bands_colour), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
      write (stdout, '(1x,a78)') '|   K-space path sections:                                                   |'
      if (spec_points%bands_num_spec_points == 0) then
        write (stdout, '(1x,a78)') '|     None defined                                                           |'
      else
        do loop = 1, spec_points%bands_num_spec_points, 2
          write (stdout, '(1x,a10,2x,a1,2x,3F7.3,5x,a3,2x,a1,2x,3F7.3,7x,a1)') '|    From:', &
            spec_points%bands_label(loop), (spec_points%bands_spec_points(i, loop), i=1, 3), &
            'To:', spec_points%bands_label(loop + 1), (spec_points%bands_spec_points(i, loop + 1), i=1, 3), '|'
        end do
      end if
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%kslice .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- KSLICE -----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Plot Properties along a slice in k-space  :', pw90_calcs%kslice, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Fermi level used for slice                :', fermi%energy_list(1), '|'
      write (stdout, '(1x,a46,10x,I8,13x,a1)') '|  Divisions along first kpath section       :', kpath%num_points, '|'
      if (index(kslice%task, 'fermi_lines') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot energy contours (fermi lines)        :', '       F', '|'
      endif
      if (index(kslice%task, 'curv') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot Berry curvature (sum over occ states):', '       F', '|'
      endif
      if (index(kslice%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot orbital magnetisation contribution   :', '       F', '|'
      endif
      if (index(kslice%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Plot spin Hall conductivity contribution  :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Property used to colour code the lines    :', &
        trim(kslice%fermi_lines_colour), '|'
      write (stdout, '(1x,a78)') '|  2D slice parameters (in reduced coordinates):                             |'
      write (stdout, '(1x,a14,2x,3F8.3,37x,a1)') '|     Corner: ', (kslice%corner(i), i=1, 3), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector1: ', (kslice%b1(i), i=1, 3), ' Divisions:', kslice%kmesh2d(1), '|'
      write (stdout, '(1x,a14,2x,3F8.3,10x,a12,2x,i4,9x,a1)') &
        '|    Vector2: ', (kslice%b2(i), i=1, 3), ' Divisions:', kslice%kmesh2d(1), '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%berry .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- BERRY ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Berry Phase related properties    :', pw90_calcs%berry, '|'
      if (index(berry%task, 'kubo') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Optical Conductivity and JDOS     :', '       F', '|'
      endif
      if (index(berry%task, 'ahc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Anomalous Hall Conductivity       :', '       F', '|'
      endif
      if (index(berry%task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Shift Current                     :', '       F', '|'
      endif
      if (index(berry%task, 'morb') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Orbital Magnetisation             :', '       F', '|'
      endif
      if (index(berry%task, 'shc') > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       T', '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Compute Spin Hall Conductivity            :', '       F', '|'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for optical responses     :', kubo_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency for optical responses     :', kubo_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for optical responses           :', kubo_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue for optical responses    :', berry%kubo_eigval_max, '|'
      if (index(berry%task, 'sc') > 0) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing factor for shift current         :', berry%sc_eta, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Frequency theshold for shift current      :', berry%sc_w_thr, '|'
        write (stdout, '(1x,a46,1x,a27,3x,a1)') '|  Bloch sums                                :', &
          trim(param_get_convention_type(berry%sc_phase_conv)), '|'
      end if
      if (berry%kubo_adpt_smr .eqv. adpt_smr .and. &
          berry%kubo_adpt_smr_fac == adpt_smr_fac .and. &
          berry%kubo_adpt_smr_max == adpt_smr_max &
          .and. berry%kubo_smr_fixed_en_width == smr_fixed_en_width .and. &
          smr_index == berry%kubo_smr_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        if (berry%kubo_adpt_smr) then
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive width smearing                   :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Adaptive smearing factor                  :', berry%kubo_adpt_smr_fac, '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum allowed smearing width            :', berry%kubo_adpt_smr_max, '|'
        else
          write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
          write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
            berry%kubo_smr_fixed_en_width, '|'
        endif
        write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', trim(param_get_smearing_type(berry%kubo_smr_index)), '|'
      endif
      if (kmesh(1) == berry%kmesh(1) .and. kmesh(2) == berry%kmesh(2) .and. kmesh(3) == berry%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (berry%kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            berry%kmesh(1), 'x', berry%kmesh(2), 'x', berry%kmesh(3), ' Spacing = ', berry%kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , berry%kmesh(1), 'x', berry%kmesh(2), 'x', berry%kmesh(3), '|'
        endif
      endif
      if (berry%curv_adpt_kmesh > 1) then
        write (stdout, '(1x,a46,10x,i8,13x,a1)') '|  Using an adaptive refinement mesh of size :', berry%curv_adpt_kmesh, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Threshold for adaptive refinement         :', &
          berry%curv_adpt_kmesh_thresh, '|'
      else
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    none', '|'
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%gyrotropic .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*--------------------------------- GYROTROPIC   ------------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '| Compute Gyrotropic properties              :', pw90_calcs%gyrotropic, '|'
      write (stdout, '(1x,a46,10x,a20,1x,a1)') '| gyrotropic_task                            :', gyrotropic%task, '|'
      call parameters_gyro_write_task(gyrotropic%task, '-d0', 'calculate the D tensor')
      call parameters_gyro_write_task(gyrotropic%task, '-dw', 'calculate the tildeD tensor')
      call parameters_gyro_write_task(gyrotropic%task, '-c', 'calculate the C tensor')
      call parameters_gyro_write_task(gyrotropic%task, '-k', 'calculate the K tensor')
      call parameters_gyro_write_task(gyrotropic%task, '-noa', 'calculate the interbad natural optical activity')
      call parameters_gyro_write_task(gyrotropic%task, '-dos', 'calculate the density of states')

      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Lower frequency for tildeD,NOA            :', gyrotropic_freq_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper frequency                           :', gyrotropic_freq_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for frequency                   :', gyrotropic_freq_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Upper eigenvalue                          :', gyrotropic%eigval_max, '|'
      if (gyrotropic%smr_fixed_en_width == smr_fixed_en_width .and. smr_index == gyrotropic%smr_index) then
        write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
      else
        write (stdout, '(1x,a78)') '|  Using local  smearing parameters                                          |'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Fixed width smearing                      :', '       T', '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Smearing width                            :', &
        gyrotropic%smr_fixed_en_width, '|'
      write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function                         :', &
        trim(param_get_smearing_type(gyrotropic%smr_index)), '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  degen_thresh                              :', gyrotropic%degen_thresh, '|'

      if (kmesh(1) == gyrotropic%kmesh(1) .and. kmesh(2) == gyrotropic%kmesh(2) .and. kmesh(3) == gyrotropic%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      elseif (gyrotropic%kmesh_spacing > 0.0_dp) then
        write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
          gyrotropic%kmesh(1), 'x', gyrotropic%kmesh(2), 'x', gyrotropic%kmesh(3), ' Spacing = ', gyrotropic%kmesh_spacing, '|'
      else
        write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
          , gyrotropic%kmesh(1), 'x', gyrotropic%kmesh(2), 'x', gyrotropic%kmesh(3), '|'
      endif
      write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  Adaptive refinement                       :', '    not implemented', '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%boltzwann .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------------- BOLTZWANN ----------------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Boltzmann transport properties    :', pw90_calcs%boltzwann, '|'
      if (boltz%dir_num_2d > 0) then
        write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  2d structure: non-periodic dimension  :', trim(boltz_2d_dir), '|'
      else
        write (stdout, '(1x,a78)') '|  3d Structure                              :                 T             |'
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Relaxation Time (fs)                      :', boltz%relax_time, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Chemical Potential (eV)  :', boltz%mu_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Chemical Potential (eV)  :', boltz%mu_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Chemical Potential (eV)     :', boltz%mu_step, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum Value of Temperature (K)          :', boltz%temp_min, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum Value of Temperature (K)          :', boltz%temp_max, '|'
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for Temperature (K)             :', boltz%temp_step, '|'

      if (kmesh(1) == boltz%kmesh(1) .and. kmesh(2) == boltz%kmesh(2) .and. kmesh(3) == boltz%kmesh(3)) then
        write (stdout, '(1x,a78)') '|  Using global k-point set for interpolation                                |'
      else
        if (boltz%kmesh_spacing > 0.0_dp) then
          write (stdout, '(1x,a15,i4,1x,a1,i4,1x,a1,i4,16x,a11,f8.3,11x,1a)') '|  Grid size = ', &
            boltz%kmesh(1), 'x', boltz%kmesh(2), 'x', boltz%kmesh(3), ' Spacing = ', boltz%kmesh_spacing, '|'
        else
          write (stdout, '(1x,a46,2x,i4,1x,a1,i4,1x,a1,i4,13x,1a)') '|  Grid size                                 :' &
            , boltz%kmesh(1), 'x', boltz%kmesh(2), 'x', boltz%kmesh(3), '|'
        endif
      endif
      write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Step size for TDF (eV)                    :', boltz%tdf_energy_step, '|'
      write (stdout, '(1x,a25,5x,a43,4x,a1)') '|  TDF Smearing Function ', trim(param_get_smearing_type(boltz%tdf_smr_index)), '|'
      if (boltz%tdf_smr_fixed_en_width > 0.0_dp) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
          '|  TDF fixed Smearing width (eV)             :', boltz%tdf_smr_fixed_en_width, '|'
      else
        write (stdout, '(1x,a78)') '|  TDF fixed Smearing width                  :         unsmeared             |'
      endif
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute DOS at same time                  :', boltz%calc_also_dos, '|'
      if (boltz%calc_also_dos .and. param_input%iprint > 2) then
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Minimum energy range for DOS plot         :', boltz%dos_energy_min, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Maximum energy range for DOS plot         :', boltz%dos_energy_max, '|'
        write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  Energy step for DOS plot                  :', boltz%dos_energy_step, '|'
        if (boltz%dos_adpt_smr .eqv. adpt_smr .and. boltz%dos_adpt_smr_fac == adpt_smr_fac &
            .and. boltz%dos_adpt_smr_max == adpt_smr_max &
            .and. boltz%dos_smr_fixed_en_width == smr_fixed_en_width .and. smr_index == boltz%dos_smr_index) then
          write (stdout, '(1x,a78)') '|  Using global smearing parameters                                          |'
        else
          if (boltz%dos_adpt_smr) then
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Adaptive width smearing               :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Adaptive smearing factor              :', boltz%dos_adpt_smr_fac, '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') &
              '|  DOS Maximum allowed smearing width        :', boltz%dos_adpt_smr_max, '|'
          else
            write (stdout, '(1x,a46,10x,a8,13x,a1)') '|  DOS Fixed width smearing                  :', '       T', '|'
            write (stdout, '(1x,a46,10x,f8.3,13x,a1)') '|  DOS Smearing width                         :', &
              boltz%dos_smr_fixed_en_width, '|'
          endif
          write (stdout, '(1x,a21,5x,a47,4x,a1)') '|  Smearing Function ', trim(param_get_smearing_type(boltz%dos_smr_index)), '|'
        endif
      endif
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

    if (pw90_calcs%geninterp .or. param_input%iprint > 2) then
      write (stdout, '(1x,a78)') '*------------------------Generic Band Interpolation--------------------------*'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Compute Properties at given k-points      :', pw90_calcs%geninterp, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Calculate band gradients                  :', geninterp%alsofirstder, '|'
      write (stdout, '(1x,a46,10x,L8,13x,a1)') '|  Write data into a single file             :', geninterp%single_file, '|'
      write (stdout, '(1x,a78)') '*----------------------------------------------------------------------------*'
    endif

101 format(20x, a3, 2x, 3F11.6)

  end subroutine param_postw90_write

  subroutine param_pw90_dealloc
    use w90_io, only: io_error
    implicit none
    integer :: ierr

    call param_dealloc
    if (allocated(dos_data%project)) then
      deallocate (dos_data%project, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating dos_project in param_pw90_dealloc')
    endif
    if (allocated(fermi%energy_list)) then
      deallocate (fermi%energy_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating fermi_energy_list in param_pw90_dealloc')
    endif
    if (allocated(berry%kubo_freq_list)) then
      deallocate (berry%kubo_freq_list, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating kubo_freq_list in param_pw90_dealloc')
    endif
  end subroutine param_pw90_dealloc

  ! extra postw90 memory
  subroutine param_pw90_mem_estimate(mem_param, mem_bw)
    use w90_comms, only: on_root
    implicit none
    real(kind=dp), parameter :: size_log = 1.0_dp
    real(kind=dp), parameter :: size_int = 4.0_dp
    real(kind=dp), parameter :: size_real = 8.0_dp
    real(kind=dp), parameter :: size_cmplx = 16.0_dp
    real(kind=dp), intent(in) :: mem_param
    real(kind=dp), intent(inout) :: mem_bw
    integer :: NumPoints1, NumPoints2, NumPoints3, ndim
    real(kind=dp) :: TDF_exceeding_energy

    !if (ispostw90) then
    if (pw90_calcs%boltzwann) then
      if (pw90_common%spin_decomp) then
        ndim = 3
      else
        ndim = 1
      end if

      ! I set a big value to have a rough estimate
      TDF_exceeding_energy = 2._dp
      NumPoints1 = int(floor((boltz%temp_max - boltz%temp_min)/boltz%temp_step)) + 1 ! temperature array
      NumPoints2 = int(floor((boltz%mu_max - boltz%mu_min)/boltz%mu_step)) + 1  ! mu array
      NumPoints3 = int(floor((dis_data%win_max - dis_data%win_min &
                              + 2._dp*TDF_exceeding_energy)/boltz%tdf_energy_step)) + 1 ! tdfenergyarray
      mem_bw = mem_bw + NumPoints1*size_real                         !TempArray
      mem_bw = mem_bw + NumPoints1*size_real                         !KTArray
      mem_bw = mem_bw + NumPoints2*size_real                         !MuArray
      mem_bw = mem_bw + NumPoints3*size_real                         !TDFEnergyArray
      mem_bw = mem_bw + 6*NumPoints3*ndim*size_real                  !TDFArray
      mem_bw = mem_bw + 6*NumPoints3*size_real                       !IntegrandArray
      mem_bw = mem_bw + (9*4 + 6)*size_real
      !ElCondTimesSeebeckFP,ThisElCond,ElCondInverse,ThisSeebeck,ElCondTimesSeebeck
      mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !ElCond
      mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !Seebeck
      mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !ThermCond
      ! I put a upper bound here below (as if there was only 1 node), because I do not have any knowledge at this point
      ! of the number of processors, so I cannot have a correct estimate
      mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !LocalElCond
      mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !LocalSeebeck
      mem_bw = mem_bw + 6*NumPoints1*NumPoints2*size_real            !LocalThermCond

      mem_bw = mem_bw + num_wann*num_wann*size_cmplx                 !HH
      mem_bw = mem_bw + 3*num_wann*num_wann*size_cmplx               !delHH
      mem_bw = mem_bw + num_wann*num_wann*size_cmplx                 !UU
      mem_bw = mem_bw + 3*num_wann*size_real                         !del_eig
      mem_bw = mem_bw + num_wann*size_real                           !eig
      mem_bw = mem_bw + num_wann*size_real                           !levelspacing_k

      NumPoints1 = int(floor((boltz%dos_energy_max - boltz%dos_energy_min)/boltz%dos_energy_step)) + 1!dosnumpoints
      mem_bw = mem_bw + NumPoints1*size_real                         !DOS_EnergyArray
      mem_bw = mem_bw + 6*ndim*NumPoints3*size_real                  !TDF_k
      mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_k
      mem_bw = mem_bw + ndim*NumPoints1*size_real                    !DOS_all
    end if
    !end if

    if (on_root) then
      !if (ispostw90) then
      if (pw90_calcs%boltzwann) &
        write (stdout, '(1x,"|",24x,a15,f16.2,a,18x,"|")') 'BoltzWann:', (mem_param + mem_bw)/(1024**2), ' Mb'
      !end if
    endif

  end subroutine param_pw90_mem_estimate

  subroutine parameters_gyro_write_task(task, key, comment)
    use w90_io, only: stdout

    character(len=*), intent(in) :: task, key, comment
    character(len=42) :: comment1

    comment1 = comment
    if ((index(task, key) > 0) .or. (index(task, 'all') > 0)) then
      write (stdout, '(1x,a2,a42,a2,10x,a8,13x,a1)') '| ', comment1, ' :', '       T', '|'
    else
      write (stdout, '(1x,a2,a42,a2,10x,a8,13x,a1)') '| ', comment1, ' :', '       F', '|'
    endif
  end subroutine parameters_gyro_write_task

end module pw90_param_methods
