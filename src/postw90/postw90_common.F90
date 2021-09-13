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

module w90_postw90_common

!==============================================================================
!! This contains the common variables and procedures needed to set up a Wannier
!! interpolatation calculation for any physical property
!==============================================================================

  ! Should we remove this 'use w90_comms' and invoke in individual routines
  ! when needed?
  !
  !use w90_comms
  use w90_constants, only: dp

  implicit none

  private

  public :: pw90common_wanint_setup, pw90common_wanint_get_kpoint_file, pw90common_wanint_param_dist
  public :: pw90common_wanint_data_dist, pw90common_get_occ
  public :: pw90common_fourier_R_to_k, pw90common_fourier_R_to_k_new, pw90common_fourier_R_to_k_vec
  public :: kpoint_dist_type, wigner_seitz_type
  public :: pw90common_kmesh_spacing
  public :: pw90common_fourier_R_to_k_new_second_d, pw90common_fourier_R_to_k_new_second_d_TB_conv, &
            pw90common_fourier_R_to_k_vec_dadb, pw90common_fourier_R_to_k_vec_dadb_TB_conv

! AAM PROBABLY REMOVE THIS
  ! This 'save' statement could probably be ommited, since this module
  ! is USEd by the main program 'wannier_parint'
  !
  save

! AAM REMOVE THIS
  ! Default accessibility is PUBLIC
  !
!  private :: wigner_seitz
!
!  private :: kmesh_spacing_singleinteger, kmesh_spacing_mesh

  interface pw90common_kmesh_spacing
    module procedure kmesh_spacing_singleinteger
    module procedure kmesh_spacing_mesh
  end interface pw90common_kmesh_spacing

  ! Parameters describing the direct lattice points R on a
  ! Wigner-Seitz supercell
  !
  type wigner_seitz_type
    integer, allocatable       :: irvec(:, :)
    real(kind=dp), allocatable :: crvec(:, :)
    integer, allocatable       :: ndegen(:)
    integer                    :: nrpts
    integer                    :: rpt_origin
  end type wigner_seitz_type

  type kpoint_dist_type ! kpoints from file
    integer                       :: max_int_kpts_on_node, num_int_kpts
    integer, allocatable          :: num_int_kpts_on_node(:)
    real(kind=dp), allocatable    :: int_kpts(:, :), weight(:)
  end type kpoint_dist_type
  !complex(kind=dp), allocatable :: v_matrix(:, :, :)
  !real(kind=dp), public :: cell_volume

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       !
  !===========================================================!

  ! Public procedures have names starting with wanint_

  subroutine pw90common_wanint_setup(num_wann, verbose, real_lattice, mp_grid, effective_model, &
                                     ws_vec, stdout, seedname, world)
    !! Setup data ready for interpolation
    use w90_constants, only: dp !, cmplx_0
    use w90_io, only: io_error, io_file_unit
    !use w90_utility, only: utility_cart_to_frac
    use w90_param_types, only: print_output_type
    use w90_comms, only: mpirank, w90commtype, comms_bcast

    integer, intent(in) :: num_wann
    type(print_output_type), intent(in) :: verbose
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    integer, intent(in) :: mp_grid(3)
    type(wigner_seitz_type), intent(inout) :: ws_vec
    integer, intent(in) :: stdout
    logical, intent(in) :: effective_model
    character(len=50), intent(in)  :: seedname
    type(w90commtype), intent(in) :: world

    integer :: ierr, ir, file_unit, num_wann_loc
    logical :: on_root = .false.

    if (mpirank(world) == 0) on_root = .true.

    ! Find nrpts, the number of points in the Wigner-Seitz cell
    !
    if (effective_model) then
      if (on_root) then
        ! nrpts is read from file, together with num_wann
        file_unit = io_file_unit()
        open (file_unit, file=trim(seedname)//'_HH_R.dat', form='formatted', &
              status='old', err=101)
        read (file_unit, *) !header
        read (file_unit, *) num_wann_loc
        if (num_wann_loc /= num_wann) &
          call io_error('Inconsistent values of num_wann in ' &
                        //trim(seedname)//'_HH_R.dat and '//trim(seedname)//'.win', stdout, seedname)
        read (file_unit, *) ws_vec%nrpts
        close (file_unit)
      endif
      call comms_bcast(ws_vec%nrpts, 1, stdout, seedname, world)
    else
      call wigner_seitz(verbose, real_lattice, mp_grid, ws_vec, stdout, seedname, .true., world)
    endif

    ! Now can allocate several arrays
    !
    allocate (ws_vec%irvec(3, ws_vec%nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating irvec in pw90common_wanint_setup', stdout, seedname)
    ws_vec%irvec = 0
    allocate (ws_vec%crvec(3, ws_vec%nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating crvec in pw90common_wanint_setup', stdout, seedname)
    ws_vec%crvec = 0.0_dp
    allocate (ws_vec%ndegen(ws_vec%nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ndegen in pw90common_wanint_setup', stdout, seedname)
    ws_vec%ndegen = 0
    !
    ! Also rpt_origin, so that when effective_model=.true it is not
    ! passed to get_HH_R without being initialized.
    ws_vec%rpt_origin = 0

    ! If effective_model, this is done in get_HH_R
    if (.not. effective_model) then
      !
      ! Set up the lattice vectors on the Wigner-Seitz supercell
      ! where the Wannier functions live
      !
      call wigner_seitz(verbose, real_lattice, mp_grid, ws_vec, stdout, seedname, &
                        .false., world)
      !
      ! Convert from reduced to Cartesian coordinates
      !
      do ir = 1, ws_vec%nrpts
        ! Note that 'real_lattice' stores the lattice vectors as *rows*
        ws_vec%crvec(:, ir) = matmul(transpose(real_lattice), ws_vec%irvec(:, ir))
      end do
    endif

    return

101 call io_error('Error in pw90common_wanint_setup: problem opening file '// &
                  trim(seedname)//'_HH_R.dat', stdout, seedname)

  end subroutine pw90common_wanint_setup

  !===========================================================!
  subroutine pw90common_wanint_get_kpoint_file(kpoints, stdout, seedname, world)
    !===========================================================!
    !                                                           !
    !! read kpoints from kpoint.dat and distribute
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp
    use w90_io, only: io_error, io_file_unit, io_date, io_time, io_stopwatch
    use w90_comms, only: mpirank, mpisize, w90commtype, comms_send, comms_recv, comms_bcast

    ! arguments
    type(kpoint_dist_type), intent(inout) :: kpoints
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90commtype), intent(in) :: world

    ! local variables
    integer :: loop_nodes, loop_kpt, i, ierr, my_node_id, num_nodes, k_unit
    real(kind=dp) :: sum
    logical :: on_root = .false.

    my_node_id = mpirank(world)
    num_nodes = mpisize(world)

    if (my_node_id == 0) on_root = .true.

    k_unit = io_file_unit()
    if (on_root) then
      open (unit=k_unit, file='kpoint.dat', status='old', form='formatted', err=106)
      read (k_unit, *) kpoints%num_int_kpts
    end if
    call comms_bcast(kpoints%num_int_kpts, 1, stdout, seedname, world)

    allocate (kpoints%num_int_kpts_on_node(0:num_nodes - 1))
    kpoints%num_int_kpts_on_node(:) = kpoints%num_int_kpts/num_nodes
    kpoints%max_int_kpts_on_node = kpoints%num_int_kpts &
                                   - (num_nodes - 1)*(kpoints%num_int_kpts/num_nodes)
    kpoints%num_int_kpts_on_node(0) = kpoints%max_int_kpts_on_node
!    if(my_node_id < num_int_kpts- num_int_kpts_on_node*num_nodes)  num_int_kpts_on_node= num_int_kpts_on_node+1

    allocate (kpoints%int_kpts(3, kpoints%max_int_kpts_on_node), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating max_int_kpts_on_node in param_read_um', stdout, seedname)
    kpoints%int_kpts = 0.0_dp
    allocate (kpoints%weight(kpoints%max_int_kpts_on_node), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating weight in param_read_um', stdout, seedname)
    kpoints%weight = 0.0_dp

    sum = 0.0_dp
    if (on_root) then
      do loop_nodes = 1, num_nodes - 1
        do loop_kpt = 1, kpoints%num_int_kpts_on_node(loop_nodes)
          read (k_unit, *) (kpoints%int_kpts(i, loop_kpt), i=1, 3), kpoints%weight(loop_kpt)
          sum = sum + kpoints%weight(loop_kpt)
        end do

        call comms_send(kpoints%int_kpts(1, 1), 3*kpoints%num_int_kpts_on_node(loop_nodes), &
                        loop_nodes, stdout, seedname, world)
        call comms_send(kpoints%weight(1), kpoints%num_int_kpts_on_node(loop_nodes), loop_nodes, &
                        stdout, seedname, world)
      end do
      do loop_kpt = 1, kpoints%num_int_kpts_on_node(0)
        read (k_unit, *) (kpoints%int_kpts(i, loop_kpt), i=1, 3), kpoints%weight(loop_kpt)
        sum = sum + kpoints%weight(loop_kpt)
      end do
!       print*,'rsum',sum
    end if

    if (.not. on_root) then
      call comms_recv(kpoints%int_kpts(1, 1), 3*kpoints%num_int_kpts_on_node(my_node_id), 0, &
                      stdout, seedname, world)
      call comms_recv(kpoints%weight(1), kpoints%num_int_kpts_on_node(my_node_id), 0, &
                      stdout, seedname, world)
    end if

    return

106 call io_error('Error: Problem opening file kpoint.dat in pw90common_wanint_get_kpoint_file', stdout, seedname)

  end subroutine pw90common_wanint_get_kpoint_file

  !===========================================================!
  subroutine pw90common_wanint_param_dist(verbose, rs_region, kmesh_info, kpt_latt, num_kpts, &
                                          dis_window, system, fermi_energy_list, num_bands, &
                                          num_wann, eigval, mp_grid, real_lattice, &
                                          pw90_calcs, scissors_shift, effective_model, pw90_spin, &
                                          pw90_ham, kpath, kslice, dos_data, berry, spin_hall, &
                                          gyrotropic, geninterp, boltz, eig_found, stdout, &
                                          seedname, world)
    !===========================================================!
    !                                                           !
    !! distribute the parameters across processors
    !! NOTE: we only send the ones postw90 uses, not all in w90
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp !, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, io_date, io_time, &
      io_stopwatch
    use w90_comms, only: mpirank, w90commtype, comms_bcast
    use w90_param_types
    use pw90_parameters, only: pw90_calculation_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, pw90_kpath_mod_type, pw90_kslice_mod_type, pw90_dos_mod_type, &
      pw90_berry_mod_type, pw90_spin_hall_type, pw90_gyrotropic_type, pw90_geninterp_mod_type, &
      pw90_boltzwann_type

    type(print_output_type), intent(inout) :: verbose
    type(ws_region_type), intent(inout) :: rs_region
    type(w90_system_type), intent(inout) :: system
    type(kmesh_info_type), intent(inout) :: kmesh_info
    real(kind=dp), allocatable, intent(inout) :: kpt_latt(:, :)
    integer, intent(inout) :: num_kpts
    type(dis_manifold_type), intent(inout) :: dis_window
    real(kind=dp), allocatable, intent(inout) :: fermi_energy_list(:)
    integer, intent(inout) :: num_bands
    integer, intent(inout) :: num_wann
    real(kind=dp), allocatable, intent(inout) :: eigval(:, :)
    integer, intent(inout) :: mp_grid(3)
    real(kind=dp), intent(inout) :: real_lattice(3, 3)
    real(kind=dp), intent(inout) :: scissors_shift
    type(pw90_calculation_type), intent(inout) :: pw90_calcs
    type(pw90_spin_mod_type), intent(inout) :: pw90_spin
    type(pw90_band_deriv_degen_type), intent(inout) :: pw90_ham
    type(pw90_kpath_mod_type), intent(inout) :: kpath
    type(pw90_kslice_mod_type), intent(inout) :: kslice
    type(pw90_dos_mod_type), intent(inout) :: dos_data
    type(pw90_berry_mod_type), intent(inout) :: berry
    type(pw90_spin_hall_type), intent(inout) :: spin_hall
    type(pw90_gyrotropic_type), intent(inout) :: gyrotropic
    type(pw90_geninterp_mod_type), intent(inout) :: geninterp
    type(pw90_boltzwann_type), intent(inout) :: boltz
    logical, intent(inout) :: eig_found
    logical, intent(inout) :: effective_model
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90commtype), intent(in) :: world

    integer :: ierr
    integer :: iprintroot
    integer :: fermi_n
    logical :: on_root = .false.

    if (mpirank(world) == 0) on_root = .true.

    call comms_bcast(effective_model, 1, stdout, seedname, world)
    call comms_bcast(eig_found, 1, stdout, seedname, world)

    if (.not. effective_model) then
      call comms_bcast(mp_grid(1), 3, stdout, seedname, world)
      call comms_bcast(num_kpts, 1, stdout, seedname, world)
      call comms_bcast(num_bands, 1, stdout, seedname, world)
    endif
    call comms_bcast(num_wann, 1, stdout, seedname, world)
    call comms_bcast(verbose%timing_level, 1, stdout, seedname, world)

    !______________________________________
    !JJ fixme maybe? not so pretty solution to setting iprint to zero on non-root processes
    iprintroot = verbose%iprint
    verbose%iprint = 0
    call comms_bcast(verbose%iprint, 1, stdout, seedname, world)
    if (on_root) verbose%iprint = iprintroot
    !______________________________________

    call comms_bcast(rs_region%ws_distance_tol, 1, stdout, seedname, world)
    call comms_bcast(rs_region%ws_search_size(1), 3, stdout, seedname, world)
!    call comms_bcast(num_atoms,1)   ! Ivo: not used in postw90, right?
!    call comms_bcast(num_species,1) ! Ivo: not used in postw90, right?
    call comms_bcast(real_lattice(1, 1), 9, stdout, seedname, world)
    !call comms_bcast(recip_lattice(1, 1), 9, stdout, seedname, world)
    !call comms_bcast(real_metric(1, 1), 9)
    !call comms_bcast(recip_metric(1, 1), 9)
    !call comms_bcast(cell_volume, 1, stdout, seedname, world)
    call comms_bcast(dos_data%energy_step, 1, stdout, seedname, world)
    call comms_bcast(dos_data%smearing%use_adaptive, 1, stdout, seedname, world)
    call comms_bcast(dos_data%smearing%type_index, 1, stdout, seedname, world)
    call comms_bcast(dos_data%kmesh%spacing, 1, stdout, seedname, world)
    call comms_bcast(dos_data%kmesh%mesh(1), 3, stdout, seedname, world)
    call comms_bcast(dos_data%smearing%adaptive_max_width, 1, stdout, seedname, world)
    call comms_bcast(dos_data%smearing%fixed_width, 1, stdout, seedname, world)
    call comms_bcast(dos_data%smearing%adaptive_prefactor, 1, stdout, seedname, world)
    call comms_bcast(dos_data%num_project, 1, stdout, seedname, world)

    call comms_bcast(pw90_calcs%berry, 1, stdout, seedname, world)
    call comms_bcast(berry%task, len(berry%task), stdout, seedname, world)
    call comms_bcast(berry%kmesh%spacing, 1, stdout, seedname, world)
    call comms_bcast(berry%kmesh%mesh(1), 3, stdout, seedname, world)
    call comms_bcast(berry%curv_adpt_kmesh, 1, stdout, seedname, world)
    call comms_bcast(berry%curv_adpt_kmesh_thresh, 1, stdout, seedname, world)
    call comms_bcast(berry%curv_unit, len(berry%curv_unit), stdout, seedname, world)

! Tsirkin
    call comms_bcast(pw90_calcs%gyrotropic, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%task, len(gyrotropic%task), stdout, seedname, world)
    call comms_bcast(gyrotropic%kmesh%spacing, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%kmesh%mesh(1), 3, stdout, seedname, world)
    call comms_bcast(gyrotropic%smr_fixed_en_width, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%smr_index, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%eigval_max, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%nfreq, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%degen_thresh, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%num_bands, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%box(1, 1), 9, stdout, seedname, world)
    call comms_bcast(gyrotropic%box_corner(1), 3, stdout, seedname, world)
    call comms_bcast(gyrotropic%smr_max_arg, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%smr_fixed_en_width, 1, stdout, seedname, world)
    call comms_bcast(gyrotropic%smr_index, 1, stdout, seedname, world)

    call comms_bcast(system%spinors, 1, stdout, seedname, world)

    call comms_bcast(spin_hall%freq_scan, 1, stdout, seedname, world)
    call comms_bcast(spin_hall%alpha, 1, stdout, seedname, world)
    call comms_bcast(spin_hall%beta, 1, stdout, seedname, world)
    call comms_bcast(spin_hall%gamma, 1, stdout, seedname, world)
    call comms_bcast(spin_hall%bandshift, 1, stdout, seedname, world)
    call comms_bcast(spin_hall%bandshift_firstband, 1, stdout, seedname, world)
    call comms_bcast(spin_hall%bandshift_energyshift, 1, stdout, seedname, world)

    call comms_bcast(berry%kubo_smearing%use_adaptive, 1, stdout, seedname, world)
    call comms_bcast(berry%kubo_smearing%adaptive_prefactor, 1, stdout, seedname, world)
    call comms_bcast(berry%kubo_smearing%adaptive_max_width, 1, stdout, seedname, world)
    call comms_bcast(berry%kubo_smearing%fixed_width, 1, stdout, seedname, world)
    call comms_bcast(berry%kubo_smearing%type_index, 1, stdout, seedname, world)
    call comms_bcast(berry%kubo_eigval_max, 1, stdout, seedname, world)
    call comms_bcast(berry%kubo_nfreq, 1, stdout, seedname, world)
    fermi_n = 0
    if (on_root) then
      if (allocated(fermi_energy_list)) fermi_n = size(fermi_energy_list)
    endif
    call comms_bcast(fermi_n, 1, stdout, seedname, world)
    call comms_bcast(dos_data%energy_min, 1, stdout, seedname, world)
    call comms_bcast(dos_data%energy_max, 1, stdout, seedname, world)
    call comms_bcast(pw90_spin%kmesh%spacing, 1, stdout, seedname, world)
    call comms_bcast(pw90_spin%kmesh%mesh(1), 3, stdout, seedname, world)
    call comms_bcast(berry%wanint_kpoint_file, 1, stdout, seedname, world)
    call comms_bcast(dis_window%win_min, 1, stdout, seedname, world)
    call comms_bcast(dis_window%win_max, 1, stdout, seedname, world)
    call comms_bcast(berry%sc_eta, 1, stdout, seedname, world)
    call comms_bcast(berry%sc_w_thr, 1, stdout, seedname, world)
    call comms_bcast(berry%sc_phase_conv, 1, stdout, seedname, world)
! ----------------------------------------------
!
! New input variables in development
!
    !call comms_bcast(verbose%devel_flag, len(verbose%devel_flag), stdout, seedname, world)
    call comms_bcast(pw90_calcs%spin_moment, 1, stdout, seedname, world)
    call comms_bcast(pw90_spin%axis_polar, 1, stdout, seedname, world)
    call comms_bcast(pw90_spin%axis_azimuth, 1, stdout, seedname, world)
    call comms_bcast(pw90_calcs%spin_decomp, 1, stdout, seedname, world)
    call comms_bcast(pw90_ham%use_degen_pert, 1, stdout, seedname, world)
    call comms_bcast(pw90_ham%degen_thr, 1, stdout, seedname, world)
    call comms_bcast(system%num_valence_bands, 1, stdout, seedname, world)
    call comms_bcast(pw90_calcs%dos, 1, stdout, seedname, world)
    call comms_bcast(dos_data%task, len(dos_data%task), stdout, seedname, world)
    call comms_bcast(pw90_calcs%kpath, 1, stdout, seedname, world)
    call comms_bcast(kpath%task, len(kpath%task), stdout, seedname, world)
    call comms_bcast(kpath%bands_colour, len(kpath%bands_colour), stdout, seedname, world)
    call comms_bcast(pw90_calcs%kslice, 1, stdout, seedname, world)
    call comms_bcast(kslice%task, len(kslice%task), stdout, seedname, world)
    call comms_bcast(kslice%corner(1), 3, stdout, seedname, world)
    call comms_bcast(kslice%b1(1), 3, stdout, seedname, world)
    call comms_bcast(kslice%b2(1), 3, stdout, seedname, world)
    call comms_bcast(kslice%kmesh2d(1), 2, stdout, seedname, world)
    call comms_bcast(kslice%fermi_lines_colour, len(kslice%fermi_lines_colour), stdout, seedname, &
                     world)
    call comms_bcast(berry%transl_inv, 1, stdout, seedname, world)
    call comms_bcast(system%num_elec_per_state, 1, stdout, seedname, world)
    call comms_bcast(scissors_shift, 1, stdout, seedname, world)
    !
    ! Do these have to be broadcasted? (Plots done on root node only)
    !
!    call comms_bcast(bands_num_points,1)
!    call comms_bcast(bands_num_spec_points,1)
!    if(allocated(bands_spec_points)) &
!         call comms_bcast(bands_spec_points(1,1),3*bands_num_spec_points)
!    if(allocated(bands_label)) &
!         call comms_bcast(bands_label(:),len(bands_label(1))*bands_num_spec_points)
! ----------------------------------------------
    call comms_bcast(pw90_calcs%geninterp, 1, stdout, seedname, world)
    call comms_bcast(geninterp%alsofirstder, 1, stdout, seedname, world)
    call comms_bcast(geninterp%single_file, 1, stdout, seedname, world)
    ! [gp-begin, Apr 12, 2012]
    ! BoltzWann variables
    call comms_bcast(pw90_calcs%boltzwann, 1, stdout, seedname, world)
    call comms_bcast(boltz%calc_also_dos, 1, stdout, seedname, world)
    call comms_bcast(boltz%dir_num_2d, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_energy_step, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_energy_min, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_energy_max, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_smearing%use_adaptive, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_smearing%fixed_width, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_smearing%adaptive_prefactor, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_smearing%adaptive_max_width, 1, stdout, seedname, world)
    call comms_bcast(boltz%mu_min, 1, stdout, seedname, world)
    call comms_bcast(boltz%mu_max, 1, stdout, seedname, world)
    call comms_bcast(boltz%mu_step, 1, stdout, seedname, world)
    call comms_bcast(boltz%temp_min, 1, stdout, seedname, world)
    call comms_bcast(boltz%temp_max, 1, stdout, seedname, world)
    call comms_bcast(boltz%temp_step, 1, stdout, seedname, world)
    call comms_bcast(boltz%kmesh%spacing, 1, stdout, seedname, world)
    call comms_bcast(boltz%kmesh%mesh(1), 3, stdout, seedname, world)
    call comms_bcast(boltz%tdf_energy_step, 1, stdout, seedname, world)
    call comms_bcast(boltz%relax_time, 1, stdout, seedname, world)
    call comms_bcast(boltz%TDF_smr_fixed_en_width, 1, stdout, seedname, world)
    call comms_bcast(boltz%TDF_smr_index, 1, stdout, seedname, world)
    call comms_bcast(boltz%dos_smearing%type_index, 1, stdout, seedname, world)
    call comms_bcast(boltz%bandshift, 1, stdout, seedname, world)
    call comms_bcast(boltz%bandshift_firstband, 1, stdout, seedname, world)
    call comms_bcast(boltz%bandshift_energyshift, 1, stdout, seedname, world)
    ! [gp-end]
    call comms_bcast(rs_region%use_ws_distance, 1, stdout, seedname, world)

    ! These variables are different from the ones above in that they are
    ! allocatable, and in param_read they were allocated on the root node only
    !
    if (.not. on_root) then
      allocate (fermi_energy_list(fermi_n), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating fermi_energy_read in postw90_param_dist', stdout, seedname)
      allocate (berry%kubo_freq_list(berry%kubo_nfreq), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating kubo_freq_list in postw90_param_dist', stdout, seedname)

      allocate (gyrotropic%band_list(gyrotropic%num_bands), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating gyrotropic_band_list in postw90_param_dist', stdout, seedname)

      allocate (gyrotropic%freq_list(gyrotropic%nfreq), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating gyrotropic_freq_list in postw90_param_dist', stdout, seedname)

      allocate (dos_data%project(dos_data%num_project), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating dos_project in postw90_param_dist', stdout, seedname)
      if (.not. effective_model) then
        if (eig_found) then
          allocate (eigval(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating eigval in postw90_param_dist', stdout, seedname)
        end if
        allocate (kpt_latt(3, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error allocating kpt_latt in postw90_param_dist', stdout, seedname)
      endif
    end if

    if (fermi_n > 0) call comms_bcast(fermi_energy_list(1), fermi_n, stdout, seedname, world)
    call comms_bcast(gyrotropic%freq_list(1), gyrotropic%nfreq, stdout, seedname, world)
    call comms_bcast(gyrotropic%band_list(1), gyrotropic%num_bands, stdout, seedname, world)
    call comms_bcast(berry%kubo_freq_list(1), berry%kubo_nfreq, stdout, seedname, world)
    call comms_bcast(dos_data%project(1), dos_data%num_project, stdout, seedname, world)
    if (.not. effective_model) then
      if (eig_found) then
        call comms_bcast(eigval(1, 1), num_bands*num_kpts, stdout, seedname, world)
      end if
      call comms_bcast(kpt_latt(1, 1), 3*num_kpts, stdout, seedname, world)
    endif

    ! kmesh: only nntot,wb, and bk are needed to evaluate the WF matrix
    ! elements of the position operator in reciprocal space. For the
    ! extra matrix elements entering the orbital magnetization, also
    ! need nnlist. In principle could only broadcast those four variables

    if (.not. effective_model) then

      call comms_bcast(kmesh_info%nnh, 1, stdout, seedname, world)
      call comms_bcast(kmesh_info%nntot, 1, stdout, seedname, world)

      if (.not. on_root) then
        allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nnlist in pw90common_wanint_param_dist', stdout, seedname)
        allocate (kmesh_info%neigh(num_kpts, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating neigh in pw90common_wanint_param_dist', stdout, seedname)
        allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nncell in pw90common_wanint_param_dist', stdout, seedname)
        allocate (kmesh_info%wb(kmesh_info%nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating wb in pw90common_wanint_param_dist', stdout, seedname)
        allocate (kmesh_info%bka(3, kmesh_info%nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bka in pw90common_wanint_param_dist', stdout, seedname)
        allocate (kmesh_info%bk(3, kmesh_info%nntot, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bk in pw90common_wanint_param_dist', stdout, seedname)
      end if

      call comms_bcast(kmesh_info%nnlist(1, 1), num_kpts*kmesh_info%nntot, stdout, seedname, world)
      call comms_bcast(kmesh_info%neigh(1, 1), num_kpts*kmesh_info%nntot/2, stdout, seedname, world)
      call comms_bcast(kmesh_info%nncell(1, 1, 1), 3*num_kpts*kmesh_info%nntot, stdout, seedname, world)
      call comms_bcast(kmesh_info%wb(1), kmesh_info%nntot, stdout, seedname, world)
      call comms_bcast(kmesh_info%bka(1, 1), 3*kmesh_info%nntot/2, stdout, seedname, world)
      call comms_bcast(kmesh_info%bk(1, 1, 1), 3*kmesh_info%nntot*num_kpts, stdout, seedname, world)

    endif

  end subroutine pw90common_wanint_param_dist

  !===========================================================!
  subroutine pw90common_wanint_data_dist(num_wann, num_kpts, num_bands, u_matrix_opt, u_matrix, &
                                         dis_window, wann_data, scissors_shift, v_matrix, &
                                         num_valence_bands, have_disentangled, stdout, seedname, &
                                         world)
    !===========================================================!
    !                                                           !
    !! Distribute the um and chk files
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp, cmplx_0 !, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, &
      io_date, io_time, io_stopwatch
    use w90_param_types, only: dis_manifold_type, wannier_data_type
    use w90_comms, only: w90commtype, mpirank, comms_bcast

    implicit none
    integer, intent(in) :: num_wann, num_kpts, num_bands
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :), u_matrix(:, :, :)
    type(dis_manifold_type), intent(inout) :: dis_window
    type(wannier_data_type), intent(inout) :: wann_data
    complex(kind=dp), allocatable :: v_matrix(:, :, :)
    integer, intent(in) :: num_valence_bands
    logical, intent(inout) :: have_disentangled
    integer, intent(in) :: stdout
    real(kind=dp), intent(in) :: scissors_shift

    character(len=50), intent(in)  :: seedname
    type(w90commtype), intent(in) :: world

    integer :: ierr, loop_kpt, m, i, j
    logical :: on_root = .false.

    if (mpirank(world) == 0) on_root = .true.

    if (.not. on_root) then
      ! wannier_centres is allocated in param_read, so only on root node
      ! It is then read in param_read_chpkt
      ! Therefore, now we need to allocate it on all nodes, and then broadcast it
      allocate (wann_data%centres(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in pw90common_wanint_data_dist', stdout, seedname)
    end if
    call comms_bcast(wann_data%centres(1, 1), 3*num_wann, stdout, seedname, world)

    ! -------------------
    ! Ivo: added 8april11
    ! -------------------
    !
    ! Calculate the matrix that describes the combined effect of
    ! disentanglement and maximal localization. This is the combination
    ! that is most often needed for interpolation purposes
    !
    ! Allocate on all nodes
    allocate (v_matrix(num_bands, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) &
      call io_error('Error allocating v_matrix in pw90common_wanint_data_dist', stdout, seedname)
    ! u_matrix and u_matrix_opt are stored on root only
    if (on_root) then
      if (.not. have_disentangled) then
        v_matrix(1:num_wann, :, :) = u_matrix(1:num_wann, :, :)
      else
        v_matrix = cmplx_0
        do loop_kpt = 1, num_kpts
          do j = 1, num_wann
            do m = 1, dis_window%ndimwin(loop_kpt)
              do i = 1, num_wann
                v_matrix(m, j, loop_kpt) = v_matrix(m, j, loop_kpt) &
                                           + u_matrix_opt(m, i, loop_kpt)*u_matrix(i, j, loop_kpt)
              enddo
            enddo
          enddo
        enddo
      endif
      if (allocated(u_matrix_opt)) deallocate (u_matrix_opt)
      if (.not. (num_valence_bands > 0 .and. abs(scissors_shift) > 1.0e-7_dp)) then
        if (allocated(u_matrix)) deallocate (u_matrix)
      endif
    endif
    call comms_bcast(v_matrix(1, 1, 1), num_bands*num_wann*num_kpts, stdout, seedname, world)

    if (num_valence_bands > 0 .and. abs(scissors_shift) > 1.0e-7_dp) then
    if (.not. on_root .and. .not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating u_matrix in pw90common_wanint_data_dist', stdout, seedname)
    endif
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, stdout, seedname, world)
    endif

!    if (.not.on_root .and. .not.allocated(m_matrix)) then
!       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!       if (ierr/=0)&
!            call io_error('Error allocating m_matrix in pw90common_wanint_data_dist')
!    endif
!    call comms_bcast(m_matrix(1,1,1,1),num_wann*num_wann*nntot*num_kpts)

    call comms_bcast(have_disentangled, 1, stdout, seedname, world)

    if (have_disentangled) then
      if (.not. on_root) then

        ! Do we really need these 'if not allocated'? Didn't use them for
        ! eigval and kpt_latt above!

!          if (.not.allocated(u_matrix_opt)) then
!             allocate(u_matrix_opt(num_bands,num_wann,num_kpts),stat=ierr)
!             if (ierr/=0)&
!              call io_error('Error allocating u_matrix_opt in pw90common_wanint_data_dist')
!          endif

        if (.not. allocated(dis_window%lwindow)) then
          allocate (dis_window%lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating lwindow in pw90common_wanint_data_dist', stdout, seedname)
        endif

        if (.not. allocated(dis_window%ndimwin)) then
          allocate (dis_window%ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating ndimwin in pw90common_wanint_data_dist', stdout, seedname)
        endif

      end if

!       call comms_bcast(u_matrix_opt(1,1,1),num_bands*num_wann*num_kpts)
      call comms_bcast(dis_window%lwindow(1, 1), num_bands*num_kpts, stdout, seedname, world)
      call comms_bcast(dis_window%ndimwin(1), num_kpts, stdout, seedname, world)
    end if

  end subroutine pw90common_wanint_data_dist

!=======================================================================

  subroutine pw90common_get_occ(ef, eig, occ, num_wann)
    !! Compute the electronic occupancy

    use w90_constants, only: dp !,eps7
    !use w90_parameters, only: num_wann !,smear_temp
!    use w90_constants, only    : elem_charge_SI,k_B_SI

    ! Arguments
    !
    integer, intent(in) :: num_wann

    real(kind=dp), intent(in)  :: eig(num_wann)
    !! Eigenvalues
    real(kind=dp), intent(in)  :: ef
    !! Fermi level
    real(kind=dp), intent(out) :: occ(num_wann)
    !! Occupancy of states

    ! Misc/Dummy
    !
    integer       :: i
!    real(kind=dp) :: kt

    ! State occupancies
    !
!    if(smear_temp < eps7) then
    !
    ! Use a step function occupancy (T=0)
    !
    occ(:) = 0.0_dp
    do i = 1, num_wann
      if (eig(i) < ef) occ(i) = 1.0_dp
    end do
!    else
    !
    ! Use a Fermi-Dirac occupancy (T=smear_temp, in Kelvin)
    !
    ! k_B.T in electron-volts
    !
!       kt=k_B_SI*smear_temp/elem_charge_SI
!       do i=1,num_wann
!          occ(i)=1.0_dp/(exp((eig(i)-ef)/kt)+1.0_dp)
!       end do
!    end if

  end subroutine pw90common_get_occ

!=======================================================================

  function kmesh_spacing_singleinteger(num_points, recip_lattice)

    !! Set up the value of the interpolation mesh spacing, needed for
    !! adaptive smearing [see Eqs. (34-35) YWVS07]. Choose it as the largest of
    !! the three Delta_k's for each of the primitive translations b1, b2, and b3

    integer, intent(in) :: num_points
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp)       :: kmesh_spacing_singleinteger

    integer        :: i
    real(kind=dp) :: Delta_k_i(3)

    ! NOTE: The vectors b_i are stored as *rows* in recip_lattice (argh!).
    ! Hence I believe Jonathan's original code confused rows with columns
    ! when computing Delta_k, which he called 'rspace'
    ! (See my e-mail of 20Sept07)
    !
    do i = 1, 3
      Delta_k_i(i) = sqrt(dot_product(recip_lattice(i, :), recip_lattice(i, :)))/num_points
    end do
    kmesh_spacing_singleinteger = maxval(Delta_k_i)

  end function kmesh_spacing_singleinteger
  function kmesh_spacing_mesh(mesh, recip_lattice)
    !! Same as kmesh_spacing_singleinteger, but for a kmesh with three
    !! different mesh samplings along the three directions

    integer, dimension(3), intent(in) :: mesh
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp) :: kmesh_spacing_mesh

    integer        :: i
    real(kind=dp) :: Delta_k_i(3)

    do i = 1, 3
      Delta_k_i(i) = sqrt(dot_product(recip_lattice(i, :), recip_lattice(i, :)))/mesh(i)
    end do
    kmesh_spacing_mesh = maxval(Delta_k_i)

  end function kmesh_spacing_mesh
  !
  !=========================================================!
  subroutine pw90common_fourier_R_to_k(rs_region, wann_data, ws_distance, ws_vec, OO, OO_R, kpt, &
                                       real_lattice, mp_grid, alpha, num_wann, seedname, stdout)
    !=========================================================!
    !                                                         !
    !! For alpha=0:
    !! O_ij(R) --> O_ij(k) = sum_R e^{+ik.R}*O_ij(R)
    !!
    !! For alpha=1,2,3:
    !! sum_R [cmplx_i*R_alpha*e^{+ik.R}*O_ij(R)]
    !! where R_alpha is a Cartesian component of R
    !! ***REMOVE EVENTUALLY*** (replace with pw90common_fourier_R_to_k_new)

    !                                                         !
    !=========================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_param_types, only: wannier_data_type, ws_region_type
    use w90_ws_distance, only: ws_translate_dist, ws_distance_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: rs_region
    type(wannier_data_type), intent(in) :: wann_data
    type(wigner_seitz_type), intent(in) :: ws_vec
    type(ws_distance_type), intent(inout) :: ws_distance

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout
    integer, intent(in) :: alpha

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), intent(out) :: OO(:, :)

    character(len=50), intent(in)  :: seedname

    ! local variables
    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

!
!  subroutine ws_translate_dist(ws_distance_tol, ws_search_size, num_wann, &
!                               wannier_centres, real_lattice, recip_lattice, iprint, mp_grid, nrpts, &
!                               irvec, force_recompute)

    if (rs_region%use_ws_distance) then
      CALL ws_translate_dist(ws_distance, stdout, seedname, rs_region, num_wann, &
                             wann_data%centres, real_lattice, mp_grid, ws_vec%nrpts, ws_vec%irvec)
    endif

    OO(:, :) = cmplx_0
    do ir = 1, ws_vec%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (rs_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(ws_vec%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (alpha == 0) then
              OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
              OO(i, j) = OO(i, j) + cmplx_i*ws_distance%crdist(alpha, ideg, i, j, ir) &
                         *phase_fac*OO_R(i, j, ir)
            else
              stop 'wrong value of alpha in pw90common_fourier_R_to_k'
            endif
          enddo
        enddo
        enddo
      else
        ! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), ws_vec%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir), dp)
        if (alpha == 0) then
          OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
          OO(:, :) = OO(:, :) + &
                     cmplx_i*ws_vec%crvec(alpha, ir)*phase_fac*OO_R(:, :, ir)
        else
          stop 'wrong value of alpha in pw90common_fourier_R_to_k'
        endif
      endif

    enddo

  end subroutine pw90common_fourier_R_to_k

  ! ***NEW***
  !
  !=========================================================!
  subroutine pw90common_fourier_R_to_k_new(rs_region, wann_data, ws_distance, ws_vec, OO_R, kpt, &
                                           real_lattice, mp_grid, num_wann, seedname, stdout, &
                                           OO, OO_dx, OO_dy, OO_dz)
    !=======================================================!
    !                                                       !
    !! For OO:
    !! $$O_{ij}(k) = \sum_R e^{+ik.R}.O_{ij}(R)$$
    !! For $$OO_{dx,dy,dz}$$:
    !! $$\sum_R [i.R_{dx,dy,dz}.e^{+ik.R}.O_{ij}(R)]$$
    !! where R_{x,y,z} are the Cartesian components of R
    !                                                       !
    !=======================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_param_types, only: ws_region_type, wannier_data_type
    use w90_ws_distance, only: ws_translate_dist, ws_distance_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: rs_region
    type(wannier_data_type), intent(in) :: wann_data
    type(wigner_seitz_type), intent(in) :: ws_vec
    type(ws_distance_type), intent(inout) :: ws_distance

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO(:, :)
    complex(kind=dp), optional, intent(out) :: OO_dx(:, :)
    complex(kind=dp), optional, intent(out) :: OO_dy(:, :)
    complex(kind=dp), optional, intent(out) :: OO_dz(:, :)

    character(len=50), intent(in)  :: seedname

    ! local variables
    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (rs_region%use_ws_distance) CALL ws_translate_dist(ws_distance, stdout, seedname, &
                                                          rs_region, num_wann, &
                                                          wann_data%centres, real_lattice, &
                                                          mp_grid, ws_vec%nrpts, ws_vec%irvec)

    if (present(OO)) OO = cmplx_0
    if (present(OO_dx)) OO_dx = cmplx_0
    if (present(OO_dy)) OO_dy = cmplx_0
    if (present(OO_dz)) OO_dz = cmplx_0
    do ir = 1, ws_vec%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (rs_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(ws_vec%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_dx)) OO_dx(i, j) = OO_dx(i, j) + &
                                              cmplx_i*ws_distance%crdist(1, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
            if (present(OO_dy)) OO_dy(i, j) = OO_dy(i, j) + &
                                              cmplx_i*ws_distance%crdist(2, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
            if (present(OO_dz)) OO_dz(i, j) = OO_dz(i, j) + &
                                              cmplx_i*ws_distance%crdist(3, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), ws_vec%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir), dp)
        if (present(OO)) OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        if (present(OO_dx)) OO_dx(:, :) = OO_dx(:, :) + &
                                          cmplx_i*ws_vec%crvec(1, ir)*phase_fac*OO_R(:, :, ir)
        if (present(OO_dy)) OO_dy(:, :) = OO_dy(:, :) + &
                                          cmplx_i*ws_vec%crvec(2, ir)*phase_fac*OO_R(:, :, ir)
        if (present(OO_dz)) OO_dz(:, :) = OO_dz(:, :) + &
                                          cmplx_i*ws_vec%crvec(3, ir)*phase_fac*OO_R(:, :, ir)
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new

  !=========================================================!
  subroutine pw90common_fourier_R_to_k_new_second_d(kpt, OO_R, num_wann, rs_region, wann_data, &
                                                    real_lattice, mp_grid, ws_distance, ws_vec, &
                                                    stdout, seedname, OO, OO_da, OO_dadb)
    !=======================================================!
    !                                                       !
    !! For OO:
    !! $$O_{ij}(k) = \sum_R e^{+ik.R}.O_{ij}(R)$$
    !! For $$OO_{dx,dy,dz}$$:
    !! $$\sum_R [i.R_{dx,dy,dz}.e^{+ik.R}.O_{ij}(R)]$$
    !! where R_{x,y,z} are the Cartesian components of R
    !! For $$OO_{dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$-\sum_R [R_{dx1,dy1,dz1}.R_{dx2,dy2,dz2}.e^{+ik.R}.O_{ij}(R)]$$
    !! where R_{xi,yi,zi} are the Cartesian components of R
    !                                                       !
    !=======================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_param_types, only: ws_region_type, wannier_data_type
    use w90_ws_distance, only: ws_translate_dist, ws_distance_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: rs_region
    type(wannier_data_type), intent(in) :: wann_data
    type(wigner_seitz_type), intent(in) :: ws_vec
    type(ws_distance_type), intent(inout) :: ws_distance

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    character(len=50), intent(in)  :: seedname
    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO(:, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_dadb(:, :, :, :)

    ! local variables
    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (rs_region%use_ws_distance) CALL ws_translate_dist(ws_distance, stdout, seedname, &
                                                          rs_region, num_wann, &
                                                          wann_data%centres, real_lattice, &
                                                          mp_grid, ws_vec%nrpts, ws_vec%irvec)

    if (present(OO)) OO = cmplx_0
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, ws_vec%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (rs_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(ws_vec%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i*ws_distance%crdist(a, ideg, i, j, ir)* &
                                 phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) &
                                        - ws_distance%crdist(a, ideg, i, j, ir) &
                                        *ws_distance%crdist(b, ideg, i, j, ir) &
                                        *phase_fac*OO_R(i, j, ir)
                enddo
              enddo
            end if

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), ws_vec%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir), dp)
        if (present(OO)) OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        if (present(OO_da)) then
          do a = 1, 3
            OO_da(:, :, a) = OO_da(:, :, a) + cmplx_i*ws_vec%crvec(a, ir)*phase_fac*OO_R(:, :, ir)
          enddo
        endif
        if (present(OO_dadb)) then
          do a = 1, 3
            do b = 1, 3
              OO_dadb(:, :, a, b) = OO_dadb(:, :, a, b) - &
                                    ws_vec%crvec(a, ir)*ws_vec%crvec(b, ir)*phase_fac*OO_R(:, :, ir)
            enddo
          enddo
        end if
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new_second_d

  subroutine pw90common_fourier_R_to_k_new_second_d_TB_conv(kpt, OO_R, oo_a_R, num_wann, &
                                                            rs_region, wann_data, real_lattice, &
                                                            mp_grid, ws_distance, ws_vec, stdout, &
                                                            seedname, OO, OO_da, OO_dadb)
    !=======================================================!
    ! modified version of pw90common_fourier_R_to_k_new_second_d, includes wannier centres in
    ! the exponential inside the sum (so called TB convention)
    !
    !! For OO:
    !! $$O_{ij}(k) = \sum_R e^{+ik.(R+tau_ij)}.O_{ij}(R)$$
    !! For $$OO_{dx,dy,dz}$$:
    !! $$\sum_R [i.(R+tau_ij)_{dx,dy,dz}.e^{+ik.(R+tau_ij)}.O_{ij}(R)]$$
    !! where R_{x,y,z} are the Cartesian components of R
    !! For $$OO_{dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$-\sum_R [(R+tau_ij)_{dx1,dy1,dz1}.(R+tau_ij)_{dx2,dy2,dz2}.e^{+ik.(R+tau_ij)}.O_{ij}(R)]$$
    !! where {xi,yi,zi} denote the Cartesian components and
    !  tau_ij = tau_j - tau_i, being tau_i=<0i|r|0i> the individual wannier centres
    !=======================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_param_types, only: ws_region_type, wannier_data_type
    use w90_ws_distance, only: ws_translate_dist, ws_distance_type
    use w90_utility, only: utility_cart_to_frac, utility_inverse_mat

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: rs_region
    type(wannier_data_type), intent(in) :: wann_data
    type(wigner_seitz_type), intent(in) :: ws_vec
    type(ws_distance_type), intent(inout) :: ws_distance

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: oo_a_R(:, :, :, :)
    complex(kind=dp), intent(in) :: OO_R(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO(:, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_dadb(:, :, :, :)

    character(len=50), intent(in)  :: seedname

    ! local variables
    real(kind=dp) :: inv_lattice(3, 3)
    integer :: ir, i, j, ideg, a, b
    real(kind=dp) :: rdotk
    real(kind=dp) :: local_wannier_centres(3, num_wann), wannier_centres_frac(3, num_wann)
    real(kind=dp) :: r_sum(3)
    complex(kind=dp) :: phase_fac

    r_sum = 0.d0

    if (rs_region%use_ws_distance) CALL ws_translate_dist(ws_distance, stdout, seedname, &
                                                          rs_region, num_wann, &
                                                          wann_data%centres, real_lattice, &
                                                          mp_grid, ws_vec%nrpts, ws_vec%irvec)

    ! calculate wannier centres in cartesian
    local_wannier_centres(:, :) = 0.d0
    do j = 1, num_wann
      do ir = 1, ws_vec%nrpts
        if ((ws_vec%irvec(1, ir) .eq. 0) .and. (ws_vec%irvec(2, ir) .eq. 0) .and. &
            (ws_vec%irvec(3, ir) .eq. 0)) then
          local_wannier_centres(1, j) = real(oo_a_R(j, j, ir, 1))
          local_wannier_centres(2, j) = real(oo_a_R(j, j, ir, 2))
          local_wannier_centres(3, j) = real(oo_a_R(j, j, ir, 3))
        endif
      enddo
    enddo
    ! rotate wannier centres from cartesian to fractional coordinates
    wannier_centres_frac(:, :) = 0.d0
    call utility_inverse_mat(real_lattice, inv_lattice)
    do ir = 1, num_wann
      call utility_cart_to_frac(local_wannier_centres(:, ir), wannier_centres_frac(:, ir), &
                                inv_lattice)
    enddo

    if (present(OO)) OO = cmplx_0
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, ws_vec%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (rs_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir) + &
                                                   wannier_centres_frac(:, j) - wannier_centres_frac(:, i), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(ws_vec%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i* &
                                 (ws_distance%crdist(a, ideg, i, j, ir) &
                                  + local_wannier_centres(a, j) &
                                  - local_wannier_centres(a, i))*phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) - &
                                        (ws_distance%crdist(a, ideg, i, j, ir) &
                                         + local_wannier_centres(a, j) &
                                         - local_wannier_centres(a, i)) &
                                        *(ws_distance%crdist(b, ideg, i, j, ir) &
                                          + local_wannier_centres(b, j) - &
                                          local_wannier_centres(b, i))*phase_fac*OO_R(i, j, ir)
                enddo
              enddo
            end if

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        do j = 1, num_wann
          do i = 1, num_wann
            r_sum(:) = real(ws_vec%irvec(:, ir)) &
                       + wannier_centres_frac(:, j) - wannier_centres_frac(:, i)
            rdotk = twopi*dot_product(kpt(:), r_sum(:))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i* &
                                 (ws_vec%crvec(a, ir) + local_wannier_centres(a, j) - &
                                  local_wannier_centres(a, i))*phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = &
                    OO_dadb(i, j, a, b) - &
                    (ws_vec%crvec(a, ir) + local_wannier_centres(a, j) - local_wannier_centres(a, i))* &
                    (ws_vec%crvec(b, ir) + local_wannier_centres(b, j) - local_wannier_centres(b, i))* &
                    phase_fac*OO_R(i, j, ir)
                enddo
              enddo
            end if
          enddo
        enddo
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new_second_d_TB_conv

  ! ***NEW***
  !
  !=========================================================!
  subroutine pw90common_fourier_R_to_k_vec(rs_region, wann_data, ws_distance, ws_vec, OO_R, kpt, &
                                           real_lattice, mp_grid, num_wann, seedname, stdout, &
                                           OO_true, OO_pseudo)
    !====================================================================!
    !                                                                    !
    !! For OO_true (true vector):
    !! $${\vec O}_{ij}(k) = \sum_R e^{+ik.R} {\vec O}_{ij}(R)$$
    !                                                                    !
    !====================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_param_types, only: ws_region_type, wannier_data_type
    use w90_ws_distance, only: ws_translate_dist, ws_distance_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: rs_region
    type(wannier_data_type), intent(in) :: wann_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: ws_vec

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in)  :: OO_R(:, :, :, :)
    complex(kind=dp), optional, intent(out) :: OO_true(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_pseudo(:, :, :)

    character(len=50), intent(in) :: seedname

    ! local variables
    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (rs_region%use_ws_distance) CALL ws_translate_dist(ws_distance, stdout, seedname, &
                                                          rs_region, num_wann, &
                                                          wann_data%centres, real_lattice, &
                                                          mp_grid, ws_vec%nrpts, ws_vec%irvec)

    if (present(OO_true)) OO_true = cmplx_0
    if (present(OO_pseudo)) OO_pseudo = cmplx_0
    do ir = 1, ws_vec%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (rs_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir) &
                                                               *ws_distance%ndeg(i, j, ir), dp)
            if (present(OO_true)) then
              OO_true(i, j, 1) = OO_true(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_true(i, j, 2) = OO_true(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_true(i, j, 3) = OO_true(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
            if (present(OO_pseudo)) then
              OO_pseudo(i, j, 1) = OO_pseudo(i, j, 1) &
                                   + cmplx_i*ws_distance%crdist(2, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 3) &
                                   - cmplx_i*ws_distance%crdist(3, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 2)
              OO_pseudo(i, j, 2) = OO_pseudo(i, j, 2) &
                                   + cmplx_i*ws_distance%crdist(3, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 1) &
                                   - cmplx_i*ws_distance%crdist(1, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 3)
              OO_pseudo(i, j, 3) = OO_pseudo(i, j, 3) &
                                   + cmplx_i*ws_distance%crdist(1, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 2) &
                                   - cmplx_i*ws_distance%crdist(2, ideg, i, j, ir) &
                                   *phase_fac*OO_R(i, j, ir, 1)
            endif
          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), ws_vec%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir), dp)
        if (present(OO_true)) then
          OO_true(:, :, 1) = OO_true(:, :, 1) + phase_fac*OO_R(:, :, ir, 1)
          OO_true(:, :, 2) = OO_true(:, :, 2) + phase_fac*OO_R(:, :, ir, 2)
          OO_true(:, :, 3) = OO_true(:, :, 3) + phase_fac*OO_R(:, :, ir, 3)
        endif
        if (present(OO_pseudo)) then
          OO_pseudo(:, :, 1) = OO_pseudo(:, :, 1) &
                               + cmplx_i*ws_vec%crvec(2, ir)*phase_fac*OO_R(:, :, ir, 3) &
                               - cmplx_i*ws_vec%crvec(3, ir)*phase_fac*OO_R(:, :, ir, 2)
          OO_pseudo(:, :, 2) = OO_pseudo(:, :, 2) &
                               + cmplx_i*ws_vec%crvec(3, ir)*phase_fac*OO_R(:, :, ir, 1) &
                               - cmplx_i*ws_vec%crvec(1, ir)*phase_fac*OO_R(:, :, ir, 3)
          OO_pseudo(:, :, 3) = OO_pseudo(:, :, 3) &
                               + cmplx_i*ws_vec%crvec(1, ir)*phase_fac*OO_R(:, :, ir, 2) &
                               - cmplx_i*ws_vec%crvec(2, ir)*phase_fac*OO_R(:, :, ir, 1)
        endif
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec

  !=========================================================!
  subroutine pw90common_fourier_R_to_k_vec_dadb(rs_region, wann_data, ws_distance, ws_vec, OO_R, &
                                                kpt, real_lattice, mp_grid, num_wann, seedname, &
                                                stdout, OO_da, OO_dadb)
    !====================================================================!
    !                                                                    !
    !! For $$OO_{ij;dx,dy,dz}$$:
    !! $$O_{ij;dx,dy,dz}(k) = \sum_R e^{+ik.R} O_{ij;dx,dy,dz}(R)$$
    !! For $$OO_{ij;dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$O_{ij;dx1,dy1,dz1;dx2,dy2,dz2}(k) = \sum_R e^{+ik.R} i.R_{dx2,dy2,dz2}
    !!                                       .O_{ij;dx1,dy1,dz1}(R)$$
    !                                                                    !
    !====================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_param_types, only: ws_region_type, wannier_data_type
    use w90_ws_distance, only: ws_translate_dist, ws_distance_type

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: rs_region
    type(wannier_data_type), intent(in) :: wann_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: ws_vec

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in)  :: OO_R(:, :, :, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out) :: OO_dadb(:, :, :, :)

    character(len=50), intent(in) :: seedname

    ! local variables
    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (rs_region%use_ws_distance) CALL ws_translate_dist(ws_distance, stdout, seedname, &
                                                          rs_region, num_wann, &
                                                          wann_data%centres, real_lattice, &
                                                          mp_grid, ws_vec%nrpts, ws_vec%irvec)

    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, ws_vec%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (rs_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(ws_vec%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO_da)) then
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) &
                                        + cmplx_i*ws_distance%crdist(b, ideg, i, j, ir) &
                                        *phase_fac*OO_R(i, j, ir, a)
                enddo
              enddo
            endif

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), ws_vec%irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir), dp)
        if (present(OO_da)) then
          OO_da(:, :, 1) = OO_da(:, :, 1) + phase_fac*OO_R(:, :, ir, 1)
          OO_da(:, :, 2) = OO_da(:, :, 2) + phase_fac*OO_R(:, :, ir, 2)
          OO_da(:, :, 3) = OO_da(:, :, 3) + phase_fac*OO_R(:, :, ir, 3)
        endif
        if (present(OO_dadb)) then
          do a = 1, 3
            do b = 1, 3
              OO_dadb(:, :, a, b) = OO_dadb(:, :, a, b) &
                                    + cmplx_i*ws_vec%crvec(b, ir)*phase_fac*OO_R(:, :, ir, a)
            enddo
          enddo
        endif
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec_dadb

  !=========================================================!
  subroutine pw90common_fourier_R_to_k_vec_dadb_TB_conv(rs_region, wann_data, ws_distance, &
                                                        ws_vec, OO_R, kpt, real_lattice, &
                                                        mp_grid, num_wann, seedname, stdout, &
                                                        OO_da, OO_dadb)
    !====================================================================!
    !                                                                    !
    ! modified version of pw90common_fourier_R_to_k_vec_dadb, includes wannier centres in
    ! the exponential inside the sum (so called TB convention)
    !
    !! For $$OO_{ij;dx,dy,dz}$$:
    !! $$O_{ij;dx,dy,dz}(k) = \sum_R e^{+ik.(R+tau_ij)} O_{ij;dx,dy,dz}(R)$$
    !! For $$OO_{ij;dx1,dy1,dz1;dx2,dy2,dz2}$$:
    !! $$O_{ij;dx1,dy1,dz1;dx2,dy2,dz2}(k) = \sum_R e^{+ik.(R+tau_ij)} i.(R+tau_ij)_{dx2,dy2,dz2}
    !!                                       .O_{ij;dx1,dy1,dz1}(R)$$
    ! with tau_ij = tau_j - tau_i, being tau_i=<0i|r|0i> the individual wannier centres
    !                                                                    !
    !====================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_param_types, only: ws_region_type, wannier_data_type
    use w90_ws_distance, only: ws_translate_dist, ws_distance_type
    use w90_utility, only: utility_cart_to_frac, utility_inverse_mat

    implicit none

    ! arguments
    type(ws_region_type), intent(in) :: rs_region
    type(wannier_data_type), intent(in) :: wann_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: ws_vec

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)

    complex(kind=dp), intent(in) :: OO_R(:, :, :, :)
    complex(kind=dp), optional, intent(out) :: OO_da(:, :, :)
    complex(kind=dp), optional, intent(out)   :: OO_dadb(:, :, :, :)

    character(len=50), intent(in) :: seedname

    ! local variables
    real(kind=dp)    :: inv_lattice(3, 3)
    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac
    real(kind=dp)    :: local_wannier_centres(3, num_wann), wannier_centres_frac(3, num_wann)
    real(kind=dp)                                                 :: r_sum(3)

    r_sum = 0.d0

    if (rs_region%use_ws_distance) CALL ws_translate_dist(ws_distance, stdout, seedname, &
                                                          rs_region, num_wann, &
                                                          wann_data%centres, real_lattice, &
                                                          mp_grid, ws_vec%nrpts, ws_vec%irvec)

    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0

    ! calculate wannier centres in cartesian
    local_wannier_centres(:, :) = 0.d0
    do j = 1, num_wann
      do ir = 1, ws_vec%nrpts
        if ((ws_vec%irvec(1, ir) .eq. 0) .and. (ws_vec%irvec(2, ir) .eq. 0) &
            .and. (ws_vec%irvec(3, ir) .eq. 0)) then
          local_wannier_centres(1, j) = real(OO_R(j, j, ir, 1))
          local_wannier_centres(2, j) = real(OO_R(j, j, ir, 2))
          local_wannier_centres(3, j) = real(OO_R(j, j, ir, 3))
        endif
      enddo
    enddo
    ! rotate wannier centres from cartesian to fractional coordinates
    wannier_centres_frac(:, :) = 0.d0
    call utility_inverse_mat(real_lattice, inv_lattice)
    do ir = 1, num_wann
      call utility_cart_to_frac(local_wannier_centres(:, ir), wannier_centres_frac(:, ir), &
                                inv_lattice)
    enddo

!    print *, 'wannier_centres_frac'
!    do ir = 1,num_wann
!     print *, wannier_centres_frac(:,ir)
!    enddo
!    stop
!
!    print *, 'crvec'
!    do ir = 1,nrpts
!     print *, crvec(:,ir)
!    enddo
!    stop
!    print *, 'wannier_centres'
!    do ir = 1,num_wann
!     print *, wannier_centres(:,ir)
!    enddo
!    stop

    do ir = 1, ws_vec%nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (rs_region%use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(ws_distance%irdist(:, ideg, i, j, ir) + &
                                                   wannier_centres_frac(:, j) &
                                                   - wannier_centres_frac(:, i), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                        /real(ws_vec%ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            if (present(OO_da)) then
              ! if we are at the origin and at the same band, then the
              ! matrix element is zero in this convention
              if ((ws_vec%irvec(1, ir) .eq. 0) .and. (ws_vec%irvec(2, ir) .eq. 0) .and. &
                  (ws_vec%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
                cycle
              else
                OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
                OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
                OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
              endif
            endif
            if (present(OO_dadb)) then
              ! same skip as before
              if ((ws_vec%irvec(1, ir) .eq. 0) .and. (ws_vec%irvec(2, ir) .eq. 0) .and. &
                  (ws_vec%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
                cycle
              else
                do a = 1, 3
                  do b = 1, 3
                    OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + cmplx_i* &
                                          (ws_distance%crdist(b, ideg, i, j, ir) &
                                           + local_wannier_centres(b, j) &
                                           - local_wannier_centres(b, i)) &
                                          *phase_fac*OO_R(i, j, ir, a)
                  enddo
                enddo
              endif
            endif

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        do j = 1, num_wann
        do i = 1, num_wann
          r_sum(:) = real(ws_vec%irvec(:, ir)) &
                     + wannier_centres_frac(:, j) - wannier_centres_frac(:, i)
          rdotk = twopi*dot_product(kpt(:), r_sum(:))
          phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_vec%ndegen(ir), dp)
          if (present(OO_da)) then
            ! if we are at the origin and at the same band, then the
            ! matrix element is zero in this convention
            if ((ws_vec%irvec(1, ir) .eq. 0) .and. (ws_vec%irvec(2, ir) .eq. 0) .and. &
                (ws_vec%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*(OO_R(i, j, ir, 1) &
                                                           - local_wannier_centres(1, j))
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*(OO_R(i, j, ir, 2) &
                                                           - local_wannier_centres(2, j))
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*(OO_R(i, j, ir, 3) &
                                                           - local_wannier_centres(3, j))
!            print *, 'OO_R(i,j,ir,1)', OO_R(i,j,ir,1)
!            print *, 'local_wannier_centres(1,j)', local_wannier_centres(1,j)
!            print *, 'OO_R(i,j,ir,2)', OO_R(i,j,ir,2)
!            print *, 'local_wannier_centres(2,j)', local_wannier_centres(2,j)
              cycle
            else
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
          endif
          if (present(OO_dadb)) then
            ! same skip as before
            if ((ws_vec%irvec(1, ir) .eq. 0) .and. (ws_vec%irvec(2, ir) .eq. 0) .and. &
                (ws_vec%irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + cmplx_i*(ws_vec%crvec(b, ir) &
                                                                       + local_wannier_centres(b, j) &
                                                                       - local_wannier_centres(b, i))*phase_fac* &
                                        (OO_R(i, j, ir, a) - local_wannier_centres(a, j))
                enddo
              enddo
!           cycle
            else
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + cmplx_i*(ws_vec%crvec(b, ir) &
                                                                       + local_wannier_centres(b, j) &
                                                                       - local_wannier_centres(b, i))*phase_fac*OO_R(i, j, ir, a)
                enddo
              enddo
            endif
          endif
        enddo
        enddo
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec_dadb_TB_conv

  !===========================================================!
  !                   PRIVATE PROCEDURES                      !
  !===========================================================!

  !================================!
  subroutine wigner_seitz(verbose, real_lattice, mp_grid, ws_vec, stdout, seedname, &
                          count_pts, world)
    !================================!
    !! Calculates a grid of lattice vectors r that fall inside (and eventually
    !! on the surface of) the Wigner-Seitz supercell centered on the
    !! origin of the Bravais superlattice with primitive translations
    !! mp_grid(1)*a_1, mp_grid(2)*a_2, and mp_grid(3)*a_3
    !==========================================================================!

    use w90_constants, only: dp
    use w90_io, only: io_error, io_stopwatch
    use w90_param_types, only: print_output_type
    use w90_utility, only: utility_metric
    use w90_comms, only: w90commtype, mpirank

    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points

    ! arguments
    type(print_output_type), intent(in) :: verbose
    type(w90commtype), intent(in) :: world
    type(wigner_seitz_type), intent(inout) :: ws_vec
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout
    logical, intent(in) :: count_pts
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    character(len=50), intent(in) :: seedname

    ! local variables
    integer       :: ndiff(3)
    real(kind=dp) :: dist(125), tot, dist_min
    integer       :: n1, n2, n3, i1, i2, i3, icnt, i, j, ir
    real(kind=dp) :: real_metric(3, 3)
    logical :: on_root = .false.
    if (mpirank(world) == 0) on_root = .true.

    if (verbose%timing_level > 1 .and. on_root) &
      call io_stopwatch('postw90_common: wigner_seitz', 1, stdout, seedname)

    call utility_metric(real_lattice, real_metric)
    ! The Wannier functions live in a periodic supercell of the real space unit
    ! cell. This supercell is mp_grid(i) unit cells long along each primitive
    ! translation vector a_i of the unit cell
    !
    ! We loop over grid points r on a cell that is approx. 8 times
    ! larger than this "primitive supercell."
    !
    ! One of these points is in the W-S supercell if it is closer to R=0 than any
    ! of the other points R (where R are the translation vectors of the
    ! supercell). In practice it is sufficient to inspect only 125 R-points.

    ! In the end, nrpts contains the total number of grid points that have been
    ! found in the Wigner-Seitz cell

    ws_vec%nrpts = 0
    do n1 = -mp_grid(1), mp_grid(1)
      do n2 = -mp_grid(2), mp_grid(2)
        do n3 = -mp_grid(3), mp_grid(3)
          ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=0,
          ! or icnt=63
          icnt = 0
          do i1 = -2, 2
            do i2 = -2, 2
              do i3 = -2, 2
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1*mp_grid(1)
                ndiff(2) = n2 - i2*mp_grid(2)
                ndiff(3) = n3 - i3*mp_grid(3)
                dist(icnt) = 0.0_dp
                do i = 1, 3
                  do j = 1, 3
                    dist(icnt) = dist(icnt) + &
                                 real(ndiff(i), dp)*real_metric(i, j)*real(ndiff(j), dp)
                  enddo
                enddo
              enddo
            enddo
          enddo
          dist_min = minval(dist)
          if (abs(dist(63) - dist_min) .lt. 1.e-7_dp) then
            ws_vec%nrpts = ws_vec%nrpts + 1
            if (.not. count_pts) then
              ws_vec%ndegen(ws_vec%nrpts) = 0
              do i = 1, 125
                if (abs(dist(i) - dist_min) .lt. 1.e-7_dp) &
                  ws_vec%ndegen(ws_vec%nrpts) = ws_vec%ndegen(ws_vec%nrpts) + 1
              end do
              ws_vec%irvec(1, ws_vec%nrpts) = n1
              ws_vec%irvec(2, ws_vec%nrpts) = n2
              ws_vec%irvec(3, ws_vec%nrpts) = n3
              !
              ! Remember which grid point r is at the origin
              !
              if (n1 == 0 .and. n2 == 0 .and. n3 == 0) ws_vec%rpt_origin = ws_vec%nrpts
            endif
          end if

          !n3
        enddo
        !n2
      enddo
      !n1
    enddo
    !
    if (count_pts) then
      if (verbose%timing_level > 1 .and. on_root) &
        call io_stopwatch('postw90_common: wigner_seitz', 2, stdout, seedname)
      return
    end if

    if (verbose%iprint >= 3 .and. on_root) then
      write (stdout, '(1x,i4,a,/)') ws_vec%nrpts, &
        ' lattice points in Wigner-Seitz supercell:'
      do ir = 1, ws_vec%nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', ws_vec%irvec(1, ir), &
          ws_vec%irvec(2, ir), ws_vec%irvec(3, ir), '  degeneracy: ', ws_vec%ndegen(ir)
      enddo
    endif
    ! Check the "sum rule"
    tot = 0.0_dp
    do ir = 1, ws_vec%nrpts
      !
      ! Corrects weights in Fourier sums for R-vectors on the boundary of the
      ! W-S supercell
      !
      tot = tot + 1.0_dp/real(ws_vec%ndegen(ir), dp)
    enddo
    if (abs(tot - real(mp_grid(1)*mp_grid(2)*mp_grid(3), dp)) > 1.e-8_dp) &
      call io_error('ERROR in wigner_seitz: error in finding Wigner-Seitz points', stdout, seedname)

    if (verbose%timing_level > 1 .and. on_root) &
      call io_stopwatch('postw90_common: wigner_seitz', 2, stdout, seedname)

    return
  end subroutine wigner_seitz

end module w90_postw90_common
