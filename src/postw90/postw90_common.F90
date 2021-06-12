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
  use w90_comms
  use w90_constants, only: dp

  implicit none

  private

  public :: pw90common_wanint_setup, pw90common_wanint_get_kpoint_file, pw90common_wanint_param_dist
  public :: pw90common_wanint_data_dist, pw90common_get_occ
  public :: pw90common_fourier_R_to_k, pw90common_fourier_R_to_k_new, pw90common_fourier_R_to_k_vec
  public :: nrpts, rpt_origin, v_matrix, ndegen, irvec, crvec
  public :: num_int_kpts_on_node, int_kpts, weight
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
  integer, allocatable       :: irvec(:, :)
  real(kind=dp), allocatable :: crvec(:, :)
  integer, allocatable       :: ndegen(:)
  integer                    :: nrpts
  integer                    :: rpt_origin

  integer                       :: max_int_kpts_on_node, num_int_kpts
  integer, allocatable          :: num_int_kpts_on_node(:)
  real(kind=dp), allocatable    :: int_kpts(:, :), weight(:)
  complex(kind=dp), allocatable :: v_matrix(:, :, :)

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       !
  !===========================================================!

  ! Public procedures have names starting with wanint_

  subroutine pw90common_wanint_setup
    !! Setup data ready for interpolation
    use w90_constants, only: dp, cmplx_0
    use w90_io, only: io_error, io_file_unit, stdout, seedname
    use w90_utility, only: utility_cart_to_frac
    use w90_parameters, only: real_lattice, effective_model, num_wann

    integer        :: ierr, i, j, k, ikpt, ir, file_unit, num_wann_loc

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
                        //trim(seedname)//'_HH_R.dat and '//trim(seedname)//'.win')
        read (file_unit, *) nrpts
        close (file_unit)
      endif
      call comms_bcast(nrpts, 1)
    else
      call wigner_seitz(count_pts=.true.)
    endif

    ! Now can allocate several arrays
    !
    allocate (irvec(3, nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating irvec in pw90common_wanint_setup')
    irvec = 0
    allocate (crvec(3, nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating crvec in pw90common_wanint_setup')
    crvec = 0.0_dp
    allocate (ndegen(nrpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ndegen in pw90common_wanint_setup')
    ndegen = 0
    !
    ! Also rpt_origin, so that when effective_model=.true it is not
    ! passed to get_HH_R without being initialized.
    rpt_origin = 0

    ! If effective_model, this is done in get_HH_R
    if (.not. effective_model) then
      !
      ! Set up the lattice vectors on the Wigner-Seitz supercell
      ! where the Wannier functions live
      !
      call wigner_seitz(count_pts=.false.)
      !
      ! Convert from reduced to Cartesian coordinates
      !
      do ir = 1, nrpts
        ! Note that 'real_lattice' stores the lattice vectors as *rows*
        crvec(:, ir) = matmul(transpose(real_lattice), irvec(:, ir))
      end do
    endif

    return

101 call io_error('Error in pw90common_wanint_setup: problem opening file '// &
                  trim(seedname)//'_HH_R.dat')

  end subroutine pw90common_wanint_setup

  !===========================================================!
  subroutine pw90common_wanint_get_kpoint_file
    !===========================================================!
    !                                                           !
    !! read kpoints from kpoint.dat and distribute
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp
    use w90_io, only: io_error, io_file_unit, &
      io_date, io_time, io_stopwatch

    integer       :: k_unit
    integer       :: loop_nodes, loop_kpt, i, ierr
    real(kind=dp) :: sum

    k_unit = io_file_unit()
    if (on_root) then
      open (unit=k_unit, file='kpoint.dat', status='old', form='formatted', err=106)
      read (k_unit, *) num_int_kpts
    end if
    call comms_bcast(num_int_kpts, 1)

    allocate (num_int_kpts_on_node(0:num_nodes - 1))
    num_int_kpts_on_node(:) = num_int_kpts/num_nodes
    max_int_kpts_on_node = num_int_kpts - (num_nodes - 1)*(num_int_kpts/num_nodes)
    num_int_kpts_on_node(0) = max_int_kpts_on_node
!    if(my_node_id < num_int_kpts- num_int_kpts_on_node*num_nodes)  num_int_kpts_on_node= num_int_kpts_on_node+1

    allocate (int_kpts(3, max_int_kpts_on_node), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating max_int_kpts_on_node in param_read_um')
    int_kpts = 0.0_dp
    allocate (weight(max_int_kpts_on_node), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating weight in param_read_um')
    weight = 0.0_dp

    sum = 0.0_dp
    if (on_root) then
      do loop_nodes = 1, num_nodes - 1
        do loop_kpt = 1, num_int_kpts_on_node(loop_nodes)
          read (k_unit, *) (int_kpts(i, loop_kpt), i=1, 3), weight(loop_kpt)
          sum = sum + weight(loop_kpt)
        end do

        call comms_send(int_kpts(1, 1), 3*num_int_kpts_on_node(loop_nodes), loop_nodes)
        call comms_send(weight(1), num_int_kpts_on_node(loop_nodes), loop_nodes)

      end do
      do loop_kpt = 1, num_int_kpts_on_node(0)
        read (k_unit, *) (int_kpts(i, loop_kpt), i=1, 3), weight(loop_kpt)
        sum = sum + weight(loop_kpt)
      end do
!       print*,'rsum',sum
    end if

    if (.not. on_root) then
      call comms_recv(int_kpts(1, 1), 3*num_int_kpts_on_node(my_node_id), root_id)
      call comms_recv(weight(1), num_int_kpts_on_node(my_node_id), root_id)

    end if

    return

106 call io_error('Error: Problem opening file kpoint.dat in pw90common_wanint_get_kpoint_file')

  end subroutine pw90common_wanint_get_kpoint_file

  !===========================================================!
  subroutine pw90common_wanint_param_dist
    !===========================================================!
    !                                                           !
    !! distribute the parameters across processors
    !! NOTE: we only send the ones postw90 uses, not all in w90
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, io_date, io_time, &
      io_stopwatch
    use w90_parameters

    integer :: ierr

    call comms_bcast(effective_model, 1)
    call comms_bcast(eig_found, 1)

    if (.not. effective_model) then
      call comms_bcast(mp_grid(1), 3)
      call comms_bcast(num_kpts, 1)
      call comms_bcast(num_bands, 1)
    endif
    call comms_bcast(num_wann, 1)
    call comms_bcast(timing_level, 1)
    call comms_bcast(iprint, 1)
    call comms_bcast(ws_distance_tol, 1)
    call comms_bcast(ws_search_size(1), 3)
!    call comms_bcast(num_atoms,1)   ! Ivo: not used in postw90, right?
!    call comms_bcast(num_species,1) ! Ivo: not used in postw90, right?
    call comms_bcast(real_lattice(1, 1), 9)
    call comms_bcast(recip_lattice(1, 1), 9)
    call comms_bcast(real_metric(1, 1), 9)
    call comms_bcast(recip_metric(1, 1), 9)
    call comms_bcast(cell_volume, 1)
    call comms_bcast(dos_energy_step, 1)
    call comms_bcast(dos_adpt_smr, 1)
    call comms_bcast(dos_smr_index, 1)
    call comms_bcast(dos_kmesh_spacing, 1)
    call comms_bcast(dos_kmesh(1), 3)
    call comms_bcast(dos_adpt_smr_max, 1)
    call comms_bcast(dos_smr_fixed_en_width, 1)
    call comms_bcast(dos_adpt_smr_fac, 1)
    call comms_bcast(num_dos_project, 1)

    call comms_bcast(berry, 1)
    call comms_bcast(berry_task, len(berry_task))
    call comms_bcast(berry_kmesh_spacing, 1)
    call comms_bcast(berry_kmesh(1), 3)
    call comms_bcast(berry_curv_adpt_kmesh, 1)
    call comms_bcast(berry_curv_adpt_kmesh_thresh, 1)
    call comms_bcast(berry_curv_unit, len(berry_curv_unit))

! Tsirkin
    call comms_bcast(gyrotropic, 1)
    call comms_bcast(gyrotropic_task, len(gyrotropic_task))
    call comms_bcast(gyrotropic_kmesh_spacing, 1)
    call comms_bcast(gyrotropic_kmesh(1), 3)
    call comms_bcast(gyrotropic_smr_fixed_en_width, 1)
    call comms_bcast(gyrotropic_smr_index, 1)
    call comms_bcast(gyrotropic_eigval_max, 1)
    call comms_bcast(gyrotropic_nfreq, 1)
    call comms_bcast(gyrotropic_degen_thresh, 1)
    call comms_bcast(gyrotropic_num_bands, 1)
    call comms_bcast(gyrotropic_box(1, 1), 9)
    call comms_bcast(gyrotropic_box_corner(1), 3)
    call comms_bcast(gyrotropic_smr_max_arg, 1)
    call comms_bcast(gyrotropic_smr_fixed_en_width, 1)
    call comms_bcast(gyrotropic_smr_index, 1)

    call comms_bcast(spinors, 1)

    call comms_bcast(shc_freq_scan, 1)
    call comms_bcast(shc_alpha, 1)
    call comms_bcast(shc_beta, 1)
    call comms_bcast(shc_gamma, 1)
    call comms_bcast(shc_bandshift, 1)
    call comms_bcast(shc_bandshift_firstband, 1)
    call comms_bcast(shc_bandshift_energyshift, 1)

    call comms_bcast(kubo_adpt_smr, 1)
    call comms_bcast(kubo_adpt_smr_fac, 1)
    call comms_bcast(kubo_adpt_smr_max, 1)
    call comms_bcast(kubo_smr_fixed_en_width, 1)
    call comms_bcast(kubo_smr_index, 1)
    call comms_bcast(kubo_eigval_max, 1)
    call comms_bcast(kubo_nfreq, 1)
    call comms_bcast(nfermi, 1)
    call comms_bcast(dos_energy_min, 1)
    call comms_bcast(dos_energy_max, 1)
    call comms_bcast(spin_kmesh_spacing, 1)
    call comms_bcast(spin_kmesh(1), 3)
    call comms_bcast(wanint_kpoint_file, 1)
    call comms_bcast(dis_win_min, 1)
    call comms_bcast(dis_win_max, 1)
    call comms_bcast(sc_eta, 1)
    call comms_bcast(sc_w_thr, 1)
    call comms_bcast(sc_phase_conv, 1)
    call comms_bcast(sc_use_eta_corr, 1)
! ----------------------------------------------
!
! New input variables in development
!
    call comms_bcast(devel_flag, len(devel_flag))
    call comms_bcast(spin_moment, 1)
    call comms_bcast(spin_axis_polar, 1)
    call comms_bcast(spin_axis_azimuth, 1)
    call comms_bcast(spin_decomp, 1)
    call comms_bcast(use_degen_pert, 1)
    call comms_bcast(degen_thr, 1)
    call comms_bcast(num_valence_bands, 1)
    call comms_bcast(dos, 1)
    call comms_bcast(dos_task, len(dos_task))
    call comms_bcast(kpath, 1)
    call comms_bcast(kpath_task, len(kpath_task))
    call comms_bcast(kpath_bands_colour, len(kpath_bands_colour))
    call comms_bcast(kslice, 1)
    call comms_bcast(kslice_task, len(kslice_task))
    call comms_bcast(kslice_corner(1), 3)
    call comms_bcast(kslice_b1(1), 3)
    call comms_bcast(kslice_b2(1), 3)
    call comms_bcast(kslice_2dkmesh(1), 2)
    call comms_bcast(kslice_fermi_lines_colour, len(kslice_fermi_lines_colour))
    call comms_bcast(transl_inv, 1)
    call comms_bcast(num_elec_per_state, 1)
    call comms_bcast(scissors_shift, 1)
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
    call comms_bcast(geninterp, 1)
    call comms_bcast(geninterp_alsofirstder, 1)
    call comms_bcast(geninterp_single_file, 1)
    ! [gp-begin, Apr 12, 2012]
    ! BoltzWann variables
    call comms_bcast(boltzwann, 1)
    call comms_bcast(boltz_calc_also_dos, 1)
    call comms_bcast(boltz_2d_dir_num, 1)
    call comms_bcast(boltz_dos_energy_step, 1)
    call comms_bcast(boltz_dos_energy_min, 1)
    call comms_bcast(boltz_dos_energy_max, 1)
    call comms_bcast(boltz_dos_adpt_smr, 1)
    call comms_bcast(boltz_dos_smr_fixed_en_width, 1)
    call comms_bcast(boltz_dos_adpt_smr_fac, 1)
    call comms_bcast(boltz_dos_adpt_smr_max, 1)
    call comms_bcast(boltz_mu_min, 1)
    call comms_bcast(boltz_mu_max, 1)
    call comms_bcast(boltz_mu_step, 1)
    call comms_bcast(boltz_temp_min, 1)
    call comms_bcast(boltz_temp_max, 1)
    call comms_bcast(boltz_temp_step, 1)
    call comms_bcast(boltz_kmesh_spacing, 1)
    call comms_bcast(boltz_kmesh(1), 3)
    call comms_bcast(boltz_tdf_energy_step, 1)
    call comms_bcast(boltz_relax_time, 1)
    call comms_bcast(boltz_TDF_smr_fixed_en_width, 1)
    call comms_bcast(boltz_TDF_smr_index, 1)
    call comms_bcast(boltz_dos_smr_index, 1)
    call comms_bcast(boltz_bandshift, 1)
    call comms_bcast(boltz_bandshift_firstband, 1)
    call comms_bcast(boltz_bandshift_energyshift, 1)
    ! [gp-end]
    call comms_bcast(use_ws_distance, 1)

    ! These variables are different from the ones above in that they are
    ! allocatable, and in param_read they were allocated on the root node only
    !
    if (.not. on_root) then
      allocate (fermi_energy_list(nfermi), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating fermi_energy_read in postw90_param_dist')
      allocate (kubo_freq_list(kubo_nfreq), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating kubo_freq_list in postw90_param_dist')

      allocate (gyrotropic_band_list(gyrotropic_num_bands), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating gyrotropic_band_list in postw90_param_dist')

      allocate (gyrotropic_freq_list(gyrotropic_nfreq), stat=ierr)
      if (ierr /= 0) call io_error( &
        'Error allocating gyrotropic_freq_list in postw90_param_dist')

      allocate (dos_project(num_dos_project), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating dos_project in postw90_param_dist')
      if (.not. effective_model) then
        if (eig_found) then
          allocate (eigval(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating eigval in postw90_param_dist')
        end if
        allocate (kpt_latt(3, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error allocating kpt_latt in postw90_param_dist')
      endif
    end if
    if (nfermi > 0) call comms_bcast(fermi_energy_list(1), nfermi)
    call comms_bcast(gyrotropic_freq_list(1), gyrotropic_nfreq)
    call comms_bcast(gyrotropic_band_list(1), gyrotropic_num_bands)
    call comms_bcast(kubo_freq_list(1), kubo_nfreq)
    call comms_bcast(dos_project(1), num_dos_project)
    if (.not. effective_model) then
      if (eig_found) then
        call comms_bcast(eigval(1, 1), num_bands*num_kpts)
      end if
      call comms_bcast(kpt_latt(1, 1), 3*num_kpts)
    endif

    ! kmesh: only nntot,wb, and bk are needed to evaluate the WF matrix
    ! elements of the position operator in reciprocal space. For the
    ! extra matrix elements entering the orbital magnetization, also
    ! need nnlist. In principle could only broadcast those four variables

    if (.not. effective_model) then

      call comms_bcast(nnh, 1)
      call comms_bcast(nntot, 1)

      if (.not. on_root) then
        allocate (nnlist(num_kpts, nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nnlist in pw90common_wanint_param_dist')
        allocate (neigh(num_kpts, nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating neigh in pw90common_wanint_param_dist')
        allocate (nncell(3, num_kpts, nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating nncell in pw90common_wanint_param_dist')
        allocate (wb(nntot), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating wb in pw90common_wanint_param_dist')
        allocate (bka(3, nntot/2), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bka in pw90common_wanint_param_dist')
        allocate (bk(3, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) &
          call io_error('Error in allocating bk in pw90common_wanint_param_dist')
      end if

      call comms_bcast(nnlist(1, 1), num_kpts*nntot)
      call comms_bcast(neigh(1, 1), num_kpts*nntot/2)
      call comms_bcast(nncell(1, 1, 1), 3*num_kpts*nntot)
      call comms_bcast(wb(1), nntot)
      call comms_bcast(bka(1, 1), 3*nntot/2)
      call comms_bcast(bk(1, 1, 1), 3*nntot*num_kpts)

    endif

  end subroutine pw90common_wanint_param_dist

  !===========================================================!
  subroutine pw90common_wanint_data_dist
    !===========================================================!
    !                                                           !
    !! Distribute the um and chk files
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, io_file_unit, &
      io_date, io_time, io_stopwatch
    use w90_parameters, only: num_wann, num_kpts, num_bands, have_disentangled, &
      u_matrix_opt, u_matrix, m_matrix, &
      ndimwin, lwindow, nntot, wannier_centres, &
      num_valence_bands, scissors_shift

    implicit none

    integer :: ierr, loop_kpt, m, i, j

    if (.not. on_root) then
      ! wannier_centres is allocated in param_read, so only on root node
      ! It is then read in param_read_chpkt
      ! Therefore, now we need to allocate it on all nodes, and then broadcast it
      allocate (wannier_centres(3, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating wannier_centres in pw90common_wanint_data_dist')
    end if
    call comms_bcast(wannier_centres(1, 1), 3*num_wann)

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
      call io_error('Error allocating v_matrix in pw90common_wanint_data_dist')
    ! u_matrix and u_matrix_opt are stored on root only
    if (on_root) then
      if (.not. have_disentangled) then
        v_matrix = u_matrix
      else
        v_matrix = cmplx_0
        do loop_kpt = 1, num_kpts
          do j = 1, num_wann
            do m = 1, ndimwin(loop_kpt)
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
    call comms_bcast(v_matrix(1, 1, 1), num_bands*num_wann*num_kpts)

    if (num_valence_bands > 0 .and. abs(scissors_shift) > 1.0e-7_dp) then
    if (.not. on_root .and. .not. allocated(u_matrix)) then
      allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) &
        call io_error('Error allocating u_matrix in pw90common_wanint_data_dist')
    endif
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts)
    endif

!    if (.not.on_root .and. .not.allocated(m_matrix)) then
!       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!       if (ierr/=0)&
!            call io_error('Error allocating m_matrix in pw90common_wanint_data_dist')
!    endif
!    call comms_bcast(m_matrix(1,1,1,1),num_wann*num_wann*nntot*num_kpts)

    call comms_bcast(have_disentangled, 1)

    if (have_disentangled) then
      if (.not. on_root) then

        ! Do we really need these 'if not allocated'? Didn't use them for
        ! eigval and kpt_latt above!

!          if (.not.allocated(u_matrix_opt)) then
!             allocate(u_matrix_opt(num_bands,num_wann,num_kpts),stat=ierr)
!             if (ierr/=0)&
!              call io_error('Error allocating u_matrix_opt in pw90common_wanint_data_dist')
!          endif

        if (.not. allocated(lwindow)) then
          allocate (lwindow(num_bands, num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating lwindow in pw90common_wanint_data_dist')
        endif

        if (.not. allocated(ndimwin)) then
          allocate (ndimwin(num_kpts), stat=ierr)
          if (ierr /= 0) &
            call io_error('Error allocating ndimwin in pw90common_wanint_data_dist')
        endif

      end if

!       call comms_bcast(u_matrix_opt(1,1,1),num_bands*num_wann*num_kpts)
      call comms_bcast(lwindow(1, 1), num_bands*num_kpts)
      call comms_bcast(ndimwin(1), num_kpts)
    end if

  end subroutine pw90common_wanint_data_dist

!=======================================================================

  subroutine pw90common_get_occ(eig, occ, ef)
    !! Compute the electronic occupancy

    use w90_constants, only: dp !,eps7
    use w90_parameters, only: num_wann !,smear_temp
!    use w90_constants, only    : elem_charge_SI,k_B_SI

    ! Arguments
    !
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

  function kmesh_spacing_singleinteger(num_points)

    !! Set up the value of the interpolation mesh spacing, needed for
    !! adaptive smearing [see Eqs. (34-35) YWVS07]. Choose it as the largest of
    !! the three Delta_k's for each of the primitive translations b1, b2, and b3

    use w90_parameters, only: recip_lattice

    integer, intent(in) :: num_points
    real(kind=dp)       :: kmesh_spacing_singleinteger

    integer        :: i
    real(kind=dp) :: Delta_k_i(3)

    ! NOTE: The vectors b_i are stored as *rows* in recip_lattice (argh!).
    ! Hence I believe Jonathan's original code confused rows with columns
    ! when computing Delta_k, which he called 'rspace'
    ! (See my e-mail of 20Sept07)
    !
    do i = 1, 3
      Delta_k_i(i) = sqrt(dot_product(recip_lattice(i, :), recip_lattice(i, :))) &
                     /num_points
    end do
    kmesh_spacing_singleinteger = maxval(Delta_k_i)

  end function kmesh_spacing_singleinteger
  function kmesh_spacing_mesh(mesh)
    !! Same as kmesh_spacing_singleinteger, but for a kmesh with three
    !! different mesh samplings along the three directions
    use w90_parameters, only: recip_lattice

    integer, dimension(3), intent(in) :: mesh
    real(kind=dp)                     :: kmesh_spacing_mesh

    integer        :: i
    real(kind=dp) :: Delta_k_i(3)

    do i = 1, 3
      Delta_k_i(i) = sqrt(dot_product(recip_lattice(i, :), recip_lattice(i, :))) &
                     /mesh(i)
    end do
    kmesh_spacing_mesh = maxval(Delta_k_i)

  end function kmesh_spacing_mesh
  !
  !=========================================================!
  subroutine pw90common_fourier_R_to_k(kpt, OO_R, OO, alpha)
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
    use w90_parameters, only: num_kpts, kpt_latt, num_wann, use_ws_distance
    use w90_ws_distance, only: irdist_ws, crdist_ws, &
      wdist_ndeg, ws_translate_dist

    implicit none

    ! Arguments
    !
    real(kind=dp)                                   :: kpt(3)
    complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R
    complex(kind=dp), dimension(:, :), intent(out)   :: OO
    integer                                         :: alpha

    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (use_ws_distance) CALL ws_translate_dist(nrpts, irvec)

    OO(:, :) = cmplx_0
    do ir = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, wdist_ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(irdist_ws(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir)*wdist_ndeg(i, j, ir), dp)
            if (alpha == 0) then
              OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
              OO(i, j) = OO(i, j) + cmplx_i*crdist_ws(alpha, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir)
            else
              stop 'wrong value of alpha in pw90common_fourier_R_to_k'
            endif
          enddo
        enddo
        enddo
      else
        ! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
        if (alpha == 0) then
          OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
          OO(:, :) = OO(:, :) + &
                     cmplx_i*crvec(alpha, ir)*phase_fac*OO_R(:, :, ir)
        else
          stop 'wrong value of alpha in pw90common_fourier_R_to_k'
        endif
      endif

    enddo

  end subroutine pw90common_fourier_R_to_k

  ! ***NEW***
  !
  !=========================================================!
  subroutine pw90common_fourier_R_to_k_new(kpt, OO_R, OO, OO_dx, OO_dy, OO_dz)
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
    use w90_parameters, only: timing_level, num_kpts, kpt_latt, num_wann, use_ws_distance
    use w90_ws_distance, only: irdist_ws, crdist_ws, wdist_ndeg, ws_translate_dist

    implicit none

    ! Arguments
    !
    real(kind=dp)                                             :: kpt(3)
    complex(kind=dp), dimension(:, :, :), intent(in)            :: OO_R
    complex(kind=dp), optional, dimension(:, :), intent(out)   :: OO
    complex(kind=dp), optional, dimension(:, :), intent(out)   :: OO_dx
    complex(kind=dp), optional, dimension(:, :), intent(out)   :: OO_dy
    complex(kind=dp), optional, dimension(:, :), intent(out)   :: OO_dz

    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (use_ws_distance) CALL ws_translate_dist(nrpts, irvec)

    if (present(OO)) OO = cmplx_0
    if (present(OO_dx)) OO_dx = cmplx_0
    if (present(OO_dy)) OO_dy = cmplx_0
    if (present(OO_dz)) OO_dz = cmplx_0
    do ir = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, wdist_ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(irdist_ws(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir)*wdist_ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_dx)) OO_dx(i, j) = OO_dx(i, j) + &
                                              cmplx_i*crdist_ws(1, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
            if (present(OO_dy)) OO_dy(i, j) = OO_dy(i, j) + &
                                              cmplx_i*crdist_ws(2, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
            if (present(OO_dz)) OO_dz(i, j) = OO_dz(i, j) + &
                                              cmplx_i*crdist_ws(3, ideg, i, j, ir)* &
                                              phase_fac*OO_R(i, j, ir)
          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
        if (present(OO)) OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        if (present(OO_dx)) OO_dx(:, :) = OO_dx(:, :) + &
                                          cmplx_i*crvec(1, ir)*phase_fac*OO_R(:, :, ir)
        if (present(OO_dy)) OO_dy(:, :) = OO_dy(:, :) + &
                                          cmplx_i*crvec(2, ir)*phase_fac*OO_R(:, :, ir)
        if (present(OO_dz)) OO_dz(:, :) = OO_dz(:, :) + &
                                          cmplx_i*crvec(3, ir)*phase_fac*OO_R(:, :, ir)
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new

  !=========================================================!
  subroutine pw90common_fourier_R_to_k_new_second_d(kpt, OO_R, OO, OO_da, OO_dadb)
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
    use w90_parameters, only: timing_level, num_kpts, kpt_latt, num_wann, use_ws_distance
    use w90_ws_distance, only: irdist_ws, crdist_ws, wdist_ndeg, ws_translate_dist

    implicit none

    ! Arguments
    !
    real(kind=dp)                                                 :: kpt(3)
    complex(kind=dp), dimension(:, :, :), intent(in)                :: OO_R
    complex(kind=dp), optional, dimension(:, :), intent(out)       :: OO
    complex(kind=dp), optional, dimension(:, :, :), intent(out)     :: OO_da
    complex(kind=dp), optional, dimension(:, :, :, :), intent(out)   :: OO_dadb

    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (use_ws_distance) CALL ws_translate_dist(nrpts, irvec)

    if (present(OO)) OO = cmplx_0
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, wdist_ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(irdist_ws(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir)*wdist_ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i*crdist_ws(a, ideg, i, j, ir)* &
                                 phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) - crdist_ws(a, ideg, i, j, ir)* &
                                        crdist_ws(b, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir)
                enddo
              enddo
            end if

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
        if (present(OO)) OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
        if (present(OO_da)) then
          do a = 1, 3
            OO_da(:, :, a) = OO_da(:, :, a) + cmplx_i*crvec(a, ir)*phase_fac*OO_R(:, :, ir)
          enddo
        endif
        if (present(OO_dadb)) then
          do a = 1, 3
            do b = 1, 3
              OO_dadb(:, :, a, b) = OO_dadb(:, :, a, b) - &
                                    crvec(a, ir)*crvec(b, ir)*phase_fac*OO_R(:, :, ir)
            enddo
          enddo
        end if
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_new_second_d

  subroutine pw90common_fourier_R_to_k_new_second_d_TB_conv(kpt, OO_R, oo_a_R, OO, OO_da, OO_dadb)
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
    use w90_parameters, only: timing_level, num_kpts, kpt_latt, num_wann, &
      use_ws_distance, wannier_centres, recip_lattice
    use w90_ws_distance, only: irdist_ws, crdist_ws, wdist_ndeg, ws_translate_dist
    use w90_utility, only: utility_cart_to_frac

    implicit none

    ! Arguments
    !
    real(kind=dp)                                                 :: kpt(3)
    complex(kind=dp), dimension(:, :, :), intent(in)                :: OO_R
    complex(kind=dp), optional, dimension(:, :), intent(out)       :: OO
    complex(kind=dp), optional, dimension(:, :, :), intent(out)     :: OO_da
    complex(kind=dp), optional, dimension(:, :, :, :), intent(out)   :: OO_dadb
    complex(kind=dp), dimension(:, :, :, :), intent(in)     :: oo_a_R

    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac
    real(kind=dp)    :: local_wannier_centres(3, num_wann), wannier_centres_frac(3, num_wann)
    real(kind=dp)                                                 :: r_sum(3)

    r_sum = 0.d0

    if (use_ws_distance) CALL ws_translate_dist(nrpts, irvec)

    ! calculate wannier centres in cartesian
    local_wannier_centres(:, :) = 0.d0
    do j = 1, num_wann
      do ir = 1, nrpts
        if ((irvec(1, ir) .eq. 0) .and. (irvec(2, ir) .eq. 0) .and. (irvec(3, ir) .eq. 0)) then
          local_wannier_centres(1, j) = real(oo_a_R(j, j, ir, 1))
          local_wannier_centres(2, j) = real(oo_a_R(j, j, ir, 2))
          local_wannier_centres(3, j) = real(oo_a_R(j, j, ir, 3))
        endif
      enddo
    enddo
    ! rotate wannier centres from cartesian to fractional coordinates
    wannier_centres_frac(:, :) = 0.d0
    do ir = 1, num_wann
      call utility_cart_to_frac(local_wannier_centres(:, ir), wannier_centres_frac(:, ir), recip_lattice)
    enddo

    if (present(OO)) OO = cmplx_0
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, wdist_ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(irdist_ws(:, ideg, i, j, ir) + &
                                                   wannier_centres_frac(:, j) - wannier_centres_frac(:, i), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir)*wdist_ndeg(i, j, ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i* &
                                 (crdist_ws(a, ideg, i, j, ir) + local_wannier_centres(a, j) - &
                                  local_wannier_centres(a, i))*phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) - &
                                        (crdist_ws(a, ideg, i, j, ir) + local_wannier_centres(a, j) - &
                                         local_wannier_centres(a, i))* &
                                        (crdist_ws(b, ideg, i, j, ir) + local_wannier_centres(b, j) - &
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
            r_sum(:) = real(irvec(:, ir)) + wannier_centres_frac(:, j) - wannier_centres_frac(:, i)
            rdotk = twopi*dot_product(kpt(:), r_sum(:))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
            if (present(OO)) OO(i, j) = OO(i, j) + phase_fac*OO_R(i, j, ir)
            if (present(OO_da)) then
              do a = 1, 3
                OO_da(i, j, a) = OO_da(i, j, a) + cmplx_i* &
                                 (crvec(a, ir) + local_wannier_centres(a, j) - &
                                  local_wannier_centres(a, i))*phase_fac*OO_R(i, j, ir)
              enddo
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = &
                    OO_dadb(i, j, a, b) - &
                    (crvec(a, ir) + local_wannier_centres(a, j) - local_wannier_centres(a, i))* &
                    (crvec(b, ir) + local_wannier_centres(b, j) - local_wannier_centres(b, i))* &
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
  subroutine pw90common_fourier_R_to_k_vec(kpt, OO_R, OO_true, OO_pseudo)
    !====================================================================!
    !                                                                    !
    !! For OO_true (true vector):
    !! $${\vec O}_{ij}(k) = \sum_R e^{+ik.R} {\vec O}_{ij}(R)$$
    !                                                                    !
    !====================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_parameters, only: num_kpts, kpt_latt, num_wann, use_ws_distance
    use w90_ws_distance, only: irdist_ws, crdist_ws, wdist_ndeg, ws_translate_dist

    implicit none

    ! Arguments
    !
    real(kind=dp)                                     :: kpt(3)
    complex(kind=dp), dimension(:, :, :, :), intent(in)  :: OO_R
    complex(kind=dp), optional, dimension(:, :, :), intent(out)   :: OO_true
    complex(kind=dp), optional, dimension(:, :, :), intent(out)   :: OO_pseudo

    integer          :: ir, i, j, ideg
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (use_ws_distance) CALL ws_translate_dist(nrpts, irvec)
    if (present(OO_true)) OO_true = cmplx_0
    if (present(OO_pseudo)) OO_pseudo = cmplx_0
    do ir = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, wdist_ndeg(i, j, ir)
            rdotk = twopi*dot_product(kpt(:), real(irdist_ws(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir)*wdist_ndeg(i, j, ir), dp)
            if (present(OO_true)) then
              OO_true(i, j, 1) = OO_true(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_true(i, j, 2) = OO_true(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_true(i, j, 3) = OO_true(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
            if (present(OO_pseudo)) then
              OO_pseudo(i, j, 1) = OO_pseudo(i, j, 1) &
                                   + cmplx_i*crdist_ws(2, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir, 3) &
                                   - cmplx_i*crdist_ws(3, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir, 2)
              OO_pseudo(i, j, 2) = OO_pseudo(i, j, 2) &
                                   + cmplx_i*crdist_ws(3, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir, 1) &
                                   - cmplx_i*crdist_ws(1, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir, 3)
              OO_pseudo(i, j, 3) = OO_pseudo(i, j, 3) &
                                   + cmplx_i*crdist_ws(1, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir, 2) &
                                   - cmplx_i*crdist_ws(2, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir, 1)
            endif
          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
        if (present(OO_true)) then
          OO_true(:, :, 1) = OO_true(:, :, 1) + phase_fac*OO_R(:, :, ir, 1)
          OO_true(:, :, 2) = OO_true(:, :, 2) + phase_fac*OO_R(:, :, ir, 2)
          OO_true(:, :, 3) = OO_true(:, :, 3) + phase_fac*OO_R(:, :, ir, 3)
        endif
        if (present(OO_pseudo)) then
          OO_pseudo(:, :, 1) = OO_pseudo(:, :, 1) &
                               + cmplx_i*crvec(2, ir)*phase_fac*OO_R(:, :, ir, 3) &
                               - cmplx_i*crvec(3, ir)*phase_fac*OO_R(:, :, ir, 2)
          OO_pseudo(:, :, 2) = OO_pseudo(:, :, 2) &
                               + cmplx_i*crvec(3, ir)*phase_fac*OO_R(:, :, ir, 1) &
                               - cmplx_i*crvec(1, ir)*phase_fac*OO_R(:, :, ir, 3)
          OO_pseudo(:, :, 3) = OO_pseudo(:, :, 3) &
                               + cmplx_i*crvec(1, ir)*phase_fac*OO_R(:, :, ir, 2) &
                               - cmplx_i*crvec(2, ir)*phase_fac*OO_R(:, :, ir, 1)
        endif
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec

  !=========================================================!
  subroutine pw90common_fourier_R_to_k_vec_dadb(kpt, OO_R, OO_da, OO_dadb)
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
    use w90_parameters, only: num_kpts, kpt_latt, num_wann, use_ws_distance
    use w90_ws_distance, only: irdist_ws, crdist_ws, wdist_ndeg, ws_translate_dist

    implicit none

    ! Arguments
    !
    real(kind=dp)                                     :: kpt(3)
    complex(kind=dp), dimension(:, :, :, :), intent(in)  :: OO_R
    complex(kind=dp), optional, dimension(:, :, :), intent(out)     :: OO_da
    complex(kind=dp), optional, dimension(:, :, :, :), intent(out)   :: OO_dadb

    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac

    if (use_ws_distance) CALL ws_translate_dist(nrpts, irvec)
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0
    do ir = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, wdist_ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(irdist_ws(:, ideg, i, j, ir), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir)*wdist_ndeg(i, j, ir), dp)
            if (present(OO_da)) then
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
            endif
            if (present(OO_dadb)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + &
                                        cmplx_i*crdist_ws(b, ideg, i, j, ir)*phase_fac*OO_R(i, j, ir, a)
                enddo
              enddo
            endif

          enddo
        enddo
        enddo
      else
! [lp] Original code, without IJ-dependent shift:
        rdotk = twopi*dot_product(kpt(:), irvec(:, ir))
        phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
        if (present(OO_da)) then
          OO_da(:, :, 1) = OO_da(:, :, 1) + phase_fac*OO_R(:, :, ir, 1)
          OO_da(:, :, 2) = OO_da(:, :, 2) + phase_fac*OO_R(:, :, ir, 2)
          OO_da(:, :, 3) = OO_da(:, :, 3) + phase_fac*OO_R(:, :, ir, 3)
        endif
        if (present(OO_dadb)) then
          do a = 1, 3
            do b = 1, 3
              OO_dadb(:, :, a, b) = OO_dadb(:, :, a, b) + cmplx_i*crvec(b, ir)*phase_fac*OO_R(:, :, ir, a)
            enddo
          enddo
        endif
      endif
    enddo

  end subroutine pw90common_fourier_R_to_k_vec_dadb

  !=========================================================!
  subroutine pw90common_fourier_R_to_k_vec_dadb_TB_conv(kpt, OO_R, OO_da, OO_dadb)
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
    use w90_parameters, only: num_kpts, kpt_latt, num_wann, use_ws_distance, &
      wannier_centres, recip_lattice
    use w90_ws_distance, only: irdist_ws, crdist_ws, wdist_ndeg, ws_translate_dist
    use w90_utility, only: utility_cart_to_frac

    implicit none

    ! Arguments
    !
    real(kind=dp)                                     :: kpt(3)
    complex(kind=dp), dimension(:, :, :, :), intent(in)  :: OO_R
    complex(kind=dp), optional, dimension(:, :, :), intent(out)     :: OO_da
    complex(kind=dp), optional, dimension(:, :, :, :), intent(out)   :: OO_dadb

    integer          :: ir, i, j, ideg, a, b
    real(kind=dp)    :: rdotk
    complex(kind=dp) :: phase_fac
    real(kind=dp)    :: local_wannier_centres(3, num_wann), wannier_centres_frac(3, num_wann)
    real(kind=dp)                                                 :: r_sum(3)

    r_sum = 0.d0

    if (use_ws_distance) CALL ws_translate_dist(nrpts, irvec)
    if (present(OO_da)) OO_da = cmplx_0
    if (present(OO_dadb)) OO_dadb = cmplx_0

    ! calculate wannier centres in cartesian
    local_wannier_centres(:, :) = 0.d0
    do j = 1, num_wann
      do ir = 1, nrpts
        if ((irvec(1, ir) .eq. 0) .and. (irvec(2, ir) .eq. 0) .and. (irvec(3, ir) .eq. 0)) then
          local_wannier_centres(1, j) = real(OO_R(j, j, ir, 1))
          local_wannier_centres(2, j) = real(OO_R(j, j, ir, 2))
          local_wannier_centres(3, j) = real(OO_R(j, j, ir, 3))
        endif
      enddo
    enddo
    ! rotate wannier centres from cartesian to fractional coordinates
    wannier_centres_frac(:, :) = 0.d0
    do ir = 1, num_wann
      call utility_cart_to_frac(local_wannier_centres(:, ir), wannier_centres_frac(:, ir), recip_lattice)
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

    do ir = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
      if (use_ws_distance) then
        do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, wdist_ndeg(i, j, ir)

            rdotk = twopi*dot_product(kpt(:), real(irdist_ws(:, ideg, i, j, ir) + &
                                                   wannier_centres_frac(:, j) - wannier_centres_frac(:, i), dp))
            phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir)*wdist_ndeg(i, j, ir), dp)
            if (present(OO_da)) then
              ! if we are at the origin and at the same band, then the
              ! matrix element is zero in this convention
              if ((irvec(1, ir) .eq. 0) .and. (irvec(2, ir) .eq. 0) .and. (irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
                cycle
              else
                OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*OO_R(i, j, ir, 1)
                OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*OO_R(i, j, ir, 2)
                OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*OO_R(i, j, ir, 3)
              endif
            endif
            if (present(OO_dadb)) then
              ! same skip as before
              if ((irvec(1, ir) .eq. 0) .and. (irvec(2, ir) .eq. 0) .and. (irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
                cycle
              else
                do a = 1, 3
                  do b = 1, 3
                    OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + cmplx_i* &
                                          (crdist_ws(b, ideg, i, j, ir) + local_wannier_centres(b, j) - &
                                           local_wannier_centres(b, i))*phase_fac*OO_R(i, j, ir, a)
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
          r_sum(:) = real(irvec(:, ir)) + wannier_centres_frac(:, j) - wannier_centres_frac(:, i)
          rdotk = twopi*dot_product(kpt(:), r_sum(:))
          phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
          if (present(OO_da)) then
            ! if we are at the origin and at the same band, then the
            ! matrix element is zero in this convention
            if ((irvec(1, ir) .eq. 0) .and. (irvec(2, ir) .eq. 0) .and. (irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
              OO_da(i, j, 1) = OO_da(i, j, 1) + phase_fac*(OO_R(i, j, ir, 1) - local_wannier_centres(1, j))
              OO_da(i, j, 2) = OO_da(i, j, 2) + phase_fac*(OO_R(i, j, ir, 2) - local_wannier_centres(2, j))
              OO_da(i, j, 3) = OO_da(i, j, 3) + phase_fac*(OO_R(i, j, ir, 3) - local_wannier_centres(3, j))
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
            if ((irvec(1, ir) .eq. 0) .and. (irvec(2, ir) .eq. 0) .and. (irvec(3, ir) .eq. 0) .and. (i .eq. j)) then
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + cmplx_i*(crvec(b, ir) + local_wannier_centres(b, j) - &
                                                                       local_wannier_centres(b, i))*phase_fac* &
                                        (OO_R(i, j, ir, a) - local_wannier_centres(a, j))
                enddo
              enddo
!           cycle
            else
              do a = 1, 3
                do b = 1, 3
                  OO_dadb(i, j, a, b) = OO_dadb(i, j, a, b) + cmplx_i*(crvec(b, ir) + local_wannier_centres(b, j) - &
                                                                       local_wannier_centres(b, i))*phase_fac*OO_R(i, j, ir, a)
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
  subroutine wigner_seitz(count_pts)
    !================================!
    !! Calculates a grid of lattice vectors r that fall inside (and eventually
    !! on the surface of) the Wigner-Seitz supercell centered on the
    !! origin of the Bravais superlattice with primitive translations
    !! mp_grid(1)*a_1, mp_grid(2)*a_2, and mp_grid(3)*a_3
    !==========================================================================!

    use w90_constants, only: eps8, dp
    use w90_io, only: stdout, io_error, io_stopwatch
    use w90_parameters, only: iprint, mp_grid, real_metric, timing_level, &
      ws_search_size, ws_distance_tol

    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points

    implicit none

    logical, intent(in) :: count_pts

    integer       :: ndiff(3)
    real(kind=dp) :: tot, dist_min
    real(kind=dp), allocatable :: dist(:)
    integer       :: n1, n2, n3, i1, i2, i3, icnt, i, j, ir
    integer       :: ierr, dist_dim

    if (timing_level > 1 .and. on_root) &
      call io_stopwatch('postw90_common: wigner_seitz', 1)

    dist_dim = 1
    do i = 1, 3
      dist_dim = dist_dim*((ws_search_size(i) + 1)*2 + 1)
    end do
    allocate (dist(dist_dim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating dist in wigner_seitz')

    ! The Wannier functions live in a periodic supercell of the real space unit
    ! cell. This supercell is mp_grid(i) unit cells long along each primitive
    ! translation vector a_i of the unit cell
    !
    ! We loop over grid points r on a unit cell that is (2*ws_search_size+1)**3 times
    ! larger than this primitive supercell.
    !
    ! One of these points is in the W-S supercell if it is closer to R=0 than any
    ! of the other points R (where R are the translation vectors of the
    ! supercell).

    ! In the end, nrpts contains the total number of grid points that have been
    ! found in the Wigner-Seitz cell

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
                    dist(icnt) = dist(icnt) + &
                                 real(ndiff(i), dp)*real_metric(i, j)*real(ndiff(j), dp)
                  enddo
                enddo
              enddo
            enddo
          enddo
          dist_min = minval(dist)
          if (abs(dist((dist_dim + 1)/2) - dist_min) .lt. ws_distance_tol**2) then
            nrpts = nrpts + 1
            if (.not. count_pts) then
              ndegen(nrpts) = 0
              do i = 1, dist_dim
                if (abs(dist(i) - dist_min) .lt. ws_distance_tol**2) &
                  ndegen(nrpts) = ndegen(nrpts) + 1
              end do
              irvec(1, nrpts) = n1
              irvec(2, nrpts) = n2
              irvec(3, nrpts) = n3
              !
              ! Remember which grid point r is at the origin
              !
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
    if (ierr /= 0) call io_error('Error in deallocating dist wigner_seitz')

    if (count_pts) then
      if (timing_level > 1 .and. on_root) &
        call io_stopwatch('postw90_common: wigner_seitz', 2)
      return
    end if

    ! Check the "sum rule"
    tot = 0.0_dp
    do i = 1, nrpts
      tot = tot + 1.0_dp/real(ndegen(i), dp)
    enddo

    if (iprint >= 3 .and. on_root) then
      write (stdout, '(1x,i4,a,/)') nrpts, &
        ' lattice points in Wigner-Seitz supercell:'
      do ir = 1, nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', irvec(1, ir), &
          irvec(2, ir), irvec(3, ir), '  degeneracy: ', ndegen(ir)
      enddo
      write (stdout, '(1x,a,f12.3)') ' tot = ', tot
      write (stdout, '(1x,a,i12)') ' mp_grid product = ', mp_grid(1)*mp_grid(2)*mp_grid(3)
    endif
    if (abs(tot - real(mp_grid(1)*mp_grid(2)*mp_grid(3), dp)) > eps8) then
      call io_error('ERROR in wigner_seitz: error in finding Wigner-Seitz points')

    endif

    if (timing_level > 1 .and. on_root) &
      call io_stopwatch('postw90_common: wigner_seitz', 2)

    return
  end subroutine wigner_seitz

end module w90_postw90_common
