!
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
!  w90_plot: plot various data                               !
!                                                            !
!------------------------------------------------------------!

module w90_plot_mod

  !! This module handles various plots

  use w90_comms, only: comms_reduce, w90_comm_type, mpisize, mpirank, comms_array_split

  implicit none

  private

  public :: plot_main

contains

  !================================================!

  subroutine plot_main(atom_data, band_plot, dis_manifold, fermi_energy_list, fermi_surface_plot, &
                       ham_logical, kmesh_info, kpt_latt, output_file, wvfn_read, real_space_ham, &
                       kpoint_path, print_output, wannier_data, wannier_plot, ws_region, &
                       w90_calculation, ham_k, ham_r, m_matrix, u_matrix, u_matrix_opt, eigval, &
                       real_lattice, wannier_centres_translated, bohr, irvec, mp_grid, ndegen, &
                       shift_vec, nrpts, num_bands, num_kpts, num_wann, rpt_origin, &
                       transport_mode, have_disentangled, lsitesymmetry, w90_system, seedname, &
                       stdout, timer, dist_k, error, comm)
    !================================================!
    !
    !! Main plotting routine
    !
    !================================================!

    use w90_constants, only: eps6, dp
    use w90_hamiltonian, only: hamiltonian_get_hr, hamiltonian_write_hr, hamiltonian_setup, &
      hamiltonian_write_tb
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_info_type, wannier_data_type, atom_data_type, dis_manifold_type, &
      kpoint_path_type, print_output_type, ws_region_type, ws_distance_type, timer_list_type, &
      w90_system_type
    use w90_utility, only: utility_recip_lattice_base
    use w90_wannier90_types, only: w90_calculation_type, wvfn_read_type, output_file_type, &
      fermi_surface_plot_type, band_plot_type, wannier_plot_type, real_space_ham_type, &
      ham_logical_type
    use w90_ws_distance, only: ws_translate_dist, ws_write_vec
    use w90_error, only: w90_error_type

    implicit none

    ! arguments
    type(real_space_ham_type), intent(inout) :: real_space_ham
    type(w90_error_type), allocatable, intent(out) :: error
    type(timer_list_type), intent(inout) :: timer
    type(ham_logical_type), intent(inout) :: ham_logical

    complex(kind=dp), intent(inout), allocatable :: ham_r(:, :, :)
    complex(kind=dp), intent(inout), allocatable :: ham_k(:, :, :)

    real(kind=dp), intent(inout), allocatable :: wannier_centres_translated(:, :)

    integer, intent(inout), allocatable :: irvec(:, :)
    integer, intent(inout), allocatable :: ndegen(:)
    integer, intent(inout), allocatable :: shift_vec(:, :)
    integer, intent(inout)              :: nrpts
    integer, intent(inout)              :: rpt_origin

    type(atom_data_type), intent(in) :: atom_data
    type(band_plot_type), intent(in) :: band_plot
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(fermi_surface_plot_type), intent(in) :: fermi_surface_plot
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(kpoint_path_type), intent(in) :: kpoint_path
    type(output_file_type), intent(in) :: output_file
    type(print_output_type), intent(in) :: print_output
    type(w90_calculation_type), intent(in) :: w90_calculation
    type(w90_comm_type), intent(in) :: comm
    type(w90_system_type), intent(in) :: w90_system
    type(wannier_data_type), intent(in) :: wannier_data
    type(wannier_plot_type), intent(in) :: wannier_plot
    type(ws_region_type), intent(in) :: ws_region
    type(wvfn_read_type), intent(in) :: wvfn_read

    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)

    real(kind=dp), intent(in), allocatable :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: dist_k(:)

    character(len=20), intent(in) :: transport_mode
    character(len=50), intent(in) :: seedname

    logical, intent(in) :: have_disentangled
    logical, intent(in) :: lsitesymmetry

    ! local variables
    type(ws_distance_type) :: ws_distance
    real(kind=dp) :: recip_lattice(3, 3), volume
    integer :: nkp, bands_num_spec_points, my_node_id, i
    logical :: have_gamma
    logical :: on_root = .false.

    my_node_id = mpirank(comm)

    if (my_node_id == 0) on_root = .true.

    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)

    ! setup RS cluster for calculations that need it
    if (output_file%write_hr .or. &
        output_file%write_hr_diag .or. &
        output_file%write_r2mn .or. &
        output_file%write_rmn .or. &
        output_file%write_tb .or. &
        w90_calculation%bands_plot .or. &
        w90_calculation%fermi_surface_plot) then

      call hamiltonian_setup(ham_logical, print_output, ws_region, w90_calculation, ham_k, ham_r, &
                             real_lattice, wannier_centres_translated, irvec, mp_grid, ndegen, &
                             num_kpts, num_wann, nrpts, rpt_origin, band_plot%mode, stdout, &
                             timer, error, transport_mode, comm)
      if (allocated(error)) return
    endif

    ! setup RS Hamilton eqn for calculations that need it
    if (output_file%write_hr .or. &
        output_file%write_hr_diag .or. &
        output_file%write_tb .or. &
        w90_calculation%bands_plot .or. &
        w90_calculation%fermi_surface_plot) then

      ! Check if the kmesh includes the gamma point
      have_gamma = .false.
      do nkp = 1, num_kpts
        if (all(abs(kpt_latt(:, nkp)) < eps6)) have_gamma = .true.
      end do
      if (.not. have_gamma) then
        write (stdout, '(1x,a)') '!!!! Kpoint grid does not include Gamma. '// &
          ' Interpolation may be incorrect. !!!!'
      endif

      ! Transform Hamiltonian to WF basis
      call hamiltonian_get_hr(atom_data, dis_manifold, ham_logical, real_space_ham, print_output, &
                              ham_k, ham_r, u_matrix, u_matrix_opt, eigval, kpt_latt, &
                              real_lattice, wannier_data%centres, wannier_centres_translated, &
                              irvec, shift_vec, nrpts, num_bands, num_kpts, num_wann, &
                              have_disentangled, stdout, timer, error, lsitesymmetry, comm)

      if (allocated(error)) return
    endif

    if (on_root) then
      if (print_output%timing_level > 0) call io_stopwatch_start('plot: main', timer)

      if (print_output%iprint > 0) then
        ! Print the header only if there is something to plot
        if (output_file%write_hr .or. &
            output_file%write_r2mn .or. &
            output_file%write_rmn .or. &
            output_file%write_tb .or. &
            output_file%write_u_matrices .or. &
            w90_calculation%bands_plot .or. &
            w90_calculation%fermi_surface_plot .or. &
            w90_calculation%wannier_plot) then
          write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
          write (stdout, '(1x,a)') '|                               PLOTTING                                    |'
          write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
          write (stdout, *)
        endif
      endif

      if (w90_calculation%fermi_surface_plot) then
        ! plot_fermi_surface can be trivially parallelised--fixme
        call plot_fermi_surface(fermi_energy_list, recip_lattice, fermi_surface_plot, num_wann, &
                                ham_r, irvec, ndegen, nrpts, print_output%timing_level, stdout, &
                                seedname, timer, error, comm)
        if (allocated(error)) return
      endif

      if (output_file%write_hr .or. output_file%write_tb) then
        call ws_translate_dist(ws_distance, ws_region, num_wann, &
                               wannier_data%centres, real_lattice, mp_grid, nrpts, irvec, &
                               error, comm, force_recompute=.false.)
        if (allocated(error)) return

        call ws_write_vec(ws_distance, nrpts, irvec, num_wann, ws_region%use_ws_distance, &
                          seedname, error, comm)
        if (allocated(error)) return
      end if

      ! calculate and write projection of WFs on original bands in outer window
      ! only meaningful in disentanglement case (and otherwise lwindow is not available)
      if (output_file%write_proj .and. num_bands > num_wann) then
        call plot_calc_projection(num_bands, num_wann, num_kpts, u_matrix_opt, eigval, &
                                  dis_manifold%lwindow, print_output%timing_level, &
                                  print_output%iprint, stdout, timer)
        if (allocated(error)) return
      endif

      if (output_file%write_bvec) then
        call plot_bvec(kmesh_info, num_kpts, seedname, error, comm)
        if (allocated(error)) return
      endif

      if (output_file%write_hr) then
        ! this is a trivial matrix write; no need to parallelize
        call hamiltonian_write_hr(ham_r, irvec, ndegen, nrpts, num_wann, &
                                  print_output%timing_level, seedname, timer, error, comm)
        if (allocated(error)) return
      endif

      if (output_file%write_hr_diag) then
        if (print_output%iprint > 0) then
          write (stdout, *)
          write (stdout, '(1x,a)') 'On-site Hamiltonian matrix elements'
          write (stdout, '(3x,a)') '  n        <0n|H|0n> (eV)'
          write (stdout, '(3x,a)') '-------------------------'
          do i = 1, num_wann
            write (stdout, '(3x,i3,5x,f12.6)') i, real(ham_r(i, i, rpt_origin), kind=dp)
          enddo
          write (stdout, *)
        endif
      endif

      if (output_file%write_u_matrices) then
        call plot_u_matrices(u_matrix_opt, u_matrix, kpt_latt, dis_manifold, have_disentangled, &
                             num_wann, num_kpts, num_bands, seedname)
      endif

      if (output_file%write_vdw_data) then
        ! aam: write data required for vdW utility
        call plot_write_vdw_data(num_wann, wannier_data, real_lattice, u_matrix, u_matrix_opt, &
                                 have_disentangled, w90_system, error, comm, stdout, seedname)
        if (allocated(error)) return
      endif

      if (output_file%write_xyz) then
        call plot_write_xyz(real_space_ham%translate_home_cell, num_wann, wannier_data%centres, &
                            real_lattice, atom_data, print_output, error, comm, stdout, seedname)
        if (allocated(error)) return
      endif
    end if !on_root

    if (w90_calculation%bands_plot) then
      bands_num_spec_points = 0
      if (allocated(kpoint_path%labels)) bands_num_spec_points = size(kpoint_path%labels)

      call plot_interpolate_bands(mp_grid, real_lattice, band_plot, kpoint_path, real_space_ham, &
                                  ws_region, print_output, recip_lattice, num_wann, wannier_data, &
                                  ham_r, irvec, ndegen, nrpts, wannier_centres_translated, &
                                  ws_distance, bands_num_spec_points, stdout, seedname, timer, &
                                  error, comm)
      if (allocated(error)) return
    endif

    if (output_file%svd_omega) then
      call plot_svd_omega_i(num_wann, num_kpts, kmesh_info, m_matrix, print_output, timer, dist_k, &
                            error, comm, stdout)
      if (allocated(error)) return
    endif

    if (output_file%write_rmn) then
      ! parallel write_rmn
      call plot_write_rmn(kmesh_info, m_matrix, kpt_latt, irvec, nrpts, num_kpts, num_wann, &
                          seedname, dist_k, error, comm)
      if (allocated(error)) return
    endif

    if (output_file%write_r2mn) then
      ! write matrix elements <m|r^2|n> to file
      call plot_write_r2mn(num_kpts, num_wann, kmesh_info, m_matrix, seedname, dist_k, error, comm)
      if (allocated(error)) return
    endif

    if (output_file%write_tb) then
      call hamiltonian_write_tb(kmesh_info, ham_r, m_matrix, kpt_latt, real_lattice, irvec, &
                                ndegen, nrpts, num_kpts, num_wann, print_output%timing_level, &
                                seedname, timer, dist_k, error, comm)
      if (allocated(error)) return
    endif

    if (w90_calculation%wannier_plot) then
      call plot_wannier(wannier_plot, wvfn_read, wannier_data, print_output, u_matrix_opt, &
                        dis_manifold, real_lattice, atom_data, kpt_latt, u_matrix, num_kpts, &
                        num_bands, num_wann, have_disentangled, w90_system%spinors, bohr, stdout, &
                        seedname, timer, dist_k, error, comm)
      if (allocated(error)) return
    endif

    if (on_root .and. print_output%timing_level > 0) call io_stopwatch_stop('plot: main', timer)
  end subroutine plot_main

  !-----------------------------------------------------------------
  !----------------- Private Routines ------------------------------
  !-----------------------------------------------------------------

  !================================================!
  subroutine plot_interpolate_bands(mp_grid, real_lattice, band_plot, kpoint_path, real_space_ham, &
                                    ws_region, print_output, recip_lattice, num_wann, &
                                    wannier_data, ham_r, irvec, ndegen, nrpts, &
                                    wannier_centres_translated, ws_distance, &
                                    bands_num_spec_points, stdout, seedname, timer, error, comm)
    !================================================!
    !                                            !
    !! Plots the interpolated band structure
    !                                            !
    !================================================!

    use w90_constants, only: dp, cmplx_0, twopi
    use w90_io, only: io_time, io_stopwatch_start, io_stopwatch_stop
    use w90_ws_distance, only: ws_translate_dist
    use w90_utility, only: utility_metric
    use w90_types, only: wannier_data_type, kpoint_path_type, print_output_type, ws_region_type, &
      ws_distance_type, timer_list_type
    use w90_wannier90_types, only: band_plot_type, real_space_ham_type
    use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
      set_error_warn

    implicit none

    ! arguments
    type(band_plot_type), intent(in) :: band_plot
    type(kpoint_path_type), intent(in) :: kpoint_path
    type(print_output_type), intent(in) :: print_output
    type(real_space_ham_type), intent(in) :: real_space_ham
    type(wannier_data_type), intent(in) :: wannier_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(ws_region_type), intent(in) :: ws_region
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(inout) :: nrpts
    integer, intent(in) :: ndegen(:)
    integer, intent(in) :: irvec(:, :)
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann
    integer, intent(in) :: bands_num_spec_points
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) :: wannier_centres_translated(:, :)

    complex(kind=dp), intent(in) :: ham_r(:, :, :)
    character(len=50), intent(in)  :: seedname

    ! local variables
    integer, allocatable :: irvec_cut(:, :)
    integer              :: irvec_max(3)
    integer              :: nrpts_cut
    integer, allocatable :: iwork(:), ifail(:)
    integer              :: info, i, j
    integer              :: irpt, nfound, loop_kpt, counter
    integer              :: loop_spts, total_pts, loop_i, nkp, ideg
    integer              :: num_paths, num_spts, ierr
    integer              :: bndunit, gnuunit, loop_w, loop_p
    integer, allocatable :: kpath_pts(:)
    integer, allocatable :: idx_special_points(:)
    integer, allocatable :: label_idx_special_points(:)

    real(kind=dp)              :: rdotk, vec(3), emin, emax, time0
    real(kind=dp), allocatable :: kpath_len(:)
    real(kind=dp), allocatable :: rwork(:)
    real(kind=dp), allocatable :: xval(:)
    real(kind=dp), allocatable :: eig_int(:, :), plot_kpoint(:, :)
    real(kind=dp), allocatable :: bands_proj(:, :)
    real(kind=dp), allocatable :: xval_special_points(:)
    real(kind=dp)              :: recip_metric(3, 3)

    complex(kind=dp)              :: fac
    complex(kind=dp), allocatable :: ham_r_cut(:, :, :)
    complex(kind=dp), allocatable :: ham_pack(:)
    complex(kind=dp), allocatable :: ham_kprm(:, :)
    complex(kind=dp), allocatable :: U_int(:, :)
    complex(kind=dp), allocatable :: cwork(:)

    logical, allocatable :: kpath_print_first_point(:)

    character(len=20), allocatable :: glabel(:)
    character(len=20), allocatable :: xlabel(:)
    character(len=20), allocatable :: ctemp(:)

    ! mpi variables
    integer ::  my_node_id, num_nodes
    logical ::  on_root
    integer, allocatable :: counts(:)
    integer, allocatable :: displs(:)

    num_nodes = mpisize(comm)
    my_node_id = mpirank(comm)

    on_root = .false.
    if (my_node_id == 0) on_root = .true.

    allocate (counts(0:num_nodes - 1))
    allocate (displs(0:num_nodes - 1))
    !
    if (on_root) then
      if (print_output%timing_level > 1) then
        call io_stopwatch_start('plot: interpolate_bands', timer)
      endif
      !
      time0 = io_time()

      write (stdout, *)
      write (stdout, '(1x,a)') 'Calculating interpolated band-structure'
      write (stdout, *)
    endif ! on_root
    call utility_metric(recip_lattice, recip_metric)
    !
    allocate (ham_pack((num_wann*(num_wann + 1))/2), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ham_pack in plot_interpolate_bands', comm)
      return
    endif
    allocate (ham_kprm(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ham_kprm in plot_interpolate_bands', comm)
      return
    endif
    allocate (U_int(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating U_int in plot_interpolate_bands', comm)
      return
    endif
    allocate (cwork(2*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cwork in plot_interpolate_bands', comm)
      return
    endif
    allocate (rwork(7*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating rwork in plot_interpolate_bands', comm)
      return
    endif
    allocate (iwork(5*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating iwork in plot_interpolate_bands', comm)
      return
    endif
    allocate (ifail(num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ifail in plot_interpolate_bands', comm)
      return
    endif

    if (.not. kpoint_path%bands_kpt_explicit) then
      allocate (kpath_len(bands_num_spec_points/2), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating kpath_len in plot_interpolate_bands', comm)
        return
      endif
      allocate (kpath_pts(bands_num_spec_points/2), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating kpath_pts in plot_interpolate_bands', comm)
        return
      endif
      allocate (kpath_print_first_point(bands_num_spec_points/2), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating kpath_print_first_point in plot_interpolate_bands', comm)
        return
      endif
      allocate (idx_special_points(bands_num_spec_points), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating idx_special_points in plot_interpolate_bands', comm)
        return
      endif
      allocate (xval_special_points(bands_num_spec_points), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating xval_special_points in plot_interpolate_bands', comm)
        return
      endif
      idx_special_points = -1
      xval_special_points = -1._dp
    end if
    !
    ! Work out how many points in the total path and the positions of the special points
    !
    if (kpoint_path%bands_kpt_explicit) then
      total_pts = size(kpoint_path%bands_kpt_frac, 2)
      ! Count the total number of special points
      num_spts = 0
      do i = 1, total_pts
        do j = 1, bands_num_spec_points
          if (sum((kpoint_path%bands_kpt_frac(:, i) - kpoint_path%points(:, j))**2) <= 1.e-6) then
            num_spts = num_spts + 1
            exit
          end if
        end do
      end do
    else
      num_paths = bands_num_spec_points/2

      kpath_print_first_point = .false.

      ! Loop over paths, set to False print_first_point if the starting point
      ! is the same as the ending point of the previous path.
      ! I skip the first path for which I always want to print the first point.
      kpath_print_first_point(1) = .true.
      do i = 2, num_paths
        ! If either the coordinates are different or the label is different, compute again the point
        ! (it will end up at the same x coordinate)
        if ((SUM((kpoint_path%points(:, (i - 1)*2) - &
                  kpoint_path%points(:, (i - 1)*2 + 1))**2) > 1.e-6) .or. &
            (TRIM(kpoint_path%labels((i - 1)*2)) .ne. &
             TRIM(kpoint_path%labels((i - 1)*2 + 1)))) then
          kpath_print_first_point(i) = .true.
        end if
      enddo

      ! Count the total number of special points
      num_spts = num_paths
      do i = 1, num_paths
        if (kpath_print_first_point(i)) num_spts = num_spts + 1
      end do

      do loop_spts = 1, num_paths
        vec = kpoint_path%points(:, 2*loop_spts) - kpoint_path%points(:, 2*loop_spts - 1)
        kpath_len(loop_spts) = sqrt(dot_product(vec, (matmul(recip_metric, vec))))
        if (loop_spts == 1) then
          kpath_pts(loop_spts) = kpoint_path%num_points_first_segment
        else
          kpath_pts(loop_spts) = nint(real(kpoint_path%num_points_first_segment, dp) &
                                      *kpath_len(loop_spts)/kpath_len(1))
          ! At least 1 point
          !if (kpath_pts(loop_spts) .eq. 0) kpath_pts(loop_spts) = 1
        end if
      end do
      total_pts = sum(kpath_pts)
      do i = 1, num_paths
        if (kpath_print_first_point(i)) total_pts = total_pts + 1
      end do
    end if

    allocate (plot_kpoint(3, total_pts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating plot_kpoint in plot_interpolate_bands', comm)
      return
    endif
    allocate (xval(total_pts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating xval in plot_interpolate_bands', comm)
      return
    endif
    allocate (eig_int(num_wann, total_pts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating eig_int in plot_interpolate_bands', comm)
      return
    endif
    allocate (bands_proj(num_wann, total_pts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating bands_proj in plot_interpolate_bands', comm)
      return
    endif
    allocate (glabel(num_spts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating glabel in plot_interpolate_bands', comm)
      return
    endif
    allocate (xlabel(num_spts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating xlabel in plot_interpolate_bands', comm)
      return
    endif
    allocate (ctemp(bands_num_spec_points), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ctemp in plot_interpolate_bands', comm)
      return
    endif
    eig_int = 0.0_dp; bands_proj = 0.0_dp
    !
    ! Find the position of each kpoint in the path
    !
    if (kpoint_path%bands_kpt_explicit) then
      allocate (idx_special_points(num_spts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating idx_special_points in plot_interpolate_bands', comm)
        return
      endif
      allocate (xval_special_points(num_spts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating xval_special_points in plot_interpolate_bands', comm)
        return
      endif
      idx_special_points = -1
      xval_special_points = -1._dp
      allocate (label_idx_special_points(num_spts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating label_idx_special_points in plot_interpolate_bands', comm)
        return
      endif
      plot_kpoint(:, :) = kpoint_path%bands_kpt_frac(:, :)
      xval = 0.0_dp
      counter = 0
      do i = 1, total_pts
        if (i > 1) then
          vec = plot_kpoint(:, i) - plot_kpoint(:, i - 1)
          xval(i) = xval(i - 1) + sqrt(dot_product(vec, (matmul(recip_metric, vec))))
        end if
        do j = 1, bands_num_spec_points
          if (sum((kpoint_path%bands_kpt_frac(:, i) - kpoint_path%points(:, j))**2) <= 1.e-6) then
            counter = counter + 1
            idx_special_points(counter) = i
            label_idx_special_points(counter) = j
            xval_special_points(counter) = xval(i)
            if (counter > 1) then
              if (idx_special_points(counter) == idx_special_points(counter - 1) + 1) then
                ! If the two points are consecutive, it means that the x coordinate should be the same
                xval(i) = xval(i - 1)
                xval_special_points(counter) = xval(i)
              end if
            end if
            exit
          end if
        end do
      end do

    else
      counter = 0
      do loop_spts = 1, num_paths
        if (kpath_print_first_point(loop_spts)) then
          counter = counter + 1
          if (counter == 1) then
            xval(counter) = 0.0_dp
          else
            ! If we are printing the first point in a path,
            ! It means that the coordinate did not change (otherwise
            ! we would not be printing it). Therefore I do not move
            ! on the x axis, there was a jump in the path here.
            xval(counter) = xval(counter - 1)
          endif
          plot_kpoint(:, counter) = kpoint_path%points(:, 2*loop_spts - 1)

          idx_special_points(2*loop_spts - 1) = counter
          xval_special_points(2*loop_spts - 1) = xval(counter)
        end if

        ! This is looping on all points but the first (1 is the first point
        ! after the first in the path)
        do loop_i = 1, kpath_pts(loop_spts)
          counter = counter + 1
          ! Set xval, the x position on the path of the current path
          if (counter == 1) then
            ! This case should never happen but I keep it in for "safety"
            xval(counter) = 0.0_dp
          else
            xval(counter) = xval(counter - 1) + kpath_len(loop_spts)/real(kpath_pts(loop_spts), dp)
          endif
          plot_kpoint(:, counter) = kpoint_path%points(:, 2*loop_spts - 1) + &
                                    (kpoint_path%points(:, 2*loop_spts) &
                                     - kpoint_path%points(:, 2*loop_spts - 1))* &
                                    (real(loop_i, dp)/real(kpath_pts(loop_spts), dp))
        end do
        idx_special_points(2*loop_spts) = counter
        xval_special_points(2*loop_spts) = xval(counter)
      end do
      !xval(total_pts)=sum(kpath_len)
      plot_kpoint(:, total_pts) = kpoint_path%points(:, bands_num_spec_points)
    end if
    !
    ! Write out the kpoints in the path
    !
    if (on_root) then
      open (newunit=bndunit, file=trim(seedname)//'_band.kpt', form='formatted')
      write (bndunit, *) total_pts
      do loop_spts = 1, total_pts
        write (bndunit, '(3f12.6,3x,a)') (plot_kpoint(loop_i, loop_spts), loop_i=1, 3), "1.0"
      end do
      close (bndunit)
      !
      ! Write out information on high-symmetry points in the path
      !
      open (newunit=bndunit, file=trim(seedname)//'_band.labelinfo.dat', form='formatted')
      if (kpoint_path%bands_kpt_explicit) then
        do loop_spts = 1, num_spts
          write (bndunit, '(a,3x,I10,3x,4f18.10)') &
            kpoint_path%labels(label_idx_special_points(loop_spts)), &
            idx_special_points(loop_spts), &
            xval_special_points(loop_spts), &
            (plot_kpoint(loop_i, idx_special_points(loop_spts)), loop_i=1, 3)
        end do
      else
        do loop_spts = 1, bands_num_spec_points
          if ((MOD(loop_spts, 2) .eq. 1) .and. &
              (kpath_print_first_point((loop_spts + 1)/2) .eqv. .false.)) cycle
          write (bndunit, '(a,3x,I10,3x,4f18.10)') &
            kpoint_path%labels(loop_spts), &
            idx_special_points(loop_spts), &
            xval_special_points(loop_spts), &
            (plot_kpoint(loop_i, idx_special_points(loop_spts)), loop_i=1, 3)
        end do
      end if
      close (bndunit)
    endif ! on_root
    !
    ! Cut H matrix in real-space
    !
    if (index(band_plot%mode, 'cut') .ne. 0) then
      call plot_cut_hr(real_space_ham, real_lattice, mp_grid, num_wann, &
                       wannier_centres_translated, stdout, error)
      if (allocated(error)) return
    endif
    !
    ! Interpolate the Hamiltonian at each kpoint
    !
    if (ws_region%use_ws_distance) then
      if (index(band_plot%mode, 's-k') .ne. 0) then
        call ws_translate_dist(ws_distance, ws_region, num_wann, &
                               wannier_data%centres, real_lattice, mp_grid, nrpts, &
                               irvec, error, comm, force_recompute=.true.)
        if (allocated(error)) return

      elseif (index(band_plot%mode, 'cut') .ne. 0) then
        call ws_translate_dist(ws_distance, ws_region, num_wann, &
                               wannier_data%centres, real_lattice, mp_grid, nrpts_cut, &
                               irvec_cut, error, comm, force_recompute=.true.)
        if (allocated(error)) return

      else
        call set_error_warn(error, 'Error in plot_interpolate bands: value of bands_plot_mode not recognised', comm)
        return
      endif
    endif
    if (on_root .and. print_output%timing_level > 2) then
      call io_stopwatch_start('plot: interpolate_bands: loop_kpoints', timer)
    endif
    ! [lp] the s-k and cut codes are very similar when use_ws_distance is used, a complete
    !      merge after this point is not impossible
    call comms_array_split(total_pts, counts, displs, comm)
    ! Don't worry about serial run, it is OK to call comms_array_split!
    do loop_kpt = displs(my_node_id) + 1, displs(my_node_id) + counts(my_node_id)
      ! do loop_kpt = 1, total_pts
      ham_kprm = cmplx_0
      !
      if (index(band_plot%mode, 's-k') .ne. 0) then
        do irpt = 1, nrpts
          ! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
          if (ws_region%use_ws_distance) then
            do j = 1, num_wann
            do i = 1, num_wann
              do ideg = 1, ws_distance%ndeg(i, j, irpt)
                rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), &
                                          real(ws_distance%irdist(:, ideg, i, j, irpt), dp))
                fac = cmplx(cos(rdotk), sin(rdotk), dp) &
                      /real(ndegen(irpt)*ws_distance%ndeg(i, j, irpt), dp)
                ham_kprm(i, j) = ham_kprm(i, j) + fac*ham_r(i, j, irpt)
              enddo
            enddo
            enddo
          else
            ! [lp] Original code, without IJ-dependent shift:
            rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), irvec(:, irpt))
            fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(irpt), dp)
            ham_kprm = ham_kprm + fac*ham_r(:, :, irpt)
          endif
        end do
        ! end of s-k mode
      elseif (index(band_plot%mode, 'cut') .ne. 0) then
        do irpt = 1, nrpts_cut
          ! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
          if (ws_region%use_ws_distance) then
            do j = 1, num_wann
            do i = 1, num_wann
              do ideg = 1, ws_distance%ndeg(i, j, irpt)
                rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), &
                                          real(ws_distance%irdist(:, ideg, i, j, irpt), dp))
                fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ws_distance%ndeg(i, j, irpt), dp)
                ham_kprm(i, j) = ham_kprm(i, j) + fac*ham_r_cut(i, j, irpt)
              enddo
            enddo
            enddo
            ! [lp] Original code, without IJ-dependent shift:
          else
            rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), irvec_cut(:, irpt))
            !~[aam] check divide by ndegen?
            fac = cmplx(cos(rdotk), sin(rdotk), dp)
            ham_kprm = ham_kprm + fac*ham_r_cut(:, :, irpt)
          endif ! end of use_ws_distance
        end do
      endif ! end of "cut" mode
      !
      ! Diagonalise H_k (->basis of eigenstates)
      do j = 1, num_wann
        do i = 1, j
          ham_pack(i + ((j - 1)*j)/2) = ham_kprm(i, j)
        enddo
      enddo
      call ZHPEVX('V', 'A', 'U', num_wann, ham_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                  nfound, eig_int(1, loop_kpt), U_int, num_wann, cwork, rwork, iwork, ifail, info)
      if (info < 0) then
        write (stdout, '(a,i3,a)') 'THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
        call set_error_fatal(error, 'Error in plot_interpolate_bands', comm)
        return
      endif
      if (info > 0) then
        write (stdout, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
        call set_error_warn(error, 'Error in plot_interpolate_bands', comm)
        return
      endif
      ! Compute projection onto WF if requested
      if (allocated(band_plot%project)) then
      do loop_w = 1, num_wann
        do loop_p = 1, num_wann
          if (any(band_plot%project == loop_p)) then
            bands_proj(loop_w, loop_kpt) = bands_proj(loop_w, loop_kpt) + &
                                           abs(U_int(loop_p, loop_w))**2
          end if
        end do
      end do
      end if
      !
    end do

    call comms_reduce(eig_int(1, 1), num_wann*total_pts, 'SUM', error, comm)
    call comms_reduce(bands_proj(1, 1), num_wann*total_pts, 'SUM', error, comm)

    if (on_root .and. print_output%timing_level > 2) then
      call io_stopwatch_stop('plot: interpolate_bands: loop_kpoints', timer)
    endif
    !
    ! Interpolation Finished!
    ! Now we write plotting files
    !
    if (on_root) then
      emin = minval(eig_int) - 1.0_dp
      emax = maxval(eig_int) + 1.0_dp

      if (index(band_plot%format, 'gnu') > 0) then
        call plot_interpolate_gnuplot(band_plot, kpoint_path, bands_num_spec_points, num_wann)
      endif
      if (index(band_plot%format, 'xmgr') > 0) then
        call plot_interpolate_xmgrace(kpoint_path, bands_num_spec_points, num_wann, error)
      endif
      write (stdout, '(1x,a,f11.3,a)') &
        'Time to calculate interpolated band structure ', io_time() - time0, ' (sec)'
      write (stdout, *)
    endif

    if (allocated(ham_r_cut)) then
      deallocate (ham_r_cut, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ham_r_cut in plot_interpolate_bands', comm)
        return
      endif
    endif
    if (allocated(irvec_cut)) then
      deallocate (irvec_cut, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating irvec_cut in plot_interpolate_bands', comm)
        return
      endif
    endif

    if (print_output%timing_level > 1) call io_stopwatch_stop('plot: interpolate_bands', timer)

    if (allocated(idx_special_points)) then
      deallocate (idx_special_points, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating idx_special_points in &
            &plot_interpolate_bands', comm)
        return
      endif
    endif
    if (allocated(xval_special_points)) then
      deallocate (xval_special_points, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating xval_special_points in &
            &plot_interpolate_bands', comm)
        return
      endif
    endif
    if (allocated(label_idx_special_points)) then
      deallocate (label_idx_special_points, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating label_idx_special_points in &
            &plot_interpolate_bands', comm)
        return
      endif
    endif

  contains

    !================================================!
    subroutine plot_cut_hr(real_space_ham, real_lattice, mp_grid, num_wann, &
                           wannier_centres_translated, stdout, error)
      !================================================!
      !
      !!  In real-space picture, ham_r(j,i,k) is an interaction between
      !!  j_th WF at 0 and i_th WF at the lattice point translated
      !!  by matmul(real_lattice(:,:),irvec(:,k))
      !!  We truncate Hamiltonian matrix when
      !!   1) |  r_i(0) - r_j (R) | > dist_cutoff
      !!   2) |  ham_r(i,j,k)     | < hr_cutoff
      !!  while the condition 1) is essential to get a meaningful band structure,
      !!    ( dist_cutoff must be smaller than the shortest distance from
      !!      the center of W-S supercell to the points at the cell boundaries )
      !!  the condition 2) is optional.
      !!
      !!  limitation: when bands_plot_dim .ne. 3
      !!      one_dim_vec must be parallel to one of the cartesian axis
      !!      and perpendicular to the other two primitive lattice vectors
      !================================================!

      use w90_constants, only: dp, cmplx_0, eps8
      use w90_wannier90_types, only: band_plot_type, real_space_ham_type
      use w90_error, only: w90_error_type, set_error_alloc, set_error_warn

      implicit none

      ! arguments
      type(real_space_ham_type), intent(in) :: real_space_ham
      type(w90_error_type), allocatable, intent(out) :: error

      real(kind=dp), intent(in) :: real_lattice(3, 3)
      real(kind=dp), intent(in) :: wannier_centres_translated(:, :)

      integer, intent(in) :: mp_grid(3)
      integer, intent(in) :: num_wann
      integer, intent(in) :: stdout

      ! local variables
      integer :: nrpts_tmp
      integer :: one_dim_vec, two_dim_vec(2)
      integer :: i, j, n1, n2, n3, i1, i2, i3
      real(kind=dp), allocatable :: ham_r_tmp(:, :)
      real(kind=dp), allocatable :: shift_vec(:, :)
      real(kind=dp) :: dist_ij_vec(3)
      real(kind=dp) :: dist_vec(3)
      real(kind=dp) :: dist

      allocate (ham_r_tmp(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating ham_r_tmp in plot_cut_hr', comm)
        return
      endif

      irvec_max = maxval(irvec, DIM=2) + 1

      if (real_space_ham%system_dim .ne. 3) then
        ! Find one_dim_vec which is parallel to one_dim_dir
        ! two_dim_vec - the other two lattice vectors
        ! Along the confined directions, take only irvec=0
        j = 0
        do i = 1, 3
          if (abs(abs(real_lattice(real_space_ham%one_dim_dir, i)) &
                  - sqrt(dot_product(real_lattice(:, i), real_lattice(:, i)))) .lt. eps8) then
            one_dim_vec = i
            j = j + 1
          end if
        end do
        if (j .ne. 1) then
          call set_error_warn(error, 'Error: 1-d lattice vector not defined in plot_cut_hr', comm)
          return
        endif
        j = 0
        do i = 1, 3
          if (i .ne. one_dim_vec) then
            j = j + 1
            two_dim_vec(j) = i
          end if
        end do
        if (real_space_ham%system_dim .eq. 1) then
          irvec_max(two_dim_vec(1)) = 0
          irvec_max(two_dim_vec(2)) = 0
        end if
        if (real_space_ham%system_dim .eq. 2) irvec_max(one_dim_vec) = 0
      end if

      nrpts_cut = (2*irvec_max(1) + 1)*(2*irvec_max(2) + 1)*(2*irvec_max(3) + 1)
      allocate (ham_r_cut(num_wann, num_wann, nrpts_cut), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating ham_r_cut in plot_cut_hr', comm)
        return
      endif
      allocate (irvec_cut(3, nrpts_cut), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating irvec_cut in plot_cut_hr', comm)
        return
      endif
      allocate (shift_vec(3, nrpts_cut), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating shift_vec in plot_cut_hr', comm)
        return
      endif

      nrpts_tmp = 0
      do n1 = -irvec_max(1), irvec_max(1)
        do n2 = -irvec_max(2), irvec_max(2)
          loop_n3: do n3 = -irvec_max(3), irvec_max(3)
            do irpt = 1, nrpts
              i1 = mod(n1 - irvec(1, irpt), mp_grid(1))
              i2 = mod(n2 - irvec(2, irpt), mp_grid(2))
              i3 = mod(n3 - irvec(3, irpt), mp_grid(3))
              if (i1 .eq. 0 .and. i2 .eq. 0 .and. i3 .eq. 0) then
                nrpts_tmp = nrpts_tmp + 1
                ham_r_cut(:, :, nrpts_tmp) = ham_r(:, :, irpt)
                irvec_cut(1, nrpts_tmp) = n1
                irvec_cut(2, nrpts_tmp) = n2
                irvec_cut(3, nrpts_tmp) = n3
                cycle loop_n3
              end if
            end do
          end do loop_n3
        end do
      end do

      if (nrpts_tmp .ne. nrpts_cut) then
        write (stdout, '(a)') 'FAILED TO EXPAND ham_r'
        call set_error_warn(error, 'Error in plot_cut_hr', comm)
        return
      end if

      ! AAM: 29/10/2009 Bug fix thanks to Dr Shujun Hu, NIMS, Japan.
      do irpt = 1, nrpts_cut
        ! line below is incorrect for non-orthorhombic cells
        !shift_vec(:,irpt) = matmul(real_lattice(:,:),real(irvec_cut(:,irpt),dp))
        ! line below is the same as calculating
        ! matmul(transpose(real_lattice(:,:)),irvec_cut(:,irpt))
        shift_vec(:, irpt) = matmul(real(irvec_cut(:, irpt), dp), real_lattice(:, :))
      end do

      ! note: dist_cutoff_mode does not necessarily follow bands_plot_dim
      ! e.g. for 1-d system (bands_plot_dim=1) we can still apply 3-d dist_cutoff (dist_cutoff_mode=three_dim)
      if (index(real_space_ham%dist_cutoff_mode, 'one_dim') > 0) then
        do i = 1, num_wann
          do j = 1, num_wann
            dist_ij_vec(real_space_ham%one_dim_dir) = &
              wannier_centres_translated(real_space_ham%one_dim_dir, i) - &
              wannier_centres_translated(real_space_ham%one_dim_dir, j)
            do irpt = 1, nrpts_cut
              dist_vec(real_space_ham%one_dim_dir) = dist_ij_vec(real_space_ham%one_dim_dir) + &
                                                     shift_vec(real_space_ham%one_dim_dir, irpt)
              dist = abs(dist_vec(real_space_ham%one_dim_dir))
              if (dist .gt. real_space_ham%dist_cutoff) &
                ham_r_cut(j, i, irpt) = cmplx_0
            end do
          end do
        end do
      else if (index(real_space_ham%dist_cutoff_mode, 'two_dim') > 0) then
        do i = 1, num_wann
          do j = 1, num_wann
            dist_ij_vec(:) = wannier_centres_translated(:, i) - wannier_centres_translated(:, j)
            do irpt = 1, nrpts_cut
              dist_vec(:) = dist_ij_vec(:) + shift_vec(:, irpt)
              dist_vec(real_space_ham%one_dim_dir) = 0.0_dp
              dist = sqrt(dot_product(dist_vec, dist_vec))
              if (dist .gt. real_space_ham%dist_cutoff) &
                ham_r_cut(j, i, irpt) = cmplx_0
            end do
          end do
        end do
      else
        do i = 1, num_wann
          do j = 1, num_wann
            dist_ij_vec(:) = wannier_centres_translated(:, i) - wannier_centres_translated(:, j)
            do irpt = 1, nrpts_cut
              dist_vec(:) = dist_ij_vec(:) + shift_vec(:, irpt)
              dist = sqrt(dot_product(dist_vec, dist_vec))
              if (dist .gt. real_space_ham%dist_cutoff) &
                ham_r_cut(j, i, irpt) = cmplx_0
            end do
          end do
        end do
      end if

      do irpt = 1, nrpts_cut
        do i = 1, num_wann
          do j = 1, num_wann
            if (abs(ham_r_cut(j, i, irpt)) .lt. real_space_ham%hr_cutoff) &
              ham_r_cut(j, i, irpt) = cmplx_0
          end do
        end do
      end do

      write (stdout, '(/1x,a78)') repeat('-', 78)
      write (stdout, '(1x,4x,a)') &
        'Maximum absolute value of Real-space Hamiltonian at each lattice point'
      write (stdout, '(1x,8x,a62)') repeat('-', 62)
      write (stdout, '(1x,11x,a,11x,a)') 'Lattice point R', 'Max |H_ij(R)|'
      !  output maximum ham_r_cut at each lattice point
      do irpt = 1, nrpts_cut
        ham_r_tmp(:, :) = abs(ham_r_cut(:, :, irpt))
        write (stdout, '(1x,9x,3I5,9x,F12.6)') irvec_cut(:, irpt), maxval(ham_r_tmp)
      end do
      !
      return

    end subroutine plot_cut_hr

    !================================================!
    subroutine plot_interpolate_gnuplot(band_plot, kpoint_path, bands_num_spec_points, num_wann)
      !================================================!
      !
      !! Plots the interpolated band structure in gnuplot format
      !
      !================================================!

      use w90_constants, only: dp
      use w90_types, only: kpoint_path_type
      use w90_wannier90_types, only: band_plot_type

      implicit none

      ! arguments
      type(band_plot_type), intent(in) :: band_plot
      type(kpoint_path_type), intent(in) :: kpoint_path
      integer, intent(in) :: num_wann, bands_num_spec_points

      open (newunit=bndunit, file=trim(seedname)//'_band.dat', form='formatted')
      open (newunit=gnuunit, file=trim(seedname)//'_band.gnu', form='formatted')
      !
      ! Gnuplot format
      !
      do i = 1, num_wann
        do nkp = 1, total_pts
          if (allocated(band_plot%project)) then
            write (bndunit, '(3E16.8)') xval(nkp), eig_int(i, nkp), bands_proj(i, nkp)
          else
            write (bndunit, '(2E16.8)') xval(nkp), eig_int(i, nkp)
          end if
        enddo
        write (bndunit, *) ' '
      enddo
      close (bndunit)
      ! Axis labels
      if (kpoint_path%bands_kpt_explicit) then
        do i = 1, num_spts
          if (i > 1) then
            if (idx_special_points(i) .eq. idx_special_points(i - 1) + 1 .and. &
                kpoint_path%labels(label_idx_special_points(i)) .ne. kpoint_path%labels(label_idx_special_points(i - 1))) then
              ! If two different point indeces are consecutive, the label should be combined
              glabel(i) = TRIM(kpoint_path%labels(label_idx_special_points(i - 1)))//'|'// &
                          TRIM(kpoint_path%labels(label_idx_special_points(i)))
            end if
          else
            glabel(i) = TRIM(kpoint_path%labels(label_idx_special_points(i)))
          end if
        end do
      else
        glabel(1) = TRIM(kpoint_path%labels(1))
        do i = 2, num_paths
          if (kpoint_path%labels(2*(i - 1)) /= kpoint_path%labels(2*(i - 1) + 1)) then
            glabel(i) = TRIM(kpoint_path%labels(2*(i - 1)))//'|'// &
                        TRIM(kpoint_path%labels(2*(i - 1) + 1))
          else
            glabel(i) = TRIM(kpoint_path%labels(2*(i - 1)))
          end if
        end do
        glabel(num_paths + 1) = TRIM(kpoint_path%labels(2*num_paths))
      end if
      ! gnu file
      write (gnuunit, 701) xval(total_pts), emin, emax
      if (kpoint_path%bands_kpt_explicit) then
        do i = 1, num_spts
          write (gnuunit, 705) xval(idx_special_points(i)), emin, xval(idx_special_points(i)), emax
        end do
        write (gnuunit, 702, advance="no") TRIM(glabel(1)), 0.0_dp, &
          (TRIM(glabel(i)), xval(idx_special_points(i)), i=2, num_spts - 1)
        write (gnuunit, 703) TRIM(glabel(num_spts)), xval(idx_special_points(num_spts))
      else
        do i = 1, num_paths - 1
          write (gnuunit, 705) sum(kpath_len(1:i)), emin, sum(kpath_len(1:i)), emax
        enddo
        write (gnuunit, 702, advance="no") TRIM(glabel(1)), 0.0_dp, &
          (TRIM(glabel(i + 1)), sum(kpath_len(1:i)), i=1, bands_num_spec_points/2 - 1)
        write (gnuunit, 703) TRIM(glabel(1 + bands_num_spec_points/2)), sum(kpath_len(:))
      end if
      write (gnuunit, *) 'plot ', '"'//trim(seedname)//'_band.dat', '"'
      close (gnuunit)

      if (allocated(band_plot%project)) then
        open (newunit=gnuunit, file=trim(seedname)//'_band_proj.gnu', form='formatted')
        write (gnuunit, '(a)') '#File to plot a colour-mapped Bandstructure'
        write (gnuunit, '(a)') 'set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )'
        write (gnuunit, '(a)') 'unset ztics'
        write (gnuunit, '(a)') 'unset key'
        write (gnuunit, '(a)') '# can make pointsize smaller (~0.5). Too small and nothing is &
          &printed'
        write (gnuunit, '(a)') 'set pointsize 0.8'
        write (gnuunit, '(a)') 'set grid xtics'
        write (gnuunit, '(a)') 'set view 0,0'
        write (gnuunit, '(a,f9.5,a)') 'set xrange [0:', xval(total_pts), ']'
        write (gnuunit, '(a,f9.5,a,f9.5,a)') 'set yrange [', emin, ':', emax, ']'
        if (kpoint_path%bands_kpt_explicit) then
          write (gnuunit, 702, advance="no") TRIM(glabel(1)), 0.0_dp, &
            (TRIM(glabel(i)), xval(idx_special_points(i)), i=2, num_spts - 1)
          write (gnuunit, 703) TRIM(glabel(num_spts)), xval(idx_special_points(num_spts))
        else
          write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
            (glabel(i + 1), sum(kpath_len(1:i)), i=1, bands_num_spec_points/2 - 1)
          write (gnuunit, 703) glabel(1 + bands_num_spec_points/2), sum(kpath_len(:))
        end if

        write (gnuunit, '(a,a,a,a)') 'splot ', '"'//trim(seedname)//'_band.dat', '"', &
          ' u 1:2:3 w p pt 13 palette'
        write (gnuunit, '(a)') '#use the next lines to make a nice figure for a paper'
        write (gnuunit, '(a)') '#set term postscript enhanced eps color lw 0.5 dl 0.5'
        write (gnuunit, '(a)') '#set pointsize 0.275'
      end if
      !
701   format('set style data dots', /, 'set nokey', /, 'set xrange [0:', F8.5, ']', /, &
             'set yrange [', F9.5, ' :', F9.5, ']')
702   format('set xtics (', :20('"', A, '" ', F8.5, ','))
703   format(A, '" ', F8.5, ')')
705   format('set arrow from ', F8.5, ',', F10.5, ' to ', F8.5, ',', F10.5, ' nohead')

    end subroutine plot_interpolate_gnuplot

    !================================================!
    subroutine plot_interpolate_xmgrace(kpoint_path, bands_num_spec_points, num_wann, error)
      !================================================!
      !
      !! Plots the interpolated band structure in Xmgrace format
      !
      !================================================!

      use w90_io, only: io_date
      use w90_types, only: kpoint_path_type
      use w90_error, only: w90_error_type, set_error_warn

      implicit none

      type(kpoint_path_type), intent(in) :: kpoint_path
      type(w90_error_type), allocatable, intent(out) :: error

      integer, intent(in) :: num_wann, bands_num_spec_points
      character(len=9) :: cdate, ctime

      if (kpoint_path%bands_kpt_explicit) then
        call set_error_warn(error, "bands_kpt_explicit not implemented with xmgrace format", comm)
      end if

      call io_date(cdate, ctime)

      ! Axis labels

      ! Switch any G to Gamma

      do i = 1, bands_num_spec_points
        if (kpoint_path%labels(i) == 'G') then
          ctemp(i) = '\xG\0'
        else
          ctemp(i) = kpoint_path%labels(i)
        end if
      end do

      xlabel(1) = ' '//trim(ctemp(1))//' '
      do i = 2, num_paths
        if (ctemp(2*(i - 1)) /= ctemp(2*(i - 1) + 1)) then
          xlabel(i) = trim(ctemp(2*(i - 1)))//'|'//trim(ctemp(2*(i - 1) + 1))
        else
          xlabel(i) = ctemp(2*(i - 1))
        end if
      end do
      xlabel(num_paths + 1) = ctemp(bands_num_spec_points)

      open (newunit=gnuunit, file=trim(seedname)//'_band.agr', form='formatted')
      !
      ! Xmgrace format
      !
      write (gnuunit, '(a)') '# Grace project file                      '
      write (gnuunit, '(a)') '# written using Wannier90 www.wannier.org '
      write (gnuunit, '(a)') '@version 50113                            '
      write (gnuunit, '(a)') '@page size 792, 612                       '
      write (gnuunit, '(a)') '@page scroll 5%                           '
      write (gnuunit, '(a)') '@page inout 5%                            '
      write (gnuunit, '(a)') '@link page off                            '
      write (gnuunit, '(a)') '@timestamp def "'//cdate//' at '//ctime//'" '
      write (gnuunit, '(a)') '@with g0'
      write (gnuunit, '(a)') '@    world xmin 0.00'
      write (gnuunit, '(a,f10.5)') '@    world xmax ', xval(total_pts)
      write (gnuunit, '(a,f10.5)') '@    world ymin ', emin
      write (gnuunit, '(a,f10.5)') '@    world ymax ', emax
      write (gnuunit, '(a)') '@default linewidth 1.5'
      write (gnuunit, '(a)') '@    xaxis  tick on'
      write (gnuunit, '(a)') '@    xaxis  tick major 1'
      write (gnuunit, '(a)') '@    xaxis  tick major color 1'
      write (gnuunit, '(a)') '@    xaxis  tick major linestyle 3'
      write (gnuunit, '(a)') '@    xaxis  tick major grid on'
      write (gnuunit, '(a)') '@    xaxis  tick spec type both'
      write (gnuunit, '(a,i0)') '@    xaxis  tick spec ', 1 + bands_num_spec_points/2
      write (gnuunit, '(a)') '@    xaxis  tick major 0, 0'
      do i = 1, bands_num_spec_points/2
        write (gnuunit, '(a,i0,a,a)') '@    xaxis  ticklabel ', i - 1, ',', '"'// &
          trim(adjustl(xlabel(i)))//'"'
        write (gnuunit, '(a,i0,a,f10.5)') '@    xaxis  tick major ', i, ' , ', sum(kpath_len(1:i))
      end do
      write (gnuunit, '(a,i0,a)') '@    xaxis  ticklabel ', bands_num_spec_points/2 &
        , ',"'//trim(adjustl(xlabel(1 + bands_num_spec_points/2)))//'"'
      write (gnuunit, '(a)') '@    xaxis  ticklabel char size 1.500000'
      write (gnuunit, '(a)') '@    yaxis  tick major 10'
      write (gnuunit, '(a)') '@    yaxis  label "Band Energy (eV)"'
      write (gnuunit, '(a)') '@    yaxis  label char size 1.500000'
      write (gnuunit, '(a)') '@    yaxis  ticklabel char size 1.500000'
      do i = 1, num_wann
        write (gnuunit, '(a,i0,a)') '@    s', i - 1, ' line color 1'
      end do
      do i = 1, num_wann
        write (gnuunit, '(a,i0)') '@target G0.S', i - 1
        write (gnuunit, '(a)') '@type xy'
        do nkp = 1, total_pts
          write (gnuunit, '(2E16.8)') xval(nkp), eig_int(i, nkp)
        end do
        write (gnuunit, '(a,i0)') '&'
      end do

    end subroutine plot_interpolate_xmgrace

  end subroutine plot_interpolate_bands

  !================================================!
  subroutine plot_fermi_surface(fermi_energy_list, recip_lattice, fermi_surface_plot, num_wann, &
                                ham_r, irvec, ndegen, nrpts, timing_level, stdout, seedname, &
                                timer, error, comm)
    !================================================!
    !
    !!  Prepares a Xcrysden bxsf file to view the fermi surface
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, twopi
    use w90_io, only: io_date, io_time, io_stopwatch_start, io_stopwatch_stop
    use w90_wannier90_types, only: fermi_surface_plot_type
    use w90_error, only: w90_error_type, set_error_alloc, set_error_fatal, set_error_warn
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(fermi_surface_plot_type), intent(in)   :: fermi_surface_plot
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    complex(kind=dp), intent(in) :: ham_r(:, :, :)
    character(len=50), intent(in)  :: seedname
    real(kind=dp), allocatable, intent(in)      :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    integer, intent(in) :: irvec(:, :)
    integer, intent(in) :: ndegen(:)
    integer, intent(in) :: nrpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level

    ! local variables
    complex(kind=dp), allocatable :: ham_pack(:)
    complex(kind=dp)   :: fac
    complex(kind=dp), allocatable :: ham_kprm(:, :)
    complex(kind=dp), allocatable :: U_int(:, :)
    complex(kind=dp), allocatable :: cwork(:)
    real(kind=dp), allocatable    :: rwork(:)
    real(kind=dp), allocatable  :: eig_int(:, :)
    real(kind=dp)      :: rdotk, time0
    integer, allocatable :: iwork(:), ifail(:)
    integer              :: loop_x, loop_y, loop_z, INFO, ikp, i, j, ierr
    integer              :: irpt, nfound, npts_plot, loop_kpt, bxsf_unit
    integer              :: fermi_n
    character(len=9)     :: cdate, ctime

    if (timing_level > 1) call io_stopwatch_start('plot: fermi_surface', timer)
    time0 = io_time()
    write (stdout, *)
    write (stdout, '(1x,a)') 'Calculating Fermi surface'
    write (stdout, *)

    fermi_n = 0
    if (allocated(fermi_energy_list)) fermi_n = size(fermi_energy_list)
    if (fermi_n > 1) then
      call set_error_alloc(error, "Error in plot: nfermi>1. Set the fermi level " &
                           //"using the input parameter 'fermi_level'", comm)
      return
    endif

    allocate (ham_pack((num_wann*(num_wann + 1))/2), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ham_pack plot_fermi_surface', comm)
      return
    endif
    allocate (ham_kprm(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ham_kprm plot_fermi_surface', comm)
      return
    endif
    allocate (U_int(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating U_int in plot_fermi_surface', comm)
      return
    endif
    allocate (cwork(2*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cwork in plot_fermi_surface', comm)
      return
    endif
    allocate (rwork(7*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating rwork in plot_fermi_surface', comm)
      return
    endif
    allocate (iwork(5*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating iwork in plot_fermi_surface', comm)
      return
    endif
    allocate (ifail(num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ifail in plot_fermi_surface', comm)
      return
    endif

    npts_plot = (fermi_surface_plot%num_points + 1)**3
    allocate (eig_int(num_wann, npts_plot), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating eig_int in plot_fermi_surface', comm)
      return
    endif
    eig_int = 0.0_dp
    U_int = (0.0_dp, 0.0_dp)

    ikp = 0
    do loop_x = 1, fermi_surface_plot%num_points + 1
      do loop_y = 1, fermi_surface_plot%num_points + 1
        do loop_z = 1, fermi_surface_plot%num_points + 1
          ikp = ikp + 1

          ham_kprm = cmplx_0
          do irpt = 1, nrpts
            rdotk = twopi*real((loop_x - 1)*irvec(1, irpt) + &
                               (loop_y - 1)*irvec(2, irpt) + (loop_z - 1)* &
                               irvec(3, irpt), dp)/real(fermi_surface_plot%num_points, dp)
            fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(irpt), dp)
            ham_kprm = ham_kprm + fac*ham_r(:, :, irpt)
          end do
          ! Diagonalise H_k (->basis of eigenstates)
          do j = 1, num_wann
            do i = 1, j
              ham_pack(i + ((j - 1)*j)/2) = ham_kprm(i, j)
            enddo
          enddo
          call ZHPEVX('N', 'A', 'U', num_wann, ham_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                      nfound, eig_int(1, ikp), U_int, num_wann, cwork, rwork, iwork, ifail, info)
          if (info < 0) then
            write (stdout, '(a,i3,a)') 'THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
            call set_error_fatal(error, 'Error in plot_fermi_surface', comm)
            return
          endif
          if (info > 0) then
            write (stdout, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
            call set_error_warn(error, 'Error in plot_fermi_surface', comm)
            return
          endif
        end do
      end do
    end do

    call io_date(cdate, ctime)
    open (newunit=bxsf_unit, FILE=trim(seedname)//'.bxsf', STATUS='UNKNOWN', FORM='FORMATTED')
    write (bxsf_unit, *) ' BEGIN_INFO'
    write (bxsf_unit, *) '      #'
    write (bxsf_unit, *) '      # this is a Band-XCRYSDEN-Structure-File'
    write (bxsf_unit, *) '      # for Fermi Surface Visualisation'
    write (bxsf_unit, *) '      #'
    write (bxsf_unit, *) '      # Generated by the Wannier90 code http://www.wannier.org'
    write (bxsf_unit, *) '      # On ', cdate, ' at ', ctime
    write (bxsf_unit, *) '      #'
    write (bxsf_unit, *) '      Fermi Energy:', fermi_energy_list(1)
    write (bxsf_unit, *) ' END_INFO'
    write (bxsf_unit, *)
    write (bxsf_unit, *) ' BEGIN_BLOCK_BANDGRID_3D'
    write (bxsf_unit, *) 'from_wannier_code'
    write (bxsf_unit, *) ' BEGIN_BANDGRID_3D_fermi'
    write (bxsf_unit, *) num_wann
    write (bxsf_unit, *) fermi_surface_plot%num_points + 1, fermi_surface_plot%num_points + 1, &
      fermi_surface_plot%num_points + 1
    write (bxsf_unit, *) '0.0 0.0 0.0'
    write (bxsf_unit, *) (recip_lattice(1, i), i=1, 3)
    write (bxsf_unit, *) (recip_lattice(2, i), i=1, 3)
    write (bxsf_unit, *) (recip_lattice(3, i), i=1, 3)
    do i = 1, num_wann
      write (bxsf_unit, *) 'BAND: ', i
      do loop_kpt = 1, npts_plot
        write (bxsf_unit, '(2E16.8)') eig_int(i, loop_kpt)
      enddo
    enddo
    write (bxsf_unit, *) 'END_BANDGRID_3D'
    write (bxsf_unit, *) ' END_BLOCK_BANDGRID_3D'
    close (bxsf_unit)

    write (stdout, '(1x,a,f11.3,a)') 'Time to calculate interpolated Fermi surface ', &
      io_time() - time0, ' (sec)'
    write (stdout, *)

    if (timing_level > 1) call io_stopwatch_stop('plot: fermi_surface', timer)

    return

  end subroutine plot_fermi_surface

  !================================================!
  subroutine plot_wannier(wannier_plot, wvfn_read, wannier_data, print_output, u_matrix_opt, &
                          dis_manifold, real_lattice, atom_data, kpt_latt, u_matrix, num_kpts, &
                          num_bands, num_wann, have_disentangled, spinors, bohr, stdout, seedname, &
                          timer, dist_k, error, comm)
    !================================================!
    !! Plot the WF in Xcrysden format
    !! based on code written by Michel Posternak
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi, cmplx_1
    use w90_io, only: io_date, io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: wannier_data_type, atom_data_type, dis_manifold_type, print_output_type, &
      timer_list_type
    use w90_wannier90_types, only: wvfn_read_type, wannier_plot_type
    use w90_comms, only: w90_comm_type
    use w90_error, only: w90_error_type, set_error_alloc, set_error_file, set_error_file, &
      set_error_warn

    implicit none

    ! arguments
    type(atom_data_type), intent(in) :: atom_data
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(print_output_type), intent(in) :: print_output
    type(wannier_data_type), intent(in) :: wannier_data
    type(wannier_plot_type), intent(in) :: wannier_plot
    type(wvfn_read_type), intent(in) :: wvfn_read
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)

    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    integer, intent(in) :: dist_k(:)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    logical, intent(in) :: have_disentangled
    logical, intent(in) :: spinors

    character(len=50), intent(in)  :: seedname

    ! local variables
    real(kind=dp) :: scalfac, tmax, tmaxx, x_0ang, y_0ang, z_0ang
    real(kind=dp) :: fxcry(3), dirl(3, 3), w_real, w_imag, ratmax, ratio
    real(kind=dp) :: upspinor, dnspinor, upphase, dnphase

    complex(kind=dp), allocatable :: wann_func(:, :, :, :)
    complex(kind=dp), allocatable :: r_wvfn(:, :)
    complex(kind=dp), allocatable :: r_wvfn_tmp(:, :)
    complex(kind=dp), allocatable :: wann_func_nc(:, :, :, :, :) ! add the spinor dim.
    complex(kind=dp), allocatable :: r_wvfn_nc(:, :, :) ! add the spinor dim.
    complex(kind=dp), allocatable :: r_wvfn_tmp_nc(:, :, :) ! add the spinor dim.
    complex(kind=dp) :: catmp, wmod

    logical :: have_file, on_root

    integer :: num_nodes, my_node_id
    integer :: i, j, nsp, nat, nbnd, counter, ierr
    integer :: loop_kpt, ik, ix, iy, iz, nk, ngx, ngy, ngz, nxx, nyy, nzz
    integer :: loop_b, nx, ny, nz, npoint, file_unit, loop_w, num_inc
    integer :: wann_plot_num

    character(len=11) :: wfnname
    character(len=60) :: wanxsf, wancube
    character(len=9)  :: cdate, ctime
    logical           :: inc_band(num_bands)

    num_nodes = mpisize(comm)
    my_node_id = mpirank(comm)

    on_root = .false.
    if (my_node_id == 0) on_root = .true.

    !
    if (print_output%timing_level > 1) call io_stopwatch_start('plot: wannier', timer)
    !
    associate (ngs=>wannier_plot%supercell)
      !
      if (.not. spinors) then
        write (wfnname, 200) 1, wvfn_read%spin_channel
      else
        write (wfnname, 199) 1
      endif
      inquire (file=wfnname, exist=have_file)
      if (.not. have_file) then
        call set_error_file(error, 'plot_wannier: file '//wfnname//' not found', comm)
        return
      endif

      if (wvfn_read%formatted) then
        open (newunit=file_unit, file=wfnname, form='formatted')
        read (file_unit, *) ngx, ngy, ngz, nk, nbnd
      else
        open (newunit=file_unit, file=wfnname, form='unformatted')
        read (file_unit) ngx, ngy, ngz, nk, nbnd
      end if
      close (file_unit)

200   format('UNK', i5.5, '.', i1)
199   format('UNK', i5.5, '.', 'NC')

      if (allocated(wannier_plot%list)) then
        wann_plot_num = size(wannier_plot%list)
      else
        wann_plot_num = 0
      endif
      allocate (wann_func(-((ngs(1))/2)*ngx:((ngs(1) + 1)/2)*ngx - 1, &
                          -((ngs(2))/2)*ngy:((ngs(2) + 1)/2)*ngy - 1, &
                          -((ngs(3))/2)*ngz:((ngs(3) + 1)/2)*ngz - 1, wann_plot_num), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating wann_func in plot_wannier', comm)
        return
      endif
      wann_func = cmplx_0
      if (spinors) then
        allocate (wann_func_nc(-((ngs(1))/2)*ngx:((ngs(1) + 1)/2)*ngx - 1, &
                               -((ngs(2))/2)*ngy:((ngs(2) + 1)/2)*ngy - 1, &
                               -((ngs(3))/2)*ngz:((ngs(3) + 1)/2)*ngz - 1, 2, wann_plot_num), &
                  stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating wann_func_nc in plot_wannier', comm)
          return
        endif
        wann_func_nc = cmplx_0
      endif
      if (.not. spinors) then
        if (have_disentangled) then
          allocate (r_wvfn_tmp(ngx*ngy*ngz, maxval(dis_manifold%ndimwin)), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error in allocating r_wvfn_tmp in plot_wannier', comm)
            return
          endif
        end if
        allocate (r_wvfn(ngx*ngy*ngz, num_wann), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating r_wvfn in plot_wannier', comm)
          return
        endif
      else
        if (have_disentangled) then
          allocate (r_wvfn_tmp_nc(ngx*ngy*ngz, maxval(dis_manifold%ndimwin), 2), stat=ierr)
          if (ierr /= 0) then
            call set_error_alloc(error, 'Error in allocating r_wvfn_tmp_nc in plot_wannier', comm)
            return
          endif
        end if
        allocate (r_wvfn_nc(ngx*ngy*ngz, num_wann, 2), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating r_wvfn_nc in plot_wannier', comm)
          return
        endif
      endif

      call io_date(cdate, ctime)
      do loop_kpt = 1, num_kpts
        if (dist_k(loop_kpt) /= my_node_id) cycle

        inc_band = .true.
        num_inc = num_wann
        if (have_disentangled) then
          inc_band(:) = dis_manifold%lwindow(:, loop_kpt)
          num_inc = dis_manifold%ndimwin(loop_kpt)
        end if

        if (.not. spinors) then
          write (wfnname, 200) loop_kpt, wvfn_read%spin_channel
        else
          write (wfnname, 199) loop_kpt
        endif
        if (wvfn_read%formatted) then
          open (newunit=file_unit, file=wfnname, form='formatted')
          read (file_unit, *) ix, iy, iz, ik, nbnd
        else
          open (newunit=file_unit, file=wfnname, form='unformatted')
          read (file_unit) ix, iy, iz, ik, nbnd
        end if

        if ((ix /= ngx) .or. (iy /= ngy) .or. (iz /= ngz) .or. (ik /= loop_kpt)) then
          write (stdout, '(1x,a,a)') 'WARNING: mismatch in file', trim(wfnname)
          write (stdout, '(1x,5(a6,I5))') '   ix=', ix, '   iy=', iy, '   iz=', iz, '   ik=', ik, ' nbnd=', nbnd
          write (stdout, '(1x,5(a6,I5))') '  ngx=', ngx, '  ngy=', ngy, '  ngz=', ngz, '  kpt=', loop_kpt, 'bands=', num_bands
          call set_error_file(error, 'plot_wannier', comm)
          return
        end if

        if (have_disentangled) then
          counter = 1
          do loop_b = 1, num_bands
            if (counter > num_inc) exit
            if (wvfn_read%formatted) then
              do nx = 1, ngx*ngy*ngz
                read (file_unit, *) w_real, w_imag
                if (.not. spinors) then
                  r_wvfn_tmp(nx, counter) = cmplx(w_real, w_imag, kind=dp)
                else
                  r_wvfn_tmp_nc(nx, counter, 1) = cmplx(w_real, w_imag, kind=dp) ! up-spinor
                endif
              end do
              if (spinors) then
                do nx = 1, ngx*ngy*ngz
                  read (file_unit, *) w_real, w_imag
                  r_wvfn_tmp_nc(nx, counter, 2) = cmplx(w_real, w_imag, kind=dp) ! down-spinor
                end do
              endif
            else
              if (.not. spinors) then
                read (file_unit) (r_wvfn_tmp(nx, counter), nx=1, ngx*ngy*ngz)
              else
                read (file_unit) (r_wvfn_tmp_nc(nx, counter, 1), nx=1, ngx*ngy*ngz) ! up-spinor
                read (file_unit) (r_wvfn_tmp_nc(nx, counter, 2), nx=1, ngx*ngy*ngz) ! down-spinor
              endif
            end if
            if (inc_band(loop_b)) counter = counter + 1
          end do
        else
          do loop_b = 1, num_bands
            if (wvfn_read%formatted) then
              do nx = 1, ngx*ngy*ngz
                read (file_unit, *) w_real, w_imag
                if (.not. spinors) then
                  r_wvfn(nx, loop_b) = cmplx(w_real, w_imag, kind=dp)
                else
                  r_wvfn_nc(nx, loop_b, 1) = cmplx(w_real, w_imag, kind=dp) ! up-spinor
                endif
              end do
              if (spinors) then
                do nx = 1, ngx*ngy*ngz
                  read (file_unit, *) w_real, w_imag
                  r_wvfn_nc(nx, loop_b, 2) = cmplx(w_real, w_imag, kind=dp) ! down-spinor
                end do
              endif
            else
              if (.not. spinors) then
                read (file_unit) (r_wvfn(nx, loop_b), nx=1, ngx*ngy*ngz)
              else
                read (file_unit) (r_wvfn_nc(nx, loop_b, 1), nx=1, ngx*ngy*ngz) ! up-spinor
                read (file_unit) (r_wvfn_nc(nx, loop_b, 2), nx=1, ngx*ngy*ngz) ! down-spinor
              endif
            end if
          end do
        end if

        close (file_unit)

        if (have_disentangled) then
          if (.not. spinors) then
            r_wvfn = cmplx_0
            do loop_w = 1, num_wann
              do loop_b = 1, num_inc
                r_wvfn(:, loop_w) = r_wvfn(:, loop_w) + &
                                    u_matrix_opt(loop_b, loop_w, loop_kpt)*r_wvfn_tmp(:, loop_b)
              end do
            end do
          else
            r_wvfn_nc = cmplx_0
            do loop_w = 1, num_wann
              do loop_b = 1, num_inc
                call zaxpy(ngx*ngy*ngz, u_matrix_opt(loop_b, loop_w, loop_kpt), r_wvfn_tmp_nc(1, loop_b, 1), 1, & ! up-spinor
                           r_wvfn_nc(1, loop_w, 1), 1)
                call zaxpy(ngx*ngy*ngz, u_matrix_opt(loop_b, loop_w, loop_kpt), r_wvfn_tmp_nc(1, loop_b, 2), 1, & ! down-spinor
                           r_wvfn_nc(1, loop_w, 2), 1)
              end do
            end do
          endif
        end if

        ! nxx, nyy, nzz span a parallelogram in the real space mesh, of side
        ! 2*nphir, and centered around the maximum of phi_i, nphimx(i, 1 2 3)
        !
        ! nx ny nz are the nxx nyy nzz brought back to the unit cell in
        ! which u_nk(r)=cptwrb(r,n)  is represented
        !
        ! There is a big performance improvement in looping over num_wann
        ! in the inner loop. This is poor memory access for wann_func and
        ! but the reduced number of operations wins out.

        do nzz = -((ngs(3))/2)*ngz, ((ngs(3) + 1)/2)*ngz - 1
          nz = mod(nzz, ngz)
          if (nz .lt. 1) nz = nz + ngz
          do nyy = -((ngs(2))/2)*ngy, ((ngs(2) + 1)/2)*ngy - 1
            ny = mod(nyy, ngy)
            if (ny .lt. 1) ny = ny + ngy
            do nxx = -((ngs(1))/2)*ngx, ((ngs(1) + 1)/2)*ngx - 1
              nx = mod(nxx, ngx)
              if (nx .lt. 1) nx = nx + ngx

              scalfac = kpt_latt(1, loop_kpt)*real(nxx - 1, dp)/real(ngx, dp) + &
                        kpt_latt(2, loop_kpt)*real(nyy - 1, dp)/real(ngy, dp) + &
                        kpt_latt(3, loop_kpt)*real(nzz - 1, dp)/real(ngz, dp)
              npoint = nx + (ny - 1)*ngx + (nz - 1)*ngy*ngx
              catmp = exp(twopi*cmplx_i*scalfac)
              do loop_b = 1, num_wann
                do loop_w = 1, wann_plot_num
                  if (.not. spinors) then
                    wann_func(nxx, nyy, nzz, loop_w) = wann_func(nxx, nyy, nzz, loop_w) + &
                                                       u_matrix(loop_b, wannier_plot%list(loop_w), loop_kpt)* &
                                                       r_wvfn(npoint, loop_b)*catmp
                  else
                    wann_func_nc(nxx, nyy, nzz, 1, loop_w) = &
                      wann_func_nc(nxx, nyy, nzz, 1, loop_w) + & ! up-spinor
                      u_matrix(loop_b, wannier_plot%list(loop_w), loop_kpt)*r_wvfn_nc(npoint, loop_b, 1)*catmp
                    wann_func_nc(nxx, nyy, nzz, 2, loop_w) = &
                      wann_func_nc(nxx, nyy, nzz, 2, loop_w) + & ! down-spinor
                      u_matrix(loop_b, wannier_plot%list(loop_w), loop_kpt)*r_wvfn_nc(npoint, loop_b, 2)*catmp
                    if (loop_b == num_wann) then ! last loop
                      upspinor = real(wann_func_nc(nxx, nyy, nzz, 1, loop_w)* &
                                      conjg(wann_func_nc(nxx, nyy, nzz, 1, loop_w)), dp)
                      dnspinor = real(wann_func_nc(nxx, nyy, nzz, 2, loop_w)* &
                                      conjg(wann_func_nc(nxx, nyy, nzz, 2, loop_w)), dp)
                      if (wannier_plot%spinor_phase) then
                        upphase = sign(1.0_dp, real(wann_func_nc(nxx, nyy, nzz, 1, loop_w), dp))
                        dnphase = sign(1.0_dp, real(wann_func_nc(nxx, nyy, nzz, 2, loop_w), dp))
                      else
                        upphase = 1.0_dp; dnphase = 1.0_dp
                      endif
                      select case (wannier_plot%spinor_mode)
                      case ('total')
                        wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(upspinor + dnspinor), 0.0_dp, dp)
                      case ('up')
                        wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(upspinor), 0.0_dp, dp)*upphase
                      case ('down')
                        wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(dnspinor), 0.0_dp, dp)*dnphase
                      case default
                        call set_error_file(error, 'plot_wannier: Invalid wannier_plot_spinor_mode '&
                            &//trim(wannier_plot%spinor_mode), comm)
                        return
                      end select
                      wann_func(nxx, nyy, nzz, loop_w) = &
                        wann_func(nxx, nyy, nzz, loop_w)/real(num_kpts, dp)
                    endif
                  endif
                end do
              end do
            end do
          end do
        end do

      end do !loop over kpoints

      if (spinors) then
        call comms_reduce(wann_func_nc(-((ngs(1))/2)*ngx, -((ngs(2))/2)*ngy, -((ngs(3))/2)*ngz, 1, 1), &
                          size(wann_func_nc), 'SUM', error, comm)
      else
        call comms_reduce(wann_func(-((ngs(1))/2)*ngx, -((ngs(2))/2)*ngy, -((ngs(3))/2)*ngz, 1), &
                          size(wann_func), 'SUM', error, comm)
      endif
      if (allocated(error)) return

      if (on_root) then
        if (spinors) then
          do nzz = -((ngs(3))/2)*ngz, ((ngs(3) + 1)/2)*ngz - 1
            do nyy = -((ngs(2))/2)*ngy, ((ngs(2) + 1)/2)*ngy - 1
              do nxx = -((ngs(1))/2)*ngx, ((ngs(1) + 1)/2)*ngx - 1
                do loop_w = 1, wann_plot_num
                  upspinor = real(wann_func_nc(nxx, nyy, nzz, 1, loop_w)* &
                                  conjg(wann_func_nc(nxx, nyy, nzz, 1, loop_w)), dp)
                  dnspinor = real(wann_func_nc(nxx, nyy, nzz, 2, loop_w)* &
                                  conjg(wann_func_nc(nxx, nyy, nzz, 2, loop_w)), dp)
                  if (wannier_plot%spinor_phase) then
                    upphase = sign(1.0_dp, real(wann_func_nc(nxx, nyy, nzz, 1, loop_w), dp))
                    dnphase = sign(1.0_dp, real(wann_func_nc(nxx, nyy, nzz, 2, loop_w), dp))
                  else
                    upphase = 1.0_dp; dnphase = 1.0_dp
                  endif
                  select case (wannier_plot%spinor_mode)
                  case ('total')
                    wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(upspinor + dnspinor), 0.0_dp, dp)
                  case ('up')
                    wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(upspinor), 0.0_dp, dp)*upphase
                  case ('down')
                    wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(dnspinor), 0.0_dp, dp)*dnphase
                  case default
                    call set_error_file(error, 'plot_wannier: Invalid wannier_plot_spinor_mode ' &
                                        //trim(wannier_plot%spinor_mode), comm)
                    return
                  end select
                  wann_func(nxx, nyy, nzz, loop_w) = wann_func(nxx, nyy, nzz, loop_w)/real(num_kpts, dp)
                end do
              end do
            end do
          end do
        endif
      endif

      if (on_root) then
        if (.not. spinors) then !!!!! For spinor Wannier functions, the steps below are not necessary.
          ! fix the global phase by setting the wannier to
          ! be real at the point where it has max. modulus

          do loop_w = 1, wann_plot_num
            tmaxx = 0.0
            wmod = cmplx_1
            do nzz = -((ngs(3))/2)*ngz, ((ngs(3) + 1)/2)*ngz - 1
              do nyy = -((ngs(2))/2)*ngy, ((ngs(2) + 1)/2)*ngy - 1
                do nxx = -((ngs(1))/2)*ngx, ((ngs(1) + 1)/2)*ngx - 1
                  wann_func(nxx, nyy, nzz, loop_w) = wann_func(nxx, nyy, nzz, loop_w)/real(num_kpts, dp)
                  tmax = real(wann_func(nxx, nyy, nzz, loop_w)* &
                              conjg(wann_func(nxx, nyy, nzz, loop_w)), dp)
                  if (tmax > tmaxx) then
                    tmaxx = tmax
                    wmod = wann_func(nxx, nyy, nzz, loop_w)
                  end if
                end do
              end do
            end do
            wmod = wmod/sqrt(real(wmod)**2 + aimag(wmod)**2)
            wann_func(:, :, :, loop_w) = wann_func(:, :, :, loop_w)/wmod
            write (stdout, '(6x,a,i4,7x,a,f11.6,SP,f11.6,"i")') "Wannier Function Num: ", &
              wannier_plot%list(loop_w), "Phase Factor = ", 1/wmod
          end do
          write (stdout, *) ''
          !
          ! Check the 'reality' of the WF
          !
          do loop_w = 1, wann_plot_num
            ratmax = 0.0_dp
            do nzz = -((ngs(3))/2)*ngz, ((ngs(3) + 1)/2)*ngz - 1
              do nyy = -((ngs(2))/2)*ngy, ((ngs(2) + 1)/2)*ngy - 1
                do nxx = -((ngs(1))/2)*ngx, ((ngs(1) + 1)/2)*ngx - 1
                  if (abs(real(wann_func(nxx, nyy, nzz, loop_w), dp)) >= 0.01_dp) then
                    ratio = abs(aimag(wann_func(nxx, nyy, nzz, loop_w)))/ &
                            abs(real(wann_func(nxx, nyy, nzz, loop_w), dp))
                    ratmax = max(ratmax, ratio)
                  end if
                end do
              end do
            end do
            write (stdout, '(6x,a,i4,7x,a,f11.6)') 'Wannier Function Num: ', wannier_plot%list(loop_w), &
              'Maximum Im/Re Ratio = ', ratmax
          end do
        endif !!!!!
        write (stdout, *) ' '
        if (wannier_plot%format .eq. 'xcrysden') then
          call internal_xsf_format()
        elseif (wannier_plot%format .eq. 'cube') then
          call internal_cube_format(atom_data, wannier_data, real_lattice, bohr, error)
          if (allocated(error)) return
        else
          call set_error_warn(error, 'wannier_plot_format not recognised in wannier_plot', comm)
          return
        endif

        if (print_output%timing_level > 1) call io_stopwatch_stop('plot: wannier', timer)
      end if !on_root

    end associate

    return

  contains

    !================================================!
    subroutine internal_cube_format(atom_data, wannier_data, real_lattice, bohr, error)
      !================================================!
      !
      !! Write WFs in Gaussian cube format.
      !
      !================================================!

      use w90_utility, only: utility_translate_home, utility_cart_to_frac, utility_frac_to_cart, &
        utility_inverse_mat, utility_recip_lattice_base
      use w90_types, only: wannier_data_type, atom_data_type
      use w90_wannier90_types, only: wvfn_read_type
      use w90_error, only: w90_error_type, set_error_alloc, set_error_warn, set_error_dealloc

      implicit none

      type(wannier_data_type), intent(in) :: wannier_data
      type(atom_data_type), intent(in) :: atom_data
      type(w90_error_type), allocatable, intent(out) :: error
      real(kind=dp), intent(in) :: bohr

      real(kind=dp), intent(in) :: real_lattice(3, 3)

      real(kind=dp), allocatable :: wann_cube(:, :, :)
      real(kind=dp) :: inv_lattice(3, 3), recip_lattice(3, 3), pos_frac(3), volume
      real(kind=dp) :: rstart(3), rend(3), rlength(3), orig(3), dgrid(3)
      real(kind=dp) :: moda(3), modb(3)
      real(kind=dp) :: val_Q
      real(kind=dp) :: comf(3), wcf(3), diff(3), difc(3), dist
      integer :: ierr, iname, max_elements
      integer :: isp, iat, nzz, nyy, nxx, loop_w, qxx, qyy, qzz, wann_index
      integer :: istart(3), iend(3), ilength(3)
      integer :: ixx, iyy, izz
      integer :: irdiff(3), icount
      integer, allocatable :: atomic_Z(:)
      logical :: lmol, lcrys
      character(len=2), dimension(109) :: periodic_table = (/ &
           & 'H ', 'He', &
           & 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
           & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
           & 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
           & 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
           & 'Cs', 'Ba', &
           & 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
           & 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
           & 'Fr', 'Ra', &
           & 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', &
           & 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt'/)

      associate (ngs=>wannier_plot%supercell)

        allocate (atomic_Z(atom_data%num_species), stat=ierr)
        if (ierr .ne. 0) then
          call set_error_alloc(error, 'Error: allocating atomic_Z in wannier_plot', comm)
          return
        endif

        call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
        lmol = .false.
        lcrys = .false.
        if (index(wannier_plot%mode, 'mol') > 0) lmol = .true.      ! molecule mode
        if (index(wannier_plot%mode, 'crys') > 0) lcrys = .true.    ! crystal mode

        val_Q = 1.0_dp ! dummy value for cube file

        ! Assign atomic numbers to species
        max_elements = size(periodic_table)
        do isp = 1, atom_data%num_species
          do iname = 1, max_elements
            if (atom_data%symbol(isp) .eq. periodic_table(iname)) then
              atomic_Z(isp) = iname
              exit
            endif
          enddo
        end do

202     format(a, '_', i5.5, '.cube')

        ! Lengths of real and reciprocal lattice vectors
        do i = 1, 3
          moda(i) = sqrt(real_lattice(i, 1)*real_lattice(i, 1) &
                         + real_lattice(i, 2)*real_lattice(i, 2) &
                         + real_lattice(i, 3)*real_lattice(i, 3))
          modb(i) = sqrt(recip_lattice(i, 1)*recip_lattice(i, 1) &
                         + recip_lattice(i, 2)*recip_lattice(i, 2) &
                         + recip_lattice(i, 3)*recip_lattice(i, 3))
        enddo

        ! Grid spacing in each lattice direction
        dgrid(1) = moda(1)/ngx; dgrid(2) = moda(2)/ngy; dgrid(3) = moda(3)/ngz

        ! Find "centre of mass" of atomic positions (in fractional coordinates)
        call utility_inverse_mat(real_lattice, inv_lattice)
        comf(:) = 0.0_dp
        do isp = 1, atom_data%num_species
          do iat = 1, atom_data%species_num(isp)
            call utility_cart_to_frac(atom_data%pos_cart(:, iat, isp), pos_frac, inv_lattice)
            comf(:) = comf(:) + pos_frac(:)
          enddo
        enddo
        comf(:) = comf(:)/atom_data%num_atoms

        ! Loop over WFs
        do loop_w = 1, wann_plot_num

          wann_index = wannier_plot%list(loop_w)
          write (wancube, 202) trim(seedname), wann_index

          ! Find start and end of cube wrt simulation (home) cell origin
          do i = 1, 3
            ! ... in terms of distance along each lattice vector direction i
            rstart(i) = (wannier_data%centres(1, wann_index)*recip_lattice(i, 1) &
                         + wannier_data%centres(2, wann_index)*recip_lattice(i, 2) &
                         + wannier_data%centres(3, wann_index)*recip_lattice(i, 3))*moda(i)/twopi &
                        - twopi*wannier_plot%radius/(moda(i)*modb(i))
            rend(i) = (wannier_data%centres(1, wann_index)*recip_lattice(i, 1) &
                       + wannier_data%centres(2, wann_index)*recip_lattice(i, 2) &
                       + wannier_data%centres(3, wann_index)*recip_lattice(i, 3))*moda(i)/twopi &
                      + twopi*wannier_plot%radius/(moda(i)*modb(i))
          enddo

          rlength(:) = rend(:) - rstart(:)
          ilength(:) = ceiling(rlength(:)/dgrid(:))

          ! ... in terms of integer gridpoints along each lattice vector direction i
          istart(:) = floor(rstart(:)/dgrid(:)) + 1
          iend(:) = istart(:) + ilength(:) - 1

          ! Origin of cube wrt simulation (home) cell in Cartesian co-ordinates
          do i = 1, 3
            orig(i) = real(istart(1) - 1, dp)*dgrid(1)*real_lattice(1, i)/moda(1) &
                      + real(istart(2) - 1, dp)*dgrid(2)*real_lattice(2, i)/moda(2) &
                      + real(istart(3) - 1, dp)*dgrid(3)*real_lattice(3, i)/moda(3)
          enddo

          ! Debugging
          if (print_output%iprint > 3) then
            write (stdout, '(a,i12)') 'loop_w  =', loop_w
            write (stdout, '(a,3f12.6)') 'comf    =', (comf(i), i=1, 3)
            write (stdout, '(a,3i12)') 'ngi     =', ngx, ngy, ngz
            write (stdout, '(a,3f12.6)') 'dgrid   =', (dgrid(i), i=1, 3)
            write (stdout, '(a,3f12.6)') 'rstart  =', (rstart(i), i=1, 3)
            write (stdout, '(a,3f12.6)') 'rend    =', (rend(i), i=1, 3)
            write (stdout, '(a,3f12.6)') 'rlength =', (rlength(i), i=1, 3)
            write (stdout, '(a,3i12)') 'istart  =', (istart(i), i=1, 3)
            write (stdout, '(a,3i12)') 'iend    =', (iend(i), i=1, 3)
            write (stdout, '(a,3i12)') 'ilength =', (ilength(i), i=1, 3)
            write (stdout, '(a,3f12.6)') 'orig    =', (orig(i), i=1, 3)
            write (stdout, '(a,3f12.6)') 'wann_cen=', (wannier_data%centres(i, wann_index), i=1, 3)
          endif

          allocate (wann_cube(1:ilength(1), 1:ilength(2), 1:ilength(3)), stat=ierr)
          if (ierr .ne. 0) then
            call set_error_alloc(error, 'Error: allocating wann_cube in wannier_plot', comm)
            return
          endif

          ! initialise
          wann_cube = 0.0_dp

          do nzz = 1, ilength(3)
            qzz = nzz + istart(3) - 1
            izz = int((abs(qzz) - 1)/ngz)
!            if (qzz.lt.-ngz) qzz=qzz+izz*ngz
!            if (qzz.gt.(ngs(3)-1)*ngz-1) then
            if (qzz .lt. (-((ngs(3))/2)*ngz)) qzz = qzz + izz*ngz
            if (qzz .gt. ((ngs(3) + 1)/2)*ngz - 1) then
              write (stdout, *) 'Error plotting WF cube. Try one of the following:'
              write (stdout, *) '   (1) increase wannier_plot_supercell;'
              write (stdout, *) '   (2) decrease wannier_plot_radius;'
              write (stdout, *) '   (3) set wannier_plot_format=xcrysden'
              call set_error_warn(error, 'Error plotting WF cube.', comm)
              return
            endif
            do nyy = 1, ilength(2)
              qyy = nyy + istart(2) - 1
              iyy = int((abs(qyy) - 1)/ngy)
!               if (qyy.lt.-ngy) qyy=qyy+iyy*ngy
!               if (qyy.gt.(ngs(2)-1)*ngy-1) then
              if (qyy .lt. (-((ngs(2))/2)*ngy)) qyy = qyy + iyy*ngy
              if (qyy .gt. ((ngs(2) + 1)/2)*ngy - 1) then
                write (stdout, *) 'Error plotting WF cube. Try one of the following:'
                write (stdout, *) '   (1) increase wannier_plot_supercell;'
                write (stdout, *) '   (2) decrease wannier_plot_radius;'
                write (stdout, *) '   (3) set wannier_plot_format=xcrysden'
                call set_error_warn(error, 'Error plotting WF cube.', comm)
                return
              endif
              do nxx = 1, ilength(1)
                qxx = nxx + istart(1) - 1
                ixx = int((abs(qxx) - 1)/ngx)
!                  if (qxx.lt.-ngx) qxx=qxx+ixx*ngx
!                  if (qxx.gt.(ngs(1)-1)*ngx-1) then
                if (qxx .lt. (-((ngs(1))/2)*ngx)) qxx = qxx + ixx*ngx
                if (qxx .gt. ((ngs(1) + 1)/2)*ngx - 1) then
                  write (stdout, *) 'Error plotting WF cube. Try one of the following:'
                  write (stdout, *) '   (1) increase wannier_plot_supercell;'
                  write (stdout, *) '   (2) decrease wannier_plot_radius;'
                  write (stdout, *) '   (3) set wannier_plot_format=xcrysden'
                  call set_error_warn(error, 'Error plotting WF cube.', comm)
                  return
                endif
                wann_cube(nxx, nyy, nzz) = real(wann_func(qxx, qyy, qzz, loop_w), dp)
              enddo
            enddo
          enddo

          ! WF centre in fractional coordinates
          call utility_cart_to_frac(wannier_data%centres(:, wann_index), wcf(:), inv_lattice)

          ! The vector (in fractional coordinates) from WF centre to "centre of mass"
          diff(:) = comf(:) - wcf(:)

          ! Corresponding nearest cell vector
          irdiff(:) = nint(diff(:))

          if (print_output%iprint > 3) then
            write (stdout, '(a,3f12.6)') 'wcf     =', (wcf(i), i=1, 3)
            write (stdout, '(a,3f12.6)') 'diff    =', (diff(i), i=1, 3)
            write (stdout, '(a,3i12)') 'irdiff  =', (irdiff(i), i=1, 3)
          endif

          if (lmol) then ! In "molecule mode" translate origin of cube to bring it in coincidence with the atomic positions
            orig(:) = orig(:) + real(irdiff(1), kind=dp)*real_lattice(1, :) &
                      + real(irdiff(2), kind=dp)*real_lattice(2, :) &
                      + real(irdiff(3), kind=dp)*real_lattice(3, :)
            if (print_output%iprint > 3) write (stdout, '(a,3f12.6,/)') 'orig-new=', (orig(i), i=1, 3)
          else ! In "crystal mode" count number of atoms within a given radius of wannier centre
            icount = 0
            do isp = 1, atom_data%num_species
              do iat = 1, atom_data%species_num(isp)
                call utility_cart_to_frac(atom_data%pos_cart(:, iat, isp), pos_frac, inv_lattice)
                do nzz = -ngs(3)/2, (ngs(3) + 1)/2
                  do nyy = -ngs(2)/2, (ngs(2) + 1)/2
                    do nxx = -ngs(1)/2, (ngs(1) + 1)/2
                      diff(:) = pos_frac(:) - wcf(:) &
                                + (/real(nxx, kind=dp), real(nyy, kind=dp), real(nzz, kind=dp)/)
                      call utility_frac_to_cart(diff, difc, real_lattice)
                      dist = sqrt(difc(1)*difc(1) + difc(2)*difc(2) + difc(3)*difc(3))
                      if (dist .le. (wannier_plot%scale*wannier_plot%radius)) then
                        icount = icount + 1
                      endif
                    enddo
                  enddo
                enddo
              enddo ! iat
            enddo ! isp
            if (print_output%iprint > 3) write (stdout, '(a,i12)') 'icount  =', icount
          endif

          ! Write cube file (everything in Bohr)
          open (newunit=file_unit, file=trim(wancube), form='formatted', status='unknown')
          ! First two lines are comments
          write (file_unit, *) '     Generated by Wannier90 code http://www.wannier.org'
          write (file_unit, *) '     On ', cdate, ' at ', ctime
          ! Number of atoms, origin of cube (Cartesians) wrt simulation (home) cell
          if (lmol) then
            write (file_unit, '(i4,3f13.5)') atom_data%num_atoms, orig(1)/bohr, orig(2)/bohr, orig(3)/bohr
          else
            write (file_unit, '(i4,3f13.5)') icount, orig(1)/bohr, orig(2)/bohr, orig(3)/bohr
          endif
          ! Number of grid points in each direction, lattice vector
          write (file_unit, '(i4,3f13.5)') ilength(1), real_lattice(1, 1)/(real(ngx, dp)*bohr), &
            real_lattice(1, 2)/(real(ngx, dp)*bohr), real_lattice(1, 3)/(real(ngx, dp)*bohr)
          write (file_unit, '(i4,3f13.5)') ilength(2), real_lattice(2, 1)/(real(ngy, dp)*bohr), &
            real_lattice(2, 2)/(real(ngy, dp)*bohr), real_lattice(2, 3)/(real(ngy, dp)*bohr)
          write (file_unit, '(i4,3f13.5)') ilength(3), real_lattice(3, 1)/(real(ngz, dp)*bohr), &
            real_lattice(3, 2)/(real(ngz, dp)*bohr), real_lattice(3, 3)/(real(ngz, dp)*bohr)

          ! Atomic number, valence charge, position of atom
!         do isp=1,num_species
!            do iat=1,atoms_species_num(isp)
!               write(file_unit,'(i4,4f13.5)') atomic_Z(isp), val_Q, (atoms_pos_cart(i,iat,isp)/bohr,i=1,3)
!            end do
!         end do

          do isp = 1, atom_data%num_species
            do iat = 1, atom_data%species_num(isp)
              if (lmol) then ! In "molecule mode", write atomic coordinates as they appear in input file
                write (file_unit, '(i4,4f13.5)') atomic_Z(isp), val_Q, (atom_data%pos_cart(i, iat, isp)/bohr, i=1, 3)
              else           ! In "crystal mode", write atoms in supercell within a given radius of Wannier centre
                call utility_cart_to_frac(atom_data%pos_cart(:, iat, isp), pos_frac, inv_lattice)
                do nzz = -ngs(3)/2, (ngs(3) + 1)/2
                  do nyy = -ngs(2)/2, (ngs(2) + 1)/2
                    do nxx = -ngs(1)/2, (ngs(1) + 1)/2
                      diff(:) = pos_frac(:) - wcf(:) &
                                + (/real(nxx, kind=dp), real(nyy, kind=dp), real(nzz, kind=dp)/)
                      call utility_frac_to_cart(diff, difc, real_lattice)
                      dist = sqrt(difc(1)*difc(1) + difc(2)*difc(2) + difc(3)*difc(3))
                      if (dist .le. (wannier_plot%scale*wannier_plot%radius)) then
                        diff(:) = pos_frac(:) &
                                  + (/real(nxx, kind=dp), real(nyy, kind=dp), real(nzz, kind=dp)/)
                        call utility_frac_to_cart(diff, difc, real_lattice)
                        write (file_unit, '(i4,4f13.5)') atomic_Z(isp), val_Q, (difc(i)/bohr, i=1, 3)
                      endif
                    enddo
                  enddo
                enddo
              endif
            enddo ! iat
          enddo ! isp

          ! Volumetric data in batches of 6 values per line, 'z'-direction first.
          do nxx = 1, ilength(1)
            do nyy = 1, ilength(2)
              do nzz = 1, ilength(3)
                write (file_unit, '(E13.5)', advance='no') wann_cube(nxx, nyy, nzz)
                if ((mod(nzz, 6) .eq. 0) .or. (nzz .eq. ilength(3))) write (file_unit, '(a)') ''
              enddo
            enddo
          enddo

          deallocate (wann_cube, stat=ierr)
          if (ierr .ne. 0) then
            call set_error_dealloc(error, 'Error: deallocating wann_cube in wannier_plot', comm)
            return
          endif

        end do

        deallocate (atomic_Z, stat=ierr)
        if (ierr .ne. 0) then
          call set_error_dealloc(error, 'Error: deallocating atomic_Z in wannier_plot', comm)
          return
        endif

      end associate

      return

    end subroutine internal_cube_format

    subroutine internal_xsf_format()

      implicit none

201   format(a, '_', i5.5, '.xsf')

      associate (ngs=>wannier_plot%supercell)

        ! this is to create the WF...xsf output, to be read by XCrySDen
        ! (coordinates + isosurfaces)

        x_0ang = -real(((ngs(1))/2)*ngx + 1, dp)/real(ngx, dp)*real_lattice(1, 1) - &
                 real(((ngs(2))/2)*ngy + 1, dp)/real(ngy, dp)*real_lattice(2, 1) - &
                 real(((ngs(3))/2)*ngz + 1, dp)/real(ngz, dp)*real_lattice(3, 1)
        y_0ang = -real(((ngs(1))/2)*ngx + 1, dp)/real(ngx, dp)*real_lattice(1, 2) - &
                 real(((ngs(2))/2)*ngy + 1, dp)/real(ngy, dp)*real_lattice(2, 2) - &
                 real(((ngs(3))/2)*ngz + 1, dp)/real(ngz, dp)*real_lattice(3, 2)
        z_0ang = -real(((ngs(1))/2)*ngx + 1, dp)/real(ngx, dp)*real_lattice(1, 3) - &
                 real(((ngs(2))/2)*ngy + 1, dp)/real(ngy, dp)*real_lattice(2, 3) - &
                 real(((ngs(3))/2)*ngz + 1, dp)/real(ngz, dp)*real_lattice(3, 3)

        fxcry(1) = real(ngs(1)*ngx - 1, dp)/real(ngx, dp)
        fxcry(2) = real(ngs(2)*ngy - 1, dp)/real(ngy, dp)
        fxcry(3) = real(ngs(3)*ngz - 1, dp)/real(ngz, dp)
        do j = 1, 3
          dirl(:, j) = fxcry(:)*real_lattice(:, j)
        end do

        do loop_b = 1, wann_plot_num

          write (wanxsf, 201) trim(seedname), wannier_plot%list(loop_b)

          open (newunit=file_unit, file=trim(wanxsf), form='formatted', status='unknown')
          write (file_unit, *) '      #'
          write (file_unit, *) '      # Generated by the Wannier90 code http://www.wannier.org'
          write (file_unit, *) '      # On ', cdate, ' at ', ctime
          write (file_unit, *) '      #'
          ! should pass this into the code
          if (index(wannier_plot%mode, 'mol') > 0) then
            write (file_unit, '("ATOMS")')
          else
            write (file_unit, '("CRYSTAL")')
            write (file_unit, '("PRIMVEC")')
            write (file_unit, '(3f12.7)') real_lattice(1, 1), real_lattice(1, 2), real_lattice(1, 3)
            write (file_unit, '(3f12.7)') real_lattice(2, 1), real_lattice(2, 2), real_lattice(2, 3)
            write (file_unit, '(3f12.7)') real_lattice(3, 1), real_lattice(3, 2), real_lattice(3, 3)
            write (file_unit, '("CONVVEC")')
            write (file_unit, '(3f12.7)') real_lattice(1, 1), real_lattice(1, 2), real_lattice(1, 3)
            write (file_unit, '(3f12.7)') real_lattice(2, 1), real_lattice(2, 2), real_lattice(2, 3)
            write (file_unit, '(3f12.7)') real_lattice(3, 1), real_lattice(3, 2), real_lattice(3, 3)
            write (file_unit, '("PRIMCOORD")')
            write (file_unit, '(i6,"  1")') atom_data%num_atoms
          endif
          do nsp = 1, atom_data%num_species
            do nat = 1, atom_data%species_num(nsp)
              write (file_unit, '(a2,3x,3f12.7)') atom_data%symbol(nsp), (atom_data%pos_cart(i, nat, nsp), i=1, 3)
            end do
          end do

          write (file_unit, '(/)')
          write (file_unit, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
          write (file_unit, '(3i6)') ngs(1)*ngx, ngs(2)*ngy, ngs(3)*ngz
          write (file_unit, '(3f12.6)') x_0ang, y_0ang, z_0ang
          write (file_unit, '(3f12.7)') dirl(1, 1), dirl(1, 2), dirl(1, 3)
          write (file_unit, '(3f12.7)') dirl(2, 1), dirl(2, 2), dirl(2, 3)
          write (file_unit, '(3f12.7)') dirl(3, 1), dirl(3, 2), dirl(3, 3)
          write (file_unit, '(6e13.5)') &
            (((real(wann_func(nx, ny, nz, loop_b)), nx=-((ngs(1))/2)*ngx, ((ngs(1) + 1)/2)*ngx - 1), &
              ny=-((ngs(2))/2)*ngy, ((ngs(2) + 1)/2)*ngy - 1), nz=-((ngs(3))/2)*ngz, ((ngs(3) + 1)/2)*ngz - 1)
          write (file_unit, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
          close (file_unit)

        end do

      end associate
      return
    end subroutine internal_xsf_format

  end subroutine plot_wannier

  !================================================!
  subroutine plot_u_matrices(u_matrix_opt, u_matrix, kpt_latt, dis_manifold, &
                             have_disentangled, num_wann, num_kpts, num_bands, seedname)
    !================================================!
    !
    !! Plot u_matrix and u_matrix_opt to textfiles in readable format
    !
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_time, io_date
    use w90_types, only: dis_manifold_type

    implicit none

    character(len=50), intent(in)  :: seedname
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    logical, intent(in) :: have_disentangled
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(dis_manifold_type), intent(in) :: dis_manifold

    character(len=33)  :: header
    character(len=9)   :: cdate, ctime
    complex(kind=dp), allocatable :: utmp(:, :) ! re-indexed u_matrix_opt for printout
    integer :: matunit, i, j, nkp, ioff, nbw

    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime

    open (newunit=matunit, file=trim(seedname)//'_u.mat', form='formatted')

    write (matunit, *) header
    write (matunit, *) num_kpts, num_wann, num_wann

    do nkp = 1, num_kpts
      write (matunit, *)
      write (matunit, '(f15.10,sp,f15.10,sp,f15.10)') kpt_latt(:, nkp)
      write (matunit, '(f15.10,sp,f15.10)') ((u_matrix(i, j, nkp), i=1, num_wann), j=1, num_wann)
    end do
    close (matunit)

    if (have_disentangled) then
      allocate (utmp(num_bands, num_wann))

      open (newunit=matunit, file=trim(seedname)//'_u_dis.mat', form='formatted')
      write (matunit, *) header
      write (matunit, *) num_kpts, num_wann, num_bands
      do nkp = 1, num_kpts
        utmp = 0.d0
        ioff = dis_manifold%nfirstwin(nkp)
        nbw = dis_manifold%ndimwin(nkp)
        utmp(ioff:ioff + nbw - 1, :) = u_matrix_opt(1:nbw, :, nkp)
        write (matunit, *)
        write (matunit, '(f15.10,sp,f15.10,sp,f15.10)') kpt_latt(:, nkp)
        !write (matunit, '(f15.10,sp,f15.10)') ((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann)
        write (matunit, '(f15.10,sp,f15.10)') ((utmp(i, j), i=1, num_bands), j=1, num_wann)
      end do
      close (matunit)
      deallocate (utmp)
    endif

  end subroutine plot_u_matrices

  !================================================!
  subroutine plot_bvec(kmesh_info, num_kpts, seedname, error, comm)
    !================================================!
    !! June 2018: RM and SP
    !! Write to file the matrix elements of bvector and their weights
    !! This is used by EPW to compute the velocity.
    !! You need "write_bvec = .true." in your wannier input
    !!
    !================================================!

    use w90_io, only: io_date
    use w90_constants, only: dp
    use w90_types, only: kmesh_info_type
    use w90_error, only: w90_error_type, set_error_file

    implicit none

    type(kmesh_info_type), intent(in) :: kmesh_info
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer :: nkp, nn, file_unit
    character(len=33) :: header
    character(len=9)  :: cdate, ctime

    integer, intent(in) :: num_kpts
    character(len=50), intent(in)  :: seedname

    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime

    open (newunit=file_unit, file=trim(seedname)//'.bvec', form='formatted', status='unknown', err=101)
    write (file_unit, *) header ! Date and time
    write (file_unit, *) num_kpts, kmesh_info%nntot
    do nkp = 1, num_kpts
      do nn = 1, kmesh_info%nntot
        write (file_unit, '(4F20.14)') kmesh_info%bk(:, nn, nkp), kmesh_info%wb(nn)
      enddo
    enddo
    close (file_unit)

    return

101 call set_error_file(error, 'Error: plot_bvec: problem opening file '//trim(seedname)//'.bvec', comm)
    return

  end subroutine plot_bvec

  !================================================!
  subroutine plot_write_rmn(kmesh_info, m_matrix, kpt_latt, irvec, nrpts, num_kpts, &
                            num_wann, seedname, dist_k, error, comm)
    !================================================!
    !
    !! Write out the matrix elements of r
    !
    !================================================!

    use w90_comms, only: comms_reduce, w90_comm_type, mpisize, mpirank
    use w90_constants, only: twopi, cmplx_i, dp
    use w90_error, only: w90_error_type, set_error_file
    use w90_io, only: io_date
    use w90_types, only: kmesh_info_type

    implicit none

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: nrpts
    integer, intent(in) :: irvec(:, :)
    integer, intent(in) :: num_wann
    integer, intent(in) :: num_kpts
    integer, intent(in) :: dist_k(:) ! MPI k-point distribution
    real(kind=dp), intent(in)     :: kpt_latt(:, :)
    complex(kind=dp), intent(in)  :: m_matrix(:, :, :, :)
    character(len=50), intent(in) :: seedname

    ! local variables
    integer :: loop_rpt, m, n, nkp, ind, nn, file_unit, ierr
    integer :: my_node_id, nkp_rank
    ! nkp_rank is the rank-local kpoint index for m_matrix decomposition
    real(kind=dp) :: rdotk
    complex(kind=dp) :: fac
    complex(kind=dp) :: position(3)
    character(len=33) :: header
    character(len=9)  :: cdate, ctime
    logical :: on_root = .false.

    my_node_id = mpirank(comm)

    if (my_node_id == 0) on_root = .true.

    if (on_root) then
      open (newunit=file_unit, file=trim(seedname)//'_r.dat', form='formatted', status='unknown', &
            iostat=ierr)
      if (ierr /= 0) then
        call set_error_file(error, 'Error: hamiltonian_write_rmn: problem opening file '//trim(seedname)//'_r', comm)
        return
      endif

      call io_date(cdate, ctime)
      header = 'written on '//cdate//' at '//ctime
      write (file_unit, *) header ! Date and time
      write (file_unit, *) num_wann
      write (file_unit, *) nrpts
    endif

    do loop_rpt = 1, nrpts
      do m = 1, num_wann
        do n = 1, num_wann

          position(:) = 0._dp
          nkp_rank = 1
          do nkp = 1, num_kpts
            if (dist_k(nkp) /= my_node_id) cycle

            rdotk = twopi*dot_product(kpt_latt(:, nkp), real(irvec(:, loop_rpt), dp))
            fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
            do ind = 1, 3
              do nn = 1, kmesh_info%nntot
                if (m .eq. n) then
                  ! For loop_rpt==rpt_origin, this reduces to
                  ! Eq.(32) of Marzari and Vanderbilt PRB 56,
                  ! 12847 (1997). Otherwise, is is Eq.(44)
                  ! Wang, Yates, Souza and Vanderbilt PRB 74,
                  ! 195118 (2006), modified according to
                  ! Eqs.(27,29) of Marzari and Vanderbilt
                  position(ind) = position(ind) - kmesh_info%wb(nn)*kmesh_info%bk(ind, nn, nkp) &
                                  *aimag(log(m_matrix(n, m, nn, nkp_rank)))*fac
                else
                  ! Eq.(44) Wang, Yates, Souza and Vanderbilt PRB 74, 195118 (2006)
                  position(ind) = position(ind) + cmplx_i*kmesh_info%wb(nn) &
                                  *kmesh_info%bk(ind, nn, nkp)*m_matrix(n, m, nn, nkp_rank)*fac
                endif
              end do
            end do
            nkp_rank = nkp_rank + 1
          end do ! global k list
          call comms_reduce(position(1), 3, 'SUM', error, comm)
          if (on_root) write (file_unit, '(5I5,6F12.6)') irvec(:, loop_rpt), n, m, position(:)
        end do
      end do
    end do

    if (on_root) close (file_unit)
  end subroutine plot_write_rmn

  !================================================!
  subroutine plot_write_r2mn(num_kpts, num_wann, kmesh_info, m_matrix, seedname, dist_k, error, comm)
    !================================================!
    !
    ! Write seedname.r2mn file
    !
    !================================================!

    use w90_comms, only: w90_comm_type
    use w90_constants, only: dp
    use w90_error, only: w90_error_type, set_error_file
    use w90_types, only: kmesh_info_type

    implicit none

    type(kmesh_info_type), intent(in) :: kmesh_info
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_kpts, num_wann
    integer, intent(in) :: dist_k(:) ! MPI k-point distribution
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    character(len=50), intent(in)  :: seedname

    integer :: r2mnunit, nw1, nw2, nkp, nn, ierr
    integer :: nkp_rank, my_node_id
    real(kind=dp) :: r2ave_mn, delta
    logical :: on_root = .false.

    ! note that here I use formulas analogue to Eq. 23, and not to the
    ! shift-invariant Eq. 32 .

    my_node_id = mpirank(comm)

    if (my_node_id == 0) on_root = .true.

    if (on_root) then
      open (newunit=r2mnunit, file=trim(seedname)//'.r2mn', form='formatted', iostat=ierr)
      if (ierr /= 0) then
        call set_error_file(error, 'Error opening file '//trim(seedname)//'.r2mn in plot_write_r2mn', comm)
        return
      endif
    endif

    do nw1 = 1, num_wann
      do nw2 = 1, num_wann
        r2ave_mn = 0.0_dp
        delta = 0.0_dp
        if (nw1 .eq. nw2) delta = 1.0_dp
        nkp_rank = 1
        do nkp = 1, num_kpts
          if (dist_k(nkp) /= my_node_id) cycle

          do nn = 1, kmesh_info%nntot
            ! [GP-begin, Apr13, 2012: corrected sign inside "real"]
            r2ave_mn = r2ave_mn + kmesh_info%wb(nn)* &
                       (2.0_dp*delta - real(m_matrix(nw1, nw2, nn, nkp_rank) + &
                                            conjg(m_matrix(nw2, nw1, nn, nkp_rank)), dp))
          enddo
          nkp_rank = nkp_rank + 1
        enddo ! global k list
        call comms_reduce(r2ave_mn, 1, 'SUM', error, comm)
        r2ave_mn = r2ave_mn/real(num_kpts, dp)
        if (on_root) write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn
      enddo
    enddo
    if (on_root) close (r2mnunit)
  end subroutine plot_write_r2mn

  !================================================!
  subroutine plot_calc_projection(num_bands, num_wann, num_kpts, u_matrix_opt, eigval, lwindow, &
                                  timing_level, iprint, stdout, timer)
    !================================================!
    !
    ! Calculates and writes the projection of each Wannier function
    ! on the original bands within the outer window.
    !
    !================================================!

    use w90_comms, only: w90_comm_type
    use w90_constants, only: dp
    use w90_error, only: w90_error_type
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level, iprint
    logical, intent(in) :: lwindow(:, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    real(kind=dp), intent(in) :: eigval(:, :)

    ! local variables
    integer :: nw, nb, nkp, counter
    real(kind=dp) :: summ

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch_start('wann: calc_projection', timer)

    if (iprint > 0) then
      write (stdout, '(/1x,a78)') repeat('-', 78)
      write (stdout, '(1x,9x,a)') &
        'Projection of Bands in Outer Window on all Wannier Functions'
      write (stdout, '(1x,8x,62a)') repeat('-', 62)
      write (stdout, '(1x,16x,a)') '   Kpt  Band      Eigval      |Projection|^2'
      write (stdout, '(1x,16x,a47)') repeat('-', 47)
    endif

    do nkp = 1, num_kpts
      counter = 0
      do nb = 1, num_bands
        if (lwindow(nb, nkp)) then
          counter = counter + 1
          summ = 0.0_dp
          do nw = 1, num_wann
            summ = summ + abs(u_matrix_opt(counter, nw, nkp))**2
          enddo
          if (iprint > 0) write (stdout, '(1x,16x,i5,1x,i5,1x,f14.6,2x,f14.8)') &
            nkp, nb, eigval(nb, nkp), summ
        endif
      enddo
    enddo
    if (iprint > 0) write (stdout, '(1x,a78/)') repeat('-', 78)

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch_stop('wann: calc_projection', timer)

    return

  end subroutine plot_calc_projection

  !================================================!
  subroutine plot_write_vdw_data(num_wann, wannier_data, real_lattice, u_matrix, u_matrix_opt, &
                                 have_disentangled, w90_system, error, comm, stdout, seedname)
    !================================================!
    !
    ! Write a file with Wannier centres, spreads and occupations for
    ! post-processing computation of vdW C6 coeffients.
    !
    ! Based on code written by Lampros Andrinopoulos.
    !================================================!

    use w90_constants, only: cmplx_0
    use w90_constants, only: dp
    use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_input
    use w90_io, only: io_date
    use w90_types, only: wannier_data_type, w90_system_type
    use w90_utility, only: utility_translate_home

    implicit none

    type(wannier_data_type), intent(in) :: wannier_data
    type(w90_system_type), intent(in) :: w90_system
    type(w90_error_type), allocatable, intent(out)  :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    integer, intent(in) :: stdout
    logical, intent(in) :: have_disentangled
    character(len=50), intent(in)  :: seedname

    integer          :: iw, vdw_unit, r, s, k, m, ierr, ndim
    real(kind=dp)    :: wc(3, num_wann)
    real(kind=dp)    :: ws(num_wann)
    complex(kind=dp), allocatable :: f_w(:, :), v_matrix(:, :) !f_w2(:,:)

    wc = wannier_data%centres
    ws = wannier_data%spreads

    ! translate Wannier centres to the home unit cell
    do iw = 1, num_wann
      call utility_translate_home(wc(:, iw), real_lattice)
    enddo

    allocate (f_w(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating f_w in plot_write_vdw_data', comm)
      return
    endif

!~    ! aam: remove f_w2 at end
!~    allocate(f_w2(num_wann, num_wann),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating f_w2 in plot_write_vdw_data')

    if (have_disentangled) then

      ! dimension of occupied subspace
      if (w90_system%num_valence_bands .le. 0) then
        call set_error_input(error, 'Please set num_valence_bands in seedname.win', comm)
        return
      endif

      ndim = w90_system%num_valence_bands

      allocate (v_matrix(ndim, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating V_matrix in plot_write_vdw_data', comm)
        return
      endif

      ! aam: initialise
      f_w(:, :) = cmplx_0
      v_matrix(:, :) = cmplx_0
!~       f_w2(:,:) = cmplx_0

      ! aam: IN THE END ONLY NEED DIAGONAL PART, SO COULD SIMPLIFY...
      ! aam: calculate V = U_opt . U
      do s = 1, num_wann
        do k = 1, ndim
          do m = 1, num_wann
            v_matrix(k, s) = v_matrix(k, s) + u_matrix_opt(k, m, 1)*u_matrix(m, s, 1)
          enddo
        enddo
      enddo

      ! aam: calculate f = V^dagger . V
      do r = 1, num_wann
        do s = 1, num_wann
          do k = 1, ndim
            f_w(r, s) = f_w(r, s) + v_matrix(k, s)*conjg(v_matrix(k, r))
          enddo
        enddo
      enddo

!~       ! original formulation
!~       do r=1,num_wann
!~          do s=1,num_wann
!~             do nkp=1,num_kpts
!~                do k=1,ndimfroz(nkp)
!~                   do m=1,num_wann
!~                      do l=1,num_wann
!~                         f_w2(r,s) = f_w2(r,s) + &
!~                              u_matrix_opt(k,m,nkp) * u_matrix(m,s,nkp) * &
!~                              conjg(u_matrix_opt(k,l,nkp)) * conjg(u_matrix(l,r,nkp))
!~                      end do
!~                   end do
!~                end do
!~             end do
!~          end do
!~       end do

!~       ! test equivalence
!~       do r=1,num_wann
!~          do s=1,num_wann
!~             if (abs(real(f_w(r,s),dp)-real(f_w2(r,s),dp)).gt.eps6) then
!~                write(*,'(i6,i6,f16.10,f16.10)') r,s,real(f_w(r,s),dp),real(f_w2(r,s),dp)
!~             endif
!~             if (abs(aimag(f_w(r,s))-aimag(f_w2(r,s))).gt.eps6) then
!~                write(*,'(a,i6,i6,f16.10,f16.10)') 'Im: ',r,s,aimag(f_w(r,s)),aimag(f_w2(r,s))
!~             endif
!~          enddo
!~       enddo

    else
      ! for valence only, all occupancies are unity
      f_w(:, :) = 1.0_dp
    endif

    ! aam: write the seedname.vdw file directly here
    open (newunit=vdw_unit, file=trim(seedname)//'.vdw', action='write')
    if (have_disentangled) then
      write (vdw_unit, '(a)') 'disentangle T'
    else
      write (vdw_unit, '(a)') 'disentangle F'
    endif
    write (vdw_unit, '(a)') 'amalgamate F'
    write (vdw_unit, '(a,i3)') 'degeneracy', w90_system%num_elec_per_state
    write (vdw_unit, '(a)') 'num_frag 2'
    write (vdw_unit, '(a)') 'num_wann'
    write (vdw_unit, '(i3,1x,i3)') num_wann/2, num_wann/2
    write (vdw_unit, '(a)') 'tol_occ 0.9'
    write (vdw_unit, '(a)') 'pxyz'
    write (vdw_unit, '(a)') 'F F F'
    write (vdw_unit, '(a)') 'F F F'
    write (vdw_unit, '(a)') 'tol_dist 0.05'
    write (vdw_unit, '(a)') 'centres_spreads_occ'
    write (vdw_unit, '(a)') 'ang'
    do iw = 1, num_wann
      write (vdw_unit, '(4(f13.10,1x),1x,f11.8)') wc(1:3, iw), ws(iw), real(f_w(iw, iw))
    end do
    close (vdw_unit)

    write (stdout, '(/a/)') ' vdW data written to file '//trim(seedname)//'.vdw'

    if (have_disentangled) then
      deallocate (v_matrix, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating v_matrix in plot_write_vdw_data', comm)
        return
      endif
    endif

    deallocate (f_w, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating f_w in plot_write_vdw_data', comm)
      return
    endif

    return

  end subroutine plot_write_vdw_data

  !================================================!
  subroutine plot_svd_omega_i(num_wann, num_kpts, kmesh_info, m_matrix, print_output, timer, &
                              dist_k, error, comm, stdout)
    !================================================!

    use w90_comms, only: w90_comm_type, mpirank, comms_allreduce
    use w90_constants, only: dp, cmplx_0
    use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_info_type, print_output_type, timer_list_type

    implicit none

    ! arguments
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    integer, intent(in) :: dist_k(:)
    integer, intent(in) :: num_wann, num_kpts
    integer, intent(in) :: stdout
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    ! local variables
    complex(kind=dp), allocatable :: cv1(:, :), cv2(:, :)
    complex(kind=dp), allocatable :: cw1(:), cw2(:)
    complex(kind=dp), allocatable :: cpad1(:)
    real(kind=dp), allocatable :: singvd(:)

    integer :: ierr, info, nkrank
    integer :: nkp, nn, nb, na, ind
    real(kind=dp) :: omt1, omt2, omt3

    nkrank = count(dist_k == mpirank(comm)) ! number k this rank, for dimensioning

    if (print_output%timing_level > 1 .and. print_output%iprint > 0) then
      call io_stopwatch_start('wann: svd_omega_i', timer)
    endif

    allocate (cw1(10*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cw1 in plot_svd_omega_i', comm)
      return
    endif
    allocate (cw2(10*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cw2 in plot_svd_omega_i', comm)
      return
    endif
    allocate (cv1(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cv1 in plot_svd_omega_i', comm)
      return
    endif
    allocate (cv2(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cv2 in plot_svd_omega_i', comm)
      return
    endif
    allocate (singvd(num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating singvd in plot_svd_omega_i', comm)
      return
    endif
    allocate (cpad1(num_wann*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cpad1 in plot_svd_omega_i', comm)
      return
    endif

    cw1 = cmplx_0; cw2 = cmplx_0; cv1 = cmplx_0; cv2 = cmplx_0; cpad1 = cmplx_0
    singvd = 0.0_dp

    ! singular value decomposition
    omt1 = 0.0_dp; omt2 = 0.0_dp; omt3 = 0.0_dp
    do nkp = 1, nkrank
      do nn = 1, kmesh_info%nntot
        ind = 1
        do nb = 1, num_wann
          do na = 1, num_wann
            cpad1(ind) = m_matrix(na, nb, nn, nkp)
            ind = ind + 1
          enddo
        enddo
        call zgesvd('A', 'A', num_wann, num_wann, cpad1, num_wann, singvd, cv1, &
                    num_wann, cv2, num_wann, cw1, 10*num_wann, cw2, info)
        if (info .ne. 0) then
          call set_error_fatal(error, 'ERROR: Singular value decomp. zgesvd failed', comm)
          return
        endif

        do nb = 1, num_wann
          omt1 = omt1 + kmesh_info%wb(nn)*(1.0_dp - singvd(nb)**2)
          omt2 = omt2 - kmesh_info%wb(nn)*(2.0_dp*log(singvd(nb)))
          omt3 = omt3 + kmesh_info%wb(nn)*(acos(singvd(nb))**2)
        enddo
      enddo
    enddo
    call comms_allreduce(omt1, 1, 'SUM', error, comm)
    call comms_allreduce(omt2, 1, 'SUM', error, comm)
    call comms_allreduce(omt3, 1, 'SUM', error, comm)
    omt1 = omt1/real(num_kpts, dp)
    omt2 = omt2/real(num_kpts, dp)
    omt3 = omt3/real(num_kpts, dp)
    if (print_output%iprint > 0) then
      write (stdout, *) ' '
      write (stdout, '(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ', &
        omt1*print_output%lenconfac**2, '('//trim(print_output%length_unit)//'^2)'
      write (stdout, '(2x,a,f15.9,1x,a)') '                 -2log s = ', &
        omt2*print_output%lenconfac**2, '('//trim(print_output%length_unit)//'^2)'
      write (stdout, '(2x,a,f15.9,1x,a)') '                  acos^2 = ', &
        omt3*print_output%lenconfac**2, '('//trim(print_output%length_unit)//'^2)'
    endif

    deallocate (cpad1, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cpad1 in plot_svd_omega_i', comm)
      return
    endif
    deallocate (singvd, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating singvd in plot_svd_omega_i', comm)
      return
    endif
    deallocate (cv2, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cv2 in plot_svd_omega_i', comm)
      return
    endif
    deallocate (cv1, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cv1 in plot_svd_omega_i', comm)
      return
    endif
    deallocate (cw2, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cw2 in plot_svd_omega_i', comm)
      return
    endif
    deallocate (cw1, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cw1 in plot_svd_omega_i', comm)
      return
    endif

    if (print_output%timing_level > 1 .and. print_output%iprint > 0) then
      call io_stopwatch_stop('wann: svd_omega_i', timer)
    endif
  end subroutine plot_svd_omega_i

  !================================================!
  subroutine plot_write_xyz(translate_home_cell, num_wann, wannier_centres, real_lattice, &
                            atom_data, print_output, error, comm, stdout, seedname)
    !================================================!
    !
    ! Write xyz file with Wannier centres
    !
    !================================================!

    use w90_io, only: io_date
    use w90_constants, only: dp, cmplx_0
    use w90_utility, only: utility_translate_home
    use w90_error, only: w90_error_type, set_error_file
    use w90_types, only: atom_data_type, print_output_type

    implicit none

    type(atom_data_type), intent(in) :: atom_data
    type(print_output_type), intent(in) :: print_output
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    logical, intent(in) :: translate_home_cell
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: wannier_centres(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: iw, ind, xyz_unit, nsp, nat, ierr
    character(len=9) :: cdate, ctime
    real(kind=dp) :: wc(3, num_wann)

    wc = wannier_centres

    if (translate_home_cell) then
      do iw = 1, num_wann
        call utility_translate_home(wc(:, iw), real_lattice)
      enddo
    endif

    if (print_output%iprint > 2) then
      write (stdout, '(1x,a)') 'Final centres (translated to home cell for writing xyz file)'
      do iw = 1, num_wann
        write (stdout, '(2x, "WF centre", i5, 2x, "(", f10.6, ",", f10.6, ",", f10.6, " )")') &
          iw, (wc(ind, iw)*print_output%lenconfac, ind=1, 3)
      end do
      write (stdout, '(1x,a78)') repeat('-', 78)
      write (stdout, *)
    endif

    open (newunit=xyz_unit, file=trim(seedname)//'_centres.xyz', form='formatted', iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error opening file '//trim(seedname)//'_centres.xyz in wann_write_xyz', comm)
      return
    endif
    write (xyz_unit, '(i6)') num_wann + atom_data%num_atoms
    call io_date(cdate, ctime)
    write (xyz_unit, *) 'Wannier centres, written by Wannier90 on'//cdate//' at '//ctime
    do iw = 1, num_wann
      write (xyz_unit, '("X",6x,3(f14.8,3x))') (wc(ind, iw), ind=1, 3)
    end do
    do nsp = 1, atom_data%num_species
      do nat = 1, atom_data%species_num(nsp)
        write (xyz_unit, '(a2,5x,3(f14.8,3x))') atom_data%symbol(nsp), atom_data%pos_cart(:, nat, nsp)
      end do
    end do
    close (xyz_unit)

    write (stdout, '(/a)') ' Wannier centres written to file '//trim(seedname)//'_centres.xyz'

    return

  end subroutine plot_write_xyz

end module w90_plot_mod
