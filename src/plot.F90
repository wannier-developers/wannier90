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

module w90_plot
  !! This module handles various plots

  implicit none

  private
  public :: plot_main

contains

  !============================================!
  subroutine plot_main()
    !! Main plotting routine
    !============================================!

    use w90_constants, only: eps6
    use w90_io, only: stdout, io_stopwatch
    use w90_parameters, only: num_kpts, bands_plot, dos_plot, &
      kpt_latt, fermi_surface_plot, &
      wannier_plot, timing_level, write_bvec, &
      write_hr, write_rmn, write_tb, write_u_matrices
    use w90_hamiltonian, only: hamiltonian_get_hr, hamiltonian_write_hr, &
      hamiltonian_setup, hamiltonian_write_rmn, &
      hamiltonian_write_tb, nrpts, irvec
    use w90_ws_distance, only: done_ws_distance, ws_translate_dist, &
      ws_write_vec

    implicit none

    integer :: nkp
    logical :: have_gamma

    if (timing_level > 0) call io_stopwatch('plot: main', 1)

    ! Print the header only if there is something to plot
    if (bands_plot .or. dos_plot .or. fermi_surface_plot .or. write_hr .or. &
        wannier_plot .or. write_u_matrices .or. write_tb) then
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '|                               PLOTTING                                    |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, *)
    end if

    if (bands_plot .or. dos_plot .or. fermi_surface_plot .or. write_hr .or. &
        write_tb) then
      ! Check if the kmesh includes the gamma point
      have_gamma = .false.
      do nkp = 1, num_kpts
        if (all(abs(kpt_latt(:, nkp)) < eps6)) have_gamma = .true.
      end do
      if (.not. have_gamma) &
           write (stdout, '(1x,a)') '!!!! Kpoint grid does not include Gamma. '// &
           & ' Interpolation may be incorrect. !!!!'
      ! Transform Hamiltonian to WF basis
      !
      call hamiltonian_setup()
      !
      call hamiltonian_get_hr()
      !
      if (bands_plot) call plot_interpolate_bands
      !
      if (fermi_surface_plot) call plot_fermi_surface
      !
      if (write_hr) call hamiltonian_write_hr()
      !
      if (write_rmn) call hamiltonian_write_rmn()
      !
      if (write_tb) call hamiltonian_write_tb()
      !
      if (write_hr .or. write_rmn .or. write_tb) then
        if (.not. done_ws_distance) call ws_translate_dist(nrpts, irvec)
        call ws_write_vec(nrpts, irvec)
      end if
    end if

    if (wannier_plot) call plot_wannier

    if (write_bvec) call plot_bvec

    if (write_u_matrices) call plot_u_matrices

    if (timing_level > 0) call io_stopwatch('plot: main', 2)

  end subroutine plot_main

  !-----------------------------------!
  !----------------- Private Routines ------------------------------
  !-----------------------------------!

  !============================================!
  subroutine plot_interpolate_bands
    !============================================!
    !                                            !
    !! Plots the interpolated band structure
    !                                            !
    !============================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, stdout, io_file_unit, seedname, &
      io_time, io_stopwatch
    use w90_parameters, only: num_wann, bands_num_points, recip_metric, &
      bands_num_spec_points, timing_level, &
      bands_spec_points, bands_label, bands_plot_format, &
      bands_plot_mode, num_bands_project, bands_plot_project, &
      use_ws_distance
    use w90_hamiltonian, only: irvec, nrpts, ndegen, ham_r
    use w90_ws_distance, only: irdist_ws, wdist_ndeg, &
      ws_translate_dist

    implicit none

    complex(kind=dp), allocatable  :: ham_r_cut(:, :, :)
    complex(kind=dp), allocatable  :: ham_pack(:)
    complex(kind=dp)              :: fac
    complex(kind=dp), allocatable  :: ham_kprm(:, :)
    complex(kind=dp), allocatable  :: U_int(:, :)
    complex(kind=dp), allocatable  :: cwork(:)
    real(kind=dp), allocatable     :: rwork(:)
    real(kind=dp)      :: kpath_len(bands_num_spec_points/2)
    integer            :: kpath_pts(bands_num_spec_points/2)
    logical            :: kpath_print_first_point(bands_num_spec_points/2)
    real(kind=dp), allocatable :: xval(:)
    real(kind=dp), allocatable :: eig_int(:, :), plot_kpoint(:, :)
    real(kind=dp), allocatable :: bands_proj(:, :)
    real(kind=dp)        :: rdotk, vec(3), emin, emax, time0
    integer, allocatable :: irvec_cut(:, :)
    integer              :: irvec_max(3)
    integer              :: nrpts_cut
    integer, allocatable :: iwork(:), ifail(:)
    integer              :: info, i, j
    integer              :: irpt, nfound, loop_kpt, counter
    integer              :: loop_spts, total_pts, loop_i, nkp, ideg
    integer              :: num_paths, num_spts, ierr
    integer              :: bndunit, gnuunit, loop_w, loop_p
    character(len=20), allocatable   :: glabel(:)
    character(len=10), allocatable  :: xlabel(:)
    character(len=10), allocatable  :: ctemp(:)
    integer, allocatable :: idx_special_points(:)
    real(kind=dp), allocatable :: xval_special_points(:)

    !
    if (timing_level > 1) call io_stopwatch('plot: interpolate_bands', 1)
    !
    time0 = io_time()
    write (stdout, *)
    write (stdout, '(1x,a)') 'Calculating interpolated band-structure'
    write (stdout, *)
    !
    allocate (ham_pack((num_wann*(num_wann + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_pack in plot_interpolate_bands')
    allocate (ham_kprm(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_kprm in plot_interpolate_bands')
    allocate (U_int(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating U_int in plot_interpolate_bands')
    allocate (cwork(2*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwork in plot_interpolate_bands')
    allocate (rwork(7*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rwork in plot_interpolate_bands')
    allocate (iwork(5*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating iwork in plot_interpolate_bands')
    allocate (ifail(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ifail in plot_interpolate_bands')

    allocate (idx_special_points(bands_num_spec_points), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating idx_special_points in plot_interpolate_bands')
    allocate (xval_special_points(bands_num_spec_points), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating xval_special_points in plot_interpolate_bands')
    idx_special_points = -1
    xval_special_points = -1._dp
    !
    ! Work out how many points in the total path and the positions of the special points
    !
    num_paths = bands_num_spec_points/2

    kpath_print_first_point = .false.

    ! Loop over paths, set to False print_first_point if the starting point
    ! is the same as the ending point of the previous path.
    ! I skip the first path for which I always want to print the first point.
    kpath_print_first_point(1) = .true.
    do i = 2, num_paths
      ! If either the coordinates are different or the label is different, compute again the point
      ! (it will end up at the same x coordinate)
      if ((SUM((bands_spec_points(:, (i - 1)*2) - bands_spec_points(:, (i - 1)*2 + 1))**2) > 1.e-6) .or. &
          (TRIM(bands_label((i - 1)*2)) .ne. TRIM(bands_label((i - 1)*2 + 1)))) then
        kpath_print_first_point(i) = .true.
      end if
    enddo

    ! Count the total number of special points
    num_spts = num_paths
    do i = 1, num_paths
      if (kpath_print_first_point(i)) num_spts = num_spts + 1
    end do

    do loop_spts = 1, num_paths
      vec = bands_spec_points(:, 2*loop_spts) - bands_spec_points(:, 2*loop_spts - 1)
      kpath_len(loop_spts) = sqrt(dot_product(vec, (matmul(recip_metric, vec))))
      if (loop_spts == 1) then
        kpath_pts(loop_spts) = bands_num_points
      else
        kpath_pts(loop_spts) = nint(real(bands_num_points, dp)*kpath_len(loop_spts)/kpath_len(1))
        ! At least 1 point
        !if (kpath_pts(loop_spts) .eq. 0) kpath_pts(loop_spts) = 1
      end if
    end do
    total_pts = sum(kpath_pts)
    do i = 1, num_paths
      if (kpath_print_first_point(i)) total_pts = total_pts + 1
    end do

    allocate (plot_kpoint(3, total_pts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating plot_kpoint in plot_interpolate_bands')
    allocate (xval(total_pts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating xval in plot_interpolate_bands')
    allocate (eig_int(num_wann, total_pts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating eig_int in plot_interpolate_bands')
    allocate (bands_proj(num_wann, total_pts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating bands_proj in plot_interpolate_bands')
    allocate (glabel(num_spts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating num_spts in plot_interpolate_bands')
    allocate (xlabel(num_spts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating xlabel in plot_interpolate_bands')
    allocate (ctemp(bands_num_spec_points), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ctemp in plot_interpolate_bands')
    eig_int = 0.0_dp; bands_proj = 0.0_dp
    !
    ! Find the position of each kpoint in the path
    !
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
        plot_kpoint(:, counter) = bands_spec_points(:, 2*loop_spts - 1)

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
        plot_kpoint(:, counter) = bands_spec_points(:, 2*loop_spts - 1) + &
                                  (bands_spec_points(:, 2*loop_spts) - bands_spec_points(:, 2*loop_spts - 1))* &
                                  (real(loop_i, dp)/real(kpath_pts(loop_spts), dp))
      end do
      idx_special_points(2*loop_spts) = counter
      xval_special_points(2*loop_spts) = xval(counter)
    end do
    !xval(total_pts)=sum(kpath_len)
    plot_kpoint(:, total_pts) = bands_spec_points(:, bands_num_spec_points)
    !
    ! Write out the kpoints in the path
    !
    bndunit = io_file_unit()
    open (bndunit, file=trim(seedname)//'_band.kpt', form='formatted')
    write (bndunit, *) total_pts
    do loop_spts = 1, total_pts
      write (bndunit, '(3f12.6,3x,a)') (plot_kpoint(loop_i, loop_spts), loop_i=1, 3), "1.0"
    end do
    close (bndunit)
    !
    ! Write out information on high-symmetry points in the path
    !
    bndunit = io_file_unit()
    open (bndunit, file=trim(seedname)//'_band.labelinfo.dat', form='formatted')
    do loop_spts = 1, bands_num_spec_points
      if ((MOD(loop_spts, 2) .eq. 1) .and. (kpath_print_first_point((loop_spts + 1)/2) .eqv. .false.)) cycle
      write (bndunit, '(a,3x,I10,3x,4f18.10)') &
        bands_label(loop_spts), &
        idx_special_points(loop_spts), &
        xval_special_points(loop_spts), &
        (plot_kpoint(loop_i, idx_special_points(loop_spts)), loop_i=1, 3)
    end do
    close (bndunit)
    !
    ! Cut H matrix in real-space
    !
    if (index(bands_plot_mode, 'cut') .ne. 0) call plot_cut_hr()
    !
    ! Interpolate the Hamiltonian at each kpoint
    !
    if (use_ws_distance) then
      if (index(bands_plot_mode, 's-k') .ne. 0) then
        call ws_translate_dist(nrpts, irvec, force_recompute=.true.)
      elseif (index(bands_plot_mode, 'cut') .ne. 0) then
        call ws_translate_dist(nrpts_cut, irvec_cut, force_recompute=.true.)
      else
        call io_error('Error in plot_interpolate bands: value of bands_plot_mode not recognised')
      endif
    endif

    ! [lp] the s-k and cut codes are very similar when use_ws_distance is used, a complete
    !      mercge after this point is not impossible
    do loop_kpt = 1, total_pts
      ham_kprm = cmplx_0
      !
      if (index(bands_plot_mode, 's-k') .ne. 0) then
        do irpt = 1, nrpts
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
          if (use_ws_distance) then
            do j = 1, num_wann
            do i = 1, num_wann
              do ideg = 1, wdist_ndeg(i, j, irpt)
                rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), &
                                          real(irdist_ws(:, ideg, i, j, irpt), dp))
                fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(irpt)*wdist_ndeg(i, j, irpt), dp)
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
      elseif (index(bands_plot_mode, 'cut') .ne. 0) then
        do irpt = 1, nrpts_cut
! [lp] Shift the WF to have the minimum distance IJ, see also ws_distance.F90
          if (use_ws_distance) then
            do j = 1, num_wann
            do i = 1, num_wann
              do ideg = 1, wdist_ndeg(i, j, irpt)
                rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), &
                                          real(irdist_ws(:, ideg, i, j, irpt), dp))
                fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(wdist_ndeg(i, j, irpt), dp)
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
        call io_error('Error in plot_interpolate_bands')
      endif
      if (info > 0) then
        write (stdout, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
        call io_error('Error in plot_interpolate_bands')
      endif
      ! Compute projection onto WF if requested
      if (num_bands_project > 0) then
      do loop_w = 1, num_wann
        do loop_p = 1, num_wann
          if (any(bands_plot_project == loop_p)) then
            bands_proj(loop_w, loop_kpt) = bands_proj(loop_w, loop_kpt) + abs(U_int(loop_p, loop_w))**2
          end if
        end do
      end do
      end if
      !
    end do
    !
    ! Interpolation Finished!
    ! Now we write plotting files
    !
    emin = minval(eig_int) - 1.0_dp
    emax = maxval(eig_int) + 1.0_dp

    if (index(bands_plot_format, 'gnu') > 0) call plot_interpolate_gnuplot
    if (index(bands_plot_format, 'xmgr') > 0) call plot_interpolate_xmgrace

    write (stdout, '(1x,a,f11.3,a)') &
      'Time to calculate interpolated band structure ', io_time() - time0, ' (sec)'
    write (stdout, *)

    if (allocated(ham_r_cut)) deallocate (ham_r_cut, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ham_r_cut in plot_interpolate_bands')
    if (allocated(irvec_cut)) deallocate (irvec_cut, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating irvec_cut in plot_interpolate_bands')
    !
    if (timing_level > 1) call io_stopwatch('plot: interpolate_bands', 2)
    !
    if (allocated(idx_special_points)) deallocate (idx_special_points, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating idx_special_points in plot_interpolate_bands')
    if (allocated(xval_special_points)) deallocate (xval_special_points, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating xval_special_points in plot_interpolate_bands')

  contains

    !============================================!
    subroutine plot_cut_hr()
      !============================================!
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
      !============================================!

      use w90_constants, only: dp, cmplx_0, eps8
      use w90_io, only: io_error, stdout
      use w90_parameters, only: num_wann, mp_grid, real_lattice, &
        one_dim_dir, bands_plot_dim, &
        hr_cutoff, dist_cutoff, dist_cutoff_mode
      use w90_hamiltonian, only: wannier_centres_translated

      implicit none
      !
      integer :: nrpts_tmp
      integer :: one_dim_vec, two_dim_vec(2)
      integer :: i, j, n1, n2, n3, i1, i2, i3
      real(kind=dp), allocatable :: ham_r_tmp(:, :)
      real(kind=dp), allocatable :: shift_vec(:, :)
      real(kind=dp) :: dist_ij_vec(3)
      real(kind=dp) :: dist_vec(3)
      real(kind=dp) :: dist

      allocate (ham_r_tmp(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating ham_r_tmp in plot_cut_hr')

      irvec_max = maxval(irvec, DIM=2) + 1

      if (bands_plot_dim .ne. 3) then
        ! Find one_dim_vec which is parallel to one_dim_dir
        ! two_dim_vec - the other two lattice vectors
        ! Along the confined directions, take only irvec=0
        j = 0
        do i = 1, 3
          if (abs(abs(real_lattice(one_dim_dir, i)) &
                  - sqrt(dot_product(real_lattice(:, i), real_lattice(:, i)))) .lt. eps8) then
            one_dim_vec = i
            j = j + 1
          end if
        end do
        if (j .ne. 1) call io_error('Error: 1-d lattice vector not defined in plot_cut_hr')
        j = 0
        do i = 1, 3
          if (i .ne. one_dim_vec) then
            j = j + 1
            two_dim_vec(j) = i
          end if
        end do
        if (bands_plot_dim .eq. 1) then
          irvec_max(two_dim_vec(1)) = 0
          irvec_max(two_dim_vec(2)) = 0
        end if
        if (bands_plot_dim .eq. 2) irvec_max(one_dim_vec) = 0
      end if

      nrpts_cut = (2*irvec_max(1) + 1)*(2*irvec_max(2) + 1)*(2*irvec_max(3) + 1)
      allocate (ham_r_cut(num_wann, num_wann, nrpts_cut), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating ham_r_cut in plot_cut_hr')
      allocate (irvec_cut(3, nrpts_cut), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating irvec_cut in plot_cut_hr')
      allocate (shift_vec(3, nrpts_cut), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating shift_vec in plot_cut_hr')

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
        call io_error('Error in plot_cut_hr')
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
      if (index(dist_cutoff_mode, 'one_dim') > 0) then
        do i = 1, num_wann
          do j = 1, num_wann
            dist_ij_vec(one_dim_dir) = &
              wannier_centres_translated(one_dim_dir, i) - wannier_centres_translated(one_dim_dir, j)
            do irpt = 1, nrpts_cut
              dist_vec(one_dim_dir) = dist_ij_vec(one_dim_dir) + shift_vec(one_dim_dir, irpt)
              dist = abs(dist_vec(one_dim_dir))
              if (dist .gt. dist_cutoff) &
                ham_r_cut(j, i, irpt) = cmplx_0
            end do
          end do
        end do
      else if (index(dist_cutoff_mode, 'two_dim') > 0) then
        do i = 1, num_wann
          do j = 1, num_wann
            dist_ij_vec(:) = wannier_centres_translated(:, i) - wannier_centres_translated(:, j)
            do irpt = 1, nrpts_cut
              dist_vec(:) = dist_ij_vec(:) + shift_vec(:, irpt)
              dist_vec(one_dim_dir) = 0.0_dp
              dist = sqrt(dot_product(dist_vec, dist_vec))
              if (dist .gt. dist_cutoff) &
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
              if (dist .gt. dist_cutoff) &
                ham_r_cut(j, i, irpt) = cmplx_0
            end do
          end do
        end do
      end if

      do irpt = 1, nrpts_cut
        do i = 1, num_wann
          do j = 1, num_wann
            if (abs(ham_r_cut(j, i, irpt)) .lt. hr_cutoff) &
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

    !============================================!
    subroutine plot_interpolate_gnuplot
      !============================================!
      !                                            !
      !! Plots the interpolated band structure in gnuplot format
      !============================================!

      use w90_constants, only: dp
      use w90_io, only: io_file_unit, seedname
      use w90_parameters, only: num_wann, bands_num_spec_points, &
        bands_label, num_bands_project

      implicit none

      !
      bndunit = io_file_unit()
      open (bndunit, file=trim(seedname)//'_band.dat', form='formatted')
      gnuunit = io_file_unit()
      open (gnuunit, file=trim(seedname)//'_band.gnu', form='formatted')
      !
      ! Gnuplot format
      !
      do i = 1, num_wann
        do nkp = 1, total_pts
          if (num_bands_project > 0) then
            write (bndunit, '(3E16.8)') xval(nkp), eig_int(i, nkp), bands_proj(i, nkp)
          else
            write (bndunit, '(2E16.8)') xval(nkp), eig_int(i, nkp)
          end if
        enddo
        write (bndunit, *) ' '
      enddo
      close (bndunit)
      ! Axis labels
      glabel(1) = TRIM(bands_label(1))
      do i = 2, num_paths
        if (bands_label(2*(i - 1)) /= bands_label(2*(i - 1) + 1)) then
          glabel(i) = TRIM(bands_label(2*(i - 1)))//'|'//TRIM(bands_label(2*(i - 1) + 1))
        else
          glabel(i) = TRIM(bands_label(2*(i - 1)))
        end if
      end do
      glabel(num_paths + 1) = TRIM(bands_label(2*num_paths))
      ! gnu file
      write (gnuunit, 701) xval(total_pts), emin, emax
      do i = 1, num_paths - 1
        write (gnuunit, 705) sum(kpath_len(1:i)), emin, sum(kpath_len(1:i)), emax
      enddo
      write (gnuunit, 702, advance="no") TRIM(glabel(1)), 0.0_dp, &
        (TRIM(glabel(i + 1)), sum(kpath_len(1:i)), i=1, bands_num_spec_points/2 - 1)
      write (gnuunit, 703) TRIM(glabel(1 + bands_num_spec_points/2)), sum(kpath_len(:))
      write (gnuunit, *) 'plot ', '"'//trim(seedname)//'_band.dat', '"'
      close (gnuunit)

      if (num_bands_project > 0) then
        gnuunit = io_file_unit()
        open (gnuunit, file=trim(seedname)//'_band_proj.gnu', form='formatted')
        write (gnuunit, '(a)') '#File to plot a colour-mapped Bandstructure'
        write (gnuunit, '(a)') 'set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )'
        write (gnuunit, '(a)') 'unset ztics'
        write (gnuunit, '(a)') 'unset key'
        write (gnuunit, '(a)') '# can make pointsize smaller (~0.5). Too small and nothing is printed'
        write (gnuunit, '(a)') 'set pointsize 0.8'
        write (gnuunit, '(a)') 'set grid xtics'
        write (gnuunit, '(a)') 'set view 0,0'
        write (gnuunit, '(a,f9.5,a)') 'set xrange [0:', xval(total_pts), ']'
        write (gnuunit, '(a,f9.5,a,f9.5,a)') 'set yrange [', emin, ':', emax, ']'
        write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
          (glabel(i + 1), sum(kpath_len(1:i)), i=1, bands_num_spec_points/2 - 1)
        write (gnuunit, 703) glabel(1 + bands_num_spec_points/2), sum(kpath_len(:))

        write (gnuunit, '(a,a,a,a)') 'splot ', '"'//trim(seedname)//'_band.dat', '"', ' u 1:2:3 w p pt 13 palette'
        write (gnuunit, '(a)') '#use the next lines to make a nice figure for a paper'
        write (gnuunit, '(a)') '#set term postscript enhanced eps color lw 0.5 dl 0.5'
        write (gnuunit, '(a)') '#set pointsize 0.275'
      end if
      !
701   format('set style data dots', /, 'set nokey', /, 'set xrange [0:', F8.5, ']', /, 'set yrange [', F9.5, ' :', F9.5, ']')
702   format('set xtics (', :20('"', A, '" ', F8.5, ','))
703   format(A, '" ', F8.5, ')')
705   format('set arrow from ', F8.5, ',', F10.5, ' to ', F8.5, ',', F10.5, ' nohead')

    end subroutine plot_interpolate_gnuplot

    subroutine plot_interpolate_xmgrace
      !============================================!
      !                                            !
      !! Plots the interpolated band structure in Xmgrace format
      !============================================!

      use w90_io, only: io_file_unit, seedname, io_date
      use w90_parameters, only: num_wann, bands_num_spec_points

      implicit none

      character(len=9) :: cdate, ctime

      call io_date(cdate, ctime)

      ! Axis labels

      ! Switch any G to Gamma

      do i = 1, bands_num_spec_points
        if (bands_label(i) == 'G') then
          ctemp(i) = '\xG\0'
        else
          ctemp(i) = bands_label(i)
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

      gnuunit = io_file_unit()
      open (gnuunit, file=trim(seedname)//'_band.agr', form='formatted')
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
        write (gnuunit, '(a,i0,a,a)') '@    xaxis  ticklabel ', i - 1, ',', '"'//trim(adjustl(xlabel(i)))//'"'
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

  !===========================================================!
  subroutine plot_fermi_surface
    !===========================================================!
    !                                                           !
    !!  Prepares a Xcrysden bxsf file to view the fermi surface
    !                                                           !
    !===========================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
    use w90_io, only: io_error, stdout, io_file_unit, seedname, &
      io_date, io_time, io_stopwatch
    use w90_parameters, only: num_wann, fermi_surface_num_points, timing_level, &
      recip_lattice, nfermi, fermi_energy_list
    use w90_hamiltonian, only: irvec, nrpts, ndegen, ham_r

    implicit none

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
    character(len=9)     :: cdate, ctime
    !
    if (timing_level > 1) call io_stopwatch('plot: fermi_surface', 1)
    time0 = io_time()
    write (stdout, *)
    write (stdout, '(1x,a)') 'Calculating Fermi surface'
    write (stdout, *)
    !
    if (nfermi > 1) call io_error("Error in plot: nfermi>1. Set the fermi level " &
                                  //"using the input parameter 'fermi_level'")
    !
    allocate (ham_pack((num_wann*(num_wann + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_pack plot_fermi_surface')
    allocate (ham_kprm(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ham_kprm plot_fermi_surface')
    allocate (U_int(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating U_int in plot_fermi_surface')
    allocate (cwork(2*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwork in plot_fermi_surface')
    allocate (rwork(7*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rwork in plot_fermi_surface')
    allocate (iwork(5*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating iwork in plot_fermi_surface')
    allocate (ifail(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ifail in plot_fermi_surface')
    !
    npts_plot = (fermi_surface_num_points + 1)**3
    allocate (eig_int(num_wann, npts_plot), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating eig_int in plot_fermi_surface')
    eig_int = 0.0_dp
    U_int = (0.0_dp, 0.0_dp)
    !
    ikp = 0
    do loop_x = 1, fermi_surface_num_points + 1
      do loop_y = 1, fermi_surface_num_points + 1
        do loop_z = 1, fermi_surface_num_points + 1
          ikp = ikp + 1

          ham_kprm = cmplx_0
          do irpt = 1, nrpts
            rdotk = twopi*real((loop_x - 1)*irvec(1, irpt) + &
                               (loop_y - 1)*irvec(2, irpt) + (loop_z - 1)* &
                               irvec(3, irpt), dp)/real(fermi_surface_num_points, dp)
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
            call io_error('Error in plot_fermi_surface')
          endif
          if (info > 0) then
            write (stdout, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
            call io_error('Error in plot_fermi_surface')
          endif
        end do
      end do
    end do

    call io_date(cdate, ctime)
    bxsf_unit = io_file_unit()
    open (bxsf_unit, FILE=trim(seedname)//'.bxsf', STATUS='UNKNOWN', FORM='FORMATTED')
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
    write (bxsf_unit, *) fermi_surface_num_points + 1, fermi_surface_num_points + 1, fermi_surface_num_points + 1
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

    write (stdout, '(1x,a,f11.3,a)') 'Time to calculate interpolated Fermi surface ', io_time() - time0, ' (sec)'
    write (stdout, *)
    !
    if (timing_level > 1) call io_stopwatch('plot: fermi_surface', 2)
    !
    return

  end subroutine plot_fermi_surface

  !============================================!
  subroutine plot_wannier
    !============================================!
    !                                            !
    !! Plot the WF in Xcrysden format
    !! based on code written by Michel Posternak
    !                                            !
    !============================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi, cmplx_1
    use w90_io, only: io_error, stdout, io_file_unit, seedname, &
      io_date, io_stopwatch
    use w90_parameters, only: num_wann, num_bands, num_kpts, u_matrix, spin, &
      ngs => wannier_plot_supercell, kpt_latt, num_species, atoms_species_num, &
      atoms_symbol, atoms_pos_cart, num_atoms, real_lattice, have_disentangled, &
      ndimwin, lwindow, u_matrix_opt, num_wannier_plot, wannier_plot_list, &
      wannier_plot_mode, wvfn_formatted, timing_level, wannier_plot_format, &
      spinors, wannier_plot_spinor_mode, wannier_plot_spinor_phase

    implicit none

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

    logical :: have_file
    integer :: i, j, nsp, nat, nbnd, counter, ierr
    integer :: loop_kpt, ik, ix, iy, iz, nk, ngx, ngy, ngz, nxx, nyy, nzz
    integer :: loop_b, nx, ny, nz, npoint, file_unit, loop_w, num_inc
    integer :: ispinor
    character(len=11) :: wfnname
    character(len=60) :: wanxsf, wancube
    character(len=9)  :: cdate, ctime
    logical           :: inc_band(num_bands)
    !
    if (timing_level > 1) call io_stopwatch('plot: wannier', 1)
    !
    if (.not. spinors) then
      write (wfnname, 200) 1, spin
    else
      write (wfnname, 199) 1
    endif
    inquire (file=wfnname, exist=have_file)
    if (.not. have_file) call io_error('plot_wannier: file '//wfnname//' not found')

    file_unit = io_file_unit()
    if (wvfn_formatted) then
      open (unit=file_unit, file=wfnname, form='formatted')
      read (file_unit, *) ngx, ngy, ngz, nk, nbnd
    else
      open (unit=file_unit, file=wfnname, form='unformatted')
      read (file_unit) ngx, ngy, ngz, nk, nbnd
    end if
    close (file_unit)

200 format('UNK', i5.5, '.', i1)
199 format('UNK', i5.5, '.', 'NC')

    allocate (wann_func(-((ngs(1))/2)*ngx:((ngs(1) + 1)/2)*ngx - 1, &
                        -((ngs(2))/2)*ngy:((ngs(2) + 1)/2)*ngy - 1, &
                        -((ngs(3))/2)*ngz:((ngs(3) + 1)/2)*ngz - 1, num_wannier_plot), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating wann_func in plot_wannier')
    wann_func = cmplx_0
    if (spinors) then
      allocate (wann_func_nc(-((ngs(1))/2)*ngx:((ngs(1) + 1)/2)*ngx - 1, &
                             -((ngs(2))/2)*ngy:((ngs(2) + 1)/2)*ngy - 1, &
                             -((ngs(3))/2)*ngz:((ngs(3) + 1)/2)*ngz - 1, 2, num_wannier_plot), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating wann_func_nc in plot_wannier')
      wann_func_nc = cmplx_0
    endif
    if (.not. spinors) then
      if (have_disentangled) then
        allocate (r_wvfn_tmp(ngx*ngy*ngz, maxval(ndimwin)), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating r_wvfn_tmp in plot_wannier')
      end if
      allocate (r_wvfn(ngx*ngy*ngz, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating r_wvfn in plot_wannier')
    else
      if (have_disentangled) then
        allocate (r_wvfn_tmp_nc(ngx*ngy*ngz, maxval(ndimwin), 2), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating r_wvfn_tmp_nc in plot_wannier')
      end if
      allocate (r_wvfn_nc(ngx*ngy*ngz, num_wann, 2), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating r_wvfn_nc in plot_wannier')
    endif

    call io_date(cdate, ctime)
    do loop_kpt = 1, num_kpts

      inc_band = .true.
      num_inc = num_wann
      if (have_disentangled) then
        inc_band(:) = lwindow(:, loop_kpt)
        num_inc = ndimwin(loop_kpt)
      end if

      if (.not. spinors) then
        write (wfnname, 200) loop_kpt, spin
      else
        write (wfnname, 199) loop_kpt
      endif
      file_unit = io_file_unit()
      if (wvfn_formatted) then
        open (unit=file_unit, file=wfnname, form='formatted')
        read (file_unit, *) ix, iy, iz, ik, nbnd
      else
        open (unit=file_unit, file=wfnname, form='unformatted')
        read (file_unit) ix, iy, iz, ik, nbnd
      end if

      if ((ix /= ngx) .or. (iy /= ngy) .or. (iz /= ngz) .or. (ik /= loop_kpt)) then
        write (stdout, '(1x,a,a)') 'WARNING: mismatch in file', trim(wfnname)
        write (stdout, '(1x,5(a6,I5))') '   ix=', ix, '   iy=', iy, '   iz=', iz, '   ik=', ik, ' nbnd=', nbnd
        write (stdout, '(1x,5(a6,I5))') '  ngx=', ngx, '  ngy=', ngy, '  ngz=', ngz, '  kpt=', loop_kpt, 'bands=', num_bands
        call io_error('plot_wannier')
      end if

      if (have_disentangled) then
        counter = 1
        do loop_b = 1, num_bands
          if (counter > num_inc) exit
          if (wvfn_formatted) then
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
          if (wvfn_formatted) then
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
              do loop_w = 1, num_wannier_plot
                if (.not. spinors) then
                  wann_func(nxx, nyy, nzz, loop_w) = &
                    wann_func(nxx, nyy, nzz, loop_w) + &
                    u_matrix(loop_b, wannier_plot_list(loop_w), loop_kpt)*r_wvfn(npoint, loop_b)*catmp
                else
                  wann_func_nc(nxx, nyy, nzz, 1, loop_w) = &
                    wann_func_nc(nxx, nyy, nzz, 1, loop_w) + & ! up-spinor
                    u_matrix(loop_b, wannier_plot_list(loop_w), loop_kpt)*r_wvfn_nc(npoint, loop_b, 1)*catmp
                  wann_func_nc(nxx, nyy, nzz, 2, loop_w) = &
                    wann_func_nc(nxx, nyy, nzz, 2, loop_w) + & ! down-spinor
                    u_matrix(loop_b, wannier_plot_list(loop_w), loop_kpt)*r_wvfn_nc(npoint, loop_b, 2)*catmp
                  if (loop_b == num_wann) then ! last loop
                    upspinor = real(wann_func_nc(nxx, nyy, nzz, 1, loop_w)* &
                                    conjg(wann_func_nc(nxx, nyy, nzz, 1, loop_w)), dp)
                    dnspinor = real(wann_func_nc(nxx, nyy, nzz, 2, loop_w)* &
                                    conjg(wann_func_nc(nxx, nyy, nzz, 2, loop_w)), dp)
                    if (wannier_plot_spinor_phase) then
                      upphase = sign(1.0_dp, real(wann_func_nc(nxx, nyy, nzz, 1, loop_w), dp))
                      dnphase = sign(1.0_dp, real(wann_func_nc(nxx, nyy, nzz, 2, loop_w), dp))
                    else
                      upphase = 1.0_dp; dnphase = 1.0_dp
                    endif
                    select case (wannier_plot_spinor_mode)
                    case ('total')
                      wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(upspinor + dnspinor), 0.0_dp, dp)
                    case ('up')
                      wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(upspinor), 0.0_dp, dp)*upphase
                    case ('down')
                      wann_func(nxx, nyy, nzz, loop_w) = cmplx(sqrt(dnspinor), 0.0_dp, dp)*dnphase
                    case default
                      call io_error('plot_wannier: Invalid wannier_plot_spinor_mode '//trim(wannier_plot_spinor_mode))
                    end select
                    wann_func(nxx, nyy, nzz, loop_w) = wann_func(nxx, nyy, nzz, loop_w)/real(num_kpts, dp)
                  endif
                endif
              end do
            end do
          end do
        end do

      end do

    end do !loop over kpoints

    if (.not. spinors) then !!!!! For spinor Wannier functions, the steps below are not necessary.
      ! fix the global phase by setting the wannier to
      ! be real at the point where it has max. modulus

      do loop_w = 1, num_wannier_plot
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
      end do
      !
      ! Check the 'reality' of the WF
      !
      do loop_w = 1, num_wannier_plot
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
        write (stdout, '(6x,a,i4,7x,a,f11.6)') 'Wannier Function Num: ', wannier_plot_list(loop_w), &
          'Maximum Im/Re Ratio = ', ratmax
      end do
    endif !!!!!
    write (stdout, *) ' '
    if (wannier_plot_format .eq. 'xcrysden') then
      call internal_xsf_format()
    elseif (wannier_plot_format .eq. 'cube') then
      call internal_cube_format()
    else
      call io_error('wannier_plot_format not recognised in wannier_plot')
    endif

    if (timing_level > 1) call io_stopwatch('plot: wannier', 2)

    return

  contains

    !============================================!
    subroutine internal_cube_format()
      !============================================!
      !                                            !
      !! Write WFs in Gaussian cube format.
      !                                            !
      !============================================!

      use w90_constants, only: bohr
      use w90_parameters, only: recip_lattice, iprint, &
        wannier_plot_radius, wannier_centres, atoms_symbol, &
        wannier_plot_scale, atoms_pos_frac, num_atoms
      use w90_utility, only: utility_translate_home, &
        utility_cart_to_frac, utility_frac_to_cart

      implicit none

      real(kind=dp), allocatable :: wann_cube(:, :, :)
      real(kind=dp) :: rstart(3), rend(3), rlength(3), orig(3), dgrid(3)
      real(kind=dp) :: moda(3), modb(3)
      real(kind=dp) :: val_Q
      real(kind=dp) :: comf(3), wcf(3), diff(3), difc(3), dist
      integer :: ierr, iname, max_elements, iw
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

      allocate (atomic_Z(num_species), stat=ierr)
      if (ierr .ne. 0) call io_error('Error: allocating atomic_Z in wannier_plot')

      lmol = .false.
      lcrys = .false.
      if (index(wannier_plot_mode, 'mol') > 0) lmol = .true.      ! molecule mode
      if (index(wannier_plot_mode, 'crys') > 0) lcrys = .true.    ! crystal mode

      val_Q = 1.0_dp ! dummy value for cube file

      ! Assign atomic numbers to species
      max_elements = size(periodic_table)
      do isp = 1, num_species
        do iname = 1, max_elements
          if (atoms_symbol(isp) .eq. periodic_table(iname)) then
            atomic_Z(isp) = iname
            exit
          endif
        enddo
      end do

202   format(a, '_', i5.5, '.cube')

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
      comf(:) = 0.0_dp
      do isp = 1, num_species
        do iat = 1, atoms_species_num(isp)
          comf(:) = comf(:) + atoms_pos_frac(:, iat, isp)
        enddo
      enddo
      comf(:) = comf(:)/num_atoms

      ! Loop over WFs
      do loop_w = 1, num_wannier_plot

        wann_index = wannier_plot_list(loop_w)
        write (wancube, 202) trim(seedname), wann_index

        ! Find start and end of cube wrt simulation (home) cell origin
        do i = 1, 3
          ! ... in terms of distance along each lattice vector direction i
          rstart(i) = (wannier_centres(1, wann_index)*recip_lattice(i, 1) &
                       + wannier_centres(2, wann_index)*recip_lattice(i, 2) &
                       + wannier_centres(3, wann_index)*recip_lattice(i, 3))*moda(i)/twopi &
                      - twopi*wannier_plot_radius/(moda(i)*modb(i))
          rend(i) = (wannier_centres(1, wann_index)*recip_lattice(i, 1) &
                     + wannier_centres(2, wann_index)*recip_lattice(i, 2) &
                     + wannier_centres(3, wann_index)*recip_lattice(i, 3))*moda(i)/twopi &
                    + twopi*wannier_plot_radius/(moda(i)*modb(i))
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
        if (iprint > 3) then
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
          write (stdout, '(a,3f12.6)') 'wann_cen=', (wannier_centres(i, wann_index), i=1, 3)
        endif

        allocate (wann_cube(1:ilength(1), 1:ilength(2), 1:ilength(3)), stat=ierr)
        if (ierr .ne. 0) call io_error('Error: allocating wann_cube in wannier_plot')

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
            call io_error('Error plotting WF cube.')
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
              call io_error('Error plotting WF cube.')
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
                call io_error('Error plotting WF cube.')
              endif
              wann_cube(nxx, nyy, nzz) = real(wann_func(qxx, qyy, qzz, loop_w), dp)
            enddo
          enddo
        enddo

        ! WF centre in fractional coordinates
        call utility_cart_to_frac(wannier_centres(:, wann_index), wcf(:), recip_lattice)

        ! The vector (in fractional coordinates) from WF centre to "centre of mass"
        diff(:) = comf(:) - wcf(:)

        ! Corresponding nearest cell vector
        irdiff(:) = nint(diff(:))

        if (iprint > 3) then
          write (stdout, '(a,3f12.6)') 'wcf     =', (wcf(i), i=1, 3)
          write (stdout, '(a,3f12.6)') 'diff    =', (diff(i), i=1, 3)
          write (stdout, '(a,3i12)') 'irdiff  =', (irdiff(i), i=1, 3)
        endif

        if (lmol) then ! In "molecule mode" translate origin of cube to bring it in coincidence with the atomic positions
          orig(:) = orig(:) + real(irdiff(1), kind=dp)*real_lattice(1, :) &
                    + real(irdiff(2), kind=dp)*real_lattice(2, :) &
                    + real(irdiff(3), kind=dp)*real_lattice(3, :)
          if (iprint > 3) write (stdout, '(a,3f12.6,/)') 'orig-new=', (orig(i), i=1, 3)
        else ! In "crystal mode" count number of atoms within a given radius of wannier centre
          icount = 0
          do isp = 1, num_species
            do iat = 1, atoms_species_num(isp)
              do nzz = -ngs(3)/2, (ngs(3) + 1)/2
                do nyy = -ngs(2)/2, (ngs(2) + 1)/2
                  do nxx = -ngs(1)/2, (ngs(1) + 1)/2
                    diff(:) = atoms_pos_frac(:, iat, isp) - wcf(:) &
                              + (/real(nxx, kind=dp), real(nyy, kind=dp), real(nzz, kind=dp)/)
                    call utility_frac_to_cart(diff, difc, real_lattice)
                    dist = sqrt(difc(1)*difc(1) + difc(2)*difc(2) + difc(3)*difc(3))
                    if (dist .le. (wannier_plot_scale*wannier_plot_radius)) then
                      icount = icount + 1
                    endif
                  enddo
                enddo
              enddo
            enddo ! iat
          enddo ! isp
          if (iprint > 3) write (stdout, '(a,i12)') 'icount  =', icount
        endif

        ! Write cube file (everything in Bohr)
        file_unit = io_file_unit()
        open (unit=file_unit, file=trim(wancube), form='formatted', status='unknown')
        ! First two lines are comments
        write (file_unit, *) '     Generated by Wannier90 code http://www.wannier.org'
        write (file_unit, *) '     On ', cdate, ' at ', ctime
        ! Number of atoms, origin of cube (Cartesians) wrt simulation (home) cell
        if (lmol) then
          write (file_unit, '(i4,3f13.5)') num_atoms, orig(1)/bohr, orig(2)/bohr, orig(3)/bohr
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

        do isp = 1, num_species
          do iat = 1, atoms_species_num(isp)
            if (lmol) then ! In "molecule mode", write atomic coordinates as they appear in input file
              write (file_unit, '(i4,4f13.5)') atomic_Z(isp), val_Q, (atoms_pos_cart(i, iat, isp)/bohr, i=1, 3)
            else           ! In "crystal mode", write atoms in supercell within a given radius of Wannier centre
              do nzz = -ngs(3)/2, (ngs(3) + 1)/2
                do nyy = -ngs(2)/2, (ngs(2) + 1)/2
                  do nxx = -ngs(1)/2, (ngs(1) + 1)/2
                    diff(:) = atoms_pos_frac(:, iat, isp) - wcf(:) &
                              + (/real(nxx, kind=dp), real(nyy, kind=dp), real(nzz, kind=dp)/)
                    call utility_frac_to_cart(diff, difc, real_lattice)
                    dist = sqrt(difc(1)*difc(1) + difc(2)*difc(2) + difc(3)*difc(3))
                    if (dist .le. (wannier_plot_scale*wannier_plot_radius)) then
                      diff(:) = atoms_pos_frac(:, iat, isp) &
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
        if (ierr .ne. 0) call io_error('Error: deallocating wann_cube in wannier_plot')

      end do

      deallocate (atomic_Z, stat=ierr)
      if (ierr .ne. 0) call io_error('Error: deallocating atomic_Z in wannier_plot')

      return

    end subroutine internal_cube_format

    subroutine internal_xsf_format()

      implicit none

201   format(a, '_', i5.5, '.xsf')

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

      do loop_b = 1, num_wannier_plot

        write (wanxsf, 201) trim(seedname), wannier_plot_list(loop_b)

        file_unit = io_file_unit()
        open (unit=file_unit, file=trim(wanxsf), form='formatted', status='unknown')
        write (file_unit, *) '      #'
        write (file_unit, *) '      # Generated by the Wannier90 code http://www.wannier.org'
        write (file_unit, *) '      # On ', cdate, ' at ', ctime
        write (file_unit, *) '      #'
        ! should pass this into the code
        if (index(wannier_plot_mode, 'mol') > 0) then
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
          write (file_unit, '(i6,"  1")') num_atoms
        endif
        do nsp = 1, num_species
          do nat = 1, atoms_species_num(nsp)
            write (file_unit, '(a2,3x,3f12.7)') atoms_symbol(nsp), (atoms_pos_cart(i, nat, nsp), i=1, 3)
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

      return

    end subroutine internal_xsf_format

  end subroutine plot_wannier

  !============================================!
  subroutine plot_u_matrices
    !============================================!
    !                                            !
    !! Plot u_matrix and u_matrix_opt to textfiles in readable format
    !                                            !
    !============================================!

    use w90_parameters, only: num_bands, num_kpts, num_wann, have_disentangled, &
      kpt_latt, u_matrix, u_matrix_opt
    use w90_io, only: io_error, stdout, io_file_unit, seedname, &
      io_time, io_stopwatch, io_date

    implicit none
    integer             :: matunit
    integer             :: i, j, nkp
    character(len=33)  :: header
    character(len=9)   :: cdate, ctime

    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime

    matunit = io_file_unit()
    open (matunit, file=trim(seedname)//'_u.mat', form='formatted')

    write (matunit, *) header
    write (matunit, *) num_kpts, num_wann, num_wann

    do nkp = 1, num_kpts
      write (matunit, *)
      write (matunit, '(f15.10,sp,f15.10,sp,f15.10)') kpt_latt(:, nkp)
      write (matunit, '(f15.10,sp,f15.10)') ((u_matrix(i, j, nkp), i=1, num_wann), j=1, num_wann)
    end do
    close (matunit)

    if (have_disentangled) then
      matunit = io_file_unit()
      open (matunit, file=trim(seedname)//'_u_dis.mat', form='formatted')
      write (matunit, *) header
      write (matunit, *) num_kpts, num_wann, num_bands
      do nkp = 1, num_kpts
        write (matunit, *)
        write (matunit, '(f15.10,sp,f15.10,sp,f15.10)') kpt_latt(:, nkp)
        write (matunit, '(f15.10,sp,f15.10)') ((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann)
      end do
      close (matunit)
    endif

  end subroutine plot_u_matrices

  !============================================!
  subroutine plot_bvec()
    !!
    !! June 2018: RM and SP
    !! Write to file the matrix elements of bvector and their weights
    !! This is used by EPW to compute the velocity.
    !! You need "write_bvec = .true." in your wannier input
    !!
    !============================================!
    use w90_parameters, only: wb, bk, num_kpts, nntot
    use w90_io, only: io_error, io_file_unit, seedname, io_date
    !
    implicit none
    !
    integer            :: nkp, nn, file_unit
    character(len=33) :: header
    character(len=9)  :: cdate, ctime
    !
    file_unit = io_file_unit()
    call io_date(cdate, ctime)
    header = 'written on '//cdate//' at '//ctime
    !
    open (file_unit, file=trim(seedname)//'.bvec', form='formatted', status='unknown', err=101)
    write (file_unit, *) header ! Date and time
    write (file_unit, *) num_kpts, nntot
    do nkp = 1, num_kpts
      do nn = 1, nntot
        write (file_unit, '(4F14.8)') bk(:, nn, nkp), wb(nn)
      enddo
    enddo
    close (file_unit)
    !
    return
    !
101 call io_error('Error: plot_bvec: problem opening file '//trim(seedname)//'.bvec')

  end subroutine plot_bvec

end module w90_plot

