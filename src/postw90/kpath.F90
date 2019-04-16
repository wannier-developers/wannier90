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

module w90_kpath

  !! Calculates quantities along a specified k-path:
  !!
  !!  - Energy bands (eventually colored by the spin)
  !!  - Energy bands (colored by Berry curvature-like
  !!                  term of spin Hall conductivity)
  !!
  !!  - (Berry curvature)x(-1) summed over occupied bands
  !!  - Berry curvature-like term of spin Hall conductivity
  !!    summed over occupied bands
  !!
  !!  - Integrand of orbital magnetization Morb=LCtil+ICtil

  use w90_constants, only: dp

  implicit none

  private

  public :: k_path

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       !
  !===========================================================!

  subroutine k_path
    !! Main routine

    use w90_comms
    use w90_constants, only: dp, cmplx_0, cmplx_i, twopi, eps8
    use w90_io, only: io_error, io_file_unit, seedname, &
      io_time, io_stopwatch, stdout
    use w90_utility, only: utility_diagonalize
    use w90_postw90_common, only: pw90common_fourier_R_to_k
    use w90_parameters, only: num_wann, kpath_task, &
      bands_num_spec_points, bands_label, &
      kpath_bands_colour, nfermi, fermi_energy_list, &
      berry_curv_unit, shc_alpha, shc_beta, shc_gamma, kubo_adpt_smr
    use w90_get_oper, only: get_HH_R, HH_R, get_AA_R, get_BB_R, get_CC_R, &
      get_FF_R, get_SS_R, get_SHC_R
    use w90_spin, only: spin_get_nk
    use w90_berry, only: berry_get_imf_klist, berry_get_imfgh_klist, &
      berry_get_shc_klist
    use w90_constants, only: bohr

    integer, dimension(0:num_nodes - 1) :: counts, displs

    integer           :: i, j, n, num_paths, num_spts, loop_path, loop_kpt, &
                         total_pts, counter, loop_i, dataunit, gnuunit, pyunit, &
                         my_num_pts
    real(kind=dp)     :: ymin, ymax, kpt(3), spn_k(num_wann), &
                         imf_k_list(3, 3, nfermi), img_k_list(3, 3, nfermi), &
                         imh_k_list(3, 3, nfermi), Morb_k(3, 3), &
                         range, zmin, zmax
    real(kind=dp)     :: shc_k_band(num_wann), shc_k_fermi(nfermi)
    real(kind=dp), allocatable, dimension(:) :: kpath_len
    logical           :: plot_bands, plot_curv, plot_morb, plot_shc
    character(len=120) :: file_name

    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: UU(:, :)
    real(kind=dp), allocatable    :: xval(:), eig(:, :), my_eig(:, :), &
                                     curv(:, :), my_curv(:, :), &
                                     morb(:, :), my_morb(:, :), &
                                     color(:, :), my_color(:, :), &
                                     plot_kpoint(:, :), my_plot_kpoint(:, :), &
                                     shc(:), my_shc(:)
    character(len=3), allocatable  :: glabel(:)

    ! Everything is done on the root node (not worthwhile parallelizing)
    ! However, we still have to read and distribute the data if we
    ! are in parallel. So calls to get_oper are done on all nodes at the moment
    !
    plot_bands = index(kpath_task, 'bands') > 0
    plot_curv = index(kpath_task, 'curv') > 0
    plot_morb = index(kpath_task, 'morb') > 0
    plot_shc = index(kpath_task, 'shc') > 0

    if (on_root) then
      if (plot_shc .or. (plot_bands .and. kpath_bands_colour == 'shc')) then
        ! not allowed to use adpt smr, since adpt smr needs berry_kmesh,
        ! see line 1837 of berry.F90
        if (kubo_adpt_smr) call io_error( &
          'Error: Must use fixed smearing when plotting spin Hall conductivity')
      end if
      if (plot_shc) then
        if (nfermi == 0) then
          call io_error('Error: must specify Fermi energy')
        else if (nfermi /= 1) then
          call io_error('Error: kpath plot only accept one Fermi energy, ' &
                        //'use fermi_energy instead of fermi_energy_min')
        end if
      end if
    end if

    call k_path_print_info(plot_bands, plot_curv, plot_morb, plot_shc)

    ! Set up the needed Wannier matrix elements
    call get_HH_R
    if (plot_curv .or. plot_morb) call get_AA_R
    if (plot_morb) then
      call get_BB_R
      call get_CC_R
    endif

    if (plot_shc .or. (plot_bands .and. kpath_bands_colour == 'shc')) then
      call get_AA_R
      call get_SS_R
      call get_SHC_R
    endif

    if (plot_bands .and. kpath_bands_colour == 'spin') call get_SS_R

    if (on_root) then
      ! Determine the number of k-points (total_pts) as well as
      ! their reciprocal-lattice coordinates long the path (plot_kpoint)
      ! and their associated horizontal coordinate for the plot (xval)
      call k_path_get_points(num_paths, kpath_len, total_pts, xval, plot_kpoint)
      ! (paths)
      num_spts = num_paths + 1 ! number of path endpoints (special pts)
    else
      ! Dummy allocation for making scatterv work
      allocate (plot_kpoint(1, 1))
    end if

    ! Broadcast number of k-points on the path
    call comms_bcast(total_pts, 1)

    ! Partition set of k-points into junks
    call comms_array_split(total_pts, counts, displs); 
    !kpt_lo = displs(my_node_id)+1
    !kpt_hi = displs(my_node_id)+counts(my_node_id)
    my_num_pts = counts(my_node_id)

    ! Distribute coordinates
    allocate (my_plot_kpoint(3, my_num_pts))
    call comms_scatterv(my_plot_kpoint, 3*my_num_pts, &
                        plot_kpoint, 3*counts, 3*displs)

    ! Value of the vertical coordinate in the actual plots: energy bands
    !
    if (plot_bands) then
      allocate (HH(num_wann, num_wann))
      allocate (UU(num_wann, num_wann))
      allocate (my_eig(num_wann, my_num_pts))
      if (kpath_bands_colour /= 'none') allocate (my_color(num_wann, my_num_pts))
    end if

    ! Value of the vertical coordinate in the actual plots
    !
    if (plot_curv) allocate (my_curv(my_num_pts, 3))
    if (plot_morb) allocate (my_morb(my_num_pts, 3))
    if (plot_shc) allocate (my_shc(my_num_pts))

    ! Loop over local junk of k-points on the path and evaluate the requested quantities
    !
    do loop_kpt = 1, my_num_pts
      kpt(:) = my_plot_kpoint(:, loop_kpt)

      if (plot_bands) then
        call pw90common_fourier_R_to_k(kpt, HH_R, HH, 0)
        call utility_diagonalize(HH, num_wann, my_eig(:, loop_kpt), UU)
        !
        ! Color-code energy bands with the spin projection along the
        ! chosen spin quantization axis
        !
        if (kpath_bands_colour == 'spin') then
          call spin_get_nk(kpt, spn_k)
          my_color(:, loop_kpt) = spn_k(:)
          !
          ! The following is needed to prevent bands from disappearing
          ! when the magnitude of the Wannier interpolated spn_k (very
          ! slightly) exceeds 1.0 (e.g. in bcc Fe along N--G--H)
          !
          do n = 1, num_wann
            if (my_color(n, loop_kpt) > 1.0_dp - eps8) then
              my_color(n, loop_kpt) = 1.0_dp - eps8
            else if (my_color(n, loop_kpt) < -1.0_dp + eps8) then
              my_color(n, loop_kpt) = -1.0_dp + eps8
            end if
          end do
        else if (kpath_bands_colour == 'shc') then
          call berry_get_shc_klist(kpt, shc_k_band=shc_k_band)
          my_color(:, loop_kpt) = shc_k_band
        end if
      end if

      if (plot_morb) then
        call berry_get_imfgh_klist(kpt, imf_k_list, img_k_list, imh_k_list)
        Morb_k = img_k_list(:, :, 1) + imh_k_list(:, :, 1) &
                 - 2.0_dp*fermi_energy_list(1)*imf_k_list(:, :, 1)
        Morb_k = -Morb_k/2.0_dp ! differs by -1/2 from Eq.97 LVTS12
        my_morb(loop_kpt, 1) = sum(Morb_k(:, 1))
        my_morb(loop_kpt, 2) = sum(Morb_k(:, 2))
        my_morb(loop_kpt, 3) = sum(Morb_k(:, 3))
      end if

      if (plot_curv) then
        if (.not. plot_morb) then
          call berry_get_imf_klist(kpt, imf_k_list)
        end if
        my_curv(loop_kpt, 1) = sum(imf_k_list(:, 1, 1))
        my_curv(loop_kpt, 2) = sum(imf_k_list(:, 2, 1))
        my_curv(loop_kpt, 3) = sum(imf_k_list(:, 3, 1))
      end if

      if (plot_shc) then
        call berry_get_shc_klist(kpt, shc_k_fermi=shc_k_fermi)
        my_shc(loop_kpt) = shc_k_fermi(1)
      end if
    end do !loop_kpt

    ! Send results to root process
    if (plot_bands) then
      allocate (eig(num_wann, total_pts))
      call comms_gatherv(my_eig, num_wann*my_num_pts, &
                         eig, num_wann*counts, num_wann*displs)
      if (kpath_bands_colour /= 'none') then
        allocate (color(num_wann, total_pts))
        call comms_gatherv(my_color, num_wann*my_num_pts, &
                           color, num_wann*counts, num_wann*displs)
      end if
    end if

    if (plot_curv) then
      allocate (curv(total_pts, 3))
      do i = 1, 3
        call comms_gatherv(my_curv(:, i), my_num_pts, &
                           curv(:, i), counts, displs)
      end do
    end if

    if (plot_morb) then
      allocate (morb(total_pts, 3))
      do i = 1, 3
        call comms_gatherv(my_morb(:, i), my_num_pts, &
                           morb(:, i), counts, displs)
      end do
    end if

    if (plot_shc) then
      allocate (shc(total_pts))
      call comms_gatherv(my_shc, my_num_pts, shc, counts, displs)
    end if

    if (on_root) then
      num_spts = num_paths + 1
      allocate (glabel(num_spts))

      if (plot_bands) then
        !
        ! Write out the kpoints in the path in a format that can be inserted
        ! directly in the pwscf input file (the '1.0_dp' in the second column is
        ! a k-point weight, expected by pwscf)
        !
        dataunit = io_file_unit()
        open (dataunit, file=trim(seedname)//'-path.kpt', form='formatted')
        write (dataunit, *) total_pts
        do loop_kpt = 1, total_pts
          write (dataunit, '(3f12.6,3x,f4.1)') &
            (plot_kpoint(loop_i, loop_kpt), loop_i=1, 3), 1.0_dp
        end do
        close (dataunit)
      end if
      if (plot_curv .and. berry_curv_unit == 'bohr2') curv = curv/bohr**2

      if (plot_bands .and. kpath_bands_colour == 'shc') then
        if (berry_curv_unit == 'bohr2') color = color/bohr**2
      end if
      if (plot_shc) then
        if (berry_curv_unit == 'bohr2') shc = shc/bohr**2
      end if

      ! Axis labels
      !
      glabel(1) = ' '//bands_label(1)//' '
      do i = 2, num_paths
        if (bands_label(2*(i - 1)) /= bands_label(2*(i - 1) + 1)) then
          glabel(i) = bands_label(2*(i - 1))//'/'//bands_label(2*(i - 1) + 1)
        else
          glabel(i) = ' '//bands_label(2*(i - 1))//' '
        end if
      end do
      glabel(num_spts) = ' '//bands_label(bands_num_spec_points)//' '

      ! Now write the plotting files

      write (stdout, '(/,1x,a)') 'Output files:'

      if (plot_bands) then
        file_name = trim(seedname)//'-path.kpt'
        write (stdout, '(/,3x,a)') file_name
        !
        ! Data file
        !
        dataunit = io_file_unit()
        file_name = trim(seedname)//'-bands.dat'
        write (stdout, '(/,3x,a)') file_name
        open (dataunit, file=file_name, form='formatted')
        do i = 1, num_wann
          do loop_kpt = 1, total_pts
            if (kpath_bands_colour == 'none') then
              write (dataunit, '(2E16.8)') xval(loop_kpt), eig(i, loop_kpt)
            else
              write (dataunit, '(3E16.8)') xval(loop_kpt), &
                eig(i, loop_kpt), color(i, loop_kpt)
            end if
          end do
          write (dataunit, *) ' '
        end do
        close (dataunit)
      end if

      if (plot_bands .and. .not. plot_curv .and. .not. plot_morb &
          .and. .not. plot_shc) then
        !
        ! Gnuplot script
        !
        ymin = minval(eig) - 1.0_dp
        ymax = maxval(eig) + 1.0_dp
        gnuunit = io_file_unit()
        file_name = trim(seedname)//'-bands.gnu'
        write (stdout, '(/,3x,a)') file_name
        open (gnuunit, file=file_name, form='formatted')
        do i = 1, num_paths - 1
          write (gnuunit, 705) sum(kpath_len(1:i)), ymin, &
            sum(kpath_len(1:i)), ymax
        end do
        if (kpath_bands_colour == 'none') then
          write (gnuunit, 701) xval(total_pts), ymin, ymax
          write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
            (glabel(i + 1), sum(kpath_len(1:i)), i=1, num_paths - 1)
          write (gnuunit, 703) glabel(1 + num_paths), sum(kpath_len(:))
          write (gnuunit, *) 'plot ', '"'//trim(seedname)//'-bands.dat', '"'
        else if (kpath_bands_colour == 'spin') then
          !
          ! Only works with gnuplot v4.2 and higher
          !
          write (gnuunit, 706) xval(total_pts), ymin, ymax
          write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
            (glabel(i + 1), sum(kpath_len(1:i)), i=1, num_paths - 1)
          write (gnuunit, 703) glabel(1 + num_paths), sum(kpath_len(:))
          write (gnuunit, *) &
            'set palette defined (-1 "blue", 0 "green", 1 "red")'
          write (gnuunit, *) 'set pm3d map'
          write (gnuunit, *) 'set zrange [-1:1]'
          write (gnuunit, *) 'splot ', '"'//trim(seedname)//'-bands.dat', &
            '" with dots palette'
        else if (kpath_bands_colour == 'shc') then
          !
          ! Only works with gnuplot v4.2 and higher
          !
          write (gnuunit, 706) xval(total_pts), ymin, ymax
          write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
            (glabel(i + 1), sum(kpath_len(1:i)), i=1, num_paths - 1)
          write (gnuunit, 703) glabel(1 + num_paths), sum(kpath_len(:))
          write (gnuunit, *) &
            'sgnlog10(x) = abs(x) > 10.0 ? log10(abs(x))*sgn(x) : x/10.0'
          ! merge: Fortran ternary operator: similar to ? : in C
          zmin = minval(color)
          zmin = merge(sign(log10(abs(zmin)), zmin), zmin/10.0_dp, abs(zmin) > 10.0_dp)
          zmax = maxval(color)
          zmax = merge(sign(log10(abs(zmax)), zmax), zmax/10.0_dp, abs(zmax) > 10.0_dp)
          write (gnuunit, *) &
            'set palette defined (', zmin, ' "blue", 0 "green", ', zmax, ' "red")'
          write (gnuunit, *) 'set pm3d map'
          write (gnuunit, *) 'set zrange [', zmin, ':', zmax, ']'
          write (gnuunit, *) 'splot ', '"'//trim(seedname)//'-bands.dat', &
            '" u 1:2:(sgnlog10($3)) with dots palette'
        end if
        close (gnuunit)
        !
        ! python script
        !
        pyunit = io_file_unit()
        file_name = trim(seedname)//'-bands.py'
        write (stdout, '(/,3x,a)') file_name
        open (pyunit, file=file_name, form='formatted')
        write (pyunit, '(a)') 'import pylab as pl'
        write (pyunit, '(a)') 'import numpy as np'
        write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
          "-bands.dat')"
        write (pyunit, '(a)') "x=data[:,0]"
        write (pyunit, '(a)') "y=data[:,1]"
        if (kpath_bands_colour == 'spin' &
            .or. kpath_bands_colour == 'shc') write (pyunit, '(a)') "z=data[:,2]"
        if (kpath_bands_colour == 'shc') write (pyunit, '(a)') &
          "z=np.array([np.log10(abs(elem))*np.sign(elem) " &
          //"if abs(elem)>10 else elem/10.0 for elem in z])"
        write (pyunit, '(a)') "tick_labels=[]"
        write (pyunit, '(a)') "tick_locs=[]"
        do j = 1, num_spts
          if (trim(glabel(j)) == ' G') then
            write (pyunit, '(a)') "tick_labels.append('$\Gamma$')"
          else
            write (pyunit, '(a)') "tick_labels.append('"//trim(glabel(j)) &
              //"'.strip())"
          end if
          if (j == 1) then
            write (pyunit, '(a,F12.6,a)') "tick_locs.append(0)"
          else
            write (pyunit, '(a,F12.6,a)') "tick_locs.append(", &
              sum(kpath_len(1:j - 1)), ")"
          end if
        end do
        if (kpath_bands_colour == 'none') then
          write (pyunit, '(a)') "pl.scatter(x,y,color='k',marker='+',s=0.1)"
        else if (kpath_bands_colour == 'spin' .or. &
                 kpath_bands_colour == 'shc') then
          write (pyunit, '(a)') &
            "pl.scatter(x,y,c=z,marker='+',s=1,cmap=pl.cm.jet)"
        end if
        write (pyunit, '(a)') "pl.xlim([0,max(x)])"
        write (pyunit, '(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y))," &
          //"max(y)+0.025*(max(y)-min(y))])"
        write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
        write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
        write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
          //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
          //"linestyle='-',linewidth=0.5)"
        write (pyunit, '(a)') "pl.ylabel('Energy [eV]')"
        if (kpath_bands_colour == 'spin' .or. kpath_bands_colour == 'shc') then
          write (pyunit, '(a)') &
            "pl.axes().set_aspect(aspect=0.65*max(x)/(max(y)-min(y)))"
          write (pyunit, '(a)') "pl.colorbar(shrink=0.7)"
        end if
        write (pyunit, '(a)') "outfile = '"//trim(seedname)//"-bands.pdf'"
        write (pyunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
        write (pyunit, '(a)') "pl.show()"

      end if ! plot_bands .and. .not.plot_curv .and. .not.plot_morb .and. .not. plot_shc

      if (plot_curv) then
        ! It is conventional to plot the negative curvature
        curv = -curv
        dataunit = io_file_unit()
        file_name = trim(seedname)//'-curv.dat'
        write (stdout, '(/,3x,a)') file_name
        open (dataunit, file=file_name, form='formatted')
        do loop_kpt = 1, total_pts
          write (dataunit, '(4E16.8)') xval(loop_kpt), &
            curv(loop_kpt, :)
        end do
        write (dataunit, *) ' '
        close (dataunit)
      end if

      if (plot_curv .and. .not. plot_bands) then

        do i = 1, 3
          !
          ! gnuplot script
          !
          gnuunit = io_file_unit()
          file_name = trim(seedname)//'-curv_'//achar(119 + i)//'.gnu'
          write (stdout, '(/,3x,a)') file_name
          open (gnuunit, file=file_name, form='formatted')
          ymin = minval(curv(:, i))
          ymax = maxval(curv(:, i))
          range = ymax - ymin
          ymin = ymin - 0.02_dp*range
          ymax = ymax + 0.02_dp*range
          write (gnuunit, 707) xval(total_pts), ymin, ymax
          do j = 1, num_paths - 1
            write (gnuunit, 705) sum(kpath_len(1:j)), ymin, &
              sum(kpath_len(1:j)), ymax
          end do
          write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
            (glabel(j + 1), sum(kpath_len(1:j)), j=1, num_paths - 1)
          write (gnuunit, 703) glabel(1 + num_paths), sum(kpath_len(:))
          write (gnuunit, *) &
            'plot ', '"'//trim(seedname)//'-curv.dat', '" u 1:'//achar(49 + i)
          close (gnuunit)
          !
          ! python script
          !
          pyunit = io_file_unit()
          file_name = trim(seedname)//'-curv_'//achar(119 + i)//'.py'
          write (stdout, '(/,3x,a)') file_name
          open (pyunit, file=file_name, form='formatted')
          write (pyunit, '(a)') 'import pylab as pl'
          write (pyunit, '(a)') 'import numpy as np'
          write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
            "-curv.dat')"
          write (pyunit, '(a)') "x=data[:,0]"
          write (pyunit, '(a)') "y=data[:,"//achar(48 + i)//"]"
          write (pyunit, '(a)') "tick_labels=[]"
          write (pyunit, '(a)') "tick_locs=[]"
          do j = 1, num_spts
            if (trim(glabel(j)) == ' G') then
              write (pyunit, '(a)') "tick_labels.append('$\Gamma$')"
            else
              write (pyunit, '(a)') "tick_labels.append('" &
                //trim(glabel(j))//"'.strip())"
            end if
            if (j == 1) then
              write (pyunit, '(a,F12.6,a)') "tick_locs.append(0)"
            else
              write (pyunit, '(a,F12.6,a)') "tick_locs.append(", &
                sum(kpath_len(1:j - 1)), ")"
            end if
          end do
          write (pyunit, '(a)') "pl.plot(x,y,color='k')"
          write (pyunit, '(a)') "pl.xlim([0,max(x)])"
          write (pyunit, '(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y))," &
            //"max(y)+0.025*(max(y)-min(y))])"
          write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
          write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
          write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
            //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
            //"linestyle='-',linewidth=0.5)"
          if (berry_curv_unit == 'ang2') then
            write (pyunit, '(a)') "pl.ylabel('$-\Omega_"//achar(119 + i) &
              //"(\mathbf{k})$  [ $\AA^2$ ]')"
          else if (berry_curv_unit == 'bohr2') then
            write (pyunit, '(a)') "pl.ylabel('$-\Omega_"//achar(119 + i) &
              //"(\mathbf{k})$  [ bohr$^2$ ]')"
          end if
          write (pyunit, '(a)') "outfile = '"//trim(seedname)// &
            "-curv_"//achar(119 + i)//".pdf'"
          write (pyunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
          write (pyunit, '(a)') "pl.show()"
        end do

      end if ! plot_curv .and. .not.plot_bands

      if (plot_morb) then
        dataunit = io_file_unit()
        file_name = trim(seedname)//'-morb.dat'
        write (stdout, '(/,3x,a)') file_name
        open (dataunit, file=file_name, form='formatted')
        do loop_kpt = 1, total_pts
          write (dataunit, '(4E16.8)') xval(loop_kpt), morb(loop_kpt, :)
        end do
        write (dataunit, *) ' '
        close (dataunit)
      end if

      if (plot_morb .and. .not. plot_bands) then
        do i = 1, 3
          !
          ! gnuplot script
          !
          gnuunit = io_file_unit()
          file_name = trim(seedname)//'-morb_'//achar(119 + i)//'.gnu'
          write (stdout, '(/,3x,a)') file_name
          open (gnuunit, file=file_name, form='formatted')
          ymin = minval(morb(:, i))
          ymax = maxval(morb(:, i))
          range = ymax - ymin
          ymin = ymin - 0.02_dp*range
          ymax = ymax + 0.02_dp*range
          write (gnuunit, 707) xval(total_pts), ymin, ymax
          do j = 1, num_paths - 1
            write (gnuunit, 705) sum(kpath_len(1:j)), ymin, &
              sum(kpath_len(1:j)), ymax
          end do
          write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
            (glabel(j + 1), sum(kpath_len(1:j)), j=1, num_paths - 1)
          write (gnuunit, 703) glabel(1 + num_paths), sum(kpath_len(:))
          write (gnuunit, *) &
            'plot ', '"'//trim(seedname)//'-morb.dat', '" u 1:'//achar(49 + i)
          close (gnuunit)
          !
          ! python script
          !
          pyunit = io_file_unit()
          file_name = trim(seedname)//'-morb_'//achar(119 + i)//'.py'
          write (stdout, '(/,3x,a)') file_name
          open (pyunit, file=file_name, form='formatted')
          write (pyunit, '(a)') 'import pylab as pl'
          write (pyunit, '(a)') 'import numpy as np'
          write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
            "-morb.dat')"
          write (pyunit, '(a)') "x=data[:,0]"
          write (pyunit, '(a)') "y=data[:,"//achar(48 + i)//"]"
          write (pyunit, '(a)') "tick_labels=[]"
          write (pyunit, '(a)') "tick_locs=[]"
          do j = 1, num_spts
            if (trim(glabel(j)) == ' G') then
              write (pyunit, '(a)') "tick_labels.append('$\Gamma$')"
            else
              write (pyunit, '(a)') &
                "tick_labels.append('"//trim(glabel(j))//"'.strip())"
            end if
            if (j == 1) then
              write (pyunit, '(a,F12.6,a)') "tick_locs.append(0)"
            else
              write (pyunit, '(a,F12.6,a)') "tick_locs.append(", &
                sum(kpath_len(1:j - 1)), ")"
            end if
          end do
          write (pyunit, '(a)') "pl.plot(x,y,color='k')"
          write (pyunit, '(a)') "pl.xlim([0,max(x)])"
          write (pyunit, '(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y))," &
            //"max(y)+0.025*(max(y)-min(y))])"
          write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
          write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
          write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
            //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
            //"linestyle='-',linewidth=0.5)"
          write (pyunit, '(a)') "pl.ylabel(r'$M^{\rm{orb}}_z(\mathbf{k})$" &
            //"  [ Ry$\cdot\AA^2$ ]')"
          write (pyunit, '(a)') "outfile = '"//trim(seedname)// &
            "-morb_"//achar(119 + i)//".pdf'"
          write (pyunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
          write (pyunit, '(a)') "pl.show()"
        end do

      end if ! plot_morb .and. .not.plot_bands

      if (plot_shc) then
        dataunit = io_file_unit()
        file_name = trim(seedname)//'-shc.dat'
        write (stdout, '(/,3x,a)') file_name
        open (dataunit, file=file_name, form='formatted')
        do loop_kpt = 1, total_pts
          write (dataunit, '(2E16.8)') xval(loop_kpt), &
            shc(loop_kpt)
        end do
        write (dataunit, *) ' '
        close (dataunit)
      end if

      if (plot_shc .and. .not. plot_bands) then
        !
        ! gnuplot script
        !
        gnuunit = io_file_unit()
        file_name = trim(seedname)//'-shc'//'.gnu'
        write (stdout, '(/,3x,a)') file_name
        open (gnuunit, file=file_name, form='formatted')
        ymin = minval(shc(:))
        ymax = maxval(shc(:))
        range = ymax - ymin
        ymin = ymin - 0.02_dp*range
        ymax = ymax + 0.02_dp*range
        write (gnuunit, 707) xval(total_pts), ymin, ymax
        do j = 1, num_paths - 1
          write (gnuunit, 705) sum(kpath_len(1:j)), ymin, &
            sum(kpath_len(1:j)), ymax
        end do
        write (gnuunit, 702, advance="no") glabel(1), 0.0_dp, &
          (glabel(j + 1), sum(kpath_len(1:j)), j=1, num_paths - 1)
        write (gnuunit, 703) glabel(1 + num_paths), sum(kpath_len(:))
        write (gnuunit, *) &
          'plot ', '"'//trim(seedname)//'-shc.dat', '" u 1:2'
        close (gnuunit)
        !
        ! python script
        !
        pyunit = io_file_unit()
        file_name = trim(seedname)//'-shc'//'.py'
        write (stdout, '(/,3x,a)') file_name
        open (pyunit, file=file_name, form='formatted')
        write (pyunit, '(a)') "# uncomment these two lines if you are " &
          //"running in non-GUI environment"
        write (pyunit, '(a)') "#import matplotlib"
        write (pyunit, '(a)') "#matplotlib.use('Agg')"
        write (pyunit, '(a)') 'import pylab as pl'
        write (pyunit, '(a)') 'import numpy as np'
        write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
          "-shc.dat')"
        write (pyunit, '(a)') "x=data[:,0]"
        write (pyunit, '(a)') "y=data[:,1]"
        write (pyunit, '(a)') "tick_labels=[]"
        write (pyunit, '(a)') "tick_locs=[]"
        do j = 1, num_spts
          if (trim(glabel(j)) == ' G') then
            write (pyunit, '(a)') "tick_labels.append('$\Gamma$')"
          else
            write (pyunit, '(a)') "tick_labels.append('" &
              //trim(glabel(j))//"'.strip())"
          end if
          if (j == 1) then
            write (pyunit, '(a,F12.6,a)') "tick_locs.append(0)"
          else
            write (pyunit, '(a,F12.6,a)') "tick_locs.append(", &
              sum(kpath_len(1:j - 1)), ")"
          end if
        end do
        write (pyunit, '(a)') "pl.plot(x,y,color='k')"
        write (pyunit, '(a)') "pl.xlim([0,max(x)])"
        write (pyunit, '(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y))," &
          //"max(y)+0.025*(max(y)-min(y))])"
        write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
        write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
        write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
          //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
          //"linestyle='-',linewidth=0.5)"
        if (berry_curv_unit == 'ang2') then
          write (pyunit, '(a)') "pl.ylabel('$\Omega_{" &
            //achar(119 + shc_alpha)//achar(119 + shc_beta) &
            //"}^{spin"//achar(119 + shc_gamma) &
            //"}(\mathbf{k})$  [ $\AA^2$ ]')"
        else if (berry_curv_unit == 'bohr2') then
          write (pyunit, '(a)') "pl.ylabel('$\Omega_{" &
            //achar(119 + shc_alpha)//achar(119 + shc_beta) &
            //"}^{spin"//achar(119 + shc_gamma) &
            //"}(\mathbf{k})$  [ bohr$^2$ ]')"
        end if
        write (pyunit, '(a)') "outfile = '"//trim(seedname)// &
          "-shc"//".pdf'"
        write (pyunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
        write (pyunit, '(a)') "pl.show()"

      end if ! plot_shc .and. .not.plot_bands

      if (plot_bands .and. plot_shc) then
        !
        ! python script
        !
        pyunit = io_file_unit()
        file_name = trim(seedname)//'-bands+shc'//'.py'
        write (stdout, '(/,3x,a)') file_name
        open (pyunit, file=file_name, form='formatted')
        write (pyunit, '(a)') "# uncomment these two lines if you are " &
          //"running in non-GUI environment"
        write (pyunit, '(a)') "#import matplotlib"
        write (pyunit, '(a)') "#matplotlib.use('Agg')"
        write (pyunit, '(a)') 'import pylab as pl'
        write (pyunit, '(a)') 'import numpy as np'
        write (pyunit, '(a)') 'from matplotlib.gridspec import GridSpec'
        write (pyunit, '(a)') "tick_labels=[]"
        write (pyunit, '(a)') "tick_locs=[]"
        do j = 1, num_spts
          if (trim(glabel(j)) == ' G') then
            write (pyunit, '(a)') "tick_labels.append('$\Gamma$')"
          else
            write (pyunit, '(a)') "tick_labels.append('"//trim(glabel(j)) &
              //"'.strip())"
          end if
          if (j == 1) then
            write (pyunit, '(a,F12.6,a)') "tick_locs.append(0)"
          else
            write (pyunit, '(a,F12.6,a)') "tick_locs.append(", &
              sum(kpath_len(1:j - 1)), ")"
          end if
        end do
        write (pyunit, '(a)') "fig = pl.figure()"
        write (pyunit, '(a)') "gs = GridSpec(2,52,hspace=0.00,wspace=1)"
        !
        ! upper panel (energy bands)
        !
        write (pyunit, '(a)') "axes1 = pl.subplot(gs[0, :-2])"
        write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
          "-bands.dat')"
        write (pyunit, '(a)') "x=data[:,0]"
        write (pyunit, '(a,F12.6)') "y=data[:,1]-", fermi_energy_list(1)
        if (kpath_bands_colour == 'spin' .or. kpath_bands_colour == 'shc') &
          write (pyunit, '(a)') "z=data[:,2]"
        if (kpath_bands_colour == 'shc') &
          write (pyunit, '(a)') "z=np.array([np.log10(abs(elem))*np.sign(elem) " &
          //"if abs(elem)>10 else elem/10.0 for elem in z])"
        if (kpath_bands_colour == 'none') then
          write (pyunit, '(a)') "pl.scatter(x,y,color='k',marker='+',s=0.1)"
        else if (kpath_bands_colour == 'spin' &
                 .or. kpath_bands_colour == 'shc') then
          write (pyunit, '(a)') &
            "pl.scatter(x,y,c=z,marker='+',s=1,cmap=pl.cm.jet)"
        end if
        write (pyunit, '(a)') "pl.xlim([0,max(x)])"
        write (pyunit, '(a)') "pl.ylim([-0.65,0.65]) # Adjust this range as needed"
        write (pyunit, '(a)') "pl.plot([tick_locs[0],tick_locs[-1]],[0,0]," &
          //"color='black',linestyle='--',linewidth=0.5)"
        write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
        write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
        write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
          //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
          //"linestyle='-',linewidth=0.5)"
        write (pyunit, '(a)') "pl.ylabel('Energy$-$E$_F$ [eV]')"
        write (pyunit, '(a)') "pl.tick_params(axis='x'," &
          //"which='both',bottom='off',top='off',labelbottom='off')"
        !
        ! color bar for upper panel
        !
        write (pyunit, '(a)') "# Now adding the colorbar"
        write (pyunit, '(a)') "cbaxes = pl.subplot(gs[:, -2:])"
        write (pyunit, '(a)') "cb = pl.colorbar(cax = cbaxes," &
          //"orientation='vertical') "
        write (pyunit, '(a)') "#cblim = int(min(abs(max(z)),abs(min(z))))"
        write (pyunit, '(a)') "#cb.set_ticks([-cblim,0,cblim])"
        write (pyunit, '(a)') "#cblim = [min(z),0,max(z)]"
        write (pyunit, '(a)') "#cb.set_ticks([min(z),0,max(z)])"
        write (pyunit, '(a)') "#cb.set_ticklabels(['-','0','+'])"
        !
        ! lower panel (SHC)
        !
        write (pyunit, '(a)') "axes2 = pl.subplot(gs[1, :-2])"
        write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
          "-shc.dat')"
        write (pyunit, '(a)') "x=data[:,0]"
        write (pyunit, '(a)') "y=data[:,1]"
        write (pyunit, '(a)') &
                "y=np.array([np.log10(abs(elem))*np.sign(elem) &
                        &if abs(elem)>10 else elem/10.0 for elem in y])"
        write (pyunit, '(a)') "pl.plot(x,y,color='k')"
        write (pyunit, '(a)') "pl.xlim([0,max(x)])"
        write (pyunit, '(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y))," &
          //"max(y)+0.025*(max(y)-min(y))])"
        write (pyunit, '(a)') "pl.plot([0,max(x)],[0,0],color='black'," &
          //"linestyle='--',linewidth=0.5)"
        write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
        write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
        write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
          //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
          //"linestyle='-',linewidth=0.5)"
        if (berry_curv_unit == 'ang2') then
          write (pyunit, '(a)') "pl.ylabel('$log_{10}|\Omega_{" &
            //achar(119 + shc_alpha)//achar(119 + shc_beta) &
            //"}^{spin"//achar(119 + shc_gamma) &
            //"}(\mathbf{k})|$  [ $\AA^2$ ]')"
        else if (berry_curv_unit == 'bohr2') then
          write (pyunit, '(a)') "pl.ylabel('$\Omega_{" &
            //achar(119 + shc_alpha)//achar(119 + shc_beta) &
            //"}^{spin"//achar(119 + shc_gamma) &
            //"}(\mathbf{k})$  [ bohr$^2$ ]')"
        end if
        write (pyunit, '(a)') "outfile = '"//trim(seedname)// &
          "-bands+shc"//".pdf'"
        write (pyunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
        write (pyunit, '(a)') "pl.show()"

      end if ! plot_bands .and. plot_shc

      if (plot_bands .and. (plot_curv .or. plot_morb)) then
        !
        ! python script
        !
        do i = 1, 3
          pyunit = io_file_unit()
          if (plot_curv) then
            file_name = trim(seedname)//'-bands+curv_'//achar(119 + i)//'.py'
          else if (plot_morb) then
            file_name = trim(seedname)//'-bands+morb_'//achar(119 + i)//'.py'
          end if
          write (stdout, '(/,3x,a)') file_name
          open (pyunit, file=file_name, form='formatted')
          write (pyunit, '(a)') 'import pylab as pl'
          write (pyunit, '(a)') 'import numpy as np'
          write (pyunit, '(a)') 'from matplotlib.gridspec import GridSpec'
          write (pyunit, '(a)') "tick_labels=[]"
          write (pyunit, '(a)') "tick_locs=[]"
          do j = 1, num_spts
            if (trim(glabel(j)) == ' G') then
              write (pyunit, '(a)') "tick_labels.append('$\Gamma$')"
            else
              write (pyunit, '(a)') "tick_labels.append('"//trim(glabel(j)) &
                //"'.strip())"
            end if
            if (j == 1) then
              write (pyunit, '(a,F12.6,a)') "tick_locs.append(0)"
            else
              write (pyunit, '(a,F12.6,a)') "tick_locs.append(", &
                sum(kpath_len(1:j - 1)), ")"
            end if
          end do
          write (pyunit, '(a)') "fig = pl.figure()"
          write (pyunit, '(a)') "gs = GridSpec(2, 1,hspace=0.00)"
          !
          ! upper panel (energy bands)
          !
          write (pyunit, '(a)') "axes1 = pl.subplot(gs[0, 0:])"
          write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
            "-bands.dat')"
          write (pyunit, '(a)') "x=data[:,0]"
          write (pyunit, '(a,F12.6)') "y=data[:,1]-", fermi_energy_list(1)
          if (kpath_bands_colour == 'spin') write (pyunit, '(a)') "z=data[:,2]"
          if (kpath_bands_colour == 'shc') then
            write (pyunit, '(a)') "z=data[:,2]"
            write (pyunit, '(a)') "z=np.array([np.log10(abs(elem))*np.sign(elem) " &
              //"if abs(elem)>10 else elem/10.0 for elem in z])"
          end if
          if (kpath_bands_colour == 'none') then
            write (pyunit, '(a)') "pl.scatter(x,y,color='k',marker='+',s=0.1)"
          else if (kpath_bands_colour == 'spin' .or. kpath_bands_colour == 'shc') then
            write (pyunit, '(a)') &
              "pl.scatter(x,y,c=z,marker='+',s=1,cmap=pl.cm.jet)"
          end if
          write (pyunit, '(a)') "pl.xlim([0,max(x)])"
          write (pyunit, '(a)') "pl.ylim([-0.65,0.65]) # Adjust this range as needed"
          write (pyunit, '(a)') "pl.plot([tick_locs[0],tick_locs[-1]],[0,0]," &
            //"color='black',linestyle='--',linewidth=0.5)"
          write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
          write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
          write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
            //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
            //"linestyle='-',linewidth=0.5)"
          write (pyunit, '(a)') "pl.ylabel('Energy$-$E$_F$ [eV]')"
          write (pyunit, '(a)') "pl.tick_params(axis='x'," &
            //"which='both',bottom='off',top='off',labelbottom='off')"
          !
          ! lower panel (curvature or orbital magnetization)
          !
          write (pyunit, '(a)') "axes2 = pl.subplot(gs[1, 0:])"
          if (plot_curv) then
            write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
              "-curv.dat')"
          else if (plot_morb) then
            write (pyunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
              "-morb.dat')"
          end if
          write (pyunit, '(a)') "x=data[:,0]"
          write (pyunit, '(a)') "y=data[:,"//achar(48 + i)//"]"
          write (pyunit, '(a)') "pl.plot(x,y,color='k')"
          write (pyunit, '(a)') "pl.xlim([0,max(x)])"
          write (pyunit, '(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y))," &
            //"max(y)+0.025*(max(y)-min(y))])"
          write (pyunit, '(a)') "pl.xticks(tick_locs,tick_labels)"
          write (pyunit, '(a)') "for n in range(1,len(tick_locs)):"
          write (pyunit, '(a)') "   pl.plot([tick_locs[n],tick_locs[n]]," &
            //"[pl.ylim()[0],pl.ylim()[1]],color='gray'," &
            //"linestyle='-',linewidth=0.5)"
          if (plot_curv) then
            if (berry_curv_unit == 'ang2') then
              write (pyunit, '(a)') "pl.ylabel('$-\Omega_"//achar(119 + i) &
                //"(\mathbf{k})$  [ $\AA^2$ ]')"
            else if (berry_curv_unit == 'bohr2') then
              write (pyunit, '(a)') "pl.ylabel('$-\Omega_"//achar(119 + i) &
                //"(\mathbf{k})$  [ bohr$^2$ ]')"
            end if
            write (pyunit, '(a)') "outfile = '"//trim(seedname)// &
              "-bands+curv_"//achar(119 + i)//".pdf'"
          else if (plot_morb) then
            write (pyunit, '(a)') "pl.ylabel(r'$M^{\rm{orb}}_z(\mathbf{k})$" &
              //"  [ Ry$\cdot\AA^2$ ]')"
            write (pyunit, '(a)') "outfile = '"//trim(seedname)// &
              "-morb_"//achar(119 + i)//".pdf'"
          end if
          write (pyunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
          write (pyunit, '(a)') "pl.show()"
        end do

      end if ! plot_bands .and. plot_curv

    end if ! on_root

701 format('set style data dots', /, 'unset key', /, &
           'set xrange [0:', F8.5, ']', /, 'set yrange [', F16.8, ' :', F16.8, ']')
702 format('set xtics (', :20('"', A3, '" ', F8.5, ','))
703 format(A3, '" ', F8.5, ')')
704 format('set palette defined (', F8.5, ' "red", 0 "green", ', F8.5, ' "blue")')
705 format('set arrow from ', F16.8, ',', F16.8, ' to ', F16.8, ',', F16.8, ' nohead')
706 format('unset key', /, &
           'set xrange [0:', F9.5, ']', /, 'set yrange [', F16.8, ' :', F16.8, ']')
707 format('set style data lines', /, 'set nokey', /, &
           'set xrange [0:', F8.5, ']', /, 'set yrange [', F16.8, ' :', F16.8, ']')

  end subroutine k_path

  !===========================================================!
  !                   PRIVATE PROCEDURES                      !
  !===========================================================!
  subroutine k_path_print_info(plot_bands, plot_curv, plot_morb, plot_shc)

    use w90_comms, only: on_root
    use w90_parameters, only: kpath_bands_colour, berry_curv_unit, nfermi
    use w90_io, only: stdout, io_error

    logical, intent(in)      :: plot_bands, plot_curv, plot_morb, plot_shc

    if (on_root) then
      write (stdout, '(/,/,1x,a)') &
        'Properties calculated in module  k p a t h'
      write (stdout, '(1x,a)') &
        '------------------------------------------'

      if (plot_bands) then
        select case (kpath_bands_colour)
        case ("none")
          write (stdout, '(/,3x,a)') '* Energy bands in eV'
        case ("spin")
          write (stdout, '(/,3x,a)') '* Energy bands in eV, coloured by spin'
        case ("shc")
          write (stdout, '(/,3x,a)') '* Energy bands in eV, coloured by SHC'
        end select
      end if
      if (plot_curv) then
        if (berry_curv_unit == 'ang2') then
          write (stdout, '(/,3x,a)') '* Negative Berry curvature in Ang^2'
        else if (berry_curv_unit == 'bohr2') then
          write (stdout, '(/,3x,a)') '* Negative Berry curvature in Bohr^2'
        endif
        if (nfermi /= 1) call io_error( &
          'Must specify one Fermi level when kpath_task=curv')
      end if
      if (plot_morb) then
        write (stdout, '(/,3x,a)') &
          '* Orbital magnetization k-space integrand in eV.Ang^2'
        if (nfermi /= 1) call io_error( &
          'Must specify one Fermi level when kpath_task=morb')
      end if
      if (plot_shc) then
        if (berry_curv_unit == 'ang2') then
          write (stdout, '(/,3x,a)') '* Berry curvature-like term for' &
            //' spin Hall conductivity in Ang^2'
        else if (berry_curv_unit == 'bohr2') then
          write (stdout, '(/,3x,a)') '* Berry curvature-like term for' &
            //' spin Hall conductivity in Bohr^2'
        end if
        if (nfermi /= 1) call io_error( &
          'Must specify one Fermi level when kpath_task=shc')
      end if
    end if ! on_root

  end subroutine

  !===================================================================!
  subroutine k_path_get_points(num_paths, kpath_len, total_pts, xval, plot_kpoint)
    !===================================================================!
    ! Determine the number of k-points (total_pts) as well as           !
    ! their reciprocal-lattice coordinates long the path (plot_kpoint)  !
    ! and their associated horizontal coordinate for the plot (xval)    !
    !===================================================================!

    use w90_parameters, only: kpath_num_points, &
      bands_num_spec_points, &
      bands_spec_points, &
      recip_metric

    integer, intent(out)    :: num_paths, total_pts
    real(kind=dp), allocatable, dimension(:), intent(out)   :: kpath_len, xval
    real(kind=dp), allocatable, dimension(:, :), intent(out) :: plot_kpoint

    integer :: counter, loop_path, loop_i

    integer, allocatable, dimension(:) :: kpath_pts
    real(kind=dp)                      :: vec(3)

    ! Work out how many points there are in the total path, and the
    ! positions of the special points
    !
    num_paths = bands_num_spec_points/2 ! number of straight line segments
    allocate (kpath_pts(num_paths))
    allocate (kpath_len(num_paths))
    do loop_path = 1, num_paths
      vec = bands_spec_points(:, 2*loop_path) &
            - bands_spec_points(:, 2*loop_path - 1)
      kpath_len(loop_path) = &
        sqrt(dot_product(vec, (matmul(recip_metric, vec))))
      !
      ! kpath_pts(loop_path) is the number of points in path number
      ! loop_path (all segments have the same density of points)
      !
      if (loop_path == 1) then
        kpath_pts(loop_path) = kpath_num_points
      else
        kpath_pts(loop_path) = nint(real(kpath_num_points, dp) &
                                    *kpath_len(loop_path)/kpath_len(1))
      end if
    end do
    total_pts = sum(kpath_pts) + 1

    ! Reciprocal-lattice coordinates of the k-points along the path
    !
    allocate (plot_kpoint(3, total_pts))

    ! Value of the horizontal coordinate in the actual plots (units of
    ! distance in k-space)
    !
    allocate (xval(total_pts))

    ! Find the position of each kpoint along the path
    !
    counter = 0
    do loop_path = 1, num_paths
      do loop_i = 1, kpath_pts(loop_path)
        counter = counter + 1
        if (counter == 1) then
          xval(counter) = 0.0_dp
        else
          xval(counter) = xval(counter - 1) &
                          + kpath_len(loop_path)/real(kpath_pts(loop_path), dp)
        end if
        plot_kpoint(:, counter) = bands_spec_points(:, 2*loop_path - 1) &
                                  + (bands_spec_points(:, 2*loop_path) &
                                     - bands_spec_points(:, 2*loop_path - 1) &
                                     ) &
                                  *(real(loop_i - 1, dp)/real(kpath_pts(loop_path), dp))
      end do
    end do
    !
    ! Last point
    !
    xval(total_pts) = sum(kpath_len)
    plot_kpoint(:, total_pts) = bands_spec_points(:, bands_num_spec_points)

  end subroutine

end module w90_kpath
