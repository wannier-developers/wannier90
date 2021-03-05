module w90_lib

  public wannier_setup, wannier_run

contains

!-*- mode: F90 -*-!
!------------------------------------------------------------!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Please cite                                                !
!                                                            !
!  [ref] "Wannier90 as a community code:                     !
!        new features and applications",                     !
!        G. Pizzi et al.,  J. Phys. Cond. Matt. 32,          !
!        165902 (2020).                                      !
!        http://doi.org/10.1088/1361-648X/ab51ff             !
!                                                            !
! in any publications arising from the use of this code.     !
!                                                            !
! Wannier90 is based on Wannier77, written by N. Marzari,    !
! I. Souza and D. Vanderbilt. For the method please cite     !
!                                                            !
! [ref] N. Marzari and D. Vanderbilt,                        !
!       Phys. Rev. B 56 12847 (1997)                         !
!       http://dx.doi.org/10.1103/PhysRevB.56.12847          !
!                                                            !
! [ref] I. Souza, N. Marzari and D. Vanderbilt,              !
!       Phys. Rev. B 65 035109 (2001)                        !
!       http://dx.doi.org/10.1103/PhysRevB.65.035109         !
!                                                            !
! [ref] N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza,    !
!       D. Vanderbilt, "Maximally localized Wannier          !
!       functions: theory and applications",                 !
!       Rev. Mod. Phys. 84, 1419 (2012)                      !
!       http://dx.doi.org/10.1103/RevModPhys.84.1419         !
!                                                            !
! For a full list of authors and contributors, please        !
! see the README file in the root directory of the           !
! distribution.                                              !
!                                                            !
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

  subroutine wannier_setup(seed__name, mp_grid_loc, num_kpts_loc, &
                           real_lattice_loc, recip_lattice_loc, kpt_latt_loc, &
                           num_bands_tot, num_atoms_loc, atom_symbols_loc, &
                           atoms_cart_loc, gamma_only_loc, spinors_loc, &
                           nntot_loc, nnlist_loc, nncell_loc, num_bands_loc, &
                           num_wann_loc, proj_site_loc, proj_l_loc, proj_m_loc, &
                           proj_radial_loc, proj_z_loc, proj_x_loc, &
                           proj_zona_loc, exclude_bands_loc, proj_s_loc, &
                           proj_s_qaxis_loc)
    use w90_constants
    use w90_parameters
    use w90_sitesym
    use w90_io
    use w90_kmesh
    use w90_comms, only: comms_setup_vars, on_root, num_nodes, my_node_id

    use w90_param_methods, only: param_write_header, param_lib_set_atoms, &
      param_read, param_write, param_dealloc, param_dist

    implicit none

    integer, intent(in) :: num_atoms_loc
    integer, intent(in) :: num_bands_tot
    integer, intent(in) :: num_kpts_loc
    integer, intent(inout) :: nntot_loc  !JJ, was "out", but don't set yet
    integer, intent(inout) :: num_bands_loc  !JJ, was "out", but don't set yet
    integer, intent(inout) :: num_wann_loc  !JJ, was "out", but don't set yet
    logical, intent(in) :: gamma_only_loc
    logical, intent(in) :: spinors_loc

    integer, dimension(3), intent(in) :: mp_grid_loc
    integer, dimension(3, num_kpts_loc, num_nnmax), intent(inout) :: nncell_loc  !JJ, was "out", but don't set yet
    integer, dimension(num_bands_tot), intent(inout) :: exclude_bands_loc  !JJ, was "out", but don't set yet
    integer, dimension(num_bands_tot), intent(inout) :: proj_l_loc  !JJ, was "out", but don't set yet
    integer, dimension(num_bands_tot), intent(inout) :: proj_m_loc  !JJ, was "out", but don't set yet
    integer, dimension(num_bands_tot), intent(inout) :: proj_radial_loc  !JJ, was "out", but don't set yet
    ! JJ fixme we need a proper interface
    !integer, dimension(num_bands_tot), optional, intent(inout) :: proj_s_loc !JJ, was "out", but don't set yet
    integer, dimension(num_bands_tot), intent(inout) :: proj_s_loc !JJ, was "out", but don't set yet
    integer, dimension(num_kpts_loc, num_nnmax), intent(inout) :: nnlist_loc !JJ, was "out", but don't set yet

    real(kind=dp), dimension(3, 3), intent(in) :: real_lattice_loc
    real(kind=dp), dimension(3, 3), intent(in) :: recip_lattice_loc
    real(kind=dp), dimension(3, num_atoms_loc), intent(in) :: atoms_cart_loc
    real(kind=dp), dimension(3, num_bands_tot), intent(inout) :: proj_site_loc !JJ, was "out", but don't set yet
    real(kind=dp), dimension(3, num_bands_tot), intent(inout) :: proj_x_loc !JJ, was "out", but don't set yet
    real(kind=dp), dimension(3, num_bands_tot), intent(inout) :: proj_z_loc !JJ, was "out", but don't set yet
    !real(kind=dp), dimension(3, num_bands_tot), optional, intent(inout) :: proj_s_qaxis_loc !JJ, was "out", but don't set yet
    real(kind=dp), dimension(3, num_bands_tot), intent(inout) :: proj_s_qaxis_loc !JJ, was "out", but don't set yet
    real(kind=dp), dimension(3, num_kpts_loc), intent(in) :: kpt_latt_loc
    real(kind=dp), dimension(num_bands_tot), intent(inout) :: proj_zona_loc !JJ, was "out", but don't set yet

    character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
    character(len=*), intent(in) :: seed__name

    integer :: ierr
    real(kind=dp) time1, time2
    character(len=9) :: stat, pos
    logical :: wout_found

    time2 = io_time()

    call comms_setup_vars

    if (on_root) then

      call param_read(driver, w90_calcs, pp_calc, param_input, param_plot, &
                      param_wannierise, lsitesymmetry, symmetrize_eps, &
                      wann_data, param_hamil, kmesh_data, kmesh_info, k_points, &
                      num_kpts, dis_data, fermi_surface_data, fermi, tran, atoms, &
                      num_bands, num_wann, eigval, mp_grid, num_proj, select_proj, &
                      real_lattice, recip_lattice, spec_points, pw90_calcs, &
                      postw90_oper, pw90_common, pw90_spin, pw90_ham, kpath, &
                      kslice, dos_data, berry, spin_hall, gyrotropic, geninterp, &
                      boltz, eig_found, .true.)

      seedname = trim(adjustl(seed__name))
      if (driver%restart .eq. ' ') then
        stat = 'replace'
        pos = 'rewind'
      else
        inquire (file=trim(seedname)//'.wout', exist=wout_found)
        if (wout_found) then
          stat = 'old'
        else
          stat = 'replace'
        endif
        pos = 'append'
      endif

      stdout = io_file_unit()
      open (unit=stdout, file=trim(seedname)//'.wout', status=trim(stat), position=trim(pos))

      call param_write_header()

      if (num_nodes == 1) then
#ifdef MPI
        write (stdout, '(/,1x,a)') 'Running in serial (with parallel executable)'
#else
        write (stdout, '(/,1x,a)') 'Running in serial (with serial executable)'
#endif
      else
        write (stdout, '(/,1x,a,i3,a/)') 'Running in parallel on ', num_nodes, ' CPUs'
      endif

      !%%%% only time arguments are referenced
      ! copy local data into module variables
      if (driver%library) then
        mp_grid = mp_grid_loc
        num_kpts = num_kpts_loc
        real_lattice = real_lattice_loc
        recip_lattice = recip_lattice_loc
        allocate (k_points%kpt_latt(3, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating kpt_latt in wannier_setup')
        k_points%kpt_latt = kpt_latt_loc
        atoms%num_atoms = num_atoms_loc
        call param_lib_set_atoms(atoms, atom_symbols_loc, atoms_cart_loc, recip_lattice)
        param_input%gamma_only = gamma_only_loc
        param_input%spinors = spinors_loc
        ! GP: at this point we don't know yet the number of excluded bands...
        num_bands = num_bands_tot
      endif

      call param_write(driver, w90_calcs, param_input, param_plot, &
                       param_wannierise, lsitesymmetry, symmetrize_eps, &
                       wann_data, param_hamil, kmesh_data, k_points, num_kpts, &
                       dis_data, fermi_surface_data, fermi, tran, atoms, &
                       num_bands, num_wann, mp_grid, num_proj, select_proj, &
                       real_lattice, recip_lattice, spec_points, pw90_calcs)

      time1 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time2, ' (sec)'

      if (.not. driver%explicit_nnkpts) then
        call kmesh_get(recip_lattice, k_points%kpt_cart, &
                       param_input%timing_level, kmesh_info%nncell, &
                       kmesh_info%neigh, kmesh_info%nnlist, kmesh_info%nntot, &
                       kmesh_data%shell_list, param_input%devel_flag, &
                       param_input%iprint, param_input%lenconfac, &
                       kmesh_data%tol, num_kpts, kmesh_data%search_shells, &
                       param_input%gamma_only, kmesh_info%nnh, kmesh_info%wbtot, &
                       kmesh_data%skip_B1_tests, kmesh_info%bk, kmesh_info%bka, &
                       kmesh_info%wb, kmesh_data%num_shells, &
                       param_input%length_unit)
        time2 = io_time()
        write (stdout, '(1x,a25,f11.3,a)') 'Time to get kmesh        ', time2 - time1, ' (sec)'
      endif

      if (driver%postproc_setup) then
        call kmesh_write(recip_lattice, param_input%timing_level, &
                         kmesh_info%nncell, kmesh_info%nnlist, kmesh_info%nntot, &
                         num_kpts, kmesh_data%input_proj%l, num_proj, &
                         kmesh_data%input_proj_site, param_input%spinors, &
                         k_points%kpt_latt, real_lattice, pp_calc%only_A, &
                         kmesh_data%input_proj%zona, kmesh_data%input_proj%x, &
                         kmesh_data%input_proj%z, kmesh_data%input_proj%radial, &
                         kmesh_data%input_proj%m, param_input%exclude_bands, &
                         param_input%num_exclude_bands, &
                         kmesh_data%auto_projections, &
                         kmesh_data%input_proj%s_qaxis, kmesh_data%input_proj%s)

        time1 = io_time()
        write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', time1 - time2, ' (sec)'
        write (stdout, '(/a)') ' '//trim(seedname)//'.nnkp written.'

        call kmesh_dealloc(kmesh_info%nncell, kmesh_info%neigh, kmesh_info%nnlist, &
                           kmesh_info%bk, kmesh_info%bka, kmesh_info%wb)

        call param_dealloc(driver, param_input, param_plot, param_wannierise, &
                           wann_data, kmesh_data, k_points, dis_data, fermi, &
                           atoms, eigval, spec_points, dos_data, berry)

        write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
        stop !JJ seems a bit abrupt?
      endif

    endif ! on_root

    ! We now distribute the parameters to the other nodes
    call param_dist(driver, w90_calcs, pp_calc, param_input, param_plot, &
                    param_wannierise, lsitesymmetry, symmetrize_eps, wann_data, &
                    param_hamil, kmesh_data, kmesh_info, k_points, num_kpts, &
                    dis_data, fermi_surface_data, fermi, tran, atoms, num_bands, &
                    num_wann, eigval, mp_grid, num_proj, real_lattice, &
                    recip_lattice, pw90_calcs, postw90_oper, pw90_common, &
                    pw90_spin, pw90_ham, kpath, kslice, dos_data, berry, &
                    spin_hall, gyrotropic, geninterp, boltz, eig_found)

    if (param_input%gamma_only .and. num_nodes > 1) &
      call io_error('Gamma point branch is serial only at the moment')

  end subroutine wannier_setup

  subroutine wannier_run(seed__name, mp_grid_loc, num_kpts_loc, &
                         real_lattice_loc, recip_lattice_loc, kpt_latt_loc, &
                         num_bands_loc, num_wann_loc, nntot_loc, num_atoms_loc, &
                         atom_symbols_loc, atoms_cart_loc, gamma_only_loc, &
                         m_matrix_loc, a_matrix_loc, eigenvalues_loc, &
                         u_matrix_loc, u_matrix_opt_loc, lwindow_loc, &
                         wann_centres_loc, wann_spreads_loc, spread_loc)

    !! For more information, check a (minimal) example of how it can be used
    !! in the folder test-suite/library-mode-test/test_library.F90

    use w90_constants
    use w90_parameters
    use w90_io
    use w90_hamiltonian
    use w90_kmesh
    use w90_disentangle
    use w90_overlap
    use w90_wannierise
    use w90_plot
    use w90_transport
    use w90_comms, only: my_node_id, num_nodes, comms_array_split, comms_scatterv, &
      on_root

    use w90_param_methods, only: param_lib_set_atoms, param_read, param_write, &
      param_write_chkpt, param_dealloc

    implicit none

    ! passed variables
    integer, intent(in) :: nntot_loc
    integer, intent(in) :: num_atoms_loc
    integer, intent(in) :: num_bands_loc
    integer, intent(in) :: num_kpts_loc
    integer, intent(in) :: num_wann_loc
    integer, dimension(3), intent(in) :: mp_grid_loc

    logical, intent(in) :: gamma_only_loc
    logical, dimension(num_bands_loc, num_kpts_loc), optional, intent(out) :: lwindow_loc !JJ optional

    real(kind=dp), dimension(3, 3), intent(in) :: real_lattice_loc
    real(kind=dp), dimension(3, 3), intent(in) :: recip_lattice_loc
    real(kind=dp), dimension(3, num_atoms_loc), intent(in) :: atoms_cart_loc
    real(kind=dp), dimension(3, num_kpts_loc), intent(in) :: kpt_latt_loc
    real(kind=dp), dimension(3, num_wann_loc), optional, intent(out) :: wann_centres_loc !JJ optional
    real(kind=dp), dimension(3), optional, intent(out) :: spread_loc !JJ optional
    real(kind=dp), dimension(num_bands_loc, num_kpts_loc), intent(in) :: eigenvalues_loc
    real(kind=dp), dimension(num_wann_loc), optional, intent(out) :: wann_spreads_loc !JJ optional

    complex(kind=dp), dimension(num_bands_loc, num_bands_loc, nntot_loc, num_kpts_loc), intent(in) :: m_matrix_loc
    complex(kind=dp), dimension(num_bands_loc, num_wann_loc, num_kpts_loc), intent(in) :: a_matrix_loc
    complex(kind=dp), dimension(num_bands_loc, num_wann_loc, num_kpts_loc), optional, intent(out) :: u_matrix_opt_loc !JJ optional
    complex(kind=dp), dimension(num_wann_loc, num_wann_loc, num_kpts_loc), intent(out) :: u_matrix_loc

    character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
    character(len=*), intent(in) :: seed__name

    ! local variables
    integer :: ierr, loop_k, loop_w
    integer :: nkp, nn, n, m
    integer :: nrpts
    integer :: rpt_origin
    integer, allocatable :: irvec(:, :)
    integer, allocatable :: ndegen(:)
    integer, allocatable :: shift_vec(:, :)
    integer, dimension(0:num_nodes - 1) :: counts ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: displs ! Needed to split an array on different nodes

    logical :: wout_found

    real(kind=dp), allocatable :: wannier_centres_translated(:, :)
    real(kind=dp) time0, time1, time2

    complex(kind=dp), allocatable :: ham_k(:, :, :)
    complex(kind=dp), allocatable :: ham_r(:, :, :)

    character(len=9) :: stat, pos

    type(ham_logical) :: hmlg

    time0 = io_time()

    seedname = trim(adjustl(seed__name))
    inquire (file=trim(seedname)//'.wout', exist=wout_found)
    if (wout_found) then
      stat = 'old'
    else
      stat = 'replace'
    endif
    pos = 'append'

    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.wout', status=trim(stat), position=trim(pos))

    ! copy local data into module variables
    ! if in library mode, then override the data obtained with param_read()
    if (driver%library) then
      num_bands = num_bands_loc
      mp_grid = mp_grid_loc
      num_kpts = num_kpts_loc
      real_lattice = real_lattice_loc
      recip_lattice = recip_lattice_loc
      allocate (k_points%kpt_latt(3, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpt_latt in wannier_setup')
      k_points%kpt_latt = kpt_latt_loc
      allocate (eigval(num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eigval in wannier_setup')
      eigval = eigenvalues_loc
      atoms%num_atoms = num_atoms_loc
      param_input%gamma_only = gamma_only_loc
      call param_lib_set_atoms(atoms, atom_symbols_loc, atoms_cart_loc, recip_lattice)

      call comms_array_split(num_kpts, counts, displs)
      call overlap_allocate(u_matrix, m_matrix_local, m_matrix, u_matrix_opt, a_matrix, m_matrix_orig_local, &
                            m_matrix_orig, param_input%timing_level, kmesh_info%nntot, num_kpts, num_wann, num_bands, &
                            w90_calcs%disentanglement)

      if (w90_calcs%disentanglement) then
        m_matrix_orig = m_matrix_loc
        a_matrix = a_matrix_loc
        u_matrix_opt = cmplx_0
        u_matrix = cmplx_0
      else
        m_matrix = m_matrix_loc
        u_matrix = a_matrix_loc
      endif

      ! IMPORTANT NOTE: _loc are variables local to this function, passed in as variables
      ! Instead, _local are variables local to the MPI process.
      if (w90_calcs%disentanglement) then
        call comms_scatterv(m_matrix_orig_local, num_bands*num_bands*kmesh_info%nntot*counts(my_node_id), &
                            m_matrix_orig, num_bands*num_bands*kmesh_info%nntot*counts, num_bands*num_bands*kmesh_info%nntot*displs)
      else
        call comms_scatterv(m_matrix_local, num_wann*num_wann*kmesh_info%nntot*counts(my_node_id), &
                            m_matrix, num_wann*num_wann*kmesh_info%nntot*counts, num_wann*num_wann*kmesh_info%nntot*displs)
      endif
    endif

! JJ these lines are sadly orphaned
!~  ! Check Mmn(k,b) is symmetric in m and n for gamma_only case
!~  if (gamma_only) call overlap_check_m_symmetry()

    if (w90_calcs%disentanglement) then
      time1 = io_time()
      param_input%have_disentangled = .false.
      call dis_main(num_kpts, kmesh_info%nntot, num_wann, num_bands, &
                    dis_data%spheres_num, dis_data%num_iter, dis_data%spheres_first_wann, &
                    dis_data%conv_window, param_input%timing_level, num_nodes, my_node_id, &
                    param_input%optimisation, param_input%iprint, kmesh_info%nnlist, &
                    dis_data%ndimwin, dis_data%win_min, dis_data%win_max, dis_data%froz_min, &
                    dis_data%froz_max, dis_data%mix_ratio, dis_data%conv_tol, kmesh_info%wbtot, &
                    param_input%lenconfac, param_input%omega_invariant, eigval, recip_lattice, &
                    k_points%kpt_latt, dis_data%spheres, kmesh_info%wb, param_input%devel_flag, &
                    param_input%length_unit, lsitesymmetry, param_input%gamma_only, on_root, &
                    dis_data%frozen_states, dis_data%lwindow, u_matrix, u_matrix_opt, m_matrix, &
                    m_matrix_local, m_matrix_orig, m_matrix_orig_local, a_matrix, sym)

      param_input%have_disentangled = .true.

      call param_write_chkpt('postdis', param_input, wann_data, kmesh_info, &
                             k_points, num_kpts, dis_data, num_bands, &
                             num_wann, u_matrix, u_matrix_opt, m_matrix, &
                             mp_grid, real_lattice, recip_lattice)
      time2 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, ' (sec)'
    else
! JJ, when is this intended?  What does it do??
!    if (param_input%gamma_only) then
!      call overlap_project_gamma(kmesh_info%nntot, m_matrix, u_matrix, param_input%timing_level, num_wann)  !lp note this not called by wannier_prog.F90
!    else
!      call overlap_project(m_matrix_local, kmesh_info%nnlist, kmesh_info%nntot, m_matrix, u_matrix, &
!                           param_input%timing_level, num_kpts, num_wann, num_bands, lsitesymmetry, sym) !lp note this not called by wannier_prog.F90
!    endif
!    time1 = io_time()
!    write (stdout, '(1x,a25,f11.3,a)') 'Time to project overlaps ', time1 - time2, ' (sec)'
    end if

    !there could be a w90_calcs%wannierise then these calls would be more consistent
    time1 = io_time()
    if (param_input%gamma_only) then
      call wann_main_gamma(num_wann, param_wannierise, kmesh_info, &
                           param_input, u_matrix, m_matrix, num_kpts, &
                           real_lattice, wann_data, num_bands, u_matrix_opt, &
                           eigval, dis_data%lwindow, recip_lattice, atoms, &
                           k_points, dis_data, mp_grid, stdout)
    else
      call wann_main(num_wann, param_wannierise, kmesh_info, param_input, &
                     u_matrix, m_matrix, num_kpts, real_lattice, num_proj, &
                     wann_data, k_points, num_bands, u_matrix_opt, &
                     eigval, dis_data, recip_lattice, atoms, &
                     lsitesymmetry, stdout, mp_grid, w90_calcs, &
                     tran%mode, param_hamil, sym, ham_r, irvec, &
                     shift_vec, ndegen, nrpts, rpt_origin, &
                     wannier_centres_translated, hmlg, ham_k)
    endif
    time2 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time2 - time1, ' (sec)'

    call param_write_chkpt('postwann', param_input, wann_data, kmesh_info, &
                           k_points, num_kpts, dis_data, num_bands, &
                           num_wann, u_matrix, u_matrix_opt, m_matrix, &
                           mp_grid, real_lattice, recip_lattice)

    ! fixme, write_u_matrices should trigger wannier_plot by some other means than this test here JJ
    ! fixme,  also write_bvec
    if (w90_calcs%wannier_plot .or. w90_calcs%bands_plot .or. w90_calcs%fermi_surface_plot &
        .or. w90_calcs%write_hr .or. param_plot%write_u_matrices .or. param_plot%write_bvec) then
      call plot_main(num_kpts, w90_calcs%bands_plot, dos_plot, k_points%kpt_latt, &
                     w90_calcs%fermi_surface_plot, w90_calcs%wannier_plot, &
                     param_input%timing_level, param_plot%write_bvec, w90_calcs%write_hr, &
                     param_plot%write_rmn, param_plot%write_tb, param_plot%write_u_matrices, &
                     real_lattice, num_wann, kmesh_info%wb, kmesh_info%bk, m_matrix, &
                     kmesh_info%nntot, recip_lattice, wann_data%centres, atoms%num_atoms, &
                     atoms%pos_cart, param_hamil%translation_centre_frac, &
                     param_hamil%automatic_translation, atoms%num_species, atoms%species_num, &
                     param_input%lenconfac, param_input%have_disentangled, dis_data%ndimwin, &
                     dis_data%lwindow, u_matrix_opt, eigval, u_matrix, lsitesymmetry, &
                     num_bands, param_input%ws_distance_tol, param_input%ws_search_size, &
                     mp_grid, tran%mode, param_input%bands_plot_mode, &
                     w90_calcs%transport, param_input%iprint, param_plot%wannier_plot_radius, &
                     param_plot%wannier_plot_scale, atoms%pos_frac, &
                     param_plot%wannier_plot_spinor_phase, param_plot%wannier_plot_spinor_mode, &
                     param_input%spinors, param_plot%wannier_plot_format, &
                     param_plot%wvfn_formatted, param_plot%wannier_plot_mode, &
                     param_plot%wannier_plot_list, param_plot%num_wannier_plot, atoms%symbol, &
                     param_plot%spin, param_plot%wannier_plot_supercell, fermi%energy_list, fermi%n, &
                     fermi_surface_data%num_points, param_input%one_dim_dir, &
                     param_plot%bands_plot_dim, param_input%hr_cutoff, param_input%dist_cutoff, &
                     param_input%dist_cutoff_mode, param_input%use_ws_distance, &
                     param_plot%bands_plot_project, param_plot%num_bands_project, &
                     param_plot%bands_plot_format, spec_points%bands_label, &
                     spec_points%bands_spec_points, spec_points%bands_num_spec_points, &
                     param_plot%bands_num_points, ham_r, irvec, shift_vec, ndegen, nrpts, &
                     rpt_origin, wannier_centres_translated, hmlg, ham_k)

      time1 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
    end if

    time2 = io_time()
    if (w90_calcs%transport) then
      call tran_main(tran%mode, tran%read_ht, param_input%timing_level, &
                     w90_calcs%write_hr, param_input%write_xyz, num_wann, real_lattice, &
                     recip_lattice, wann_data%centres, atoms%num_atoms, w90_calcs%bands_plot, &
                     param_input%iprint, param_hamil%translation_centre_frac, &
                     param_hamil%automatic_translation, atoms%num_species, atoms%species_num, &
                     param_input%lenconfac, param_input%have_disentangled, dis_data%ndimwin, &
                     dis_data%lwindow, u_matrix_opt, k_points%kpt_latt, eigval, u_matrix, &
                     lsitesymmetry, num_bands, num_kpts, atoms%pos_cart, &
                     param_input%ws_distance_tol, param_input%ws_search_size, mp_grid, &
                     param_input%bands_plot_mode, w90_calcs%transport, param_input%dist_cutoff_hc, &
                     param_input%dist_cutoff, param_input%dist_cutoff_mode, tran%num_bandc, &
                     tran%num_cc, tran%num_rr, tran%num_lc, tran%num_cr, tran%write_ht, &
                     fermi%energy_list, fermi%n, k_points%kpt_cart, tran%num_ll, tran%num_cell_ll, &
                     tran%easy_fix, atoms%symbol, wann_data%spreads, tran%group_threshold, &
                     param_input%one_dim_dir, tran%use_same_lead, tran%energy_step, tran%win_min, &
                     tran%win_max, tran%num_bb, param_input%length_unit, param_input%hr_cutoff, &
                     ham_r, irvec, shift_vec, ndegen, nrpts, rpt_origin, &
                     wannier_centres_translated, hmlg, ham_k)

      call tran_dealloc()
      time1 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
    end if

    ! Now we zero all of the local output data, then copy in the data
    ! from the parameters module

    u_matrix_loc = u_matrix
    if (present(u_matrix_opt_loc) .and. present(lwindow_loc)) then
    if (w90_calcs%disentanglement) then
      u_matrix_opt_loc = u_matrix_opt
      lwindow_loc = dis_data%lwindow
    else
      u_matrix_opt_loc = cmplx_0
      do loop_k = 1, num_kpts
        do loop_w = 1, num_wann
          u_matrix_opt_loc(loop_w, loop_w, loop_k) = cmplx_1
        end do
      end do
      lwindow_loc = .true.
    end if
    end if

    if (present(wann_centres_loc)) wann_centres_loc = wann_data%centres
    if (present(wann_spreads_loc)) wann_spreads_loc = wann_data%spreads
    if (present(spread_loc)) then
      spread_loc(1) = param_wannierise%omega_total
      spread_loc(2) = param_input%omega_invariant   !JJ maybe mv omg_inv to param_wann?
      spread_loc(3) = param_wannierise%omega_tilde
    endif
!  call hamiltonian_dealloc(ham_r, irvec, ndegen, wannier_centres_translated, &
!                           hmlg, ham_k)
!  call overlap_dealloc(m_matrix_orig_local, m_matrix_local, u_matrix_opt, &
!                       a_matrix, m_matrix_orig, m_matrix, u_matrix)
!  call kmesh_dealloc(kmesh_info%nncell, kmesh_info%neigh, kmesh_info%nnlist, kmesh_info%bk, kmesh_info%bka, kmesh_info%wb)
!  call param_dealloc(driver, param_input, param_plot, param_wannierise, &
!                     wann_data, kmesh_data, k_points, dis_data, fermi, &
!                     atoms, eigval, spec_points, dos_data, berry)
!
    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time() - time0, ' (sec)'

    if (param_input%timing_level > 0) call io_print_timings()

    write (stdout, *)
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'
    close (stdout)

  end subroutine wannier_run

end module
