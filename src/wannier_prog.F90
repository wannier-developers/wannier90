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

program wannier
  !! The main Wannier90 program

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
  use w90_comms, only: on_root, num_nodes, comms_setup, comms_end, comms_bcast, my_node_id
  use w90_sitesym !YN:

  use w90_param_methods, only: param_read, param_write, param_postw90_write, param_dealloc, &
    param_write_header, param_write_chkpt, param_read_chkpt, param_lib_set_atoms, &
    param_memory_estimate, param_get_smearing_type, param_get_convention_type, &
    param_dist, param_chkpt_dist

  implicit none

! integer, save :: rpt_origin
  !! index of R=0
! integer, save :: nrpts
  !! number of Wigner-Seitz grid points
! integer, save, allocatable :: irvec(:, :)
  !!  The irpt-th Wigner-Seitz grid point has components
  !! irvec(1:3,irpt) in the basis of the lattice vectors
! integer, save, allocatable :: shift_vec(:, :)
! integer, save, allocatable :: ndegen(:)
  !! Weight of the irpt-th point is 1/ndegen(irpt)
! real(kind=dp), save, allocatable :: wannier_centres_translated(:, :)
! complex(kind=dp), save, allocatable :: ham_r(:, :, :)
  !! Hamiltonian matrix in WF representation

  integer :: rpt_origin
  !! index of R=0
  integer :: nrpts
  !! number of Wigner-Seitz grid points
  integer, allocatable :: irvec(:, :)
  !!  The irpt-th Wigner-Seitz grid point has components
  !! irvec(1:3,irpt) in the basis of the lattice vectors
  integer, allocatable :: shift_vec(:, :)
  integer, allocatable :: ndegen(:)
  !! Weight of the irpt-th point is 1/ndegen(irpt)
  real(kind=dp), allocatable :: wannier_centres_translated(:, :)
  complex(kind=dp), allocatable :: ham_r(:, :, :)
  !! Hamiltonian matrix in WF representation

  complex(kind=dp), allocatable :: ham_k(:, :, :)

  real(kind=dp) time0, time1, time2
  character(len=9) :: stat, pos, cdate, ctime
  logical :: wout_found, dryrun
  integer :: len_seedname
  character(len=50) :: prog

  logical :: library !JJ
  type(sitesym_data) :: sym !JJ

! logical :: ham_have_setup = .false.
! logical :: have_translated = .false.
! logical :: use_translation = .false.
  type(ham_logical) :: hmlg

  call comms_setup

  library = .false.

  time0 = io_time()

  if (on_root) then
    prog = 'wannier90'
    call io_commandline(prog, dryrun)
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1)
  call comms_bcast(seedname, len_seedname)
  call comms_bcast(dryrun, 1)

  if (on_root) then
    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.werr')
    call io_date(cdate, ctime)
    write (stdout, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

    call param_read
    close (stdout, status='delete')

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
      write (stdout, '(/,1x,a,i3,a/)') &
        'Running in parallel on ', num_nodes, ' CPUs'
    endif
    call param_write()

    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

    if (.not. driver%explicit_nnkpts) call kmesh_get(recip_lattice, &
                                                     k_points%kpt_cart, param_input%timing_level, kmesh_info%nncell, &
                                                     kmesh_info%neigh, kmesh_info%nnlist, kmesh_info%nntot, kmesh_data%shell_list, &
                                                     param_input%devel_flag, param_input%iprint, param_input%lenconfac, &
                                                     kmesh_data%tol, num_kpts, kmesh_data%search_shells, param_input%gamma_only, &
                                                     kmesh_info%nnh, kmesh_info%wbtot, kmesh_data%skip_B1_tests, kmesh_info%bk, &
                                                     kmesh_info%bka, kmesh_info%wb, kmesh_data%num_shells, param_input%length_unit)
    time2 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') &
      'Time to get kmesh        ', time2 - time1, ' (sec)'

    call param_memory_estimate
  end if

  if (dryrun) then
    if (on_root) then
      write (stdout, *) ' '
      write (stdout, *) '                       ==============================='
      write (stdout, *) '                                   DRYRUN             '
      write (stdout, *) '                       No problems found with win file'
      write (stdout, *) '                       ==============================='
    endif
    stop
  endif

  ! We now distribute the parameters to the other nodes
  call param_dist
  if (param_input%gamma_only .and. num_nodes > 1) &
    call io_error('Gamma point branch is serial only at the moment')

  if (w90_calcs%transport .and. tran%read_ht) goto 3003

  ! Sort out restarts
  if (driver%restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) call param_read_chkpt()
    call param_chkpt_dist
    if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
    if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

    select case (driver%restart)
    case ('default')    ! continue from where last checkpoint was written
      if (on_root) write (stdout, '(/1x,a)', advance='no') 'Resuming a previous Wannier90 calculation '
      if (driver%checkpoint .eq. 'postdis') then
        if (on_root) write (stdout, '(a/)') 'from wannierisation ...'
        goto 1001         ! go to wann_main
      elseif (driver%checkpoint .eq. 'postwann') then
        if (on_root) write (stdout, '(a/)') 'from plotting ...'
        goto 2002         ! go to plot_main
      else
        if (on_root) write (stdout, '(/a/)')
        call io_error('Value of checkpoint not recognised in wann_prog')
      endif
    case ('wannierise') ! continue from wann_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
      goto 1001
    case ('plot')       ! continue from plot_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
      goto 2002
    case ('transport')   ! continue from tran_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from transport routines ...'
      goto 3003
    case default        ! for completeness... (it is already trapped in param_read)
      call io_error('Value of restart not recognised in wann_prog')
    end select
  endif

  if (driver%postproc_setup) then
    if (on_root) call kmesh_write(recip_lattice, param_input%timing_level, &
                                  kmesh_info%nncell, kmesh_info%nnlist, kmesh_info%nntot, num_kpts, &
                                  kmesh_data%input_proj%l, num_proj, kmesh_data%input_proj_site, &
                                  param_input%spinors, k_points%kpt_latt, real_lattice, pp_calc%only_A, &
                                  kmesh_data%input_proj%zona, kmesh_data%input_proj%x, kmesh_data%input_proj%z, &
                                  kmesh_data%input_proj%radial, kmesh_data%input_proj%m, &
                                  param_input%exclude_bands, param_input%num_exclude_bands, &
                                  kmesh_data%auto_projections, kmesh_data%input_proj%s_qaxis, &
                                  kmesh_data%input_proj%s)
    call kmesh_dealloc(kmesh_info%nncell, kmesh_info%neigh, kmesh_info%nnlist, kmesh_info%bk, kmesh_info%bka, kmesh_info%wb)
    call param_dealloc()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (on_root) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
    call comms_end
    stop
  endif

  if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
  if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

  call overlap_allocate(u_matrix, m_matrix_local, m_matrix, u_matrix_opt, &
                        a_matrix, m_matrix_orig_local, m_matrix_orig, param_input%timing_level, &
                        kmesh_info%nntot, num_kpts, num_wann, num_bands, &
                        w90_calcs%disentanglement)
  call overlap_read(lsitesymmetry, m_matrix_orig_local, m_matrix_local, &
                    param_input%gamma_only, w90_calcs%use_bloch_phases, w90_calcs%cp_pp, &
                    u_matrix_opt, m_matrix_orig, param_input%timing_level, a_matrix, m_matrix, &
                    u_matrix, select_proj%proj2wann_map, &
                    select_proj%lselproj, num_proj, kmesh_info%nnlist, kmesh_info%nncell, &
                    kmesh_info%nntot, num_kpts, num_wann, num_bands, &
                    w90_calcs%disentanglement, sym)

  time1 = io_time()
  if (on_root) write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, ' (sec)'

  param_input%have_disentangled = .false.

  if (w90_calcs%disentanglement) then
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
    time2 = io_time()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, ' (sec)'
  endif

  if (on_root) call param_write_chkpt('postdis')
!~  call param_write_um

1001 time2 = io_time()

  if (.not. param_input%gamma_only) then
    call wann_main(num_wann, param_wannierise, kmesh_info, param_input, &
                   u_matrix, m_matrix, num_kpts, real_lattice, num_proj, &
                   wann_data, k_points%kpt_latt, num_bands, u_matrix_opt, &
                   eigval, dis_data, recip_lattice, atoms, &
                   lsitesymmetry, stdout, mp_grid, w90_calcs, &
                   tran%mode, param_hamil, sym, ham_r, irvec, &
                   shift_vec, ndegen, nrpts, rpt_origin, &
                   wannier_centres_translated, hmlg, ham_k)
  else
    call wann_main_gamma(num_wann, param_wannierise, kmesh_info, &
                         param_input, u_matrix, m_matrix, num_kpts, &
                         real_lattice, wann_data, num_bands, u_matrix_opt, &
                         eigval, dis_data%lwindow, recip_lattice, atoms, stdout)
  end if

  time1 = io_time()
  if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time1 - time2, ' (sec)'

  if (on_root) call param_write_chkpt('postwann')

2002 continue
  if (on_root) then
    ! I call the routine always; the if statements to decide if/what
    ! to plot are inside the function
    time2 = io_time()
    !
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
    !
    time1 = io_time()
    ! Now time is always printed, even if no plotting is done/required, but
    ! it shouldn't be a problem.
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
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
      time1 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (tran%read_ht) goto 4004
    end if
  endif

  call tran_dealloc()
  call hamiltonian_dealloc(ham_r, irvec, ndegen, wannier_centres_translated, &
                           hmlg, ham_k)
  call overlap_dealloc(m_matrix_orig_local, m_matrix_local, u_matrix_opt, &
                       a_matrix, m_matrix_orig, m_matrix, u_matrix)
  call kmesh_dealloc(kmesh_info%nncell, kmesh_info%neigh, kmesh_info%nnlist, kmesh_info%bk, kmesh_info%bka, kmesh_info%wb)
  call param_dealloc()
  if (lsitesymmetry) call sitesym_dealloc(sym) !YN:

4004 continue

  if (on_root) then
    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time(), ' (sec)'

    if (param_input%timing_level > 0) call io_print_timings()

    write (stdout, *)
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'

    close (stdout)
  endif

  call comms_end

end program wannier

