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

    if (restart .eq. ' ') then
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

    if (.not. explicit_nnkpts) call kmesh_get(recip_lattice, kpt_cart, timing_level, nncell, neigh, &
                                              nnlist, nntot, shell_list, devel_flag, iprint, lenconfac, &
                                              kmesh_tol, num_kpts, search_shells, gamma_only, nnh, wbtot, &
                                              skip_B1_tests, bk, bka, wb, num_shells, length_unit)
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
  if (gamma_only .and. num_nodes > 1) &
    call io_error('Gamma point branch is serial only at the moment')

  if (transport .and. tran_read_ht) goto 3003

  ! Sort out restarts
  if (restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) call param_read_chkpt()
    call param_chkpt_dist
    if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
    if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

    select case (restart)
    case ('default')    ! continue from where last checkpoint was written
      if (on_root) write (stdout, '(/1x,a)', advance='no') 'Resuming a previous Wannier90 calculation '
      if (checkpoint .eq. 'postdis') then
        if (on_root) write (stdout, '(a/)') 'from wannierisation ...'
        goto 1001         ! go to wann_main
      elseif (checkpoint .eq. 'postwann') then
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

  if (postproc_setup) then
    if (on_root) call kmesh_write(recip_lattice, timing_level, nncell, nnlist, nntot, num_kpts, &
                                  input_proj_l, num_proj, input_proj_site, spinors, kpt_latt, real_lattice, &
                                  calc_only_A, input_proj_zona, input_proj_x, input_proj_z, input_proj_radial, &
                                  input_proj_m, exclude_bands, num_exclude_bands, auto_projections, &
                                  input_proj_s_qaxis, input_proj_s)
    call kmesh_dealloc(nncell, neigh, nnlist, bk, bka, wb)
    call param_dealloc()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (on_root) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
    call comms_end
    stop
  endif

  if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
  if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

  call overlap_allocate(u_matrix, m_matrix_local, m_matrix, u_matrix_opt, a_matrix, m_matrix_orig_local, &
                        m_matrix_orig, timing_level, nntot, num_kpts, num_wann, num_bands, disentanglement)
  call overlap_read(lsitesymmetry, m_matrix_orig_local, m_matrix_local, gamma_only, use_bloch_phases, &
                    cp_pp, u_matrix_opt, m_matrix_orig, timing_level, a_matrix, m_matrix, u_matrix, &
                    devel_flag, proj2wann_map, lselproj, num_proj, nnlist, nncell, nntot, num_kpts, &
                    num_wann, num_bands, disentanglement, sym)

  time1 = io_time()
  if (on_root) write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, ' (sec)'

  have_disentangled = .false.

  if (disentanglement) then
    if (.not. allocated(dis_spheres)) then ! this is *sometimes* true for non-root mpi tasks (depending on input options)
      allocate (dis_spheres(1, 1)) !JJ temporary workaround to avoid runtime check failure
    endif

    !call dis_main()
    call dis_main(num_kpts, nntot, num_wann, num_bands, dis_spheres_num, &
                  dis_num_iter, dis_spheres_first_wann, dis_conv_window, timing_level, &
                  num_nodes, my_node_id, optimisation, iprint, nnlist, ndimwin, dis_win_min, &
                  dis_win_max, dis_froz_min, dis_froz_max, dis_mix_ratio, dis_conv_tol, &
                  wbtot, lenconfac, omega_invariant, eigval, recip_lattice, kpt_latt, &
                  dis_spheres, wb, devel_flag, length_unit, lsitesymmetry, gamma_only, &
                  on_root, frozen_states, lwindow, u_matrix, u_matrix_opt, m_matrix, &
                  m_matrix_local, m_matrix_orig, m_matrix_orig_local, a_matrix, sym)

    have_disentangled = .true.
    time2 = io_time()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, ' (sec)'
  endif

  if (on_root) call param_write_chkpt('postdis')
!~  call param_write_um

1001 time2 = io_time()

  ! JJ hack to workaround mpi_scatterv requirement that all arrays are valid *for all mpi procs*
  ! m_matrix* usually alloc'd in overlaps.F90, but not always
  if (.not. allocated(m_matrix)) allocate (m_matrix(1, 1, 1, 1)) !JJ temporary workaround to avoid runtime check failure

  if (.not. gamma_only) then
    call wann_main(num_wann, num_cg_steps, num_iter, nnlist, nntot, &
                   wbtot, u_matrix, m_matrix, num_kpts, iprint, &
                   num_print_cycles, num_dump_cycles, omega_invariant, &
                   length_unit, lenconfac, proj_site, real_lattice, &
                   write_r2mn, guiding_centres, num_guide_cycles, &
                   num_no_guide_iter, timing_level, trial_step, precond, &
                   fixed_step, lfixstep, write_proj, have_disentangled, &
                   conv_tol, num_proj, conv_window, conv_noise_amp, &
                   conv_noise_num, wannier_centres, write_xyz, &
                   wannier_spreads, omega_total, omega_tilde, &
                   optimisation, write_vdw_data, write_hr_diag, kpt_latt, &
                   bk, ccentres_cart, slwf_num, selective_loc, &
                   slwf_constrain, slwf_lambda, neigh, nnh, bka, &
                   num_bands, u_matrix_opt, eigval, lwindow, wb, &
                   translate_home_cell, recip_lattice, num_atoms, &
                   atoms_symbol, atoms_pos_cart, num_species, &
                   atoms_species_num, num_valence_bands, &
                   num_elec_per_state, lsitesymmetry, stdout, &
                   ws_distance_tol, ws_search_size, real_metric, mp_grid, &
                   transport_mode, bands_plot_mode, transport, bands_plot, &
                   translation_centre_frac, automatic_translation, ndimwin, sym, &
                   ham_r, irvec, shift_vec, ndegen, nrpts, rpt_origin, &
                   wannier_centres_translated, &
                   hmlg, ham_k)

  else
    call wann_main_gamma(num_wann, num_iter, wb, nntot, u_matrix, m_matrix, &
                         num_kpts, iprint, num_print_cycles, &
                         num_dump_cycles, omega_invariant, length_unit, &
                         lenconfac, proj_site, real_lattice, write_r2mn, &
                         guiding_centres, num_guide_cycles, &
                         num_no_guide_iter, timing_level, write_proj, &
                         have_disentangled, conv_tol, conv_window, &
                         wannier_centres, write_xyz, wannier_spreads, &
                         omega_total, omega_tilde, write_vdw_data, neigh, &
                         nnh, bk, bka, num_bands, u_matrix_opt, eigval, &
                         lwindow, wbtot, translate_home_cell, &
                         recip_lattice, num_atoms, atoms_symbol, &
                         atoms_pos_cart, num_species, atoms_species_num, &
                         num_valence_bands, num_elec_per_state, stdout)
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
    call plot_main(num_kpts, bands_plot, dos_plot, kpt_latt, fermi_surface_plot, wannier_plot, &
                   timing_level, write_bvec, write_hr, write_rmn, write_tb, write_u_matrices, &
                   real_lattice, num_wann, wb, bk, m_matrix, nntot, recip_lattice, wannier_centres, &
                   num_atoms, atoms_pos_cart, translation_centre_frac, automatic_translation, &
                   num_species, atoms_species_num, lenconfac, have_disentangled, ndimwin, lwindow, &
                   u_matrix_opt, eigval, u_matrix, lsitesymmetry, num_bands, ws_distance_tol, &
                   ws_search_size, real_metric, mp_grid, transport_mode, bands_plot_mode, transport, &
                   iprint, wannier_plot_radius, wannier_plot_scale, atoms_pos_frac, &
                   wannier_plot_spinor_phase, wannier_plot_spinor_mode, spinors, wannier_plot_format, &
                   wvfn_formatted, wannier_plot_mode, wannier_plot_list, num_wannier_plot, &
                   atoms_symbol, spin, wannier_plot_supercell, fermi_energy_list, nfermi, &
                   fermi_surface_num_points, one_dim_dir, bands_plot_dim, hr_cutoff, dist_cutoff, &
                   dist_cutoff_mode, use_ws_distance, bands_plot_project, num_bands_project, &
                   bands_plot_format, bands_label, bands_spec_points, bands_num_spec_points, &
                   recip_metric, bands_num_points, ham_r, irvec, shift_vec, ndegen, nrpts, rpt_origin, &
                   wannier_centres_translated, hmlg, ham_k)
    !
    time1 = io_time()
    ! Now time is always printed, even if no plotting is done/required, but
    ! it shouldn't be a problem.
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
    time2 = io_time()
    if (transport) then
      call tran_main(transport_mode, tran_read_ht, timing_level, write_hr, write_xyz, num_wann, &
                     real_lattice, recip_lattice, wannier_centres, num_atoms, bands_plot, iprint, &
                     translation_centre_frac, automatic_translation, num_species, atoms_species_num, &
                     lenconfac, have_disentangled, ndimwin, lwindow, u_matrix_opt, kpt_latt, &
                     eigval, u_matrix, lsitesymmetry, num_bands, num_kpts, atoms_pos_cart, &
                     ws_distance_tol, ws_search_size, real_metric, mp_grid, bands_plot_mode, transport, &
                     dist_cutoff_hc, dist_cutoff, dist_cutoff_mode, tran_num_bandc, tran_num_cc, &
                     tran_num_rr, tran_num_lc, tran_num_cr, tran_write_ht, fermi_energy_list, nfermi, &
                     kpt_cart, tran_num_ll, tran_num_cell_ll, tran_easy_fix, atoms_symbol, &
                     wannier_spreads, tran_group_threshold, one_dim_dir, tran_use_same_lead, &
                     tran_energy_step, tran_win_min, tran_win_max, tran_num_bb, length_unit, hr_cutoff, &
                     ham_r, irvec, shift_vec, ndegen, nrpts, rpt_origin, wannier_centres_translated, &
                     hmlg, ham_k)
      time1 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (tran_read_ht) goto 4004
    end if
  endif

  call tran_dealloc()
  call hamiltonian_dealloc(ham_r, irvec, ndegen, wannier_centres_translated, &
                           hmlg, ham_k)
  call overlap_dealloc(m_matrix_orig_local, m_matrix_local, u_matrix_opt, &
                       a_matrix, m_matrix_orig, m_matrix, u_matrix)
  call kmesh_dealloc(nncell, neigh, nnlist, bk, bka, wb)
  call param_dealloc()
  if (lsitesymmetry) call sitesym_dealloc(sym) !YN:

4004 continue

  if (on_root) then
    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time(), ' (sec)'

    if (timing_level > 0) call io_print_timings()

    write (stdout, *)
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'

    close (stdout)
  endif

  call comms_end

end program wannier

