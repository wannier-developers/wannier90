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
                         real_lattice_loc, recip_lattice_loc, kpt_latt_loc, num_bands_tot, &
                         num_atoms_loc, atom_symbols_loc, atoms_cart_loc, gamma_only_loc, spinors_loc, &
                         nntot_loc, nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc, &
                         proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc, proj_z_loc, &
                         proj_x_loc, proj_zona_loc, exclude_bands_loc, proj_s_loc, proj_s_qaxis_loc)

  !! This routine should be called first from a code calling the library
  !! mode to setup all the variables.
  !! NOTE! The library mode currently works ONLY in serial (when called from
  !! a parallel code, make sure to run it only on 1 MPI process)
  !!
  !! For more information, check a (minimal) example of how it can be used
  !! in the folder test-suite/library-mode-test/test_library.F90

  use w90_constants
  use w90_parameters
  use w90_sitesym
  use w90_io
  use w90_kmesh
  use w90_comms, only: comms_setup_vars

  implicit none

  character(len=*), intent(in) :: seed__name
  integer, dimension(3), intent(in) :: mp_grid_loc
  integer, intent(in) :: num_kpts_loc
  real(kind=dp), dimension(3, 3), intent(in) :: real_lattice_loc
  real(kind=dp), dimension(3, 3), intent(in) :: recip_lattice_loc
  real(kind=dp), dimension(3, num_kpts_loc), intent(in) :: kpt_latt_loc
  integer, intent(in) :: num_bands_tot
  integer, intent(in) :: num_atoms_loc
  character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
  real(kind=dp), dimension(3, num_atoms_loc), intent(in) :: atoms_cart_loc
  logical, intent(in) :: gamma_only_loc
  logical, intent(in) :: spinors_loc
  integer, intent(out) :: nntot_loc
  integer, dimension(num_kpts_loc, num_nnmax), intent(out) :: nnlist_loc
  integer, dimension(3, num_kpts_loc, num_nnmax), intent(out) :: nncell_loc
  integer, intent(out) :: num_bands_loc
  integer, intent(out) :: num_wann_loc
  real(kind=dp), dimension(3, num_bands_tot), intent(out) :: proj_site_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_l_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_m_loc
  integer, dimension(num_bands_tot), intent(out) :: proj_radial_loc
  real(kind=dp), dimension(3, num_bands_tot), intent(out) :: proj_z_loc
  real(kind=dp), dimension(3, num_bands_tot), intent(out) :: proj_x_loc
  real(kind=dp), dimension(num_bands_tot), intent(out) :: proj_zona_loc
  integer, dimension(num_bands_tot), intent(out) :: exclude_bands_loc
  integer, dimension(num_bands_tot), optional, intent(out) :: proj_s_loc
  real(kind=dp), dimension(3, num_bands_tot), optional, intent(out) :: proj_s_qaxis_loc

  real(kind=dp) time0, time1, time2
  character(len=9) :: stat, pos, cdate, ctime
  integer :: ierr
  logical :: wout_found

  time0 = io_time()

  call comms_setup_vars

  library = .true.
!  seedname="wannier"
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

  call param_write_header()

  write (stdout, '(/a/)') ' Wannier90 is running in LIBRARY MODE'
  write (stdout, '(a/)') ' Setting up k-point neighbours...'

  ! copy local data into module variables
  mp_grid = mp_grid_loc
  num_kpts = num_kpts_loc
  real_lattice = real_lattice_loc
  recip_lattice = recip_lattice_loc
  allocate (kpt_latt(3, num_kpts), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating kpt_latt in wannier_setup')
  kpt_latt = kpt_latt_loc
  num_atoms = num_atoms_loc
  call param_lib_set_atoms(atom_symbols_loc, atoms_cart_loc)
  gamma_only = gamma_only_loc
  spinors = spinors_loc

  ! GP: at this point we don't know yet the number of excluded bands...
  num_bands = num_bands_tot
  library_param_read_first_pass = .true.
  call param_read()
  ! Following calls will all NOT be first_pass, and I need to pass
  ! directly num_bands, that is already set internally now to num_bands = num_bands_tot - num_exclude_bands
  library_param_read_first_pass = .false.
  ! set cell_volume as it is written to output in param_write
  cell_volume = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - real_lattice(3, 2)*real_lattice(2, 3)) + &
                real_lattice(1, 2)*(real_lattice(2, 3)*real_lattice(3, 1) - real_lattice(3, 3)*real_lattice(2, 1)) + &
                real_lattice(1, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - real_lattice(3, 1)*real_lattice(2, 2))
  call param_write()

  time1 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

  if (.not. explicit_nnkpts) call kmesh_get(recip_lattice, kpt_cart, timing_level, nncell, neigh, &
                                            nnlist, nntot, shell_list, devel_flag, iprint, lenconfac, &
                                            kmesh_tol, num_kpts, search_shells, gamma_only, nnh, wbtot, &
                                            skip_B1_tests, bk, bka, wb, num_shells, length_unit)

  ! Now we zero all of the local output data, then copy in the data
  ! from the parameters module

  nntot_loc = 0
  nnlist_loc = 0
  nncell_loc = 0
  proj_site_loc = 0.0_dp
  proj_l_loc = 0
  proj_m_loc = 0
  proj_z_loc = 0.0_dp
  proj_x_loc = 0.0_dp
  proj_radial_loc = 0
  proj_zona_loc = 0.0_dp
  exclude_bands_loc = 0

  nntot_loc = nntot
  nnlist_loc(:, 1:nntot) = nnlist(:, 1:nntot)
  nncell_loc(:, :, 1:nntot) = nncell(:, :, 1:nntot)
  num_bands_loc = num_bands
  num_wann_loc = num_wann
  if (allocated(proj_site)) then
    proj_site_loc(:, 1:num_proj) = proj_site(:, 1:num_proj)
    proj_l_loc(1:num_proj) = proj_l(1:num_proj)
    proj_m_loc(1:num_proj) = proj_m(1:num_proj)
    proj_z_loc(:, 1:num_proj) = proj_z(:, 1:num_proj)
    proj_x_loc(:, 1:num_proj) = proj_x(:, 1:num_proj)
    proj_radial_loc(1:num_proj) = proj_radial(1:num_proj)
    proj_zona_loc(1:num_proj) = proj_zona(1:num_proj)
    if (allocated(proj_s) .and. present(proj_s_loc) .and. present(proj_s_qaxis_loc)) then
      proj_s_loc(1:num_proj) = proj_s(1:num_proj)
      proj_s_qaxis_loc(:, 1:num_proj) = proj_s_qaxis(:, 1:num_proj)
    end if
  endif
  if (allocated(exclude_bands)) then
    exclude_bands_loc(1:num_exclude_bands) = exclude_bands(1:num_exclude_bands)
  end if

  if (postproc_setup) then
    call kmesh_write(recip_lattice, timing_level, nncell, nnlist, nntot, num_kpts, &
                     input_proj_l, num_proj, input_proj_site, spinors, kpt_latt, real_lattice, &
                     calc_only_A, input_proj_zona, input_proj_x, input_proj_z, input_proj_radial, &
                     input_proj_m, exclude_bands, num_exclude_bands, auto_projections, &
                     input_proj_s_qaxis, input_proj_s)
    write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    write (stdout, '(/a)') ' '//trim(seedname)//'.nnkp written.'
  endif

  call kmesh_dealloc(nncell, neigh, nnlist, bk, bka, wb)
  call param_dealloc()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'

  write (stdout, '(/a/)') ' Finished setting up k-point neighbours.'

  call io_date(cdate, ctime)

  write (stdout, '(2a)') ' Exiting wannier_setup in wannier90 ', ctime

  close (stdout)

end subroutine wannier_setup

subroutine wannier_run(seed__name, mp_grid_loc, num_kpts_loc, &
                       real_lattice_loc, recip_lattice_loc, kpt_latt_loc, num_bands_loc, &
                       num_wann_loc, nntot_loc, num_atoms_loc, atom_symbols_loc, &
                       atoms_cart_loc, gamma_only_loc, M_matrix_loc, A_matrix_loc, eigenvalues_loc, &
                       U_matrix_loc, U_matrix_opt_loc, lwindow_loc, wann_centres_loc, &
                       wann_spreads_loc, spread_loc)

  !! This routine should be called after wannier_setup from a code calling
  !! the library mode to actually run the Wannier code.
  !!
  !! NOTE! The library mode currently works ONLY in serial.
  !! When called from an external code, wannier90 needs to be compiled
  !! in sequential and wannier_run called with 1 MPI process.
  !!
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
  use w90_comms, only: my_node_id, num_nodes, &
    comms_array_split, comms_scatterv, on_root

  implicit none

  character(len=*), intent(in) :: seed__name
  integer, dimension(3), intent(in) :: mp_grid_loc
  integer, intent(in) :: num_kpts_loc
  real(kind=dp), dimension(3, 3), intent(in) :: real_lattice_loc
  real(kind=dp), dimension(3, 3), intent(in) :: recip_lattice_loc
  real(kind=dp), dimension(3, num_kpts_loc), intent(in) :: kpt_latt_loc
  integer, intent(in) :: num_bands_loc
  integer, intent(in) :: num_wann_loc
  integer, intent(in) :: nntot_loc
  integer, intent(in) :: num_atoms_loc
  character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
  real(kind=dp), dimension(3, num_atoms_loc), intent(in) :: atoms_cart_loc
  logical, intent(in) :: gamma_only_loc
  complex(kind=dp), dimension(num_bands_loc, num_bands_loc, nntot_loc, num_kpts_loc), intent(in) :: M_matrix_loc
  complex(kind=dp), dimension(num_bands_loc, num_wann_loc, num_kpts_loc), intent(in) :: A_matrix_loc
  real(kind=dp), dimension(num_bands_loc, num_kpts_loc), intent(in) :: eigenvalues_loc
  complex(kind=dp), dimension(num_wann_loc, num_wann_loc, num_kpts_loc), intent(out) :: U_matrix_loc
  complex(kind=dp), dimension(num_bands_loc, num_wann_loc, num_kpts_loc), optional, intent(out) :: U_matrix_opt_loc
  logical, dimension(num_bands_loc, num_kpts_loc), optional, intent(out) :: lwindow_loc
  real(kind=dp), dimension(3, num_wann_loc), optional, intent(out) :: wann_centres_loc
  real(kind=dp), dimension(num_wann_loc), optional, intent(out) :: wann_spreads_loc
  real(kind=dp), dimension(3), optional, intent(out) :: spread_loc

  real(kind=dp) time0, time1, time2
  character(len=9) :: stat, pos, cdate, ctime
  integer :: ierr, loop_k, loop_w
  logical :: wout_found

  integer :: nkp, nn, n, m

! Needed to split an array on different nodes
  integer, dimension(0:num_nodes - 1) :: counts
  integer, dimension(0:num_nodes - 1) :: displs

  complex(kind=dp), allocatable :: ham_r(:, :, :)
  integer, allocatable :: irvec(:, :)
  integer, allocatable :: shift_vec(:, :)
  integer, allocatable :: ndegen(:)
  integer :: rpt_origin
  real(kind=dp), allocatable :: wannier_centres_translated(:, :)
  complex(kind=dp), allocatable :: ham_k(:, :, :)
  integer :: nrpts

  type(sitesym_data) :: sym
  type(ham_logical) :: hmlg

  time0 = io_time()

  library = .true.
!  seedname="wannier"
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

  call io_date(cdate, ctime)

  write (stdout, '(/,2a,/)') ' Resuming Wannier90 at ', ctime

!  call param_write_header

  ! copy local data into module variables
  num_bands = num_bands_loc
  mp_grid = mp_grid_loc
  num_kpts = num_kpts_loc
  real_lattice = real_lattice_loc
  recip_lattice = recip_lattice_loc
  allocate (kpt_latt(3, num_kpts), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating kpt_latt in wannier_setup')
  kpt_latt = kpt_latt_loc
  allocate (eigval(num_bands, num_kpts), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating eigval in wannier_setup')
  eigval = eigenvalues_loc
  num_atoms = num_atoms_loc
  gamma_only = gamma_only_loc

  call param_lib_set_atoms(atom_symbols_loc, atoms_cart_loc)

  call param_read()

  call param_write()

  time1 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

  call kmesh_get(recip_lattice, kpt_cart, timing_level, nncell, neigh, &
                 nnlist, nntot, shell_list, devel_flag, iprint, lenconfac, &
                 kmesh_tol, num_kpts, search_shells, gamma_only, nnh, wbtot, &
                 skip_B1_tests, bk, bka, wb, num_shells, length_unit)

  time2 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to get kmesh        ', time2 - time1, ' (sec)'

  call comms_array_split(num_kpts, counts, displs)
  call overlap_allocate(u_matrix, m_matrix_local, m_matrix, u_matrix_opt, a_matrix, m_matrix_orig_local, &
                        m_matrix_orig, timing_level, nntot, num_kpts, num_wann, num_bands, disentanglement)

  if (disentanglement) then
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
  if (disentanglement) then
    call comms_scatterv(m_matrix_orig_local, num_bands*num_bands*nntot*counts(my_node_id), &
                        m_matrix_orig, num_bands*num_bands*nntot*counts, num_bands*num_bands*nntot*displs)
  else
    call comms_scatterv(m_matrix_local, num_wann*num_wann*nntot*counts(my_node_id), &
                        m_matrix, num_wann*num_wann*nntot*counts, num_wann*num_wann*nntot*displs)
  endif

!~  ! Check Mmn(k,b) is symmetric in m and n for gamma_only case
!~  if (gamma_only) call overlap_check_m_symmetry()

  if (disentanglement) then
    have_disentangled = .false.
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
    call param_write_chkpt('postdis')
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle      ', time1 - time2, ' (sec)'
  else
    if (gamma_only) then
      call overlap_project_gamma(nntot, m_matrix, u_matrix, timing_level, num_wann)  !lp note this not called by wannier_prog.F90
    else
      call overlap_project(m_matrix_local, nnlist, nntot, m_matrix, u_matrix, &
                           timing_level, num_kpts, num_wann, num_bands, lsitesymmetry, sym) !lp note this not called by wannier_prog.F90
    endif
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to project overlaps ', time1 - time2, ' (sec)'
  end if

  if (gamma_only) then
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
  else
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
                   translation_centre_frac, automatic_translation, ndimwin, &
                   sym, ham_r, irvec, shift_vec, ndegen, nrpts, rpt_origin, &
                   wannier_centres_translated, hmlg, ham_k)

  endif

  call param_write_chkpt('postwann')

  time2 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time2 - time1, ' (sec)'

  if (wannier_plot .or. bands_plot .or. fermi_surface_plot .or. write_hr) then
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
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  end if

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
  end if

  ! Now we zero all of the local output data, then copy in the data
  ! from the parameters module

  u_matrix_loc = u_matrix
  if (present(u_matrix_opt_loc) .and. present(lwindow_loc)) then
  if (disentanglement) then
    u_matrix_opt_loc = u_matrix_opt
    lwindow_loc = lwindow
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

  if (present(wann_centres_loc)) wann_centres_loc = wannier_centres
  if (present(wann_spreads_loc)) wann_spreads_loc = wannier_spreads
  if (present(spread_loc)) then
    spread_loc(1) = omega_total
    spread_loc(2) = omega_invariant
    spread_loc(3) = omega_tilde
  endif
  call hamiltonian_dealloc(ham_r, irvec, ndegen, wannier_centres_translated, &
                           hmlg, ham_k)
  call overlap_dealloc(m_matrix_orig_local, m_matrix_local, u_matrix_opt, &
                       a_matrix, m_matrix_orig, m_matrix, u_matrix)
  call kmesh_dealloc(nncell, neigh, nnlist, bk, bka, wb)
  call param_dealloc()

  write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time() - time0, ' (sec)'

  if (timing_level > 0) call io_print_timings()

  write (stdout, *)
  write (stdout, '(1x,a)') 'All done: wannier90 exiting'
  close (stdout)

end subroutine wannier_run
