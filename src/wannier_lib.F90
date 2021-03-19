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
  use wannier_param_data
  use wannier_methods, only: param_read, param_write, proj
  use w90_sitesym
  use w90_io
  use w90_kmesh
  use w90_comms, only: comms_setup_vars

  use w90_param_methods, only: param_write_header, param_lib_set_atoms

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

  real(kind=dp) time0, time1
  character(len=9) :: stat, pos, cdate, ctime
  integer :: ierr
  logical :: wout_found

  time0 = io_time()

  call comms_setup_vars

  !library = .true.
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
  allocate (k_points%kpt_latt(3, num_kpts), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating kpt_latt in wannier_setup')
  k_points%kpt_latt = kpt_latt_loc
  atoms%num_atoms = num_atoms_loc
  call param_lib_set_atoms(atoms, atom_symbols_loc, atoms_cart_loc, recip_lattice)
  param_input%gamma_only = gamma_only_loc
  param_input%spinors = spinors_loc

  ! GP: at this point we don't know yet the number of excluded bands...
  num_bands = num_bands_tot
  !library_param_read_first_pass = .true.
  call param_read(driver, w90_calcs, pp_calc, param_input, param_plot, &
                  param_wannierise, lsitesymmetry, symmetrize_eps, &
                  wann_data, param_hamil, kmesh_data, kmesh_info, &
                  k_points, num_kpts, dis_data, fermi_surface_data, &
                  fermi, tran, atoms, num_bands, num_wann, eigval, &
                  mp_grid, num_proj, select_proj, real_lattice, &
                  recip_lattice, spec_points, eig_found, .true., .true.)
  ! Following calls will all NOT be first_pass, and I need to pass
  ! directly num_bands, that is already set internally now to num_bands = num_bands_tot - num_exclude_bands
  !library_param_read_first_pass = .false.

  call param_write(driver, w90_calcs, param_input, param_plot, &
                   param_wannierise, lsitesymmetry, symmetrize_eps, &
                   wann_data, param_hamil, kmesh_data, k_points, &
                   num_kpts, dis_data, fermi_surface_data, fermi, tran, &
                   atoms, num_bands, num_wann, mp_grid, num_proj, &
                   select_proj, real_lattice, recip_lattice, spec_points)

  time1 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

  if (.not. driver%explicit_nnkpts) call kmesh_get(recip_lattice, k_points%kpt_cart, param_input, &
                                                   kmesh_info, kmesh_data, num_kpts)

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

  nntot_loc = kmesh_info%nntot
  nnlist_loc(:, 1:kmesh_info%nntot) = kmesh_info%nnlist(:, 1:kmesh_info%nntot)
  nncell_loc(:, :, 1:kmesh_info%nntot) = kmesh_info%nncell(:, :, 1:kmesh_info%nntot)
  num_bands_loc = num_bands
  num_wann_loc = num_wann
  if (allocated(param_wannierise%proj_site)) then
    proj_site_loc(:, 1:num_proj) = param_wannierise%proj_site(:, 1:num_proj)
    proj_l_loc(1:num_proj) = proj%l(1:num_proj)
    proj_m_loc(1:num_proj) = proj%m(1:num_proj)
    proj_z_loc(:, 1:num_proj) = proj%z(:, 1:num_proj)
    proj_x_loc(:, 1:num_proj) = proj%x(:, 1:num_proj)
    proj_radial_loc(1:num_proj) = proj%radial(1:num_proj)
    proj_zona_loc(1:num_proj) = proj%zona(1:num_proj)
    if (allocated(proj%s) .and. present(proj_s_loc) .and. present(proj_s_qaxis_loc)) then
      proj_s_loc(1:num_proj) = proj%s(1:num_proj)
      proj_s_qaxis_loc(:, 1:num_proj) = proj%s_qaxis(:, 1:num_proj)
    end if
  endif
  if (allocated(param_input%exclude_bands)) then
    exclude_bands_loc(1:param_input%num_exclude_bands) = &
      param_input%exclude_bands(1:param_input%num_exclude_bands)
  end if

  if (driver%postproc_setup) then
    call kmesh_write(recip_lattice, param_input, kmesh_info, num_kpts, kmesh_data, num_proj, &
                     k_points%kpt_latt, real_lattice, pp_calc%only_A)
    write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    write (stdout, '(/a)') ' '//trim(seedname)//'.nnkp written.'
  endif

  call kmesh_dealloc(kmesh_info)
  call param_w90_dealloc(param_input, param_plot, param_wannierise, &
                         wann_data, kmesh_data, k_points, dis_data, fermi, &
                         atoms, eigval, spec_points)
  write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'

  write (stdout, '(/a/)') ' Finished setting up k-point neighbours.'

  call io_date(cdate, ctime)

  write (stdout, '(2a)') ' Exiting wannier_setup in wannier90 ', ctime

  close (stdout)

end subroutine wannier_setup

! setup uses the same arguments as run and some extra ones:
!subroutine wannier_setup(
!                         spinors_loc, &
!                          nnlist_loc, nncell_loc
!                         proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc, proj_z_loc, &
!                         proj_x_loc, proj_zona_loc, exclude_bands_loc, proj_s_loc, proj_s_qaxis_loc)
!
!

subroutine wannier_run(seed__name, mp_grid_loc, num_kpts_loc, real_lattice_loc, recip_lattice_loc, &
                       kpt_latt_loc, num_bands_loc, num_wann_loc, nntot_loc, num_atoms_loc, &
                       atom_symbols_loc, atoms_cart_loc, gamma_only_loc, M_matrix_loc, &
                       A_matrix_loc, eigenvalues_loc, U_matrix_loc, U_matrix_opt_loc, lwindow_loc, &
                       wann_centres_loc, wann_spreads_loc, spread_loc)

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
  use wannier_param_data
  use wannier_methods, only: param_read, param_write, param_write_chkpt, &
    param_w90_dealloc
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

  use w90_param_methods, only: param_lib_set_atoms

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
  complex(kind=dp), dimension(num_bands_loc, num_bands_loc, nntot_loc, num_kpts_loc), &
    intent(in) :: M_matrix_loc
  complex(kind=dp), dimension(num_bands_loc, num_wann_loc, num_kpts_loc), intent(in) :: A_matrix_loc
  real(kind=dp), dimension(num_bands_loc, num_kpts_loc), intent(in) :: eigenvalues_loc
  complex(kind=dp), dimension(num_wann_loc, num_wann_loc, num_kpts_loc), intent(out) :: U_matrix_loc
  complex(kind=dp), dimension(num_bands_loc, num_wann_loc, num_kpts_loc), optional, &
    intent(out) :: U_matrix_opt_loc
  logical, dimension(num_bands_loc, num_kpts_loc), optional, intent(out) :: lwindow_loc
  real(kind=dp), dimension(3, num_wann_loc), optional, intent(out) :: wann_centres_loc
  real(kind=dp), dimension(num_wann_loc), optional, intent(out) :: wann_spreads_loc
  real(kind=dp), dimension(3), optional, intent(out) :: spread_loc

  real(kind=dp) time0, time1, time2
  character(len=9) :: stat, pos, cdate, ctime
  integer :: ierr, loop_k, loop_w
  logical :: wout_found

  !integer :: nkp, nn, n, m

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

  !driver%library = .true.
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
  allocate (k_points%kpt_latt(3, num_kpts), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating kpt_latt in wannier_setup')
  k_points%kpt_latt = kpt_latt_loc
  allocate (eigval(num_bands, num_kpts), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating eigval in wannier_setup')
  eigval = eigenvalues_loc
  atoms%num_atoms = num_atoms_loc
  param_input%gamma_only = gamma_only_loc

  call param_lib_set_atoms(atoms, atom_symbols_loc, atoms_cart_loc, recip_lattice)

  call param_read(driver, w90_calcs, pp_calc, param_input, param_plot, &
                  param_wannierise, lsitesymmetry, symmetrize_eps, &
                  wann_data, param_hamil, kmesh_data, kmesh_info, &
                  k_points, num_kpts, dis_data, fermi_surface_data, &
                  fermi, tran, atoms, num_bands, num_wann, eigval, &
                  mp_grid, num_proj, select_proj, real_lattice, &
                  recip_lattice, spec_points, eig_found, .true., .false.)

  call param_write(driver, w90_calcs, param_input, param_plot, &
                   param_wannierise, lsitesymmetry, symmetrize_eps, &
                   wann_data, param_hamil, kmesh_data, k_points, &
                   num_kpts, dis_data, fermi_surface_data, fermi, tran, &
                   atoms, num_bands, num_wann, mp_grid, num_proj, &
                   select_proj, real_lattice, recip_lattice, spec_points)

  time1 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

  call kmesh_get(recip_lattice, k_points%kpt_cart, param_input, kmesh_info, kmesh_data, num_kpts)

  time2 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to get kmesh        ', time2 - time1, ' (sec)'

  call comms_array_split(num_kpts, counts, displs)
  call overlap_allocate(u_matrix, m_matrix_local, m_matrix, u_matrix_opt, a_matrix, &
                        m_matrix_orig_local, m_matrix_orig, param_input%timing_level, &
                        kmesh_info%nntot, num_kpts, num_wann, num_bands, w90_calcs%disentanglement)

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
    call comms_scatterv(m_matrix_orig_local, &
                        num_bands*num_bands*kmesh_info%nntot*counts(my_node_id), m_matrix_orig, &
                        num_bands*num_bands*kmesh_info%nntot*counts, &
                        num_bands*num_bands*kmesh_info%nntot*displs)
  else
    call comms_scatterv(m_matrix_local, num_wann*num_wann*kmesh_info%nntot*counts(my_node_id), &
                        m_matrix, num_wann*num_wann*kmesh_info%nntot*counts, num_wann*num_wann* &
                        kmesh_info%nntot*displs)
  endif

!~  ! Check Mmn(k,b) is symmetric in m and n for gamma_only case
!~  if (gamma_only) call overlap_check_m_symmetry()

  if (w90_calcs%disentanglement) then
    param_input%have_disentangled = .false.

    call dis_main(num_bands, num_kpts, num_wann, recip_lattice, eigval, a_matrix, m_matrix, &
                  m_matrix_local, m_matrix_orig, m_matrix_orig_local, u_matrix, u_matrix_opt, &
                  dis_data, kmesh_info, k_points, param_input, num_nodes, my_node_id, on_root, &
                  lsitesymmetry, sym)

    param_input%have_disentangled = .true.
    call param_write_chkpt('postdis', param_input, wann_data, kmesh_info, k_points, num_kpts, &
                           dis_data, num_bands, num_wann, u_matrix, u_matrix_opt, m_matrix, &
                           mp_grid, real_lattice, recip_lattice)
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle      ', time1 - time2, ' (sec)'
  else
    if (param_input%gamma_only) then
      call overlap_project_gamma(kmesh_info%nntot, m_matrix, u_matrix, param_input%timing_level, &
                                 num_wann)  !lp note this not called by wannier_prog.F90
    else
      call overlap_project(m_matrix_local, kmesh_info%nnlist, kmesh_info%nntot, m_matrix, &
                           u_matrix, param_input%timing_level, num_kpts, num_wann, num_bands, &
                           lsitesymmetry, sym) !lp note this not called by wannier_prog.F90
    endif
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to project overlaps ', time1 - time2, ' (sec)'
  end if

  if (param_input%gamma_only) then
    call wann_main_gamma(num_wann, param_wannierise, kmesh_info, param_input, u_matrix, m_matrix, &
                         num_kpts, real_lattice, wann_data, num_bands, u_matrix_opt, eigval, &
                         dis_data%lwindow, recip_lattice, atoms, k_points, dis_data, mp_grid, &
                         stdout)
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

  call param_write_chkpt('postwann', param_input, wann_data, kmesh_info, k_points, num_kpts, &
                         dis_data, num_bands, num_wann, u_matrix, u_matrix_opt, m_matrix, mp_grid, &
                         real_lattice, recip_lattice)

  time2 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time2 - time1, ' (sec)'

  if (w90_calcs%wannier_plot .or. w90_calcs%bands_plot .or. w90_calcs%fermi_surface_plot .or. w90_calcs%write_hr) then
    call plot_main(num_kpts, w90_calcs, k_points, param_input, &
                   param_plot, real_lattice, num_wann, kmesh_info, m_matrix, &
                   recip_lattice, wann_data, atoms, param_hamil, dis_data, &
                   u_matrix_opt, eigval, u_matrix, lsitesymmetry, num_bands, &
                   mp_grid, tran%mode, fermi, fermi_surface_data, spec_points, &
                   ham_r, irvec, shift_vec, ndegen, nrpts, rpt_origin, &
                   wannier_centres_translated, hmlg, ham_k)

    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  end if

  time2 = io_time()
  if (w90_calcs%transport) then
    call tran_main(tran, param_input, w90_calcs, num_wann, real_lattice, recip_lattice, wann_data, &
                   atoms, param_hamil, dis_data, u_matrix_opt, k_points, eigval, u_matrix, &
                   lsitesymmetry, num_bands, num_kpts, mp_grid, fermi, ham_r, irvec, shift_vec, &
                   ndegen, nrpts, rpt_origin, wannier_centres_translated, hmlg, ham_k)

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
  call hamiltonian_dealloc(ham_r, irvec, ndegen, wannier_centres_translated, &
                           hmlg, ham_k)
  call overlap_dealloc(m_matrix_orig_local, m_matrix_local, u_matrix_opt, &
                       a_matrix, m_matrix_orig, m_matrix, u_matrix)
  call kmesh_dealloc(kmesh_info)
  call param_w90_dealloc(param_input, param_plot, param_wannierise, &
                         wann_data, kmesh_data, k_points, dis_data, &
                         atoms, eigval, spec_points)

  write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time() - time0, ' (sec)'

  if (param_input%timing_level > 0) call io_print_timings()

  write (stdout, *)
  write (stdout, '(1x,a)') 'All done: wannier90 exiting'
  close (stdout)

end subroutine wannier_run
