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
  use w90_param_types
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

  use w90_lib

  use w90_param_methods, only: param_read, param_write, param_dealloc, &
    param_write_header, param_write_chkpt, param_read_chkpt, &
    param_memory_estimate, param_dist, param_chkpt_dist

  use w90_parameters

  implicit none

  ! data from parameters module
!  type(w90_calculation_type) :: w90_calcs
!  ! Are we running postw90?
!  !logical :: ispostw90 = .false.
!  type(pw90_calculation_type) :: pw90_calcs
!  type(param_driver_type) :: driver
!  type(postproc_type) :: pp_calc
!  type(parameter_input_type) :: param_input
!  ! only in parameters and postw90_common !BGS localise somehow
!  logical :: eig_found
!  type(param_plot_type) :: param_plot
!  type(postw90_oper_type) :: postw90_oper
!  type(param_wannierise_type) :: param_wannierise
!  ! RS: symmetry-adapted Wannier functions
!  logical :: lsitesymmetry = .false.
!  real(kind=dp) :: symmetrize_eps = 1.d-3
!  type(wannier_data_type) :: wann_data
!  type(param_hamiltonian_type) :: param_hamil
!  type(param_kmesh_type) :: kmesh_data
!  type(kmesh_info_type) :: kmesh_info
!  type(k_point_type) :: k_points
!  integer :: num_kpts !BGS put in k_point_type?
!  type(postw90_common_type) :: pw90_common
!  type(postw90_spin_type) :: pw90_spin
!  type(postw90_ham_type) :: pw90_ham
!  type(disentangle_type) :: dis_data
!  type(fermi_surface_type) :: fermi_surface_data
!  type(kpath_type) :: kpath
!  type(kslice_type) :: kslice
!  type(dos_plot_type) :: dos_data
!  type(berry_type) :: berry
!  type(spin_hall_type) :: spin_hall
!  type(gyrotropic_type) :: gyrotropic
!  type(fermi_data_type) :: fermi
!  type(geninterp_type) :: geninterp
!  type(boltzwann_type) :: boltz
!  type(transport_type) :: tran
!  type(atom_data_type) :: atoms

!  integer :: num_bands
  !! Number of bands

!  integer :: num_wann
  !! number of wannier functions

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states
  !BGS a_matrix, m_matrix in disentangle and overlap
!  complex(kind=dp), allocatable :: a_matrix(:, :, :)
!  complex(kind=dp), allocatable :: m_matrix_orig(:, :, :, :)
!  complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)
!  !BGS disentangle, hamiltonian, a wannierise print, and postw90/get_oper
!  real(kind=dp), allocatable :: eigval(:, :)

  !BGS need to sort these further, u_matrix in lots of places
  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states
!  complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)

  ! u_matrix gives the unitary rotations from the optimal subspace to the
  ! optimally smooth states.
  ! m_matrix we store here, becuase it is needed for restart of wannierise
!  complex(kind=dp), allocatable :: u_matrix(:, :, :)
  ! disentangle, hamiltonain, overlap and wannierise
!  complex(kind=dp), allocatable :: m_matrix(:, :, :, :)
  !BGS is disentangle and overlap
!  complex(kind=dp), allocatable :: m_matrix_local(:, :, :, :)

!  integer :: mp_grid(3)
  !! Dimensions of the Monkhorst-Pack grid

!  integer :: num_proj
  !BGS used by stuff in driver/kmesh/wannier - keep separate or duplicate?

!  type(select_projection_type) :: select_proj

!  real(kind=dp) :: real_lattice(3, 3)

  !parameters derived from input
!  real(kind=dp) :: recip_lattice(3, 3)

!  type(special_kpoints_type) :: spec_points
  ! end data from parameters module

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

  type(ham_logical) :: hmlg

  ! fixme
  !JJ data that will *not be referenced*; for calling library only
  !1/ setup
  integer :: num_atoms_loc
  integer :: num_bands_tot
  integer :: num_kpts_loc
  integer :: nntot_loc
  integer :: num_bands_loc
  integer :: num_wann_loc
  logical :: gamma_only_loc
  logical :: spinors_loc

  integer :: mp_grid_loc(3)
  integer, allocatable :: nncell_loc(:, :, :)
  integer, allocatable :: exclude_bands_loc(:)
  integer, allocatable :: proj_l_loc(:)
  integer, allocatable :: proj_m_loc(:)
  integer, allocatable :: proj_radial_loc(:)
  integer, allocatable :: proj_s_loc(:)
  integer, allocatable :: nnlist_loc(:, :)

  real(kind=dp) :: real_lattice_loc(3, 3)
  real(kind=dp) :: recip_lattice_loc(3, 3)
  real(kind=dp), allocatable :: atoms_cart_loc(:, :)
  real(kind=dp), allocatable :: proj_site_loc(:, :)
  real(kind=dp), allocatable :: proj_x_loc(:, :)
  real(kind=dp), allocatable :: proj_z_loc(:, :)
  real(kind=dp), allocatable :: proj_s_qaxis_loc(:, :)
  real(kind=dp), allocatable :: kpt_latt_loc(:, :)
  real(kind=dp), allocatable :: proj_zona_loc(:)

  character, allocatable :: atom_symbols_loc(:, :)

  !2/ run
  complex(kind=dp), allocatable :: m_matrix_loc(:, :, :, :)
  complex(kind=dp), allocatable :: a_matrix_loc(:, :, :)
  complex(kind=dp), allocatable :: u_matrix_loc(:, :, :)
  complex(kind=dp), allocatable :: u_matrix_opt_loc(:, :, :)
  logical, allocatable :: lwindow_loc(:, :)
  real(kind=dp), allocatable :: eigenvalues_loc(:, :)
  real(kind=dp), allocatable :: wann_centres_loc(:, :)
  real(kind=dp), allocatable :: wann_spreads_loc(:)
  real(kind=dp) :: spread_loc(3)
  ! end of vars for lib calls

  logical :: overlap = .true.

  call comms_setup

  driver%library = .false.

  time0 = io_time()

  if (on_root) then
    prog = 'wannier90'
    call io_commandline(prog, dryrun)
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1)
  call comms_bcast(seedname, len_seedname)
  call comms_bcast(dryrun, 1)

  ! JJ >>snip<<

  call wannier_setup(seedname, mp_grid, num_kpts_loc, &
                     real_lattice_loc, recip_lattice_loc, kpt_latt_loc, num_bands_tot, &
                     num_atoms_loc, atom_symbols_loc, atoms_cart_loc, gamma_only_loc, spinors_loc, &
                     nntot_loc, nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc, &
                     proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc, proj_z_loc, &
                     proj_x_loc, proj_zona_loc, exclude_bands_loc, proj_s_loc, proj_s_qaxis_loc)

  ! JJ >>snip<<

  if (on_root) then
    call param_memory_estimate(w90_calcs, param_input, param_wannierise, &
                               kmesh_data, kmesh_info, num_kpts, dis_data, &
                               atoms, num_bands, num_wann, num_proj, &
                               pw90_calcs, pw90_common, boltz, .false.)
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

  ! JJ >>snip<<

  if (w90_calcs%transport .and. tran%read_ht) then
    ! skip to tran_main; replaces goto 3003
    ! seems unnecessary?
    w90_calcs%disentanglement = .false.
    w90_calcs%wannierise = .false.
    overlap = .false.
  endif

  ! Sort out restarts
  if (driver%restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) then
      call param_read_chkpt(driver, param_input, wann_data, kmesh_info, &
                            k_points, num_kpts, dis_data, num_bands, num_wann, &
                            u_matrix, u_matrix_opt, m_matrix, mp_grid, &
                            real_lattice, recip_lattice, .false.)
    endif
    call param_chkpt_dist(driver, param_input, wann_data, num_kpts, dis_data, &
                          num_bands, num_wann, u_matrix, u_matrix_opt)

    if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
    if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

    select case (driver%restart)
    case ('default')    ! continue from where last checkpoint was written
      if (on_root) write (stdout, '(/1x,a)', advance='no') 'Resuming a previous Wannier90 calculation '

      if (driver%checkpoint .eq. 'postdis') then
        if (on_root) write (stdout, '(a/)') 'from wannierisation ...'

        w90_calcs%disentanglement = .false.
        overlap = .false.

      elseif (driver%checkpoint .eq. 'postwann') then
        if (on_root) write (stdout, '(a/)') 'from plotting ...'
        w90_calcs%disentanglement = .false.
        w90_calcs%wannierise = .false.
        overlap = .false.

      else
        if (on_root) write (stdout, '(/a/)')
        call io_error('Value of checkpoint not recognised in wann_prog')

      endif

    case ('wannierise') ! continue from wann_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
      w90_calcs%disentanglement = .false.
      w90_calcs%wannierise = .true.
      overlap = .false.

    case ('plot')       ! continue from plot_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
      ! fixme, the logic for choosing plot and trans in wannier_run needs fixing
      w90_calcs%disentanglement = .false.
      w90_calcs%wannierise = .false.
      overlap = .false.

    case ('transport')   ! continue from tran_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from transport routines ...'
      ! fixme, the logic for choosing plot and trans in wannier_run needs fixing
      w90_calcs%disentanglement = .false.
      w90_calcs%wannierise = .false.
      overlap = .false.

    case default        ! for completeness... (it is already trapped in param_read)
      call io_error('Value of restart not recognised in wann_prog')
    end select
  endif

  if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
  if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

  if (overlap) then
    time2 = io_time()
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
  endif

  ! JJ >>snip<<

  call wannier_run(seedname, mp_grid_loc, num_kpts_loc, real_lattice_loc, recip_lattice_loc, &
                   kpt_latt_loc, num_bands_loc, num_wann_loc, nntot_loc, num_atoms_loc, &
                   atom_symbols_loc, atoms_cart_loc, gamma_only_loc, m_matrix_loc, a_matrix_loc, &
                   eigenvalues_loc, u_matrix_loc, u_matrix_opt_loc=u_matrix_opt_loc, &
                   lwindow_loc=lwindow_loc, wann_centres_loc=wann_centres_loc, &
                   wann_spreads_loc=wann_spreads_loc, spread_loc=spread_loc)

  if (lsitesymmetry) call sitesym_dealloc(sym)

  if (param_input%timing_level > 0) call io_print_timings()

  write (stdout, *)
  write (stdout, '(1x,a)') 'All done: wannier90 exiting'
  close (stdout)

  call comms_end()
end program wannier
