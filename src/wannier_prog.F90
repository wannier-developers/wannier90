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

  use w90_constants, only: w90_physical_constants, dp
  use w90_param_types
  use w90_io
  use w90_hamiltonian
  use w90_kmesh
  use w90_disentangle
  use w90_overlap
  use w90_wannierise
  use w90_plot
  use w90_transport
  use w90_comms
  use w90_sitesym !YN:

  use w90_param_methods, only: param_write_header, param_read_chkpt, param_chkpt_dist
  use wannier_param_types
  use wannier_methods, only: param_read, param_w90_dealloc, param_write, &
    param_dist, param_memory_estimate, param_write_chkpt, w90_extra_io_type

#ifdef MPI
#  if !(defined(MPI08) || defined(MPI90) || defined(MPIH))
#    error "You need to define which MPI interface you are using"
#  endif
#endif

#ifdef MPI08
  use mpi_f08 ! use f08 interface if possible
#endif
#ifdef MPI90
  use mpi ! next best, use fortran90 interface
#endif

  implicit none

#ifdef MPIH
  include 'mpif.h' ! worst case, use legacy interface
#endif

  type(w90_physical_constants) :: physics
  ! data from parameters module
  type(w90_calculation_type) :: w90_calcs
  ! Are we running postw90?
  !logical :: ispostw90 = .false.
  character(len=20) :: checkpoint
  type(print_output_type) :: verbose
  type(w90_system_type) :: system
  type(real_space_ham_type) :: rs_region
  type(exclude_bands_type) :: excluded_bands
  logical :: eig_found
  type(param_plot_type) :: param_plot
  type(band_plot_type) :: band_plot
  type(wannier_plot_type) :: wann_plot
  type(param_wannierise_type) :: param_wannierise
  logical :: lsitesymmetry = .false. ! RS: symmetry-adapted Wannier functions
  real(kind=dp) :: symmetrize_eps = 1.d-3
  type(wannier_data_type) :: wann_data
  type(param_hamiltonian_type) :: param_hamil
  type(kmesh_input_type) :: kmesh_data
  type(proj_input_type) :: input_proj
  type(kmesh_info_type) :: kmesh_info
  type(k_point_type) :: k_points
  integer :: num_kpts
  type(disentangle_type) :: dis_data
  type(disentangle_manifold_type) :: dis_window
  type(fermi_surface_type) :: fermi_surface_data
  type(fermi_data_type) :: fermi
  type(transport_type) :: tran
  type(atom_data_type) :: atoms

  integer :: num_bands
  !! Number of bands

  integer :: num_wann
  !! number of wannier functions

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states
  !BGS a_matrix, m_matrix in disentangle and overlap
  complex(kind=dp), allocatable :: a_matrix(:, :, :)
  complex(kind=dp), allocatable :: m_matrix_orig(:, :, :, :)
  complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)
  real(kind=dp), allocatable :: eigval(:, :)

  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states
  complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)

  ! u_matrix gives the unitary rotations from the optimal subspace to the
  ! optimally smooth states.
  ! m_matrix we store here, becuase it is needed for restart of wannierise
  complex(kind=dp), allocatable :: u_matrix(:, :, :)
  ! disentangle, hamiltonain, overlap and wannierise
  complex(kind=dp), allocatable :: m_matrix(:, :, :, :)
  complex(kind=dp), allocatable :: m_matrix_local(:, :, :, :)

  integer :: mp_grid(3)
  !! Dimensions of the Monkhorst-Pack grid

  integer :: num_proj

  type(select_projection_type) :: select_proj

  real(kind=dp) :: real_lattice(3, 3)

  !parameters derived from input
  real(kind=dp) :: recip_lattice(3, 3)

  type(special_kpoints_type) :: spec_points
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
  integer :: stdout
  character(len=50) :: prog
  character(len=50) :: seedname
  logical :: on_root = .false.
  integer :: num_nodes, my_node_id, ierr

  type(w90_extra_io_type) :: write_data
  ! was in driver, only used by wannier_lib
  type(proj_input_type) :: proj
  !Projections
  logical :: lhasproj

  logical :: use_bloch_phases, cp_pp, calc_only_A
  logical :: gamma_only
  logical :: have_disentangled

  type(sitesym_data) :: sym
  type(ham_logical) :: hmlg
  type(w90commtype) :: w90comm

#ifdef MPI
  w90comm%comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  if (ierr .ne. 0) call io_error('MPI initialisation error', stdout, seedname)  ! JJ, fixme, what are stdout, seedname here?  unassigned!
#endif

  num_nodes = mpisize(w90comm)
  my_node_id = mpirank(w90comm)
  if (my_node_id == 0) on_root = .true.

  time0 = io_time()

  if (on_root) then
    prog = 'wannier90'
    call io_commandline(prog, dryrun, seedname)
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1, stdout, seedname, w90comm)
  call comms_bcast(seedname, len_seedname, stdout, seedname, w90comm)
  call comms_bcast(dryrun, 1, stdout, seedname, w90comm)

  if (on_root) then
    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.werr')
    call io_date(cdate, ctime)
    write (stdout, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

    call param_read(atoms, band_plot, dis_data, dis_window, excluded_bands, fermi, &
                    fermi_surface_data, kmesh_data, kmesh_info, k_points, param_hamil, &
                    param_plot, param_wannierise, proj, input_proj, rs_region, select_proj, &
                    spec_points, system, tran, verbose, wann_data, wann_plot, write_data, &
                    w90_calcs, eigval, real_lattice, recip_lattice, physics%bohr, symmetrize_eps, &
                    mp_grid, num_bands, num_kpts, num_proj, num_wann, eig_found, calc_only_A, &
                    cp_pp, gamma_only, lhasproj, .false., .false., lsitesymmetry, &
                    use_bloch_phases, seedname, stdout)
    have_disentangled = .false.
    close (stdout, status='delete')

    if (w90_calcs%restart .eq. ' ') then
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
    call param_write_header(physics%bohr_version_str, physics%constants_version_str1, &
                            physics%constants_version_str2, stdout)
    if (num_nodes == 1) then
#ifdef MPI
      write (stdout, '(/,1x,a)') 'Running in serial (with parallel executable)'
#else
      write (stdout, '(/,1x,a)') 'Running in serial (with serial executable)'
#endif
    else
      write (stdout, '(/,1x,a,i3,a/)') 'Running in parallel on ', num_nodes, ' CPUs'
    endif
    call param_write(atoms, band_plot, dis_data, fermi, fermi_surface_data, k_points, &
                     param_hamil, param_plot, param_wannierise, proj, input_proj, &
                     rs_region, select_proj, spec_points, tran, verbose, wann_data, wann_plot, &
                     write_data, w90_calcs, real_lattice, recip_lattice, symmetrize_eps, mp_grid, &
                     num_bands, num_kpts, num_proj, num_wann, cp_pp, gamma_only, lsitesymmetry, &
                     system%spinors, use_bloch_phases, stdout)
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

    if (.not. kmesh_info%explicit_nnkpts) call kmesh_get(kmesh_data, kmesh_info, verbose, &
                                                         k_points%kpt_cart, recip_lattice, num_kpts, &
                                                         gamma_only, seedname, stdout)
    time2 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to get kmesh        ', time2 - time1, ' (sec)'

    call param_memory_estimate(atoms, kmesh_info, param_wannierise, input_proj, verbose, &
                               w90_calcs, num_bands, num_kpts, num_proj, num_wann, gamma_only, &
                               stdout)
  end if !on_root

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
  call param_dist(atoms, band_plot, dis_data, dis_window, excluded_bands, fermi, &
                  fermi_surface_data, kmesh_data, kmesh_info, k_points, param_hamil, param_plot, &
                  param_wannierise, input_proj, rs_region, system, tran, verbose, wann_data, &
                  wann_plot, w90_calcs, eigval, real_lattice, recip_lattice, symmetrize_eps, &
                  mp_grid, num_bands, num_kpts, num_proj, num_wann, eig_found, cp_pp, gamma_only, &
                  have_disentangled, lhasproj, lsitesymmetry, use_bloch_phases, seedname, &
                  stdout, w90comm)
  if (gamma_only .and. num_nodes > 1) &
    call io_error('Gamma point branch is serial only at the moment', stdout, seedname)

  if (w90_calcs%transport .and. tran%read_ht) goto 3003

  ! Sort out restarts
  if (w90_calcs%restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) then
      call param_read_chkpt(dis_window, excluded_bands, kmesh_info, k_points, wann_data, m_matrix, &
                            u_matrix, u_matrix_opt, real_lattice, recip_lattice, &
                            param_wannierise%omega%invariant, mp_grid, num_bands, num_kpts, &
                            num_wann, checkpoint, have_disentangled, .false., &
                            seedname, stdout)
    endif
    call param_chkpt_dist(dis_window, wann_data, u_matrix, u_matrix_opt, &
                          param_wannierise%omega%invariant, num_bands, num_kpts, num_wann, &
                          checkpoint, have_disentangled, seedname, stdout, w90comm)
    if (lsitesymmetry) call sitesym_read(sym, num_bands, num_kpts, num_wann, seedname, stdout)  ! update this to read on root and bcast - JRY
    if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

    select case (w90_calcs%restart)
    case ('default')    ! continue from where last checkpoint was written
      if (on_root) write (stdout, '(/1x,a)', advance='no') &
        'Resuming a previous Wannier90 calculation '
      if (checkpoint .eq. 'postdis') then
        if (on_root) write (stdout, '(a/)') 'from wannierisation ...'
        goto 1001         ! go to wann_main
      elseif (checkpoint .eq. 'postwann') then
        if (on_root) write (stdout, '(a/)') 'from plotting ...'
        goto 2002         ! go to plot_main
      else
        if (on_root) write (stdout, '(/a/)')
        call io_error('Value of checkpoint not recognised in wann_prog', stdout, seedname)
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
      call io_error('Value of restart not recognised in wann_prog', stdout, seedname)
    end select
  endif

  if (w90_calcs%postproc_setup) then
    if (on_root) call kmesh_write(excluded_bands, kmesh_info, input_proj, verbose, &
                                  k_points%kpt_latt, real_lattice, recip_lattice, num_kpts, &
                                  num_proj, calc_only_A, system%spinors, seedname, stdout)
    call kmesh_dealloc(kmesh_info, stdout, seedname)
    call param_w90_dealloc(atoms, band_plot, dis_data, dis_window, excluded_bands, kmesh_data, &
                           k_points, param_wannierise, proj, input_proj, spec_points, &
                           wann_data, wann_plot, write_data, eigval, seedname, stdout)
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (on_root) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
    call comms_end
    stop
  endif

  if (lsitesymmetry) call sitesym_read(sym, num_bands, num_kpts, num_wann, seedname, stdout) ! update this to read on root and bcast - JRY
  if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

  call overlap_allocate(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, &
                        u_matrix, u_matrix_opt, kmesh_info%nntot, num_bands, num_kpts, num_wann, &
                        verbose%timing_level, w90_calcs%disentanglement, seedname, stdout, &
                        w90comm)
  call overlap_read(kmesh_info, select_proj, sym, w90_calcs, a_matrix, m_matrix, m_matrix_local, &
                    m_matrix_orig, m_matrix_orig_local, u_matrix, u_matrix_opt, num_bands, &
                    num_kpts, num_proj, num_wann, verbose%timing_level, cp_pp, &
                    gamma_only, lsitesymmetry, use_bloch_phases, seedname, stdout, w90comm)
  time1 = io_time()
  if (on_root) write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, &
    ' (sec)'

  have_disentangled = .false.

  if (w90_calcs%disentanglement) then

    call dis_main(dis_data, dis_window, kmesh_info, k_points, sym, verbose, a_matrix, m_matrix, &
                  m_matrix_local, m_matrix_orig, m_matrix_orig_local, u_matrix, u_matrix_opt, &
                  eigval, recip_lattice, param_wannierise%omega%invariant, num_bands, num_kpts, &
                  num_wann, gamma_only, lsitesymmetry, stdout, seedname, w90comm)
    have_disentangled = .true.
    time2 = io_time()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, &
      ' (sec)'
  endif

  if (on_root) then
    call param_write_chkpt('postdis', excluded_bands, wann_data, kmesh_info, &
                           k_points, num_kpts, dis_window, num_bands, num_wann, u_matrix, &
                           u_matrix_opt, m_matrix, mp_grid, real_lattice, recip_lattice, &
                           param_wannierise%omega%invariant, have_disentangled, stdout, seedname)
  endif
!~  call param_write_um

1001 time2 = io_time()

  ! JJ hack to workaround mpi_scatterv requirement that all arrays are valid *for all mpi procs*
  ! m_matrix* usually alloc'd in overlaps.F90, but not always
  if (.not. allocated(m_matrix)) allocate (m_matrix(1, 1, 1, 1)) !JJ temporary workaround to avoid runtime check failure

  if (.not. gamma_only) then
    call wann_main(atoms, dis_window, excluded_bands, hmlg, kmesh_info, k_points, param_hamil, &
                   param_wannierise, rs_region, sym, system, verbose, wann_data, w90_calcs, ham_k, &
                   ham_r, m_matrix, u_matrix, u_matrix_opt, eigval, real_lattice, recip_lattice, &
                   wannier_centres_translated, irvec, mp_grid, ndegen, shift_vec, nrpts, &
                   num_bands, num_kpts, num_proj, num_wann, rpt_origin, band_plot%plot_mode, &
                   tran%mode, have_disentangled, lsitesymmetry, seedname, stdout, w90comm)
  else
    call wann_main_gamma(atoms, dis_window, excluded_bands, kmesh_info, k_points, &
                         param_wannierise, system, verbose, wann_data, w90_calcs, m_matrix, &
                         u_matrix, u_matrix_opt, eigval, real_lattice, recip_lattice, mp_grid, &
                         num_bands, num_kpts, num_wann, have_disentangled, seedname, stdout, &
                         w90comm)
  end if

  time1 = io_time()
  if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time1 - time2, &
    ' (sec)'

  if (on_root) then
    call param_write_chkpt('postwann', excluded_bands, wann_data, kmesh_info, &
                           k_points, num_kpts, dis_window, num_bands, num_wann, u_matrix, &
                           u_matrix_opt, m_matrix, mp_grid, real_lattice, recip_lattice, &
                           param_wannierise%omega%invariant, have_disentangled, stdout, seedname)
  endif

2002 continue
  if (on_root) then
    ! I call the routine always; the if statements to decide if/what to plot are inside the function
    time2 = io_time()

    call plot_main(atoms, band_plot, dis_window, fermi, fermi_surface_data, hmlg, kmesh_info, &
                   k_points, param_hamil, param_plot, rs_region, spec_points, verbose, wann_data, &
                   wann_plot, w90_calcs, ham_k, ham_r, m_matrix, u_matrix, u_matrix_opt, eigval, &
                   real_lattice, recip_lattice, wannier_centres_translated, physics%bohr, irvec, &
                   mp_grid, ndegen, shift_vec, nrpts, num_bands, num_kpts, num_wann, rpt_origin, &
                   tran%mode, have_disentangled, lsitesymmetry, system%spinors, seedname, stdout)
    time1 = io_time()

    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
    if (w90_calcs%transport) then
      time2 = io_time()

      call tran_main(atoms, dis_window, fermi, hmlg, k_points, param_hamil, &
                     rs_region, tran, verbose, wann_data, w90_calcs, ham_k, ham_r, u_matrix, &
                     u_matrix_opt, eigval, real_lattice, recip_lattice, &
                     wannier_centres_translated, irvec, mp_grid, ndegen, shift_vec, nrpts, &
                     num_bands, num_kpts, num_wann, rpt_origin, band_plot%plot_mode, &
                     have_disentangled, lsitesymmetry, seedname, stdout)
      time1 = io_time()

      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (tran%read_ht) goto 4004
    end if
  endif

  call tran_dealloc(stdout, seedname)
  call hamiltonian_dealloc(hmlg, ham_k, ham_r, wannier_centres_translated, irvec, ndegen, &
                           stdout, seedname)
  call overlap_dealloc(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, &
                       u_matrix, u_matrix_opt, seedname, stdout, w90comm)
  call kmesh_dealloc(kmesh_info, stdout, seedname)
  call param_w90_dealloc(atoms, band_plot, dis_data, dis_window, excluded_bands, kmesh_data, &
                         k_points, param_wannierise, proj, input_proj, spec_points, &
                         wann_data, wann_plot, write_data, eigval, seedname, stdout)
  if (lsitesymmetry) call sitesym_dealloc(sym, stdout, seedname) !YN:

4004 continue

  if (on_root) then
    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time(), ' (sec)'

    if (verbose%timing_level > 0) call io_print_timings(stdout)

    write (stdout, *)
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'

    close (stdout)
  endif

  call comms_end

end program wannier

