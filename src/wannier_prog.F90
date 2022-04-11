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

  use w90_constants, only: w90_physical_constants_type, dp
  use w90_error
  use w90_types
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

  use w90_readwrite, only: w90_readwrite_write_header, w90_readwrite_read_chkpt, &
    w90_readwrite_chkpt_dist, w90_readwrite_in_file, w90_readwrite_clean_infile, &
    w90_readwrite_read_final_alloc, w90_readwrite_uppercase
  use w90_wannier90_types
  use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, &
    w90_wannier90_readwrite_w90_dealloc, &
    w90_wannier90_readwrite_write, w90_wannier90_readwrite_dist, &
    w90_wannier90_readwrite_memory_estimate, &
    w90_wannier90_readwrite_write_chkpt, w90_extra_io_type

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

  type(atom_data_type) :: atom_data
  type(band_plot_type) :: band_plot
  type(dis_control_type) :: dis_control
  type(dis_manifold_type) :: dis_manifold
  type(dis_spheres_type) :: dis_spheres
  type(fermi_surface_plot_type) :: fermi_surface_plot
  type(ham_logical_type) :: ham_logical
  type(kmesh_info_type) :: kmesh_info
  type(kmesh_input_type) :: kmesh_input
  type(kpoint_path_type) :: kpoint_path
  type(output_file_type) :: output_file
  type(print_output_type) :: print_output
  type(proj_input_type) :: input_proj
  type(real_space_ham_type) :: real_space_ham
  type(select_projection_type) :: select_projection
  type(sitesym_type) :: sitesym
  type(transport_type) :: transport
  type(w90_calculation_type) :: w90_calculation
  type(w90comm_type) :: comm
  type(w90_physical_constants_type) :: physics
  type(w90_system_type) :: w90_system
  type(wann_control_type) :: wann_control
  type(wannier_data_type) :: wannier_data
  type(wannier_plot_type) :: wannier_plot
  type(wann_omega_type) :: omega
  type(ws_region_type) :: ws_region
  type(wvfn_read_type) :: wvfn_read

  integer, allocatable :: exclude_bands(:)
  integer :: mp_grid(3)  !! Dimensions of the Monkhorst-Pack grid
  integer :: num_bands   !! Number of bands
  integer :: num_exclude_bands
  integer :: num_kpts
  integer :: num_proj
  integer :: num_wann    !! number of wannier functions
  integer :: optimisation

  real(kind=dp), allocatable :: kpt_latt(:, :) !! kpoints in lattice vecs
  real(kind=dp), allocatable :: eigval(:, :)
  real(kind=dp), allocatable :: fermi_energy_list(:)
  real(kind=dp) :: real_lattice(3, 3)
  real(kind=dp) :: recip_lattice(3, 3)

  ! Are we running postw90?
  !logical :: ispostw90 = .false.

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states
  !BGS a_matrix, m_matrix in disentangle and overlap
  complex(kind=dp), allocatable :: a_matrix(:, :, :)
  complex(kind=dp), allocatable :: m_matrix_orig(:, :, :, :)
  complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)
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

  integer :: rpt_origin
  !! index of R=0
  integer :: nrpts
  !! number of Wigner-Seitz grid points
  integer, allocatable :: irvec(:, :)
  !! The irpt-th Wigner-Seitz grid point has components
  !! irvec(1:3,irpt) in the basis of the lattice vectors
  integer, allocatable :: shift_vec(:, :)
  integer, allocatable :: ndegen(:)
  !! Weight of the irpt-th point is 1/ndegen(irpt)
  real(kind=dp), allocatable :: wannier_centres_translated(:, :)
  complex(kind=dp), allocatable :: ham_r(:, :, :)
  !! Hamiltonian matrix in WF representation
  complex(kind=dp), allocatable :: ham_k(:, :, :)

  type(proj_input_type) :: proj_input ! Projections
  type(w90_extra_io_type) :: w90_extra_io

  real(kind=dp) time0, time1, time2

  integer :: len_seedname
  integer :: num_nodes, my_node_id, ierr
  integer :: stdout, stderr

  logical :: eig_found
  logical :: gamma_only
  logical :: have_disentangled, disentanglement
  logical :: lhasproj
  logical :: lsitesymmetry = .false. ! RS: symmetry-adapted Wannier functions
  logical :: on_root
  logical :: use_bloch_phases, cp_pp, calc_only_A
  logical :: wout_found, dryrun

  character(len=20) :: checkpoint
  character(len=50) :: prog
  character(len=50) :: seedname
  character(len=9) :: stat, pos, cdate, ctime

  type(timer_list_type) :: timer
  type(w90_error_type), allocatable :: error

  stdout = 6 ! stdout stream
  stderr = 0 ! stderr stream
  prog = "wannier90"
  seedname = "wannier"

#ifdef MPI
  comm%comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  if (ierr .ne. 0) then
    call set_error_mpi(error, 'MPI initialisation error', comm)
    call prterr(error, stdout, stderr, comm)
  endif
#endif

  num_nodes = mpisize(comm)
  my_node_id = mpirank(comm)
  on_root = .false.
  if (my_node_id == 0) on_root = .true.

  ! global inits (non-type based) from readwrite files
  optimisation = 3
  num_wann = -99
  num_proj = 0
  lhasproj = .false.
  !lsitesymmetry = .false.
  cp_pp = .false.
  calc_only_A = .false.
  use_bloch_phases = .false.
  eig_found = .false.
  ! end global inits

  time0 = io_time()

  if (on_root) then
    call io_commandline(prog, dryrun, seedname)
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1, error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)
  call comms_bcast(seedname, len_seedname, error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)
  call comms_bcast(dryrun, 1, error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  if (on_root) then
    open (newunit=stderr, file=trim(seedname)//'.werr') ! redefining from unit 0
    call io_date(cdate, ctime)
    write (stderr, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

    ! read infile to in_data structure
    call w90_readwrite_in_file(seedname, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)

    call w90_wannier90_readwrite_read(atom_data, band_plot, dis_control, dis_spheres, &
                                      dis_manifold, exclude_bands, fermi_energy_list, &
                                      fermi_surface_plot, kmesh_input, kmesh_info, kpt_latt, &
                                      output_file, wvfn_read, wann_control, omega, proj_input, &
                                      input_proj, real_space_ham, select_projection, kpoint_path, &
                                      w90_system, transport, print_output, wannier_data, &
                                      wannier_plot, w90_extra_io, ws_region, w90_calculation, &
                                      eigval, real_lattice, physics%bohr, sitesym%symmetrize_eps, &
                                      mp_grid, num_bands, num_kpts, num_proj, num_wann, &
                                      optimisation, eig_found, calc_only_A, cp_pp, gamma_only, &
                                      lhasproj, .false., .false., lsitesymmetry, use_bloch_phases, &
                                      seedname, stderr, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)

    call w90_readwrite_clean_infile(stdout, seedname, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)

    disentanglement = (num_bands > num_wann)

    if (.not. (w90_calculation%transport .and. transport%read_ht)) then
      ! For aesthetic purposes, convert some things to uppercase
      call w90_readwrite_uppercase(atom_data, kpoint_path, print_output%length_unit)
      ! Initialise
      omega%total = -999.0_dp
      omega%tilde = -999.0_dp
      omega%invariant = -999.0_dp
      call w90_readwrite_read_final_alloc(disentanglement, dis_manifold, wannier_data, num_wann, &
                                          num_bands, num_kpts, error, comm)
      if (allocated(error)) call prterr(error, stdout, stderr, comm)
    endif

    have_disentangled = .false.

    if (w90_calculation%restart .eq. ' ') then
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

    open (newunit=stdout, file=trim(seedname)//'.wout', status=trim(stat), position=trim(pos)) ! redefining from unit 6

    call w90_readwrite_write_header(physics%bohr_version_str, physics%constants_version_str1, &
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
    call w90_wannier90_readwrite_write(atom_data, band_plot, dis_control, dis_spheres, &
                                       fermi_energy_list, fermi_surface_plot, kpt_latt, &
                                       output_file, wvfn_read, wann_control, proj_input, &
                                       input_proj, real_space_ham, select_projection, kpoint_path, &
                                       transport, print_output, wannier_data, wannier_plot, &
                                       w90_extra_io, w90_calculation, real_lattice, &
                                       sitesym%symmetrize_eps, mp_grid, num_bands, num_kpts, &
                                       num_proj, num_wann, optimisation, cp_pp, gamma_only, &
                                       lsitesymmetry, w90_system%spinors, use_bloch_phases, &
                                       stdout)

    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

    if (.not. kmesh_info%explicit_nnkpts) call kmesh_get(kmesh_input, kmesh_info, print_output, &
                                                         kpt_latt, real_lattice, &
                                                         num_kpts, gamma_only, stdout, &
                                                         timer, error, comm)
    if (allocated(error)) then
      call prterr(error, stdout, stderr, comm)
    end if

    time2 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to get kmesh        ', time2 - time1, ' (sec)'

    call w90_wannier90_readwrite_memory_estimate(atom_data, kmesh_info, wann_control, input_proj, &
                                                 print_output, num_bands, num_kpts, num_proj, &
                                                 num_wann, optimisation, gamma_only, stdout)
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
  call w90_wannier90_readwrite_dist(atom_data, band_plot, dis_control, dis_spheres, dis_manifold, &
                                    exclude_bands, fermi_energy_list, fermi_surface_plot, &
                                    kmesh_input, kmesh_info, kpt_latt, output_file, wvfn_read, &
                                    wann_control, omega, input_proj, real_space_ham, &
                                    w90_system, transport, print_output, wannier_data, &
                                    wannier_plot, ws_region, w90_calculation, eigval, &
                                    real_lattice, sitesym%symmetrize_eps, mp_grid, &
                                    kpoint_path%num_points_first_segment, num_bands, num_kpts, &
                                    num_proj, num_wann, optimisation, eig_found, cp_pp, &
                                    gamma_only, have_disentangled, lhasproj, lsitesymmetry, &
                                    use_bloch_phases, error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  if (gamma_only .and. num_nodes > 1) then
    call set_error_fatal(error, 'Gamma point branch is serial only at the moment', comm)
    call prterr(error, stdout, stderr, comm)
  endif

  disentanglement = (num_bands > num_wann) ! for non-root

  if (w90_calculation%transport .and. transport%read_ht) goto 3003

  ! Sort out restarts
  if (w90_calculation%restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) then
      num_exclude_bands = 0
      if (allocated(exclude_bands)) num_exclude_bands = size(exclude_bands)
      call w90_readwrite_read_chkpt(dis_manifold, exclude_bands, kmesh_info, kpt_latt, &
                                    wannier_data, m_matrix, u_matrix, u_matrix_opt, real_lattice, &
                                    omega%invariant, mp_grid, num_bands, num_exclude_bands, &
                                    num_kpts, num_wann, checkpoint, have_disentangled, .false., &
                                    seedname, stdout, error, comm)
      if (allocated(error)) call prterr(error, stdout, stderr, comm)
    endif
    call w90_readwrite_chkpt_dist(dis_manifold, wannier_data, u_matrix, u_matrix_opt, &
                                  omega%invariant, num_bands, num_kpts, num_wann, checkpoint, &
                                  have_disentangled, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)

    if (lsitesymmetry) then
      call sitesym_read(sitesym, num_bands, num_kpts, num_wann, seedname, error, comm)
      ! update this to read on root and bcast - JRY
      if (allocated(error)) call prterr(error, stdout, stderr, comm)
    endif

    select case (w90_calculation%restart)
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
        call set_error_input(error, 'Value of checkpoint not recognised in wann_prog', comm)
        call prterr(error, stdout, stderr, comm)
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
    case default        ! for completeness... (it is already trapped in w90_wannier90_readwrite_read)
      call set_error_input(error, 'Value of restart not recognised in wann_prog', comm)
      call prterr(error, stdout, stderr, comm)
    end select
  endif

  if (w90_calculation%postproc_setup) then
    if (on_root) call kmesh_write(exclude_bands, kmesh_info, input_proj, print_output, kpt_latt, &
                                  real_lattice, num_kpts, num_proj, calc_only_A, &
                                  w90_system%spinors, seedname, timer)

    call kmesh_dealloc(kmesh_info, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)
    call w90_wannier90_readwrite_w90_dealloc(atom_data, band_plot, dis_spheres, dis_manifold, &
                                             exclude_bands, kmesh_input, kpt_latt, wann_control, &
                                             proj_input, input_proj, select_projection, &
                                             kpoint_path, wannier_data, wannier_plot, &
                                             w90_extra_io, eigval, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (on_root) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
#ifdef MPI
    call mpi_finalize(ierr)
#endif
    stop
  endif

  if (lsitesymmetry) then
    call sitesym_read(sitesym, num_bands, num_kpts, num_wann, seedname, error, comm)
    ! update this to read on root and bcast - JRY
    if (allocated(error)) call prterr(error, stdout, stderr, comm)
  endif

  call overlap_allocate(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, &
                        u_matrix, u_matrix_opt, kmesh_info%nntot, num_bands, num_kpts, num_wann, &
                        print_output%timing_level, timer, error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  call overlap_read(kmesh_info, select_projection, sitesym, a_matrix, m_matrix, m_matrix_local, &
                    m_matrix_orig, m_matrix_orig_local, u_matrix, u_matrix_opt, num_bands, &
                    num_kpts, num_proj, num_wann, print_output%timing_level, cp_pp, &
                    gamma_only, lsitesymmetry, use_bloch_phases, seedname, stdout, timer, error, &
                    comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  time1 = io_time()
  if (on_root) write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, &
    ' (sec)'

  have_disentangled = .false.

  if (disentanglement) then

    call dis_main(dis_control, dis_spheres, dis_manifold, kmesh_info, kpt_latt, sitesym, &
                  print_output, a_matrix, m_matrix, m_matrix_local, m_matrix_orig, &
                  m_matrix_orig_local, u_matrix, u_matrix_opt, eigval, real_lattice, &
                  omega%invariant, num_bands, num_kpts, num_wann, optimisation, gamma_only, &
                  lsitesymmetry, stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)

    have_disentangled = .true.
    time2 = io_time()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, &
      ' (sec)'
  endif

  if (on_root) then
    call w90_wannier90_readwrite_write_chkpt('postdis', exclude_bands, wannier_data, kmesh_info, &
                                             kpt_latt, num_kpts, dis_manifold, num_bands, &
                                             num_wann, u_matrix, u_matrix_opt, m_matrix, mp_grid, &
                                             real_lattice, omega%invariant, have_disentangled, &
                                             stdout, seedname)
  endif

1001 time2 = io_time()

  ! JJ workaround mpi_scatterv requirement that all arrays are valid *for all mpi procs*
  ! m_matrix* usually alloc'd in overlaps.F90, but not invariably, need to check here
  if (.not. allocated(m_matrix)) allocate (m_matrix(1, 1, 1, 1))

  if (.not. gamma_only) then
    call wann_main(atom_data, dis_manifold, exclude_bands, ham_logical, kmesh_info, kpt_latt, &
                   output_file, real_space_ham, wann_control, omega, sitesym, w90_system, &
                   print_output, wannier_data, ws_region, w90_calculation, ham_k, ham_r, m_matrix, &
                   u_matrix, u_matrix_opt, eigval, real_lattice, wannier_centres_translated, &
                   irvec, mp_grid, ndegen, shift_vec, nrpts, num_bands, num_kpts, num_proj, &
                   num_wann, optimisation, rpt_origin, band_plot%mode, transport%mode, &
                   have_disentangled, lsitesymmetry, seedname, stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)

  else
    call wann_main_gamma(atom_data, dis_manifold, exclude_bands, kmesh_info, kpt_latt, &
                         output_file, wann_control, omega, w90_system, print_output, wannier_data, &
                         m_matrix, u_matrix, u_matrix_opt, eigval, real_lattice, mp_grid, &
                         num_bands, num_kpts, num_wann, have_disentangled, &
                         real_space_ham%translate_home_cell, seedname, stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)

  end if

  time1 = io_time()
  if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time1 - time2, &
    ' (sec)'

  if (on_root) then
    call w90_wannier90_readwrite_write_chkpt('postwann', exclude_bands, wannier_data, kmesh_info, &
                                             kpt_latt, num_kpts, dis_manifold, num_bands, &
                                             num_wann, u_matrix, u_matrix_opt, m_matrix, mp_grid, &
                                             real_lattice, omega%invariant, have_disentangled, &
                                             stdout, seedname)
  endif

2002 continue
  if (on_root) then
    time2 = io_time()
  end if

  ! I call the routine always; the if statements to decide if/what to plot are inside the function
  call plot_main(atom_data, band_plot, dis_manifold, fermi_energy_list, fermi_surface_plot, &
                 ham_logical, kmesh_info, kpt_latt, output_file, wvfn_read, real_space_ham, &
                 kpoint_path, print_output, wannier_data, wannier_plot, ws_region, &
                 w90_calculation, ham_k, ham_r, m_matrix, u_matrix, u_matrix_opt, eigval, &
                 real_lattice, wannier_centres_translated, physics%bohr, irvec, mp_grid, ndegen, &
                 shift_vec, nrpts, num_bands, num_kpts, num_wann, rpt_origin, transport%mode, &
                 have_disentangled, lsitesymmetry, w90_system%spinors, seedname, stdout, timer, &
                 error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  if (on_root) then
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
    if (w90_calculation%transport) then
      time2 = io_time()

      call tran_main(atom_data, dis_manifold, fermi_energy_list, ham_logical, kpt_latt, &
                     output_file, real_space_ham, transport, print_output, wannier_data, &
                     ws_region, w90_calculation, ham_k, ham_r, u_matrix, u_matrix_opt, eigval, &
                     real_lattice, wannier_centres_translated, irvec, mp_grid, ndegen, shift_vec, &
                     nrpts, num_bands, num_kpts, num_wann, rpt_origin, band_plot%mode, &
                     have_disentangled, lsitesymmetry, seedname, stdout, timer, error, comm)
      if (allocated(error)) call prterr(error, stdout, stderr, comm)

      time1 = io_time()

      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (transport%read_ht) goto 4004
    end if
  endif

  call hamiltonian_dealloc(ham_logical, ham_k, ham_r, wannier_centres_translated, irvec, ndegen, &
                           error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  call overlap_dealloc(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, &
                       u_matrix, u_matrix_opt, error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  call kmesh_dealloc(kmesh_info, error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)
  call w90_wannier90_readwrite_w90_dealloc(atom_data, band_plot, dis_spheres, dis_manifold, &
                                           exclude_bands, kmesh_input, kpt_latt, wann_control, &
                                           proj_input, input_proj, select_projection, kpoint_path, &
                                           wannier_data, wannier_plot, w90_extra_io, eigval, &
                                           error, comm)
  if (allocated(error)) call prterr(error, stdout, stderr, comm)

  if (lsitesymmetry) then
    call sitesym_dealloc(sitesym, error, comm)
    if (allocated(error)) call prterr(error, stdout, stderr, comm)
  endif

4004 continue

  if (on_root) then
    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time(), ' (sec)'

    if (print_output%timing_level > 0) call io_print_timings(timer, stdout)

    write (stdout, *)
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'

    close (stdout)
    close (stderr, status='delete') ! this should not be unit 0
  endif

#ifdef MPI
  call mpi_finalize(ierr)
#endif

contains

  subroutine prterr(error, stdout, stderr, comm)
    use w90_error_base, only: code_remote
    use w90_comms, only: comms_no_sync_send, comms_no_sync_recv

    type(w90_error_type), allocatable, intent(in) :: error
    integer, intent(in) :: stderr, stdout
    type(w90comm_type), intent(in) :: comm

    type(w90_error_type), allocatable :: le ! unchecked error state for calls made in this routine
    integer :: ie ! global error value to be returned
    integer :: je ! error value on remote ranks
    integer :: j ! rank index
    integer :: failrank ! lowest rank reporting an error
    character(len=128) :: mesg ! only print 128 chars of error

    ie = 0
    mesg = 'not set'

    if (mpirank(comm) == 0) then
      ! fixme, report all failing ranks
      do j = mpisize(comm) - 1, 1, -1
        call comms_no_sync_recv(je, 1, j, le, comm)

        if (je /= code_remote .and. je /= 0) then
          failrank = j
          ie = je
          call comms_no_sync_recv(mesg, 128, j, le, comm)
        endif
      enddo
      ! if the error is on rank0
      if (error%code /= code_remote .and. error%code /= 0) then
        failrank = 0
        ie = error%code
        mesg = error%message
      endif

      !if (ie == 0) write (stderr, *) "logic error" ! to arrive here requires this

      write (stderr, *) 'Exiting.......'
      write (stderr, '(1x,a)') trim(mesg)
      write (stderr, '(1x,a,i0,a)') '(rank: ', failrank, ')'
      write (stdout, '(1x,a)') ' error encountered; check error .werr error log'

    else ! non 0 ranks
      je = error%code
      call comms_no_sync_send(je, 1, 0, le, comm)
      if (je /= code_remote .and. je /= 0) then
        mesg = error%message
        call comms_no_sync_send(mesg, 128, 0, le, comm)
      endif
    endif

#ifdef MPI
    call mpi_finalize(je) ! je overwritten here
#endif
    !call exit(ie) ! return true fail code (gnu extension)
    stop
  end subroutine prterr

end program wannier

