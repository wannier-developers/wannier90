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

  use w90_param_methods, only: param_write_header, param_read_chkpt, param_chkpt_dist
  use wannier_param_types
  use wannier_methods, only: param_read, param_w90_dealloc, param_write, &
    param_dist, param_memory_estimate, param_write_chkpt
#ifdef MPI
  use mpi_f08
#endif

  implicit none

  ! data from parameters module
  type(w90_calculation_type) :: w90_calcs
  ! Are we running postw90?
  !logical :: ispostw90 = .false.
  type(param_driver_type) :: driver
  type(postproc_type) :: pp_calc
  type(parameter_input_type) :: param_input
  logical :: eig_found
  type(param_plot_type) :: param_plot
  type(param_wannierise_type) :: param_wannierise
  logical :: lsitesymmetry = .false. ! RS: symmetry-adapted Wannier functions
  real(kind=dp) :: symmetrize_eps = 1.d-3
  type(wannier_data_type) :: wann_data
  type(param_hamiltonian_type) :: param_hamil
  type(param_kmesh_type) :: kmesh_data
  type(kmesh_info_type) :: kmesh_info
  type(k_point_type) :: k_points
  integer :: num_kpts
  type(disentangle_type) :: dis_data
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
  character(len=50) :: prog

  logical :: library !JJ
  type(sitesym_data) :: sym !JJ

! logical :: ham_have_setup = .false.
! logical :: have_translated = .false.
! logical :: use_translation = .false.
  type(ham_logical) :: hmlg

#ifdef MPI
  type(mpi_comm) :: comm = MPI_COMM_WORLD ! standalone executable case
#else
  integer :: comm
#endif

  call comms_setup(comm)

  library = .false.

  time0 = io_time()

  if (on_root) then
    prog = 'wannier90'
    call io_commandline(prog, dryrun)
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1, comm)
  call comms_bcast(seedname, len_seedname, comm)
  call comms_bcast(dryrun, 1, comm)

  if (on_root) then
    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.werr')
    call io_date(cdate, ctime)
    write (stdout, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

    call param_read(driver, w90_calcs, pp_calc, param_input, param_plot, &
                    param_wannierise, lsitesymmetry, symmetrize_eps, &
                    wann_data, param_hamil, kmesh_data, kmesh_info, &
                    k_points, num_kpts, dis_data, fermi_surface_data, &
                    fermi, tran, atoms, num_bands, num_wann, eigval, &
                    mp_grid, num_proj, select_proj, real_lattice, &
                    recip_lattice, spec_points, eig_found, .false., .false.)
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
      write (stdout, '(/,1x,a,i3,a/)') 'Running in parallel on ', num_nodes, ' CPUs'
    endif
    call param_write(driver, w90_calcs, param_input, param_plot, &
                     param_wannierise, lsitesymmetry, symmetrize_eps, &
                     wann_data, param_hamil, kmesh_data, k_points, num_kpts, &
                     dis_data, fermi_surface_data, fermi, tran, atoms, &
                     num_bands, num_wann, mp_grid, num_proj, select_proj, &
                     real_lattice, recip_lattice, spec_points)

    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

    if (.not. driver%explicit_nnkpts) call kmesh_get(recip_lattice, k_points%kpt_cart, &
                                                     param_input, kmesh_info, kmesh_data, num_kpts)
    time2 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to get kmesh        ', time2 - time1, ' (sec)'

    call param_memory_estimate(w90_calcs, param_input, param_wannierise, &
                               kmesh_data, kmesh_info, num_kpts, &
                               atoms, num_bands, num_wann, num_proj)
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
  call param_dist(driver, w90_calcs, pp_calc, param_input, param_plot, param_wannierise, &
                  lsitesymmetry, symmetrize_eps, wann_data, param_hamil, kmesh_data, kmesh_info, &
                  k_points, num_kpts, dis_data, fermi_surface_data, fermi, tran, atoms, num_bands, &
                  num_wann, eigval, mp_grid, num_proj, real_lattice, recip_lattice, eig_found, comm)
  if (param_input%gamma_only .and. num_nodes > 1) &
    call io_error('Gamma point branch is serial only at the moment')

  if (w90_calcs%transport .and. tran%read_ht) goto 3003

  ! Sort out restarts
  if (driver%restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) then
      call param_read_chkpt(driver%checkpoint, param_input, wann_data, kmesh_info, &
                            k_points, num_kpts, dis_data, num_bands, num_wann, &
                            u_matrix, u_matrix_opt, m_matrix, mp_grid, &
                            real_lattice, recip_lattice, .false.)
    endif
    call param_chkpt_dist(driver%checkpoint, param_input, wann_data, num_kpts, dis_data, &
                          num_bands, num_wann, u_matrix, u_matrix_opt, comm)

    if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
    if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

    select case (driver%restart)
    case ('default')    ! continue from where last checkpoint was written
      if (on_root) write (stdout, '(/1x,a)', advance='no') &
        'Resuming a previous Wannier90 calculation '
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
    if (on_root) call kmesh_write(recip_lattice, param_input, kmesh_info, num_kpts, kmesh_data, &
                                  num_proj, k_points%kpt_latt, real_lattice, pp_calc%only_A)
    call kmesh_dealloc(kmesh_info)
    call param_w90_dealloc(param_input, param_plot, param_wannierise, &
                           wann_data, kmesh_data, k_points, dis_data, &
                           atoms, eigval, spec_points)
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (on_root) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
    call comms_end
    stop
  endif

  if (lsitesymmetry) call sitesym_read(num_bands, num_wann, num_kpts, sym)   ! update this to read on root and bcast - JRY
  if (lsitesymmetry) sym%symmetrize_eps = symmetrize_eps ! for the time being, copy value from w90_parameters  (JJ)

  call overlap_allocate(u_matrix, m_matrix_local, m_matrix, u_matrix_opt, a_matrix, &
                        m_matrix_orig_local, m_matrix_orig, param_input%timing_level, &
                        kmesh_info%nntot, num_kpts, num_wann, num_bands, w90_calcs%disentanglement)

  call overlap_read(lsitesymmetry, m_matrix_orig_local, m_matrix_local, param_input, w90_calcs, &
                    u_matrix_opt, m_matrix_orig, a_matrix, m_matrix, u_matrix, select_proj, &
                    num_proj, kmesh_info, num_kpts, num_wann, num_bands, sym, comm)

  time1 = io_time()
  if (on_root) write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, &
    ' (sec)'

  param_input%have_disentangled = .false.

  if (w90_calcs%disentanglement) then

    call dis_main(num_bands, num_kpts, num_wann, recip_lattice, eigval, a_matrix, m_matrix, &
                  m_matrix_local, m_matrix_orig, m_matrix_orig_local, u_matrix, u_matrix_opt, &
                  dis_data, kmesh_info, k_points, param_input, num_nodes, my_node_id, on_root, &
                  lsitesymmetry, sym, comm)

    param_input%have_disentangled = .true.
    time2 = io_time()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, &
      ' (sec)'
  endif

  if (on_root) then
    call param_write_chkpt('postdis', param_input, wann_data, kmesh_info, k_points, num_kpts, &
                           dis_data, num_bands, num_wann, u_matrix, u_matrix_opt, m_matrix, &
                           mp_grid, real_lattice, recip_lattice)
  endif
!~  call param_write_um

1001 time2 = io_time()

  ! JJ hack to workaround mpi_scatterv requirement that all arrays are valid *for all mpi procs*
  ! m_matrix* usually alloc'd in overlaps.F90, but not always
  if (.not. allocated(m_matrix)) allocate (m_matrix(1, 1, 1, 1)) !JJ temporary workaround to avoid runtime check failure

  if (.not. param_input%gamma_only) then
    call wann_main(num_wann, param_wannierise, kmesh_info, param_input, u_matrix, m_matrix, &
                   num_kpts, real_lattice, num_proj, wann_data, k_points, num_bands, u_matrix_opt, &
                   eigval, dis_data, recip_lattice, atoms, lsitesymmetry, stdout, mp_grid, &
                   w90_calcs, tran%mode, param_hamil, sym, ham_r, irvec, shift_vec, ndegen, nrpts, &
                   rpt_origin, wannier_centres_translated, hmlg, ham_k, comm)
  else
    call wann_main_gamma(num_wann, param_wannierise, kmesh_info, param_input, u_matrix, m_matrix, &
                         num_kpts, real_lattice, wann_data, num_bands, u_matrix_opt, eigval, &
                         dis_data%lwindow, recip_lattice, atoms, k_points, dis_data, mp_grid, &
                         stdout, comm)
  end if

  time1 = io_time()
  if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time1 - time2, &
    ' (sec)'

  if (on_root) then
    call param_write_chkpt('postwann', param_input, wann_data, kmesh_info, k_points, num_kpts, &
                           dis_data, num_bands, num_wann, u_matrix, u_matrix_opt, m_matrix, &
                           mp_grid, real_lattice, recip_lattice)
  endif

2002 continue
  if (on_root) then
    ! I call the routine always; the if statements to decide if/what to plot are inside the function
    time2 = io_time()

    call plot_main(num_kpts, w90_calcs, k_points, param_input, param_plot, real_lattice, &
                   num_wann, kmesh_info, m_matrix, recip_lattice, wann_data, atoms, param_hamil, &
                   dis_data, u_matrix_opt, eigval, u_matrix, lsitesymmetry, num_bands, mp_grid, &
                   tran%mode, fermi, fermi_surface_data, spec_points, ham_r, irvec, shift_vec, ndegen, &
                   nrpts, rpt_origin, wannier_centres_translated, hmlg, ham_k)
    time1 = io_time()

    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
    if (w90_calcs%transport) then
      time2 = io_time()

      call tran_main(tran, param_input, w90_calcs, num_wann, real_lattice, recip_lattice, &
                     wann_data, atoms, param_hamil, dis_data, u_matrix_opt, k_points, eigval, &
                     u_matrix, lsitesymmetry, num_bands, num_kpts, mp_grid, fermi, ham_r, irvec, &
                     shift_vec, ndegen, nrpts, rpt_origin, wannier_centres_translated, hmlg, ham_k)
      time1 = io_time()

      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (tran%read_ht) goto 4004
    end if
  endif

  call tran_dealloc()
  call hamiltonian_dealloc(ham_r, irvec, ndegen, wannier_centres_translated, hmlg, ham_k)
  call overlap_dealloc(m_matrix_orig_local, m_matrix_local, u_matrix_opt, a_matrix, &
                       m_matrix_orig, m_matrix, u_matrix)
  call kmesh_dealloc(kmesh_info)
  call param_w90_dealloc(param_input, param_plot, param_wannierise, &
                         wann_data, kmesh_data, k_points, dis_data, &
                         atoms, eigval, spec_points)
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

