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

#ifdef MPI08
  use mpi_f08
#endif
#ifdef MPI90
  use mpi
#endif

  use w90_library
  use w90_library_extra ! for input_reader_special, overlaps, etc

  use w90_comms, only: w90_comm_type, comms_sync_error
  use w90_io, only: io_commandline, io_date, io_time, prterr
  use w90_sitesym, only: sitesym_read
  use w90_error, only: w90_error_type, set_error_input

  implicit none

  integer, parameter :: dp = kind(0.d0)

  character(len=:), allocatable :: seedname, progname, cpstatus
  character(len=:), pointer :: restart
  complex(kind=dp), allocatable :: m_matrix_loc(:, :, :, :)
  complex(kind=dp), allocatable :: u_matrix(:, :, :)
  complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)
  real(kind=dp), allocatable :: eigval(:, :)
  integer, allocatable :: dist_k(:)
  integer :: mpisize, rank, ierr, nkl
  integer, pointer :: nb, nk, nw, nn
  integer :: stdout, stderr
  logical, pointer :: pp
  logical :: ld, lovlp, ldsnt, lwann, lplot, ltran, need_eigvals
  type(lib_common_type), target :: common_data
  type(w90_error_type), allocatable :: error
  character(len=9) :: cdate, ctime

  pp => common_data%w90_calculation%postproc_setup
  restart => common_data%w90_calculation%restart
  nw => common_data%num_wann
  nb => common_data%num_bands
  nk => common_data%num_kpts
  nn => common_data%kmesh_info%nntot

  progname = 'wannier90' ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=91442
  call io_commandline(progname, ld, pp, seedname)

#ifdef MPI
  call mpi_init(ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: mpi_init() returned an error!'
    stop
  endif
  call mpi_comm_rank(mpi_comm_world, rank, ierr) ! the type of comm_world depends on interface used
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: mpi_comm_rank() returned an error!'
    stop
  endif
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: mpi_comm_size() returned an error!'
    stop
  endif
  call w90_set_comm(common_data, mpi_comm_world)
#else
  rank = 0
  mpisize = 1
#endif

  ! open main output file
  if (rank == 0) open (newunit=stdout, file=seedname//'.wout', status="replace")

  ! open main error file
  ! call w90_get_fortran_stderr(stderr) !alternative for terminal output
  if (rank == 0) open (newunit=stderr, file=seedname//'.werr', status="replace")

  call io_date(cdate, ctime)
  if (rank == 0) write (stderr, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

  ! read key parameters from .win file
  call input_reader_special(common_data, seedname, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! read all remaining parameters from .win file
  call w90_input_reader(common_data, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! write useful info (includes jazzy header info)
  call w90_print_info(common_data, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! special branch for writing nnkp file
  ! exit immediately after writing the nnkp file
  if (pp) then
    call write_kmesh(common_data, stdout, stderr, ierr) ! only active on rank 0
    if (ierr /= 0) stop
    if (rank == 0) close (unit=stderr, status='delete')
    if (rank == 0) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (rank == 0) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
#ifdef MPI
    call mpi_finalize(ierr)
#endif
    stop
  endif

  ! test mpi error handling using "unlucky" input token
  if (rank == -common_data%print_output%timing_level) then
    call set_error_input(error, 'received unlucky_rank', common_data%comm)
  else
    ! this is necessary since non-root may never enter an mpi collective if root has exited here
    call comms_sync_error(common_data%comm, error, 0)
  endif
  if (allocated(error)) then ! applies (is t) for all ranks now
    call prterr(error, ierr, stdout, stderr, common_data%comm)
#ifdef MPI
    call mpi_finalize(ierr) ! let's be nice
#endif
    stop
  endif
  ! end unlucky code

  ! setup kpoint distribution
  allocate (dist_k(nk), stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to allocate dist_k array!'
    stop
  endif
  ! get a basic k-point/rank distribution
  call w90_distribute_kpts(common_data, nk, mpisize, dist_k, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! copy distribution to library
  call set_kpoint_distribution(common_data, dist_k, stdout, stderr, ierr)
  if (ierr /= 0) stop

  ! setup SAWF data
  if (common_data%lsitesymmetry) then
    call sitesym_read(common_data%sitesym, nb, nk, nw, seedname, error, common_data%comm) ! (not a library call)
    if (allocated(error)) then
      write (stderr, *) 'Wannier90: failed to setup symmetry!'
      deallocate (error)
      stop
    endif
  endif

  call w90_get_nn(common_data, nn, stdout, stderr, ierr)
  nkl = count(dist_k == rank) ! number of kpoints this rank
  !write (*, *) 'rank, nw, nb, nk, nn, nk(rank): ', rank, nw, nb, nk, nn, nkl

  allocate (m_matrix_loc(nb, nb, nn, nkl), stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to allocate m_matrix_loc!'
    stop
  endif
  call w90_set_m_local(common_data, m_matrix_loc)  ! we don't need global m

  allocate (u_matrix(nw, nw, nk), stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to allocate u_matrix!'
    stop
  endif
  call w90_set_u_matrix(common_data, u_matrix)

  allocate (u_matrix_opt(nb, nw, nk), stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to allocate u_matrix_opt!'
    stop
  endif
  call w90_set_u_opt(common_data, u_matrix_opt)

! restart system
  lovlp = .true.
  ldsnt = .true.
  lwann = .true.
  lplot = .true.
  ltran = .false.

  if (restart == '') then
    if (rank == 0) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else
    cpstatus = ''
    call read_chkpt(common_data, cpstatus, stdout, stderr, ierr)
    if (ierr /= 0) stop

    if (restart == 'wannierise' .or. (restart == 'default' .and. cpstatus == 'postdis')) then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .true.
      lplot = .true.
      ltran = .false.
    elseif (restart == 'plot' .or. (restart == 'default' .and. cpstatus == 'postwann')) then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .false.
      lplot = .true.
      ltran = .false.
    elseif (restart == 'transport') then
      if (rank == 0) write (stdout, '(1x,a/)') 'Restarting Wannier90 from transport routines ...'
      lovlp = .false.
      ldsnt = .false.
      lwann = .false.
      lplot = .false.
      ltran = .true.
      !else
      ! illegitimate restart choice, should declaim the acceptable choices
    endif
  endif
  ltran = (ltran .or. common_data%w90_calculation%transport)
  ldsnt = (ldsnt .and. (nw < nb)) ! disentanglement only needed if space reduced

  ! circumstances where eigenvalues are needed are a little overcomplicated
  need_eigvals = .false.
  need_eigvals = common_data%w90_calculation%bands_plot
  need_eigvals = (need_eigvals .or. common_data%w90_calculation%fermi_surface_plot)
  need_eigvals = (need_eigvals .or. common_data%output_file%write_hr)
  need_eigvals = (need_eigvals .or. common_data%output_file%write_tb)
  need_eigvals = (need_eigvals .or. ldsnt) ! disentanglement anyway requires evals

  if (need_eigvals) then
    allocate (eigval(nb, nk), stat=ierr)
    if (ierr /= 0) then
      write (stderr, *) 'Wannier90: failed to allocate eigval array!'
      stop
    endif
    call read_eigvals(common_data, eigval, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call w90_set_eigval(common_data, eigval)
  endif

  ! ends setup

  if (lovlp) then
    call overlaps(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (ldsnt) then
    call w90_disentangle(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(common_data, 'postdis', stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lwann) then
    call w90_project_overlap(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call w90_wannierise(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
    call write_chkpt(common_data, 'postwann', stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (lplot) then
    call w90_plot(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  if (ltran) then
    call w90_transport(common_data, stdout, stderr, ierr)
    if (ierr /= 0) stop
  endif

  ! cleanup

  if (need_eigvals) then
    deallocate (eigval, stat=ierr)
    if (ierr /= 0) then
      write (stderr, *) 'Wannier90: failed to deallocate eigval array!'
      stop
    endif
  endif
  deallocate (dist_k, stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to deallocate dist_k array!'
    stop
  endif
  deallocate (m_matrix_loc, stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to deallocate m_matrix_loc!'
    stop
  endif
  deallocate (u_matrix, stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to deallocate u_matrix!'
    stop
  endif
  deallocate (u_matrix_opt, stat=ierr)
  if (ierr /= 0) then
    write (stderr, *) 'Wannier90: failed to deallocate u_matrix_opt!'
    stop
  endif

  call print_times(common_data, stdout)

  if (rank == 0) then
    close (unit=stderr, status='delete')
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'
    close (unit=stdout)
  endif

#ifdef MPI
  call mpi_finalize(ierr)
#endif
end program wannier
