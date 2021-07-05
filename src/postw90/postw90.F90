!-*- mode: F90 -*-!
!------------------------------------------------------------!
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

program postw90
  !! The postw90 program
  use w90_constants, only: dp, eps6, pw90_physical_constants
  use w90_param_types
  use pw90_parameters
  use w90_param_methods, only: param_read_chkpt, param_write_header
  use pw90_param_methods
  use w90_io
  use w90_kmesh
  use w90_comms, only: comms_end, comms_bcast, comms_barrier, w90commtype, mpirank, mpisize
  use w90_postw90_common, only: pw90common_wanint_setup, pw90common_wanint_get_kpoint_file, &
    pw90common_wanint_param_dist, pw90common_wanint_data_dist, kpoint_dist_type, wigner_seitz_type

  ! These modules deal with the interpolation of specific physical properties
  !
  use w90_dos
  use w90_berry, only: berry_main
  use w90_gyrotropic
  use w90_spin
  use w90_kpath
  use w90_kslice

  use w90_boltzwann
  use w90_geninterp

  use w90_ws_distance, only: ws_distance_type

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

  type(pw90_physical_constants) :: physics
  integer       :: nkp, len_seedname
  integer       :: stdout
  character(len=50) :: seedname
  logical       :: have_gamma
  real(kind=dp) :: time0, time1, time2
  character(len=9) :: stat, pos
  logical :: wpout_found, werr_found, dryrun
  character(len=50) :: prog
  character(len=20) :: checkpoint
  ! this is a dummy that is not used in postw90, DO NOT use
  complex(kind=dp), allocatable :: m_matrix(:, :, :, :)

!
  complex(kind=dp), allocatable :: HH_R(:, :, :) !  <0n|r|Rm>
  !! $$\langle 0n | H | Rm \rangle$$

  complex(kind=dp), allocatable :: AA_R(:, :, :, :) ! <0n|r|Rm>
  !! $$\langle 0n |  \hat{r} | Rm \rangle$$

  complex(kind=dp), allocatable :: BB_R(:, :, :, :) ! <0|H(r-R)|R>
  !! $$\langle 0n | H(\hat{r}-R) | Rm \rangle$$

  complex(kind=dp), allocatable :: CC_R(:, :, :, :, :) ! <0|r_alpha.H(r-R)_beta|R>
  !! $$\langle 0n | \hat{r}_{\alpha}.H.(\hat{r}- R)_{\beta} | Rm \rangle$$

  complex(kind=dp), allocatable :: FF_R(:, :, :, :, :) ! <0|r_alpha.(r-R)_beta|R>
  !! $$\langle 0n | \hat{r}_{\alpha}.(\hat{r}-R)_{\beta} | Rm \rangle$$

  complex(kind=dp), allocatable :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>
  !! $$\langle 0n | \sigma_{x,y,z} | Rm \rangle$$

  complex(kind=dp), allocatable :: SR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.(r-R)_alpha|Rm>
  !! $$\langle 0n | \sigma_{x,y,z}.(\hat{r}-R)_{\alpha}  | Rm \rangle$$

  complex(kind=dp), allocatable :: SHR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.H.(r-R)_alpha|Rm>
  !! $$\langle 0n | \sigma_{x,y,z}.H.(\hat{r}-R)_{\alpha}  | Rm \rangle$$

  complex(kind=dp), allocatable :: SH_R(:, :, :, :) ! <0n|sigma_x,y,z.H|Rm>
  !! $$\langle 0n | \sigma_{x,y,z}.H  | Rm \rangle$$
!
  ! w90_parameters stuff
  type(parameter_input_type) :: param_input
  type(wannier_data_type) :: wann_data
  type(param_kmesh_type) :: kmesh_data
  type(kmesh_info_type) :: kmesh_info
  type(k_point_type) :: k_points
  integer :: num_kpts !BGS put in k_point_type?
  type(disentangle_manifold_type) :: dis_window
  type(fermi_data_type) :: fermi
  type(atom_data_type) :: atoms

  integer :: num_bands
  !! Number of bands

  integer :: num_wann
  !! number of wannier functions

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states
  !BGS disentangle, hamiltonian, a wannierise print, and postw90/get_oper
  real(kind=dp), allocatable :: eigval(:, :)

  ! u_matrix_opt in postw90 only for generation of v_matrix
  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states
  complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)

  ! u_matrix gives the unitary rotations from the optimal subspace to the
  ! optimally smooth states.
  ! m_matrix we store here, becuase it is needed for restart of wannierise
  complex(kind=dp), allocatable :: u_matrix(:, :, :)

  integer :: mp_grid(3)
  !! Dimensions of the Monkhorst-Pack grid

  integer :: num_proj

  real(kind=dp) :: real_lattice(3, 3)

  !parameters derived from input
  real(kind=dp) :: recip_lattice(3, 3)

  type(special_kpoints_type) :: spec_points
  ! end w90_parameters
  ! data from pw90_parameters
  type(pw90_calculation_type), save :: pw90_calcs

  logical, save :: eig_found ! used to control broadcast of eigval

  type(postw90_oper_type), save :: postw90_oper
  type(postw90_common_type), save :: pw90_common
  type(postw90_spin_type), save :: pw90_spin
  type(postw90_ham_type), save :: pw90_ham
  type(kpath_type), save :: kpath
  type(kslice_type), save :: kslice

  ! module  d o s
  ! No need to save 'dos_plot', only used here (introduced 'dos_task')
  logical          :: dos_plot

  type(dos_plot_type), save :: dos_data
  type(berry_type), save :: berry
  type(spin_hall_type), save :: spin_hall
  type(gyrotropic_type), save :: gyrotropic
  type(geninterp_type), save :: geninterp
  type(boltzwann_type), save :: boltz
  ! end pw90_parameters
  ! from postw90_common
  complex(kind=dp), allocatable :: v_matrix(:, :, :)
  type(wigner_seitz_type) :: ws_vec
  type(kpoint_dist_type) :: kpt_dist

  ! local vars
  integer :: my_node_id, num_nodes, ierr
  logical :: on_root = .false.
  type(w90commtype) :: comm
  type(ws_distance_type) :: ws_distance
  type(pw90_extra_io_type) :: write_data
  real(kind=dp) :: omega_invariant

#ifdef MPI
  comm%comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  if (ierr .ne. 0) call io_error('MPI initialisation error', stdout, seedname)  ! JJ, fixme, what are stdout, seedname here?  unassigned!
#endif

  my_node_id = mpirank(comm)
  num_nodes = mpisize(comm)
  if (my_node_id == 0) on_root = .true.

  if (on_root) then
    time0 = io_time()
    prog = 'postw90'
    call io_commandline(prog, dryrun, seedname)
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1, stdout, seedname, comm)
  call comms_bcast(seedname, len_seedname, stdout, seedname, comm)
  call comms_bcast(dryrun, 1, stdout, seedname, comm)

  if (on_root) then
    ! If an error file (generated by postw90) exists, I delete it
    ! Note: I do it only for the error file generated from the root node;
    ! If error files generated by other nodes exist, I don't do anything for them
    inquire (file=trim(seedname)//'.node_00000.werr', exist=werr_found)
    if (werr_found) then
      stdout = io_file_unit()
      open (unit=stdout, file=trim(seedname)//'.node_00000.werr', status='old', position='append')
      close (stdout, status='delete')
    end if

    inquire (file=trim(seedname)//'.wpout', exist=wpout_found)
    if (wpout_found) then
      stat = 'old'
    else
      stat = 'replace'
    endif
    pos = 'append'

    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.wpout', status=trim(stat), position=trim(pos))

    call param_write_header(physics%bohr_version_str, physics%constants_version_str1, &
                            physics%constants_version_str2, stdout)
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
  end if

  ! Read onto the root node all the input parameters from seendame.win,
  ! as well as the energy eigenvalues on the ab-initio q-mesh from seedname.eig
  !
  if (on_root) then
    call param_postw90_read(param_input, wann_data, kmesh_data, k_points, num_kpts, &
                            dis_window, fermi, atoms, num_bands, num_wann, eigval, &
                            mp_grid, real_lattice, recip_lattice, spec_points, &
                            pw90_calcs, postw90_oper, pw90_common, pw90_spin, &
                            pw90_ham, kpath, kslice, dos_data, berry, &
                            spin_hall, gyrotropic, geninterp, boltz, eig_found, write_data, &
                            physics%bohr, stdout, seedname)
    call param_postw90_write(param_input, fermi, atoms, num_wann, &
                             real_lattice, recip_lattice, spec_points, &
                             pw90_calcs, postw90_oper, pw90_common, &
                             pw90_spin, kpath, kslice, dos_data, berry, &
                             gyrotropic, geninterp, boltz, write_data, stdout)
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') &
      'Time to read parameters  ', time1 - time0, ' (sec)'

    if (.not. pw90_common%effective_model) then
      ! Check if the q-mesh includes the gamma point
      !
      have_gamma = .false.
      do nkp = 1, num_kpts
        if (all(abs(k_points%kpt_latt(:, nkp)) < eps6)) have_gamma = .true.
      end do
      if (.not. have_gamma) write (stdout, '(1x,a)') &
        'Ab-initio does not include Gamma. Interpolation may be incorrect!!!'
      !
      ! Need nntot, wb, and bk to evaluate WF matrix elements of
      ! the position operator in reciprocal space. Also need
      ! nnlist to compute the additional matrix elements entering
      ! the orbital magnetization
      !
      call kmesh_get(kmesh_data, kmesh_info, param_input, k_points%kpt_cart, recip_lattice, &
                     num_kpts, seedname, stdout)
      time2 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') &
        'Time to get kmesh        ', time2 - time1, ' (sec)'
    endif

    ! GP, May 10, 2012: for the moment I leave this commented
    ! since we need first to tune that routine so that it doesn't
    ! print the memory information related to wannier90.x.
    ! Note that the code for the memory estimation for the
    ! Boltzwann routine is already there.
    !       call param_memory_estimate
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

  ! We now distribute a subset of the parameters to the other nodes
  !
  call pw90common_wanint_param_dist(param_input, kmesh_info, k_points, num_kpts, dis_window, &
                                    fermi, num_bands, num_wann, eigval, mp_grid, &
                                    real_lattice, recip_lattice, pw90_calcs, &
                                    pw90_common, pw90_spin, pw90_ham, kpath, kslice, &
                                    dos_data, berry, spin_hall, gyrotropic, geninterp, &
                                    boltz, eig_found, stdout, seedname, comm)

  if (.not. pw90_common%effective_model) then
    !
    ! Read files seedname.chk (overlap matrices, unitary matrices for
    ! both disentanglement and maximal localization, etc.)
    !
    if (on_root) then
      call param_read_chkpt(dis_window, kmesh_info, k_points, param_input, wann_data, m_matrix, &
                            u_matrix, u_matrix_opt, real_lattice, recip_lattice, omega_invariant, &
                            mp_grid, num_bands, num_kpts, num_wann, checkpoint, .true., seedname, &
                            stdout)
    endif
    !
    ! Distribute the information in the um and chk files to the other nodes
    !
    ! Ivo: For interpolation purposes we do not need u_matrix_opt and
    !      u_matrix separately, only their product v_matrix, and this
    !      is what is distributed now
    !
    call pw90common_wanint_data_dist(num_wann, num_kpts, num_bands, u_matrix_opt, u_matrix, &
                                     dis_window, param_input, wann_data, pw90_common, &
                                     v_matrix, stdout, seedname, comm)
    !
  end if

  ! Read list of k-points in irreducible BZ and their weights
  !
  ! Should this be done on root node only?
  !
  if (berry%wanint_kpoint_file) call pw90common_wanint_get_kpoint_file(kpt_dist, stdout, &
                                                                       seedname, comm)

  ! Setup a number of common variables for all interpolation tasks
  !
  call pw90common_wanint_setup(num_wann, param_input, real_lattice, mp_grid, pw90_common, &
                               ws_vec, stdout, seedname, comm)

  if (on_root) then
    time1 = io_time()
    write (stdout, '(/1x,a25,f11.3,a)') &
      'Time to read and process .chk    ', time1 - time2, ' (sec)'
  endif
  !
  ! Now perform one or more of the following tasks

  ! ---------------------------------------------------------------
  ! Density of states calculated using a uniform interpolation mesh
  ! ---------------------------------------------------------------
  !
  if (pw90_calcs%dos .and. index(dos_data%task, 'dos_plot') > 0) then
    call dos_main(berry, dis_window, dos_data, kpt_dist, k_points, param_input, pw90_common, &
                  pw90_ham, postw90_oper, pw90_spin, wann_data, ws_distance, ws_vec, HH_R, SS_R, &
                  u_matrix, v_matrix, eigval, real_lattice, recip_lattice, mp_grid, num_bands, &
                  num_kpts, num_wann, seedname, stdout, comm)
  endif

! find_fermi_level commented for the moment in dos.F90
!  if(dos .and. index(dos_task,'find_fermi_energy')>0) call find_fermi_level

  ! --------------------------------------------------------------------
  ! Bands, Berry curvature, or orbital magnetization plot along a k-path
  ! --------------------------------------------------------------------
  if (pw90_calcs%kpath) then
    call k_path(berry, dis_window, fermi, kmesh_info, kpath, k_points, param_input, postw90_oper, &
                pw90_common, pw90_ham, pw90_spin, spec_points, spin_hall, wann_data, ws_vec, &
                ws_distance, AA_R, BB_R, CC_R, HH_R, SH_R, SHR_R, SR_R, SS_R, v_matrix, u_matrix, &
                physics%bohr, eigval, real_lattice, recip_lattice, mp_grid, num_wann, num_bands, &
                num_kpts, seedname, stdout, comm)
  end if

  ! ---------------------------------------------------------------------------
  ! Bands, Berry curvature, or orbital magnetization plot on a slice in k-space
  ! ---------------------------------------------------------------------------
  !
  if (pw90_calcs%kslice) then

    call k_slice(berry, dis_window, fermi, kmesh_info, k_points, kslice, param_input, pw90_common, &
                 pw90_ham, postw90_oper, pw90_spin, spin_hall, wann_data, ws_distance, ws_vec, &
                 AA_R, BB_R, CC_R, HH_R, SH_R, SHR_R, SR_R, SS_R, v_matrix, u_matrix, &
                 physics%bohr, eigval, real_lattice, recip_lattice, mp_grid, num_bands, num_kpts, &
                 num_wann, seedname, stdout, comm)
  end if

  ! --------------------
  ! Spin magnetic moment
  ! --------------------
  !
  if (pw90_calcs%spin_moment) then
    call spin_get_moment(dis_window, fermi, kpt_dist, k_points, param_input, pw90_common, &
                         postw90_oper, pw90_spin, wann_data, ws_distance, ws_vec, HH_R, SS_R, &
                         u_matrix, v_matrix, eigval, real_lattice, recip_lattice, mp_grid, &
                         num_wann, num_bands, num_kpts, berry%wanint_kpoint_file, seedname, &
                         stdout, comm)
  end if

  ! -------------------------------------------------------------------
  ! dc Anomalous Hall conductivity and eventually (if 'mcd' string also
  ! present in addition to 'ahe', e.g., 'ahe+mcd') dichroic optical
  ! conductivity, both calculated on the same (adaptively-refined) mesh
  ! -------------------------------------------------------------------
  !
  ! ---------------------------------------------------------------
  ! Absorptive dichroic optical conductivity & JDOS on uniform mesh
  ! ---------------------------------------------------------------
  !
  ! -----------------------------------------------------------------
  ! Absorptive ordinary optical conductivity & JDOS on a uniform mesh
  ! -----------------------------------------------------------------
  !
  ! -----------------------------------------------------------------
  ! Orbital magnetization
  ! -----------------------------------------------------------------
  !
  if (pw90_calcs%berry) then
    call berry_main(berry, dis_window, fermi, kmesh_info, kpt_dist, k_points, param_input, &
                    pw90_common, pw90_ham, postw90_oper, pw90_spin, physics, spin_hall, wann_data, &
                    ws_distance, ws_vec, AA_R, BB_R, CC_R, HH_R, SH_R, SHR_R, SR_R, SS_R, &
                    u_matrix, v_matrix, eigval, real_lattice, recip_lattice, mp_grid, num_wann, &
                    num_kpts, num_bands, seedname, stdout, comm)
  end if
  ! -----------------------------------------------------------------
  ! Boltzmann transport coefficients (BoltzWann module)
  ! -----------------------------------------------------------------
  !
  if (on_root) then
    time1 = io_time()
  endif

  if (pw90_calcs%geninterp) then
    call geninterp_main(dis_window, geninterp, k_points, param_input, pw90_common, pw90_ham, &
                        wann_data, ws_distance, ws_vec, HH_R, v_matrix, u_matrix, eigval, &
                        real_lattice, recip_lattice, mp_grid, num_bands, num_kpts, num_wann, &
                        seedname, stdout, comm)
  end if

  if (pw90_calcs%boltzwann) then
    call boltzwann_main(boltz, dis_window, dos_data, k_points, param_input, pw90_common, pw90_ham, &
                        postw90_oper, pw90_spin, physics, wann_data, ws_distance, ws_vec, HH_R, &
                        SS_R, v_matrix, u_matrix, eigval, real_lattice, recip_lattice, mp_grid, &
                        num_wann, num_bands, num_kpts, seedname, stdout, comm)
  end if

  if (pw90_calcs%gyrotropic) then
    call gyrotropic_main(berry, dis_window, fermi, gyrotropic, kmesh_info, k_points, param_input, &
                         pw90_common, pw90_ham, postw90_oper, physics, wann_data, ws_vec, &
                         ws_distance, AA_R, BB_R, CC_R, HH_R, SS_R, u_matrix, v_matrix, &
                         eigval, real_lattice, recip_lattice, mp_grid, num_bands, num_kpts, &
                         num_wann, seedname, stdout, comm)
  endif

  if (on_root .and. pw90_calcs%boltzwann) then
    time2 = io_time()
    write (stdout, '(/1x,a,f11.3,a)') &
      'Time for BoltzWann (Boltzmann transport) ', time2 - time1, ' (sec)'
  endif

  ! I put a barrier here before calling the final time printing routines,
  ! just to be sure that all processes have arrived here.
  call comms_barrier(comm)
  if (on_root) then
    write (stdout, '(/,1x,a25,f11.3,a)') &
      'Total Execution Time     ', io_time(), ' (sec)'
    if (param_input%timing_level > 0) call io_print_timings(stdout)
    write (stdout, *)
    write (stdout, '(/,1x,a)') 'All done: postw90 exiting'
    close (stdout)
  end if

  call comms_end

end program postw90

