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

module w90_libv1_types
  ! global type instances (of w90_types) to allow legacy library

  use w90_types
  use w90_error

  implicit none

  public

  integer, allocatable, save :: exclude_bands(:)
  integer, save :: num_bands !! Number of bands
  integer, save :: num_kpts
  integer, save :: num_wann !! number of wannier functions
  integer, save :: optimisation

  type(atom_data_type), save :: atoms
  type(dis_manifold_type), save :: dis_window
  type(kmesh_info_type), save :: kmesh_info
  type(kmesh_input_type), save :: kmesh_data
  type(print_output_type), save :: verbose
  type(proj_input_type), save :: input_proj
  type(w90_system_type), save :: system
  type(wannier_data_type), save :: wann_data
  type(ws_region_type) :: ws_region

  real(kind=dp), allocatable, save :: fermi_energy_list(:)
  real(kind=dp), allocatable, save :: kpt_latt(:, :) !! kpoints in lattice vecs

  logical, save :: cp_pp, calc_only_A
  logical, save :: gamma_only
  logical, save :: have_disentangled
  logical, save :: use_bloch_phases

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states
  !BGS disentangle, hamiltonian, a wannierise print, and postw90/get_oper
  real(kind=dp), allocatable, save :: eigval(:, :)

  !BGS u_matrix_opt in postw90 only for generation of v_matrix
  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states
  complex(kind=dp), allocatable, save :: u_matrix_opt(:, :, :)

  ! u_matrix gives the unitary rotations from the optimal subspace to the
  ! optimally smooth states.
  ! m_matrix we store here, becuase it is needed for restart of wannierise
  complex(kind=dp), allocatable, save :: u_matrix(:, :, :)

  integer, save :: mp_grid(3)
  !! Dimensions of the Monkhorst-Pack grid

  integer, save :: num_proj
  !BGS used by stuff in driver/kmesh/wannier - keep separate or duplicate?

  real(kind=dp), save :: real_lattice(3, 3)

  !parameters derived from input
  !real(kind=dp), save :: recip_lattice(3, 3)

  type(kpoint_path_type), save :: spec_points

contains

  subroutine prterr(error, stdout)
    type(w90_error_type), intent(in) :: error
    integer, intent(in) :: stdout
    !write (stdout, *) "ERROR CODE: ", error%code
    !write (stdout, *) "ERROR CONDITION: ", error%message
    ! rep old behaviour for a moment while testing...
    write (stdout, *) 'Exiting.......'
    write (stdout, *) error%message
    stop
    !call exit(error%code)
  end subroutine prterr
end module w90_libv1_types

module w90_wannier90_libv1_types
  ! global type instances (of w90_wannier90_types) to allow legacy library

  use w90_constants, only: dp, maxlen
  !use w90_io, only: maxlen
  use w90_wannier90_types

  implicit none

  public

  type(w90_calculation_type), save :: w90_calcs
  type(output_file_type), save :: out_files
  type(real_space_ham_type) :: rs_region
  type(wvfn_read_type), save :: plot
  type(band_plot_type), save :: band_plot
  type(wannier_plot_type), save :: wann_plot
  type(dis_control_type), save :: dis_data
  type(dis_spheres_type), save :: dis_spheres
  type(wann_control_type), save :: wannierise
  type(wann_omega_type), save :: wann_omega
  ! RS: symmetry-adapted Wannier functions
  logical, save :: lsitesymmetry = .false.
  real(kind=dp), save :: symmetrize_eps = 1.d-3
  !type(hamiltonian_type), save :: hamiltonian
  type(fermi_surface_plot_type), save :: fermi_surface_data
  type(transport_type), save :: tran
  type(select_projection_type), save :: select_proj

  logical, save :: eig_found

  ! a_matrix, m_matrix in disentangle and overlap
  complex(kind=dp), allocatable, save :: a_matrix(:, :, :)
  complex(kind=dp), allocatable, save :: m_matrix_orig(:, :, :, :)
  complex(kind=dp), allocatable, save :: m_matrix_orig_local(:, :, :, :)
  ! disentangle, hamiltonian, overlap and wannierise
  complex(kind=dp), allocatable, save :: m_matrix(:, :, :, :)
  ! in disentangle and overlap
  complex(kind=dp), allocatable, save :: m_matrix_local(:, :, :, :)

end module w90_wannier90_libv1_types

!================================================!
subroutine wannier_setup(seed__name, mp_grid_loc, num_kpts_loc, real_lattice_loc, &
                         recip_lattice_loc, kpt_latt_loc, num_bands_tot, num_atoms_loc, &
                         atom_symbols_loc, atoms_cart_loc, gamma_only_loc, spinors_loc, &
                         nntot_loc, nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc, &
                         proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc, proj_z_loc, &
                         proj_x_loc, proj_zona_loc, exclude_bands_loc, proj_s_loc, proj_s_qaxis_loc)
  !================================================!
  !! This routine should be called first from a code calling the library
  !! mode to setup all the variables.
  !! NOTE! The library mode currently works ONLY in serial (when called from
  !! a parallel code, make sure to run it only on 1 MPI process)
  !!
  !! For more information, check a (minimal) example of how it can be used
  !! in the folder test-suite/library-mode-test/test_library.F90
  !================================================!

  use w90_comms, only: w90comm_type
  use w90_constants, only: w90_physical_constants_type, dp
  use w90_io
  use w90_kmesh
  use w90_libv1_types
  use w90_readwrite, only: w90_readwrite_write_header, w90_readwrite_lib_set_atoms
  use w90_sitesym
  use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_wannier90_readwrite_write, &
    w90_wannier90_readwrite_w90_dealloc, w90_extra_io_type
  use w90_wannier90_libv1_types
  use w90_error, only: w90_error_type

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

  type(w90_physical_constants_type) :: physics
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
  type(timer_list_type) :: timer
  type(w90_error_type), allocatable :: error
  type(w90comm_type) :: comm

  type(w90_extra_io_type) :: write_data
  ! was in driver, only used by wannier_lib
  type(proj_input_type) :: proj
  !Projections
  logical :: lhasproj

  real(kind=dp) time0, time1
  character(len=9) :: stat, pos, cdate, ctime
  integer :: ierr
  integer :: stdout
  character(len=50)  :: seedname
  logical :: wout_found
  logical :: disentanglement
  logical :: mpiinitalready

#ifdef MPI
  call mpi_initialized(mpiinitalready, ierr)
  if (.not. mpiinitalready) then
    call set_error_mpi(error, 'mpi_init must be called before wannier_setup() when libwannier is compiled with MPI support', comm)
  endif
  if (allocated(error)) call prterr(error, stdout)
#endif

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

  call w90_readwrite_write_header(physics%bohr_version_str, physics%constants_version_str1, &
                                  physics%constants_version_str2, stdout)

  write (stdout, '(/a/)') ' Wannier90 is running in LIBRARY MODE'
  write (stdout, '(a/)') ' Setting up k-point neighbours...'

  ! copy local data into module variables
  mp_grid = mp_grid_loc
  num_kpts = num_kpts_loc
  real_lattice = real_lattice_loc
  !recip_lattice = recip_lattice_loc
  allocate (kpt_latt(3, num_kpts), stat=ierr)
  if (ierr /= 0) call set_error_alloc(error, 'Error allocating kpt_latt in wannier_setup', comm)
  if (allocated(error)) call prterr(error, stdout)
  kpt_latt = kpt_latt_loc
  atoms%num_atoms = num_atoms_loc
  call w90_readwrite_lib_set_atoms(atoms, atom_symbols_loc, atoms_cart_loc, real_lattice, &
                                   error, comm)
  if (allocated(error)) call prterr(error, stdout)
  gamma_only = gamma_only_loc
  system%spinors = spinors_loc

  ! GP: at this point we don't know yet the number of excluded bands...
  num_bands = num_bands_tot
  !library_w90_wannier90_readwrite_read_first_pass = .true.
  call w90_wannier90_readwrite_read(atoms, band_plot, dis_data, dis_spheres, dis_window, exclude_bands, fermi_energy_list, &
                                    fermi_surface_data, kmesh_data, kmesh_info, kpt_latt, out_files, &
                                    plot, wannierise, wann_omega, proj, input_proj, rs_region, select_proj, &
                                    spec_points, system, tran, verbose, wann_data, wann_plot, write_data, ws_region, &
                                    w90_calcs, eigval, real_lattice, physics%bohr, symmetrize_eps, mp_grid, &
                                    num_bands, num_kpts, num_proj, num_wann, optimisation, eig_found, calc_only_A, &
                                    cp_pp, gamma_only, lhasproj, .true., .true., lsitesymmetry, use_bloch_phases, &
                                    seedname, stdout, error, comm)
  if (allocated(error)) call prterr(error, stdout)
  have_disentangled = .false.
  disentanglement = (num_bands > num_wann)
  ! Following calls will all NOT be first_pass, and I need to pass
  ! directly num_bands, that is already set internally now to num_bands = num_bands_tot - num_exclude_bands
  !library_w90_wannier90_readwrite_read_first_pass = .false.

  call w90_wannier90_readwrite_write(atoms, band_plot, dis_data, dis_spheres, fermi_energy_list, &
                                     fermi_surface_data, kpt_latt, out_files, plot, wannierise, &
                                     proj, input_proj, rs_region, select_proj, spec_points, tran, &
                                     verbose, wann_data, wann_plot, write_data, w90_calcs, &
                                     real_lattice, symmetrize_eps, mp_grid, num_bands, num_kpts, &
                                     num_proj, num_wann, optimisation, cp_pp, gamma_only, &
                                     lsitesymmetry, system%spinors, use_bloch_phases, stdout)
  time1 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

  if (.not. kmesh_info%explicit_nnkpts) then
    call kmesh_get(kmesh_data, kmesh_info, verbose, kpt_latt, real_lattice, num_kpts, gamma_only, &
                   stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout)
  endif
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
  if (allocated(wannierise%guiding_centres%centres)) then
    proj_site_loc(:, 1:num_proj) = wannierise%guiding_centres%centres(:, 1:num_proj)
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
  if (allocated(exclude_bands)) then
    exclude_bands_loc(1:size(exclude_bands)) = exclude_bands(1:size(exclude_bands))
  end if

  if (w90_calcs%postproc_setup) then
    call kmesh_write(exclude_bands, kmesh_info, input_proj, verbose, kpt_latt, real_lattice, &
                     num_kpts, num_proj, calc_only_A, system%spinors, seedname, timer)
    write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    write (stdout, '(/a)') ' '//trim(seedname)//'.nnkp written.'
  endif

  call kmesh_dealloc(kmesh_info, error, comm)
  if (allocated(error)) call prterr(error, stdout)

  call w90_wannier90_readwrite_w90_dealloc(atoms, band_plot, dis_spheres, dis_window, &
                                           exclude_bands, kmesh_data, kpt_latt, wannierise, proj, &
                                           input_proj, select_proj, spec_points, wann_data, &
                                           wann_plot, write_data, eigval, error, comm)
  if (allocated(error)) call prterr(error, stdout)
  write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'

  write (stdout, '(/a/)') ' Finished setting up k-point neighbours.'

  call io_date(cdate, ctime)

  write (stdout, '(2a)') ' Exiting wannier_setup in wannier90 ', ctime

  close (stdout)

end subroutine wannier_setup

!================================================!
subroutine wannier_run(seed__name, mp_grid_loc, num_kpts_loc, real_lattice_loc, recip_lattice_loc, &
                       kpt_latt_loc, num_bands_loc, num_wann_loc, nntot_loc, num_atoms_loc, &
                       atom_symbols_loc, atoms_cart_loc, gamma_only_loc, M_matrix_loc, A_matrix_loc, &
                       eigenvalues_loc, U_matrix_loc, U_matrix_opt_loc, lwindow_loc, &
                       wann_centres_loc, wann_spreads_loc, spread_loc)

  !================================================!
  !! This routine should be called after wannier_setup from a code calling
  !! the library mode to actually run the Wannier code.
  !!
  !! NOTE! The library mode currently works ONLY in serial.
  !! When called from an external code, wannier90 needs to be compiled
  !! in sequential and wannier_run called with 1 MPI process.
  !!
  !! For more information, check a (minimal) example of how it can be used
  !! in the folder test-suite/library-mode-test/test_library.F90
  !================================================!

  use w90_constants, only: w90_physical_constants_type, dp
  use w90_libv1_types
  use w90_wannier90_libv1_types
  use w90_wannier90_readwrite, only: w90_wannier90_readwrite_read, w90_wannier90_readwrite_write, &
    w90_wannier90_readwrite_write_chkpt, w90_wannier90_readwrite_w90_dealloc, w90_extra_io_type
  use w90_io
  use w90_hamiltonian
  use w90_kmesh
  use w90_disentangle
  use w90_overlap
  use w90_wannierise
  use w90_plot
  use w90_transport
  use w90_comms, only: comms_array_split, comms_scatterv, w90comm_type, &
    mpisize, mpirank
  use w90_readwrite, only: w90_readwrite_lib_set_atoms
  use w90_error

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

  type(w90_physical_constants_type) :: physics
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
  integer :: stdout
  character(len=50) :: seedname
  logical :: wout_found

  complex(kind=dp), allocatable :: ham_r(:, :, :)
  integer, allocatable :: irvec(:, :)
  integer, allocatable :: shift_vec(:, :)
  integer, allocatable :: ndegen(:)
  integer :: rpt_origin
  real(kind=dp), allocatable :: wannier_centres_translated(:, :)
  complex(kind=dp), allocatable :: ham_k(:, :, :)
  integer :: nrpts

  type(sitesym_type) :: sym
  type(ham_logical_type) :: hmlg

  type(w90_extra_io_type) :: write_data
  type(proj_input_type) :: proj
  !Projections
  logical :: lhasproj

! Needed to split an array on different nodes (not done here)
  integer, allocatable :: counts(:)
  integer, allocatable :: displs(:)
  integer :: num_nodes, my_node_id
  type(w90comm_type) :: comm
  logical :: disentanglement
  logical :: mpiinitalready

  type(timer_list_type) :: timer
  type(w90_error_type), allocatable :: error

  ! CORRECT ONLY FOR SERIAL CASE!!!
  ! THESE LIBRARY ROUTINES ARE OBSOLETE
  !
  ! depending on how the rest of the code is compiled (w or w/o MPI)
  ! the various functions may or may not require a legitimate
  ! MPI communicator
  ! because the old "version 1" library interface does not support
  ! providing a communicator, we open mpi_comm_world here and pass it
  !
  ! it is expected that this library is only invoked with one process
  ! even when compiled with MPI
  ! use with more than one process is not supported/tested
  ! --JJ 13Aug21
  !
#ifdef MPI
  call mpi_initialized(mpiinitalready, ierr)
  if (.not. mpiinitalready) then
    call set_error_mpi(error, 'mpi_init must be called before wannier_run() when libwannier is compiled with MPI support', comm)
  endif
  if (allocated(error)) call prterr(error, stdout)
  comm%comm = MPI_COMM_WORLD
  num_nodes = 1
  my_node_id = 0
#else
  num_nodes = 1
  my_node_id = 0
#endif

  allocate (counts(0:num_nodes - 1)); 
  allocate (displs(0:num_nodes - 1)); 
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

!  call w90_readwrite_write_header

  num_bands = num_bands_loc
  mp_grid = mp_grid_loc
  num_kpts = num_kpts_loc
  real_lattice = real_lattice_loc
  allocate (kpt_latt(3, num_kpts), stat=ierr)
  if (ierr /= 0) call set_error_alloc(error, 'Error allocating kpt_latt in wannier_setup', comm)
  if (allocated(error)) call prterr(error, stdout)

  kpt_latt = kpt_latt_loc
  allocate (eigval(num_bands, num_kpts), stat=ierr)
  if (ierr /= 0) call set_error_alloc(error, 'Error allocating eigval in wannier_setup', comm)
  if (allocated(error)) call prterr(error, stdout)
  eigval = eigenvalues_loc
  atoms%num_atoms = num_atoms_loc
  gamma_only = gamma_only_loc

  call w90_readwrite_lib_set_atoms(atoms, atom_symbols_loc, atoms_cart_loc, real_lattice, &
                                   error, comm)
  if (allocated(error)) call prterr(error, stdout)

  call w90_wannier90_readwrite_read(atoms, band_plot, dis_data, dis_spheres, dis_window, exclude_bands, &
                                    fermi_energy_list, fermi_surface_data, kmesh_data, kmesh_info, kpt_latt, &
                                    out_files, plot, wannierise, wann_omega, proj, input_proj, rs_region, &
                                    select_proj, spec_points, system, tran, verbose, wann_data, wann_plot, &
                                    write_data, ws_region, w90_calcs, eigval, real_lattice, physics%bohr, &
                                    symmetrize_eps, mp_grid, num_bands, num_kpts, num_proj, num_wann, optimisation, &
                                    eig_found, calc_only_A, cp_pp, gamma_only, lhasproj, .true., .false., &
                                    lsitesymmetry, use_bloch_phases, seedname, stdout, error, comm)
  if (allocated(error)) call prterr(error, stdout)
  have_disentangled = .false.
  disentanglement = (num_bands > num_wann)
  call w90_wannier90_readwrite_write(atoms, band_plot, dis_data, dis_spheres, fermi_energy_list, &
                                     fermi_surface_data, kpt_latt, out_files, plot, wannierise, &
                                     proj, input_proj, rs_region, select_proj, spec_points, tran, &
                                     verbose, wann_data, wann_plot, write_data, w90_calcs, &
                                     real_lattice, symmetrize_eps, mp_grid, num_bands, num_kpts, &
                                     num_proj, num_wann, optimisation, cp_pp, gamma_only, &
                                     lsitesymmetry, system%spinors, use_bloch_phases, stdout)
  time1 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

  call kmesh_get(kmesh_data, kmesh_info, verbose, kpt_latt, real_lattice, num_kpts, gamma_only, &
                 stdout, timer, error, comm)
  if (allocated(error)) call prterr(error, stdout)

  time2 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time to get kmesh        ', time2 - time1, ' (sec)'

  call comms_array_split(num_kpts, counts, displs, comm)
  call overlap_allocate(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, &
                        u_matrix, u_matrix_opt, kmesh_info%nntot, num_bands, num_kpts, num_wann, &
                        verbose%timing_level, timer, error, comm)
  if (allocated(error)) call prterr(error, stdout)
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
    call comms_scatterv(m_matrix_orig_local, &
                        num_bands*num_bands*kmesh_info%nntot*counts(my_node_id), m_matrix_orig, &
                        num_bands*num_bands*kmesh_info%nntot*counts, &
                        num_bands*num_bands*kmesh_info%nntot*displs, error, comm)
    if (allocated(error)) call prterr(error, stdout)
  else
    call comms_scatterv(m_matrix_local, num_wann*num_wann*kmesh_info%nntot*counts(my_node_id), &
                        m_matrix, num_wann*num_wann*kmesh_info%nntot*counts, num_wann*num_wann* &
                        kmesh_info%nntot*displs, error, comm)
    if (allocated(error)) call prterr(error, stdout)
  endif

!~  ! Check Mmn(k,b) is symmetric in m and n for gamma_only case
!~  if (gamma_only) call overlap_check_m_symmetry()

  if (disentanglement) then
    have_disentangled = .false.

    call dis_main(dis_data, dis_spheres, dis_window, kmesh_info, kpt_latt, sym, verbose, a_matrix, &
                  m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, u_matrix, &
                  u_matrix_opt, eigval, real_lattice, wann_omega%invariant, &
                  num_bands, num_kpts, num_wann, optimisation, gamma_only, lsitesymmetry, &
                  stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout)
    have_disentangled = .true.
    call w90_wannier90_readwrite_write_chkpt('postdis', exclude_bands, wann_data, kmesh_info, &
                                             kpt_latt, num_kpts, dis_window, num_bands, num_wann, u_matrix, &
                                             u_matrix_opt, m_matrix, mp_grid, real_lattice, &
                                             wann_omega%invariant, have_disentangled, stdout, seedname)

    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle      ', time1 - time2, ' (sec)'
  else
    if (gamma_only) then
      call overlap_project_gamma(m_matrix, u_matrix, kmesh_info%nntot, num_wann, &
                                 verbose%timing_level, stdout, timer, error, comm)
      if (allocated(error)) call prterr(error, stdout)
    else
      call overlap_project(sym, m_matrix, m_matrix_local, u_matrix, kmesh_info%nnlist, &
                           kmesh_info%nntot, num_bands, num_kpts, num_wann, &
                           verbose%timing_level, lsitesymmetry, stdout, timer, error, comm)
      if (allocated(error)) call prterr(error, stdout)
    endif
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to project overlaps ', time1 - time2, ' (sec)'
  end if

  if (gamma_only) then
    call wann_main_gamma(atoms, dis_window, exclude_bands, kmesh_info, kpt_latt, out_files, &
                         wannierise, wann_omega, system, verbose, wann_data, m_matrix, &
                         u_matrix, u_matrix_opt, eigval, real_lattice, mp_grid, &
                         num_bands, num_kpts, num_wann, have_disentangled, &
                         rs_region%translate_home_cell, seedname, stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout)
  else
    call wann_main(atoms, dis_window, exclude_bands, hmlg, kmesh_info, kpt_latt, out_files, &
                   rs_region, wannierise, wann_omega, sym, system, verbose, wann_data, &
                   ws_region, w90_calcs, ham_k, ham_r, m_matrix, u_matrix, u_matrix_opt, eigval, &
                   real_lattice, wannier_centres_translated, irvec, mp_grid, ndegen, shift_vec, &
                   nrpts, num_bands, num_kpts, num_proj, num_wann, optimisation, rpt_origin, &
                   band_plot%mode, tran%mode, have_disentangled, lsitesymmetry, &
                   seedname, stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout)
  endif

  call w90_wannier90_readwrite_write_chkpt('postwann', exclude_bands, wann_data, kmesh_info, kpt_latt, &
                                           num_kpts, dis_window, num_bands, num_wann, u_matrix, u_matrix_opt, &
                                           m_matrix, mp_grid, real_lattice, &
                                           wann_omega%invariant, have_disentangled, stdout, seedname)

  time2 = io_time()
  write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time2 - time1, ' (sec)'

  if (w90_calcs%wannier_plot .or. w90_calcs%bands_plot .or. w90_calcs%fermi_surface_plot .or. out_files%write_hr) then
    call plot_main(atoms, band_plot, dis_window, fermi_energy_list, fermi_surface_data, hmlg, &
                   kmesh_info, kpt_latt, out_files, plot, rs_region, spec_points, &
                   verbose, wann_data, wann_plot, ws_region, w90_calcs, ham_k, ham_r, m_matrix, &
                   u_matrix, u_matrix_opt, eigval, real_lattice, &
                   wannier_centres_translated, physics%bohr, irvec, mp_grid, ndegen, shift_vec, &
                   nrpts, num_bands, num_kpts, num_wann, rpt_origin, tran%mode, have_disentangled, &
                   lsitesymmetry, system%spinors, seedname, stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout)
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  end if

  time2 = io_time()
  if (w90_calcs%transport) then
    call tran_main(atoms, dis_window, fermi_energy_list, hmlg, kpt_latt, out_files, rs_region, &
                   tran, verbose, wann_data, ws_region, w90_calcs, ham_k, ham_r, u_matrix, &
                   u_matrix_opt, eigval, real_lattice, wannier_centres_translated, &
                   irvec, mp_grid, ndegen, shift_vec, nrpts, num_bands, num_kpts, num_wann, &
                   rpt_origin, band_plot%mode, have_disentangled, lsitesymmetry, seedname, &
                   stdout, timer, error, comm)
    if (allocated(error)) call prterr(error, stdout)
    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
  end if

  ! Now we zero all of the local output data, then copy in the data
  ! from the parameters module

  u_matrix_loc = u_matrix
  if (present(u_matrix_opt_loc) .and. present(lwindow_loc)) then
  if (disentanglement) then
    u_matrix_opt_loc = u_matrix_opt
    lwindow_loc = dis_window%lwindow
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
    spread_loc(1) = wann_omega%total
    spread_loc(2) = wann_omega%invariant
    spread_loc(3) = wann_omega%tilde
  endif
  call hamiltonian_dealloc(hmlg, ham_k, ham_r, wannier_centres_translated, irvec, ndegen, &
                           error, comm)
  if (allocated(error)) call prterr(error, stdout)

  call overlap_dealloc(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, m_matrix_orig_local, &
                       u_matrix, u_matrix_opt, error, comm)
  if (allocated(error)) call prterr(error, stdout)
  call kmesh_dealloc(kmesh_info, error, comm)
  if (allocated(error)) call prterr(error, stdout)
  call w90_wannier90_readwrite_w90_dealloc(atoms, band_plot, dis_spheres, dis_window, &
                                           exclude_bands, kmesh_data, kpt_latt, wannierise, proj, &
                                           input_proj, select_proj, spec_points, wann_data, &
                                           wann_plot, write_data, eigval, error, comm)
  if (allocated(error)) call prterr(error, stdout)
  write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time() - time0, ' (sec)'

  if (verbose%timing_level > 0) call io_print_timings(timer, stdout)

  write (stdout, *)
  write (stdout, '(1x,a)') 'All done: wannier90 exiting'
  close (stdout)

end subroutine wannier_run
