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

module w90_wannierise
  !! Main routines for the minimisation of the spread

  use w90_constants, only: dp

  implicit none

  private

  public :: wann_main
  public :: wann_main_gamma  ![ysl]

  type localisation_vars
    !! Contributions to the spread
    real(kind=dp) :: om_i
    !! Gauge Invarient
    real(kind=dp) :: om_d
    !! Diagonal
    real(kind=dp) :: om_od
    !! Off-diagonal
    real(kind=dp) :: om_tot
    !! Total
    real(kind=dp) :: om_iod
    !! Combined I-OD term for selective localization
    real(kind=dp) :: om_nu
    !! Lagrange multiplier term due to constrained centres
    !! real(kind=dp) :: om_c
    !! Total spread functional with constraints
!~     real(kind=dp) :: om_1
!~     real(kind=dp) :: om_2
!~     real(kind=dp) :: om_3
  end type localisation_vars

contains

  !==================================================================!
  subroutine wann_main(atoms, dis_window, excluded_bands, hmlg, kmesh_info, k_points, out_files, &
                       param_hamil, wannierise, sym, system, verbose, wann_data, &
                       ws_region, w90_calcs, ham_k, ham_r, m_matrix, u_matrix, u_matrix_opt, &
                       eigval, real_lattice, recip_lattice, wannier_centres_translated, irvec, &
                       mp_grid, ndegen, shift_vec, nrpts, num_bands, num_kpts, num_proj, num_wann, &
                       rpt_origin, bands_plot_mode, transport_mode, have_disentangled, &
                       lsitesymmetry, seedname, stdout, comm)
    !==================================================================!
    !                                                                  !
    !! Calculate the Unitary Rotations to give Maximally Localised Wannier Functions
    !                                                                  !
    !===================================================================
    use w90_constants, only: dp, cmplx_1, cmplx_0, twopi, cmplx_i
    use w90_io, only: io_error, io_wallclocktime, io_stopwatch, io_file_unit
    use wannier_param_types, only: wannierise_type, output_file_type, &
      w90_calculation_type, param_hamiltonian_type
    use w90_param_types, only: kmesh_info_type, print_output_type, &
      wannier_data_type, atom_data_type, k_point_type, disentangle_manifold_type, w90_system_type, &
      exclude_bands_type, ws_region_type
    use wannier_methods, only: param_write_chkpt
    use w90_utility, only: utility_frac_to_cart, utility_zgemm
    use w90_sitesym, only: sitesym_symmetrize_gradient, sitesym_data
    use w90_comms, only: mpisize, mpirank, comms_gatherv, comms_bcast, &
      comms_scatterv, comms_array_split, w90commtype

    !ivo
    use w90_hamiltonian, only: hamiltonian_setup, hamiltonian_get_hr, ham_logical

    implicit none

    ! passed variables
    type(atom_data_type), intent(in) :: atoms
    type(disentangle_manifold_type), intent(in) :: dis_window
    type(ham_logical), intent(inout) :: hmlg
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(k_point_type), intent(in) :: k_points
    type(w90_system_type), intent(in) :: system
    type(ws_region_type), intent(in) :: ws_region
    type(exclude_bands_type), intent(in) :: excluded_bands
    type(print_output_type), intent(in) :: verbose
    type(output_file_type), intent(in) :: out_files
    type(param_hamiltonian_type), intent(inout) :: param_hamil
    type(wannierise_type), intent(inout) :: wannierise
    type(sitesym_data), intent(in) :: sym
    type(w90_calculation_type), intent(in) :: w90_calcs
    type(w90commtype), intent(in) :: comm
    type(wannier_data_type), intent(inout) :: wann_data

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_proj
    integer, intent(in) :: num_wann
    integer, intent(inout), allocatable :: irvec(:, :)
    integer, intent(inout), allocatable :: ndegen(:)
    integer, intent(inout), allocatable :: shift_vec(:, :)
    integer, intent(inout) :: nrpts
    integer, intent(inout) :: rpt_origin
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(inout), allocatable :: wannier_centres_translated(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)

    complex(kind=dp), intent(inout), allocatable :: ham_k(:, :, :)
    complex(kind=dp), intent(inout), allocatable :: ham_r(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)

    logical, intent(in) :: lsitesymmetry
    logical, intent(in) :: have_disentangled

    character(len=*), intent(in) :: bands_plot_mode
    character(len=*), intent(in) :: transport_mode
    character(len=50), intent(in) :: seedname

    ! local variables
    type(localisation_vars) :: old_spread
    type(localisation_vars) :: wann_spread
    type(localisation_vars) :: trial_spread

    ! Data to avoid large allocation within iteration loop
    real(kind=dp), allocatable  :: rnkb(:, :, :)
    real(kind=dp), allocatable  :: rnkb_loc(:, :, :)
    real(kind=dp), allocatable  :: ln_tmp(:, :, :)
    real(kind=dp), allocatable  :: ln_tmp_loc(:, :, :)

    ! for MPI
    complex(kind=dp), allocatable  :: u_matrix_loc(:, :, :)
    complex(kind=dp), allocatable  :: m_matrix_loc(:, :, :, :)
    complex(kind=dp), allocatable  :: cdq_loc(:, :, :) ! the only large array sent from process to process in the main loop
    complex(kind=dp), allocatable  :: cdodq_loc(:, :, :)
    integer, allocatable  :: counts(:)
    integer, allocatable  :: displs(:)

    logical :: first_pass
    !! Used to trigger the calculation of the invarient spread we only need to do this on entering wann_main (_gamma)
    real(kind=dp) :: lambda_loc
    ! end of wannierise module data

    ! guiding centres
    real(kind=dp), allocatable :: rguide(:, :)
    integer :: irguide

    ! local arrays used and passed in subroutines
    complex(kind=dp), allocatable :: csheet(:, :, :)
    complex(kind=dp), allocatable :: cdodq(:, :, :)
    complex(kind=dp), allocatable :: cdodq_r(:, :, :)
    complex(kind=dp), allocatable :: k_to_r(:, :)
    complex(kind=dp), allocatable :: cdodq_precond(:, :, :)
    complex(kind=dp), allocatable :: cdodq_precond_loc(:, :, :)
    real(kind=dp), allocatable :: sheet(:, :, :)
    real(kind=dp), allocatable :: rave(:, :), r2ave(:), rave2(:)
    !real(kind=dp), dimension(3) :: rvec_cart

    !local arrays not passed into subroutines
    complex(kind=dp), allocatable  :: cwschur1(:), cwschur2(:)
    complex(kind=dp), allocatable  :: cwschur3(:), cwschur4(:)
    complex(kind=dp), allocatable  :: cdq(:, :, :)!,cdqkeep(:,:,:)
    ! cdqkeep is replaced by cdqkeep_loc
    complex(kind=dp), allocatable  :: cdqkeep_loc(:, :, :)
    complex(kind=dp), allocatable  :: cz(:, :)
    complex(kind=dp), allocatable  :: cmtmp(:, :), tmp_cdq(:, :)
    ! complex(kind=dp), allocatable  :: m0(:,:,:,:),u0(:,:,:)
    ! m0 and u0 are replaced by m0_loc and u0_loc
    complex(kind=dp), allocatable  :: m0_loc(:, :, :, :), u0_loc(:, :, :)
    complex(kind=dp), allocatable  :: cwork(:)
    real(kind=dp), allocatable  :: evals(:)
    real(kind=dp), allocatable  :: rwork(:)

    real(kind=dp) :: doda0
    real(kind=dp) :: falphamin, alphamin
    real(kind=dp) :: gcfac, gcnorm1, gcnorm0
    integer       :: i, n, iter, ind, ierr, iw, ncg, nkp, nkp_loc !, nn
    logical       :: lprint, ldump, lquad
    real(kind=dp), allocatable :: history(:)
    real(kind=dp)              :: save_spread
    logical                    :: lconverged, lrandom, lfirst
    integer                    :: conv_count, noise_count, page_unit
    complex(kind=dp) :: rdotk !, fac
    !real(kind=dp) :: alpha_precond
    integer :: irpt, loop_kpt
    !logical :: cconverged
    !real(kind=dp) :: glpar, cvalue_new
    real(kind=dp), allocatable :: rnr0n2(:)

    ! pllel setup
    logical :: on_root = .false.
    integer :: num_nodes, my_node_id

    num_nodes = mpisize(comm)
    my_node_id = mpirank(comm)
    if (my_node_id == 0) on_root = .true.

    if (verbose%timing_level > 0 .and. verbose%iprint > 0) call io_stopwatch('wann: main', 1, stdout, seedname)

    first_pass = .true.

    ! Allocate stuff

    allocate (history(wannierise%control%conv_window), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating history in wann_main', stdout, seedname)

    ! module data
!    if(optimisation>0) then
!       allocate(  m0 (num_wann, num_wann, nntot, num_kpts),stat=ierr)
!    end if
!    if (ierr/=0) call io_error('Error in allocating m0 in wann_main')
!    allocate(  u0 (num_wann, num_wann, num_kpts),stat=ierr)
!    if (ierr/=0) call io_error('Error in allocating u0 in wann_main')
    allocate (rnkb(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rnkb in wann_main', stdout, seedname)
    allocate (ln_tmp(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ln_tmp in wann_main', stdout, seedname)
    if (wannierise%constrain%selective_loc) then
      allocate (rnr0n2(wannierise%constrain%slwf_num), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating rnr0n2 in wann_main', stdout, seedname)
    end if

    rnkb = 0.0_dp

    ! sub vars passed into other subs
    allocate (csheet(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating csheet in wann_main', stdout, seedname)
    allocate (cdodq(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cdodq in wann_main', stdout, seedname)
    allocate (sheet(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating sheet in wann_main', stdout, seedname)
    allocate (rave(3, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rave in wann_main', stdout, seedname)
    allocate (r2ave(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating r2ave in wann_main', stdout, seedname)
    allocate (rave2(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rave2 in wann_main', stdout, seedname)
    allocate (rguide(3, num_wann))
    if (ierr /= 0) call io_error('Error in allocating rguide in wann_main', stdout, seedname)

    if (wannierise%control%precond) then
      call hamiltonian_setup(hmlg, verbose, ws_region, w90_calcs, ham_k, ham_r, &
                             real_lattice, wannier_centres_translated, irvec, mp_grid, ndegen, &
                             num_kpts, num_wann, nrpts, rpt_origin, bands_plot_mode, stdout, &
                             seedname, transport_mode)
      allocate (cdodq_r(num_wann, num_wann, nrpts), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating cdodq_r in wann_main', stdout, seedname)
      allocate (cdodq_precond(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating cdodq_precond in wann_main', stdout, seedname)

      ! this method of computing the preconditioning is much more efficient, but requires more RAM
      if (verbose%optimisation >= 3) then
        allocate (k_to_r(num_kpts, nrpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating k_to_r in wann_main', stdout, seedname)

        do irpt = 1, nrpts
          do loop_kpt = 1, num_kpts
            rdotk = twopi*dot_product(k_points%kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
            k_to_r(loop_kpt, irpt) = exp(-cmplx_i*rdotk)
          enddo
        enddo
      end if
    end if

    csheet = cmplx_1; cdodq = cmplx_0
    sheet = 0.0_dp; rave = 0.0_dp; r2ave = 0.0_dp; rave2 = 0.0_dp; rguide = 0.0_dp

    ! sub vars not passed into other subs
    allocate (cwschur1(num_wann), cwschur2(10*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwshur1 in wann_main', stdout, seedname)
    allocate (cwschur3(num_wann), cwschur4(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwshur3 in wann_main', stdout, seedname)
    allocate (cdq(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cdq in wann_main', stdout, seedname)

    ! for MPI
    if (allocated(counts)) deallocate (counts)
    allocate (counts(0:num_nodes - 1), stat=ierr)
    if (ierr /= 0) then
      call io_error('Error in allocating counts in wann_main', stdout, seedname)
    end if

    if (allocated(displs)) deallocate (displs)
    allocate (displs(0:num_nodes - 1), stat=ierr)
    if (ierr /= 0) then
      call io_error('Error in allocating displs in wann_main', stdout, seedname)
    end if
    call comms_array_split(num_kpts, counts, displs, comm)
    allocate (rnkb_loc(num_wann, kmesh_info%nntot, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rnkb_loc in wann_main', stdout, seedname)
    allocate (ln_tmp_loc(num_wann, kmesh_info%nntot, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ln_tmp_loc in wann_main', stdout, seedname)
    allocate (u_matrix_loc(num_wann, num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating u_matrix_loc in wann_main', stdout, seedname)
    allocate (m_matrix_loc(num_wann, num_wann, kmesh_info%nntot, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating m_matrix_loc in wann_main', stdout, seedname)
!    allocate( m_matrix_1b  (num_wann, num_wann, num_kpts),stat=ierr )
!    if (ierr/=0) call io_error('Error in allocating m_matrix_1b in wann_main')
!    allocate( m_matrix_1b_loc  (num_wann, num_wann, max(1,counts(my_node_id))),stat=ierr )
!    if (ierr/=0) call io_error('Error in allocating m_matrix_1b_loc in wann_main')
    if (wannierise%control%precond) then
      allocate (cdodq_precond_loc(num_wann, num_wann, max(1, counts(my_node_id))), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating cdodq_precond_loc in wann_main', stdout, seedname)
    end if
    ! initialize local u and m matrices with global ones
    do nkp_loc = 1, counts(my_node_id)
      nkp = nkp_loc + displs(my_node_id)
!       m_matrix_loc (:,:,:, nkp_loc) = &
!           m_matrix (:,:,:, nkp)
      u_matrix_loc(:, :, nkp_loc) = &
        u_matrix(:, :, nkp)
    end do
    call comms_scatterv(m_matrix_loc, num_wann*num_wann*kmesh_info%nntot*counts(my_node_id), &
                        m_matrix, num_wann*num_wann*kmesh_info%nntot*counts, &
                        num_wann*num_wann*kmesh_info%nntot*displs, stdout, seedname, comm)

    allocate (cdq_loc(num_wann, num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cdq_loc in wann_main', stdout, seedname)
    allocate (cdodq_loc(num_wann, num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cdodq_loc in wann_main', stdout, seedname)
    allocate (cdqkeep_loc(num_wann, num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cdqkeep_loc in wann_main', stdout, seedname)
    if (verbose%optimisation > 0) then
      allocate (m0_loc(num_wann, num_wann, kmesh_info%nntot, max(1, counts(my_node_id))), stat=ierr)
    end if
    if (ierr /= 0) call io_error('Error in allocating m0_loc in wann_main', stdout, seedname)
    allocate (u0_loc(num_wann, num_wann, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating u0_loc in wann_main', stdout, seedname)

    allocate (cz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cz in wann_main', stdout, seedname)
    allocate (cmtmp(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cmtmp in wann_main', stdout, seedname)
    allocate (tmp_cdq(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating tmp_cdq in wann_main', stdout, seedname)
    allocate (evals(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating evals in wann_main', stdout, seedname)
    allocate (cwork(4*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwork in wann_main', stdout, seedname)
    allocate (rwork(3*num_wann - 2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rwork in wann_main', stdout, seedname)

    cwschur1 = cmplx_0; cwschur2 = cmplx_0; cwschur3 = cmplx_0; cwschur4 = cmplx_0
    cdq = cmplx_0; cz = cmplx_0; cmtmp = cmplx_0; cdqkeep_loc = cmplx_0; cdq_loc = cmplx_0; ! buff=cmplx_0;

    gcnorm1 = 0.0_dp; gcnorm0 = 0.0_dp

    ! initialise rguide to projection centres (Cartesians in units of Ang)
    if (wannierise%control%guiding_centres) then
      do n = 1, num_proj
        call utility_frac_to_cart(wannierise%proj_site(:, n), rguide(:, n), real_lattice)
      enddo
!       if(spinors) then ! not needed with new changes to spinor proj 2013 JRY
!          do n=1,num_proj
!             call utility_frac_to_cart(proj_site(:,n),rguide(:,n+num_proj),real_lattice)
!          enddo
!       end if
    end if

    if (verbose%iprint > 0) then
      write (stdout, *)
      write (stdout, '(1x,a)') '*------------------------------- WANNIERISE ---------------------------------*'
      write (stdout, '(1x,a)') '+--------------------------------------------------------------------+<-- CONV'
      if (verbose%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '| Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV'
      else
        write (stdout, '(1x,a)') '| Iter  Delta Spread     RMS Gradient      Spread (Bohr^2)     Time  |<-- CONV'
      endif
      write (stdout, '(1x,a)') '+--------------------------------------------------------------------+<-- CONV'
      write (stdout, *)
    endif

    irguide = 0
    if (wannierise%control%guiding_centres .and. (wannierise%control%num_no_guide_iter .le. 0)) then
      call wann_phases(csheet, sheet, rguide, irguide, num_wann, kmesh_info, num_kpts, m_matrix, &
                       .false., counts, displs, m_matrix_loc, rnkb, verbose%timing_level, &
                       stdout, seedname, verbose%iprint, comm)
      irguide = 1
    endif

    ! constrained centres part
    lambda_loc = 0.0_dp
    if (wannierise%constrain%selective_loc .and. wannierise%constrain%slwf_constrain) then
      lambda_loc = wannierise%constrain%slwf_lambda
    end if

    ! calculate initial centers and spread
    call wann_omega(csheet, sheet, rave, r2ave, rave2, wann_spread, num_wann, kmesh_info, &
                    num_kpts, verbose, wannierise%constrain, &
                    wannierise%omega%invariant, counts, displs, ln_tmp_loc, m_matrix_loc, &
                    lambda_loc, first_pass, stdout, seedname, comm)

    ! public variables
    if (.not. wannierise%constrain%selective_loc) then
      wannierise%omega%total = wann_spread%om_tot
      wannierise%omega%invariant = wann_spread%om_i
      wannierise%omega%tilde = wann_spread%om_d + wann_spread%om_od
    else
      wannierise%omega%total = wann_spread%om_tot
      ! omega_invariant = wann_spread%om_iod
      ! omega_tilde = wann_spread%om_d + wann_spread%om_nu
    end if

    ! public arrays of Wannier centres and spreads
    wann_data%centres = rave
    wann_data%spreads = r2ave - rave2

    if (wannierise%control%lfixstep) lquad = .false.
    ncg = 0
    iter = 0
    old_spread%om_tot = 0.0_dp

    ! print initial state
    if (verbose%iprint > 0) then
      write (stdout, '(1x,a78)') repeat('-', 78)
      write (stdout, '(1x,a)') 'Initial State'
      do iw = 1, num_wann
        write (stdout, 1000) iw, (rave(ind, iw)*verbose%lenconfac, ind=1, 3), &
          (r2ave(iw) - rave2(iw))*verbose%lenconfac**2
      end do
      write (stdout, 1001) (sum(rave(ind, :))*verbose%lenconfac, ind=1, 3), &
        (sum(r2ave) - sum(rave2))*verbose%lenconfac**2
      write (stdout, *)
      if (wannierise%constrain%selective_loc .and. wannierise%constrain%slwf_constrain) then
        write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') iter, &
          (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
          sqrt(abs(gcnorm1))*verbose%lenconfac, &
          wann_spread%om_tot*verbose%lenconfac**2, io_wallclocktime(), '<-- CONV'
        write (stdout, '(7x,a,F15.7,a,F15.7,a,F15.7,a,F15.7,a)') &
          'O_D=', wann_spread%om_d*verbose%lenconfac**2, &
          ' O_IOD=', (wann_spread%om_iod + wann_spread%om_nu)*verbose%lenconfac**2, &
          ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
        write (stdout, '(1x,a78)') repeat('-', 78)
      elseif (wannierise%constrain%selective_loc .and. .not. wannierise%constrain%slwf_constrain) then
        write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') iter, &
          (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
          sqrt(abs(gcnorm1))*verbose%lenconfac, &
          wann_spread%om_tot*verbose%lenconfac**2, io_wallclocktime(), '<-- CONV'
        write (stdout, '(7x,a,F15.7,a,F15.7,a,F15.7,a)') &
          'O_D=', wann_spread%om_d*verbose%lenconfac**2, &
          ' O_IOD=', wann_spread%om_iod*verbose%lenconfac**2, &
          ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
        write (stdout, '(1x,a78)') repeat('-', 78)
      else
        write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') iter, &
          (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
          sqrt(abs(gcnorm1))*verbose%lenconfac, &
          wann_spread%om_tot*verbose%lenconfac**2, io_wallclocktime(), '<-- CONV'
        write (stdout, '(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
          'O_D=', wann_spread%om_d*verbose%lenconfac**2, ' O_OD=', &
          wann_spread%om_od*verbose%lenconfac**2, &
          ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
        write (stdout, '(1x,a78)') repeat('-', 78)
      end if
    endif

    lconverged = .false.; lfirst = .true.; lrandom = .false.
    conv_count = 0; noise_count = 0

    if (.not. wannierise%control%lfixstep .and. verbose%optimisation <= 0) then
      page_unit = io_file_unit()
      open (unit=page_unit, status='scratch', form='unformatted')
    endif

    ! main iteration loop
    do iter = 1, wannierise%control%num_iter

      lprint = .false.
      if ((mod(iter, wannierise%control%num_print_cycles) .eq. 0) .or. (iter .eq. 1) &
          .or. (iter .eq. wannierise%control%num_iter)) lprint = .true.

      ldump = .false.
      if ((wannierise%control%num_dump_cycles .gt. 0) .and. &
          (mod(iter, wannierise%control%num_dump_cycles) .eq. 0)) ldump = .true.

      if (lprint .and. verbose%iprint > 0) write (stdout, '(1x,a,i6)') 'Cycle: ', iter

      if (wannierise%control%guiding_centres .and. (iter .gt. wannierise%control%num_no_guide_iter) &
          .and. (mod(iter, wannierise%control%num_guide_cycles) .eq. 0)) then
        call wann_phases(csheet, sheet, rguide, irguide, num_wann, kmesh_info, num_kpts, m_matrix, &
                         .false., counts, displs, m_matrix_loc, rnkb, verbose%timing_level, &
                         stdout, seedname, verbose%iprint, comm)
        irguide = 1
      endif

      ! calculate gradient of omega

      if (lsitesymmetry .or. wannierise%control%precond) then
        call wann_domega(csheet, sheet, rave, num_wann, kmesh_info, num_kpts, &
                         wannierise%constrain, lsitesymmetry, counts, displs, ln_tmp_loc, &
                         m_matrix_loc, rnkb_loc, cdodq_loc, lambda_loc, verbose%timing_level, &
                         stdout, seedname, sym, comm, verbose%iprint, cdodq)
      else
        call wann_domega(csheet, sheet, rave, num_wann, kmesh_info, num_kpts, &
                         wannierise%constrain, lsitesymmetry, counts, displs, ln_tmp_loc, &
                         m_matrix_loc, rnkb_loc, cdodq_loc, lambda_loc, verbose%timing_level, &
                         stdout, seedname, sym, comm, verbose%iprint)
      endif

      if (lprint .and. verbose%iprint > 2) &
        write (stdout, *) ' LINE --> Iteration                     :', iter

      ! calculate search direction (cdq)
      if (wannierise%control%precond) then
        call precond_search_direction(cdodq, cdodq_r, cdodq_precond, cdodq_precond_loc, &
                                      k_to_r, wann_spread, num_wann, num_kpts, &
                                      k_points%kpt_latt, real_lattice, nrpts, irvec, ndegen, &
                                      counts, displs, stdout)
      endif
      call internal_search_direction(cdodq_precond_loc, cdqkeep_loc, iter, lprint, lrandom, &
                                     noise_count, ncg, gcfac, gcnorm0, gcnorm1, doda0, &
                                     wannierise%control, num_wann, &
                                     kmesh_info%wbtot, cdq_loc, cdodq_loc, counts, stdout)
      if (lsitesymmetry) call sitesym_symmetrize_gradient(sym, cdq, 2, num_kpts, num_wann) !RS:

      ! save search direction
      cdqkeep_loc(:, :, :) = cdq_loc(:, :, :)

      ! check whether we're doing fixed step lengths
      if (wannierise%control%lfixstep) then

        alphamin = wannierise%control%fixed_step

        ! or a parabolic line search
      else

        ! take trial step
        cdq_loc(:, :, :) = cdqkeep_loc(:, :, :)*(wannierise%control%trial_step/(4.0_dp*kmesh_info%wbtot))

        ! store original U and M before rotating
        u0_loc = u_matrix_loc

        if (verbose%optimisation <= 0) then
!             write(page_unit)   m_matrix
          write (page_unit) m_matrix_loc
          rewind (page_unit)
        else
          m0_loc = m_matrix_loc
        endif

        ! update U and M
        call internal_new_u_and_m(cdq, cmtmp, tmp_cdq, cwork, rwork, evals, cwschur1, cwschur2, &
                                  cwschur3, cwschur4, cz, num_wann, num_kpts, kmesh_info, &
                                  lsitesymmetry, counts, displs, cdq_loc, u_matrix_loc, &
                                  m_matrix_loc, verbose%timing_level, stdout, sym, comm)

        ! calculate spread at trial step
        call wann_omega(csheet, sheet, rave, r2ave, rave2, trial_spread, num_wann, kmesh_info, &
                        num_kpts, verbose, wannierise%constrain, &
                        wannierise%omega%invariant, counts, displs, ln_tmp_loc, &
                        m_matrix_loc, lambda_loc, first_pass, stdout, seedname, comm)

        ! Calculate optimal step (alphamin)
        call internal_optimal_step(wann_spread, trial_spread, doda0, alphamin, falphamin, lquad, &
                                   lprint, wannierise%control%trial_step, stdout)
      endif

      ! print line search information
      if (lprint .and. verbose%iprint > 2) then
        write (stdout, *) ' LINE --> Spread at initial point       :', &
          wann_spread%om_tot*verbose%lenconfac**2
        if (.not. wannierise%control%lfixstep) &
          write (stdout, *) ' LINE --> Spread at trial step          :', &
          trial_spread%om_tot*verbose%lenconfac**2
        write (stdout, *) ' LINE --> Slope along search direction  :', &
          doda0*verbose%lenconfac**2
        write (stdout, *) ' LINE --> ||SD gradient||^2             :', &
          gcnorm1*verbose%lenconfac**2
        if (.not. wannierise%control%lfixstep) then
          write (stdout, *) ' LINE --> Trial step length             :', wannierise%control%trial_step
          if (lquad) then
            write (stdout, *) ' LINE --> Optimal parabolic step length :', alphamin
            write (stdout, *) ' LINE --> Spread at predicted minimum   :', &
              falphamin*verbose%lenconfac**2
          endif
        else
          write (stdout, *) ' LINE --> Fixed step length             :', wannierise%control%fixed_step
        endif
        write (stdout, *) ' LINE --> CG coefficient                :', gcfac
      endif

      ! if taking a fixed step or if parabolic line search was successful
      if (wannierise%control%lfixstep .or. lquad) then

        ! take optimal step
        cdq_loc(:, :, :) = cdqkeep_loc(:, :, :)*(alphamin/(4.0_dp*kmesh_info%wbtot))

        ! if doing a line search then restore original U and M before rotating
        if (.not. wannierise%control%lfixstep) then
          u_matrix_loc = u0_loc
          if (verbose%optimisation <= 0) then
!                read(page_unit)  m_matrix
            read (page_unit) m_matrix_loc
            rewind (page_unit)
          else
            m_matrix_loc = m0_loc
          endif
        endif

        ! update U and M
        call internal_new_u_and_m(cdq, cmtmp, tmp_cdq, cwork, rwork, evals, cwschur1, cwschur2, &
                                  cwschur3, cwschur4, cz, num_wann, num_kpts, kmesh_info, &
                                  lsitesymmetry, counts, displs, cdq_loc, u_matrix_loc, &
                                  m_matrix_loc, verbose%timing_level, stdout, sym, comm)

        call wann_spread_copy(wann_spread, old_spread)

        ! calculate the new centers and spread
        call wann_omega(csheet, sheet, rave, r2ave, rave2, wann_spread, num_wann, kmesh_info, &
                        num_kpts, verbose, wannierise%constrain, &
                        wannierise%omega%invariant, counts, displs, ln_tmp_loc, &
                        m_matrix_loc, lambda_loc, first_pass, stdout, seedname, comm)

        ! parabolic line search was unsuccessful, use trial step already taken
      else

        call wann_spread_copy(wann_spread, old_spread)
        call wann_spread_copy(trial_spread, wann_spread)

      endif

      ! print the new centers and spreads
      if (lprint .and. verbose%iprint > 0) then
        do iw = 1, num_wann
          write (stdout, 1000) iw, (rave(ind, iw)*verbose%lenconfac, ind=1, 3), &
            (r2ave(iw) - rave2(iw))*verbose%lenconfac**2
        end do
        write (stdout, 1001) (sum(rave(ind, :))*verbose%lenconfac, ind=1, 3), &
          (sum(r2ave) - sum(rave2))*verbose%lenconfac**2
        write (stdout, *)
        if (wannierise%constrain%selective_loc .and. wannierise%constrain%slwf_constrain) then
          write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
            iter, (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
            sqrt(abs(gcnorm1))*verbose%lenconfac, &
            wann_spread%om_tot*verbose%lenconfac**2, io_wallclocktime(), '<-- CONV'
          write (stdout, '(7x,a,F15.7,a,F15.7,a,F15.7,a)') &
            'O_IOD=', (wann_spread%om_iod + wann_spread%om_nu)*verbose%lenconfac**2, &
            ' O_D=', wann_spread%om_d*verbose%lenconfac**2, &
            ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
          write (stdout, '(a,E15.7,a,E15.7,a,E15.7,a)') &
            'Delta: O_IOD=', ((wann_spread%om_iod + wann_spread%om_nu) - &
                              (old_spread%om_iod + wann_spread%om_nu))*verbose%lenconfac**2, &
            ' O_D=', (wann_spread%om_d - old_spread%om_d)*verbose%lenconfac**2, &
            ' O_TOT=', (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, ' <-- DLTA'
          write (stdout, '(1x,a78)') repeat('-', 78)
        elseif (wannierise%constrain%selective_loc .and. .not. wannierise%constrain%slwf_constrain) then
          write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
            iter, (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
            sqrt(abs(gcnorm1))*verbose%lenconfac, &
            wann_spread%om_tot*verbose%lenconfac**2, io_wallclocktime(), '<-- CONV'
          write (stdout, '(7x,a,F15.7,a,F15.7,a,F15.7,a)') &
            'O_IOD=', wann_spread%om_iod*verbose%lenconfac**2, &
            ' O_D=', wann_spread%om_d*verbose%lenconfac**2, &
            ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
          write (stdout, '(a,E15.7,a,E15.7,a,E15.7,a)') &
            'Delta: O_IOD=', (wann_spread%om_iod - old_spread%om_iod)*verbose%lenconfac**2, &
            ' O_D=', (wann_spread%om_d - old_spread%om_d)*verbose%lenconfac**2, &
            ' O_TOT=', (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, ' <-- DLTA'
          write (stdout, '(1x,a78)') repeat('-', 78)
        else
          write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
            iter, (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
            sqrt(abs(gcnorm1))*verbose%lenconfac, &
            wann_spread%om_tot*verbose%lenconfac**2, io_wallclocktime(), '<-- CONV'
          write (stdout, '(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
            'O_D=', wann_spread%om_d*verbose%lenconfac**2, &
            ' O_OD=', wann_spread%om_od*verbose%lenconfac**2, &
            ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
          write (stdout, '(1x,a,E15.7,a,E15.7,a,E15.7,a)') &
            'Delta: O_D=', (wann_spread%om_d - old_spread%om_d)*verbose%lenconfac**2, &
            ' O_OD=', (wann_spread%om_od - old_spread%om_od)*verbose%lenconfac**2, &
            ' O_TOT=', (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, ' <-- DLTA'
          write (stdout, '(1x,a78)') repeat('-', 78)
        end if
      end if

      ! Public array of Wannier centres and spreads
      wann_data%centres = rave
      wann_data%spreads = r2ave - rave2

      ! Public variables
      if (.not. wannierise%constrain%selective_loc) then
        wannierise%omega%total = wann_spread%om_tot
        wannierise%omega%tilde = wann_spread%om_d + wann_spread%om_od
      else
        wannierise%omega%total = wann_spread%om_tot
        !omega_tilde = wann_spread%om_d + wann_spread%om_nu
      end if

      if (ldump) then
        ! Before calling param_write_chkpt, I need to gather on the root node
        ! the u_matrix from the u_matrix_loc. No need to broadcast it since
        ! it's printed by the root node only
        call comms_gatherv(u_matrix_loc, num_wann*num_wann*counts(my_node_id), &
                           u_matrix, num_wann*num_wann*counts, num_wann*num_wann*displs, stdout, &
                           seedname, comm)
        ! I also transfer the M matrix
        call comms_gatherv(m_matrix_loc, num_wann*num_wann*kmesh_info%nntot*counts(my_node_id), &
                           m_matrix, num_wann*num_wann*kmesh_info%nntot*counts, &
                           num_wann*num_wann*kmesh_info%nntot*displs, stdout, seedname, comm)
        if (on_root) call param_write_chkpt('postdis', excluded_bands, wann_data, kmesh_info, &
                                            k_points, num_kpts, dis_window, num_bands, num_wann, &
                                            u_matrix, u_matrix_opt, m_matrix, mp_grid, &
                                            real_lattice, recip_lattice, &
                                            wannierise%omega%invariant, have_disentangled, &
                                            stdout, seedname)
      endif

      if (wannierise%control%conv_window .gt. 1) then
        call internal_test_convergence(old_spread, wann_spread, history, save_spread, iter, &
                                       conv_count, noise_count, lconverged, lrandom, lfirst, &
                                       wannierise%control, stdout)
      endif

      if (lconverged) then
        write (stdout, '(/13x,a,es10.3,a,i2,a)') &
          '<<<     Delta <', wannierise%control%conv_tol, &
          '  over ', wannierise%control%conv_window, ' iterations     >>>'
        write (stdout, '(13x,a/)') '<<< Wannierisation convergence criteria satisfied >>>'
        exit
      endif

    enddo
    ! end of the minimization loop

    ! the m matrix is sent by piece to avoid huge arrays
    ! But, I want to reduce the memory usage as much as possible.
!    do nn = 1, nntot
!      m_matrix_1b_loc=m_matrix_loc(:,:,nn,:)
!      call comms_gatherv(m_matrix_1b_loc,num_wann*num_wann*counts(my_node_id),&
!                 m_matrix_1b,num_wann*num_wann*counts,num_wann*num_wann*displs)
!      call comms_bcast(m_matrix_1b(1,1,1),num_wann*num_wann*num_kpts)
!      m_matrix(:,:,nn,:)=m_matrix_1b(:,:,:)
!    end do!nn
    call comms_gatherv(m_matrix_loc, num_wann*num_wann*kmesh_info%nntot*counts(my_node_id), &
                       m_matrix, num_wann*num_wann*kmesh_info%nntot*counts, &
                       num_wann*num_wann*kmesh_info%nntot*displs, stdout, seedname, comm)

    ! send u matrix
    call comms_gatherv(u_matrix_loc, num_wann*num_wann*counts(my_node_id), &
                       u_matrix, num_wann*num_wann*counts, num_wann*num_wann*displs, stdout, &
                       seedname, comm)
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, stdout, seedname, comm)

    ! Evaluate the penalty functional
    if (wannierise%constrain%selective_loc .and. wannierise%constrain%slwf_constrain) then
      rnr0n2 = 0.0_dp
      do iw = 1, wannierise%constrain%slwf_num
        rnr0n2(iw) = (wann_data%centres(1, iw) - wannierise%constrain%ccentres_cart(iw, 1))**2 &
                     + (wann_data%centres(2, iw) - wannierise%constrain%ccentres_cart(iw, 2))**2 &
                     + (wann_data%centres(3, iw) - wannierise%constrain%ccentres_cart(iw, 3))**2
      end do
    end if

    if (verbose%iprint > 0) then
      write (stdout, '(1x,a)') 'Final State'
      do iw = 1, num_wann
        write (stdout, 1000) iw, (rave(ind, iw)*verbose%lenconfac, ind=1, 3), &
          (r2ave(iw) - rave2(iw))*verbose%lenconfac**2
      end do
      write (stdout, 1001) (sum(rave(ind, :))*verbose%lenconfac, ind=1, 3), &
        (sum(r2ave) - sum(rave2))*verbose%lenconfac**2
      write (stdout, *)
      if (wannierise%constrain%selective_loc .and. wannierise%constrain%slwf_constrain) then
        write (stdout, '(3x,a21,a,f15.9)') '     Spreads ('//trim(verbose%length_unit)//'^2)', &
          '       Omega IOD_C   = ', (wann_spread%om_iod + wann_spread%om_nu)*verbose%lenconfac**2
        write (stdout, '(3x,a,f15.9)') '     ================       Omega D       = ', &
          wann_spread%om_d*verbose%lenconfac**2
        write (stdout, '(3x,a,f15.9)') '                            Omega Rest    = ', &
          (sum(r2ave) - sum(rave2) + wann_spread%om_tot)*verbose%lenconfac**2
        write (stdout, '(3x,a,f15.9)') '                            Penalty func  = ', &
          sum(rnr0n2(:))
        write (stdout, '(3x,a21,a,f15.9)') 'Final Spread ('//trim(verbose%length_unit)//'^2)', &
          '       Omega Total_C = ', wann_spread%om_tot*verbose%lenconfac**2
        write (stdout, '(1x,a78)') repeat('-', 78)
      else if (wannierise%constrain%selective_loc .and. .not. wannierise%constrain%slwf_constrain) then
        write (stdout, '(3x,a21,a,f15.9)') '     Spreads ('//trim(verbose%length_unit)//'^2)', &
          '       Omega IOD    = ', wann_spread%om_iod*verbose%lenconfac**2
        write (stdout, '(3x,a,f15.9)') '     ================       Omega D      = ', &
          wann_spread%om_d*verbose%lenconfac**2
        write (stdout, '(3x,a,f15.9)') '                            Omega Rest   = ', &
          (sum(r2ave) - sum(rave2) + wann_spread%om_tot)*verbose%lenconfac**2
        write (stdout, '(3x,a21,a,f15.9)') 'Final Spread ('//trim(verbose%length_unit)//'^2)', &
          '       Omega Total  = ', wann_spread%om_tot*verbose%lenconfac**2
        write (stdout, '(1x,a78)') repeat('-', 78)
      else
        write (stdout, '(3x,a21,a,f15.9)') '     Spreads ('//trim(verbose%length_unit)//'^2)', &
          '       Omega I      = ', wann_spread%om_i*verbose%lenconfac**2
        write (stdout, '(3x,a,f15.9)') '     ================       Omega D      = ', &
          wann_spread%om_d*verbose%lenconfac**2
        write (stdout, '(3x,a,f15.9)') '                            Omega OD     = ', &
          wann_spread%om_od*verbose%lenconfac**2
        write (stdout, '(3x,a21,a,f15.9)') 'Final Spread ('//trim(verbose%length_unit)//'^2)', &
          '       Omega Total  = ', wann_spread%om_tot*verbose%lenconfac**2
        write (stdout, '(1x,a78)') repeat('-', 78)
      end if
    endif

    if (out_files%write_xyz .and. on_root) then
      call wann_write_xyz(out_files%translate_home_cell, num_wann, wann_data%centres, &
                          real_lattice, recip_lattice, atoms, verbose, stdout, seedname)
    endif

    if (out_files%write_hr_diag) then
      call hamiltonian_setup(hmlg, verbose, ws_region, w90_calcs, ham_k, ham_r, &
                             real_lattice, wannier_centres_translated, irvec, mp_grid, ndegen, &
                             num_kpts, num_wann, nrpts, rpt_origin, bands_plot_mode, stdout, &
                             seedname, transport_mode)
      call hamiltonian_get_hr(atoms, dis_window, hmlg, param_hamil, verbose, ham_k, ham_r, &
                              u_matrix, u_matrix_opt, eigval, k_points%kpt_latt, real_lattice, &
                              recip_lattice, wann_data%centres, wannier_centres_translated, irvec, &
                              shift_vec, nrpts, num_bands, num_kpts, num_wann, have_disentangled, &
                              stdout, seedname, lsitesymmetry)
      if (verbose%iprint > 0) then
        write (stdout, *)
        write (stdout, '(1x,a)') 'On-site Hamiltonian matrix elements'
        write (stdout, '(3x,a)') '  n        <0n|H|0n> (eV)'
        write (stdout, '(3x,a)') '-------------------------'
        do i = 1, num_wann
          write (stdout, '(3x,i3,5x,f12.6)') i, real(ham_r(i, i, rpt_origin), kind=dp)
        enddo
        write (stdout, *)
      endif
    endif

    if (wannierise%control%guiding_centres) then
      call wann_phases(csheet, sheet, rguide, irguide, num_wann, kmesh_info, num_kpts, m_matrix, &
                       .false., counts, displs, m_matrix_loc, rnkb, verbose%timing_level, &
                       stdout, seedname, verbose%iprint, comm)
    endif

    ! unitarity is checked
!~    call internal_check_unitarity()
    call wann_check_unitarity(num_kpts, num_wann, u_matrix, verbose%timing_level, &
                              verbose%iprint, stdout, seedname)

    ! write extra info regarding omega_invariant
!~    if (iprint>2) call internal_svd_omega_i()
!    if (iprint>2) call wann_svd_omega_i()
    if (verbose%iprint > 2 .and. on_root) then
      call wann_svd_omega_i(num_wann, num_kpts, kmesh_info, m_matrix, verbose, stdout, seedname)
    endif

    ! write matrix elements <m|r^2|n> to file
!~    if (write_r2mn) call internal_write_r2mn()
!    if (write_r2mn) call wann_write_r2mn()
    if (out_files%write_r2mn .and. on_root) call wann_write_r2mn(num_kpts, num_wann, kmesh_info, &
                                                                 m_matrix, stdout, seedname)

    ! calculate and write projection of WFs on original bands in outer window
    if (have_disentangled .and. out_files%write_proj) &
      call wann_calc_projection(num_bands, num_wann, num_kpts, u_matrix_opt, eigval, &
                                dis_window%lwindow, verbose%timing_level, verbose%iprint, &
                                stdout, seedname)

    ! aam: write data required for vdW utility
    if (out_files%write_vdw_data .and. on_root) then
      call wann_write_vdw_data(num_wann, wann_data, real_lattice, recip_lattice, u_matrix, &
                               u_matrix_opt, have_disentangled, system, stdout, seedname)
    endif

    ! deallocate sub vars not passed into other subs
    deallocate (rwork, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rwork in wann_main', stdout, seedname)
    deallocate (cwork, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cwork in wann_main', stdout, seedname)
    deallocate (evals, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating evals in wann_main', stdout, seedname)
    deallocate (tmp_cdq, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating tmp_cdq in wann_main', stdout, seedname)
    deallocate (cmtmp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cmtmp in wann_main', stdout, seedname)
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cz in wann_main', stdout, seedname)
    deallocate (cdq, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cdq in wann_main', stdout, seedname)

    ! for MPI
    deallocate (ln_tmp_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ln_tmp_loc in wann_main', stdout, seedname)
    deallocate (rnkb_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rnkb_loc in wann_main', stdout, seedname)
    deallocate (u_matrix_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating u_matrix_loc in wann_main', stdout, seedname)
    deallocate (m_matrix_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating m_matrix_loc in wann_main', stdout, seedname)
!    deallocate(m_matrix_1b,stat=ierr)
!    if (ierr/=0) call io_error('Error in deallocating m_matrix_1b in wann_main')
!    deallocate(m_matrix_1b_loc,stat=ierr)
!    if (ierr/=0) call io_error('Error in deallocating m_matrix_1b_loc in wann_main')
    deallocate (cdq_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cdq_loc in wann_main', stdout, seedname)
    deallocate (cdodq_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cdodq_loc in wann_main', stdout, seedname)
    deallocate (cdqkeep_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cdqkeep_loc in wann_main', stdout, seedname)

    deallocate (cwschur3, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cwschur3 in wann_main', stdout, seedname)
    deallocate (cwschur1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cwschur1 in wann_main', stdout, seedname)
    if (wannierise%control%precond) then
      if (verbose%optimisation >= 3) then
        deallocate (k_to_r, stat=ierr)
        if (ierr /= 0) call io_error('Error in deallocating k_to_r in wann_main', stdout, seedname)
      end if
      deallocate (cdodq_r, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating cdodq_r in wann_main', stdout, seedname)
      deallocate (cdodq_precond, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating cdodq_precond in wann_main', stdout, seedname)
      deallocate (cdodq_precond_loc, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating cdodq_precond_loc in wann_main', stdout, seedname)
    end if

    ! deallocate sub vars passed into other subs
    deallocate (rguide, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rguide in wann_main', stdout, seedname)
    deallocate (rave2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rave2 in wann_main', stdout, seedname)
    deallocate (rave, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rave in wann_main', stdout, seedname)
    deallocate (sheet, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating sheet in wann_main', stdout, seedname)
    deallocate (cdodq, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cdodq in wann_main', stdout, seedname)
    deallocate (csheet, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating csheet in wann_main', stdout, seedname)
    if (wannierise%constrain%selective_loc) then
      deallocate (rnr0n2, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating rnr0n2 in wann_main', stdout, seedname)
    end if
    ! deallocate module data
    deallocate (ln_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ln_tmp in wann_main', stdout, seedname)
    deallocate (rnkb, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rnkb in wann_main', stdout, seedname)

    deallocate (u0_loc, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating u0_loc in wann_main', stdout, seedname)
    if (verbose%optimisation > 0) then
      deallocate (m0_loc, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating m0_loc in wann_main', stdout, seedname)
    end if

    if (allocated(counts)) deallocate (counts)
    if (allocated(displs)) deallocate (displs)

    deallocate (history, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating history in wann_main', stdout, seedname)

    if (verbose%timing_level > 0 .and. verbose%iprint > 0) call io_stopwatch('wann: main', 2, stdout, seedname)

    return

1000 format(2x, 'WF centre and spread', &
&       i5, 2x, '(', f10.6, ',', f10.6, ',', f10.6, ' )', f15.8)

1001 format(2x, 'Sum of centres and spreads', &
&       1x, '(', f10.6, ',', f10.6, ',', f10.6, ' )', f15.8)

  contains

    !===============================================!
    subroutine internal_test_convergence(old_spread, wann_spread, history, save_spread, iter, &
                                         conv_count, noise_count, lconverged, lrandom, lfirst, &
                                         wann_control, stdout)
      !===============================================!
      !                                               !
      !! Determine whether minimisation of non-gauge
      !! invariant spread is converged
      !                                               !
      !===============================================!
      use w90_io, only: io_error
      use wannier_param_types, only: wann_control_type

      implicit none

      type(wann_control_type), intent(in) :: wann_control

      type(localisation_vars), intent(in) :: old_spread
      type(localisation_vars), intent(in) :: wann_spread
      real(kind=dp), intent(inout) :: history(:)
      real(kind=dp), intent(out) :: save_spread
      integer, intent(in) :: iter
      integer, intent(inout) :: conv_count
      integer, intent(inout) :: noise_count
      logical, intent(inout) :: lconverged, lrandom, lfirst
!     integer, intent(in) :: conv_window
!     real(kind=dp), intent(in) :: conv_tol
!     real(kind=dp), intent(in) :: conv_noise_amp
!     integer, intent(in) :: conv_noise_num
      integer, intent(in) ::stdout
      ! local
      real(kind=dp) :: delta_omega
      integer :: j, ierr
      real(kind=dp), allocatable :: temp_hist(:)

      allocate (temp_hist(wann_control%conv_window), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating temp_hist in wann_main', stdout, seedname)

      delta_omega = wann_spread%om_tot - old_spread%om_tot

      if (iter .le. wann_control%conv_window) then
        history(iter) = delta_omega
      else
        temp_hist = eoshift(history, 1, delta_omega)
        history = temp_hist
      endif

      conv_count = conv_count + 1

      if (conv_count .lt. wann_control%conv_window) then
        return
      else
!~         write(stdout,*) (history(j),j=1,conv_window)
        do j = 1, wann_control%conv_window
          if (abs(history(j)) .gt. wann_control%conv_tol) return
        enddo
      endif

      if ((wann_control%conv_noise_amp .gt. 0.0_dp) .and. &
          (noise_count .lt. wann_control%conv_noise_num)) then
        if (lfirst) then
          lfirst = .false.
          save_spread = wann_spread%om_tot
          lrandom = .true.
          conv_count = 0
        else
          if (abs(save_spread - wann_spread%om_tot) .lt. wann_control%conv_tol) then
            lconverged = .true.
            return
          else
            save_spread = wann_spread%om_tot
            lrandom = .true.
            conv_count = 0
          endif
        endif
      else
        lconverged = .true.
      endif

      if (lrandom) noise_count = noise_count + 1

      deallocate (temp_hist, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating temp_hist in wann_main', stdout, seedname)

      return

    end subroutine internal_test_convergence

    !===============================================!
    subroutine internal_random_noise(conv_noise_amp, num_wann, counts, cdq_loc, stdout)
      !===============================================!
      !                                               !
      !! Add some random noise to the search direction
      !! to help escape from local minima
      !                                               !
      !===============================================!
      use w90_constants, only: cmplx_0
      use w90_io, only: io_error
      use w90_comms, only: w90commtype

      implicit none
      real(kind=dp), intent(in) :: conv_noise_amp
      integer, intent(in) :: num_wann
      complex(kind=dp), intent(inout) :: cdq_loc(:, :, :)
      integer, intent(in) :: counts(0:)
      integer, intent(in) :: stdout
      ! local
      integer :: ikp, iw, jw, ierr
      real(kind=dp), allocatable :: noise_real(:, :), noise_imag(:, :)
      complex(kind=dp), allocatable :: cnoise(:, :)

      ! Allocate
      allocate (noise_real(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating noise_real in wann_main', stdout, seedname)
      allocate (noise_imag(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating noise_imag in wann_main', stdout, seedname)
      allocate (cnoise(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating cnoise in wann_main', stdout, seedname)

      ! Initialise
      cnoise = cmplx_0; noise_real = 0.0_dp; noise_imag = 0.0_dp

      ! cdq is a num_wann x num_wann x num_kpts anti-hermitian array
      ! to which we add a random anti-hermitian matrix

      do ikp = 1, counts(my_node_id)
        do iw = 1, num_wann
          call random_seed()
          call random_number(noise_real(:, iw))
          call random_seed()
          call random_number(noise_imag(:, iw))
        enddo
        do jw = 1, num_wann
          do iw = 1, jw
            if (iw .eq. jw) then
              cnoise(iw, jw) = cmplx(0.0_dp, noise_imag(iw, jw), dp)
            else
              cnoise(iw, jw) = cmplx(noise_real(iw, jw), noise_imag(iw, jw), dp)
            endif
            cnoise(jw, iw) = -conjg(cnoise(iw, jw))
          enddo
        enddo
        ! Add noise to search direction
        cdq_loc(:, :, ikp) = cdq_loc(:, :, ikp) + conv_noise_amp*cnoise(:, :)
      enddo

      ! Deallocate
      deallocate (cnoise, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating cnoise in wann_main', stdout, seedname)
      deallocate (noise_imag, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating noise_imag in wann_main', stdout, seedname)
      deallocate (noise_real, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating noise_real in wann_main', stdout, seedname)

      return

    end subroutine internal_random_noise

    !===============================================!
    subroutine precond_search_direction(cdodq, cdodq_r, cdodq_precond, cdodq_precond_loc, &
                                        k_to_r, wann_spread, num_wann, num_kpts, &
                                        kpt_latt, real_lattice, nrpts, irvec, ndegen, &
                                        counts, displs, stdout)
      !===============================================!
      !                                               !
      !! Calculate the conjugate gradients search
      !! direction using the Fletcher-Reeves formula:
      !!
      !!     cg_coeff = [g(i).g(i)]/[g(i-1).g(i-1)]
      !                                               !
      !===============================================!
      use w90_constants, only: cmplx_0, cmplx_1, cmplx_i, twopi
      use w90_io, only: io_stopwatch
      !use w90_comms, only: comms_allreduce, w90commtype

      implicit none

      complex(kind=dp), intent(in) :: cdodq(:, :, :)
      complex(kind=dp), intent(inout) :: cdodq_r(:, :, :)
      complex(kind=dp), intent(inout) :: cdodq_precond(:, :, :)
      complex(kind=dp), intent(inout) :: cdodq_precond_loc(:, :, :)
      !complex(kind=dp), intent(inout) :: cdqkeep_loc(:, :, :)
      ! k_to_r depends on optimisation flag
      complex(kind=dp), allocatable, intent(in) :: k_to_r(:, :)
      type(localisation_vars), intent(in) :: wann_spread
      !integer, intent(in) :: iter
      !logical, intent(in) :: lprint
      !logical, intent(inout) :: lrandom
      integer, intent(in) :: num_wann, num_kpts
      real(kind=dp), intent(in) :: kpt_latt(:, :)
      real(kind=dp), intent(in) :: real_lattice(3, 3)
      integer, intent(in) :: nrpts
      integer, intent(in) :: irvec(:, :)
      integer, intent(in) :: ndegen(:)
      integer, intent(in) :: counts(0:)
      integer, intent(in) :: displs(0:)
      integer, intent(in) :: stdout
      !type(w90commtype), intent(in) :: comm

      ! local
      complex(kind=dp), external :: zdotc
      complex(kind=dp) :: fac, rdotk
      real(kind=dp), dimension(3) :: rvec_cart
      real(kind=dp) :: alpha_precond
      integer :: irpt, loop_kpt

      if (verbose%timing_level > 1 .and. verbose%iprint > 0) &
        call io_stopwatch('wann: main: search_direction', 1, stdout, seedname)

      ! gcnorm1 = Tr[gradient . gradient] -- NB gradient is anti-Hermitian
      ! gcnorm1 = real(zdotc(num_kpts*num_wann*num_wann,cdodq,1,cdodq,1),dp)

      !if (wannierise%precond) then
      ! compute cdodq_precond

      cdodq_r(:, :, :) = 0 ! intermediary gradient in R space
      cdodq_precond(:, :, :) = 0
      cdodq_precond_loc(:, :, :) = 0
!         cdodq_precond(:,:,:) = complx_0

      ! convert to real space in cdodq_r
      ! Two algorithms: either double loop or GEMM. GEMM is much more efficient but requires more RAM
      ! Ideally, we should implement FFT-based filtering here
      if (verbose%optimisation >= 3) then
        call zgemm('N', 'N', num_wann*num_wann, nrpts, num_kpts, cmplx_1, cdodq, &
                   num_wann*num_wann, k_to_r, num_kpts, cmplx_0, cdodq_r, num_wann*num_wann)
        cdodq_r = cdodq_r/real(num_kpts, dp)
      else
        do irpt = 1, nrpts
          do loop_kpt = 1, num_kpts
            rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
            fac = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
            cdodq_r(:, :, irpt) = cdodq_r(:, :, irpt) + fac*cdodq(:, :, loop_kpt)
          enddo
        enddo
      end if

      ! filter cdodq_r in real space by 1/(1+R^2/alpha)

      ! this alpha coefficient is more or less arbitrary, and could
      ! be tweaked further: the point is to have something that has
      ! the right units, and is not too small (or the filtering is
      ! too severe) or too high (or the filtering does nothing).
      !
      ! the descent direction produced has a different magnitude
      ! than the one without preconditionner, so the values of
      ! trial_step are not consistent
      alpha_precond = 10*wann_spread%om_tot/num_wann
      do irpt = 1, nrpts
        rvec_cart = matmul(real_lattice(:, :), real(irvec(:, irpt), dp))
        cdodq_r(:, :, irpt) = cdodq_r(:, :, irpt)*1/(1 + dot_product(rvec_cart, rvec_cart)/ &
                                                     alpha_precond)
      end do

      ! go back to k space
      if (verbose%optimisation >= 3) then
        do irpt = 1, nrpts
          cdodq_r(:, :, irpt) = cdodq_r(:, :, irpt)/real(ndegen(irpt), dp)
        end do
        call zgemm('N', 'C', num_wann*num_wann, num_kpts, nrpts, cmplx_1, &
            & cdodq_r, num_wann*num_wann, k_to_r, num_kpts, cmplx_0, cdodq_precond, &
            num_wann*num_wann)
      else
        do irpt = 1, nrpts
          do loop_kpt = 1, num_kpts
            rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
            fac = exp(cmplx_i*rdotk)/real(ndegen(irpt), dp)
            cdodq_precond(:, :, loop_kpt) = cdodq_precond(:, :, loop_kpt) + &
                                            fac*cdodq_r(:, :, irpt)
          enddo
        enddo
      end if
      cdodq_precond_loc(:, :, 1:counts(my_node_id)) = &
        cdodq_precond(:, :, 1 + displs(my_node_id):displs(my_node_id) + counts(my_node_id))

      !end if

    end subroutine precond_search_direction

    !===============================================!
    subroutine internal_search_direction(cdodq_precond_loc, cdqkeep_loc, iter, lprint, lrandom, &
                                         noise_count, ncg, gcfac, gcnorm0, gcnorm1, &
                                         doda0, wann_control, num_wann, &
                                         wbtot, cdq_loc, cdodq_loc, counts, stdout)
      !===============================================!
      !                                               !
      !! Calculate the conjugate gradients search
      !! direction using the Fletcher-Reeves formula:
      !!
      !!     cg_coeff = [g(i).g(i)]/[g(i-1).g(i-1)]
      !                                               !
      !===============================================!
      !use w90_constants, only: cmplx_0, cmplx_1, cmplx_i, twopi
      use w90_io, only: io_stopwatch
      use w90_comms, only: comms_allreduce, w90commtype
      use wannier_param_types, only: wann_control_type

      implicit none

      type(wann_control_type), intent(in) :: wann_control

      complex(kind=dp), allocatable, intent(inout) :: cdodq_precond_loc(:, :, :)
      complex(kind=dp), intent(inout) :: cdqkeep_loc(:, :, :)
      integer, intent(in) :: iter
      logical, intent(in) :: lprint
      logical, intent(inout) :: lrandom
      integer, intent(in) :: noise_count
      integer, intent(inout) :: ncg
      real(kind=dp), intent(out) :: gcfac
      real(kind=dp), intent(inout) :: gcnorm0, gcnorm1
      real(kind=dp), intent(out) :: doda0
      integer, intent(in) :: num_wann
      !real(kind=dp), intent(in) :: kpt_latt(:, :)
      !real(kind=dp), intent(in) :: real_lattice(3, 3)
      real(kind=dp), intent(in) :: wbtot
      !integer, intent(in) :: nrpts
      !integer, intent(in) :: irvec(:, :)
      !integer, intent(in) :: ndegen(:)
      complex(kind=dp), intent(inout) :: cdq_loc(:, :, :)
      complex(kind=dp), intent(in) :: cdodq_loc(:, :, :)
      integer, intent(in) :: counts(0:)
      integer, intent(in) :: stdout
      ! local
      complex(kind=dp), external :: zdotc

      if ((.not. wann_control%precond) .and. verbose%timing_level > 1 .and. verbose%iprint > 0) &
        call io_stopwatch('wann: main: search_direction', 1, stdout, seedname)

      ! gcnorm1 = Tr[gradient . gradient] -- NB gradient is anti-Hermitian
      if (wann_control%precond) then
!         gcnorm1 = real(zdotc(num_kpts*num_wann*num_wann,cdodq_precond,1,cdodq,1),dp)
        gcnorm1 = real(zdotc(counts(my_node_id)*num_wann*num_wann, cdodq_precond_loc, 1, &
                             cdodq_loc, 1), dp)
      else
        gcnorm1 = real(zdotc(counts(my_node_id)*num_wann*num_wann, cdodq_loc, 1, cdodq_loc, 1), dp)
      end if
      call comms_allreduce(gcnorm1, 1, 'SUM', stdout, seedname, comm)

      ! calculate cg_coefficient
      if ((iter .eq. 1) .or. (ncg .ge. wann_control%num_cg_steps)) then
        gcfac = 0.0_dp                 ! Steepest descents
        ncg = 0
      else
        if (gcnorm0 .gt. epsilon(1.0_dp)) then
          gcfac = gcnorm1/gcnorm0     ! Fletcher-Reeves CG coefficient
          ! prevent CG coefficient from getting too large
          if (gcfac .gt. 3.0_dp) then
            if (lprint .and. verbose%iprint > 2) &
              write (stdout, *) ' LINE --> CG coeff too large. Resetting :', gcfac
            gcfac = 0.0_dp
            ncg = 0
          else
            ncg = ncg + 1
          endif
        else
          gcfac = 0.0_dp
          ncg = 0
        endif
      endif

      ! save for next iteration
      gcnorm0 = gcnorm1

      ! calculate search direction

      if (wann_control%precond) then
        cdq_loc(:, :, :) = cdodq_precond_loc(:, :, :) + cdqkeep_loc(:, :, :)*gcfac !! JRY not MPI
      else
        cdq_loc(:, :, :) = cdodq_loc(:, :, :) + cdqkeep_loc(:, :, :)*gcfac
      end if

      ! add some random noise to search direction, if required
      if (lrandom) then
        if (verbose%iprint > 0) write (stdout, '(a,i3,a,i3,a)') &
          ' [ Adding random noise to search direction. Time ', noise_count, ' / ', &
          wann_control%conv_noise_num, ' ]'
        call internal_random_noise(wann_control%conv_noise_amp, num_wann, counts, cdq_loc, &
                                   stdout)
      endif
      ! calculate gradient along search direction - Tr[gradient . search direction]
      ! NB gradient is anti-hermitian
      doda0 = -real(zdotc(counts(my_node_id)*num_wann*num_wann, cdodq_loc, 1, cdq_loc, 1), dp)

      call comms_allreduce(doda0, 1, 'SUM', stdout, seedname, comm)

      doda0 = doda0/(4.0_dp*wbtot)

      ! check search direction is not uphill
      if (doda0 .gt. 0.0_dp) then
        ! if doing a CG step then reset CG
        if (ncg .gt. 0) then
          if (lprint .and. verbose%iprint > 2 .and. verbose%iprint > 0) &
            write (stdout, *) ' LINE --> Search direction uphill: resetting CG'
          cdq_loc(:, :, :) = cdodq_loc(:, :, :)
          if (lrandom) call internal_random_noise(wann_control%conv_noise_amp, num_wann, &
                                                  counts, cdq_loc, stdout)
          ncg = 0
          gcfac = 0.0_dp
          ! re-calculate gradient along search direction
          doda0 = -real(zdotc(counts(my_node_id)*num_wann*num_wann, cdodq_loc, 1, cdq_loc, 1), dp)

          call comms_allreduce(doda0, 1, 'SUM', stdout, seedname, comm)

          doda0 = doda0/(4.0_dp*wbtot)
          ! if search direction still uphill then reverse search direction
          if (doda0 .gt. 0.0_dp) then
            if (lprint .and. verbose%iprint > 2 .and. verbose%iprint > 0) &
              write (stdout, *) ' LINE --> Search direction still uphill: reversing'
            cdq_loc(:, :, :) = -cdq_loc(:, :, :)
            doda0 = -doda0
          endif
          ! if doing a SD step then reverse search direction
        else
          if (lprint .and. verbose%iprint > 2 .and. verbose%iprint > 0) &
            write (stdout, *) ' LINE --> Search direction uphill: reversing'
          cdq_loc(:, :, :) = -cdq_loc(:, :, :)
          doda0 = -doda0
        endif
      endif

      !~     ! calculate search direction
      !~     cdq(:,:,:) = cdodq(:,:,:) + cdqkeep(:,:,:) * gcfac

      if (verbose%timing_level > 1 .and. verbose%iprint > 0) &
        call io_stopwatch('wann: main: search_direction', 2, stdout, seedname)

      lrandom = .false.

      return

    end subroutine internal_search_direction

    !===============================================!
    subroutine internal_optimal_step(wann_spread, trial_spread, doda0, alphamin, falphamin, lquad, &
                                     lprint, trial_step, stdout)
      !===============================================!
      !                                               !
      !! Calculate the optimal step length based on a
      !! parabolic line search
      !                                               !
      !===============================================!
      use w90_io, only: io_stopwatch
      use w90_comms, only: w90commtype

      implicit none

      type(localisation_vars), intent(in) :: wann_spread
      type(localisation_vars), intent(in) :: trial_spread
      real(kind=dp), intent(in) :: doda0
      real(kind=dp), intent(out) :: alphamin, falphamin
      logical, intent(out) :: lquad
      logical, intent(in) :: lprint
      real(kind=dp), intent(in) :: trial_step
!     integer, intent(in) :: iprint, timing_level
      integer, intent(in) :: stdout
      ! local
      real(kind=dp) :: fac, shift, eqa, eqb

      if (verbose%timing_level > 1 .and. verbose%iprint > 0) &
        call io_stopwatch('wann: main: optimal_step', 1, stdout, seedname)

      fac = trial_spread%om_tot - wann_spread%om_tot
      if (abs(fac) .gt. tiny(1.0_dp)) then
        fac = 1.0_dp/fac
        shift = 1.0_dp
      else
        fac = 1.0e6_dp
        shift = fac*trial_spread%om_tot - fac*wann_spread%om_tot
      endif
      eqb = fac*doda0
      eqa = shift - eqb*trial_step
      if (abs(eqa/(fac*wann_spread%om_tot)) .gt. epsilon(1.0_dp)) then
        lquad = .true.
        alphamin = -0.5_dp*eqb/eqa*(trial_step**2)
        falphamin = wann_spread%om_tot &
                    - 0.25_dp*eqb*eqb/(fac*eqa)*(trial_step**2)
      else
        if (lprint .and. verbose%iprint > 2) write (stdout, *) &
          ' LINE --> Parabolic line search unstable: using trial step'
        lquad = .false.
        alphamin = trial_step
        falphamin = trial_spread%om_tot
      endif

      if (doda0*alphamin .gt. 0.0_dp) then
        if (lprint .and. verbose%iprint > 2) write (stdout, *) &
          ' LINE --> Line search unstable : using trial step'
        lquad = .false.
        alphamin = trial_step
        falphamin = trial_spread%om_tot
      endif

      if (verbose%timing_level > 1 .and. verbose%iprint > 0) &
        call io_stopwatch('wann: main: optimal_step', 2, stdout, seedname)

      return

    end subroutine internal_optimal_step

    !===============================================!
    subroutine internal_new_u_and_m(cdq, cmtmp, tmp_cdq, cwork, rwork, evals, cwschur1, cwschur2, &
                                    cwschur3, cwschur4, cz, num_wann, num_kpts, kmesh_info, &
                                    lsitesymmetry, counts, displs, cdq_loc, u_matrix_loc, &
                                    m_matrix_loc, timing_level, stdout, sym, comm)
      !===============================================!
      !                                               !
      !! Update U and M matrices after a trial step
      !                                               !
      !===============================================!
      use w90_constants, only: cmplx_i
      use w90_sitesym, only: sitesym_symmetrize_rotation, sitesym_data
      use w90_io, only: io_stopwatch, io_error
      use w90_comms, only: comms_bcast, comms_gatherv, w90commtype
      use w90_utility, only: utility_zgemm
      use w90_param_types, only: kmesh_info_type

      implicit none

      type(kmesh_info_type), intent(in) :: kmesh_info

      type(sitesym_data), intent(in) :: sym
      complex(kind=dp), intent(inout) :: cdq(:, :, :)
      complex(kind=dp), intent(inout) :: cmtmp(:, :), tmp_cdq(:, :) ! really just local?
      complex(kind=dp), intent(inout) :: cwork(:)
      real(kind=dp), intent(inout) :: evals(:)
      real(kind=dp), intent(inout) :: rwork(:)
      complex(kind=dp), intent(inout) :: cwschur1(:), cwschur2(:)
      complex(kind=dp), intent(inout) :: cwschur3(:), cwschur4(:)
      complex(kind=dp), intent(inout) :: cz(:, :)
      integer, intent(in) :: num_wann, num_kpts
      logical, intent(in) :: lsitesymmetry
      integer, intent(in) :: counts(0:)
      integer, intent(in) :: displs(0:)
      complex(kind=dp), intent(inout) :: cdq_loc(:, :, :)
      complex(kind=dp), intent(inout) :: u_matrix_loc(:, :, :)
      complex(kind=dp), intent(inout) :: m_matrix_loc(:, :, :, :)
      integer, intent(in) :: timing_level
      integer, intent(in) :: stdout
      type(w90commtype), intent(in) :: comm

      ! local vars
      integer :: i, nkp, nn, nkp2, nsdim, nkp_loc, info
      logical :: ltmp
      integer :: my_node_id !absence of this was not caught by ifort (!!) why -- fixme JJ

      my_node_id = mpirank(comm)

      if (timing_level > 1 .and. verbose%iprint > 0) call io_stopwatch('wann: main: u_and_m', 1, stdout, seedname)

      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        if (lsitesymmetry) then                !YN: RS:
          if (sym%ir2ik(sym%ik2ir(nkp)) .ne. nkp) cycle !YN: RS:
        end if                                 !YN: RS:
        ! cdq(nkp) is anti-Hermitian; tmp_cdq = i*cdq  is Hermitian
        tmp_cdq(:, :) = cmplx_i*cdq_loc(:, :, nkp_loc)
        ! Hermitian matrix eigen-solver
        call zheev('V', 'U', num_wann, tmp_cdq, num_wann, evals, cwork, 4*num_wann, rwork, info)
        if (info .ne. 0) then
          if (verbose%iprint > 0) write (stdout, *) &
            'wann_main: ZHEEV in internal_new_u_and_m failed, info= ', info
          if (verbose%iprint > 0) write (stdout, *) '           trying Schur decomposition instead'
!!$            call io_error('wann_main: problem in ZHEEV in internal_new_u_and_m')
          tmp_cdq(:, :) = cdq_loc(:, :, nkp_loc)
          call zgees('V', 'N', ltmp, num_wann, tmp_cdq, num_wann, nsdim, &
                     cwschur1, cz, num_wann, cwschur2, 10*num_wann, cwschur3, &
                     cwschur4, info)
          if (info .ne. 0) then
            if (verbose%iprint > 0) write (stdout, *) 'wann_main: SCHUR failed, info= ', info
            call io_error('wann_main: problem computing schur form 1', stdout, seedname)
          endif
          do i = 1, num_wann
            tmp_cdq(:, i) = cz(:, i)*exp(cwschur1(i))
          enddo
          ! cmtmp   = tmp_cdq . cz^{dagger}
          call utility_zgemm(cmtmp, tmp_cdq, 'N', cz, 'C', num_wann)
          cdq_loc(:, :, nkp_loc) = cmtmp(:, :)
        else
          do i = 1, num_wann
            cmtmp(:, i) = tmp_cdq(:, i)*exp(-cmplx_i*evals(i))
          enddo
          ! cdq(nkp)   = cmtmp . tmp_cdq^{dagger}
          call utility_zgemm(cdq_loc(:, :, nkp_loc), cmtmp, 'N', tmp_cdq, 'C', num_wann)
        endif
      enddo

      ! each process communicates its result to other processes
      ! it would be enough to copy only next neighbors
      call comms_gatherv(cdq_loc, num_wann*num_wann*counts(my_node_id), cdq, &
                         num_wann*num_wann*counts, num_wann*num_wann*displs, stdout, seedname, comm)
      call comms_bcast(cdq(1, 1, 1), num_wann*num_wann*num_kpts, stdout, seedname, comm)

!!$      do nkp = 1, num_kpts
!!$         tmp_cdq(:,:) = cdq(:,:,nkp)
!!$         call zgees ('V', 'N', ltmp, num_wann, tmp_cdq, num_wann, nsdim, &
!!$              cwschur1, cz, num_wann, cwschur2, 10 * num_wann, cwschur3, &
!!$              cwschur4, info)
!!$         if (info.ne.0) then
!!$            write(stdout,*) 'SCHUR: ', info
!!$            call io_error('wann_main: problem computing schur form 1')
!!$         endif
!!$         do i=1,num_wann
!!$            tmp_cdq(:,i) = cz(:,i) * exp(cwschur1(i))
!!$         enddo
!!$         ! cmtmp   = tmp_cdq . cz^{dagger}
!!$         call utility_zgemm(cmtmp,tmp_cdq,'N',cz,'C',num_wann)
!!$         cdq(:,:,nkp)=cmtmp(:,:)
!!$      enddo

      if (lsitesymmetry) then
        call sitesym_symmetrize_rotation(sym, cdq, num_kpts, num_wann, seedname, stdout) !RS: calculate cdq(Rk) from k

        cdq_loc(:, :, 1:counts(my_node_id)) = cdq(:, :, 1 + displs(my_node_id):displs(my_node_id) &
                                                  + counts(my_node_id))
      endif

      ! the orbitals are rotated
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        ! cmtmp = U(k) . cdq(k)
        call utility_zgemm(cmtmp, u_matrix_loc(:, :, nkp_loc), 'N', cdq_loc(:, :, nkp_loc), 'N', &
                           num_wann)
        u_matrix_loc(:, :, nkp_loc) = cmtmp(:, :)
      enddo

      ! and the M_ij are updated
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        do nn = 1, kmesh_info%nntot
          nkp2 = kmesh_info%nnlist(nkp, nn)
          ! tmp_cdq = cdq^{dagger} . M
          call utility_zgemm(tmp_cdq, cdq(:, :, nkp), 'C', m_matrix_loc(:, :, nn, nkp_loc), 'N', &
                             num_wann)
          ! cmtmp = tmp_cdq . cdq
          call utility_zgemm(cmtmp, tmp_cdq, 'N', cdq(:, :, nkp2), 'N', num_wann)
          m_matrix_loc(:, :, nn, nkp_loc) = cmtmp(:, :)
        enddo
      enddo

      if (timing_level > 1) call io_stopwatch('wann: main: u_and_m', 2, stdout, seedname)

      return

    end subroutine internal_new_u_and_m

!~    !========================================!
!~    subroutine internal_check_unitarity()
!~    !========================================!
!~
!~      implicit none
!~
!~      integer :: nkp,i,j,m
!~      complex(kind=dp) :: ctmp1,ctmp2
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',1)
!~
!~      do nkp = 1, num_kpts
!~         do i = 1, num_wann
!~            do j = 1, num_wann
!~               ctmp1 = cmplx_0
!~               ctmp2 = cmplx_0
!~               do m = 1, num_wann
!~                  ctmp1 = ctmp1 + u_matrix (i, m, nkp) * conjg (u_matrix (j, m, nkp) )
!~                  ctmp2 = ctmp2 + u_matrix (m, j, nkp) * conjg (u_matrix (m, i, nkp) )
!~               enddo
!~               if ( (i.eq.j) .and. (abs (ctmp1 - cmplx_1 ) .gt. eps5) ) &
!~                    then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp1
!~                  call io_error('wann_main: unitariety error 1')
!~               endif
!~               if ( (i.eq.j) .and. (abs (ctmp2 - cmplx_1 ) .gt. eps5) ) &
!~                    then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp2
!~                  call io_error('wann_main: unitariety error 2')
!~               endif
!~               if ( (i.ne.j) .and. (abs (ctmp1) .gt. eps5) ) then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp1
!~                  call io_error('wann_main: unitariety error 3')
!~               endif
!~               if ( (i.ne.j) .and. (abs (ctmp2) .gt. eps5) ) then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp2
!~                  call io_error('wann_main: unitariety error 4')
!~               endif
!~            enddo
!~         enddo
!~      enddo
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',2)
!~
!~      return
!~
!~    end subroutine internal_check_unitarity

!~    !========================================!
!~    subroutine internal_write_r2mn()
!~    !========================================!
!~    !                                        !
!~    ! Write seedname.r2mn file               !
!~    !                                        !
!~    !========================================!
!~      use w90_io, only: seedname,io_file_unit,io_error
!~
!~      implicit none
!~
!~      integer :: r2mnunit,nw1,nw2,nkp,nn
!~      real(kind=dp) :: r2ave_mn,delta
!~
!~      ! note that here I use formulas analogue to Eq. 23, and not to the
!~      ! shift-invariant Eq. 32 .
!~      r2mnunit=io_file_unit()
!~      open(r2mnunit,file=trim(seedname)//'.r2mn',form='formatted',err=158)
!~      do nw1 = 1, num_wann
!~         do nw2 = 1, num_wann
!~            r2ave_mn = 0.0_dp
!~            delta = 0.0_dp
!~            if (nw1.eq.nw2) delta = 1.0_dp
!~            do nkp = 1, num_kpts
!~               do nn = 1, nntot
!~                  r2ave_mn = r2ave_mn + wb(nn) * &
!~                       ! [GP-begin, Apr13, 2012: corrected sign inside "real"]
!~                       ( 2.0_dp * delta - real(m_matrix(nw1,nw2,nn,nkp) + &
!~                       conjg(m_matrix(nw2,nw1,nn,nkp)),kind=dp) )
!~                       ! [GP-end]
!~               enddo
!~            enddo
!~            r2ave_mn = r2ave_mn / real(num_kpts,dp)
!~            write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn
!~         enddo
!~      enddo
!~      close(r2mnunit)
!~
!~      return
!~
!~158   call io_error('Error opening file '//trim(seedname)//'.r2mn in wann_main')
!~
!~    end subroutine internal_write_r2mn

!~    !========================================!
!~    subroutine internal_svd_omega_i()
!~    !========================================!
!~
!~      implicit none
!~
!~      complex(kind=dp), allocatable  :: cv1(:,:),cv2(:,:)
!~      complex(kind=dp), allocatable  :: cw1(:),cw2(:)
!~      complex(kind=dp), allocatable  :: cpad1 (:)
!~      real(kind=dp),    allocatable  :: singvd (:)
!~
!~      integer :: nkp,nn,nb,na,ind
!~      real(kind=dp) :: omt1,omt2,omt3
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',1)
!~
!~      allocate( cw1 (10 * num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cw1 in wann_main')
!~      allocate( cw2 (10 * num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cw2 in wann_main')
!~      allocate( cv1 (num_wann, num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cv1 in wann_main')
!~      allocate( cv2 (num_wann, num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cv2 in wann_main')
!~      allocate( singvd (num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating singvd in wann_main')
!~      allocate( cpad1 (num_wann * num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cpad1 in wann_main')
!~
!~      cw1=cmplx_0; cw2=cmplx_0; cv1=cmplx_0; cv2=cmplx_0; cpad1=cmplx_0
!~      singvd=0.0_dp
!~
!~      ! singular value decomposition
!~      omt1 = 0.0_dp ; omt2 = 0.0_dp ; omt3 = 0.0_dp
!~      do nkp = 1, num_kpts
!~         do nn = 1, nntot
!~            ind = 1
!~            do nb = 1, num_wann
!~               do na = 1, num_wann
!~                  cpad1 (ind) = m_matrix (na, nb, nn, nkp)
!~                  ind = ind+1
!~               enddo
!~            enddo
!~            call zgesvd ('A', 'A', num_wann, num_wann, cpad1, num_wann, singvd, cv1, &
!~                 num_wann, cv2, num_wann, cw1, 10 * num_wann, cw2, info)
!~            if (info.ne.0) then
!~               call io_error('ERROR: Singular value decomp. zgesvd failed')
!~            endif
!~
!~            do nb = 1, num_wann
!~               omt1 = omt1 + wb(nn) * (1.0_dp - singvd (nb) **2)
!~               omt2 = omt2 - wb(nn) * (2.0_dp * log (singvd (nb) ) )
!~               omt3 = omt3 + wb(nn) * (acos (singvd (nb) ) **2)
!~            enddo
!~         enddo
!~      enddo
!~      omt1 = omt1 / real(num_kpts,dp)
!~      omt2 = omt2 / real(num_kpts,dp)
!~      omt3 = omt3 / real(num_kpts,dp)
!~      write ( stdout , * ) ' '
!~      write(stdout,'(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ',&
!~           omt1*lenconfac**2,'('//trim(length_unit)//'^2)'
!~      write(stdout,'(2x,a,f15.9,1x,a)') '                 -2log s = ',&
!~           omt2*lenconfac**2,'('//trim(length_unit)//'^2)'
!~      write(stdout,'(2x,a,f15.9,1x,a)') '                  acos^2 = ',&
!~           omt3*lenconfac**2,'('//trim(length_unit)//'^2)'
!~
!~      deallocate(cpad1,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cpad1 in wann_main')
!~      deallocate(singvd,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating singvd in wann_main')
!~      deallocate(cv2,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cv2 in wann_main')
!~      deallocate(cv1,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cv1 in wann_main')
!~      deallocate(cw2,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cw2 in wann_main')
!~      deallocate(cw1,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cw1 in wann_main')
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',2)
!~
!~      return
!~
!~    end subroutine internal_svd_omega_i

  end subroutine wann_main

  !==================================================================!

  subroutine wann_phases(csheet, sheet, rguide, irguide, num_wann, kmesh_info, num_kpts, m_matrix, &
                         gamma_only, counts, displs, m_matrix_loc, rnkb, timing_level, stdout, &
                         seedname, iprint, comm, m_w)
    !==================================================================!
    !! Uses guiding centres to pick phases which give a
    !! consistent choice of branch cut for the spread definition
    !                                                                  !
    !===================================================================
    use w90_constants, only: eps6, cmplx_0, cmplx_i
    use w90_io, only: io_stopwatch
    use w90_utility, only: utility_inv3
    use w90_comms, only: comms_allreduce, w90commtype, mpirank
    use w90_param_types, only: kmesh_info_type

    implicit none

    ! passed variables
    type(w90commtype), intent(in) :: comm
    type(kmesh_info_type), intent(in) :: kmesh_info

    integer, intent(in) :: timing_level
    integer, intent(in) :: stdout
    integer, intent(in) :: num_wann
    integer, intent(in) :: num_kpts
    integer, intent(in) :: irguide !! Zero if first call to this routine
    integer, intent(in) :: iprint
    integer, intent(in) :: displs(0:)
    integer, intent(in) :: counts(0:)

    real(kind=dp), intent(out) :: sheet(:, :, :) !! Choice of branch cut
    real(kind=dp), intent(out) :: rnkb(:, :, :)
    real(kind=dp), intent(inout) :: rguide(:, :) !! Guiding centres
    real(kind=dp), intent(in), optional :: m_w(:, :, :)

    complex(kind=dp), intent(out) :: csheet(:, :, :) !! Choice of phase
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    complex(kind=dp), allocatable, intent(in) :: m_matrix_loc(:, :, :, :)

    character(len=50), intent(in)  :: seedname

    logical, intent(in) :: gamma_only

    !local
    complex(kind=dp) :: csum(kmesh_info%nnh)
    real(kind=dp)    ::  xx(kmesh_info%nnh)
    real(kind=dp)    :: smat(3, 3), svec(3), sinv(3, 3)
    real(kind=dp)    :: xx0, det, brn
    complex(kind=dp) :: csumt
    integer :: loop_wann, na, nkp, i, j, nn, ind, m, nkp_loc
    integer :: my_node_id

    my_node_id = mpirank(comm)

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: phases', 1, stdout, seedname)

    csum = cmplx_0; xx = 0.0_dp

    ! report problem to solve
    ! for each band, csum is determined and then its appropriate
    ! guiding center rguide(3,nwann)

    do loop_wann = 1, num_wann

      if (.not. present(m_w)) then
        ! get average phase for each unique bk direction
        if (gamma_only) then
          do na = 1, kmesh_info%nnh
            csum(na) = cmplx_0
            do nkp_loc = 1, counts(my_node_id)
              nkp = nkp_loc + displs(my_node_id)
              nn = kmesh_info%neigh(nkp, na)
              csum(na) = csum(na) + m_matrix(loop_wann, loop_wann, nn, nkp_loc)
            enddo
          enddo
        else
          do na = 1, kmesh_info%nnh
            csum(na) = cmplx_0
            do nkp_loc = 1, counts(my_node_id)
              nkp = nkp_loc + displs(my_node_id)
              nn = kmesh_info%neigh(nkp, na)
              csum(na) = csum(na) + m_matrix_loc(loop_wann, loop_wann, nn, nkp_loc)
            enddo
          enddo
        endif

      else

        do na = 1, kmesh_info%nnh
          csum(na) = cmplx_0
          do nkp_loc = 1, counts(my_node_id)
            nkp = nkp_loc + displs(my_node_id)
            nn = kmesh_info%neigh(nkp, na)
            csum(na) = csum(na) &
                       + cmplx(m_w(loop_wann, loop_wann, 2*nn - 1), m_w(loop_wann, loop_wann, 2*nn), dp)
          enddo
        enddo

      end if

      call comms_allreduce(csum(1), kmesh_info%nnh, 'SUM', stdout, seedname, comm)

      ! now analyze that information to get good guess at
      ! wannier center
      !      write(*,*)
      !      do na=1,nnh
      !       write(*,'a,3f10.5,a,2f10.5)')
      !    &    ' bka=',(bka(j,na),j=1,3),'  csum=',csum(na)
      !      end do
      ! problem is to find a real-space 3-vector rguide such that
      !   phase of csum(nn) ~= phase of exp[ -i bka(nn) dot rguide ]
      ! or, letting
      !   xx(nn) = - Im ln csum(nn)  (modulo 2*pi)
      ! then
      !   bka(nn) dot rguide ~= xx(nn)
      !
      ! we take an arbitrary branch cut for first three xx(nn)
      ! and determine rguide from these; then for each additional bka
      ! vector, we first determine the most consistent branch cut,
      ! and then update rguide
      !
      ! rguide is obtained by minimizing
      !   sum_nn [ bka(nn) dot rguide - xx(nn) ] ^2
      ! or, setting the derivative with respect to rcenter to zero,
      !   sum_i smat(j,i) * rguide(i,nwann) = svec(j)
      ! where
      !   smat(j,i) = sum_nn bka(j,nn) * bka(i,nn)
      !   svec(j)   = sum_nn bka(j,nn) * xx(nn)
      ! initialize smat and svec

      smat = 0.0_dp
      svec = 0.0_dp

      do nn = 1, kmesh_info%nnh
        if (nn .le. 3) then
          !         obtain xx with arbitrary branch cut choice
          xx(nn) = -aimag(log(csum(nn)))
        else
          !         obtain xx with branch cut choice guided by rguide
          xx0 = 0.0_dp
          do j = 1, 3
            xx0 = xx0 + kmesh_info%bka(j, nn)*rguide(j, loop_wann)
          enddo
          !         xx0 is expected value for xx
!             csumt = exp (ci * xx0)
          csumt = exp(cmplx_i*xx0)
          !         csumt has opposite of expected phase of csum(nn)
          xx(nn) = xx0 - aimag(log(csum(nn)*csumt))
        endif

        !       write(*,'(a,i5,3f7.3,2f10.5)') 'nn, bka, xx, mag =',
        !    1    nn,(bka(j,nn),j=1,3),xx(nn),abs(csum(nn))/float(num_kpts)
        !       update smat and svec
        do j = 1, 3
          do i = 1, 3
            smat(j, i) = smat(j, i) + kmesh_info%bka(j, nn)*kmesh_info%bka(i, nn)
          enddo
          svec(j) = svec(j) + kmesh_info%bka(j, nn)*xx(nn)
        enddo

        if (nn .ge. 3) then
          !         determine rguide
          call utility_inv3(smat, sinv, det)
          !         the inverse of smat is sinv/det
          if (abs(det) .gt. eps6) then
            !          to check that the first nn bka vectors are not
            !          linearly dependent - this is a change from original code
            if (irguide .ne. 0) then
              do j = 1, 3
                rguide(j, loop_wann) = 0.0_dp
                do i = 1, 3
                  rguide(j, loop_wann) = rguide(j, loop_wann) + sinv(j, i) &
                                         *svec(i)/det
                enddo
              enddo
            endif
          endif
        endif

      enddo

    enddo

    !     obtain branch cut choice guided by rguid
    sheet = 0.0_dp
    do nkp = 1, num_kpts
      do nn = 1, kmesh_info%nntot
        do loop_wann = 1, num_wann
          ! sheet (loop_wann, nn, nkp) = 0.d0
          do j = 1, 3
            sheet(loop_wann, nn, nkp) = sheet(loop_wann, nn, nkp) &
                                        + kmesh_info%bk(j, nn, nkp)*rguide(j, loop_wann)
          enddo
          ! csheet (loop_wann, nn, nkp) = exp (ci * sheet (loop_wann, nn, nkp) )
        enddo
      enddo
    enddo
    csheet = exp(cmplx_i*sheet)

    ! now check that we picked the proper sheet for the log
    ! of m_matrix. criterion: q_n^{k,b}=Im(ln(M_nn^{k,b})) + b \cdot r_n are
    ! circa 0 for a good solution, circa multiples of 2 pi  for a bad one.
    ! I use the guiding center, instead of r_n, to understand which could be
    ! right sheet

    rnkb = 0.0_dp
    do nkp = 1, num_kpts
      do nn = 1, kmesh_info%nntot
        do m = 1, num_wann
          !           rnkb (m, nn, nkp) = 0.0_dp
          brn = 0.0_dp
          do ind = 1, 3
            brn = brn + kmesh_info%bk(ind, nn, nkp)*rguide(ind, m)
          enddo
          rnkb(m, nn, nkp) = rnkb(m, nn, nkp) + brn
        enddo
      enddo
    enddo
!    write ( stdout , * ) ' '
!    write ( stdout , * ) ' PHASES ARE SET USING THE GUIDING CENTERS'
!    write ( stdout , * ) ' '
!    do nkp = 1, num_kpts
!       do n = 1, num_wann
!          do nn = 1, nntot
!             pherr = aimag(log(csheet(n,nn,nkp)*m_matrix(n,n,nn,nkp))) &
!                  - sheet(n,nn,nkp)+rnkb(n,nn,nkp)-aimag(log(m_matrix(n,n,nn,nkp)))
!          enddo
!       enddo
!    enddo

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: phases', 2, stdout, seedname)

    return

  end subroutine wann_phases

  !==================================================================!
  subroutine wann_omega(csheet, sheet, rave, r2ave, rave2, wann_spread, num_wann, kmesh_info, &
                        num_kpts, verbose, wann_constrain, omega_invariant, counts, displs, &
                        ln_tmp_loc, m_matrix_loc, lambda_loc, first_pass, stdout, seedname, comm)
    !==================================================================!
    !                                                                  !
    !!   Calculate the Wannier Function spread                         !
    !                                                                  !
    ! Modified by Valerio Vitale for the SLWF+C method (PRB 90, 165125)!
    ! Jun 2018, based on previous work by Charles T. Johnson and       !
    ! Radu Miron at Implerial College London
    !===================================================================
    use w90_io, only: io_stopwatch
    use w90_comms, only: comms_allreduce, w90commtype, mpirank
    use w90_param_types, only: kmesh_info_type, print_output_type
    use wannier_param_types, only: wann_localise_type

    implicit none

    type(wann_localise_type), intent(in) :: wann_constrain
    real(kind=dp), intent(in) :: omega_invariant

    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: verbose

    complex(kind=dp), intent(in)  :: csheet(:, :, :)
    real(kind=dp), intent(in)  :: sheet(:, :, :)
    real(kind=dp), intent(out) :: rave(:, :)
    real(kind=dp), intent(out) :: r2ave(:)
    real(kind=dp), intent(out) :: rave2(:)
    type(localisation_vars), intent(out)  :: wann_spread

    ! from w90_parameters
    integer, intent(in) :: num_wann
!   integer, intent(in) :: nntot
!   real(kind=dp), intent(in) :: wb(:)
!   real(kind=dp), intent(in) :: bk(:, :, :)
    integer, intent(in) :: num_kpts
!   real(kind=dp), intent(in) :: omega_invariant
!   logical, intent(in) :: selective_loc
!   logical, intent(in) :: slwf_constrain
!   integer, intent(in) :: slwf_num
!   real(kind=dp), intent(in) :: ccentres_cart(:, :)
!   integer, intent(in) :: timing_level
    ! end of parameters
    integer, intent(in) :: counts(0:), displs(0:)
    real(kind=dp), intent(inout) :: ln_tmp_loc(:, :, :)
    complex(kind=dp), intent(in) :: m_matrix_loc(:, :, :, :)
    real(kind=dp), intent(in) :: lambda_loc
    logical, intent(inout) :: first_pass
    integer, intent(in) :: stdout
    type(w90commtype), intent(in) :: comm
    character(len=50), intent(in)  :: seedname

    !local variables
    real(kind=dp) :: summ, mnn2
    real(kind=dp) :: brn
    integer :: ind, nkp, nn, m, n, iw, nkp_loc
    integer :: my_node_id

    my_node_id = mpirank(comm)

    if (verbose%timing_level > 1 .and. verbose%iprint > 0) call io_stopwatch('wann: omega', 1, stdout, seedname)

    do nkp_loc = 1, counts(my_node_id)
      nkp = nkp_loc + displs(my_node_id)
      do nn = 1, kmesh_info%nntot
        do n = 1, num_wann
          ! Note that this ln_tmp is defined differently wrt the one in wann_domega
          ln_tmp_loc(n, nn, nkp_loc) = (aimag(log(csheet(n, nn, nkp) &
                                                  *m_matrix_loc(n, n, nn, nkp_loc))) - sheet(n, nn, nkp))
        end do
      end do
    end do

    rave = 0.0_dp
    do iw = 1, num_wann
      do ind = 1, 3
        do nkp_loc = 1, counts(my_node_id)
          nkp = nkp_loc + displs(my_node_id)
          do nn = 1, kmesh_info%nntot
            rave(ind, iw) = rave(ind, iw) + kmesh_info%wb(nn)*kmesh_info%bk(ind, nn, nkp) &
                            *ln_tmp_loc(iw, nn, nkp_loc)
          enddo
        enddo
      enddo
    enddo

    call comms_allreduce(rave(1, 1), num_wann*3, 'SUM', stdout, seedname, comm)

    rave = -rave/real(num_kpts, dp)

    rave2 = 0.0_dp
    do iw = 1, num_wann
      rave2(iw) = sum(rave(:, iw)*rave(:, iw))
    enddo

    ! aam: is this useful?
!~    rtot=0.0_dp
!~    do ind = 1, 3
!~       do loop_wann = 1, num_wann
!~          rtot (ind) = rtot (ind) + rave (ind, loop_wann)
!~       enddo
!~    enddo

    r2ave = 0.0_dp
    do iw = 1, num_wann
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        do nn = 1, kmesh_info%nntot
          mnn2 = real(m_matrix_loc(iw, iw, nn, nkp_loc) &
                      *conjg(m_matrix_loc(iw, iw, nn, nkp_loc)), kind=dp)
          r2ave(iw) = r2ave(iw) + kmesh_info%wb(nn)*(1.0_dp - mnn2 + ln_tmp_loc(iw, nn, nkp_loc)**2)
        enddo
      enddo
    enddo

    call comms_allreduce(r2ave(1), num_wann, 'SUM', stdout, seedname, comm)

    r2ave = r2ave/real(num_kpts, dp)

!~    wann_spread%om_1 = 0.0_dp
!~    do nkp = 1, num_kpts
!~       do nn = 1, nntot
!~          do loop_wann = 1, num_wann
!~             wann_spread%om_1 = wann_spread%om_1 + wb(nn) * &
!~                  ( 1.0_dp - m_matrix(loop_wann,loop_wann,nn,nkp) * &
!~                  conjg(m_matrix(loop_wann,loop_wann,nn,nkp)) )
!~          enddo
!~       enddo
!~    enddo
!~    wann_spread%om_1 = wann_spread%om_1 / real(num_kpts,dp)
!~
!~    wann_spread%om_2 = 0.0_dp
!~    do loop_wann = 1, num_wann
!~       sqim = 0.0_dp
!~       do nkp = 1, num_kpts
!~          do nn = 1, nntot
!~             sqim = sqim + wb(nn) * &
!~                  ( (aimag(log(csheet(loop_wann,nn,nkp) * &
!~                  m_matrix(loop_wann,loop_wann,nn,nkp))) - &
!~                  sheet(loop_wann,nn,nkp))**2 )
!~          enddo
!~       enddo
!~       sqim = sqim / real(num_kpts,dp)
!~       wann_spread%om_2 = wann_spread%om_2 + sqim
!~    enddo
!~
!~    wann_spread%om_3 = 0.0_dp
!~    do loop_wann = 1, num_wann
!~       bim = 0.0_dp
!~       do ind = 1, 3
!~          do nkp = 1, num_kpts
!~             do nn = 1, nntot
!~                bim(ind) = bim(ind) &
!~                     + wb(nn) * bk(ind,nn,nkp) &
!~                     * ( aimag(log(csheet(loop_wann,nn,nkp) &
!~                     * m_matrix(loop_wann,loop_wann,nn,nkp))) &
!~                     - sheet(loop_wann,nn,nkp) )
!~             enddo
!~          enddo
!~       enddo
!~       bim = bim/real(num_kpts,dp)
!~       bim2 = 0.0_dp
!~       do ind = 1, 3
!~          bim2 = bim2 + bim (ind) * bim (ind)
!~       enddo
!~       wann_spread%om_3 = wann_spread%om_3 - bim2
!~    enddo

    !jry: Either the above (om1,2,3) or the following is redundant
    !     keep it in the code base for testing

    if (wann_constrain%selective_loc) then
      wann_spread%om_iod = 0.0_dp
      do nkp_loc = 1, counts(my_node_id)
        do nn = 1, kmesh_info%nntot
          summ = 0.0_dp
          do n = 1, wann_constrain%slwf_num
            summ = summ &
                   + real(m_matrix_loc(n, n, nn, nkp_loc) &
                          *conjg(m_matrix_loc(n, n, nn, nkp_loc)), kind=dp)
            if (wann_constrain%slwf_constrain) then
              !! Centre constraint contribution. Zero if slwf_constrain=false
              summ = summ - lambda_loc*ln_tmp_loc(n, nn, nkp_loc)**2
            end if
          enddo
          wann_spread%om_iod = wann_spread%om_iod &
                               + kmesh_info%wb(nn)*(real(wann_constrain%slwf_num, dp) - summ)
        enddo
      enddo

      call comms_allreduce(wann_spread%om_iod, 1, 'SUM', stdout, seedname, comm)

      wann_spread%om_iod = wann_spread%om_iod/real(num_kpts, dp)

      wann_spread%om_d = 0.0_dp
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        do nn = 1, kmesh_info%nntot
          do n = 1, wann_constrain%slwf_num
            brn = sum(kmesh_info%bk(:, nn, nkp)*rave(:, n))
            wann_spread%om_d = wann_spread%om_d + (1.0_dp - lambda_loc)*kmesh_info%wb(nn) &
                               *(ln_tmp_loc(n, nn, nkp_loc) + brn)**2
          enddo
        enddo
      enddo

      call comms_allreduce(wann_spread%om_d, 1, 'SUM', stdout, seedname, comm)

      wann_spread%om_d = wann_spread%om_d/real(num_kpts, dp)

      wann_spread%om_nu = 0.0_dp
      !! Contribution from constrains on centres
      if (wann_constrain%slwf_constrain) then
        do nkp_loc = 1, counts(my_node_id)
          nkp = nkp_loc + displs(my_node_id)
          do nn = 1, kmesh_info%nntot
            do n = 1, wann_constrain%slwf_num
              wann_spread%om_nu = wann_spread%om_nu + 2.0_dp*kmesh_info%wb(nn)* &
                                  ln_tmp_loc(n, nn, nkp_loc)*lambda_loc* &
                                  sum(kmesh_info%bk(:, nn, nkp)*wann_constrain%ccentres_cart(n, :))
            enddo
          enddo
        enddo

        call comms_allreduce(wann_spread%om_nu, 1, 'SUM', stdout, seedname, comm)

        wann_spread%om_nu = wann_spread%om_nu/real(num_kpts, dp)

        do n = 1, wann_constrain%slwf_num
          wann_spread%om_nu = wann_spread%om_nu &
                              + lambda_loc*sum(wann_constrain%ccentres_cart(n, :)**2)
        end do

      end if

      wann_spread%om_tot = wann_spread%om_iod + wann_spread%om_d + wann_spread%om_nu
      !! wann_spread%om_c = wann_spread%om_iod + wann_spread%om_d + wann_spread%om_nu
    else
      if (first_pass) then
        wann_spread%om_i = 0.0_dp
        nkp = nkp_loc + displs(my_node_id)
        do nkp_loc = 1, counts(my_node_id)
          do nn = 1, kmesh_info%nntot
            summ = 0.0_dp
            do m = 1, num_wann
              do n = 1, num_wann
                summ = summ &
                       + real(m_matrix_loc(n, m, nn, nkp_loc) &
                              *conjg(m_matrix_loc(n, m, nn, nkp_loc)), kind=dp)
              enddo
            enddo
            wann_spread%om_i = wann_spread%om_i &
                               + kmesh_info%wb(nn)*(real(num_wann, dp) - summ)
          enddo
        enddo

        call comms_allreduce(wann_spread%om_i, 1, 'SUM', stdout, seedname, comm)

        wann_spread%om_i = wann_spread%om_i/real(num_kpts, dp)
        first_pass = .false.
      else
        wann_spread%om_i = omega_invariant
      endif

      wann_spread%om_od = 0.0_dp
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        do nn = 1, kmesh_info%nntot
          do m = 1, num_wann
            do n = 1, num_wann
              if (m .ne. n) wann_spread%om_od = wann_spread%om_od &
                                                + kmesh_info%wb(nn)*real(m_matrix_loc(n, m, nn, nkp_loc) &
                                                                         *conjg(m_matrix_loc(n, m, nn, nkp_loc)), kind=dp)
            enddo
          enddo
        enddo
      enddo

      call comms_allreduce(wann_spread%om_od, 1, 'SUM', stdout, seedname, comm)

      wann_spread%om_od = wann_spread%om_od/real(num_kpts, dp)

      wann_spread%om_d = 0.0_dp
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        do nn = 1, kmesh_info%nntot
          do n = 1, num_wann
            brn = sum(kmesh_info%bk(:, nn, nkp)*rave(:, n))
            wann_spread%om_d = wann_spread%om_d + kmesh_info%wb(nn) &
                               *(ln_tmp_loc(n, nn, nkp_loc) + brn)**2
          enddo
        enddo
      enddo

      call comms_allreduce(wann_spread%om_d, 1, 'SUM', stdout, seedname, comm)

      wann_spread%om_d = wann_spread%om_d/real(num_kpts, dp)

      wann_spread%om_tot = wann_spread%om_i + wann_spread%om_d + wann_spread%om_od
    end if

    if (verbose%timing_level > 1 .and. verbose%iprint > 0) call io_stopwatch('wann: omega', 2, stdout, seedname)

    return

  end subroutine wann_omega

  !==================================================================!
  subroutine wann_domega(csheet, sheet, rave, num_wann, kmesh_info, num_kpts, wann_constrain, &
                         lsitesymmetry, counts, displs, ln_tmp_loc, m_matrix_loc, rnkb_loc, &
                         cdodq_loc, lambda_loc, timing_level, stdout, seedname, sym, comm, iprint, &
                         cdodq)
    !==================================================================!
    !                                                                  !
    !   Calculate the Gradient of the Wannier Function spread          !
    !                                                                  !
    ! Modified by Valerio Vitale for the SLWF+C method (PRB 90, 165125)!
    ! Jun 2018, based on previous work by Charles T. Johnson and       !
    ! Radu Miron at Imperial College London
    !===================================================================
    use w90_constants, only: cmplx_0
    use w90_io, only: io_stopwatch, io_error
    use w90_sitesym, only: sitesym_symmetrize_gradient, sitesym_data !RS:
    use w90_comms, only: comms_gatherv, comms_bcast, comms_allreduce, &
      w90commtype, mpirank
    use w90_param_types, only: kmesh_info_type
    use wannier_param_types, only: wann_localise_type

    implicit none

    type(kmesh_info_type), intent(in) :: kmesh_info
    type(wann_localise_type), intent(inout) :: wann_constrain

    complex(kind=dp), intent(in)  :: csheet(:, :, :)
    real(kind=dp), intent(in)  :: sheet(:, :, :)
    real(kind=dp), intent(out) :: rave(:, :)
    ! as we work on the local cdodq, returning the full cdodq array is now
    ! made optional
    complex(kind=dp), intent(out), optional :: cdodq(:, :, :)

    type(sitesym_data), intent(in) :: sym
    ! from w90_parameters
    integer, intent(in) :: num_wann
!   integer, intent(in) :: nntot
!   real(kind=dp), intent(in) :: wb(:)
!   real(kind=dp), intent(in) :: bk(:, :, :)
    integer, intent(in) :: num_kpts
!   logical, intent(in) :: selective_loc
!   logical, intent(in) :: slwf_constrain
!    integer, intent(in) :: slwf_num
!   real(kind=dp), intent(in) :: ccentres_cart(:, :)
    logical, intent(in) :: lsitesymmetry
    integer, intent(in) :: timing_level, iprint
    character(len=50), intent(in)  :: seedname

    ! end parameters
    integer, intent(in) :: counts(0:), displs(0:)
    real(kind=dp), intent(inout) :: ln_tmp_loc(:, :, :)
    complex(kind=dp), intent(in) :: m_matrix_loc(:, :, :, :)
    real(kind=dp), intent(inout) :: rnkb_loc(:, :, :)
    complex(kind=dp), intent(out) :: cdodq_loc(:, :, :)
    real(kind=dp), intent(in) :: lambda_loc
    integer, intent(in) :: stdout
    type(w90commtype), intent(in) :: comm

    ! local
    complex(kind=dp), allocatable  :: cr(:, :)
    complex(kind=dp), allocatable  :: crt(:, :)
    real(kind=dp), allocatable :: r0kb(:, :, :)
    integer :: iw, ind, nkp, nn, m, n, ierr, nkp_loc
    complex(kind=dp) :: mnn
    integer :: my_node_id

    my_node_id = mpirank(comm)

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: domega', 1, stdout, seedname)

    allocate (cr(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cr in wann_domega', stdout, seedname)
    allocate (crt(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating crt in wann_domega', stdout, seedname)
    if (wann_constrain%selective_loc .and. wann_constrain%slwf_constrain) then
      allocate (r0kb(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating r0kb in wann_domega', stdout, seedname)
    end if

    do nkp_loc = 1, counts(my_node_id)
      nkp = nkp_loc + displs(my_node_id)
      do nn = 1, kmesh_info%nntot
        do n = 1, num_wann
          ! Note that this ln_tmp is defined differently wrt the one in wann_omega
          ln_tmp_loc(n, nn, nkp_loc) = kmesh_info%wb(nn)*(aimag(log(csheet(n, nn, nkp) &
                                                                    *m_matrix_loc(n, n, nn, nkp_loc))) - sheet(n, nn, nkp))
        end do
      end do
    end do

    ! recalculate rave
    rave = 0.0_dp
    do iw = 1, num_wann
      do ind = 1, 3
        do nkp_loc = 1, counts(my_node_id)
          nkp = nkp_loc + displs(my_node_id)
          do nn = 1, kmesh_info%nntot
            rave(ind, iw) = rave(ind, iw) + kmesh_info%bk(ind, nn, nkp) &
                            *ln_tmp_loc(iw, nn, nkp_loc)
          enddo
        enddo
      enddo
    enddo
    rave = -rave/real(num_kpts, dp)

    call comms_allreduce(rave(1, 1), num_wann*3, 'SUM', stdout, seedname, comm)

    ! b.r_0n are calculated
    if (wann_constrain%selective_loc .and. wann_constrain%slwf_constrain) then
      r0kb = 0.0_dp
      do nkp_loc = 1, counts(my_node_id)
        nkp = nkp_loc + displs(my_node_id)
        do nn = 1, kmesh_info%nntot
          do n = 1, num_wann
            r0kb(n, nn, nkp_loc) = sum(kmesh_info%bk(:, nn, nkp) &
                                       *wann_constrain%ccentres_cart(n, :))
          enddo
        enddo
      enddo
    end if

    rnkb_loc = 0.0_dp
    do nkp_loc = 1, counts(my_node_id)
      nkp = nkp_loc + displs(my_node_id)
      do nn = 1, kmesh_info%nntot
        do n = 1, num_wann
          rnkb_loc(n, nn, nkp_loc) = sum(kmesh_info%bk(:, nn, nkp)*rave(:, n))
        enddo
      enddo
    enddo

    ! cd0dq(m,n,nkp) is calculated
    cdodq_loc = cmplx_0
    cr = cmplx_0
    crt = cmplx_0
    do nkp_loc = 1, counts(my_node_id)
      nkp = nkp_loc + displs(my_node_id)
      do nn = 1, kmesh_info%nntot
        do n = 1, num_wann ! R^{k,b} and R~^{k,b} have columns of zeroes for the non-objective Wannier functions
          mnn = m_matrix_loc(n, n, nn, nkp_loc)
          crt(:, n) = m_matrix_loc(:, n, nn, nkp_loc)/mnn
          cr(:, n) = m_matrix_loc(:, n, nn, nkp_loc)*conjg(mnn)
        enddo
        if (wann_constrain%selective_loc) then
          do n = 1, num_wann
            do m = 1, num_wann
              if (m <= wann_constrain%slwf_num) then
                if (n <= wann_constrain%slwf_num) then
                  ! A[R^{k,b}]=(R-Rdag)/2
                  cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                             + kmesh_info%wb(nn)*0.5_dp*(cr(m, n) - conjg(cr(n, m)))
                  cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                             - (crt(m, n)*ln_tmp_loc(n, nn, nkp_loc) &
                                                + conjg(crt(n, m)*ln_tmp_loc(m, nn, nkp_loc))) &
                                             *cmplx(0.0_dp, -0.5_dp, kind=dp)
                  cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                             - (crt(m, n)*rnkb_loc(n, nn, nkp_loc) &
                                                + conjg(crt(n, m)*rnkb_loc(m, nn, nkp_loc))) &
                                             *cmplx(0.0_dp, -0.5_dp, kind=dp)
                  if (wann_constrain%slwf_constrain) then
                    cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) + lambda_loc &
                                               *(crt(m, n)*ln_tmp_loc(n, nn, nkp_loc) &
                                                 + conjg(crt(n, m)*ln_tmp_loc(m, nn, nkp_loc))) &
                                               *cmplx(0.0_dp, -0.5_dp, kind=dp)
                    cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                               + kmesh_info%wb(nn)*lambda_loc &
                                               *(crt(m, n)*rnkb_loc(n, nn, nkp_loc) &
                                                 + conjg(crt(n, m)*rnkb_loc(m, nn, nkp_loc))) &
                                               *cmplx(0.0_dp, -0.5_dp, kind=dp)
                    cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) - lambda_loc &
                                               *(crt(m, n)*ln_tmp_loc(n, nn, nkp_loc) &
                                                 + conjg(crt(n, m))*ln_tmp_loc(m, nn, nkp_loc)) &
                                               *cmplx(0.0_dp, -0.5_dp, kind=dp)
                    cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                               - kmesh_info%wb(nn)*lambda_loc &
                                               *(r0kb(n, nn, nkp_loc)*crt(m, n) &
                                                 + r0kb(m, nn, nkp_loc)*conjg(crt(n, m))) &
                                               *cmplx(0.0_dp, -0.5_dp, kind=dp)
                  end if
                else
                  cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) - kmesh_info%wb(nn) &
                                             *0.5_dp*conjg(cr(n, m)) &
                                             - conjg(crt(n, m)*(ln_tmp_loc(m, nn, nkp_loc) &
                                                                + kmesh_info%wb(nn)*rnkb_loc(m, nn, nkp_loc))) &
                                             *cmplx(0.0_dp, -0.5_dp, kind=dp)
                  if (wann_constrain%slwf_constrain) then
                    cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) + lambda_loc &
                                               *conjg(crt(n, m)*(ln_tmp_loc(m, nn, nkp_loc) &
                                                                 + kmesh_info%wb(nn)*rnkb_loc(m, nn, nkp_loc))) &
                                               *cmplx(0.0_dp, -0.5_dp, kind=dp) &
                                               - lambda_loc*(conjg(crt(n, m)) &
                                                             *ln_tmp_loc(m, nn, nkp_loc)) &
                                               *cmplx(0.0_dp, -0.5_dp, kind=dp)
                    cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                               - kmesh_info%wb(nn)*lambda_loc &
                                               *r0kb(m, nn, nkp_loc)*conjg(crt(n, m)) &
                                               *cmplx(0.0_dp, -0.5_dp, kind=dp)
                  end if
                end if
              else if (n <= wann_constrain%slwf_num) then
                cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                           + kmesh_info%wb(nn)*cr(m, n)*0.5_dp &
                                           - crt(m, n)*(ln_tmp_loc(n, nn, nkp_loc) &
                                                        + kmesh_info%wb(nn)*rnkb_loc(n, nn, nkp_loc)) &
                                           *cmplx(0.0_dp, -0.5_dp, kind=dp)
                if (wann_constrain%slwf_constrain) then
                  cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) + lambda_loc &
                                             *crt(m, n)*(ln_tmp_loc(n, nn, nkp_loc) &
                                                         + kmesh_info%wb(nn)*rnkb_loc(n, nn, nkp_loc)) &
                                             *cmplx(0.0_dp, -0.5_dp, kind=dp) &
                                             - lambda_loc*crt(m, n)*ln_tmp_loc(n, nn, nkp_loc) &
                                             *cmplx(0.0_dp, -0.5_dp, kind=dp)
                  cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) - kmesh_info%wb(nn) &
                                             *lambda_loc &
                                             *r0kb(n, nn, nkp_loc)*crt(m, n) &
                                             *cmplx(0.0_dp, -0.5_dp, kind=dp)
                end if
              else
                cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc)
              end if
            enddo
          enddo
        else
          do n = 1, num_wann
            do m = 1, num_wann
              ! A[R^{k,b}]=(R-Rdag)/2
              cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) &
                                         + kmesh_info%wb(nn)*0.5_dp &
                                         *(cr(m, n) - conjg(cr(n, m)))
              ! -S[T^{k,b}]=-(T+Tdag)/2i ; T_mn = Rt_mn q_n
              cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) - &
                                         (crt(m, n)*ln_tmp_loc(n, nn, nkp_loc) &
                                          + conjg(crt(n, m)*ln_tmp_loc(m, nn, nkp_loc))) &
                                         *cmplx(0.0_dp, -0.5_dp, kind=dp)
              cdodq_loc(m, n, nkp_loc) = cdodq_loc(m, n, nkp_loc) - kmesh_info%wb(nn) &
                                         *(crt(m, n)*rnkb_loc(n, nn, nkp_loc) &
                                           + conjg(crt(n, m)*rnkb_loc(m, nn, nkp_loc))) &
                                         *cmplx(0.0_dp, -0.5_dp, kind=dp)
            enddo
          enddo
        end if
      enddo
    enddo
    cdodq_loc = cdodq_loc/real(num_kpts, dp)*4.0_dp

    if (present(cdodq)) then
      ! each process communicates its result to other processes
      call comms_gatherv(cdodq_loc, num_wann*num_wann*counts(my_node_id), &
                         cdodq, num_wann*num_wann*counts, num_wann*num_wann*displs, stdout, &
                         seedname, comm)
      call comms_bcast(cdodq(1, 1, 1), num_wann*num_wann*num_kpts, stdout, seedname, comm)
      if (lsitesymmetry) then
        call sitesym_symmetrize_gradient(sym, cdodq, 1, num_kpts, num_wann) !RS:
        cdodq_loc(:, :, 1:counts(my_node_id)) = cdodq(:, :, displs(my_node_id) &
                                                      + 1:displs(my_node_id) + counts(my_node_id))
      endif
    end if

    deallocate (cr, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cr in wann_domega', stdout, seedname)
    deallocate (crt, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating crt in wann_domega', stdout, seedname)

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: domega', 2, stdout, seedname)

    return

  end subroutine wann_domega

  !==================================================================!
  subroutine wann_spread_copy(orig, copy)
    !==================================================================!
    !                                                                  !
    !==================================================================!

    implicit none

    type(localisation_vars), intent(in)  :: orig
    type(localisation_vars), intent(out) :: copy

    copy%om_i = orig%om_i
    copy%om_d = orig%om_d
    copy%om_od = orig%om_od
    copy%om_tot = orig%om_tot
    copy%om_iod = orig%om_iod
    copy%om_nu = orig%om_nu
    !! copy%om_c  =  orig%om_c
!~    copy%om_1   =  orig%om_1
!~    copy%om_2   =  orig%om_2
!~    copy%om_3   =  orig%om_3

    return

  end subroutine wann_spread_copy

  !==================================================================!
  subroutine wann_calc_projection(num_bands, num_wann, num_kpts, u_matrix_opt, eigval, lwindow, &
                                  timing_level, iprint, stdout, seedname)
    !==================================================================!
    !                                                                  !
    ! Calculates and writes the projection of each Wannier function    !
    ! on the original bands within the outer window.                   !
    !                                                                  !
    !==================================================================!

    use w90_io, only: io_stopwatch
    use w90_comms, only: w90commtype

    implicit none

    ! These were in the parameter module
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann
    integer, intent(in) :: num_kpts
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    real(kind=dp), intent(in) :: eigval(:, :)
    logical, intent(in) :: lwindow(:, :)
    integer, intent(in) :: timing_level, iprint
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    ! end of vars from parameter module

    integer :: nw, nb, nkp, counter
    real(kind=dp) :: summ

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: calc_projection', 1, stdout, seedname)

    if (iprint > 0) then
      write (stdout, '(/1x,a78)') repeat('-', 78)
      write (stdout, '(1x,9x,a)') &
        'Projection of Bands in Outer Window on all Wannier Functions'
      write (stdout, '(1x,8x,62a)') repeat('-', 62)
      write (stdout, '(1x,16x,a)') '   Kpt  Band      Eigval      |Projection|^2'
      write (stdout, '(1x,16x,a47)') repeat('-', 47)
    endif

    do nkp = 1, num_kpts
      counter = 0
      do nb = 1, num_bands
        if (lwindow(nb, nkp)) then
          counter = counter + 1
          summ = 0.0_dp
          do nw = 1, num_wann
            summ = summ + abs(u_matrix_opt(counter, nw, nkp))**2
          enddo
          if (iprint > 0) write (stdout, '(1x,16x,i5,1x,i5,1x,f14.6,2x,f14.8)') &
            nkp, nb, eigval(nb, nkp), summ
        endif
      enddo
    enddo
    if (iprint > 0) write (stdout, '(1x,a78/)') repeat('-', 78)

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: calc_projection', 2, stdout, seedname)

    return

  end subroutine wann_calc_projection

  !=====================================!
  subroutine wann_write_xyz(translate_home_cell, num_wann, wannier_centres, real_lattice, &
                            recip_lattice, atoms, verbose, stdout, seedname)
    !=====================================!
    !                                     !
    ! Write xyz file with Wannier centres !
    !                                     !
    !=====================================!

!   use w90_io, only: seedname, io_file_unit, io_date
    use w90_io, only: io_file_unit, io_date
    use w90_utility, only: utility_translate_home
    use w90_param_types, only: atom_data_type, print_output_type

    implicit none

    type(atom_data_type), intent(in) :: atoms
    type(print_output_type), intent(in) :: verbose

    ! from w90_parameters
    logical, intent(in) :: translate_home_cell
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: wannier_centres(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
!   integer, intent(in) :: num_atoms
!   character(len=2), intent(in) :: atoms_symbol(:)
!   real(kind=dp), intent(in) :: atoms_pos_cart(:, :, :)
!   integer, intent(in) :: num_species
!   integer, intent(in) :: atoms_species_num(:)
!   real(kind=dp), intent(in) :: lenconfac
!   integer, intent(in) :: iprint
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    ! end parameters

    integer          :: iw, ind, xyz_unit, nsp, nat
    character(len=9) :: cdate, ctime
    real(kind=dp)    :: wc(3, num_wann)

    wc = wannier_centres

    if (translate_home_cell) then
      do iw = 1, num_wann
        call utility_translate_home(wc(:, iw), real_lattice, recip_lattice)
      enddo
    endif

    if (verbose%iprint > 2) then
      write (stdout, '(1x,a)') 'Final centres (translated to home cell for writing xyz file)'
      do iw = 1, num_wann
        write (stdout, 888) iw, (wc(ind, iw)*verbose%lenconfac, ind=1, 3)
      end do
      write (stdout, '(1x,a78)') repeat('-', 78)
      write (stdout, *)
    endif

    xyz_unit = io_file_unit()
    open (xyz_unit, file=trim(seedname)//'_centres.xyz', form='formatted')
    write (xyz_unit, '(i6)') num_wann + atoms%num_atoms
    call io_date(cdate, ctime)
    write (xyz_unit, *) 'Wannier centres, written by Wannier90 on'//cdate//' at '//ctime
    do iw = 1, num_wann
      write (xyz_unit, '("X",6x,3(f14.8,3x))') (wc(ind, iw), ind=1, 3)
    end do
    do nsp = 1, atoms%num_species
      do nat = 1, atoms%species_num(nsp)
        write (xyz_unit, '(a2,5x,3(f14.8,3x))') atoms%symbol(nsp), atoms%pos_cart(:, nat, nsp)
      end do
    end do
    close (xyz_unit)

    write (stdout, '(/a)') ' Wannier centres written to file '//trim(seedname)//'_centres.xyz'

    return

888 format(2x, 'WF centre', i5, 2x, '(', f10.6, ',', f10.6, ',', f10.6, ' )')

  end subroutine wann_write_xyz

  !=================================================================!
  subroutine wann_write_vdw_data(num_wann, wann_data, real_lattice, recip_lattice, u_matrix, &
                                 u_matrix_opt, have_disentangled, system, stdout, seedname)
    !=================================================================!
    !                                                                 !
    ! Write a file with Wannier centres, spreads and occupations for  !
    ! post-processing computation of vdW C6 coeffients.               !
    !                                                                 !
    ! Based on code written by Lampros Andrinopoulos.                 !
    !=================================================================!

!   use w90_io, only: seedname, io_file_unit, io_date, io_error
    use w90_io, only: io_file_unit, io_date, io_error
    use w90_utility, only: utility_translate_home
    use w90_constants, only: cmplx_0
    use w90_param_types, only: wannier_data_type, w90_system_type
!~    use w90_disentangle, only : ndimfroz

    implicit none

    type(wannier_data_type), intent(in) :: wann_data
    type(w90_system_type), intent(in) :: system

    ! from w90_parameters
    !logical, intent(in) :: translate_home_cell
    integer, intent(in) :: num_wann
!   real(kind=dp), intent(in) :: wannier_centres(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    !character(len=2), intent(in) :: atoms_symbol(:)
    !real(kind=dp), intent(in) :: atoms_pos_cart(:, :, :)
    !integer, intent(in) :: num_species
    !integer, intent(in) :: atoms_species_num(:)
!   real(kind=dp), intent(in) :: wannier_spreads(:)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    !integer, intent(in) :: num_kpts
!   logical, intent(in) :: have_disentangled
!   integer, intent(in) :: num_valence_bands
!   integer, intent(in) :: num_elec_per_state
    ! end parameters
    integer, intent(in) :: stdout
    logical, intent(in) :: have_disentangled
    character(len=50), intent(in)  :: seedname

    integer          :: iw, vdw_unit, r, s, k, m, ierr, ndim
    real(kind=dp)    :: wc(3, num_wann)
    real(kind=dp)    :: ws(num_wann)
    complex(kind=dp), allocatable :: f_w(:, :), v_matrix(:, :) !f_w2(:,:)

    wc = wann_data%centres
    ws = wann_data%spreads

    ! translate Wannier centres to the home unit cell
    do iw = 1, num_wann
      call utility_translate_home(wc(:, iw), real_lattice, recip_lattice)
    enddo

    allocate (f_w(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating f_w in wann_write_vdw_data', stdout, seedname)

!~    ! aam: remove f_w2 at end
!~    allocate(f_w2(num_wann, num_wann),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating f_w2 in wann_write_vdw_data')

    if (have_disentangled) then

      ! dimension of occupied subspace
      if (system%num_valence_bands .le. 0) &
        call io_error('Please set num_valence_bands in seedname.win', stdout, seedname)
      ndim = system%num_valence_bands

      allocate (v_matrix(ndim, num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating V_matrix in wann_write_vdw_data', stdout, seedname)

      ! aam: initialise
      f_w(:, :) = cmplx_0
      v_matrix(:, :) = cmplx_0
!~       f_w2(:,:) = cmplx_0

      ! aam: IN THE END ONLY NEED DIAGONAL PART, SO COULD SIMPLIFY...
      ! aam: calculate V = U_opt . U
      do s = 1, num_wann
        do k = 1, ndim
          do m = 1, num_wann
            v_matrix(k, s) = v_matrix(k, s) + u_matrix_opt(k, m, 1)*u_matrix(m, s, 1)
          enddo
        enddo
      enddo

      ! aam: calculate f = V^dagger . V
      do r = 1, num_wann
        do s = 1, num_wann
          do k = 1, ndim
            f_w(r, s) = f_w(r, s) + v_matrix(k, s)*conjg(v_matrix(k, r))
          enddo
        enddo
      enddo

!~       ! original formulation
!~       do r=1,num_wann
!~          do s=1,num_wann
!~             do nkp=1,num_kpts
!~                do k=1,ndimfroz(nkp)
!~                   do m=1,num_wann
!~                      do l=1,num_wann
!~                         f_w2(r,s) = f_w2(r,s) + &
!~                              u_matrix_opt(k,m,nkp) * u_matrix(m,s,nkp) * &
!~                              conjg(u_matrix_opt(k,l,nkp)) * conjg(u_matrix(l,r,nkp))
!~                      end do
!~                   end do
!~                end do
!~             end do
!~          end do
!~       end do

!~       ! test equivalence
!~       do r=1,num_wann
!~          do s=1,num_wann
!~             if (abs(real(f_w(r,s),dp)-real(f_w2(r,s),dp)).gt.eps6) then
!~                write(*,'(i6,i6,f16.10,f16.10)') r,s,real(f_w(r,s),dp),real(f_w2(r,s),dp)
!~             endif
!~             if (abs(aimag(f_w(r,s))-aimag(f_w2(r,s))).gt.eps6) then
!~                write(*,'(a,i6,i6,f16.10,f16.10)') 'Im: ',r,s,aimag(f_w(r,s)),aimag(f_w2(r,s))
!~             endif
!~          enddo
!~       enddo
!~       write(*,*) ' done vdw '

    else
      ! for valence only, all occupancies are unity
      f_w(:, :) = 1.0_dp
    endif

    ! aam: write the seedname.vdw file directly here
    vdw_unit = io_file_unit()
    open (unit=vdw_unit, file=trim(seedname)//'.vdw', action='write')
    if (have_disentangled) then
      write (vdw_unit, '(a)') 'disentangle T'
    else
      write (vdw_unit, '(a)') 'disentangle F'
    endif
    write (vdw_unit, '(a)') 'amalgamate F'
    write (vdw_unit, '(a,i3)') 'degeneracy', system%num_elec_per_state
    write (vdw_unit, '(a)') 'num_frag 2'
    write (vdw_unit, '(a)') 'num_wann'
    write (vdw_unit, '(i3,1x,i3)') num_wann/2, num_wann/2
    write (vdw_unit, '(a)') 'tol_occ 0.9'
    write (vdw_unit, '(a)') 'pxyz'
    write (vdw_unit, '(a)') 'F F F'
    write (vdw_unit, '(a)') 'F F F'
    write (vdw_unit, '(a)') 'tol_dist 0.05'
    write (vdw_unit, '(a)') 'centres_spreads_occ'
    write (vdw_unit, '(a)') 'ang'
    do iw = 1, num_wann
      write (vdw_unit, '(4(f13.10,1x),1x,f11.8)') wc(1:3, iw), ws(iw), real(f_w(iw, iw))
    end do
    close (vdw_unit)

    write (stdout, '(/a/)') ' vdW data written to file '//trim(seedname)//'.vdw'

    if (have_disentangled) then
      deallocate (v_matrix, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating v_matrix in wann_write_vdw_data', stdout, seedname)
    endif

!~    deallocate(f_w2,stat=ierr)
!~    if (ierr/=0) call io_error('Error in deallocating f_w2 in wann_write_vdw_data')
    deallocate (f_w, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating f_w in wann_write_vdw_data', stdout, seedname)

    return

  end subroutine wann_write_vdw_data

  !========================================!
  subroutine wann_check_unitarity(num_kpts, num_wann, u_matrix, timing_level, iprint, stdout, &
                                  seedname)
    !========================================!

    use w90_constants, only: dp, cmplx_1, cmplx_0, eps5
    use w90_io, only: io_stopwatch, io_error
    use w90_comms, only: w90commtype

    implicit none

    ! passed variables
    integer, intent(in) :: num_kpts, num_wann, timing_level, iprint, stdout
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)
    character(len=50), intent(in)  :: seedname

    ! local variables
    integer :: nkp, i, j, m
    complex(kind=dp) :: ctmp1, ctmp2

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: check_unitarity', 1, stdout, seedname)

    do nkp = 1, num_kpts
      do i = 1, num_wann
        do j = 1, num_wann
          ctmp1 = cmplx_0
          ctmp2 = cmplx_0
          do m = 1, num_wann
            ctmp1 = ctmp1 + u_matrix(i, m, nkp)*conjg(u_matrix(j, m, nkp))
            ctmp2 = ctmp2 + u_matrix(m, j, nkp)*conjg(u_matrix(m, i, nkp))
          enddo
          if ((i .eq. j) .and. (abs(ctmp1 - cmplx_1) .gt. eps5)) &
            then
            if (iprint > 0) write (stdout, *) ' ERROR: unitariety of final U', nkp, i, j, &
              ctmp1
            call io_error('wann_check_unitarity: error 1', stdout, seedname)
          endif
          if ((i .eq. j) .and. (abs(ctmp2 - cmplx_1) .gt. eps5)) &
            then
            if (iprint > 0) write (stdout, *) ' ERROR: unitariety of final U', nkp, i, j, &
              ctmp2
            call io_error('wann_check_unitarity: error 2', stdout, seedname)
          endif
          if ((i .ne. j) .and. (abs(ctmp1) .gt. eps5)) then
            if (iprint > 0) write (stdout, *) ' ERROR: unitariety of final U', nkp, i, j, &
              ctmp1
            call io_error('wann_check_unitarity: error 3', stdout, seedname)
          endif
          if ((i .ne. j) .and. (abs(ctmp2) .gt. eps5)) then
            if (iprint > 0) write (stdout, *) ' ERROR: unitariety of final U', nkp, i, j, &
              ctmp2
            call io_error('wann_check_unitarity: error 4', stdout, seedname)
          endif
        enddo
      enddo
    enddo

    if (timing_level > 1 .and. iprint > 0) call io_stopwatch('wann: check_unitarity', 2, stdout, seedname)

    return

  end subroutine wann_check_unitarity

  !========================================!
  subroutine wann_write_r2mn(num_kpts, num_wann, kmesh_info, m_matrix, stdout, seedname)
    !========================================!
    !                                        !
    ! Write seedname.r2mn file               !
    !                                        !
    !========================================!

    use w90_constants, only: dp
!   use w90_io, only: seedname, io_file_unit, io_error
    use w90_io, only: io_file_unit, io_error
    use w90_param_types, only: kmesh_info_type

    implicit none

    type(kmesh_info_type), intent(in) :: kmesh_info

    ! from parameters
!   integer, intent(in) :: num_kpts, num_wann, nntot
    integer, intent(in) :: num_kpts, num_wann
!   real(kind=dp), intent(in) :: wb(:)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: r2mnunit, nw1, nw2, nkp, nn
    real(kind=dp) :: r2ave_mn, delta

    ! note that here I use formulas analogue to Eq. 23, and not to the
    ! shift-invariant Eq. 32 .
    r2mnunit = io_file_unit()
    open (r2mnunit, file=trim(seedname)//'.r2mn', form='formatted', err=158)
    do nw1 = 1, num_wann
      do nw2 = 1, num_wann
        r2ave_mn = 0.0_dp
        delta = 0.0_dp
        if (nw1 .eq. nw2) delta = 1.0_dp
        do nkp = 1, num_kpts
          do nn = 1, kmesh_info%nntot
            r2ave_mn = r2ave_mn + kmesh_info%wb(nn)* &
                       ! [GP-begin, Apr13, 2012: corrected sign inside "real"]
                       (2.0_dp*delta - real(m_matrix(nw1, nw2, nn, nkp) + &
                                            conjg(m_matrix(nw2, nw1, nn, nkp)), kind=dp))
            ! [GP-end]
          enddo
        enddo
        r2ave_mn = r2ave_mn/real(num_kpts, dp)
        write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn
      enddo
    enddo
    close (r2mnunit)

    return

158 call io_error('Error opening file '//trim(seedname)//'.r2mn in wann_write_r2mn', stdout, seedname)

  end subroutine wann_write_r2mn

  !========================================!
  subroutine wann_svd_omega_i(num_wann, num_kpts, kmesh_info, m_matrix, verbose, stdout, seedname)
    !========================================!

    use w90_constants, only: dp, cmplx_0
    use w90_io, only: io_stopwatch, io_error
    use w90_comms, only: w90commtype
    use w90_param_types, only: kmesh_info_type, print_output_type

    implicit none

    type(print_output_type), intent(in) :: verbose
    type(kmesh_info_type), intent(in) :: kmesh_info

    ! from w90_parameters
!   integer, intent(in) :: num_wann, num_kpts, nntot
    integer, intent(in) :: num_wann, num_kpts
!   real(kind=dp), intent(in) :: wb(:)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
!   real(kind=dp), intent(in) :: lenconfac
!   character(len=*), intent(in) :: length_unit
!   integer, intent(in) :: timing_level
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    complex(kind=dp), allocatable  :: cv1(:, :), cv2(:, :)
    complex(kind=dp), allocatable  :: cw1(:), cw2(:)
    complex(kind=dp), allocatable  :: cpad1(:)
    real(kind=dp), allocatable  :: singvd(:)

    integer :: ierr, info
    integer :: nkp, nn, nb, na, ind
    real(kind=dp) :: omt1, omt2, omt3

    if (verbose%timing_level > 1 .and. verbose%iprint > 0) call io_stopwatch('wann: svd_omega_i', 1, stdout, seedname)

    allocate (cw1(10*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cw1 in wann_svd_omega_i', stdout, seedname)
    allocate (cw2(10*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cw2 in wann_svd_omega_i', stdout, seedname)
    allocate (cv1(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cv1 in wann_svd_omega_i', stdout, seedname)
    allocate (cv2(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cv2 in wann_svd_omega_i', stdout, seedname)
    allocate (singvd(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating singvd in wann_svd_omega_i', stdout, seedname)
    allocate (cpad1(num_wann*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cpad1 in wann_svd_omega_i', stdout, seedname)

    cw1 = cmplx_0; cw2 = cmplx_0; cv1 = cmplx_0; cv2 = cmplx_0; cpad1 = cmplx_0
    singvd = 0.0_dp

    ! singular value decomposition
    omt1 = 0.0_dp; omt2 = 0.0_dp; omt3 = 0.0_dp
    do nkp = 1, num_kpts
      do nn = 1, kmesh_info%nntot
        ind = 1
        do nb = 1, num_wann
          do na = 1, num_wann
            cpad1(ind) = m_matrix(na, nb, nn, nkp)
            ind = ind + 1
          enddo
        enddo
        call zgesvd('A', 'A', num_wann, num_wann, cpad1, num_wann, singvd, cv1, &
                    num_wann, cv2, num_wann, cw1, 10*num_wann, cw2, info)
        if (info .ne. 0) then
          call io_error('ERROR: Singular value decomp. zgesvd failed', stdout, seedname)
        endif

        do nb = 1, num_wann
          omt1 = omt1 + kmesh_info%wb(nn)*(1.0_dp - singvd(nb)**2)
          omt2 = omt2 - kmesh_info%wb(nn)*(2.0_dp*log(singvd(nb)))
          omt3 = omt3 + kmesh_info%wb(nn)*(acos(singvd(nb))**2)
        enddo
      enddo
    enddo
    omt1 = omt1/real(num_kpts, dp)
    omt2 = omt2/real(num_kpts, dp)
    omt3 = omt3/real(num_kpts, dp)
    if (verbose%iprint > 0) then
      write (stdout, *) ' '
      write (stdout, '(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ', &
        omt1*verbose%lenconfac**2, '('//trim(verbose%length_unit)//'^2)'
      write (stdout, '(2x,a,f15.9,1x,a)') '                 -2log s = ', &
        omt2*verbose%lenconfac**2, '('//trim(verbose%length_unit)//'^2)'
      write (stdout, '(2x,a,f15.9,1x,a)') '                  acos^2 = ', &
        omt3*verbose%lenconfac**2, '('//trim(verbose%length_unit)//'^2)'
    endif

    deallocate (cpad1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cpad1 in wann_svd_omega_i', stdout, seedname)
    deallocate (singvd, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating singvd in wann_svd_omega_i', stdout, seedname)
    deallocate (cv2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cv2 in wann_svd_omega_i', stdout, seedname)
    deallocate (cv1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cv1 in wann_svd_omega_i', stdout, seedname)
    deallocate (cw2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cw2 in wann_svd_omega_i', stdout, seedname)
    deallocate (cw1, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cw1 in wann_svd_omega_i', stdout, seedname)

    if (verbose%timing_level > 1 .and. verbose%iprint > 0) call io_stopwatch('wann: svd_omega_i', 2, stdout, seedname)

    return

  end subroutine wann_svd_omega_i

  !==================================================================!
  subroutine wann_main_gamma(atoms, dis_window, excluded_bands, kmesh_info, k_points, out_files, &
                             wannierise, system, verbose, wann_data, m_matrix, &
                             u_matrix, u_matrix_opt, eigval, real_lattice, recip_lattice, mp_grid, &
                             num_bands, num_kpts, num_wann, have_disentangled, seedname, &
                             stdout, comm)
    !==================================================================!
    !                                                                  !
    ! Calculate the Unitary Rotations to give                          !
    !            Maximally Localised Wannier Functions                 !
    !                      Gamma version                               !
    !===================================================================
    use w90_constants, only: dp, cmplx_1, cmplx_0
    use w90_io, only: io_error, io_time, io_stopwatch
    use wannier_param_types, only: wannierise_type, output_file_type
    use w90_param_types, only: kmesh_info_type, print_output_type, &
      wannier_data_type, atom_data_type, k_point_type, disentangle_manifold_type, w90_system_type, &
      exclude_bands_type
    use wannier_methods, only: param_write_chkpt
    use w90_utility, only: utility_frac_to_cart, utility_zgemm
    use w90_comms, only: w90commtype

    implicit none

    !JJ this function has not yet been pllelised

    ! passed variables
    type(wannier_data_type), intent(inout) :: wann_data
    type(w90commtype), intent(in) :: comm
    type(wannierise_type), intent(inout) :: wannierise
    type(exclude_bands_type), intent(in) :: excluded_bands
    type(w90_system_type), intent(in) :: system
    type(print_output_type), intent(in) :: verbose
    type(k_point_type), intent(in) :: k_points ! needed for write_chkpt
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(output_file_type), intent(in) :: out_files
    type(disentangle_manifold_type), intent(in) :: dis_window ! needed for write_chkpt
    type(atom_data_type), intent(in) :: atoms

    integer, intent(in) :: stdout
    integer, intent(in) :: num_wann
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_bands
    integer, intent(in) :: mp_grid(3) ! needed for write_chkpt

    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: eigval(:, :)

    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)

    logical, intent(in) :: have_disentangled
    character(len=50), intent(in) :: seedname

    ! local variables
    type(localisation_vars) :: old_spread
    type(localisation_vars) :: wann_spread

    integer :: counts(0:0)
    integer :: displs(0:0)
    real(kind=dp), allocatable :: rnkb(:, :, :)
    real(kind=dp), allocatable :: ln_tmp(:, :, :)
    complex(kind=dp), allocatable :: m_matrix_loc(:, :, :, :)
    logical :: first_pass

    ! guiding centres
    real(kind=dp), allocatable :: rguide(:, :)
    integer :: irguide

    ! local arrays used and passed in subroutines
    real(kind=dp), allocatable :: m_w(:, :, :)
    complex(kind=dp), allocatable :: csheet(:, :, :)
    real(kind=dp), allocatable :: sheet(:, :, :)
    real(kind=dp), allocatable :: rave(:, :), r2ave(:), rave2(:)

    !local arrays not passed into subroutines
    complex(kind=dp), allocatable  :: u0(:, :, :)
    complex(kind=dp), allocatable  :: uc_rot(:, :)
    real(kind=dp), allocatable  :: ur_rot(:, :)
    complex(kind=dp), allocatable  :: cz(:, :)

    real(kind=dp) :: sqwb
    integer :: i, n, nn, iter, ind, ierr, iw
    integer :: tnntot
    logical :: lprint, ldump
    real(kind=dp), allocatable :: history(:)
    logical :: lconverged

    if (verbose%timing_level > 0) call io_stopwatch('wann: main_gamma', 1, stdout, seedname)

    first_pass = .true.

    ! Allocate stuff

    allocate (history(wannierise%control%conv_window), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating history in wann_main_gamma', stdout, seedname)

!~    if (.not.allocated(ph_g)) then
!~       allocate(  ph_g(num_wann),stat=ierr )
!~       if (ierr/=0) call io_error('Error in allocating ph_g in wann_main_gamma')
!~       ph_g = cmplx_1
!~    endif

    allocate (rnkb(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rnkb in wann_main_gamma', stdout, seedname)
    allocate (ln_tmp(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ln_tmp in wann_main_gamma', stdout, seedname)

    rnkb = 0.0_dp
    tnntot = 2*kmesh_info%nntot

    ! sub vars passed into other subs
    allocate (m_w(num_wann, num_wann, tnntot), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating m_w in wann_main_gamma', stdout, seedname)
    allocate (csheet(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating csheet in wann_main_gamma', stdout, seedname)
    allocate (sheet(num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating sheet in wann_main_gamma', stdout, seedname)
    allocate (rave(3, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rave in wann_main_gamma', stdout, seedname)
    allocate (r2ave(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating r2ave in wann_main_gamma', stdout, seedname)
    allocate (rave2(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rave2 in wann_main_gamma', stdout, seedname)
    allocate (rguide(3, num_wann))
    if (ierr /= 0) call io_error('Error in allocating rguide in wann_main_gamma', stdout, seedname)

    csheet = cmplx_1
    sheet = 0.0_dp; rave = 0.0_dp; r2ave = 0.0_dp; rave2 = 0.0_dp; rguide = 0.0_dp

    ! sub vars not passed into other subs
    allocate (u0(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating u0 in wann_main_gamma', stdout, seedname)
    allocate (uc_rot(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating uc_rot in wann_main_gamma', stdout, seedname)
    allocate (ur_rot(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating ur_rot in wann_main_gamma', stdout, seedname)
    allocate (cz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cz in wann_main_gamma', stdout, seedname)

    cz = cmplx_0

    ! Set up the MPI arrays for a serial run.
    !allocate (counts(0:0), displs(0:0), stat=ierr)
    !if (ierr /= 0) call io_error('Error in allocating counts and displs in wann_main_gamma')
    counts(0) = 1; displs(0) = 0

    ! store original U before rotating
!~    ! phase factor ph_g is applied to u_matrix
!~    ! NB: ph_g is applied to u_matrix_opt if (have_disentangled)
!~    if (have_disentangled) then
!~       u0=u_matrix
!~    else
!~       do iw=1,num_wann
!~          u0(iw,:,:)= conjg(ph_g(iw))*u_matrix(iw,:,:)
!~       end do
!~    endif
    u0 = u_matrix

!~    lguide = .false.
    ! guiding centres are not neede for orthorhombic systems
    if (kmesh_info%nntot .eq. 3) wannierise%control%guiding_centres = .false.

    if (wannierise%control%guiding_centres) then
      ! initialise rguide to projection centres (Cartesians in units of Ang)
!~       if ( use_bloch_phases) then
!~          lguide = .true.
!~       else
      do n = 1, num_wann
        call utility_frac_to_cart(wannierise%proj_site(:, n), rguide(:, n), real_lattice)
      enddo
!~       endif
    endif

    write (stdout, *)
    write (stdout, '(1x,a)') '*------------------------------- WANNIERISE ---------------------------------*'
    write (stdout, '(1x,a)') '+--------------------------------------------------------------------+<-- CONV'
    if (verbose%lenconfac .eq. 1.0_dp) then
      write (stdout, '(1x,a)') '| Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV'
    else
      write (stdout, '(1x,a)') '| Iter  Delta Spread     RMS Gradient      Spread (Bohr^2)     Time  |<-- CONV'
    endif
    write (stdout, '(1x,a)') '+--------------------------------------------------------------------+<-- CONV'
    write (stdout, *)

    irguide = 0
!~    if (guiding_centres.and.(num_no_guide_iter.le.0)) then
!~       if (nntot.gt.3) call wann_phases(csheet,sheet,rguide,irguide)
!~       irguide=1
!~    endif
    if (wannierise%control%guiding_centres .and. (wannierise%control%num_no_guide_iter .le. 0)) then
      call wann_phases(csheet, sheet, rguide, irguide, num_wann, kmesh_info, num_kpts, m_matrix, &
                       .true., counts, displs, m_matrix_loc, rnkb, verbose%timing_level, &
                       stdout, seedname, verbose%iprint, comm)
      irguide = 1
    endif

    !  weight m_matrix first to reduce number of operations
    !  m_w : weighted real matrix
    do nn = 1, kmesh_info%nntot
      sqwb = sqrt(kmesh_info%wb(nn))
      m_w(:, :, 2*nn - 1) = sqwb*real(m_matrix(:, :, nn, 1), dp)
      m_w(:, :, 2*nn) = sqwb*aimag(m_matrix(:, :, nn, 1))
    end do

    ! calculate initial centers and spread
    call wann_omega_gamma(m_w, csheet, sheet, rave, r2ave, rave2, wann_spread, num_wann, &
                          kmesh_info%nntot, kmesh_info%wbtot, kmesh_info%wb, kmesh_info%bk, &
                          wannierise%omega%invariant, ln_tmp, first_pass, &
                          verbose%timing_level, stdout, seedname)

    ! public variables
    wannierise%omega%total = wann_spread%om_tot
    wannierise%omega%invariant = wann_spread%om_i
    wannierise%omega%tilde = wann_spread%om_d + wann_spread%om_od

    ! Public array of Wannier centres and spreads
    wann_data%centres = rave
    wann_data%spreads = r2ave - rave2

    iter = 0
    old_spread%om_tot = 0.0_dp

    ! print initial state
    write (stdout, '(1x,a78)') repeat('-', 78)
    write (stdout, '(1x,a)') 'Initial State'
    do iw = 1, num_wann
      write (stdout, 1000) iw, (rave(ind, iw)*verbose%lenconfac, ind=1, 3), &
        (r2ave(iw) - rave2(iw))*verbose%lenconfac**2
    end do
    write (stdout, 1001) (sum(rave(ind, :))*verbose%lenconfac, ind=1, 3), (sum(r2ave) - sum(rave2))*verbose%lenconfac**2
    write (stdout, *)
    write (stdout, '(1x,i6,2x,E12.3,19x,F18.10,3x,F8.2,2x,a)') &
      iter, (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
      wann_spread%om_tot*verbose%lenconfac**2, io_time(), '<-- CONV'
    write (stdout, '(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
      'O_D=', wann_spread%om_d*verbose%lenconfac**2, ' O_OD=', wann_spread%om_od*verbose%lenconfac**2, &
      ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
    write (stdout, '(1x,a78)') repeat('-', 78)

    lconverged = .false.

    ! initialize ur_rot
    ur_rot = 0.0_dp
    do i = 1, num_wann
      ur_rot(i, i) = 1.0_dp
    end do

    ! main iteration loop

    do iter = 1, wannierise%control%num_iter

      lprint = .false.
      if ((mod(iter, wannierise%control%num_print_cycles) .eq. 0) .or. (iter .eq. 1) &
          .or. (iter .eq. wannierise%control%num_iter)) lprint = .true.

      ldump = .false.
      if ((wannierise%control%num_dump_cycles .gt. 0) .and. &
          (mod(iter, wannierise%control%num_dump_cycles) .eq. 0)) ldump = .true.

      if (lprint .and. verbose%iprint > 0) write (stdout, '(1x,a,i6)') 'Cycle: ', iter

!~       ! initialize rguide as rave for use_bloch_phases
!~       if ( (iter.gt.num_no_guide_iter) .and. lguide ) then
!~          rguide(:,:) = rave(:,:)
!~          lguide = .false.
!~       endif
!~       if ( guiding_centres.and.(iter.gt.num_no_guide_iter) &
!~            .and.(mod(iter,num_guide_cycles).eq.0) ) then
!~          if(nntot.gt.3) call wann_phases(csheet,sheet,rguide,irguide)
!~          irguide=1
!~       endif

      if (wannierise%control%guiding_centres .and. (iter .gt. wannierise%control%num_no_guide_iter) &
          .and. (mod(iter, wannierise%control%num_guide_cycles) .eq. 0)) then
        call wann_phases(csheet, sheet, rguide, irguide, num_wann, kmesh_info, num_kpts, m_matrix, &
                         .true., counts, displs, m_matrix_loc, rnkb, verbose%timing_level, &
                         stdout, seedname, verbose%iprint, comm, m_w)
        irguide = 1
      endif

      call internal_new_u_and_m_gamma(m_w, ur_rot, tnntot, num_wann, verbose%timing_level, &
                                      stdout)

      call wann_spread_copy(wann_spread, old_spread)

      ! calculate the new centers and spread
      call wann_omega_gamma(m_w, csheet, sheet, rave, r2ave, rave2, wann_spread, num_wann, &
                            kmesh_info%nntot, kmesh_info%wbtot, kmesh_info%wb, kmesh_info%bk, &
                            wannierise%omega%invariant, ln_tmp, first_pass, &
                            verbose%timing_level, stdout, seedname)

      ! print the new centers and spreads
      if (lprint) then
        do iw = 1, num_wann
          write (stdout, 1000) iw, (rave(ind, iw)*verbose%lenconfac, ind=1, 3), &
            (r2ave(iw) - rave2(iw))*verbose%lenconfac**2
        end do
        write (stdout, 1001) (sum(rave(ind, :))*verbose%lenconfac, ind=1, 3), &
          (sum(r2ave) - sum(rave2))*verbose%lenconfac**2
        write (stdout, *)
        write (stdout, '(1x,i6,2x,E12.3,19x,F18.10,3x,F8.2,2x,a)') &
          iter, (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, &
          wann_spread%om_tot*verbose%lenconfac**2, io_time(), '<-- CONV'
        write (stdout, '(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
          'O_D=', wann_spread%om_d*verbose%lenconfac**2, &
          ' O_OD=', wann_spread%om_od*verbose%lenconfac**2, &
          ' O_TOT=', wann_spread%om_tot*verbose%lenconfac**2, ' <-- SPRD'
        write (stdout, '(1x,a,E15.7,a,E15.7,a,E15.7,a)') &
          'Delta: O_D=', (wann_spread%om_d - old_spread%om_d)*verbose%lenconfac**2, &
          ' O_OD=', (wann_spread%om_od - old_spread%om_od)*verbose%lenconfac**2, &
          ' O_TOT=', (wann_spread%om_tot - old_spread%om_tot)*verbose%lenconfac**2, ' <-- DLTA'
        write (stdout, '(1x,a78)') repeat('-', 78)
      end if

      ! Public array of Wannier centres and spreads
      wann_data%centres = rave
      wann_data%spreads = r2ave - rave2

      ! Public variables
      wannierise%omega%total = wann_spread%om_tot
      wannierise%omega%tilde = wann_spread%om_d + wann_spread%om_od

      if (ldump) then
        uc_rot(:, :) = cmplx(ur_rot(:, :), 0.0_dp, dp)
        call utility_zgemm(u_matrix, u0, 'N', uc_rot, 'N', num_wann)
        call param_write_chkpt('postdis', excluded_bands, wann_data, kmesh_info, k_points, &
                               num_kpts, dis_window, num_bands, num_wann, u_matrix, u_matrix_opt, &
                               m_matrix, mp_grid, real_lattice, recip_lattice, &
                               wannierise%omega%invariant, have_disentangled, &
                               stdout, seedname)
      endif

      if (wannierise%control%conv_window .gt. 1) then
        call internal_test_convergence_gamma(wann_spread, old_spread, history, &
                                             iter, lconverged, wannierise%control%conv_window, &
                                             wannierise%control%conv_tol, stdout)
      endif

      if (lconverged) then
        write (stdout, '(/13x,a,es10.3,a,i2,a)') &
          '<<<     Delta <', wannierise%control%conv_tol, &
          '  over ', wannierise%control%conv_window, ' iterations     >>>'
        write (stdout, '(13x,a/)') '<<< Wannierisation convergence criteria satisfied >>>'
        exit
      endif

    enddo
    ! end of the minimization loop

    ! update M
    do nn = 1, kmesh_info%nntot
      sqwb = 1.0_dp/sqrt(kmesh_info%wb(nn))
      m_matrix(:, :, nn, 1) = sqwb*cmplx(m_w(:, :, 2*nn - 1), m_w(:, :, 2*nn), dp)
    end do
    ! update U
    uc_rot(:, :) = cmplx(ur_rot(:, :), 0.0_dp, dp)
    call utility_zgemm(u_matrix, u0, 'N', uc_rot, 'N', num_wann)

    write (stdout, '(1x,a)') 'Final State'
    do iw = 1, num_wann
      write (stdout, 1000) iw, (rave(ind, iw)*verbose%lenconfac, ind=1, 3), &
        (r2ave(iw) - rave2(iw))*verbose%lenconfac**2
    end do
    write (stdout, 1001) (sum(rave(ind, :))*verbose%lenconfac, ind=1, 3), &
      (sum(r2ave) - sum(rave2))*verbose%lenconfac**2
    write (stdout, *)
    write (stdout, '(3x,a21,a,f15.9)') '     Spreads ('//trim(verbose%length_unit)//'^2)', &
      '       Omega I      = ', wann_spread%om_i*verbose%lenconfac**2
    write (stdout, '(3x,a,f15.9)') '     ================       Omega D      = ', &
      wann_spread%om_d*verbose%lenconfac**2
    write (stdout, '(3x,a,f15.9)') '                            Omega OD     = ', &
      wann_spread%om_od*verbose%lenconfac**2
    write (stdout, '(3x,a21,a,f15.9)') 'Final Spread ('//trim(verbose%length_unit)//'^2)', &
      '       Omega Total  = ', wann_spread%om_tot*verbose%lenconfac**2
    write (stdout, '(1x,a78)') repeat('-', 78)

    if (out_files%write_xyz) then
      call wann_write_xyz(out_files%translate_home_cell, num_wann, wann_data%centres, &
                          real_lattice, recip_lattice, atoms, verbose, stdout, seedname)
    endif

    if (wannierise%control%guiding_centres) then
      call wann_phases(csheet, sheet, rguide, irguide, num_wann, kmesh_info, num_kpts, m_matrix, &
                       .true., counts, displs, m_matrix_loc, rnkb, verbose%timing_level, &
                       stdout, seedname, verbose%iprint, comm)
    endif

    ! unitarity is checked
!~    call internal_check_unitarity()
    call wann_check_unitarity(num_kpts, num_wann, u_matrix, verbose%timing_level, &
                              verbose%iprint, stdout, seedname)

    ! write extra info regarding omega_invariant
!~    if (iprint>2) call internal_svd_omega_i()
    if (verbose%iprint > 2) then
      call wann_svd_omega_i(num_wann, num_kpts, kmesh_info, m_matrix, verbose, stdout, seedname)
    endif

    ! write matrix elements <m|r^2|n> to file
!~    if (write_r2mn) call internal_write_r2mn()
    if (out_files%write_r2mn) call wann_write_r2mn(num_kpts, num_wann, kmesh_info, m_matrix, &
                                                   stdout, seedname)

    ! calculate and write projection of WFs on original bands in outer window
    if (have_disentangled .and. out_files%write_proj) &
      call wann_calc_projection(num_bands, num_wann, num_kpts, u_matrix_opt, eigval, &
                                dis_window%lwindow, verbose%timing_level, verbose%iprint, &
                                stdout, seedname)

    ! aam: write data required for vdW utility
    if (out_files%write_vdw_data) then
      call wann_write_vdw_data(num_wann, wann_data, real_lattice, recip_lattice, u_matrix, &
                               u_matrix_opt, have_disentangled, system, stdout, seedname)
    endif

    ! deallocate sub vars not passed into other subs
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cz in wann_main_gamma', stdout, seedname)
    deallocate (ur_rot, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ur_rot in wann_main_gamma', stdout, seedname)
    deallocate (uc_rot, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating uc_rot in wann_main_gamma', stdout, seedname)
    deallocate (u0, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating u0 in wann_main_gamma', stdout, seedname)

    ! deallocate sub vars passed into other subs
    deallocate (rguide, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rguide in wann_main_gamma', stdout, seedname)
    deallocate (rave2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rave2 in wann_main_gamma', stdout, seedname)
    deallocate (rave, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rave in wann_main_gamma', stdout, seedname)
    deallocate (sheet, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating sheet in wann_main_gamma', stdout, seedname)
    deallocate (csheet, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating csheet in wann_main_gamma', stdout, seedname)
    deallocate (m_w, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating m_w in wann_main_gamma', stdout, seedname)

    ! deallocate module data
    deallocate (ln_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ln_tmp in wann_main_gamma', stdout, seedname)
    deallocate (rnkb, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rnkb in wann_main_gamma', stdout, seedname)

    deallocate (history, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating history in wann_main_gamma', stdout, seedname)

    if (verbose%timing_level > 0) call io_stopwatch('wann: main_gamma', 2, stdout, seedname)

    return

1000 format(2x, 'WF centre and spread', &
&       i5, 2x, '(', f10.6, ',', f10.6, ',', f10.6, ' )', f15.8)

1001 format(2x, 'Sum of centres and spreads', &
&       1x, '(', f10.6, ',', f10.6, ',', f10.6, ' )', f15.8)

  contains

    !===============================================!
    subroutine internal_new_u_and_m_gamma(m_w, ur_rot, tnntot, num_wann, &
                                          timing_level, stdout)
      !===============================================!

      use w90_constants, only: pi, eps10
      use w90_io, only: io_stopwatch

      implicit none
      real(kind=dp), intent(inout) :: m_w(:, :, :)
      real(kind=dp), intent(inout)  :: ur_rot(:, :)
      integer, intent(in) :: tnntot
      integer, intent(in) :: num_wann
      integer, intent(in) :: timing_level
      integer, intent(in) :: stdout
      ! local
      real(kind=dp) :: theta, twotheta
      real(kind=dp) :: a11, a12, a21, a22
      real(kind=dp) :: cc, ss, rtmp1, rtmp2
      real(kind=dp), parameter :: pifour = 0.25_dp*pi
      integer       :: nn, nw1, nw2, nw3

      if (timing_level > 1) call io_stopwatch('wann: main_gamma: new_u_and_m_gamma', 1, stdout, seedname)

      loop_nw1: do nw1 = 1, num_wann
      loop_nw2: do nw2 = nw1 + 1, num_wann

        a11 = 0.0_dp; a12 = 0.0_dp; a22 = 0.0_dp
        do nn = 1, tnntot
          a11 = a11 + (m_w(nw1, nw1, nn) - m_w(nw2, nw2, nn))**2
          a12 = a12 + m_w(nw1, nw2, nn)*(m_w(nw1, nw1, nn) - m_w(nw2, nw2, nn))
          a22 = a22 + m_w(nw1, nw2, nn)**2
        end do
        a12 = 2.0_dp*a12
        a22 = 4.0_dp*a22
        a21 = a22 - a11
        if (abs(a12) .gt. eps10) then
          twotheta = 0.5_dp*(a21 + sqrt(a21**2 + 4.0_dp*a12**2))/a12
          theta = 0.5_dp*atan(twotheta)
        elseif (a21 .lt. eps10) then
          theta = 0.0_dp
        else
          theta = pifour
        endif
        cc = cos(theta)
        ss = sin(theta)

        ! update M
        do nn = 1, tnntot
          ! MR
          do nw3 = 1, num_wann
            rtmp1 = m_w(nw3, nw1, nn)*cc + m_w(nw3, nw2, nn)*ss
            rtmp2 = -m_w(nw3, nw1, nn)*ss + m_w(nw3, nw2, nn)*cc
            m_w(nw3, nw1, nn) = rtmp1
            m_w(nw3, nw2, nn) = rtmp2
          end do
          ! R^+ M R
          do nw3 = 1, num_wann
            rtmp1 = cc*m_w(nw1, nw3, nn) + ss*m_w(nw2, nw3, nn)
            rtmp2 = -ss*m_w(nw1, nw3, nn) + cc*m_w(nw2, nw3, nn)
            m_w(nw1, nw3, nn) = rtmp1
            m_w(nw2, nw3, nn) = rtmp2
          end do
        end do
        ! update U : U=UR
        do nw3 = 1, num_wann
          rtmp1 = ur_rot(nw3, nw1)*cc + ur_rot(nw3, nw2)*ss
          rtmp2 = -ur_rot(nw3, nw1)*ss + ur_rot(nw3, nw2)*cc
          ur_rot(nw3, nw1) = rtmp1
          ur_rot(nw3, nw2) = rtmp2
        end do
      end do loop_nw2
      end do loop_nw1

      if (timing_level > 1) call io_stopwatch('wann: main_gamma: new_u_and_m_gamma', 2, stdout, seedname)

      return

    end subroutine internal_new_u_and_m_gamma

    !===============================================!
    subroutine internal_test_convergence_gamma(wann_spread, old_spread, history, &
                                               iter, lconverged, conv_window, conv_tol, stdout)
      !===============================================!
      !                                               !
      ! Determine whether minimisation of non-gauge-  !
      ! invariant spread is converged                 !
      !                                               !
      !===============================================!
      use w90_io, only: io_error

      implicit none
      type(localisation_vars), intent(in) :: wann_spread
      type(localisation_vars), intent(in) :: old_spread
      real(kind=dp), intent(inout) :: history(:)
      integer, intent(in) :: iter
      logical, intent(out) :: lconverged
      integer, intent(in) :: conv_window
      real(kind=dp), intent(in) :: conv_tol
      integer, intent(in) :: stdout
      ! local
      real(kind=dp) :: delta_omega
      integer :: j, ierr
      real(kind=dp), allocatable :: temp_hist(:)

      allocate (temp_hist(conv_window), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating temp_hist in wann_main', stdout, seedname)

      delta_omega = wann_spread%om_tot - old_spread%om_tot

      if (iter .le. conv_window) then
        history(iter) = delta_omega
      else
        temp_hist = eoshift(history, 1, delta_omega)
        history = temp_hist
      endif

      lconverged = .false.

      if (iter .ge. conv_window) then
!~         write(stdout,*) (history(j),j=1,conv_window)
        do j = 1, conv_window
          if (abs(history(j)) .gt. conv_tol) exit
          lconverged = .true.
        enddo
      endif

      deallocate (temp_hist, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating temp_hist in wann_main_gamma', stdout, seedname)

      return

    end subroutine internal_test_convergence_gamma

!~    !========================================!
!~    subroutine internal_check_unitarity()
!~    !========================================!
!~
!~      implicit none
!~
!~      integer :: nkp,i,j,m
!~      complex(kind=dp) :: ctmp1,ctmp2
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',1)
!~
!~      do nkp = 1, num_kpts
!~         do i = 1, num_wann
!~            do j = 1, num_wann
!~               ctmp1 = cmplx_0
!~               ctmp2 = cmplx_0
!~               do m = 1, num_wann
!~                  ctmp1 = ctmp1 + u_matrix (i, m, nkp) * conjg (u_matrix (j, m, nkp) )
!~                  ctmp2 = ctmp2 + u_matrix (m, j, nkp) * conjg (u_matrix (m, i, nkp) )
!~               enddo
!~               if ( (i.eq.j) .and. (abs (ctmp1 - cmplx_1 ) .gt. eps5) ) &
!~                    then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp1
!~                  call io_error('wann_main: unitariety error 1')
!~               endif
!~               if ( (i.eq.j) .and. (abs (ctmp2 - cmplx_1 ) .gt. eps5) ) &
!~                    then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp2
!~                  call io_error('wann_main: unitariety error 2')
!~               endif
!~               if ( (i.ne.j) .and. (abs (ctmp1) .gt. eps5) ) then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp1
!~                  call io_error('wann_main: unitariety error 3')
!~               endif
!~               if ( (i.ne.j) .and. (abs (ctmp2) .gt. eps5) ) then
!~                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!~                       ctmp2
!~                  call io_error('wann_main: unitariety error 4')
!~               endif
!~            enddo
!~         enddo
!~      enddo
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',2)
!~
!~      return
!~
!~    end subroutine internal_check_unitarity

!~    !========================================!
!~    subroutine internal_write_r2mn()
!~    !========================================!
!~    !                                        !
!~    ! Write seedname.r2mn file               !
!~    !                                        !
!~    !========================================!
!~      use w90_io, only: seedname,io_file_unit,io_error
!~
!~      implicit none
!~
!~      integer :: r2mnunit,nw1,nw2,nkp,nn
!~      real(kind=dp) :: r2ave_mn,delta
!~
!~      ! note that here I use formulas analogue to Eq. 23, and not to the
!~      ! shift-invariant Eq. 32 .
!~      r2mnunit=io_file_unit()
!~      open(r2mnunit,file=trim(seedname)//'.r2mn',form='formatted',err=158)
!~      do nw1 = 1, num_wann
!~         do nw2 = 1, num_wann
!~            r2ave_mn = 0.0_dp
!~            delta = 0.0_dp
!~            if (nw1.eq.nw2) delta = 1.0_dp
!~            do nkp = 1, num_kpts
!~               do nn = 1, nntot
!~                  r2ave_mn = r2ave_mn + wb(nn) * &
!~                       ! [GP-begin, Apr13, 2012: corrected sign inside "real"]
!~                       ( 2.0_dp * delta - real(m_matrix(nw1,nw2,nn,nkp) + &
!~                       conjg(m_matrix(nw2,nw1,nn,nkp)),kind=dp) )
!~               enddo
!~            enddo
!~            r2ave_mn = r2ave_mn / real(num_kpts,dp)
!~            write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn
!~         enddo
!~      enddo
!~      close(r2mnunit)
!~
!~      return
!~
!~158   call io_error('Error opening file '//trim(seedname)//'.r2mn in wann_main')
!~
!~    end subroutine internal_write_r2mn

!~    !========================================!
!~    subroutine internal_svd_omega_i()
!~    !========================================!
!~
!~      implicit none
!~
!~      complex(kind=dp), allocatable  :: cv1(:,:),cv2(:,:)
!~      complex(kind=dp), allocatable  :: cw1(:),cw2(:)
!~      complex(kind=dp), allocatable  :: cpad1 (:)
!~      real(kind=dp),    allocatable  :: singvd (:)
!~
!~      integer :: nkp,nn,nb,na,ind
!~      real(kind=dp) :: omt1,omt2,omt3
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',1)
!~
!~      allocate( cw1 (10 * num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cw1 in wann_main')
!~      allocate( cw2 (10 * num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cw2 in wann_main')
!~      allocate( cv1 (num_wann, num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cv1 in wann_main')
!~      allocate( cv2 (num_wann, num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cv2 in wann_main')
!~      allocate( singvd (num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating singvd in wann_main')
!~      allocate( cpad1 (num_wann * num_wann),stat=ierr  )
!~      if (ierr/=0) call io_error('Error in allocating cpad1 in wann_main')
!~
!~      cw1=cmplx_0; cw2=cmplx_0; cv1=cmplx_0; cv2=cmplx_0; cpad1=cmplx_0
!~      singvd=0.0_dp
!~
!~      ! singular value decomposition
!~      omt1 = 0.0_dp ; omt2 = 0.0_dp ; omt3 = 0.0_dp
!~      do nkp = 1, num_kpts
!~         do nn = 1, nntot
!~            ind = 1
!~            do nb = 1, num_wann
!~               do na = 1, num_wann
!~                  cpad1 (ind) = m_matrix (na, nb, nn, nkp)
!~                  ind = ind+1
!~               enddo
!~            enddo
!~            call zgesvd ('A', 'A', num_wann, num_wann, cpad1, num_wann, singvd, cv1, &
!~                 num_wann, cv2, num_wann, cw1, 10 * num_wann, cw2, info)
!~            if (info.ne.0) then
!~               call io_error('ERROR: Singular value decomp. zgesvd failed')
!~            endif
!~
!~            do nb = 1, num_wann
!~               omt1 = omt1 + wb(nn) * (1.0_dp - singvd (nb) **2)
!~               omt2 = omt2 - wb(nn) * (2.0_dp * log (singvd (nb) ) )
!~               omt3 = omt3 + wb(nn) * (acos (singvd (nb) ) **2)
!~            enddo
!~         enddo
!~      enddo
!~      omt1 = omt1 / real(num_kpts,dp)
!~      omt2 = omt2 / real(num_kpts,dp)
!~      omt3 = omt3 / real(num_kpts,dp)
!~      write ( stdout , * ) ' '
!~      write(stdout,'(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ',&
!~           omt1*lenconfac**2,'('//trim(length_unit)//'^2)'
!~      write(stdout,'(2x,a,f15.9,1x,a)') '                 -2log s = ',&
!~           omt2*lenconfac**2,'('//trim(length_unit)//'^2)'
!~      write(stdout,'(2x,a,f15.9,1x,a)') '                  acos^2 = ',&
!~           omt3*lenconfac**2,'('//trim(length_unit)//'^2)'
!~
!~      deallocate(cpad1,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cpad1 in wann_main')
!~      deallocate(singvd,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating singvd in wann_main')
!~      deallocate(cv2,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cv2 in wann_main')
!~      deallocate(cv1,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cv1 in wann_main')
!~      deallocate(cw2,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cw2 in wann_main')
!~      deallocate(cw1,stat=ierr)
!~      if (ierr/=0) call io_error('Error in deallocating cw1 in wann_main')
!~
!~      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',2)
!~
!~      return
!~
!~    end subroutine internal_svd_omega_i

  end subroutine wann_main_gamma

  !==================================================================!
  subroutine wann_omega_gamma(m_w, csheet, sheet, rave, r2ave, rave2, wann_spread, num_wann, &
                              nntot, wbtot, wb, bk, omega_invariant, ln_tmp, first_pass, &
                              timing_level, stdout, seedname)
    !==================================================================!
    !                                                                  !
    !   Calculate the Wannier Function spread                          !
    !                                                                  !
    !===================================================================
    use w90_io, only: io_error, io_stopwatch

    implicit none

    ! passed variables
    type(localisation_vars), intent(out)  :: wann_spread

    integer, intent(in) :: timing_level
    integer, intent(in) :: stdout
    integer, intent(in) :: num_wann
    integer, intent(in) :: nntot

    real(kind=dp), intent(out) :: rave2(:)
    real(kind=dp), intent(out) :: rave(:, :)
    real(kind=dp), intent(out) :: r2ave(:)
    real(kind=dp), intent(out) :: ln_tmp(:, :, :)
    real(kind=dp), intent(in) :: wbtot
    real(kind=dp), intent(in) :: wb(:)
    real(kind=dp), intent(in) :: sheet(:, :, :)
    real(kind=dp), intent(in) :: omega_invariant
    real(kind=dp), intent(in) :: m_w(:, :, :)
    real(kind=dp), intent(in) :: bk(:, :, :)

    complex(kind=dp), intent(in)  :: csheet(:, :, :)

    logical, intent(inout) :: first_pass

    character(len=50), intent(in)  :: seedname

    !local variables
    real(kind=dp) :: summ, brn
    real(kind=dp), allocatable :: m_w_nn2(:)
    integer :: ind, nn, m, n, iw, rn, cn, ierr

    if (timing_level > 1) call io_stopwatch('wann: omega_gamma', 1, stdout, seedname)

    allocate (m_w_nn2(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating m_w_nn2 in wann_omega_gamma', stdout, seedname)

    if (nntot .eq. 3) then
      do nn = 1, nntot
        rn = 2*nn - 1
        cn = 2*nn
        do n = 1, num_wann
          ln_tmp(n, nn, 1) = atan2(m_w(n, n, cn), m_w(n, n, rn))
        end do
      end do
    else
      do nn = 1, nntot
        rn = 2*nn - 1
        cn = 2*nn
        do n = 1, num_wann
          ln_tmp(n, nn, 1) = aimag(log(csheet(n, nn, 1)*cmplx(m_w(n, n, rn), m_w(n, n, cn), dp))) &
                             - sheet(n, nn, 1)
        end do
      end do
    endif

    rave = 0.0_dp
    do iw = 1, num_wann
      do ind = 1, 3
        do nn = 1, nntot
          rave(ind, iw) = rave(ind, iw) - wb(nn)*bk(ind, nn, 1) &
                          *ln_tmp(iw, nn, 1)
        enddo
      enddo
    enddo

    rave2 = 0.0_dp
    do iw = 1, num_wann
      rave2(iw) = sum(rave(:, iw)*rave(:, iw))
    enddo

    m_w_nn2 = 0.0_dp
    r2ave = wbtot
    do iw = 1, num_wann
      do nn = 1, nntot
        rn = 2*nn - 1
        cn = 2*nn
        m_w_nn2(iw) = m_w_nn2(iw) + m_w(iw, iw, rn)**2 + m_w(iw, iw, cn)**2
        r2ave(iw) = r2ave(iw) + wb(nn)*ln_tmp(iw, nn, 1)**2
      enddo
      r2ave(iw) = r2ave(iw) - m_w_nn2(iw)
    enddo

    if (first_pass) then
      summ = 0.0_dp
      do nn = 1, nntot
        rn = 2*nn - 1
        cn = 2*nn
        do m = 1, num_wann
          do n = 1, num_wann
            summ = summ + m_w(n, m, rn)**2 + m_w(n, m, cn)**2
          enddo
        enddo
      enddo
      wann_spread%om_i = wbtot*real(num_wann, dp) - summ
      first_pass = .false.
    else
      wann_spread%om_i = omega_invariant
    endif

    wann_spread%om_od = wbtot*real(num_wann, dp) - sum(m_w_nn2(:)) - wann_spread%om_i

    if (nntot .eq. 3) then
      wann_spread%om_d = 0.0_dp
    else
      wann_spread%om_d = 0.0_dp
      do nn = 1, nntot
        do n = 1, num_wann
          brn = sum(bk(:, nn, 1)*rave(:, n))
          wann_spread%om_d = wann_spread%om_d + wb(nn)*(ln_tmp(n, nn, 1) + brn)**2
        enddo
      enddo
    end if

    wann_spread%om_tot = wann_spread%om_i + wann_spread%om_d + wann_spread%om_od

    deallocate (m_w_nn2, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating m_w_nn2 in wann_omega_gamma', &
                                 stdout, seedname)

    if (timing_level > 1) call io_stopwatch('wann: omega_gamma', 2, stdout, seedname)

    return

  end subroutine wann_omega_gamma

end module w90_wannierise
