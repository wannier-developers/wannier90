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
!                                                            !
!  w90_disentangle: extract subspace from entangled bands    !
!                                                            !
!------------------------------------------------------------!

module w90_disentangle_mod
  !! This module contains the core routines to extract an optimal
  !! subspace from a set of entangled bands.

  implicit none

  public :: dis_main
  public :: setup_m_loc

contains
  !================================================!

  subroutine dis_main(dis_control, dis_spheres, dis_manifold, kmesh_info, kpt_latt, sitesym, &
                      print_output, m_matrix_orig_local, u_matrix, u_matrix_opt, eigval, &
                      real_lattice, omega_invariant, num_bands, num_kpts, num_wann, gamma_only, &
                      lsitesymmetry, stdout, timer, dist_k, error, comm)
    !================================================!
    !
    !! Main disentanglement routine
    !
    !================================================!
    use w90_comms, only: comms_bcast, w90_comm_type, mpirank
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_sitesym, only: sitesym_replace_d_matrix_band, sitesym_symmetrize_u_matrix, &
      sitesym_symmetrize_zmatrix, sitesym_dis_extract_symmetry
    use w90_types, only: dis_manifold_type, kmesh_info_type, print_output_type, timer_list_type
    use w90_utility, only: utility_recip_lattice_base
    use w90_wannier90_types, only: dis_control_type, dis_spheres_type, sitesym_type

    ! arguments
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: dist_k(:)

    logical, intent(in) :: lsitesymmetry
    logical, intent(in) :: gamma_only

    real(kind=dp), pointer, intent(in) :: eigval(:, :) ! (num_bands, num_kpts)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(inout) :: omega_invariant
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    !complex(kind=dp), intent(inout) :: a_matrix(:, :, :) ! (num_bands, num_wann, num_kpts)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :) ! (num_wann, num_wann, num_kpts)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :) ! (num_bands, num_wann, num_kpts)
    complex(kind=dp), intent(inout) :: m_matrix_orig_local(:, :, :, :) ! this is the only "m matrix" here now

    type(dis_control_type), intent(inout)  :: dis_control
    type(dis_manifold_type), intent(inout) :: dis_manifold
    type(dis_spheres_type), intent(in) :: dis_spheres
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(sitesym_type), intent(inout) :: sitesym
    type(w90_comm_type), intent(in) :: comm
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    ! internal variables
    real(kind=dp) :: recip_lattice(3, 3), volume
    integer :: nkp, nkp2, nn, j, ierr, nkp_global
    logical :: linner                         !! Is there a frozen window
    logical :: lfrozen(num_bands, num_kpts)   !! true if the i-th band inside outer window is frozen
    integer :: ndimfroz(num_kpts)             !! number of frozen bands at nkp-th k point
    integer :: indxfroz(num_bands, num_kpts)  !! number of bands inside outer window at nkp-th k point
    integer :: indxnfroz(num_bands, num_kpts) !! outer-window band index for the i-th non-frozen state
    complex(kind=dp), allocatable :: a_matrix(:, :, :) ! (num_bands, num_wann, num_kpts)
    !! (equals 1 if it is the bottom of outer window)

    real(kind=dp), allocatable :: eigval_opt(:, :)  !! At input it contains a large set of eigenvalues. At
    !! it is slimmed down to contain only those inside the energy window.

    complex(kind=dp), allocatable :: cwb(:, :), cww(:, :)

    ! pllel setup
    integer :: nkrank, ikg, ikl, my_node_id
    integer, allocatable :: global_k(:)
    logical :: on_root = .false.

    my_node_id = mpirank(comm)
    on_root = (my_node_id == 0)
    nkrank = count(dist_k == my_node_id) ! this routine must proceed also in the case of zero k-points this rank, to ensure collective communications are matched

    allocate (a_matrix(num_bands, num_wann, num_kpts), stat=ierr) ! a_matrix is local to disentangle()
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating a_matrix in dis_main', comm)
      return
    endif
    a_matrix = u_matrix_opt ! initial projections are passed to this routine via u_matrix_opt

    allocate (global_k(nkrank), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating global_k in dis_main', comm)
      return
    endif
    global_k = huge(1); ikl = 1
    do ikg = 1, num_kpts
      if (dist_k(ikg) == my_node_id) then
        global_k(ikl) = ikg
        ikl = ikl + 1
      endif
    enddo

    if (print_output%timing_level > 0) call io_stopwatch_start('dis: main', timer)

    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)

    if (print_output%iprint > 0) write (stdout, '(/1x,a)') &
      '*------------------------------- DISENTANGLE --------------------------------*'

    ! Allocate arrays
    allocate (eigval_opt(num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating eigval_opt in dis_main', comm)
      return
    endif
    eigval_opt(1:num_bands, 1:num_kpts) = eigval(1:num_bands, 1:num_kpts)

    ! Set up energy windows
    call dis_windows(dis_spheres, dis_manifold, eigval_opt, kpt_latt, recip_lattice, indxfroz, &
                     indxnfroz, ndimfroz, print_output%iprint, num_bands, num_kpts, &
                     num_wann, print_output%timing_level, lfrozen, linner, on_root, &
                     stdout, timer, error, comm)
    if (allocated(error)) return

    ! Construct the unitarized projection
    call dis_project(a_matrix, u_matrix_opt, dis_manifold%ndimwin, dis_manifold%nfirstwin, num_bands, num_kpts, &
                     num_wann, print_output%timing_level, on_root, print_output%iprint, timer, &
                     error, stdout, comm)
    if (allocated(error)) return

    ! If there is an inner window, need to modify projection procedure
    ! (Sec. III.G SMV)
    if (linner) then
      if (lsitesymmetry) then
        call set_error_fatal(error, 'in symmetry-adapted mode, frozen window not implemented yet', &
                             comm)
        return
      endif
      if (print_output%iprint > 0) write (stdout, '(3x,a)') 'Using an inner window (linner = T)'
      call dis_proj_froz(u_matrix_opt, indxfroz, ndimfroz, dis_manifold%ndimwin, &
                         print_output%iprint, num_bands, num_kpts, num_wann, &
                         print_output%timing_level, lfrozen, on_root, timer, error, stdout, comm)
      if (allocated(error)) return
    else
      if (print_output%iprint > 0) write (stdout, '(3x,a)') 'No inner window (linner = F)'
    endif

    ! Debug
    call internal_check_orthonorm(u_matrix_opt, dis_manifold%ndimwin, num_kpts, num_wann, &
                                  print_output%timing_level, on_root, timer, error, stdout, comm)
    if (allocated(error)) return

    ! Slim down the original Mmn(k,b)
    call internal_slim_m(m_matrix_orig_local, dis_manifold%ndimwin, dis_manifold%nfirstwin, kmesh_info%nnlist, &
                         kmesh_info%nntot, num_bands, print_output%timing_level, timer, dist_k, &
                         global_k, error, comm)
    if (allocated(error)) return

    dis_manifold%lwindow = .false.
    do nkp = 1, num_kpts
      do j = dis_manifold%nfirstwin(nkp), dis_manifold%nfirstwin(nkp) + dis_manifold%ndimwin(nkp) - 1
        dis_manifold%lwindow(j, nkp) = .true.
      end do
    end do

    if (lsitesymmetry) then
      call sitesym_symmetrize_u_matrix(sitesym, u_matrix_opt, num_bands, num_bands, num_kpts, &
                                       num_wann, stdout, error, comm, dis_manifold%lwindow)
      if (allocated(error)) return
    endif

    !RS: calculate initial U_{opt}(Rk) from U_{opt}(k)
    ! Extract the optimally-connected num_wann-dimensional subspaces

    if (gamma_only) then
      call dis_extract_gamma(dis_control, kmesh_info, print_output, dis_manifold, &
                             m_matrix_orig_local, u_matrix_opt, eigval_opt, omega_invariant, &
                             indxnfroz, ndimfroz, num_bands, num_kpts, num_wann, timer, error, &
                             stdout, comm)
      if (allocated(error)) return
    else
      call dis_extract(dis_control, kmesh_info, sitesym, print_output, dis_manifold, &
                       m_matrix_orig_local, u_matrix_opt, eigval_opt, omega_invariant, indxnfroz, &
                       ndimfroz, my_node_id, num_bands, num_kpts, num_wann, lsitesymmetry, timer, &
                       nkrank, global_k, error, stdout, comm)
      if (allocated(error)) return
    end if

    ! Allocate workspace
    allocate (cwb(num_wann, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cwb in dis_main', comm)
      return
    endif
    allocate (cww(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cww in dis_main', comm)
      return
    endif

    ! Find the num_wann x num_wann overlap matrices between
    ! the basis states of the optimal subspaces
    do nkp = 1, nkrank
      nkp_global = global_k(nkp)
      do nn = 1, kmesh_info%nntot
        nkp2 = kmesh_info%nnlist(nkp_global, nn)
        call zgemm('C', 'N', num_wann, dis_manifold%ndimwin(nkp2), dis_manifold%ndimwin(nkp_global), &
                   cmplx_1, u_matrix_opt(:, :, nkp_global), num_bands, &
                   m_matrix_orig_local(:, :, nn, nkp), num_bands, cmplx_0, cwb, num_wann)
        call zgemm('N', 'N', num_wann, num_wann, dis_manifold%ndimwin(nkp2), cmplx_1, cwb, &
                   num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
        m_matrix_orig_local(1:num_wann, 1:num_wann, nn, nkp) = cww(:, :)
      enddo
    enddo

    ! Find the initial u_matrix
    if (lsitesymmetry) call sitesym_replace_d_matrix_band(sitesym, num_wann) !RS: replace d_matrix_band here

    if (gamma_only) then
      call internal_find_u_gamma(a_matrix, u_matrix, u_matrix_opt, dis_manifold%ndimwin, num_wann, &
                                 print_output%timing_level, stdout, timer, error, comm)
      if (allocated(error)) return
    else
      call internal_find_u(sitesym, a_matrix, u_matrix, u_matrix_opt, dis_manifold%ndimwin, &
                           num_bands, num_kpts, num_wann, print_output%timing_level, &
                           lsitesymmetry, on_root, stdout, timer, error, comm)
      if (allocated(error)) return
    end if

    !zero the unused elements of u_matrix_opt (just in case...)
    do nkp = 1, num_kpts
      do j = 1, num_wann
        if (dis_manifold%ndimwin(nkp) < num_bands) &
          u_matrix_opt(dis_manifold%ndimwin(nkp) + 1:, j, nkp) = cmplx_0
      enddo
    enddo

    ! Deallocate workspace
    deallocate (cww, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cww in dis_main', comm)
      return
    endif
    deallocate (cwb, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cwb in dis_main', comm)
      return
    endif
    deallocate (global_k, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating global_k in dis_main', comm)
      return
    endif
    deallocate (a_matrix, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating a_matrix in dis_main', comm)
      return
    endif

    if (print_output%timing_level > 0 .and. on_root) call io_stopwatch_stop('dis: main', timer)

    return
    !================================================!
  end subroutine dis_main

  subroutine setup_m_loc(kmesh_info, print_output, m_matrix_local, m_matrix_orig_local, u_matrix, &
                         num_bands, num_kpts, num_wann, timer, dist_k, error, comm)
    !================================================!
    !
    ! map m_matrix_orig_local to m_matrix_local
    !  at entry, m_matrix_local is not allocated
    !
    !================================================!
    use w90_comms, only: w90_comm_type, mpirank
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_info_type, print_output_type, timer_list_type

    ! arguments
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(in) :: dist_k(:)

    complex(kind=dp), intent(in) :: u_matrix(:, :, :) ! (num_wann, num_wann, num_kpts) -- full array duplicated on all ranks
    complex(kind=dp), intent(in) :: m_matrix_orig_local(:, :, :, :) ! (num_bands, num_bands, nntot, num_kpts) -- only local kpts
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :) ! (num_wann, num_wann, nntot, rank_kpts) -- only local kpts

    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(w90_comm_type), intent(in) :: comm
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    ! internal variables
    complex(kind=dp), allocatable :: cwb(:, :), cww(:, :)
    integer :: nkp, nkp2, nn, ierr, nkp_global, nkrank
    integer, allocatable :: global_k(:)
    integer :: ikg, ikl, my_node_id

    if (print_output%timing_level > 1) call io_stopwatch_start('dis: splitm', timer)

    ! local-global k index mapping
    my_node_id = mpirank(comm)
    nkrank = count(dist_k == my_node_id)
    allocate (global_k(nkrank), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating global_k in dis_main', comm)
      return
    endif
    global_k = huge(1); ikl = 1
    do ikg = 1, num_kpts
      if (dist_k(ikg) == my_node_id) then
        global_k(ikl) = ikg ! global [1,num_kpts] index corresponding to local [1,nk_this_node] index
        ikl = ikl + 1
      endif
    enddo

    allocate (cwb(num_wann, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cwb in dis_main', comm)
      return
    endif
    allocate (cww(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cww in dis_main', comm)
      return
    endif

    do nkp = 1, nkrank
      nkp_global = global_k(nkp)

      do nn = 1, kmesh_info%nntot
        nkp2 = kmesh_info%nnlist(nkp_global, nn)
        call zgemm('C', 'N', num_wann, num_wann, num_wann, cmplx_1, u_matrix(:, :, nkp_global), &
                   num_wann, m_matrix_orig_local(:, :, nn, nkp), num_bands, cmplx_0, cwb, num_wann)
        call zgemm('N', 'N', num_wann, num_wann, num_wann, cmplx_1, cwb, num_wann, &
                   u_matrix(:, :, nkp2), num_wann, cmplx_0, cww, num_wann)

        m_matrix_local(1:num_wann, 1:num_wann, nn, nkp) = cww(:, :)
      enddo
    enddo
    ! at this point m_matrix_orig_local may be deassociated/deallocated

    deallocate (cwb)
    deallocate (cww)
    deallocate (global_k)
    if (print_output%timing_level > 1) call io_stopwatch_stop('dis: splitm', timer)
    return
    !================================================!
  end subroutine setup_m_loc

  subroutine internal_check_orthonorm(u_matrix_opt, ndimwin, num_kpts, num_wann, timing_level, &
                                      on_root, timer, error, stdout, comm)
    !================================================!
    !
    !! This subroutine checks that the states in the columns of the
    !! final matrix U_opt are orthonormal at every k-point, i.e.,
    !! that the matrix is unitary in the sense that
    !! conjg(U_opt).U_opt = 1  (but not  U_opt.conjg(U_opt) = 1).
    !!
    !! In particular, this checks whether the projected gaussians
    !! are indeed orthogonal to the frozen states, at those k-points
    !! where both are present in the trial subspace.
    !
    !================================================!
    use w90_comms, only: w90_comm_type
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_constants, only: eps8
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: timing_level
    integer, intent(in) :: stdout
    integer, intent(in) :: num_kpts, num_wann
    integer, intent(in) :: ndimwin(:) ! (num_kpts)

    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)
    logical, intent(in) :: on_root

    ! local variables
    integer          :: nkp, l, m, j

    complex(kind=dp) :: ctmp

    if (timing_level > 1) call io_stopwatch_start('dis: main: check_orthonorm', timer)

    do nkp = 1, num_kpts
      do l = 1, num_wann
        do m = 1, l
          ctmp = cmplx_0
          do j = 1, ndimwin(nkp)
            ctmp = ctmp + conjg(u_matrix_opt(j, m, nkp))*u_matrix_opt(j, l, nkp)
          enddo
          if (l .eq. m) then
            if (abs(ctmp - cmplx_1) .gt. eps8) then
              if (on_root) then
                write (stdout, '(3i6,2f16.12)') nkp, l, m, ctmp
                write (stdout, '(1x,a)') &
                  'The trial orbitals for disentanglement are not orthonormal'
              endif
              call set_error_fatal(error, 'Error in dis_main: orthonormal error 1', comm)
              return
            endif
          else
            if (abs(ctmp) .gt. eps8) then
              if (on_root) then
                write (stdout, '(3i6,2f16.12)') nkp, l, m, ctmp
                write (stdout, '(1x,a)') &
                  'The trial orbitals for disentanglement are not orthonormal'
              endif
              call set_error_fatal(error, 'Error in dis_main: orthonormal error 2', comm)
              return
            endif
          endif
        enddo
      enddo
    enddo

    if (timing_level > 1) call io_stopwatch_stop('dis: main: check_orthonorm', timer)

    return
    !================================================!
  end subroutine internal_check_orthonorm

  subroutine internal_slim_m(m_matrix_orig_local, ndimwin, nfirstwin, nnlist, nntot, num_bands, &
                             timing_level, timer, dist_k, global_k, error, comm)
    !================================================!
    !
    !! This subroutine slims down the original Mmn(k,b), removing
    !! rows and columns corresponding to u_nks that fall outside
    !! the outer energy window.
    !
    !================================================!
    use w90_comms, only: w90_comm_type, mpirank
    use w90_constants, only: dp, cmplx_0
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(w90_comm_type), intent(in) :: comm
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: dist_k(:), global_k(:)
    integer, intent(in) :: ndimwin(:)
    integer, intent(in) :: nfirstwin(:) ! (num_kpts) index of lowest band inside outer window at nkp-th
    integer, intent(in) :: nntot, nnlist(:, :) ! (num_kpts, nntot)
    integer, intent(in) :: num_bands
    integer, intent(in) :: timing_level

    complex(kind=dp), intent(inout) :: m_matrix_orig_local(:, :, :, :)

    ! local variables
    integer :: nkp, nkp2, nn, i, j, m, n, ierr, nkp_global, nkrank, my_node_id

    complex(kind=dp), allocatable :: cmtmp(:, :)

    if (timing_level > 1) call io_stopwatch_start('dis: main: slim_m', timer)

    my_node_id = mpirank(comm)
    nkrank = count(dist_k == my_node_id) ! number of k points this rank

    allocate (cmtmp(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cmtmp in dis_main', comm)
      return
    endif
    do nkp = 1, nkrank
      nkp_global = global_k(nkp)
      do nn = 1, nntot
        nkp2 = nnlist(nkp_global, nn)
        do j = 1, ndimwin(nkp2)
          n = nfirstwin(nkp2) + j - 1
          do i = 1, ndimwin(nkp_global)
            m = nfirstwin(nkp_global) + i - 1
            cmtmp(i, j) = m_matrix_orig_local(m, n, nn, nkp)
          enddo
        enddo
        m_matrix_orig_local(:, :, nn, nkp) = cmplx_0
        do j = 1, ndimwin(nkp2)
          do i = 1, ndimwin(nkp_global)
            m_matrix_orig_local(i, j, nn, nkp) = cmtmp(i, j)
          enddo
        enddo
      enddo
    enddo

    deallocate (cmtmp, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cmtmp in dis_main', comm)
      return
    endif

    if (timing_level > 1) call io_stopwatch_stop('dis: main: slim_m', timer)

    return
    !================================================!
  end subroutine internal_slim_m

  subroutine internal_find_u(sitesym, a_matrix, u_matrix, u_matrix_opt, ndimwin, num_bands, &
                             num_kpts, num_wann, timing_level, lsitesymmetry, on_root, stdout, &
                             timer, error, comm)
    !================================================!
    !
    !! This subroutine finds the initial guess for the square unitary
    !! rotation matrix u_matrix. The method is similar to Sec. III.D
    !! of SMV, but with square instead of rectangular matrices:
    !!
    !! First find caa, the square overlap matrix <psitilde_nk|g_m>,
    !! where psitilde is an eigenstate of the optimal subspace.
    !!
    !! Note that, contrary to what is implied in Sec. III.E of SMV,
    !! this does *not* need to be computed by brute: instead we take
    !! advantage of the previous computation of overlaps with the
    !! same projections that are used to initiate the minimization of
    !! Omega.
    !!
    !! Note: |psi> U_opt = |psitilde> and obviously
    !! <psitilde| = (U_opt)^dagger <psi|
    !
    !================================================!
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_sitesym, only: sitesym_symmetrize_u_matrix
    use w90_types, only: timer_list_type
    use w90_wannier90_types, only: sitesym_type

    implicit none

    ! arguments
    type(sitesym_type), intent(inout) :: sitesym
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: ndimwin(:) ! (num_kpts)
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level

    complex(kind=dp), intent(in) :: a_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    logical, intent(in) :: on_root, lsitesymmetry

    ! local variables
    integer :: nkp, info, ierr
    complex(kind=dp), allocatable :: caa(:, :, :)
    ! For ZGESVD
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp), allocatable :: rwork(:)

    complex(kind=dp), allocatable :: cv(:, :)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cwork(:)

    if (timing_level > 1) call io_stopwatch_start('dis: main: find_u', timer)

    ! Currently, this part is not parallelized; thus, we perform the task only on root and then broadcast the result.
    if (on_root) then
      ! Allocate arrays needed for ZGESVD
      allocate (svals(num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating svals in dis_main', comm)
        return
      endif
      allocate (rwork(5*num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating rwork in dis_main', comm)
        return
      endif
      allocate (cv(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating cv in dis_main', comm)
        return
      endif
      allocate (cz(num_wann, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating cz in dis_main', comm)
        return
      endif
      allocate (cwork(4*num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating cwork in dis_main', comm)
        return
      endif
      allocate (caa(num_wann, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating caa in dis_main', comm)
        return
      endif

      do nkp = 1, num_kpts
        if (lsitesymmetry) then
          if (sitesym%ir2ik(sitesym%ik2ir(nkp)) .ne. nkp) cycle
        endif
        call zgemm('C', 'N', num_wann, num_wann, ndimwin(nkp), cmplx_1, u_matrix_opt(:, :, nkp), &
                   num_bands, a_matrix(:, :, nkp), num_bands, cmplx_0, caa(:, :, nkp), num_wann)
        ! Singular-value decomposition
        call zgesvd('A', 'A', num_wann, num_wann, caa(:, :, nkp), num_wann, svals, cz, num_wann, &
                    cv, num_wann, cwork, 4*num_wann, rwork, info)
        if (info .ne. 0) then
          if (on_root) write (stdout, *) ' ERROR: IN ZGESVD IN dis_main'
          if (on_root) write (stdout, *) 'K-POINT NKP=', nkp, ' INFO=', info
          if (info .lt. 0) then
            if (on_root) write (stdout, *) 'THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
          endif
          call set_error_fatal(error, 'dis_main: problem in ZGESVD 1', comm)
          return
        endif
        ! u_matrix is the initial guess for the unitary rotation of the
        ! basis states given by the subroutine extract
        call zgemm('N', 'N', num_wann, num_wann, num_wann, cmplx_1, cz, num_wann, cv, num_wann, &
                   cmplx_0, u_matrix(:, :, nkp), num_wann)
      enddo
    endif
    call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, error, comm)
    if (allocated(error)) return
!      if (lsitesymmetry) call sitesym_symmetrize_u_matrix(num_wann,u_matrix) !RS:

    if (on_root) then
      ! Deallocate arrays for ZGESVD
      deallocate (caa, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating caa in dis_main', comm)
        return
      endif
      deallocate (cwork, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating cwork in dis_main', comm)
        return
      endif
      deallocate (cz, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating cz in dis_main', comm)
        return
      endif
      deallocate (cv, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating cv in dis_main', comm)
        return
      endif
      deallocate (rwork, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating rwork in dis_main', comm)
        return
      endif
      deallocate (svals, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating svals in dis_main', comm)
        return
      endif
    endif

    if (lsitesymmetry) then
      call sitesym_symmetrize_u_matrix(sitesym, u_matrix, num_bands, num_wann, num_kpts, num_wann, &
                                       stdout, error, comm)
      if (allocated(error)) return
    endif

    if (timing_level > 1) call io_stopwatch_stop('dis: main: find_u', timer)

    return
    !================================================!
  end subroutine internal_find_u

  subroutine internal_find_u_gamma(a_matrix, u_matrix, u_matrix_opt, ndimwin, num_wann, &
                                   timing_level, stdout, timer, error, comm)
    !================================================!
    !
    !! Make initial u_matrix real
    !! Must be the case when gamma_only = .true.
    !
    !================================================!
    use w90_constants, only: dp
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: ndimwin(:)
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level, num_wann

    complex(kind=dp), intent(in) :: a_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    ! local variables
    integer :: info, ierr

    real(kind=dp), allocatable :: a_matrix_r(:, :)
    real(kind=dp), allocatable :: raa(:, :)
    real(kind=dp), allocatable :: u_opt_r(:, :)
    ! for dgesvd
    real(kind=dp), allocatable :: rv(:, :)
    real(kind=dp), allocatable :: rz(:, :)
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp), allocatable :: work(:)

    if (timing_level > 1) call io_stopwatch_start('dis: main: find_u_gamma', timer)

    ! Allocate arrays needed for getting a_matrix_r
    allocate (u_opt_r(ndimwin(1), num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating u_opt_r in dis_main', comm)
      return
    endif
    allocate (a_matrix_r(ndimwin(1), num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating a_matrix_r in dis_main', comm)
      return
    endif

    ! Allocate arrays needed for dgesvd
    allocate (svals(num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating svals in dis_main', comm)
      return
    endif
    allocate (work(5*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating rwork in dis_main', comm)
      return
    endif
    allocate (rv(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cv in dis_main', comm)
      return
    endif
    allocate (rz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cz in dis_main', comm)
      return
    endif
    allocate (raa(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating raa in dis_main', comm)
      return
    endif

    u_opt_r(:, :) = real(u_matrix_opt(1:ndimwin(1), 1:num_wann, 1), dp)

    a_matrix_r(:, :) = real(a_matrix(1:ndimwin(1), 1:num_wann, 1), kind=dp)

    call dgemm('T', 'N', num_wann, num_wann, ndimwin(1), 1.0_dp, u_opt_r, ndimwin(1), a_matrix_r, &
               ndimwin(1), 0.0_dp, raa, num_wann)
    ! Singular-value decomposition
    call dgesvd('A', 'A', num_wann, num_wann, raa, num_wann, svals, rz, num_wann, rv, num_wann, &
                work, 5*num_wann, info)
    if (info .ne. 0) then
      write (stdout, *) ' ERROR: IN DGESVD IN dis_main'
      write (stdout, *) 'K-POINT = Gamma', ' INFO=', info
      if (info .lt. 0) then
        write (stdout, *) 'THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
      endif
      call set_error_fatal(error, 'dis_main: problem in DGESVD 1', comm)
      return
    endif
    ! u_matrix is the initial guess for the unitary rotation of the
    ! basis states given by the subroutine extract
    call dgemm('N', 'N', num_wann, num_wann, num_wann, 1.0_dp, rz, num_wann, rv, num_wann, 0.0_dp, &
               raa, num_wann)

    u_matrix(:, :, 1) = cmplx(raa(:, :), 0.0_dp, dp)

    ! Deallocate arrays for DGESVD
    deallocate (raa, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating raa in dis_main', comm)
      return
    endif
    deallocate (rz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating rz in dis_main', comm)
      return
    endif
    deallocate (rv, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating rv in dis_main', comm)
      return
    endif
    deallocate (work, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating work in dis_main', comm)
      return
    endif
    deallocate (svals, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating svals in dis_main', comm)
      return
    endif

    ! Deallocate arrays for a_matrix_r
    deallocate (a_matrix_r, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating a_matrix_r in dis_main', comm)
      return
    endif
    deallocate (u_opt_r, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating u_opt_r in dis_main', comm)
      return
    endif

    if (timing_level > 1) call io_stopwatch_stop('dis: main: find_u_gamma', timer)

    return
    !================================================!
  end subroutine internal_find_u_gamma

  subroutine dis_windows(dis_spheres, dis_manifold, eigval_opt, kpt_latt, recip_lattice, indxfroz, &
                         indxnfroz, ndimfroz, iprint, num_bands, num_kpts, num_wann, &
                         timing_level, lfrozen, linner, on_root, stdout, timer, error, comm)
    !================================================!
    !
    !! This subroutine selects the states that are inside the outer
    !! window (ie, the energy window out of which we fish out the
    !! optimally-connected subspace) and those that are inside the
    !! inner window (that make up the frozen manifold, and are
    !! straightfowardly included as they are). This, in practice,
    !! amounts to slimming down the original num_wann x num_wann overlap
    !! matrices, removing rows and columns that belong to u_nks that
    !! have been excluded forever, and marking (indexing) the rows and
    !! columns that correspond to frozen states.
    !!
    !! Note - in windows eigval_opt are shifted, so the lowest ones go
    !! from nfirstwin(nkp) to nfirstwin(nkp)+ndimwin(nkp)-1, and above
    !! they are set to zero.
    !
    !================================================!

    ! OUTPUT:
    !     ndimwin(nkp)   number of bands inside outer window at nkp-th k poi
    !     ndimfroz(nkp)  number of frozen bands at nkp-th k point
    !     lfrozen(i,nkp) true if the i-th band inside outer window is frozen
    !     linner         true if there is an inner window
    !     indxfroz(i,nkp) outer-window band index for the i-th frozen state
    !                     (equals 1 if it is the bottom of outer window)
    !     indxnfroz(i,nkp) outer-window band index for the i-th non-frozen s
    !                     (equals 1 if it is the bottom of outer window)
    !     nfirstwin(nkp) index of lowest band inside outer window at nkp-th
    ! MODIFIED:
    !     eigval_opt(nb,nkp) At input it contains a large set of eigenvalues. At
    !                    it is slimmed down to contain only those inside the
    !                    energy window, stored in nb=1,...,ndimwin(nkp)
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: dis_manifold_type, kmesh_info_type, print_output_type, timer_list_type
    use w90_wannier90_types, only: dis_spheres_type

    implicit none

    ! arguments
    type(dis_spheres_type), intent(in) :: dis_spheres
    type(dis_manifold_type), intent(inout) :: dis_manifold ! ndimwin alone is modified
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: iprint, timing_level
    integer, intent(in) :: stdout
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(inout) :: ndimfroz(:)
    integer, intent(inout) :: indxfroz(:, :)
    integer, intent(inout) :: indxnfroz(:, :)

    real(kind=dp), intent(in) :: kpt_latt(3, num_kpts), recip_lattice(3, 3)
    real(kind=dp), intent(inout) :: eigval_opt(:, :)

    logical, intent(in) :: on_root
    logical, intent(inout) :: linner
    logical, intent(inout) :: lfrozen(:, :)

    ! local variables
    integer :: i, j, nkp
    integer :: imin, imax, kifroz_min, kifroz_max
    real(kind=dp) :: dk(3)
    logical :: dis_ok

    if (timing_level > 1) call io_stopwatch_start('dis: windows', timer)

    linner = .false.

    if (iprint > 0) write (stdout, '(1x,a)') &
      '+----------------------------------------------------------------------------+'
    if (iprint > 0) write (stdout, '(1x,a)') &
      '|                              Energy  Windows                               |'
    if (iprint > 0) write (stdout, '(1x,a)') &
      '|                              ---------------                               |'
    if (iprint > 0) write (stdout, '(1x,a,f10.5,a,f10.5,a)') &
      '|                   Outer: ', dis_manifold%win_min, '  to ', dis_manifold%win_max, &
      '  (eV)                   |'
    if (dis_manifold%frozen_states) then
      if (iprint > 0) write (stdout, '(1x,a,f10.5,a,f10.5,a)') &
        '|                   Inner: ', dis_manifold%froz_min, '  to ', dis_manifold%froz_max, &
        '  (eV)                   |'
    else
      if (iprint > 0) write (stdout, '(1x,a)') &
        '|                   No frozen states were specified                          |'
    endif
    if (iprint > 0) write (stdout, '(1x,a)') &
      '+----------------------------------------------------------------------------+'

    do nkp = 1, num_kpts
      ! Check which eigenvalues fall within the outer window
      if ((eigval_opt(1, nkp) .gt. dis_manifold%win_max) .or. &
          (eigval_opt(num_bands, nkp) .lt. dis_manifold%win_min)) then
        if (on_root) then
          write (stdout, *) ' ERROR AT K-POINT: ', nkp
          write (stdout, *) ' ENERGY WINDOW (eV):    [', dis_manifold%win_min, ',', &
            dis_manifold%win_max, ']'
          write (stdout, *) ' EIGENVALUE RANGE (eV): [', &
            eigval_opt(1, nkp), ',', eigval_opt(num_bands, nkp), ']'
          call set_error_fatal(error, 'dis_windows: The outer energy window contains no eigenvalues', comm)
          return
        endif
      endif

      dis_manifold%ndimwin(nkp) = num_bands
      dis_manifold%nfirstwin(nkp) = 1

      ! Note: we assume that eigvals are ordered from the bottom up
      imin = 0
      do i = 1, num_bands
        if (imin .eq. 0) then
          if ((eigval_opt(i, nkp) .ge. dis_manifold%win_min) .and. &
              (eigval_opt(i, nkp) .le. dis_manifold%win_max)) imin = i
          imax = i
        endif
        if (eigval_opt(i, nkp) .le. dis_manifold%win_max) imax = i
      enddo

      dis_manifold%ndimwin(nkp) = imax - imin + 1
      dis_manifold%nfirstwin(nkp) = imin

      !~~ GS-start
      ! disentangle at the current k-point only if it is within one of the
      ! spheres centered at the k-points listed in kpt_dis
      if (dis_spheres%num .gt. 0) then
        dis_ok = .false.
        ! loop on the sphere centers
        do i = 1, dis_spheres%num
          dk = kpt_latt(:, nkp) - dis_spheres%spheres(1:3, i)
          dk = matmul(anint(dk) - dk, recip_lattice(:, :))
          ! if the current k-point is included in at least one sphere,
          ! then perform disentanglement as usual
          if (abs(dot_product(dk, dk)) .lt. dis_spheres%spheres(4, i)**2) then
            dis_ok = .true.
            exit
          endif
        enddo
        ! this kpoint is not included in any sphere: no disentaglement
        if (.not. dis_ok) then
          dis_manifold%ndimwin(nkp) = num_wann
          dis_manifold%nfirstwin(nkp) = dis_spheres%first_wann
        endif
      endif
      !~~ GS-end

      if (dis_manifold%ndimwin(nkp) .lt. num_wann) then
        if (iprint > 0) write (stdout, '(1x, a17, i4, a8, i3, a9, i3)') 'Error at k-point ', nkp, &
          ' ndimwin=', dis_manifold%ndimwin(nkp), ' num_wann=', num_wann
        call set_error_fatal(error, 'dis_windows: Energy window contains fewer states than number of target WFs', comm)
        return
      endif

      do i = 1, dis_manifold%ndimwin(nkp)
        lfrozen(i, nkp) = .false.
      enddo

      ! Check which eigenvalues fall within the inner window
      kifroz_min = 0
      kifroz_max = -1

      ! (note that the above obeys kifroz_max-kifroz_min+1=kdimfroz=0, as we w
      if (dis_manifold%frozen_states) then
        do i = imin, imax
          if (kifroz_min .eq. 0) then
            if ((eigval_opt(i, nkp) .ge. dis_manifold%froz_min) .and. &
                (eigval_opt(i, nkp) .le. dis_manifold%froz_max)) then
              ! Relative to bottom of outer window
              kifroz_min = i - imin + 1
              kifroz_max = i - imin + 1
            endif
          elseif (eigval_opt(i, nkp) .le. dis_manifold%froz_max) then
            kifroz_max = kifroz_max + 1
            ! DEBUG
            ! if(kifroz_max.ne.i-imin+1) stop 'something wrong...'
            ! ENDDEBUG
          endif
        enddo
      endif

      ndimfroz(nkp) = kifroz_max - kifroz_min + 1

      if (ndimfroz(nkp) .gt. num_wann) then
        if (on_root) write (stdout, 401) nkp, ndimfroz(nkp), num_wann
401     format(' ERROR AT K-POINT ', i4, ' THERE ARE ', i2, &
               ' BANDS INSIDE THE INNER WINDOW AND ONLY', i2, &
               ' TARGET BANDS')
        if (on_root) write (stdout, 402) (eigval_opt(i, nkp), i=imin, imax)
402     format('BANDS: (eV)', 10(F10.5, 1X))
        call set_error_fatal(error, 'dis_windows: More states in the frozen window than target WFs', comm)
        return
      endif

      if (ndimfroz(nkp) .gt. 0) linner = .true.
      ! DEBUG
      !         write(*,'(a,i4,a,i2,a,i2)') 'k point ',nkp,
      !     &    ' lowest band in outer win is # ',imin,
      !     &    '   # frozen states is ',ndimfroz(nkp)
      ! ENDDEBUG
      ! Generate index array for frozen states (those inside inner window)
      if (ndimfroz(nkp) .gt. 0) then
        do i = 1, ndimfroz(nkp)
          indxfroz(i, nkp) = kifroz_min + i - 1
          lfrozen(indxfroz(i, nkp), nkp) = .true.
        enddo
        if (indxfroz(ndimfroz(nkp), nkp) .ne. kifroz_max) then
          if (on_root) then
            write (stdout, *) ' Error at k-point ', nkp, ' frozen band #', i
            write (stdout, *) ' ndimfroz=', ndimfroz(nkp)
            write (stdout, *) ' kifroz_min=', kifroz_min
            write (stdout, *) ' kifroz_max=', kifroz_max
            write (stdout, *) ' indxfroz(i,nkp)=', indxfroz(i, nkp)
          endif
          call set_error_fatal(error, 'dis_windows: Something fishy... disentangle.F90+1065', comm)
          return
        endif
      endif

      ! Generate index array for non-frozen states
      i = 0
      do j = 1, dis_manifold%ndimwin(nkp)
        if (.not. lfrozen(j, nkp)) then
          i = i + 1
          indxnfroz(i, nkp) = j
        endif
      enddo

      if (i .ne. dis_manifold%ndimwin(nkp) - ndimfroz(nkp)) then
        if (on_root) write (stdout, *) ' Error at k-point: ', nkp
        if (on_root) write (stdout, '(3(a,i5))') ' i: ', i, ' ndimwin: ', &
          dis_manifold%ndimwin(nkp), ' ndimfroz: ', ndimfroz(nkp)
        call set_error_input(error, 'dis_windows: i .ne. (ndimwin-ndimfroz) at k-point', comm)
        return
      endif

      ! Slim down eigval vector at present k
      do i = 1, dis_manifold%ndimwin(nkp)
        j = dis_manifold%nfirstwin(nkp) + i - 1
        eigval_opt(i, nkp) = eigval_opt(j, nkp)
      enddo

      do i = dis_manifold%ndimwin(nkp) + 1, num_bands
        eigval_opt(i, nkp) = 0.0_dp
      enddo

    enddo
    ! [k-point loop (nkp)]

![ysl-b]
!~    if (gamma_only) then
!~       if (.not. allocated(ph_g)) then
!~          allocate(  ph_g(num_bands),stat=ierr )
!~          if (ierr/=0) call io_error('Error in allocating ph_g in dis_windows')
!~          ph_g = cmplx_1
!~       endif
!~       ! Apply same operation to ph_g
!~       do i = 1, ndimwin(1)
!~          j = nfirstwin(1) + i - 1
!~          ph_g(i) = ph_g(j)
!~       enddo
!~       do i = ndimwin(1) + 1, num_bands
!~          ph_g(i) = cmplx_0
!~       enddo
!~    endif
!~![ysl-e]

    if (iprint > 1) then ! iprint > 0 implies rank > 0
      write (stdout, '(1x,a)') &
        '|                        K-points with Frozen States                         |'
      write (stdout, '(1x,a)') &
        '|                        ---------------------------                         |'
      i = 0
      do nkp = 1, num_kpts
        if (ndimfroz(nkp) .gt. 0) then
          i = i + 1
          if (i .eq. 1) then
            write (stdout, '(1x,a,i6)', advance='no') '|', nkp
          else if ((i .gt. 1) .and. (i .lt. 12)) then
            write (stdout, '(i6)', advance='no') nkp
          else if (i .eq. 12) then
            write (stdout, '(i6,a)') nkp, '    |'
            i = 0
          endif
        endif
      enddo
      if (i .ne. 0) then
        do j = 1, 12 - i
          write (stdout, '(6x)', advance='no')
        enddo
        write (stdout, '(a)') '    |'
      endif
      write (stdout, '(1x,a)') &
        '+----------------------------------------------------------------------------+'
    endif

    if (iprint > 0) write (stdout, '(3x,a,i4)') 'Number of target bands to extract: ', num_wann
    if (iprint > 1) then
      write (stdout, '(4(1x,a,/),(1x,a))') &
        '+----------------------------------------------------------------------------+', &
        '|                                  Windows                                   |', &
        '|                                  -------                                   |', &
        '|               K-point      Ndimwin     Ndimfroz    Nfirstwin               |', &
        '|               ----------------------------------------------               |'

      do nkp = 1, num_kpts
        write (stdout, 403) nkp, dis_manifold%ndimwin(nkp), ndimfroz(nkp), dis_manifold%nfirstwin(nkp)
      enddo
403   format(1x, '|', 14x, i6, 7x, i6, 7x, i6, 6x, i6, 18x, '|')
      write (stdout, '(1x,a)') &
        '+----------------------------------------------------------------------------+'
    endif

    if (timing_level > 1) call io_stopwatch_stop('dis: windows', timer)

    return
    !================================================!
  end subroutine dis_windows

  subroutine dis_project(a_matrix, u_matrix_opt, ndimwin, nfirstwin, num_bands, num_kpts, &
                         num_wann, timing_level, on_root, iprint, timer, error, stdout, comm)
    !================================================!
    !
    !! Construct projections for the start of the disentanglement routine
    !!
    !! Original notes from Nicola (refers only to the square case)
    !!
    !! This subroutine calculates the transformation matrix
    !! CU = CS^(-1/2).CA, where CS = CA.CA^dagger.
    !! CS is diagonalized with a Schur factorization, to be on the safe
    !! side of numerical stability.
    !!
    !! ZGEES computes for an N-by-N complex nonsymmetric matrix Y, the
    !! eigenvalues, the Schur form T, and, optionally, the matrix of
    !! Schur vectors. This gives the Schur factorization Y = Z*T*(Z**H).
    !!
    !! Optionally, it also orders the eigenvalues on the diagonal of the
    !! Schur form so that selected eigenvalues are at the top left.
    !! The leading components of Z then form an orthonormal basis for
    !! the invariant subspace corresponding to the selected eigenvalues.
    !!
    !! A complex matrix is in Schur form if it is upper triangular.
    !!
    !! Notes from Ivo disentangling (e.g. non-square) projection
    !! (See Sec. III.D of SMV paper)
    !! Compute the ndimwin(k) x num_wann matrix cu that yields,
    !! from the ndimwin original Bloch states, the num_wann Bloch-like
    !! states with maximal projection onto the num_wann localised
    !! functions:
    !!
    !! CU = CA.CS^{-1/2}, CS = transpose(CA).CA
    !!
    !! Use the singular-calue decomposition of the matrix CA:
    !!
    !! CA = CZ.CD.CVdagger (note: zgesvd spits out CVdagger)
    !!
    !! which yields
    !!
    !! CU = CZ.CD.CD^{-1}.CVdag
    !!
    !! where CZ is ndimwin(NKP) x ndimwin(NKP) and unitary, CD is
    !! ndimwin(NKP) x num_wann and diagonal, CD^{-1} is
    !! num_wann x num_wann and diagonal, and CVdag is
    !! num_wann x num_wann and unitary.
    !!
    !================================================!
    use w90_constants, only: dp, cmplx_0, cmplx_1, eps5
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: ndimwin(:), nfirstwin(:)
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(in) :: stdout, timing_level, iprint

    complex(kind=dp), intent(inout) :: a_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    logical, intent(in) :: on_root

    ! local variables
    integer :: i, j, l, m, nkp, info, ierr

    real(kind=dp), allocatable :: svals(:)
    real(kind=dp), allocatable :: rwork(:)

    complex(kind=dp) :: ctmp2
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cvdag(:, :)
!    complex(kind=dp), allocatable :: catmpmat(:,:,:)

    if (timing_level > 1) call io_stopwatch_start('dis: project', timer)

    if (iprint > 0) write (stdout, '(/1x,a)') &
      '                  Unitarised projection of Wannier functions                  '
    if (iprint > 0) write (stdout, '(1x,a)') &
      '                  ------------------------------------------                  '
    if (iprint > 0) write (stdout, '(3x,a)') 'A_mn = <psi_m|g_n> --> S = A.A^+ --> U = S^-1/2.A'
    if (iprint > 0) write (stdout, '(3x,a)', advance='no') 'In dis_project...'

    allocate (svals(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating svals in dis_project', comm)
      return
    endif
    allocate (rwork(5*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating rwork in dis_project', comm)
      return
    endif
    allocate (cvdag(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cvdag in dis_project', comm)
      return
    endif
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cz in dis_project', comm)
      return
    endif
    allocate (cwork(4*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cwork in dis_project', comm)
      return
    endif
    ! here we slim down the ca matrix
    ! up to here num_bands(=num_bands) X num_wann(=num_wann)
!    do nkp = 1, num_kpts
!       do j = 1, num_wann
!          do i = 1, ndimwin(nkp)
!             catmpmat(i,j,nkp) = a_matrix(nfirstwin(nkp)+i-1,j,nkp)
!          enddo
!       enddo
!       do j = 1, num_wann
!          a_matrix(1:ndimwin(nkp),j,nkp) = catmpmat(1:ndimwin(nkp),j,nkp)
!       enddo
!       do j = 1, num_wann
!          a_matrix(ndimwin(nkp)+1:num_bands,j,nkp) = cmplx_0
!       enddo
!    enddo
    ! in order to reduce the memory usage, we don't use catmpmat.
    do nkp = 1, num_kpts
      if (ndimwin(nkp) .ne. num_bands) then
        do j = 1, num_wann
          do i = 1, ndimwin(nkp)
            ctmp2 = a_matrix(nfirstwin(nkp) + i - 1, j, nkp)
            a_matrix(i, j, nkp) = ctmp2
          enddo
          a_matrix(ndimwin(nkp) + 1:num_bands, j, nkp) = cmplx_0
        enddo
      endif
    enddo

    do nkp = 1, num_kpts
      ! SINGULAR VALUE DECOMPOSITION
      call zgesvd('A', 'A', ndimwin(nkp), num_wann, a_matrix(:, :, nkp), num_bands, svals, cz, &
                  num_bands, cvdag, num_bands, cwork, 4*num_bands, rwork, info)
      if (info .ne. 0) then
        if (on_root) write (stdout, *) ' ERROR: IN ZGESVD IN dis_project'
        if (on_root) write (stdout, *) ' K-POINT NKP=', nkp, ' INFO=', info
        if (info .lt. 0) then
          if (on_root) write (stdout, *) ' THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
        endif
        call set_error_fatal(error, 'dis_project: problem in ZGESVD 1', comm)
        return
      endif

      ! NOTE THAT - AT LEAST FOR LINUX MKL LAPACK - THE OUTPUT OF ZGESVD
      ! GIVES ALREADY Vdagger, SO A=Z.S.Vdagger IS ACTUALLY GIVEN BY cz.s.cvda
      !
      ! also, cu is cz.cd.cd^-1.cvdag, and the asymmetric structure of cd.cd^-
      ! allows us to say cu=cz.cvdag, but where the sum on the inner index
      ! goes only from 1 to num_wann (i.e. DO l=1,num_wann). This is because
      ! cz.cd.cd^-1 is a matrix ndimwin(nkp) X num_wann, identical to the
      ! first num_wann columns of the ndimwin(nkp) X ndimwin(nkp) matrix cz
      !
      ! same for ca: s is a ndimwin(nkp) X num_wann matrix that is
      ! zero everywhere but for its num_wann X num_wann top square part,
      ! that is diagonal. Multiplying cz by the s matrix is equivalent
      ! to moltiplying the first num_wann columns of cz, each by the correspondin
      ! diagonal element of s, that is s(L)
      ! I'm not sure why we reconstruct ca in what follows - in one explicit t
      ! [ aam: it is because a_matrix is overwritten by ZGESVD ]
      ! it seemed to be identical to the input ca (as it should be)

      u_matrix_opt(:, :, nkp) = cmplx_0
      a_matrix(:, :, nkp) = cmplx_0
      do j = 1, num_wann
        do i = 1, ndimwin(nkp)
          do l = 1, num_wann
            u_matrix_opt(i, j, nkp) = u_matrix_opt(i, j, nkp) + cz(i, l)*cvdag(l, j)
            a_matrix(i, j, nkp) = a_matrix(i, j, nkp) + cz(i, l)*svals(l)*cvdag(l, j)
          enddo
        enddo
      enddo

      !
      ! CHECK UNITARITY
      !
      ! note that cu.transpose(cu) is *NOT* an identity ndimwin(nkp) by ndimwi
      ! matrix, but transpose(cu).cu is a num_wann by num_wann identity matrix.
      ! I have once checked the former statement, now I will just leave here t
      ! for the latter (what this means is that the columns of cu are orthonor
      ! vectors).
      do i = 1, num_wann
        do j = 1, num_wann
          ctmp2 = cmplx_0
          do m = 1, ndimwin(nkp)
            ctmp2 = ctmp2 + u_matrix_opt(m, j, nkp)*conjg(u_matrix_opt(m, i, nkp))
          enddo
          if ((i .eq. j) .and. (abs(ctmp2 - cmplx_1) .gt. eps5)) then
            if (on_root) write (stdout, *) ' ERROR: unitarity of initial U'
            if (on_root) write (stdout, '(1x,a,i2)') 'nkp= ', nkp
            if (on_root) write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
            if (on_root) write (stdout, '(1x,a,f12.6,1x,f12.6)') &
              '[u_matrix_opt.transpose(u_matrix_opt)]_ij= ', &
              real(ctmp2, dp), aimag(ctmp2)
            call set_error_fatal(error, 'dis_project: Error in unitarity of initial U in dis_project', comm)
            return
          endif
          if ((i .ne. j) .and. (abs(ctmp2) .gt. eps5)) then
            if (on_root) write (stdout, *) ' ERROR: unitarity of initial U'
            if (on_root) write (stdout, '(1x,a,i2)') 'nkp= ', nkp
            if (on_root) write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
            if (on_root) write (stdout, '(1x,a,f12.6,1x,f12.6)') &
              '[u_matrix_opt.transpose(u_matrix_opt)]_ij= ', &
              real(ctmp2, dp), aimag(ctmp2)
            call set_error_fatal(error, 'dis_project: Error in unitarity of initial U in dis_project', comm)
            return
          endif
        enddo
      enddo
    enddo
    ! NKP

    deallocate (cwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cwork in dis_project', comm)
      return
    endif
    deallocate (cz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cz in dis_project', comm)
      return
    endif
    deallocate (cvdag, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cvdag in dis_project', comm)
      return
    endif
    deallocate (rwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating rwork in dis_project', comm)
      return
    endif
    deallocate (svals, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating svals in dis_project', comm)
      return
    endif

    if (iprint > 0) write (stdout, '(a)') ' done'

    if (timing_level > 1) call io_stopwatch_stop('dis: project', timer)

    return
    !================================================!
  end subroutine dis_project

  subroutine dis_proj_froz(u_matrix_opt, indxfroz, ndimfroz, ndimwin, iprint, num_bands, num_kpts, &
                           num_wann, timing_level, lfrozen, on_root, timer, error, stdout, comm)
    !================================================!
    !
    !! COMPUTES THE LEADING EIGENVECTORS OF Q_froz . P_s . Q_froz,
    !! WHERE P_s PROJECTOR OPERATOR ONTO THE SUBSPACE S OF THE PROJECTED
    !! GAUSSIANS, P_f THE PROJECTOR ONTO THE FROZEN STATES, AND
    !! Q_froz = 1 - P_froz, ALL EXP IN THE BASIS OF THE BLOCH
    !! EIGENSTATES INSIDE THE OUTER ENERGY WINDOW
    !! (See Eq. (27) in Sec. III.G of SMV)
    !
    !================================================!
    use w90_constants, only: dp, cmplx_0, cmplx_1, eps8
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: timing_level, iprint
    integer, intent(in) :: stdout
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(in) :: ndimwin(:)
    integer, intent(in) :: ndimfroz(:)
    integer, intent(in) :: indxfroz(:, :)

    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    logical, intent(in) :: on_root, lfrozen(:, :)

    ! INPUT: num_wann,ndimwin,ndimfroz,indxfroz,lfrozen
    ! MODIFIED: u_matrix_opt (At input it contains the gaussians projected onto
    !             the window states in the routine project.f. At output
    !             the entries with the second index from 1 to ndimfroz(nkp)
    !             contain the frozen (inner window) states, while those
    !             from ndimfroz(nkp)+1 to num_wann have been replaced by
    !             the new trial states outside the inner window.)

    ! *********************************************************
    ! VARIABLES USED BY LAPACK'S ZHPEVX DIAGONALIZATION ROUTINE
    ! *********************************************************
    integer, allocatable :: iwork(:)
    integer, allocatable :: ifail(:)

    real(kind=dp), allocatable :: w(:)
    real(kind=dp), allocatable :: rwork(:)

    complex(kind=dp), allocatable :: cap(:)
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)

    ! *********
    ! INTERNAL:
    ! *********
    !
    ! CP_S(M,N)      PROJECTION OPERATOR ONTO THE SUBSPACE OF THE PROJEC
    !                  GAUSSIANS, EXPRESSED IN THE BASIS OF THE ORIGINAL BL
    !                  EIGENSTATES INSIDE THE OUTER WINDOW (FOR THE PRESENT
    !                  K-POINT)
    ! CQ_FROZ(M,N)   PROJECTION OPERATOR ONTO THE SUBSPACE OF THE STATES
    !                  THE SPACE OF FROZEN STATES (BUT INSIDE THE OUTER WIN
    !                  EXPRESSED IN THE BASIS OF THE ORIGINAL BLOCH EIGENST
    !                  INSIDE THE OUTER WINDOW (FOR THE PRESENT K-POINT)
    ! CPQ(M,N)       THE MATRIX cp_s . cq_froz FOR THE PRESENT K-POINT
    ! CQPQ(M,N)      THE MATRIX cq_froz . cp_s . cq_froz FOR THE PRESENT
    !

    integer :: goods, il, iu, nkp, l, j, n, m, info, ierr
    integer :: counter, loop_f, loop_v, vmap(num_bands)
    integer :: nzero

    complex(kind=dp) :: ctmp
    complex(kind=dp), allocatable :: cp_s(:, :)
    complex(kind=dp), allocatable :: cq_froz(:, :)
    complex(kind=dp), allocatable :: cpq(:, :)
    complex(kind=dp), allocatable :: cqpq(:, :)

    logical :: take

    character(len=4) :: rep

    if (timing_level > 1) call io_stopwatch_start('dis: proj_froz', timer)

    if (iprint > 0) write (stdout, '(3x,a)', advance='no') 'In dis_proj_froz...'

    allocate (iwork(5*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating iwork in dis_proj_froz', comm)
      return
    endif
    allocate (ifail(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating ifail in dis_proj_froz', comm)
      return
    endif
    allocate (w(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating w in dis_proj_froz', comm)
      return
    endif
    allocate (rwork(7*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating rwork in dis_proj_froz', comm)
      return
    endif
    allocate (cap((num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cap in dis_proj_froz', comm)
      return
    endif
    allocate (cwork(2*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cwork in dis_proj_froz', comm)
      return
    endif
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cz in dis_proj_froz', comm)
      return
    endif
    allocate (cp_s(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cp_s in dis_proj_froz', comm)
      return
    endif
    allocate (cq_froz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cq_froz in dis_proj_froz', comm)
      return
    endif
    allocate (cpq(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cpq in dis_proj_froz', comm)
      return
    endif
    allocate (cqpq(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cqpq in dis_proj_froz', comm)
      return
    endif

    do nkp = 1, num_kpts

      ! aam: this should be done at the end, otherwise important
      !      projection info is lost
!~         ! Put the frozen states in the lowest columns of u_matrix_opt
!~         if (ndimfroz(nkp).gt.0) then
!~            do l = 1, ndimfroz(nkp)
!~               u_matrix_opt(:,l,nkp)=cmplx_0
!~               u_matrix_opt(indxfroz(l,nkp),l,nkp) = cmplx_1
!~            enddo
!~         endif

      ! If there are non-frozen states, compute the num_wann-ndimfroz(nkp) leadin
      ! eigenvectors of cqpq
      if (num_wann .gt. ndimfroz(nkp)) then
        cq_froz = cmplx_0
        cp_s = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, ndimwin(nkp)
            do l = 1, num_wann
              cp_s(m, n) = cp_s(m, n) + u_matrix_opt(m, l, nkp)*conjg(u_matrix_opt(n, l, nkp))
            enddo
          enddo
          if (.not. lfrozen(n, nkp)) cq_froz(n, n) = cmplx_1
        enddo

        cpq = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, ndimwin(nkp)
            do l = 1, ndimwin(nkp)
              cpq(m, n) = cpq(m, n) + cp_s(m, l)*cq_froz(l, n)
            enddo
          enddo
        enddo

        cqpq = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, ndimwin(nkp)
            do l = 1, ndimwin(nkp)
              cqpq(m, n) = cqpq(m, n) + cq_froz(m, l)*cpq(l, n)
            enddo
          enddo
        enddo

        ! DEBUG
        ! check hermiticity of cqpq
        do n = 1, ndimwin(nkp)
          do m = 1, n
            if (abs(cqpq(m, n) - conjg(cqpq(n, m))) .gt. eps8) then
              if (on_root) write (stdout, *) ' matrix CQPQ is not hermitian'
              if (on_root) write (stdout, *) ' k-point ', nkp
              call set_error_fatal(error, 'dis_proj_froz: error', comm)
              return
            endif
          enddo
        enddo
        ! ENDDEBUG

        cap = cmplx_0
        do n = 1, ndimwin(nkp)
          do m = 1, n
            cap(m + (n - 1)*n/2) = cqpq(m, n)
          enddo
        enddo
        il = ndimwin(nkp) - (num_wann - ndimfroz(nkp)) + 1
        iu = ndimwin(nkp)
        call zhpevx('V', 'A', 'U', ndimwin(nkp), cap, 0.0_dp, 0.0_dp, il, iu, -1.0_dp, m, w, cz, &
                    num_bands, cwork, rwork, iwork, ifail, info)

!~            write(stdout,*) 'w:'
!~            do n=1,ndimwin(nkp)
!~               write(stdout,'(f14.10)') w(n)
!~            enddo
!~            write(stdout,*) 'cz:'
!~            do n=1,ndimwin(nkp)
!~               write(stdout,'(6f12.8)') cz(n,il), cz(n,iu)
!~            enddo

        ! DEBUG
        if (info .lt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING CQPQ MATRIX'
          if (on_root) write (stdout, *) ' THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
          call set_error_fatal(error, 'dis_proj_frozen: error', comm)
          return
        elseif (info .gt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING CQPQ MATRIX'
          if (on_root) write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
          call set_error_fatal(error, 'dis_proj_frozen: error', comm)
          return
        endif
        ! ENDDEBUG

        ! DEBUG
        if (m .ne. ndimwin(nkp)) then
          if (on_root) write (stdout, *) ' *** ERROR *** in dis_proj_froz'
          if (on_root) write (stdout, *) ' Number of eigenvalues/vectors obtained is', &
            m, ' not equal to the number asked,', ndimwin(nkp)
          call set_error_fatal(error, 'dis_proj_frozen: error', comm)
          return
        endif
        ! ENDDEBUG

        ! DEBUG
        ! check that the eigenvalues are between 0 and 1
        if (iprint > 2) then
          if (on_root) write (stdout, '(/a,i3,a,i3,a,i3,a)') ' K-point ', nkp, ' ndimwin: ', &
            ndimwin(nkp), ' we want the', num_wann - ndimfroz(nkp), &
            ' leading eigenvector(s) of QPQ'
        endif
        do j = 1, ndimwin(nkp)
          if (iprint > 2 .and. on_root) write (stdout, '(a,i3,a,f16.12)') '  lambda(', j, ')=', w(j)
!~[aam]        if ( (w(j).lt.eps8).or.(w(j).gt.1.0_dp + eps8) ) then
          if ((w(j) .lt. -eps8) .or. (w(j) .gt. 1.0_dp + eps8)) then
            call set_error_fatal(error, 'dis_proj_frozen: error - Eigenvalues not between 0 and 1', comm)
            return
          endif
        enddo
        ! ENDDEBUG

        ! [ aam: sometimes the leading eigenvalues form a degenerate set that is
        !        of higher dimensionality than (num_wann - ndimfroz). May need to
        !        fix this at some point. ]

        ! For certain cases we have found that one of the required eigenvectors of cqpq
        ! has a zero eigenvalue (ie it forms a degenerate set with the frozen states).
        ! It depends on floating point math as whether this eigenvalue is positive
        ! or negative (ie +/- 1e-17). If it's positive everything is ok. If negative we
        ! can end up putting one of the frozen states into U_opt (and failing the later
        ! orthogonality check).
        ! This fix detects this situation. If applies we choose the eigenvectors by
        ! checking their orthogonality to the frozen states.
        ! === For version 1.0.1 we make this the default ===

        !if (index(devel_flag, 'no-orth-fix') == 0) then
        nzero = 0; goods = 0
        do j = ndimwin(nkp), ndimwin(nkp) - (num_wann - ndimfroz(nkp)) + 1, -1
          if (w(j) < eps8) then
            nzero = nzero + 1
          else
            goods = goods + 1
          end if
        end do
        if (nzero > 0) then
          if (iprint > 2 .and. on_root) then
            write (stdout, *) ' '
            write (stdout, '(1x,a,i0,a)') 'An eigenvalue of QPQ is close to zero at kpoint ', &
              nkp, '. Using safety check.'
            write (stdout, '(1x,a,i4,a,i4)') 'We must find ', nzero, &
              ' eigenvectors with zero eigenvalues out of a set of ', ndimwin(nkp) - goods
          endif
          !First lets put the 'good' states into vamp
          vmap = 0
          counter = 1
          do j = ndimwin(nkp), ndimwin(nkp) - goods + 1, -1
            vmap(counter) = j
            counter = counter + 1
          end do

          if (iprint > 2 .and. on_root) then
            do loop_f = 1, ndimwin(nkp)
              write (stdout, '(1x,a,i4,a,es13.6)') 'Eigenvector number', loop_f, &
                '    Eigenvalue: ', w(loop_f)
              do loop_v = 1, ndimwin(nkp)
                write (stdout, '(20x,2f12.8)') cz(loop_v, loop_f)
              end do
              write (stdout, *)
            end do
          end if

          ! We need to find nzero vectors out of the remining ndimwin(nkp)-goods vectors

          do loop_f = 1, nzero
            do loop_v = ndimwin(nkp), 1, -1 !loop backwards for efficiency only
              if (any(vmap == loop_v)) cycle
              !check to see if vector is orthogonal to frozen states in u_matrix_opt
              take = .true.
              do m = 1, ndimfroz(nkp)
                ctmp = cmplx_0
                do j = 1, ndimwin(nkp)
                  ctmp = ctmp + conjg(u_matrix_opt(j, m, nkp))*cz(j, loop_v)
                enddo
                if (abs(ctmp) .gt. eps8) then
                  take = .false.
                endif
              enddo
              if (take) then !vector is good and we add it to vmap
                vmap(goods + loop_f) = loop_v
                exit
              end if
            end do
          end do

          if (iprint > 2 .and. on_root) then
            write (rep, '(i4)') num_wann - ndimfroz(nkp)
            write (stdout, '(1x,a,'//trim(rep)//'(i0,1x))') &
              'We use the following eigenvectors: ', vmap(1:(num_wann - ndimfroz(nkp)))
          end if
          do l = 1, num_wann - ndimfroz(nkp)
            if (vmap(l) == 0) then
              call set_error_fatal(error, 'dis_proj_froz: Ortho-fix failed to find enough vectors', comm)
              return
            endif
          end do

          ! put the correct eigenvectors into u_matrix_opt, and we're all done!
          counter = 1
          do l = ndimfroz(nkp) + 1, num_wann
            u_matrix_opt(1:ndimwin(nkp), l, nkp) = cz(1:ndimwin(nkp), vmap(counter))
            counter = counter + 1
          enddo

        else ! we don't need to use the fix

          do l = ndimfroz(nkp) + 1, num_wann
            u_matrix_opt(1:ndimwin(nkp), l, nkp) = cz(1:ndimwin(nkp), il)
            il = il + 1
          enddo

          if (il - 1 .ne. iu) then
            call set_error_fatal(error, 'dis_proj_frozen: error -  il-1.ne.iu  (in ortho-fix)', comm)
            return
          endif

        end if

        ! PICK THE num_wann-nDIMFROZ(NKP) LEADING EIGENVECTORS AS TRIAL STATES
        ! and PUT THEM RIGHT AFTER THE FROZEN STATES IN u_matrix_opt
        !do l = ndimfroz(nkp) + 1, num_wann
        !  if (on_root) write (stdout, *) 'il=', il
        !  u_matrix_opt(1:ndimwin(nkp), l, nkp) = cz(1:ndimwin(nkp), il)
        !  il = il + 1
        !enddo

        ! DEBUG
        !if (il - 1 .ne. iu) then
        !  call io_error('dis_proj_frozen: error -  il-1.ne.iu', stdout, seedname)
        !endif
        ! ENDDEBUG

        !end if

      endif   ! num_wann>nDIMFROZ(NKP)

      ! Put the frozen states in the lowest columns of u_matrix_opt
      if (ndimfroz(nkp) .gt. 0) then
        do l = 1, ndimfroz(nkp)
          u_matrix_opt(:, l, nkp) = cmplx_0
          u_matrix_opt(indxfroz(l, nkp), l, nkp) = cmplx_1
        enddo
      endif
    enddo   ! NKP

    deallocate (cqpq, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cqpq in dis_proj_froz', comm)
      return
    endif
    deallocate (cpq, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cpq in dis_proj_froz', comm)
      return
    endif
    deallocate (cq_froz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cq_froz in dis_proj_froz', comm)
      return
    endif
    deallocate (cp_s, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cp_s in dis_proj_froz', comm)
      return
    endif
    deallocate (cz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cz in dis_proj_froz', comm)
      return
    endif
    deallocate (cwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cwork in dis_proj_froz', comm)
      return
    endif
    deallocate (cap, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cap in dis_project', comm)
      return
    endif
    deallocate (rwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating rwork in dis_proj_froz', comm)
      return
    endif
    deallocate (w, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating w in dis_proj_froz', comm)
      return
    endif
    deallocate (ifail, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating ifail in dis_proj_froz', comm)
      return
    endif
    deallocate (iwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating iwork in dis_proj_froz', comm)
      return
    endif

    if (iprint > 0) write (stdout, '(a)') ' done'

    if (timing_level > 1) call io_stopwatch_stop('dis: proj_froz', timer)

    return
    !================================================!
  end subroutine dis_proj_froz

  subroutine dis_extract(dis_control, kmesh_info, sitesym, print_output, dis_manifold, &
                         m_matrix_orig_local, u_matrix_opt, eigval_opt, omega_invariant, &
                         indxnfroz, ndimfroz, my_node_id, num_bands, num_kpts, num_wann, &
                         lsitesymmetry, timer, ranknk, global_k, error, stdout, comm)
    !================================================!
    !
    !! Extracts an num_wann-dimensional subspace at each k by
    !! minimizing Omega_I
    !
    !================================================!

    ! MODIFIED:
    !           u_matrix_opt (At input it contains the initial guess for the optima
    ! subspace (expressed in terms of the original states inside the window)
    ! output it contains the  states that diagonalize the hamiltonian inside
    ! optimal subspace (again expressed in terms of the original window stat
    ! Giving out states that diagonalize the hamiltonian inside the optimal
    ! subspace (instead of the eigenstates of the Z matrix) is useful for
    ! performing the Wannier interpolation of the energy bands as described
    ! Sec. III.F of SMV)
    !
    !           eigval (At input: original energy eigenvalues.
    ! At output: eigenvalues of the hamiltonian inside optimal subspace)

    ! ----------------------------------------------------------------------
    ! TO DO: The complement subspace is computed but is not saved anywhere!
    ! (Check what was done with it in original code space.f)
    ! Diagonalize Z matrix only at those k points where ndimwin>num_wann?
    ! ----------------------------------------------------------------------

    ! *******************
    ! SHELLS OF K-VECTORS
    ! *******************
    ! nshells           number of shells of k-points to be used in the
    !                   finite-difference formulas for the k-derivatives
    ! aam: wb is now wb(1:nntot) 09/04/2006
    ! wb(nkp,nnx)       weight of the nnx-th b-vector (ordered along shells
    !                   of increasing length) associated with the nkp-th k-p
    ! wbtot             sum of the weights of all b-vectors associated with
    !                   given k-point (k-point 1 is used in calculation)
    ! nnlist(nkp,nnx)   vkpt(1:3,nnlist(nkp,nnx)) is the nnx-th neighboring
    !                   k-point of the nkp-th k-point vkpt(1:3,nkp) (or its
    !                   periodic image in the "home Brillouin zone")
    ! cm(n,m,nkp,nnx)   Overlap matrix <u_nk|u_{m,k+b}>
    use w90_comms, only: comms_bcast, comms_allreduce, w90_comm_type, mpirank
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_error
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop, io_wallclocktime
    use w90_sitesym, only: sitesym_symmetrize_u_matrix, sitesym_symmetrize_zmatrix, &
      sitesym_dis_extract_symmetry
    use w90_types, only: dis_manifold_type, kmesh_info_type, print_output_type, timer_list_type
    use w90_wannier90_types, only: dis_control_type, dis_spheres_type, sitesym_type

    implicit none

    ! arguments
    type(dis_control_type), intent(in) :: dis_control
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(sitesym_type), intent(in) :: sitesym
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: my_node_id
    integer, intent(in) :: stdout
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(in) :: ndimfroz(:) ! (num_kpts)
    integer, intent(in) :: indxnfroz(:, :) ! (num_bands,num_kpts)
    integer, intent(in) :: ranknk, global_k(:)

    real(kind=dp), intent(inout) :: eigval_opt(:, :)
    real(kind=dp), intent(out) :: omega_invariant

    complex(kind=dp), intent(inout) :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    logical, intent(in) :: lsitesymmetry

    ! Internal variables
    integer :: i, j, l, m, n, nn, nkp, nkp2, info, ierr, ndimk, p
    integer :: icompflag, iter, ndiff
    real(kind=dp) :: womegai, wkomegai, womegai1, rsum, delta_womegai
    real(kind=dp), allocatable :: wkomegai1(:)
    real(kind=dp), allocatable :: history(:)
    real(kind=dp), allocatable :: rwork(:)
    real(kind=dp), allocatable :: w(:)
    real(kind=dp), allocatable :: wkomegai1_loc(:)

    complex(kind=dp), allocatable :: camp_loc(:, :, :)
    complex(kind=dp), allocatable :: u_matrix_opt_loc(:, :, :)
    complex(kind=dp), allocatable :: ceamp(:, :, :)
    complex(kind=dp), allocatable :: camp(:, :, :)
    complex(kind=dp), allocatable :: czmat_in(:, :, :)
    complex(kind=dp), allocatable :: czmat_out(:, :, :)
    ! the z-matrices are now stored in local arrays
    complex(kind=dp), allocatable :: czmat_in_loc(:, :, :)
    complex(kind=dp), allocatable :: czmat_out_loc(:, :, :)
    complex(kind=dp), allocatable :: cham(:, :, :)

    complex(kind=dp), allocatable :: cap(:)
    complex(kind=dp), allocatable :: cwb(:, :), cww(:, :), cbw(:, :)
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp) :: lambda(num_wann, num_wann) !RS:

    integer, allocatable :: ifail(:)
    integer, allocatable :: iwork(:)
    integer :: nkp_loc

    logical :: dis_converged
    logical :: on_root = .false.

    on_root = (my_node_id == 0)

    if (print_output%timing_level > 1) call io_stopwatch_start('dis: extract', timer)

    if (print_output%iprint > 0) write (stdout, '(/1x,a)') &
      '                  Extraction of optimally-connected subspace                  '
    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      '                  ------------------------------------------                  '

    allocate (cwb(num_wann, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cwb in dis_extract', comm)
      return
    endif
    allocate (cww(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cww in dis_extract', comm)
      return
    endif
    allocate (cbw(num_bands, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cbw in dis_extract', comm)
      return
    endif
    allocate (iwork(5*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating iwork in dis_extract', comm)
      return
    endif
    allocate (ifail(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating ifail in dis_extract', comm)
      return
    endif
    allocate (w(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating w in dis_extract', comm)
      return
    endif
    allocate (rwork(7*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating rwork in dis_extract', comm)
      return
    endif
    allocate (cap((num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cap in dis_extract', comm)
      return
    endif
    allocate (cwork(2*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cwork in dis_extract', comm)
      return
    endif
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cz in dis_extract', comm)
      return
    endif
    allocate (u_matrix_opt_loc(num_bands, num_wann, ranknk), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating u_matrix_opt_loc in dis_extract', comm)
      return
    endif
    allocate (wkomegai1_loc(ranknk), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wkomegai1_loc in dis_extract', comm)
      return
    endif
    allocate (czmat_in_loc(num_bands, num_bands, ranknk), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating czmat_in_loc in dis_extract', comm)
      return
    endif
    allocate (czmat_out_loc(num_bands, num_bands, ranknk), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating czmat_out_loc in dis_extract', comm)
      return
    endif
    allocate (wkomegai1(num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wkomegai1 in dis_extract', comm)
      return
    endif
    allocate (czmat_in(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating czmat_in in dis_extract', comm)
      return
    endif
    allocate (czmat_out(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating czmat_out in dis_extract', comm)
      return
    endif
    allocate (history(dis_control%conv_window), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating history in dis_extract', comm)
      return
    endif

    cwb = cmplx_0; cww = cmplx_0; cbw = cmplx_0

    ! Copy matrix elements from global U matrix to local U matrix
    do nkp_loc = 1, ranknk
      nkp = global_k(nkp_loc)
      u_matrix_opt_loc(:, :, nkp_loc) = u_matrix_opt(:, :, nkp)
    enddo

    ! ********************************************
    ! ENERGY WINDOWS AND SUBSPACES AT EACH K-POINT
    ! ********************************************
    ! num_wann             dimensionality of the subspace at each k-point
    !                   (number of Wannier functions per unit cell that we w
    ! NDIMWIN(NKP)      number of bands at the nkp-th k-point that fall
    !                   within the outer energy window
    ! NDIMFROZ(NKP)     number of frozen bands at the nkp-th k-point
    ! INDXNFROZ(I,NKP)  INDEX (BETWEEN 1 AND NDIMWIN(NKP)) OF THE I-TH NON-F
    !                   ORIGINAL BAND STATE AT THE NKP-TH K-POINT
    ! U_MATRIX_OPT(J,L,NKP)    AMPLITUDE OF THE J-TH ENERGY EIGENVECTOR INSIDE THE
    !                   ENERGY WINDOW AT THE NKP-TH K-POINT IN THE EXPANSION
    !                   THE L-TH LEADING RLAMBDA EIGENVECTOR AT THE SAME K-P
    !                   If there are M_k frozen states, they occupy the lowe
    !                   entries of the second index of u_matrix_opt, and the leadin
    !                   nabnds-M_k eigenvectors of the Z matrix occupy the
    !                   remaining slots
    ! CAMP(J,L,NKP)     SAME AS U_MATRIX_OPT, BUT FOR THE COMPLEMENT SUBSPACE INSID
    !                   ENERGY WINDOW (I.E., THE NON-LEADING RLAMBDA EIGENVE
    ! CEAMP(J,L,NKPTS)  SAME AS U_MATRIX_OPT, BUT INSTEAD OF RLAMBDA EIGENVECTOR, I
    !                   FOR THE ENERGY EIGENVECTOR OBTAINED BY DIAGONALIZING
    !                   HAMILTONIAN IN THE OPTIMIZED SUBSPACE
    ! CZMAT_IN(M,N,NKP) Z-MATRIX [Eq. (21) SMV]
    ! CZMAT_OUT(M,N,NKP) OUTPUT Z-MATRIX FROM THE PRESENT ITERATION
    ! RLAMBDA           An eigenvalue of the Z matrix
    ! womegai           Gauge-invariant Wannier spread, computed usinf all s
    !                   from current iteration
    ! wkomegai1(NKP)    Eq. (18) of SMV
    ! womegai1          Eq.(11) of SMV (like wowmegai, but neighboring state
    !                   for computing overlaps are from previous iteration.
    !                   become equal at self-consistency)
    ! alphafixe         mixing parameter for the iterative procedure
    ! nitere            total number of iterations

    ! DEBUG
    if (print_output%iprint > 2) then
      if (on_root) then
        write (stdout, '(a,/)') '  Original eigenvalues inside outer window:'
        do nkp = 1, num_kpts
          write (stdout, '(a,i3,3x,20(f9.5,1x))') '  K-point ', nkp, &
            (eigval_opt(i, nkp), i=1, dis_manifold%ndimwin(nkp))
        enddo
      endif
    endif
    ! ENDDEBUG

    ! TO DO: Check if this is the best place to initialize icompflag
    icompflag = 0

    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'
    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      '|  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS'
    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'

    dis_converged = .false.

    ! ------------------
    ! BIG ITERATION LOOP
    ! ------------------
    do iter = 1, dis_control%num_iter

      if (print_output%timing_level > 1) call io_stopwatch_start('dis: extract_1', timer)

      if (iter .eq. 1) then
        ! Initialize Z matrix at k points w/ non-frozen states
        do nkp_loc = 1, ranknk
          nkp = global_k(nkp_loc)
          if (num_wann .gt. ndimfroz(nkp)) then
            call internal_zmatrix(cbw, czmat_in_loc(:, :, nkp_loc), m_matrix_orig_local, &
                                  u_matrix_opt, kmesh_info%wb, indxnfroz, ndimfroz, &
                                  dis_manifold%ndimwin, kmesh_info%nnlist, nkp, nkp_loc, &
                                  kmesh_info%nntot, num_bands, num_wann, print_output%timing_level, &
                                  timer)
          endif
        enddo

        if (lsitesymmetry) then
          czmat_in = cmplx_0
          do nkp_loc = 1, ranknk
            nkp = global_k(nkp_loc)
            czmat_in(:, :, nkp) = czmat_in_loc(:, :, nkp_loc)
          enddo
          call comms_allreduce(czmat_in(1, 1, 1), num_bands*num_bands*num_kpts, 'SUM', error, comm)
          if (allocated(error)) return
          call sitesym_symmetrize_zmatrix(sitesym, czmat_in, num_bands, num_kpts, &
                                          dis_manifold%lwindow)
          do nkp_loc = 1, ranknk
            nkp = global_k(nkp_loc)
            czmat_in_loc(:, :, nkp_loc) = czmat_in(:, :, nkp)
          end do
        end if

      else
        ! [iter.ne.1]
        ! Update Z matrix at k points with non-frozen states, using a mixing sch
        do nkp_loc = 1, ranknk
          nkp = global_k(nkp_loc)
          if (lsitesymmetry) then
            if (sitesym%ir2ik(sitesym%ik2ir(nkp)) .ne. nkp) cycle
          endif
          if (num_wann .gt. ndimfroz(nkp)) then
            ndimk = dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, ndimk
              do j = 1, i
                czmat_in_loc(j, i, nkp_loc) = &
                  cmplx(dis_control%mix_ratio, 0.0_dp, dp)*czmat_out_loc(j, i, nkp_loc) &
                  + cmplx(1.0_dp - dis_control%mix_ratio, 0.0_dp, dp)*czmat_in_loc(j, i, nkp_loc)
                ! hermiticity
                czmat_in_loc(i, j, nkp_loc) = conjg(czmat_in_loc(j, i, nkp_loc))
              enddo
            enddo
          endif
        enddo
      endif

      if (print_output%timing_level > 1) call io_stopwatch_stop('dis: extract_1', timer)
      if (print_output%timing_level > 1) call io_stopwatch_start('dis: extract_2', timer)

      womegai1 = 0.0_dp
      ! wkomegai1 is defined by Eq. (18) of SMV.
      ! Contribution to wkomegai1 from frozen states should be calculated now
      ! every k (before updating any k), so that for iter>1 overlaps are with
      ! non-frozen neighboring states from the previous iteration
      wkomegai1 = real(num_wann, dp)*kmesh_info%wbtot

      if (lsitesymmetry) then
        do nkp = 1, sitesym%nkptirr
          wkomegai1(sitesym%ir2ik(nkp)) = wkomegai1(sitesym%ir2ik(nkp))* &
                                          sitesym%nsymmetry/count(sitesym%kptsym(:, nkp) .eq. sitesym%ir2ik(nkp))
        enddo
      endif
      do nkp_loc = 1, ranknk
        nkp = global_k(nkp_loc)
        wkomegai1_loc(nkp_loc) = wkomegai1(nkp)
      end do

      do nkp_loc = 1, ranknk
        nkp = global_k(nkp_loc)
        if (ndimfroz(nkp) .gt. 0) then
          if (lsitesymmetry) then
            call set_error_fatal(error, 'not implemented in symmetry-adapted mode', comm)
            return
          endif
          do nn = 1, kmesh_info%nntot
            nkp2 = kmesh_info%nnlist(nkp, nn)
            call zgemm('C', 'N', ndimfroz(nkp), dis_manifold%ndimwin(nkp2), dis_manifold%ndimwin(nkp), &
                       cmplx_1, u_matrix_opt(:, :, nkp), num_bands, &
                       m_matrix_orig_local(:, :, nn, nkp_loc), num_bands, cmplx_0, cwb, num_wann)
            call zgemm('N', 'N', ndimfroz(nkp), num_wann, dis_manifold%ndimwin(nkp2), cmplx_1, cwb, &
                       num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
            rsum = 0.0_dp
            do n = 1, num_wann
              do m = 1, ndimfroz(nkp)
                rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
              enddo
            enddo
            wkomegai1_loc(nkp_loc) = wkomegai1_loc(nkp_loc) - kmesh_info%wb(nn)*rsum
          enddo
        endif
      enddo

      if (print_output%timing_level > 1) call io_stopwatch_stop('dis: extract_2', timer)
      if (print_output%timing_level > 1) call io_stopwatch_start('dis: extract_3', timer)

      !! ! send chunks of wkomegai1 to root node
      !! call comms_gatherv(wkomegai1_loc, counts(my_node_id), wkomegai1, counts, displs)
      !! ! send back the whole wkomegai1 array to other nodes
      !! call comms_bcast(wkomegai1(1), num_kpts)

      ! Refine optimal subspace at k points w/ non-frozen states
      do nkp_loc = 1, ranknk
        nkp = global_k(nkp_loc)
        if (lsitesymmetry) then
          if (sitesym%ir2ik(sitesym%ik2ir(nkp)) .ne. nkp) cycle

          call sitesym_dis_extract_symmetry(sitesym, lambda, u_matrix_opt_loc(:, :, nkp_loc), &
                                            czmat_in_loc(:, :, nkp_loc), nkp, &
                                            dis_manifold%ndimwin(nkp), num_bands, num_wann, &
                                            stdout, error, comm)
          if (allocated(error)) return
          do j = 1, num_wann
            wkomegai1_loc(nkp_loc) = wkomegai1_loc(nkp_loc) - real(lambda(j, j), kind=dp)
          enddo
        else
          if (num_wann .gt. ndimfroz(nkp)) then
            ! Diagonalize Z matrix
            do j = 1, dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
              do i = 1, j
                cap(i + ((j - 1)*j)/2) = czmat_in_loc(i, j, nkp_loc)
              enddo
            enddo
            ndiff = dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
            call ZHPEVX('V', 'A', 'U', ndiff, cap, 0.0_dp, 0.0_dp, 0, 0, &
                        -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)
            if (info .lt. 0) then
              if (on_root) then
                write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING Z MATRIX'
                write (stdout, *) ' THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
              endif
              call set_error_fatal(error, ' dis_extract: error', comm)
              return
            endif
            if (info .gt. 0) then
              if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING Z MATRIX'
              if (on_root) write (stdout, *) info, ' EIGENVECTORS FAILED TO CONVERGE'
              call set_error_fatal(error, ' dis_extract: error', comm)
              return
            endif

            ! Update the optimal subspace by incorporating the num_wann-ndimfroz(nkp) l
            ! eigenvectors of the Z matrix into u_matrix_opt. Also, add contribution from
            ! non-frozen states to wkomegai1(nkp) (minus the corresponding eigenvalu
            m = ndimfroz(nkp)
            do j = dis_manifold%ndimwin(nkp) - num_wann + 1, dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
              m = m + 1
              wkomegai1_loc(nkp_loc) = wkomegai1_loc(nkp_loc) - w(j)
              u_matrix_opt_loc(1:dis_manifold%ndimwin(nkp), m, nkp_loc) = cmplx_0
              ndimk = dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
              do i = 1, ndimk
                p = indxnfroz(i, nkp)
                u_matrix_opt_loc(p, m, nkp_loc) = cz(i, j)
              enddo
            enddo
          endif
          ! [if num_wann>ndimfroz(nkp)]
        endif !RS:

        ! Now that we have contribs. from both frozen and non-frozen states to
        ! wkomegai1(nkp), add it to womegai1
        womegai1 = womegai1 + wkomegai1_loc(nkp_loc)

        !if (index(print_output%devel_flag, 'compspace') > 0) then

        ! AT THE LAST ITERATION FIND A BASIS FOR THE (NDIMWIN(NKP)-num_wann)-DIMENS
        ! COMPLEMENT SPACE

        !  if (iter .eq. dis_control%num_iter) then
        !    allocate (camp(num_bands, num_bands, num_kpts), stat=ierr)
        !    if (ierr /= 0) then
        !      call set_error_alloc(error, 'Error allocating camp in dis_extract', comm)
        !      return
        !    endif
        !    allocate (camp_loc(num_bands, num_bands, max(1, counts(my_node_id))), stat=ierr)
        !    if (ierr /= 0) then
        !      call set_error_alloc(error, 'Error allocating ucamp_loc in dis_extract', comm)
        !      return
        !    endif
        !    if (dis_manifold%ndimwin(nkp) .gt. num_wann) then
        !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
        !        if (num_wann .gt. ndimfroz(nkp)) then
        !          ! USE THE NON-LEADING EIGENVECTORS OF THE Z-MATRIX
        !          camp_loc(1:dis_manifold%ndimwin(nkp), j, nkp_loc) = cz(1:dis_manifold%ndimwin(nkp), j)
        !        else
        ! Then num_wann=NDIMFROZ(NKP)
        ! USE THE ORIGINAL NON-FROZEN BLOCH EIGENSTATES
        !          do i = 1, dis_manifold%ndimwin(nkp)
        !            camp_loc(i, j, nkp_loc) = cmplx_0
        !            if (i .eq. indxnfroz(j, nkp)) camp_loc(i, j, nkp_loc) = cmplx_1
        !          enddo
        !        endif
        !      enddo
        !    else
        !      icompflag = 1
        !    endif
        !  endif

        !end if ! index(print_output%devel_flag,'compspace')>0

      enddo
      ! [Loop over k points (nkp)]

      !! ! send chunks of wkomegai1 to root node
      !! call comms_gatherv(wkomegai1_loc, counts(my_node_id), wkomegai1, counts, displs)
      !! ! send back the whole wkomegai1 array to other nodes
      !! call comms_bcast(wkomegai1(1), num_kpts)

      call comms_allreduce(womegai1, 1, 'SUM', error, comm)
      if (allocated(error)) return

      u_matrix_opt(:, :, :) = cmplx_0
      do nkp_loc = 1, ranknk
        nkp = global_k(nkp_loc)
        u_matrix_opt(:, :, nkp) = u_matrix_opt_loc(:, :, nkp_loc)
      enddo
      call comms_allreduce(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts, 'SUM', error, comm)
      if (allocated(error)) return

      if (lsitesymmetry) then
        call sitesym_symmetrize_u_matrix(sitesym, u_matrix_opt, num_bands, num_bands, num_kpts, &
                                         num_wann, stdout, error, comm, dis_manifold%lwindow)
        if (allocated(error)) return
        do nkp_loc = 1, ranknk
          nkp = global_k(nkp_loc)
          u_matrix_opt_loc(:, :, nkp_loc) = u_matrix_opt(:, :, nkp)
        enddo
      endif

      if (print_output%timing_level > 1) call io_stopwatch_stop('dis: extract_3', timer)

      womegai1 = womegai1/real(num_kpts, dp)

      ! DEBUG
      ! Orthonormality check
      !         do nkp=1,nkpts
      !           write(*,*) ' '
      !           write(*,'(a8,i4)') 'k-point ',nkp
      !           do l=1,num_wann
      !           do m=1,l
      !             ctmp=czero
      !             do j=1,ndimwin(nkp)
      !               ctmp=ctmp+conjg(u_matrix_opt(j,m,nkp))*u_matrix_opt(j,l,nkp)
      !             enddo
      !             write(*,'(i2,2x,i2,f16.12,1x,f16.12)') l,m,ctmp
      !             if(l.eq.m) then
      !               if(abs(ctmp-cmplx(1.0d0,0.0d0)).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             else
      !               if(abs(ctmp).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             endif
      !           enddo
      !           enddo
      !         enddo
      ! ENDDEBUG

      ! Compute womegai  using the updated subspaces at all k, i.e.,
      ! replacing (i-1) by (i) in Eq. (12) SMV
      if (print_output%timing_level > 1) call io_stopwatch_start('dis: extract_4', timer)

      womegai = 0.0_dp
      do nkp_loc = 1, ranknk
        nkp = global_k(nkp_loc)
        wkomegai = 0.0_dp
        do nn = 1, kmesh_info%nntot
          nkp2 = kmesh_info%nnlist(nkp, nn) ! nkp2 here may not be local to this processor hence operate on global u
          call zgemm('C', 'N', num_wann, dis_manifold%ndimwin(nkp2), dis_manifold%ndimwin(nkp), cmplx_1, &
                     u_matrix_opt(:, :, nkp), num_bands, m_matrix_orig_local(:, :, nn, nkp_loc), &
                     num_bands, cmplx_0, cwb, num_wann)
          call zgemm('N', 'N', num_wann, num_wann, dis_manifold%ndimwin(nkp2), cmplx_1, cwb, num_wann, &
                     u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
          rsum = 0.0_dp
          do n = 1, num_wann
            do m = 1, num_wann
              rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
            enddo
          enddo
          wkomegai = wkomegai + kmesh_info%wb(nn)*rsum
        enddo
        wkomegai = real(num_wann, dp)*kmesh_info%wbtot - wkomegai
        womegai = womegai + wkomegai
      enddo

      call comms_allreduce(womegai, 1, 'SUM', error, comm)
      if (allocated(error)) return
      womegai = womegai/real(num_kpts, dp)
      delta_womegai = womegai1/womegai - 1.0_dp

      if (print_output%timing_level > 1) call io_stopwatch_stop('dis: extract_4', timer)

      if (print_output%iprint > 0) then
        write (stdout, 124) iter, womegai1*print_output%lenconfac**2, &
          womegai*print_output%lenconfac**2, delta_womegai, io_wallclocktime()
      endif
124   format(2x, i6, 3x, f14.8, 3x, f14.8, 6x, es10.3, 2x, f8.2, 4x, '<-- DIS')

      ! Construct the updated Z matrix, CZMAT_OUT, at k points w/ non-frozen s
      do nkp_loc = 1, ranknk
        nkp = global_k(nkp_loc)
        if (num_wann .gt. ndimfroz(nkp)) then
          call internal_zmatrix(cbw, czmat_out_loc(:, :, nkp_loc), m_matrix_orig_local, &
                                u_matrix_opt, kmesh_info%wb, indxnfroz, ndimfroz, &
                                dis_manifold%ndimwin, kmesh_info%nnlist, nkp, nkp_loc, &
                                kmesh_info%nntot, num_bands, num_wann, print_output%timing_level, &
                                timer)
        endif
      enddo

      if (lsitesymmetry) then
        czmat_out = cmplx_0
        do nkp_loc = 1, ranknk
          nkp = global_k(nkp_loc)
          czmat_out(:, :, nkp) = czmat_out_loc(:, :, nkp_loc)
        enddo
        call comms_allreduce(czmat_out(1, 1, 1), num_bands*num_bands*num_kpts, 'SUM', error, comm)
        if (allocated(error)) return
        call sitesym_symmetrize_zmatrix(sitesym, czmat_out, num_bands, num_kpts, dis_manifold%lwindow)
        do nkp_loc = 1, ranknk
          nkp = global_k(nkp_loc)
          czmat_out_loc(:, :, nkp_loc) = czmat_out(:, :, nkp)
        end do
      end if

      call internal_test_convergence(history, delta_womegai, dis_control%conv_tol, iter, &
                                     dis_control%conv_window, dis_converged, error, comm)
      if (allocated(error)) return

      if (dis_converged) then
        if (print_output%iprint > 0) then
          write (stdout, '(/13x,a,es10.3,a,i2,a)') '<<<      Delta <', dis_control%conv_tol, &
            '  over ', dis_control%conv_window, ' iterations     >>>'
          write (stdout, '(13x,a)') '<<< Disentanglement convergence criteria satisfied >>>'
        endif
        exit
      endif

    enddo
    ! [BIG ITERATION LOOP (iter)]

    deallocate (czmat_out, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating czmat_out in dis_extract', comm)
      return
    endif
    deallocate (czmat_in, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating czmat_in in dis_extract', comm)
      return
    endif
    deallocate (czmat_out_loc, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating czmat_out_loc in dis_extract', comm)
      return
    endif
    deallocate (czmat_in_loc, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating czmat_in_loc in dis_extract', comm)
      return
    endif

    if (on_root) then
      allocate (ceamp(num_bands, num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating ceamp in dis_extract', comm)
        return
      endif
      allocate (cham(num_bands, num_bands, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating cham in dis_extract', comm)
        return
      endif
    endif

    if (.not. dis_converged) then
      if (on_root) then
        write (stdout, '(/5x,a)') &
          '<<< Warning: Maximum number of disentanglement iterations reached >>>'
        write (stdout, '(10x,a)') '<<< Disentanglement convergence criteria not satisfied >>>'
      endif
    endif

    !if (index(print_output%devel_flag, 'compspace') > 0) then

    !  if (icompflag .eq. 1) then
    !    if (print_output%iprint > 2) then
    !      if (on_root) write (stdout, ('(/4x,a)')) &
    !        'WARNING: Complement subspace has zero dimensions at the following k-points:'
    !      i = 0
    !      if (on_root) write (stdout, '(4x)', advance='no')
    !      do nkp = 1, num_kpts
    !        if (dis_manifold%ndimwin(nkp) .eq. num_wann) then
    !          i = i + 1
    !          if (i .le. 12) then
    !            if (on_root) write (stdout, '(i6)', advance='no') nkp
    !          else
    !            i = 1
    !            if (on_root) write (stdout, '(/4x)', advance='no')
    !            if (on_root) write (stdout, '(i6)', advance='no') nkp
    !          endif
    !        endif
    !      enddo
    !    endif
    !  endif

    !endif

    ! Write the final womegai. This should remain unchanged during the
    ! subsequent minimization of Omega_tilde in wannierise.f90
    ! We store it in the checkpoint file as a sanity check
    if (print_output%iprint > 0) write (stdout, '(/8x,a,f14.8,a/)') 'Final Omega_I ', &
      womegai*print_output%lenconfac**2, ' ('//trim(print_output%length_unit)//'^2)'

    ! Set public variable omega_invariant
    omega_invariant = womegai

    ! Currently, this part is not parallelized; thus, we perform the task only on root and then broadcast the result.
    if (on_root) then
      ! DIAGONALIZE THE HAMILTONIAN WITHIN THE OPTIMIZED SUBSPACES
      do nkp = 1, num_kpts

        do j = 1, num_wann
          do i = 1, num_wann
            cham(i, j, nkp) = cmplx_0
            do l = 1, dis_manifold%ndimwin(nkp)
              cham(i, j, nkp) = cham(i, j, nkp) + conjg(u_matrix_opt(l, i, nkp)) &
                                *u_matrix_opt(l, j, nkp)*eigval_opt(l, nkp)
            enddo
          enddo
        enddo

        do j = 1, num_wann
          do i = 1, j
            cap(i + ((j - 1)*j)/2) = cham(i, j, nkp)
          enddo
        enddo

        call ZHPEVX('V', 'A', 'U', num_wann, cap, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, m, w, cz, &
                    num_bands, cwork, rwork, iwork, ifail, info)

        if (info .lt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
          if (on_root) write (stdout, *) ' THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
          call set_error_fatal(error, ' dis_extract: error', comm)
          return
        endif
        if (info .gt. 0) then
          if (on_root) write (stdout, *) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
          if (on_root) write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
          call set_error_fatal(error, ' dis_extract: error', comm)
          return
        endif

        ! Store the energy eigenvalues of the optimal subspace (used in wann_ban
        eigval_opt(1:num_wann, nkp) = w(1:num_wann)

        ! CALCULATE AMPLITUDES OF THE CORRESPONDING ENERGY EIGENVECTORS IN TERMS
        ! THE ORIGINAL ("WINDOW SPACE") ENERGY EIGENVECTORS
        do j = 1, num_wann
          do i = 1, dis_manifold%ndimwin(nkp)
            ceamp(i, j, nkp) = cmplx_0
            do l = 1, num_wann
              ceamp(i, j, nkp) = ceamp(i, j, nkp) + cz(l, j)*u_matrix_opt(i, l, nkp)
            enddo
          enddo
        enddo
        ! NKP
      enddo

      if (print_output%iprint > 2) then
        if (on_root) write (stdout, '(/,a,/)') '  Eigenvalues inside optimal subspace:'
        do nkp = 1, num_kpts
          if (on_root) write (stdout, '(a,i3,2x,20(f9.5,1x))') '  K-point ', &
            nkp, (eigval_opt(i, nkp), i=1, num_wann)
        enddo
      endif

      ! Replace u_matrix_opt by ceamp. Both span the
      ! same space, but the latter is more convenient for the purpose of obtai
      ! an optimal Fourier-interpolated band structure: see Sec. III.E of SMV.
      if (.not. lsitesymmetry) then                                                                         !YN:
        do nkp = 1, num_kpts
          do j = 1, num_wann
            u_matrix_opt(1:dis_manifold%ndimwin(nkp), j, nkp) = ceamp(1:dis_manifold%ndimwin(nkp), j, nkp)
          enddo
        enddo
        !else                                                                                                        !YN:
        ! Above is skipped as we require Uopt(Rk) to be related to Uopt(k)                                        !YN: RS:
        !write(stdout,"(a)")  &                                                                                   !YN: RS:
        !  'Note(symmetry-adapted mode): u_matrix_opt are no longer the eigenstates of the subspace Hamiltonian.' !RS:
      endif                                                                                                        !YN:
    endif
    call comms_bcast(eigval_opt(1, 1), num_bands*num_kpts, error, comm)
    if (allocated(error)) return
    call comms_bcast(u_matrix_opt(1, 1, 1), num_bands*num_wann*num_kpts, error, comm)
    if (allocated(error)) return

    !if (index(print_output%devel_flag, 'compspace') > 0) then

    !  if (icompflag .eq. 1) then
    !    if (print_output%iprint > 2) then
    !      if (on_root) then
    !        write (stdout, *) 'AT SOME K-POINT(S) COMPLEMENT SUBSPACE HAS ZERO DIMENSIONALITY'
    !        write (stdout, *) '=> DID NOT CREATE FILE COMPSPACE.DAT'
    !      endif
    !    endif
    !  else
    ! DIAGONALIZE THE HAMILTONIAN IN THE COMPLEMENT SUBSPACE, WRITE THE
    ! CORRESPONDING EIGENFUNCTIONS AND ENERGY EIGENVALUES
    !    do nkp = 1, num_kpts
    !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
    !        do i = 1, dis_manifold%ndimwin(nkp) - num_wann
    !          cham(i, j, nkp) = cmplx_0
    !          do l = 1, dis_manifold%ndimwin(nkp)
    !            cham(i, j, nkp) = cham(i, j, nkp) + conjg(camp(l, i, nkp)) &
    !                              *camp(l, j, nkp)*eigval_opt(l, nkp)
    !          enddo
    !        enddo
    !      enddo
    !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
    !        do i = 1, j
    !          cap(i + ((j - 1)*j)/2) = cham(i, j, nkp)
    !        enddo
    !      enddo
    !      ndiff = dis_manifold%ndimwin(nkp) - num_wann
    !      call ZHPEVX('V', 'A', 'U', ndiff, cap, 0.0_dp, 0.0_dp, 0, 0, &
    !                  -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)
    !      if (info .lt. 0) then
    !        if (on_root) write (stdout, *) '*** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
    !        if (on_root) write (stdout, *) 'THE ', -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
    !        call io_error(' dis_extract: error', stdout, seedname)
    !      endif
    !      if (info .gt. 0) then
    !        if (on_root) write (stdout, *) '*** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
    !        if (on_root) write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
    !        call io_error(' dis_extract: error', stdout, seedname)
    !      endif
    ! CALCULATE AMPLITUDES OF THE ENERGY EIGENVECTORS IN THE COMPLEMENT SUBS
    ! TERMS OF THE ORIGINAL ENERGY EIGENVECTORS
    !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
    !        do i = 1, dis_manifold%ndimwin(nkp)
    !          camp(i, j, nkp) = cmplx_0
    !          do l = 1, dis_manifold%ndimwin(nkp) - num_wann
!write(stdout,*) 'i=',i,'   j=',j,'   l=',l
!write(stdout,*) '           camp(i,j,nkp)=',camp(i,j,nkp)
!write(stdout,*) '           cz(l,j)=',cz(l,j)
!write(stdout,*) '           u_matrix_opt(i,l,nkp)=',u_matrix_opt(i,l,nkp)

! aam: 20/10/2006 -- the second dimension of u_matrix_opt is out of bounds (allocated as num_wann)!
! commenting this line out.
!                     camp(i,j,nkp) = camp(i,j,nkp) + cz(l,j) * u_matrix_opt(i,l,nkp)
    !         enddo
    !       enddo
    !     enddo
    !   enddo   ! [loop over k points (nkp)]

    ! endif   ! [if icompflag=1]

    !endif     ![if(index(devel_flag,'compspace')>0)]

    deallocate (history, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating history in dis_extract', comm)
      return
    endif

    if (on_root) then
      deallocate (cham, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating cham in dis_extract', comm)
        return
      endif
    endif
    if (allocated(camp)) then
      deallocate (camp, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating camp in dis_extract', comm)
        return
      endif
    end if
    if (allocated(camp_loc)) then
      deallocate (camp_loc, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating camp_loc in dis_extract', comm)
        return
      endif
    endif
    if (on_root) then
      deallocate (ceamp, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating ceamp in dis_extract', comm)
        return
      endif
    endif
    deallocate (u_matrix_opt_loc, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating u_matrix_opt_loc in dis_extract', comm)
      return
    endif
    deallocate (wkomegai1_loc, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating wkomegai1_loc in dis_extract', comm)
      return
    endif
    deallocate (wkomegai1, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating wkomegai1 in dis_extract', comm)
      return
    endif

    deallocate (cz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cz in dis_extract', comm)
      return
    endif
    deallocate (cwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cwork in dis_extract', comm)
      return
    endif
    deallocate (cap, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cap in dis_extract', comm)
      return
    endif
    deallocate (rwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating rwork in dis_extract', comm)
      return
    endif
    deallocate (w, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating w in dis_extract', comm)
      return
    endif
    deallocate (ifail, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating ifail in dis_extract', comm)
      return
    endif
    deallocate (iwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating iwork in dis_extract', comm)
      return
    endif

    deallocate (cbw, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cbw in dis_extract', comm)
      return
    endif
    deallocate (cww, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cww in dis_extract', comm)
      return
    endif
    deallocate (cwb, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cwb in dis_extract', comm)
      return
    endif

    if (print_output%iprint > 0) write (stdout, '(1x,a/)') &
      '+----------------------------------------------------------------------------+'

    if (print_output%timing_level > 1) call io_stopwatch_stop('dis: extract', timer)

    return
    !================================================!
  end subroutine dis_extract

  subroutine internal_test_convergence(history, delta_womegai, dis_conv_tol, iter, &
                                       dis_conv_window, dis_converged, error, comm)
    !================================================!
    !
    !! Check if we have converged
    !
    !================================================!
    use w90_constants, only: dp
    use w90_error
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: iter, dis_conv_window

    real(kind=dp), intent(inout) :: history(:)
    real(kind=dp), intent(in) :: delta_womegai, dis_conv_tol

    logical, intent(inout) :: dis_converged

    ! local variables
    integer :: ierr
    real(kind=dp), allocatable :: temp_hist(:)

    allocate (temp_hist(dis_conv_window), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating temp_hist in dis_extract', comm)
      return
    endif

    if (iter .le. dis_conv_window) then
      history(iter) = delta_womegai
    else
      temp_hist = eoshift(history, 1, delta_womegai)
      history = temp_hist
    endif

    dis_converged = .false.
    if (iter .ge. dis_conv_window) then
      dis_converged = all(abs(history) .lt. dis_conv_tol)
    endif

    deallocate (temp_hist, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating temp_hist in dis_extract', comm)
      return
    endif

    return
    !================================================!
  end subroutine internal_test_convergence

  subroutine internal_zmatrix(cbw, cmtrx, m_matrix_orig_local, u_matrix_opt, wb, indxnfroz, &
                              ndimfroz, ndimwin, nnlist, nkp, nkp_loc, nntot, num_bands, num_wann, &
                              timing_level, timer)
    !================================================!
    !
    !! Compute the Z-matrix
    !
    !================================================!
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer

    integer, intent(in) :: num_bands, num_wann
    integer, intent(in) :: timing_level
    integer, intent(in) :: ndimwin(:)
    integer, intent(in) :: nntot, nnlist(:, :)
    integer, intent(in) :: indxnfroz(:, :)
    integer, intent(in) :: ndimfroz(:)
    integer, intent(in) :: nkp
    integer, intent(in) :: nkp_loc

    real(kind=dp), intent(in) :: wb(:)

    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)
    complex(kind=dp), intent(in) :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), intent(inout) :: cbw(:, :)
    complex(kind=dp), intent(out) :: cmtrx(:, :)
    !! (M,N)-TH ENTRY IN THE (NDIMWIN(NKP)-NDIMFROZ(NKP)) x (NDIMWIN(NKP)-NDIMFRO
    !! HERMITIAN MATRIX AT THE NKP-TH K-POINT

    ! local variables
    integer          :: l, m, n, p, q, nn, nkp2, ndimk
    complex(kind=dp) :: csum

    if (timing_level > 1) call io_stopwatch_start('dis: extract: zmatrix', timer)

    cmtrx = cmplx_0
    ndimk = ndimwin(nkp) - ndimfroz(nkp)
    do nn = 1, nntot
      nkp2 = nnlist(nkp, nn)
      call zgemm('N', 'N', num_bands, num_wann, ndimwin(nkp2), cmplx_1, &
                 m_matrix_orig_local(:, :, nn, nkp_loc), num_bands, u_matrix_opt(:, :, nkp2), &
                 num_bands, cmplx_0, cbw, num_bands)
      do n = 1, ndimk
        q = indxnfroz(n, nkp)
        do m = 1, n
          p = indxnfroz(m, nkp)
          csum = cmplx_0
          do l = 1, num_wann
            csum = csum + cbw(p, l)*conjg(cbw(q, l))
          enddo
          cmtrx(m, n) = cmtrx(m, n) + cmplx(wb(nn), 0.0_dp, kind=dp)*csum
          cmtrx(n, m) = conjg(cmtrx(m, n))
        enddo
      enddo
    enddo

    if (timing_level > 1) call io_stopwatch_stop('dis: extract: zmatrix', timer)

    return
    !================================================!
  end subroutine internal_zmatrix

  subroutine dis_extract_gamma(dis_control, kmesh_info, print_output, dis_manifold, &
                               m_matrix_orig_local, u_matrix_opt, eigval_opt, omega_invariant, &
                               indxnfroz, ndimfroz, num_bands, num_kpts, num_wann, timer, error, &
                               stdout, comm)
    !================================================!
    !
    !! Extracts an num_wann-dimensional subspace at each k by
    !! minimizing Omega_I (Gamma point version)
    !
    !================================================!
    ! this routine is not parallelised
    use w90_comms, only: w90_comm_type
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_error
    use w90_io, only: io_time, io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: dis_manifold_type, kmesh_info_type, print_output_type, timer_list_type
    use w90_types, only: timer_list_type
    use w90_wannier90_types, only: dis_control_type

    implicit none

    ! arguments
    type(dis_control_type), intent(in) :: dis_control
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: stdout
    integer, intent(in) :: num_bands, num_kpts, num_wann
    integer, intent(in) :: ndimfroz(:)
    integer, intent(in) :: indxnfroz(:, :)

    real(kind=dp), intent(inout) :: eigval_opt(:, :)
    real(kind=dp), intent(out) :: omega_invariant

    complex(kind=dp), intent(inout) :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    ! MODIFIED:
    !           u_matrix_opt (At input it contains the initial guess for the optimal
    ! subspace (expressed in terms of the original states inside the window). At
    ! output it contains the  states that diagonalize the hamiltonian inside the
    ! optimal subspace (again expressed in terms of the original window states).
    ! Giving out states that diagonalize the hamiltonian inside the optimal
    ! subspace (instead of the eigenstates of the Z matrix) is useful for
    ! performing the Wannier interpolation of the energy bands as described in
    ! Sec. III.F of SMV)
    !
    !           eigval (At input: original energy eigenvalues.
    ! At output: eigenvalues of the hamiltonian inside optimal subspace)

    ! ----------------------------------------------------------------------
    ! TO DO: The complement subspace is computed but is not saved anywhere!
    ! (Check what was done with it in original code space.f)
    ! Diagonalize Z matrix only at those k points where ndimwin>num_wann?
    ! ----------------------------------------------------------------------

    ! *******************
    ! SHELLS OF K-VECTORS
    ! *******************
    ! nshells           number of shells of k-points to be used in the
    !                   finite-difference formulas for the k-derivatives
    ! aam: wb is now wb(1:nntot) 09/04/2006
    ! wb(nkp,nnx)       weight of the nnx-th b-vector (ordered along shells
    !                   of increasing length) associated with the nkp-th k-p
    ! wbtot             sum of the weights of all b-vectors associated with
    !                   given k-point (k-point 1 is used in calculation)
    ! nnlist(nkp,nnx)   vkpt(1:3,nnlist(nkp,nnx)) is the nnx-th neighboring
    !                   k-point of the nkp-th k-point vkpt(1:3,nkp) (or its
    !                   periodic image in the "home Brillouin zone")
    ! cm(n,m,nkp,nnx)   Overlap matrix <u_nk|u_{m,k+b}>

    ! Internal variables
    integer, allocatable :: ifail(:)
    integer, allocatable :: iwork(:)
    integer :: icompflag, iter, ndiff
    integer :: i, j, l, m, n, nn, nkp, nkp2, info, ierr, ndimk, p

    complex(kind=dp), allocatable :: camp(:, :, :)
    complex(kind=dp), allocatable :: ceamp(:, :, :)
    complex(kind=dp), allocatable :: cham(:, :, :)
    complex(kind=dp), allocatable :: cwb(:, :), cww(:, :), cbw(:, :)
    complex(kind=dp), allocatable :: cz(:, :)

    real(kind=dp), allocatable :: cap_r(:)
    real(kind=dp), allocatable :: history(:)
    real(kind=dp), allocatable :: rz(:, :)
    real(kind=dp), allocatable :: rzmat_in(:, :, :)
    real(kind=dp), allocatable :: rzmat_out(:, :, :)
    real(kind=dp), allocatable :: w(:)
    real(kind=dp), allocatable :: wkomegai1(:)
    real(kind=dp), allocatable :: work(:)
    real(kind=dp) :: womegai, wkomegai, womegai1, rsum, delta_womegai

    logical :: dis_converged

    if (print_output%timing_level > 1) call io_stopwatch_start('dis: extract', timer)

    write (stdout, '(/1x,a)') &
      '                  Extraction of optimally-connected subspace                  '
    write (stdout, '(1x,a)') &
      '                  ------------------------------------------                  '

    allocate (cwb(num_wann, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cwb in dis_extract_gamma', comm)
      return
    endif
    allocate (cww(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cww in dis_extract_gamma', comm)
      return
    endif
    allocate (cbw(num_bands, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cbw in dis_extract_gamma', comm)
      return
    endif

    cwb = cmplx_0; cww = cmplx_0; cbw = cmplx_0

    allocate (iwork(5*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating iwork in dis_extract_gamma', comm)
      return
    endif
    allocate (ifail(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating ifail in dis_extract_gamma', comm)
      return
    endif
    allocate (w(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating w in dis_extract_gamma', comm)
      return
    endif
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cz in dis_extract_gamma', comm)
      return
    endif
    allocate (work(8*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating work in dis_extract_gamma', comm)
      return
    endif
    allocate (cap_r((num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cap_r in dis_extract_gamma', comm)
      return
    endif
    allocate (rz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating rz in dis_extract_gamma', comm)
      return
    endif
    allocate (wkomegai1(num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating wkomegai1 in dis_extract_gamma', comm)
      return
    endif
    allocate (rzmat_in(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating rzmat_in indis_extract_gamma', comm)
      return
    endif
    allocate (rzmat_out(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating rzmat_out indis_extract_gamma', comm)
      return
    endif
    allocate (history(dis_control%conv_window), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating history indis_extract_gamma', comm)
      return
    endif

    ! ********************************************
    ! ENERGY WINDOWS AND SUBSPACES AT EACH K-POINT
    ! ********************************************
    ! num_wann             dimensionality of the subspace at each k-point
    !                   (number of Wannier functions per unit cell that we w
    ! NDIMWIN(NKP)      number of bands at the nkp-th k-point that fall
    !                   within the outer energy window
    ! NDIMFROZ(NKP)     number of frozen bands at the nkp-th k-point
    ! INDXNFROZ(I,NKP)  INDEX (BETWEEN 1 AND NDIMWIN(NKP)) OF THE I-TH NON-F
    !                   ORIGINAL BAND STATE AT THE NKP-TH K-POINT
    ! U_MATRIX_OPT(J,L,NKP) AMPLITUDE OF THE J-TH ENERGY EIGENVECTOR INSIDE THE
    !                   ENERGY WINDOW AT THE NKP-TH K-POINT IN THE EXPANSION OF
    !                   THE L-TH LEADING RLAMBDA EIGENVECTOR AT THE SAME K-POINT
    !                   If there are M_k frozen states, they occupy the lowest
    !                   entries of the second index of u_matrix_opt, and the leading
    !                   nbands-M_k eigenvectors of the Z matrix occupy the
    !                   remaining slots
    ! CAMP(J,L,NKP)     SAME AS U_MATRIX_OPT, BUT FOR THE COMPLEMENT SUBSPACE INSIDE THE
    !                   ENERGY WINDOW (I.E., THE NON-LEADING RLAMBDA EIGENVECTORS)
    ! CEAMP(J,L,NKPTS)  SAME AS U_MATRIX_OPT, BUT INSTEAD OF RLAMBDA EIGENVECTOR, I
    !                   FOR THE ENERGY EIGENVECTOR OBTAINED BY DIAGONALIZING
    !                   HAMILTONIAN IN THE OPTIMIZED SUBSPACE
    ! RZMAT_IN(M,N,NKP) Z-MATRIX [Eq. (21) SMV]
    ! RZMAT_OUT(M,N,NKP) OUTPUT Z-MATRIX FROM THE PRESENT ITERATION
    ! RLAMBDA           An eigenvalue of the Z matrix
    ! womegai           Gauge-invariant Wannier spread, computed usinf all s
    !                   from current iteration
    ! wkomegai1(NKP)    Eq. (18) of SMV
    ! womegai1          Eq.(11) of SMV (like wowmegai, but neighboring state
    !                   for computing overlaps are from previous iteration.
    !                   become equal at self-consistency)
    ! alphafixe         mixing parameter for the iterative procedure
    ! nitere            total number of iterations

    ! DEBUG
    if (print_output%iprint > 2) then
      write (stdout, '(a,/)') '  Original eigenvalues inside outer window:'
      do nkp = 1, num_kpts
        write (stdout, '(a,i3,3x,20(f9.5,1x))') '  K-point ', nkp, &
          (eigval_opt(i, nkp), i=1, dis_manifold%ndimwin(nkp))
      enddo
    endif
    ! ENDDEBUG

    ! TO DO: Check if this is the best place to initialize icompflag
    icompflag = 0

    write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'
    write (stdout, '(1x,a)') &
      '|  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS'
    write (stdout, '(1x,a)') &
      '+---------------------------------------------------------------------+<-- DIS'

    dis_converged = .false.

    ! ------------------
    ! BIG ITERATION LOOP
    ! ------------------
    do iter = 1, dis_control%num_iter

      if (iter .eq. 1) then
        ! Initialize Z matrix at k points w/ non-frozen states
        do nkp = 1, num_kpts
          if (num_wann .gt. ndimfroz(nkp)) then
            call internal_zmatrix_gamma(cbw, m_matrix_orig_local, u_matrix_opt, rzmat_in(:, :, nkp), &
                                        kmesh_info%wb, indxnfroz, ndimfroz, dis_manifold%ndimwin, &
                                        kmesh_info%nnlist, nkp, kmesh_info%nntot, num_bands, &
                                        num_wann, print_output%timing_level, timer)
          endif
        enddo
      else
        ! [iter.ne.1]
        ! Update Z matrix at k points with non-frozen states, using a mixing sch
        do nkp = 1, num_kpts
          if (num_wann .gt. ndimfroz(nkp)) then
            ndimk = dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, ndimk
              do j = 1, i
                rzmat_in(j, i, nkp) = &
                  dis_control%mix_ratio*rzmat_out(j, i, nkp) &
                  + (1.0_dp - dis_control%mix_ratio)*rzmat_in(j, i, nkp)
                ! hermiticity
                rzmat_in(i, j, nkp) = rzmat_in(j, i, nkp)
              enddo
            enddo
          endif
        enddo
      endif
      ! [if iter=1]

      womegai1 = 0.0_dp
      ! wkomegai1 is defined by Eq. (18) of SMV.
      ! Contribution to wkomegai1 from frozen states should be calculated now
      ! every k (before updating any k), so that for iter>1 overlaps are with
      ! non-frozen neighboring states from the previous iteration

      wkomegai1 = real(num_wann, dp)*kmesh_info%wbtot
      do nkp = 1, num_kpts
        if (ndimfroz(nkp) .gt. 0) then
          do nn = 1, kmesh_info%nntot
            nkp2 = kmesh_info%nnlist(nkp, nn)
            call zgemm('C', 'N', ndimfroz(nkp), dis_manifold%ndimwin(nkp2), dis_manifold%ndimwin(nkp), &
                       cmplx_1, u_matrix_opt(:, :, nkp), num_bands, m_matrix_orig_local(:, :, nn, nkp), &
                       num_bands, cmplx_0, cwb, num_wann)
            call zgemm('N', 'N', ndimfroz(nkp), num_wann, dis_manifold%ndimwin(nkp2), cmplx_1, &
                       cwb, num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
            rsum = 0.0_dp
            do n = 1, num_wann
              do m = 1, ndimfroz(nkp)
                rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
              enddo
            enddo
            wkomegai1(nkp) = wkomegai1(nkp) - kmesh_info%wb(nn)*rsum
          enddo
        endif
      enddo

      ! Refine optimal subspace at k points w/ non-frozen states
      do nkp = 1, num_kpts
        if (num_wann .gt. ndimfroz(nkp)) then
          ! Diagonalize Z matrix
          do j = 1, dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, j
              cap_r(i + ((j - 1)*j)/2) = rzmat_in(i, j, nkp)
            enddo
          enddo
          ndiff = dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
          call DSPEVX('V', 'A', 'U', ndiff, cap_r, 0.0_dp, 0.0_dp, 0, 0, &
                      -1.0_dp, m, w, rz, num_bands, work, iwork, ifail, info)
          if (info .lt. 0) then
            write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING Z MATRIX'
            write (stdout, *) ' THE ', -info, ' ARGUMENT OF DSPEVX HAD AN ILLEGAL VALUE'
            call set_error_fatal(error, ' dis_extract_gamma: error', comm)
            return
          endif
          if (info .gt. 0) then
            write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING Z MATRIX'
            write (stdout, *) info, ' EIGENVECTORS FAILED TO CONVERGE'
            call set_error_fatal(error, '  dis_extract_gamma: error', comm)
            return
          endif
          cz(:, :) = cmplx(rz(:, :), 0.0_dp, dp)
          !
          ! Update the optimal subspace by incorporating the num_wann-ndimfroz(nkp) l
          ! eigenvectors of the Z matrix into u_matrix_opt. Also, add contribution from
          ! non-frozen states to wkomegai1(nkp) (minus the corresponding eigenvalu
          m = ndimfroz(nkp)
          do j = dis_manifold%ndimwin(nkp) - num_wann + 1, dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
            m = m + 1
            wkomegai1(nkp) = wkomegai1(nkp) - w(j)
            u_matrix_opt(1:dis_manifold%ndimwin(nkp), m, nkp) = cmplx_0
            ndimk = dis_manifold%ndimwin(nkp) - ndimfroz(nkp)
            do i = 1, ndimk
              p = indxnfroz(i, nkp)
              u_matrix_opt(p, m, nkp) = cz(i, j)
            enddo
          enddo
        endif
        ! [if num_wann>ndimfroz(nkp)]

        ! Now that we have contribs. from both frozen and non-frozen states to
        ! wkomegai1(nkp), add it to womegai1
        womegai1 = womegai1 + wkomegai1(nkp)

        ! AT THE LAST ITERATION FIND A BASIS FOR THE (NDIMWIN(NKP)-num_wann)-DIMENS
        ! COMPLEMENT SPACE
        !if (index(print_output%devel_flag, 'compspace') > 0) then

        !  if (iter .eq. dis_control%num_iter) then
        !    allocate (camp(num_bands, num_bands, num_kpts), stat=ierr)
        !    camp = cmplx_0
        !    if (ierr /= 0) then
        !      call set_error_alloc(error, 'Error allocating camp in dis_extract_gamma', comm)
        !      return
        !    endif
        !    if (dis_manifold%ndimwin(nkp) .gt. num_wann) then
        !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
        !        if (num_wann .gt. ndimfroz(nkp)) then
        !          ! USE THE NON-LEADING EIGENVECTORS OF THE Z-MATRIX
        !          camp(1:dis_manifold%ndimwin(nkp), j, nkp) = cz(1:dis_manifold%ndimwin(nkp), j)
        !        else
        !          ! Then num_wann=NDIMFROZ(NKP)
        !          ! USE THE ORIGINAL NON-FROZEN BLOCH EIGENSTATES
        !          do i = 1, dis_manifold%ndimwin(nkp)
        !            camp(i, j, nkp) = cmplx_0
        !            if (i .eq. indxnfroz(j, nkp)) camp(i, j, nkp) = cmplx_1
        !          enddo
        !        endif
        !      enddo
        !    else
        !      icompflag = 1
        !    endif
        !  endif
        !end if

      enddo
      ! [Loop over k points (nkp)]

      womegai1 = womegai1/real(num_kpts, dp)

      ! DEBUG
      ! Orthonormality check
      !         do nkp=1,nkpts
      !           write(*,*) ' '
      !           write(*,'(a8,i4)') 'k-point ',nkp
      !           do l=1,num_wann
      !           do m=1,l
      !             ctmp=czero
      !             do j=1,ndimwin(nkp)
      !               ctmp=ctmp+conjg(u_matrix_opt(j,m,nkp))*u_matrix_opt(j,l,nkp)
      !             enddo
      !             write(*,'(i2,2x,i2,f16.12,1x,f16.12)') l,m,ctmp
      !             if(l.eq.m) then
      !               if(abs(ctmp-cmplx(1.0d0,0.0d0)).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             else
      !               if(abs(ctmp).gt.1.0e-8) then
      !                 write(*,'(a49,i4)')
      !     1           '*** ERROR *** with iterative subspace at k-point ',
      !     2           nkp
      !                 write(*,*) 'vectors in u_matrix_opt not orthonormal'
      !                 stop
      !               endif
      !             endif
      !           enddo
      !           enddo
      !         enddo
      ! ENDDEBUG

      ! Compute womegai  using the updated subspaces at all k, i.e.,
      ! replacing (i-1) by (i) in Eq. (12) SMV

      womegai = 0.0_dp
      do nkp = 1, num_kpts
        wkomegai = 0.0_dp
        do nn = 1, kmesh_info%nntot
          nkp2 = kmesh_info%nnlist(nkp, nn)
          call zgemm('C', 'N', num_wann, dis_manifold%ndimwin(nkp2), dis_manifold%ndimwin(nkp), &
                     cmplx_1, u_matrix_opt(:, :, nkp), num_bands, m_matrix_orig_local(:, :, nn, nkp), &
                     num_bands, cmplx_0, cwb, num_wann)
          call zgemm('N', 'N', num_wann, num_wann, dis_manifold%ndimwin(nkp2), cmplx_1, cwb, &
                     num_wann, u_matrix_opt(:, :, nkp2), num_bands, cmplx_0, cww, num_wann)
          rsum = 0.0_dp
          do n = 1, num_wann
            do m = 1, num_wann
              rsum = rsum + real(cww(m, n), dp)**2 + aimag(cww(m, n))**2
            enddo
          enddo
          wkomegai = wkomegai + kmesh_info%wb(nn)*rsum
        enddo
        wkomegai = real(num_wann, dp)*kmesh_info%wbtot - wkomegai
        womegai = womegai + wkomegai
      enddo
      womegai = womegai/real(num_kpts, dp)
      ! [Loop over k (nkp)]

      delta_womegai = womegai1/womegai - 1.0_dp

      write (stdout, 124) iter, womegai1*print_output%lenconfac**2, &
        womegai*print_output%lenconfac**2, delta_womegai, io_time()

124   format(2x, i6, 3x, f14.8, 3x, f14.8, 6x, es10.3, 2x, f8.2, 4x, '<-- DIS')

      ! Construct the updated Z matrix, CZMAT_OUT, at k points w/ non-frozen s
      do nkp = 1, num_kpts
        if (num_wann .gt. ndimfroz(nkp)) then
          call internal_zmatrix_gamma(cbw, m_matrix_orig_local, u_matrix_opt, rzmat_out(:, :, nkp), &
                                      kmesh_info%wb, indxnfroz, ndimfroz, dis_manifold%ndimwin, &
                                      kmesh_info%nnlist, nkp, kmesh_info%nntot, num_bands, &
                                      num_wann, print_output%timing_level, timer)
        endif
      enddo

      call internal_test_convergence(history, delta_womegai, dis_control%conv_tol, iter, &
                                     dis_control%conv_window, dis_converged, error, comm)
      if (allocated(error)) return
      if (dis_converged) then
        write (stdout, '(/13x,a,es10.3,a,i2,a)') &
          '<<<      Delta <', dis_control%conv_tol, &
          '  over ', dis_control%conv_window, ' iterations     >>>'
        write (stdout, '(13x,a)') '<<< Disentanglement convergence criteria satisfied >>>'
        exit
      endif

    enddo
    ! [BIG ITERATION LOOP (iter)]

    deallocate (rzmat_out, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating rzmat_out in dis_extract_gamma', comm)
      return
    endif
    deallocate (rzmat_in, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating rzmat_in in dis_extract_gamma', comm)
      return
    endif

    allocate (ceamp(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating ceamp in dis_extract_gamma', comm)
      return
    endif
    allocate (cham(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cham in dis_extract_gamma', comm)
      return
    endif

    if (.not. dis_converged) then
      write (stdout, '(/5x,a)') &
        '<<< Warning: Maximum number of disentanglement iterations reached >>>'
      write (stdout, '(10x,a)') '<<< Disentanglement convergence criteria not satisfied >>>'
    endif

    if (icompflag .eq. 1) then
      if (print_output%iprint > 2) then
        write (stdout, ('(/4x,a)')) &
          'WARNING: Complement subspace has zero dimensions at the following k-points:'
        i = 0
        write (stdout, '(4x)', advance='no')
        do nkp = 1, num_kpts
          if (dis_manifold%ndimwin(nkp) .eq. num_wann) then
            i = i + 1
            if (i .le. 12) then
              write (stdout, '(i6)', advance='no') nkp
            else
              i = 1
              write (stdout, '(/4x)', advance='no')
              write (stdout, '(i6)', advance='no') nkp
            endif
          endif
        enddo
      endif
    endif

    ! Write the final womegai. This should remain unchanged during the
    ! subsequent minimization of Omega_tilde in wannierise.f90
    ! We store it in the checkpoint file as a sanity check
    write (stdout, '(/8x,a,f14.8,a/)') 'Final Omega_I ', &
      womegai*print_output%lenconfac**2, ' ('//trim(print_output%length_unit)//'^2)'

    ! Set public variable omega_invariant
    omega_invariant = womegai

    ! DIAGONALIZE THE HAMILTONIAN WITHIN THE OPTIMIZED SUBSPACES
    do nkp = 1, num_kpts

      do j = 1, num_wann
        do i = 1, num_wann
          cham(i, j, nkp) = cmplx_0
          do l = 1, dis_manifold%ndimwin(nkp)
            cham(i, j, nkp) = cham(i, j, nkp) + conjg(u_matrix_opt(l, i, nkp)) &
                              *u_matrix_opt(l, j, nkp)*eigval_opt(l, nkp)
          enddo
        enddo
      enddo
      do j = 1, num_wann
        do i = 1, j
          cap_r(i + ((j - 1)*j)/2) = real(cham(i, j, nkp), dp)
        enddo
      enddo

      call DSPEVX('V', 'A', 'U', num_wann, cap_r, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                  m, w, rz, num_bands, work, iwork, ifail, info)

      if (info .lt. 0) then
        write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
        write (stdout, *) ' THE ', -info, ' ARGUMENT OF DSPEVX HAD AN ILLEGAL VALUE'
        call set_error_fatal(error, ' dis_extract_gamma: error', comm)
        return
      endif
      if (info .gt. 0) then
        write (stdout, *) ' *** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
        write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
        call set_error_fatal(error, ' dis_extract_gamma: error', comm)
        return
      endif

      cz = cmplx_0
      cz(1:num_wann, 1:num_wann) = cmplx(rz(1:num_wann, 1:num_wann), 0.0_dp, dp)

      ! Store the energy eigenvalues of the optimal subspace (used in wann_ban
      eigval_opt(1:num_wann, nkp) = w(1:num_wann)

      ! CALCULATE AMPLITUDES OF THE CORRESPONDING ENERGY EIGENVECTORS IN TERMS
      ! THE ORIGINAL ("WINDOW SPACE") ENERGY EIGENVECTORS
      do j = 1, num_wann
        do i = 1, dis_manifold%ndimwin(nkp)
          ceamp(i, j, nkp) = cmplx_0
          do l = 1, num_wann
            ceamp(i, j, nkp) = ceamp(i, j, nkp) + cz(l, j)*u_matrix_opt(i, l, nkp)
          enddo
        enddo
      enddo
      ! NKP
    enddo
    ! DEBUG
    if (print_output%iprint > 2) then
      write (stdout, '(/,a,/)') '  Eigenvalues inside optimal subspace:'
      do nkp = 1, num_kpts
        write (stdout, '(a,i3,2x,20(f9.5,1x))') '  K-point ', &
          nkp, (eigval_opt(i, nkp), i=1, num_wann)
      enddo
    endif
    ! ENDDEBUG

    ! Replace u_matrix_opt by ceamp. Both span the
    ! same space, but the latter is more convenient for the purpose of obtai
    ! an optimal Fourier-interpolated band structure: see Sec. III.E of SMV.
    do nkp = 1, num_kpts
      do j = 1, num_wann
        u_matrix_opt(1:dis_manifold%ndimwin(nkp), j, nkp) = ceamp(1:dis_manifold%ndimwin(nkp), j, nkp)
      enddo
    enddo

    ! aam: 01/05/2009: added devel_flag if statement as the complementary
    !      subspace code was causing catastrophic seg-faults
    !if (index(print_output%devel_flag, 'compspace') > 0) then

    ! The compliment subspace code needs work: jry
    !  if (icompflag .eq. 1) then
    !    if (print_output%iprint > 2) then
    !      write (stdout, *) 'AT SOME K-POINT(S) COMPLEMENT SUBSPACE HAS ZERO DIMENSIONALITY'
    !      write (stdout, *) '=> DID NOT CREATE FILE COMPSPACE.DAT'
    !    endif
    !  else
    ! DIAGONALIZE THE HAMILTONIAN IN THE COMPLEMENT SUBSPACE, WRITE THE
    ! CORRESPONDING EIGENFUNCTIONS AND ENERGY EIGENVALUES
    !    do nkp = 1, num_kpts
    !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
    !        do i = 1, dis_manifold%ndimwin(nkp) - num_wann
    !          cham(i, j, nkp) = cmplx_0
    !          do l = 1, dis_manifold%ndimwin(nkp)
    !            cham(i, j, nkp) = cham(i, j, nkp) + conjg(camp(l, i, nkp)) &
    !                              *camp(l, j, nkp)*eigval_opt(l, nkp)
    !          enddo
    !        enddo
    !      enddo
    !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
    !        do i = 1, j
    !          cap_r(i + ((j - 1)*j)/2) = real(cham(i, j, nkp), dp)
    !        enddo
    !      enddo
    !      ndiff = dis_manifold%ndimwin(nkp) - num_wann

    !     call DSPEVX('V', 'A', 'U', ndiff, cap_r, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
    !                  m, w, rz, num_bands, work, iwork, ifail, info)

    !      if (info .lt. 0) then
    !        write (stdout, *) '*** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
    !        write (stdout, *) 'THE ', -info, ' ARGUMENT OF DSPEVX HAD AN ILLEGAL VALUE'
    !        call io_error(' dis_extract_gamma: error', stdout, seedname)
    !      endif
    !      if (info .gt. 0) then
    !        write (stdout, *) '*** ERROR *** DSPEVX WHILE DIAGONALIZING HAMILTONIAN'
    !        write (stdout, *) info, 'EIGENVECTORS FAILED TO CONVERGE'
    !        call io_error(' dis_extract_gamma: error', stdout, seedname)
    !      endif

    !      cz = cmplx_0
    !      cz(1:ndiff, 1:ndiff) = cmplx(rz(1:ndiff, 1:ndiff), 0.0_dp, dp)

    ! CALCULATE AMPLITUDES OF THE ENERGY EIGENVECTORS IN THE COMPLEMENT SUBS
    ! TERMS OF THE ORIGINAL ENERGY EIGENVECTORS
    !      do j = 1, dis_manifold%ndimwin(nkp) - num_wann
    !        do i = 1, dis_manifold%ndimwin(nkp)
    !          camp(i, j, nkp) = cmplx_0
    !          do l = 1, dis_manifold%ndimwin(nkp) - num_wann
!write(stdout,*) 'i=',i,'   j=',j,'   l=',l
!write(stdout,*) '           camp(i,j,nkp)=',camp(i,j,nkp)
!write(stdout,*) '           cz(l,j)=',cz(l,j)
!write(stdout,*) '           u_matrix_opt(i,l,nkp)=',u_matrix_opt(i,l,nkp)

! aam: 20/10/2006 -- the second dimension of u_matrix_opt is out of bounds (allocated as num_wann)!
! commenting this line out.
!                     camp(i,j,nkp) = camp(i,j,nkp) + cz(l,j) * u_matrix_opt(i,l,nkp)
    !          enddo
    !        enddo
    !      enddo
    !    enddo
    ! [loop over k points (nkp)]

    !  endif
    ! [if icompflag=1]

    !endif
    ! [if index(devel_flag,'compspace')>0]

    deallocate (history, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating history in dis_extract_gamma', comm)
      return
    endif
    deallocate (cham, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cham in dis_extract_gamma', comm)
      return
    endif
    if (allocated(camp)) then
      deallocate (camp, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating camp in dis_extract_gamma', comm)
        return
      endif
    end if
    deallocate (ceamp, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating ceamp in dis_extract_gamma', comm)
      return
    endif
    deallocate (wkomegai1, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating wkomegai1 in dis_extract_gamma', comm)
      return
    endif
    deallocate (rz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating rz in dis_extract_gamma', comm)
      return
    endif
    deallocate (cap_r, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cap_r in dis_extract_gamma', comm)
      return
    endif
    deallocate (work, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating work in dis_extract_gamma', comm)
      return
    endif
    deallocate (cz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cz in dis_extract_gamma', comm)
      return
    endif
    deallocate (w, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating w in dis_extract_gamma', comm)
      return
    endif
    deallocate (ifail, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating ifail in dis_extract_gamma', comm)
      return
    endif
    deallocate (iwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating iwork in dis_extract_gamma', comm)
      return
    endif
    deallocate (cbw, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cbw in dis_extract_gamma', comm)
      return
    endif
    deallocate (cww, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cww in dis_extract_gamma', comm)
      return
    endif
    deallocate (cwb, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating cwb in dis_extract_gamma', comm)
      return
    endif

    write (stdout, '(1x,a/)') &
      '+----------------------------------------------------------------------------+'

    if (print_output%timing_level > 1) call io_stopwatch_stop('dis: extract_gamma', timer)

    return
    !================================================!
  end subroutine dis_extract_gamma

  subroutine internal_zmatrix_gamma(cbw, m_matrix_orig_local, u_matrix_opt, rmtrx, wb, indxnfroz, &
                                    ndimfroz, ndimwin, nnlist, nkp, nntot, num_bands, num_wann, &
                                    timing_level, timer)
    !================================================!
    !
    !! Compute Z-matrix (Gamma point routine)
    !
    !================================================!
    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer

    integer, intent(in) :: timing_level
    integer, intent(in) :: num_bands, num_wann, nkp, nntot
    integer, intent(in) :: ndimwin(:)
    integer, intent(in) :: nnlist(:, :)
    integer, intent(in) :: ndimfroz(:)
    integer, intent(in) :: indxnfroz(:, :)

    real(kind=dp), intent(in) :: wb(:)
    real(kind=dp), intent(out) :: rmtrx(:, :)

    complex(kind=dp), intent(in) :: cbw(:, :)
    complex(kind=dp), intent(in) :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    ! Internal variables
    integer :: l, m, n, p, q, nn, nkp2, ndimk
    complex(kind=dp) :: csum

    if (timing_level > 1) call io_stopwatch_start('dis: extract_gamma: zmatrix_gamma', timer)

    rmtrx = 0.0_dp
    ndimk = ndimwin(nkp) - ndimfroz(nkp)
    do nn = 1, nntot
      nkp2 = nnlist(nkp, nn)
      call zgemm('N', 'N', num_bands, num_wann, ndimwin(nkp2), cmplx_1, &
                 m_matrix_orig_local(:, :, nn, nkp), num_bands, u_matrix_opt(:, :, nkp2), num_bands, &
                 cmplx_0, cbw, num_bands)
      do n = 1, ndimk
        q = indxnfroz(n, nkp)
        do m = 1, n
          p = indxnfroz(m, nkp)
          csum = cmplx_0
          do l = 1, num_wann
            csum = csum + cbw(p, l)*conjg(cbw(q, l))
          enddo
          rmtrx(m, n) = rmtrx(m, n) + wb(nn)*real(csum, dp)
          rmtrx(n, m) = rmtrx(m, n)
        enddo
      enddo
    enddo

    if (timing_level > 1) call io_stopwatch_stop('dis: extract_gamma: zmatrix_gamma', timer)

    return
    !================================================!
  end subroutine internal_zmatrix_gamma

end module w90_disentangle_mod
