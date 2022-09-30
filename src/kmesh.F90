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
!  w90_kmesh: operations on BZ mesh                          !
!                                                            !
!------------------------------------------------------------!

module w90_kmesh

  !! Routines to analyse the regular k-point mesh
  !! and determine the overlaps neccessary for a finite
  !! difference representation of the spread operator.
  !! These overlaps are defined by a set of vectors (b-vectors) which
  !! connect the Bloch states.
  !! See Eq. B1 in Appendix B of Marzari and
  !!  Vanderbilt  PRB 56 12847 (1997)

  use w90_constants, only: dp
  use w90_types, only: max_shells, num_nnmax ! used for dimensioning
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
    set_error_input, set_error_file
  use w90_comms, only: w90_comm_type

  implicit none

  private

  ! Definitions of system variables set in this module
  ! nnh     ! the number of b-directions (bka)
  ! nntot   ! total number of neighbours for each k-point
  ! nnlist  ! list of neighbours for each k-point
  ! neigh
  ! nncell  ! gives BZ of each neighbour of each k-point
  ! wbtot
  ! wb      ! weights associated with neighbours of each k-point
  ! bk      ! the b-vectors that go from each k-point to its neighbours
  ! bka     ! the b-directions (not considering inversion) from
  ! 1st k-point to its neighbours

  public :: kmesh_dealloc
  public :: kmesh_get
  public :: kmesh_write

  integer, parameter :: nsupcell = 5
  !! Size of supercell (of recip cell) in which to search for k-point shells

contains

  !================================================
  subroutine kmesh_get(kmesh_input, kmesh_info, print_output, kpt_latt, real_lattice, num_kpts, &
                       gamma_only, stdout, timer, error, comm)
    !================================================
    !
    !! Main routine to calculate the b-vectors
    !
    !================================================

    use w90_constants, only: eps6
    use w90_utility, only: utility_compar, utility_recip_lattice, utility_frac_to_cart, &
      utility_cart_to_frac, utility_inverse_mat
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_info_type, kmesh_input_type, print_output_type, timer_list_type

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_kpts
    integer, intent(in) :: stdout
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    logical, intent(in) :: gamma_only

    ! local variables
    real(kind=dp), allocatable :: bvec_tmp(:, :)
    real(kind=dp), allocatable :: kpt_cart(:, :)
    real(kind=dp), allocatable :: bk_local(:, :, :)
    real(kind=dp), parameter :: eta = 99999999.0_dp    ! eta = very large
    real(kind=dp) :: dist, dnn0, dnn1, bb1, bbn, ddelta
    real(kind=dp) :: dnn(max(kmesh_input%search_shells, 6*kmesh_input%higher_order_n))
    real(kind=dp) :: recip_lattice(3, 3), volume
    real(kind=dp) :: vkpp(3), vkpp2(3)

    ! higher-order finite-difference
    integer, allocatable :: lmn(:, :) ! Order in which to search the cells (ordered in dist from origin)
    integer, allocatable :: nnlist_tmp(:, :), nncell_tmp(:, :, :) ![ysl]
    integer, allocatable :: nnshell(:, :)
    integer :: ifound, counter, na, nap, loop_s, loop_b, shell !, nbvec, bnum
    integer :: ifpos, ifneg, ierr, multi(max(kmesh_input%search_shells, 6*kmesh_input%higher_order_n))
    integer :: nlist, nkp, nkp2, l, m, n, ndnn, ndnnx, ndnntot
    integer :: nnsh, nn, nnx, loop, i, j
    integer :: num_first_shells, ndnn2, nnx2, multi_cumulative, lmn_temp(3)
    integer :: num_x((1 + kmesh_input%higher_order_n)*(1 + 2*kmesh_input%higher_order_n))
    integer :: num_y((1 + kmesh_input%higher_order_n)*(1 + 2*kmesh_input%higher_order_n))
    integer :: num_z((1 + kmesh_input%higher_order_n)*(1 + 2*kmesh_input%higher_order_n))
    real(kind=dp) :: bk_latt(3), inv_lattice(3, 3)
    real(kind=dp) :: bweight(kmesh_input%max_shells_h)
    real(kind=dp) :: wb_local(kmesh_input%num_nnmax_h)

    if (print_output%timing_level > 0) call io_stopwatch_start('kmesh: get', timer)

    allocate (bk_local(3, kmesh_input%num_nnmax_h, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating bk_local in kmesh_get', comm)
      return
    endif

    allocate (lmn(3, (2*kmesh_input%search_supcell_size + 1)**3), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating lmn in kmesh_get', comm)
      return
    endif

    allocate (nnshell(num_kpts, max(kmesh_input%search_shells, 6*kmesh_input%higher_order_n)), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating nnshell in kmesh_get', comm)
      return
    endif

    call utility_recip_lattice(real_lattice, recip_lattice, volume, error, comm)
    call utility_inverse_mat(recip_lattice, inv_lattice)
    if (print_output%iprint > 0) write (stdout, '(/1x,a)') &
      '*---------------------------------- K-MESH ----------------------------------*'

    ! Sort the cell neighbours so we loop in order of distance from the home shell
    call kmesh_supercell_sort(print_output, recip_lattice, lmn, timer)

    allocate (kpt_cart(3, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating kpt_cart in kmesh_get', comm)
      return
    endif
    do nkp = 1, num_kpts
      call utility_frac_to_cart(kpt_latt(:, nkp), kpt_cart(:, nkp), recip_lattice)
    enddo

    ! find the distance between k-point 1 and its nearest-neighbour shells
    ! if we have only one k-point, the n-neighbours are its periodic images

    dnn0 = 0.0_dp
    dnn1 = eta
    ndnntot = 0
    do nlist = 1, kmesh_input%search_shells
      do nkp = 1, num_kpts
        do loop = 1, (2*kmesh_input%search_supcell_size + 1)**3
          l = lmn(1, loop); m = lmn(2, loop); n = lmn(3, loop)
          !
          vkpp = kpt_cart(:, nkp) + matmul(lmn(:, loop), recip_lattice)
          dist = sqrt((kpt_cart(1, 1) - vkpp(1))**2 &
                      + (kpt_cart(2, 1) - vkpp(2))**2 + (kpt_cart(3, 1) - vkpp(3))**2)
          if ((dist .gt. kmesh_input%tol) .and. (dist .gt. dnn0 + kmesh_input%tol)) then
            if (dist .lt. dnn1 - kmesh_input%tol) then
              dnn1 = dist  ! found a closer shell
              counter = 0
            end if
            if (dist .gt. (dnn1 - kmesh_input%tol) .and. dist .lt. (dnn1 + kmesh_input%tol)) then
              counter = counter + 1 ! count the multiplicity of the shell
            end if
          end if
        enddo
      enddo
      if (dnn1 .lt. eta - kmesh_input%tol) ndnntot = ndnntot + 1
      dnn(nlist) = dnn1
      multi(nlist) = counter
      dnn0 = dnn1
      dnn1 = eta
    enddo

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a)') '|                    Distance to Nearest-Neighbour Shells                    |'
      write (stdout, '(1x,a)') '|                    ------------------------------------                    |'
      if (print_output%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|          Shell             Distance (Ang^-1)          Multiplicity         |'
        write (stdout, '(1x,a)') '|          -----             -----------------          ------------         |'
      else
        write (stdout, '(1x,a)') '|          Shell             Distance (Bohr^-1)         Multiplicity         |'
        write (stdout, '(1x,a)') '|          -----             ------------------         ------------         |'
      endif
      do ndnn = 1, ndnntot
        write (stdout, '(1x,a,11x,i3,17x,f10.6,19x,i4,12x,a)') '|', ndnn, &
          dnn(ndnn)/print_output%lenconfac, multi(ndnn), '|'
      enddo
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    endif

    if (print_output%iprint >= 4) then
      ! Write out all the bvectors
      if (print_output%iprint > 0) then
        write (stdout, '(1x,"|",76(" "),"|")')
        write (stdout, '(1x,a)') '|         Complete list of b-vectors and their lengths                       |'
        write (stdout, '(1x,"|",76(" "),"|")')
        write (stdout, '(1x,"+",76("-"),"+")')
      endif

      allocate (bvec_tmp(3, maxval(multi)), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating bvec_tmp in kmesh_get', comm)
        return
      endif
      bvec_tmp = 0.0_dp
      counter = 0
      do shell = 1, kmesh_input%search_shells
        call kmesh_get_bvectors(kmesh_input, print_output, bvec_tmp(:, 1:multi(shell)), kpt_cart, &
                                recip_lattice, dnn(shell), lmn, 1, multi(shell), num_kpts, &
                                timer, error, comm)
        if (allocated(error)) return

        do loop = 1, multi(shell)
          counter = counter + 1
          if (print_output%iprint > 0) write (stdout, '(a,I4,1x,a,2x,3f12.6,2x,a,2x,f12.6,a)') ' | b-vector  ', counter, ': (', &
            bvec_tmp(:, loop)/print_output%lenconfac, ')', dnn(shell)/print_output%lenconfac, '  |'
        end do
      end do
      deallocate (bvec_tmp)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating bvec_tmp in kmesh_get', comm)
        return
      endif
      if (print_output%iprint > 0) write (stdout, '(1x,"|",76(" "),"|")')
      if (print_output%iprint > 0) write (stdout, '(1x,"+",76("-"),"+")')
    end if

    ! Get the shell weights to satisfy the B1 condition
    !if (index(print_output%devel_flag, 'kmesh_degen') > 0) then
    !  call kmesh_shell_from_file(kmesh_input, print_output, bvec_inp, bweight, dnn, kpt_cart, &
    !                             recip_lattice, lmn, multi, num_kpts, seedname, stdout)
    !else
    if (kmesh_input%num_shells == 0) then
      call kmesh_shell_automatic(kmesh_input, print_output, bweight, dnn, kpt_cart, recip_lattice, &
                                 lmn, multi, num_kpts, stdout, timer, error, comm)
      if (allocated(error)) return

    elseif (kmesh_input%num_shells > 0) then
      call kmesh_shell_fixed(kmesh_input, print_output, bweight, dnn, kpt_cart, recip_lattice, &
                             lmn, multi, num_kpts, stdout, timer, error, comm)
      if (allocated(error)) return
    end if

    num_first_shells = kmesh_input%num_shells ! for convenience in printout error

    if (print_output%iprint > 0) then
      if (kmesh_input%higher_order_nearest_shells) then
        write (stdout, '(1x,a)', advance='no') '| The following shells and their multiples are used: '
      else
        write (stdout, '(1x,a)', advance='no') '| The following shells are used:                     '
      endif
      do ndnn = 1, kmesh_input%num_shells
        if (ndnn .eq. kmesh_input%num_shells) then
          write (stdout, '(i3,1x)', advance='no') kmesh_input%shell_list(ndnn)
        else
          write (stdout, '(i3,",")', advance='no') kmesh_input%shell_list(ndnn)
        endif
      enddo
      do l = 1, 6 - kmesh_input%num_shells
        write (stdout, '(4x)', advance='no')
      enddo
      if (kmesh_input%higher_order_nearest_shells) then
        write (stdout, '(20x)', advance='no')
      endif
      write (stdout, '("|")')
    endif
    !end if

    kmesh_info%nntot = 0
    do loop_s = 1, kmesh_input%num_shells
      kmesh_info%nntot = kmesh_info%nntot + multi(kmesh_input%shell_list(loop_s))
    end do

    if (kmesh_info%nntot > kmesh_input%num_nnmax_h) then
      if (print_output%iprint > 0) then
        write (stdout, '(a,i2,a)') ' **WARNING: kmesh has found >', kmesh_input%num_nnmax_h, ' nearest neighbours**'
        write (stdout, '(a)') ' '
        write (stdout, '(a)') ' This is probably caused by an error in your unit cell specification'
        write (stdout, '(a)') ' '
        write (stdout, '(a)') ' If you think this is not the problem; please send your *.win file to the '
        write (stdout, '(a)') ' wannier90 developers'
        write (stdout, '(a)') ' '
        write (stdout, '(a)') ' The problem may be caused by having accidentally degenerate shells of '
        write (stdout, '(a)') ' kpoints. The solution is then to rerun wannier90 specifying the b-vectors '
        write (stdout, '(a)') ' in each shell.  Give devel_flag=kmesh_degen in the *.win file'
        write (stdout, '(a)') ' and create a *.kshell file:'
        write (stdout, '(a)') ' '
        write (stdout, '(a)') ' $>   cat hexagonal.kshell'
        write (stdout, '(a)') ' $>   1 2'
        write (stdout, '(a)') ' $>   5 6 7 8'
        write (stdout, '(a)') ' '
        write (stdout, '(a)') ' Where each line is a new shell (so num_shells in total)'
        write (stdout, '(a)') ' The elements are the bvectors labelled according to the following '
        write (stdout, '(a)') ' list (last column is distance)'
        write (stdout, '(a)') ' '
      endif

      allocate (bvec_tmp(3, maxval(multi)), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating bvec_tmp in kmesh_get', comm)
        return
      endif
      bvec_tmp = 0.0_dp
      counter = 0
      do shell = 1, kmesh_input%search_shells
        call kmesh_get_bvectors(kmesh_input, print_output, bvec_tmp(:, 1:multi(shell)), kpt_cart, &
                                recip_lattice, dnn(shell), lmn, 1, multi(shell), num_kpts, timer, &
                                error, comm)
        if (allocated(error)) return
        do loop = 1, multi(shell)
          counter = counter + 1
          if (print_output%iprint > 0) write (stdout, '(a,I4,1x,a,2x,3f12.6,2x,a,2x,f12.6,a)') ' | b-vector  ', counter, ': (', &
            bvec_tmp(:, loop)/print_output%lenconfac, ')', dnn(shell)/print_output%lenconfac, '  |'
        end do
      end do
      if (print_output%iprint > 0) write (stdout, '(a)') ' '
      deallocate (bvec_tmp)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating bvec_tmp in kmesh_get', comm)
        return
      endif

      call set_error_fatal(error, 'kmesh_get: something wrong, found too many nearest neighbours', comm)
      return
    end if

    ! higher-order algorithm: include 2b, 3b, ..., Nb shells, and modify bweights
    if (kmesh_input%higher_order_nearest_shells) then
      write (stdout, '(a)') ' | WARNING: higher_order_nearest_shells is an experimental feature, and has   |', &
        ' | not been extensively tested.                                               |'
    else
      ! update nntot
      kmesh_info%nntot = kmesh_info%nntot*kmesh_input%higher_order_n
    endif

    allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating nnlist in kmesh_get', comm)
      return
    endif
    allocate (kmesh_info%neigh(num_kpts, kmesh_info%nntot/2), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating neigh in kmesh_get', comm)
      return
    endif
    allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating nncell in kmesh_get', comm)
      return
    endif

    allocate (kmesh_info%wb(kmesh_info%nntot), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wb in kmesh_get', comm)
      return
    endif
    allocate (kmesh_info%bka(3, kmesh_info%nntot/2), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating bka in kmesh_get', comm)
      return
    endif
    allocate (kmesh_info%bk(3, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating bk in kmesh_get', comm)
      return
    endif

    nnx = 0
    do loop_s = 1, kmesh_input%num_shells
      do loop_b = 1, multi(kmesh_input%shell_list(loop_s))
        nnx = nnx + 1
        wb_local(nnx) = bweight(loop_s)
      end do
    end do

    ! Now build up the list of nearest-neighbour shells for each k-point.
    ! nnlist(nkp,1...nnx) points to the nnx neighbours (ordered along increa
    ! shells) of the k-point nkp. nncell(i,nkp,nnth) tells us in which BZ is
    ! nnth nearest-neighbour of the k-point nkp. Construct the nnx b-vectors
    ! go from k-point nkp to each neighbour bk(1:3,nkp,1...nnx).
    ! Comment: Now we have bk(3,nntot,num_kps) 09/04/2006

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a)') '|                        Shell   # Nearest-Neighbours                        |'
      write (stdout, '(1x,a)') '|                        -----   --------------------                        |'
    endif
    !if (index(print_output%devel_flag, 'kmesh_degen') == 0) then
    !
    ! Standard routine
    !
    nnshell = 0
    do nkp = 1, num_kpts
      nnx = 0
      ok: do ndnnx = 1, kmesh_input%num_shells
        ndnn = kmesh_input%shell_list(ndnnx)
        do loop = 1, (2*kmesh_input%search_supcell_size + 1)**3
          l = lmn(1, loop); m = lmn(2, loop); n = lmn(3, loop)
          vkpp2 = matmul(lmn(:, loop), recip_lattice)
          do nkp2 = 1, num_kpts
            vkpp = vkpp2 + kpt_cart(:, nkp2)
            dist = sqrt((kpt_cart(1, nkp) - vkpp(1))**2 &
                        + (kpt_cart(2, nkp) - vkpp(2))**2 + (kpt_cart(3, nkp) - vkpp(3))**2)
            if ((dist .ge. dnn(ndnn)*(1 - kmesh_input%tol)) .and. (dist .le. dnn(ndnn)*(1 + kmesh_input%tol))) then
              nnx = nnx + 1
              nnshell(nkp, ndnn) = nnshell(nkp, ndnn) + 1
              kmesh_info%nnlist(nkp, nnx) = nkp2
              kmesh_info%nncell(1, nkp, nnx) = l
              kmesh_info%nncell(2, nkp, nnx) = m
              kmesh_info%nncell(3, nkp, nnx) = n
              bk_local(:, nnx, nkp) = vkpp(:) - kpt_cart(:, nkp)
            endif
            !if we have the right number of neighbours we can exit
            if (nnshell(nkp, ndnn) == multi(ndnn)) cycle ok
          enddo
        enddo
        ! check to see if too few neighbours here
      end do ok
    end do

    ! higher-order algorithm: include 2b, 3b, ..., Nb shells, and modify bweights
    if (.not. kmesh_input%higher_order_nearest_shells) then
      ! update num_shells, shell_list, dnn(distance to shells), multi, etc.
      call kmesh_shell_reconstruct(kmesh_input, num_kpts, multi, dnn, nnshell, bweight)

      ! update bk_local
      ! update nnlist(neighboring kpts), and nncell(G_lmn vectors of neighbors)
      do nkp = 1, num_kpts
        do nn = 2, kmesh_input%higher_order_n
          multi_cumulative = 0
          do ndnn = 1, kmesh_input%num_shells/kmesh_input%higher_order_n !first-order shells
            ! ndnn: index of shells (first order)
            ! ndnn2: index of shells (nn-th-order)
            ndnn2 = (nn - 1)*kmesh_input%num_shells/kmesh_input%higher_order_n + ndnn
            do nnx = 1 + multi_cumulative, multi(ndnn) + multi_cumulative
              counter = 0
              ! nnx: index of bvectors (first order)
              ! nnx2: index of bvectors (nn-th-order)
              nnx2 = (nn - 1)*kmesh_info%nntot/kmesh_input%higher_order_n + nnx
              bk_local(:, nnx2, nkp) = nn*bk_local(:, nnx, nkp)
              ! find nnlist and nncell
              !do loop = 1, (2*kmesh_input%search_supcell_size + 1)**3
              !  l = lmn(1, loop); m = lmn(2, loop); n = lmn(3, loop)
              !  vkpp2 = matmul(lmn(:, loop), recip_lattice) !G_lmn vector
              !  do nkp2 = 1, num_kpts
              !    vkpp = vkpp2 + kpt_cart(:, nkp2) - kpt_cart(:, nkp)
              !    ! if kp2 - kp1 == bk_local(:, nnx2, :)
              !    call utility_compar(vkpp(1), bk_local(1, nnx2, nkp), ifpos, ifneg)
              !    if (ifpos .eq. 1) then
              !      counter = counter + 1
              !      kmesh_info%nnlist(nkp, nnx2) = nkp2
              !      kmesh_info%nncell(1, nkp, nnx2) = l
              !      kmesh_info%nncell(2, nkp, nnx2) = m
              !      kmesh_info%nncell(3, nkp, nnx2) = n
              !    endif
              !  enddo
              !enddo
              ! do not search supcell
              ! find nnlist(nkp2) and nncell(lmn)
              call utility_cart_to_frac(bk_local(:, nnx2, nkp), bk_latt, inv_lattice)
              lmn_temp(1) = floor(kpt_latt(1, nkp) + bk_latt(1) + 1.e-6_dp) ! e.g. 3.999999999 is 4
              lmn_temp(2) = floor(kpt_latt(2, nkp) + bk_latt(2) + 1.e-6_dp)
              lmn_temp(3) = floor(kpt_latt(3, nkp) + bk_latt(3) + 1.e-6_dp)
              vkpp2 = matmul(lmn_temp, recip_lattice) !G_lmn vector
              vkpp = kpt_cart(:, nkp) + bk_local(:, nnx2, nkp) - vkpp2 ! k_2 = k_1 + Nb - G_lmn
              do nkp2 = 1, num_kpts
                call utility_compar(kpt_cart(:, nkp2), vkpp, ifpos, ifneg)
                if (ifpos .eq. 1) then
                  counter = counter + 1
                  kmesh_info%nnlist(nkp, nnx2) = nkp2
                  kmesh_info%nncell(1, nkp, nnx2) = lmn_temp(1)
                  kmesh_info%nncell(2, nkp, nnx2) = lmn_temp(2)
                  kmesh_info%nncell(3, nkp, nnx2) = lmn_temp(3)
                endif
              enddo
              if (counter == 0) then
                call set_error_fatal(error, 'Could not find Nb vectors', comm)
              endif
              if (counter >= 2) then
                call set_error_fatal(error, 'Error in kmesh_get, try to modify tolerance in utility_compar', comm)
              endif
            enddo
            multi_cumulative = multi_cumulative + multi(ndnn)
          enddo
        enddo
      enddo
    endif
    nnx = 0
    do loop_s = 1, kmesh_input%num_shells
      do loop_b = 1, multi(kmesh_input%shell_list(loop_s))
        nnx = nnx + 1
        wb_local(nnx) = bweight(loop_s)
      end do
    end do

    !else
    !
    ! incase we set the bvectors explicitly
    !
    !nnshell = 0
    !do nkp = 1, num_kpts
    !  nnx = 0
    !  ok2: do loop = 1, (2*kmesh_input%search_supcell_size + 1)**3
    !    l = lmn(1, loop); m = lmn(2, loop); n = lmn(3, loop)
    !    vkpp2 = matmul(lmn(:, loop), recip_lattice)
    !    do nkp2 = 1, num_kpts
    !      vkpp = vkpp2 + kpt_cart(:, nkp2)
    !      bnum = 0
    !      do ndnnx = 1, kmesh_input%num_shells
    !        do nbvec = 1, multi(ndnnx)
    !          bnum = bnum + 1
    !          kpbvec = kpt_cart(:, nkp) + bvec_inp(:, nbvec, ndnnx)
    !          dist = sqrt((kpbvec(1) - vkpp(1))**2 &
    !                      + (kpbvec(2) - vkpp(2))**2 + (kpbvec(3) - vkpp(3))**2)
    !          if (abs(dist) < kmesh_input%tol) then
    !            nnx = nnx + 1
    !            nnshell(nkp, ndnnx) = nnshell(nkp, ndnnx) + 1
    !            kmesh_info%nnlist(nkp, bnum) = nkp2
    !            kmesh_info%nncell(1, nkp, bnum) = l
    !            kmesh_info%nncell(2, nkp, bnum) = m
    !            kmesh_info%nncell(3, nkp, bnum) = n
    !            bk_local(:, bnum, nkp) = bvec_inp(:, nbvec, ndnnx)
    !          endif
    !        enddo
    !      end do
    !      if (nnx == sum(multi)) exit ok2
    !    end do
    !  enddo ok2
    ! check to see if too few neighbours here
    !end do

    !end if

    if (kmesh_input%higher_order_n .eq. 1 .or. kmesh_input%higher_order_nearest_shells) then
      do ndnnx = 1, kmesh_input%num_shells
        ndnn = kmesh_input%shell_list(ndnnx)
        if (print_output%iprint > 0) then
          write (stdout, '(1x,a,24x,i3,13x,i3,33x,a)') '|', ndnn, nnshell(1, ndnn), '|'
        endif
      end do
    else
      do ndnnx = 1, num_first_shells
        ndnn = kmesh_input%shell_list(ndnnx)
        if (print_output%iprint > 0) then
          write (stdout, '(1x,a,24x,i3,13x,i3,33x,a)') '|', ndnn, nnshell(1, ndnn), '|'
        endif
      end do
      do i = 2, kmesh_input%higher_order_n
        do ndnnx = 1, num_first_shells
          ndnn = kmesh_input%shell_list(ndnnx)
          if (print_output%iprint > 0) then
            write (stdout, '(1x,a,20x,i3,a,i2,13x,i3,33x,a)') '|', i, ' x', ndnn, nnshell(1, ndnn), '|'
          endif
        end do
      end do
    endif
    if (print_output%iprint > 0) write (stdout, '(1x,"+",76("-"),"+")')

    do nkp = 1, num_kpts
      nnx = 0
      do ndnnx = 1, kmesh_input%num_shells
        ndnn = kmesh_input%shell_list(ndnnx)
        do nnsh = 1, nnshell(nkp, ndnn)
          bb1 = 0.0_dp
          bbn = 0.0_dp
          nnx = nnx + 1
          do i = 1, 3
            bb1 = bb1 + bk_local(i, nnx, 1)*bk_local(i, nnx, 1)
            bbn = bbn + bk_local(i, nnx, nkp)*bk_local(i, nnx, nkp)
          enddo
          if (abs(sqrt(bb1) - sqrt(bbn)) .gt. kmesh_input%tol) then
            if (print_output%iprint > 0) write (stdout, '(1x,2f10.6)') bb1, bbn
            call set_error_fatal(error, 'Non-symmetric k-point neighbours!', comm)
            return
          endif
        enddo
      enddo
    enddo

    ! now check that the completeness relation is satisfied for every kpoint
    ! We know it is true for kpt=1; but we check the rest to be safe.
    ! Eq. B1 in Appendix B PRB 56 12847 (1997)

    if ((.not. kmesh_input%skip_B1_tests) .and. kmesh_input%higher_order_nearest_shells) then
      do nkp = 1, num_kpts
        do i = 1, kmesh_input%higher_order_n
          if ((.not. kmesh_input%higher_order_nearest_shells) .and. i > 1) exit
          do j = 1, (1 + i)*(1 + 2*i) ! multiset coefficient ((3, 2i)) = (1 + i)*(1 + 2*i), 3: x,y,z, 2i: num. of b
            ! separate cartesian components
            ! e.g. for i=1
            ! j=1: num_x=2, num_y=0, num_z=0
            ! j=2: num_x=1, num_y=1, num_z=0
            ! j=3: num_x=0, num_y=2, num_z=0
            ! j=4: num_x=1, num_y=0, num_z=1
            ! j=5: num_x=0, num_y=1, num_z=1
            ! j=6: num_x=0, num_y=0, num_z=2
            do l = 0, 2*i
              if ((2*i + 1)*l - l*(l - 1)/2 <= j - 1 &
                  .and. (2*i + 1)*(l + 1) - (l + 1)*l/2 > j - 1) then
                num_z(j) = l
                exit
              endif
            enddo
            num_y(j) = j - 1 - ((2*i + 1)*num_z(j) - num_z(j)*(num_z(j) - 1)/2)
            num_x(j) = 2*i - num_y(j) - num_z(j)
            ddelta = 0.0_dp
            nnx = 0
            do ndnnx = 1, kmesh_input%num_shells
              ndnn = kmesh_input%shell_list(ndnnx)
              do nnsh = 1, nnshell(1, ndnn)
                nnx = nnx + 1
                ddelta = ddelta + wb_local(nnx)*(bk_local(1, nnx, nkp)**num_x(j)) &
                         *(bk_local(2, nnx, nkp)**num_y(j))*(bk_local(3, nnx, nkp)**num_z(j))
              enddo
            enddo
            if (i .eq. 1 .and. (j .eq. 1 .or. j .eq. 3 .or. j .eq. 6)) then
              if (abs(ddelta - 1.0_dp) .gt. kmesh_input%tol) then
                if (print_output%iprint > 0) write (stdout, '(1x,3i3,f12.8)') num_x(j), num_y(j), num_z(j), ddelta
                call set_error_fatal(error, 'Eq. (B1) not satisfied in kmesh_get (1)', comm)
              endif
            else
              if (abs(ddelta) .gt. kmesh_input%tol) then
                if (print_output%iprint > 0) write (stdout, '(1x,3i3,f12.8)') num_x(j), num_y(j), num_z(j), ddelta
                call set_error_fatal(error, 'Eq. (B1) not satisfied in kmesh_get (2)', comm)
              endif
            end if
          enddo
        enddo
      enddo
    end if

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a)') '| Completeness relation is fully satisfied [Eq. (B1), PRB 56, 12847 (1997)]  |'
      if ((kmesh_input%higher_order_nearest_shells) .and. (kmesh_input%higher_order_n .gt. 1)) then
        write (stdout, '(1x,a)') '| Completeness relations for higher-order are fully satisfied                |'
      endif
      write (stdout, '(1x,"+",76("-"),"+")')
    endif

    !
    kmesh_info%wbtot = 0.0_dp
    nnx = 0
    do ndnnx = 1, kmesh_input%num_shells
      ndnn = kmesh_input%shell_list(ndnnx)
      do nnsh = 1, nnshell(1, ndnn)
        nnx = nnx + 1
        kmesh_info%wbtot = kmesh_info%wbtot + wb_local(nnx)
      enddo
    enddo

    kmesh_info%nnh = kmesh_info%nntot/2
    ! make list of bka vectors from neighbours of first k-point
    ! delete any inverse vectors as you collect them
    na = 0
    do nn = 1, kmesh_info%nntot
      ifound = 0
      if (na .ne. 0) then
        do nap = 1, na
          call utility_compar(kmesh_info%bka(1, nap), bk_local(1, nn, 1), ifpos, ifneg)
          if (ifneg .eq. 1) ifound = 1
        enddo
      endif
      if (ifound .eq. 0) then
        !         found new vector to add to set
        na = na + 1
        kmesh_info%bka(1, na) = bk_local(1, nn, 1)
        kmesh_info%bka(2, na) = bk_local(2, nn, 1)
        kmesh_info%bka(3, na) = bk_local(3, nn, 1)
      endif
    enddo
    if (na .ne. kmesh_info%nnh) then
      call set_error_fatal(error, 'Did not find right number of bk directions', comm)
      return
    endif

    if (print_output%iprint > 0) then
      if (print_output%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|                  b_k Vectors (Ang^-1) and Weights (Ang^2)                  |'
        write (stdout, '(1x,a)') '|                  ----------------------------------------                  |'
      else
        write (stdout, '(1x,a)') '|                 b_k Vectors (Bohr^-1) and Weights (Bohr^2)                 |'
        write (stdout, '(1x,a)') '|                 ------------------------------------------                 |'
      endif
      write (stdout, '(1x,a)') '|            No.         b_k(x)      b_k(y)      b_k(z)        w_b           |'
      write (stdout, '(1x,a)') '|            ---        --------------------------------     --------        |'
      do i = 1, kmesh_info%nntot
        write (stdout, '(1x,"|",11x,i3,5x,3f12.6,3x,f10.6,8x,"|")') &
          i, (bk_local(j, i, 1)/print_output%lenconfac, j=1, 3), wb_local(i)*print_output%lenconfac**2
      enddo
      write (stdout, '(1x,"+",76("-"),"+")')
      if (print_output%lenconfac .eq. 1.0_dp) then
        write (stdout, '(1x,a)') '|                           b_k Directions (Ang^-1)                          |'
        write (stdout, '(1x,a)') '|                           -----------------------                          |'
      else
        write (stdout, '(1x,a)') '|                           b_k Directions (Bohr^-1)                         |'
        write (stdout, '(1x,a)') '|                           ------------------------                         |'
      endif
      write (stdout, '(1x,a)') '|            No.           x           y           z                         |'
      write (stdout, '(1x,a)') '|            ---        --------------------------------                     |'
      do i = 1, kmesh_info%nnh
        write (stdout, '(1x,"|",11x,i3,5x,3f12.6,21x,"|")') i, (kmesh_info%bka(j, i)/print_output%lenconfac, j=1, 3)
      enddo
      write (stdout, '(1x,"+",76("-"),"+")')
      write (stdout, *) ' '
    endif

    ! find index array
    do nkp = 1, num_kpts
      do na = 1, kmesh_info%nnh
        ! first, zero the index array so we can check it gets filled
        kmesh_info%neigh(nkp, na) = 0
        ! now search through list of neighbours of this k-point
        do nn = 1, kmesh_info%nntot
          call utility_compar(kmesh_info%bka(1, na), bk_local(1, nn, nkp), ifpos, ifneg)
          if (ifpos .eq. 1) kmesh_info%neigh(nkp, na) = nn
        enddo
        ! check found
        if (kmesh_info%neigh(nkp, na) .eq. 0) then
          if (print_output%iprint > 0) write (stdout, *) ' nkp,na=', nkp, na
          call set_error_fatal(error, 'kmesh_get: failed to find neighbours for this kpoint', comm)
          return
        endif
      enddo
    enddo

    !fill in the global arrays from the local ones

    do loop = 1, kmesh_info%nntot
      kmesh_info%wb(loop) = wb_local(loop)
    end do

    do loop_s = 1, num_kpts
      do loop = 1, kmesh_info%nntot
        kmesh_info%bk(:, loop, loop_s) = bk_local(:, loop, loop_s)
      end do
    end do

![ysl-b]

    if (gamma_only) then
      ! use half of the b-vectors
      if (num_kpts .ne. 1) then
        call set_error_input(error, 'Error in kmesh_get: wrong choice of gamma_only option', comm)
        return
      endif

      ! reassign nnlist, nncell, wb, bk
      allocate (nnlist_tmp(num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating nnlist_tmp in kmesh_get', comm)
        return
      endif
      allocate (nncell_tmp(3, num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating nncell_tmp in kmesh_get', comm)
        return
      endif

      nnlist_tmp(:, :) = kmesh_info%nnlist(:, :)
      nncell_tmp(:, :, :) = kmesh_info%nncell(:, :, :)

      deallocate (kmesh_info%nnlist, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating nnlist in kmesh_get', comm)
        return
      endif
      deallocate (kmesh_info%nncell, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating nncell in kmesh_get', comm)
        return
      endif
      deallocate (kmesh_info%wb, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wb in kmesh_get', comm)
        return
      endif
      deallocate (kmesh_info%bk, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating bk in kmesh_get', comm)
        return
      endif

      kmesh_info%nntot = kmesh_info%nntot/2

      allocate (kmesh_info%nnlist(num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating nnlist in kmesh_get', comm)
        return
      endif
      allocate (kmesh_info%nncell(3, num_kpts, kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating nncell in kmesh_get', comm)
        return
      endif
      allocate (kmesh_info%wb(kmesh_info%nntot), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating wb in kmesh_get', comm)
        return
      endif
      allocate (kmesh_info%bk(3, kmesh_info%nntot, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating bk in kmesh_get', comm)
        return
      endif

      na = 0
      do nn = 1, 2*kmesh_info%nntot
        ifound = 0
        if (na .ne. 0) then
          do nap = 1, na
            call utility_compar(kmesh_info%bk(1, nap, 1), bk_local(1, nn, 1), ifpos, ifneg)
            if (ifneg .eq. 1) ifound = 1
          enddo
        endif
        if (ifound .eq. 0) then
          !         found new vector to add to set
          na = na + 1
          kmesh_info%bk(1, na, 1) = bk_local(1, nn, 1)
          kmesh_info%bk(2, na, 1) = bk_local(2, nn, 1)
          kmesh_info%bk(3, na, 1) = bk_local(3, nn, 1)
          kmesh_info%wb(na) = 2.0_dp*wb_local(nn)
          kmesh_info%nnlist(1, na) = nnlist_tmp(1, nn)
          kmesh_info%nncell(1, 1, na) = nncell_tmp(1, 1, nn)
          kmesh_info%nncell(2, 1, na) = nncell_tmp(2, 1, nn)
          kmesh_info%nncell(3, 1, na) = nncell_tmp(3, 1, nn)
          kmesh_info%neigh(1, na) = na
          ! check bk.eq.bka
          call utility_compar(kmesh_info%bk(1, na, 1), kmesh_info%bka(1, na), ifpos, ifneg)
          if (ifpos .ne. 1) then
            call set_error_input(error, 'Error in kmesh_get: bk is not identical to bka in gamma_only option', comm)
            return
          endif
        endif
      enddo

      if (na .ne. kmesh_info%nnh) then
        call set_error_fatal(error, 'Did not find right number of b-vectors in gamma_only option', comm)
        return
      endif

      if (print_output%iprint > 0) then
        write (stdout, '(1x,"+",76("-"),"+")')
        write (stdout, '(1x,a)') '|        Gamma-point: number of the b-vectors is reduced by half             |'
        write (stdout, '(1x,"+",76("-"),"+")')
        if (print_output%lenconfac .eq. 1.0_dp) then
          write (stdout, '(1x,a)') '|                  b_k Vectors (Ang^-1) and Weights (Ang^2)                  |'
          write (stdout, '(1x,a)') '|                  ----------------------------------------                  |'
        else
          write (stdout, '(1x,a)') '|                 b_k Vectors (Bohr^-1) and Weights (Bohr^2)                 |'
          write (stdout, '(1x,a)') '|                 ------------------------------------------                 |'
        endif
        write (stdout, '(1x,a)') '|            No.         b_k(x)      b_k(y)      b_k(z)        w_b           |'
        write (stdout, '(1x,a)') '|            ---        --------------------------------     --------        |'
        do i = 1, kmesh_info%nntot
          write (stdout, '(1x,"|",11x,i3,5x,3f12.6,3x,f10.6,8x,"|")') &
            i, (kmesh_info%bk(j, i, 1)/print_output%lenconfac, j=1, 3), kmesh_info%wb(i)*print_output%lenconfac**2
        enddo
        write (stdout, '(1x,"+",76("-"),"+")')
        write (stdout, *) ' '
      endif

      deallocate (nnlist_tmp, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating nnlist_tmp in kmesh_get', comm)
        return
      endif
      deallocate (nncell_tmp, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating nncell_tmp in kmesh_get', comm)
        return
      endif

    endif
![ysl-e]

    ! JJ, is use_ss_functional necessarily defined here, or must it be moved to "special"
    !if (wann_control%use_ss_functional) then
    !check allocations, please!!
    if (.not. gamma_only) then
      allocate (kmesh_info%nnord(kmesh_info%nntot, num_kpts))
      allocate (kmesh_info%nninv(kmesh_info%nntot, num_kpts))
      allocate (kmesh_info%nnrev(kmesh_info%nntot, num_kpts))
      call kmesh_bvectors_perm(kmesh_info%bk(:, :, :), kmesh_info%bk(:, :, 1), num_kpts, &
                               kmesh_info%nntot, kmesh_info%nnord, kmesh_info%nninv, &
                               kmesh_info%nnrev, error, comm)
    endif

    deallocate (kpt_cart, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating kpt_cart in kmesh_get', comm)
      return
    endif

    deallocate (bk_local, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating bk_local in kmesh_get', comm)
      return
    endif

    deallocate (lmn, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating lmn in kmesh_get', comm)
      return
    endif

    deallocate (nnshell, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating nnshell in kmesh_get', comm)
      return
    endif

    if (print_output%timing_level > 0) call io_stopwatch_stop('kmesh: get', timer)

    return

  end subroutine kmesh_get

  !================================================!
  subroutine kmesh_write(exclude_bands, kmesh_info, lauto_proj, proj, print_output, kpt_latt, &
                         real_lattice, num_kpts, num_proj, calc_only_A, spinors, seedname, timer)
    !==================================================================!
    !                                                                  !
    !! Writes nnkp file (list of overlaps needed)
    !                                                                  !
    ! Note that the format is different to (and more compact than)     !
    ! that used by the old f77 code.                                   !
    !                                                                  !
    ! The file consists of num_kpts blocks of data, one block for each !
    ! k-point of the mesh. Each block consists of nntot+1 lines,       !
    ! where nntot is the (integer) number of nearest neighbours        !
    ! belonging to k-point nkp.                                        !
    !                                                                  !
    ! The first line in each block is just nntot.                      !
    !                                                                  !
    ! The second line consists of 5 integers. The first is the k-point !
    ! nkp. The second to the fifth specify it's nearest neighbours     !
    ! k+b: the second integer points to the k-point that is the        !
    ! periodic image of k+b that we want; the last three integers give !
    ! the G-vector, in reciprocal lattice units, that brings the       !
    ! k-point specified by the second integer (which is in the first   !
    ! BZ) to the actual k+b that we need.                              !
    !                                                                  !
    ! So wannier.nnkp specifies the nearest neighbours of each         !
    ! k-point, and therefore provides the information required to      !
    ! calculate the M_mn(k,b) matrix elements -- Marzari & Vanderbilt  !
    ! PRB 56, 12847 (1997) Eq. (25) -- for each pair of band indices   !
    ! m and n.                                                         !
    !==================================================================!

    use w90_io, only: io_date, io_stopwatch_start, io_stopwatch_stop
    use w90_utility, only: utility_recip_lattice_base
    use w90_types, only: kmesh_info_type, kmesh_input_type, &
      proj_type, print_output_type, timer_list_type

    implicit none

    character(len=*), intent(in)  :: seedname
    integer, allocatable, intent(in) :: exclude_bands(:)
    integer, intent(in) :: num_kpts
    integer, intent(inout) :: num_proj
    logical, intent(in) :: calc_only_A
    logical, intent(in) :: spinors
    logical, intent(in) :: lauto_proj
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(proj_type), allocatable, intent(in) :: proj(:) ! alloc only because allocation status is tested
    type(timer_list_type), intent(inout) :: timer

    real(kind=dp) :: recip_lattice(3, 3), volume
    integer           :: i, nkp, nn, nnkpout, num_exclude_bands
    character(len=9) :: cdate, ctime

    if (print_output%timing_level > 0) call io_stopwatch_start('kmesh: write', timer)

    open (newunit=nnkpout, file=trim(seedname)//'.nnkp', form='formatted')

    ! Date and time
    call io_date(cdate, ctime)
    write (nnkpout, '(4(a),/)') '# File written on ', cdate, ' at ', ctime

    ! Calc_only_A
    write (nnkpout, '(a,l2,/)') 'calc_only_A  : ', calc_only_A

    ! Real lattice
    write (nnkpout, '(a)') 'begin real_lattice'
    write (nnkpout, '(3(f12.7))') (real_lattice(1, i), i=1, 3)
    write (nnkpout, '(3(f12.7))') (real_lattice(2, i), i=1, 3)
    write (nnkpout, '(3(f12.7))') (real_lattice(3, i), i=1, 3)
    write (nnkpout, '(a/)') 'end real_lattice'

    ! Reciprocal lattice
    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
    write (nnkpout, '(a)') 'begin recip_lattice'
    write (nnkpout, '(3f12.7)') (recip_lattice(1, i), i=1, 3)
    write (nnkpout, '(3f12.7)') (recip_lattice(2, i), i=1, 3)
    write (nnkpout, '(3f12.7)') (recip_lattice(3, i), i=1, 3)
    write (nnkpout, '(a/)') 'end recip_lattice'

    ! K-points
    write (nnkpout, '(a)') 'begin kpoints'
    write (nnkpout, '(i6)') num_kpts
    do nkp = 1, num_kpts
      write (nnkpout, '(3f14.8)') (kpt_latt(i, nkp), i=1, 3)
    enddo
    write (nnkpout, '(a/)') 'end kpoints'

    if (spinors) then
      ! Projections
      write (nnkpout, '(a)') 'begin spinor_projections'
      if (allocated(proj)) then
        write (nnkpout, '(i6)') num_proj
        do i = 1, num_proj
          write (nnkpout, '(3(f10.5,1x),2x,3i3)') &
            proj(i)%site(1), proj(i)%site(2), proj(i)%site(3), &
            proj(i)%l, proj(i)%m, proj(i)%radial
          write (nnkpout, '(2x,3f11.7,1x,3f11.7,1x,f7.2)') &
            proj(i)%z(1), proj(i)%z(2), proj(i)%z(3), &
            proj(i)%x(1), proj(i)%x(2), proj(i)%x(3), &
            proj(i)%zona
          write (nnkpout, '(2x,1i3,1x,3f11.7)') &
            proj(i)%s, &
            proj(i)%s_qaxis(1), proj(i)%s_qaxis(2), proj(i)%s_qaxis(3)
        enddo
      else
        ! No projections
        write (nnkpout, '(i6)') 0
      end if
      write (nnkpout, '(a/)') 'end spinor_projections'
    else
      ! Projections
      write (nnkpout, '(a)') 'begin projections'
      if (allocated(proj)) then
        write (nnkpout, '(i6)') num_proj
        do i = 1, num_proj
          write (nnkpout, '(3(f10.5,1x),2x,3i3)') &
            proj(i)%site(1), proj(i)%site(2), proj(i)%site(3), &
            proj(i)%l, proj(i)%m, proj(i)%radial
          write (nnkpout, '(2x,3f11.7,1x,3f11.7,1x,f7.2)') &
            proj(i)%z(1), proj(i)%z(2), proj(i)%z(3), &
            proj(i)%x(1), proj(i)%x(2), proj(i)%x(3), &
            proj(i)%zona
        enddo
      else
        ! No projections
        write (nnkpout, '(i6)') 0
      end if
      write (nnkpout, '(a/)') 'end projections'
    endif

    ! Info for automatic generation of projections
    if (lauto_proj) then
      write (nnkpout, '(a)') 'begin auto_projections'
      write (nnkpout, '(i6)') num_proj
      write (nnkpout, '(i6)') 0
      write (nnkpout, '(a/)') 'end auto_projections'
    end if

    ! Nearest neighbour k-points
    write (nnkpout, '(a)') 'begin nnkpts'
    write (nnkpout, '(i4)') kmesh_info%nntot
    do nkp = 1, num_kpts
      do nn = 1, kmesh_info%nntot
        write (nnkpout, '(2i8,3x,3i4)') &
          nkp, kmesh_info%nnlist(nkp, nn), (kmesh_info%nncell(i, nkp, nn), i=1, 3)
      end do
    end do
    write (nnkpout, '(a/)') 'end nnkpts'

    !states to exclude
    num_exclude_bands = 0
    if (allocated(exclude_bands)) num_exclude_bands = size(exclude_bands)
    write (nnkpout, '(a)') 'begin exclude_bands'
    write (nnkpout, '(i4)') num_exclude_bands
    if (num_exclude_bands > 0) then
      do i = 1, num_exclude_bands
        write (nnkpout, '(i4)') exclude_bands(i)
      end do
    endif
    write (nnkpout, '(a)') 'end exclude_bands'

    close (nnkpout)

    if (print_output%timing_level > 0) call io_stopwatch_stop('kmesh: write', timer)

    return

  end subroutine kmesh_write

  !================================================
  subroutine kmesh_dealloc(kmesh_info, error, comm)
    !================================================
    !!  Release memory from the kmesh module
    !   This routine now check to see if arrays
    !   are allocated, as there are some code
    !   paths that will not allocate on all nodes
    !================================================

    use w90_types, only: kmesh_info_type

    implicit none

    type(kmesh_info_type), intent(inout) :: kmesh_info
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    integer :: ierr

    ! Deallocate real arrays that are public
    if (allocated(kmesh_info%bk)) then
      deallocate (kmesh_info%bk, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating bk in kmesh_dealloc', comm)
        return
      endif
    endif
    if (allocated(kmesh_info%bka)) then
      deallocate (kmesh_info%bka, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating bka in kmesh_dealloc', comm)
        return
      endif
    endif
    if (allocated(kmesh_info%wb)) then
      deallocate (kmesh_info%wb, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating wb in kmesh_dealloc', comm)
        return
      endif
    end if

    ! Deallocate integer arrays that are public
    if (allocated(kmesh_info%neigh)) then
      deallocate (kmesh_info%neigh, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating neigh in kmesh_dealloc', comm)
        return
      endif
    end if
    if (allocated(kmesh_info%nncell)) then
      deallocate (kmesh_info%nncell, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating nncell in kmesh_dealloc', comm)
        return
      endif
    endif
    if (allocated(kmesh_info%nnlist)) then
      deallocate (kmesh_info%nnlist, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating nnlist in kmesh_dealloc', comm)
        return
      endif
    endif

    return

  end subroutine kmesh_dealloc

  !================================================
  subroutine kmesh_supercell_sort(print_output, recip_lattice, lmn, timer)
    !================================================
    !! We look for kpoint neighbours in a large supercell of reciprocal
    !! unit cells. Done sequentially this is very slow.
    !! Here we order the cells by the distance from the origin.
    !! Doing the search in this order gives a dramatic speed up
    !================================================

    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: print_output_type, timer_list_type

    implicit none

    type(print_output_type), intent(in) :: print_output
    integer, intent(inout) :: lmn(:, :)
    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    type(timer_list_type), intent(inout) :: timer

    integer :: counter, l, m, n, loop
    !! Order in which to search the cells (ordered in dist from origin)
    integer :: lmn_cp(3, (2*nsupcell + 1)**3), indx(1)
    real(kind=dp) :: pos(3)
    real(kind=dp) :: dist((2*nsupcell + 1)**3)
    real(kind=dp) :: dist_cp((2*nsupcell + 1)**3)

    if (print_output%timing_level > 1) call io_stopwatch_start('kmesh: supercell_sort', timer)

    counter = 1
    lmn(:, counter) = 0
    dist(counter) = 0.0_dp
    do l = -nsupcell, nsupcell
      do m = -nsupcell, nsupcell
        do n = -nsupcell, nsupcell
          if (l == 0 .and. m == 0 .and. n == 0) cycle
          counter = counter + 1
          lmn(1, counter) = l; lmn(2, counter) = m; lmn(3, counter) = n
          pos = matmul(lmn(:, counter), recip_lattice)
          dist(counter) = sqrt(dot_product(pos, pos))
        end do
      end do
    end do

    do loop = (2*nsupcell + 1)**3, 1, -1
      indx = internal_maxloc(dist, nsupcell)
      dist_cp(loop) = dist(indx(1))
      lmn_cp(:, loop) = lmn(:, indx(1))
      dist(indx(1)) = -1.0_dp
    end do

    lmn = lmn_cp
    dist = dist_cp

    if (print_output%timing_level > 1) call io_stopwatch_stop('kmesh: supercell_sort', timer)

  end subroutine kmesh_supercell_sort

  !================================================
  subroutine kmesh_get_bvectors(kmesh_input, print_output, bvector, kpt_cart, recip_lattice, &
                                shell_dist, lmn, kpt, multi, num_kpts, timer, error, comm)
    !================================================
    !
    !! Returns the b-vectors for a given shell and kpoint.
    !
    !================================================

    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_input_type, print_output_type, timer_list_type

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(kmesh_input_type), intent(in)  :: kmesh_input
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_kpts
    integer, intent(in) :: lmn(:, :)
    integer, intent(in) :: multi   ! the number of kpoints in the shell
    integer, intent(in) :: kpt     ! which kpt is our 'origin'

    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) ::kpt_cart(:, :)
    real(kind=dp), intent(in) :: shell_dist ! the bvectors
    real(kind=dp), intent(out) :: bvector(3, multi) ! the bvectors

    ! local variables
    integer :: loop, nkp2, num_bvec
    real(kind=dp) :: dist, vkpp2(3), vkpp(3)

    if (print_output%timing_level > 1) call io_stopwatch_start('kmesh: get_bvectors', timer)

    bvector = 0.0_dp

    num_bvec = 0
    ok: do loop = 1, (2*kmesh_input%search_supcell_size + 1)**3
      vkpp2 = matmul(lmn(:, loop), recip_lattice)
      do nkp2 = 1, num_kpts
        vkpp = vkpp2 + kpt_cart(:, nkp2)
        dist = sqrt((kpt_cart(1, kpt) - vkpp(1))**2 &
                    + (kpt_cart(2, kpt) - vkpp(2))**2 + (kpt_cart(3, kpt) - vkpp(3))**2)
        if ((dist .ge. shell_dist*(1.0_dp - kmesh_input%tol)) .and. dist .le. shell_dist*(1.0_dp + kmesh_input%tol)) then
          num_bvec = num_bvec + 1
          bvector(:, num_bvec) = vkpp(:) - kpt_cart(:, kpt)
        endif
        !if we have the right number of neighbours we can exit
        if (num_bvec == multi) cycle ok
      enddo
    enddo ok

    if (num_bvec < multi) then
      call set_error_fatal(error, 'kmesh_get_bvector: Not enough bvectors found', comm)
      return
    endif

    if (print_output%timing_level > 1) call io_stopwatch_stop('kmesh: get_bvectors', timer)

    return

  end subroutine kmesh_get_bvectors

  !================================================
  subroutine kmesh_shell_automatic(kmesh_input, print_output, bweight, dnn, kpt_cart, &
                                   recip_lattice, lmn, multi, num_kpts, stdout, timer, error, &
                                   comm)
    !================================================
    !! Find the correct set of shells to satisfy B1
    !!  The stratagy is:
    !!       1) Take the bvectors from the next shell
    !!       2) Reject them if they are parallel to exisiting b vectors
    !!       3) Test to see if we satisfy B1, if not add another shell and repeat
    !
    !================================================

    use w90_constants, only: eps5, eps6
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_input_type, print_output_type, timer_list_type

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_kpts
    integer, intent(in) :: stdout
    integer, intent(in) :: lmn(:, :)
    integer, intent(in) :: multi(kmesh_input%search_shells)   ! the number of kpoints in the shell

    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) :: kpt_cart(:, :)
    real(kind=dp), intent(in) :: dnn(kmesh_input%search_shells) ! the bvectors
    real(kind=dp), intent(out) :: bweight(kmesh_input%max_shells_h)

    ! local variables
    integer :: loop_bn, loop_b, loop_s, info, cur_shell, ierr, loop, shell
    logical :: lpar
    real(kind=dp), allocatable :: amat(:, :), umat(:, :), vmat(:, :), smat(:, :), tmp0(:, :)
    real(kind=dp), allocatable :: bvector(:, :, :) ! the bvectors
    real(kind=dp), allocatable :: singv(:), tmp1(:), tmp2(:), tmp3(:)
    real(kind=dp) :: delta
    real(kind=dp) :: target(kmesh_input%max_shells_aux)
    real(kind=dp) :: work((kmesh_input%max_shells_aux)*10)

    ! variables for higher-order finite-difference
    integer, dimension(:, :), allocatable :: num_x, num_y, num_z
    logical :: bsat
    integer :: loop_order, num_of_eqs, higher_order_n_local

    if (kmesh_input%higher_order_nearest_shells) then
      higher_order_n_local = kmesh_input%higher_order_n
    else
      higher_order_n_local = 1 !find 1st-order b and weights first in this subroutine
    endif
    target = 0.0_dp; target(1) = 1.0_dp; target(3) = 1.0_dp; target(6) = 1.0_dp

    if (print_output%timing_level > 1) call io_stopwatch_start('kmesh: shell_automatic', timer)
    allocate (bvector(3, maxval(multi), kmesh_input%max_shells_h), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating bvector in kmesh_shell_automatic', comm)
      return
    endif
    bvector = 0.0_dp; bweight = 0.0_dp

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a)') '| The b-vectors are chosen automatically                                     |'
    endif

    ! note allocation of kmesh_input%shell list in subroutine w90_readwrite_read_kmesh_data()
    ! kmesh_input%num_shells = 0 in same place
    bsat = .false.

    do shell = 1, kmesh_input%search_shells
      cur_shell = kmesh_input%num_shells + 1

      ! get the b vectors for the new shell
      call kmesh_get_bvectors(kmesh_input, print_output, bvector(:, 1:multi(shell), cur_shell), &
                              kpt_cart, recip_lattice, dnn(shell), lmn, 1, multi(shell), num_kpts, &
                              timer, error, comm)
      if (allocated(error)) return

      if (print_output%iprint >= 3) then
        write (stdout, '(1x,a8,1x,I2,a14,1x,I2,49x,a)') '| Shell:', shell, ' Multiplicity:', multi(shell), '|'
        do loop = 1, multi(shell)
          write (stdout, '(1x,a10,I2,1x,a1,4x,3f12.6,5x,a9,9x,a)') '| b-vector ', loop, ':', &
            bvector(:, loop, cur_shell)/print_output%lenconfac, '('//trim(print_output%length_unit)//'^-1)', '|'
        end do
      end if

      ! We check that the new shell is not parrallel to an existing shell (cosine=1)
      if (higher_order_n_local == 1) then
        lpar = .false.
        if (kmesh_input%num_shells > 0) then
          do loop_bn = 1, multi(shell)
            do loop_s = 1, kmesh_input%num_shells
              do loop_b = 1, multi(kmesh_input%shell_list(loop_s))
                delta = dot_product(bvector(:, loop_bn, cur_shell), bvector(:, loop_b, loop_s))/ &
                        sqrt(dot_product(bvector(:, loop_bn, cur_shell), bvector(:, loop_bn, cur_shell))* &
                             dot_product(bvector(:, loop_b, loop_s), bvector(:, loop_b, loop_s)))
                if (abs(abs(delta) - 1.0_dp) < eps6) lpar = .true.
              end do
            end do
          end do
        end if

        if (lpar) then
          if (print_output%iprint >= 3) then
            write (stdout, '(1x,a)') '| This shell is linearly dependent on existing shells: Trying next shell     |'
          end if
          cycle
        end if
      end if

      kmesh_input%num_shells = kmesh_input%num_shells + 1
      kmesh_input%shell_list(kmesh_input%num_shells) = shell

      allocate (tmp0(kmesh_input%max_shells_aux, kmesh_input%max_shells_aux), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating amat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (tmp1(kmesh_input%max_shells_aux), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating amat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (tmp2(kmesh_input%num_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating amat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (tmp3(kmesh_input%num_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating amat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (amat(kmesh_input%max_shells_aux, kmesh_input%num_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating amat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (umat(kmesh_input%max_shells_aux, kmesh_input%max_shells_aux), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating umat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (vmat(kmesh_input%num_shells, kmesh_input%num_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating vmat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (smat(kmesh_input%num_shells, kmesh_input%max_shells_aux), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating smat in kmesh_shell_automatic', comm)
        return
      endif
      allocate (singv(kmesh_input%num_shells), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating singv in kmesh_shell_automatic', comm)
        return
      endif
      amat(:, :) = 0.0_dp; umat(:, :) = 0.0_dp; vmat(:, :) = 0.0_dp; smat(:, :) = 0.0_dp; singv(:) = 0.0_dp

      num_of_eqs = (1 + higher_order_n_local)*(1 + 2*higher_order_n_local)
      allocate (num_x(higher_order_n_local, num_of_eqs), stat=ierr)
      if (ierr /= 0) call set_error_alloc(error, 'Error allocating num_x in kmesh_shell_automatic', comm)
      allocate (num_y(higher_order_n_local, num_of_eqs), stat=ierr)
      if (ierr /= 0) call set_error_alloc(error, 'Error allocating num_y in kmesh_shell_automatic', comm)
      allocate (num_z(higher_order_n_local, num_of_eqs), stat=ierr)
      if (ierr /= 0) call set_error_alloc(error, 'Error allocating num_z in kmesh_shell_automatic', comm)

      !find higher finite-diff weights
      ! make test suite(compare nnkp files)
      do loop_order = 1, higher_order_n_local
        call kmesh_get_amat(kmesh_input, amat, bvector, multi, loop_order, &
                            num_x(loop_order, :), num_y(loop_order, :), num_z(loop_order, :))
      enddo

      info = 0
      call dgesvd('A', 'A', kmesh_input%max_shells_aux, kmesh_input%num_shells, amat, &
                  kmesh_input%max_shells_aux, singv, umat, &
                  kmesh_input%max_shells_aux, vmat, kmesh_input%num_shells, work, kmesh_input%max_shells_aux*10, info)
      if (info < 0) then
        if (print_output%iprint > 0) then
          write (stdout, '(1x,a,1x,I1,1x,a)') 'kmesh_shell_automatic: Argument', abs(info), &
            'of dgesvd is incorrect'
        endif
        call set_error_fatal(error, 'kmesh_shell_automatic: Problem with Singular Value Decomposition', comm)
        return
      else if (info > 0) then
        call set_error_fatal(error, 'kmesh_shell_automatic: Singular Value Decomposition did not converge', comm)
        return
      end if

      if (any(abs(singv) < eps5)) then
        if (kmesh_input%num_shells == 1) then
          call set_error_fatal(error, &
                               'kmesh_shell_automatic: Singular Value Decomposition has found a very small singular value', comm)
          return
        else
          if (print_output%iprint > 0) then
            write (stdout, '(1x,a)') '| SVD found small singular value, Rejecting this shell and trying the next   |'
          endif
          bsat = .false.
          kmesh_input%num_shells = kmesh_input%num_shells - 1
          goto 200
        end if
      end if

      smat = 0.0_dp
      do loop_s = 1, kmesh_input%num_shells
        smat(loop_s, loop_s) = 1.0_dp/singv(loop_s)
      end do

      ! S. Ponce: The following below is correct but had to be unpacked because of PGI-15
      ! bweight(1:num_shells)=matmul(transpose(vmat),matmul(smat,matmul(transpose(umat),target)))
      tmp0 = transpose(umat)
      tmp1 = matmul(tmp0, target)
      tmp2 = matmul(smat, tmp1)
      tmp3 = matmul(transpose(vmat), tmp2)
      bweight(1:kmesh_input%num_shells) = tmp3

      if (print_output%iprint >= 2) then
        do loop_s = 1, kmesh_input%num_shells
          write (stdout, '(1x,a,I2,a,f12.7,5x,a8,36x,a)') '| Shell: ', loop_s, &
            ' w_b ', bweight(loop_s)*print_output%lenconfac**2, '('//trim(print_output%length_unit)//'^2)', '|'
        end do
      end if

      !check if the conditions including (B1) for finite-difference are satisfied
      bsat = .true.
      do loop_order = 1, higher_order_n_local
        call kmesh_check_condition(kmesh_input, bsat, bvector, bweight, multi, loop_order, &
                                   num_x(loop_order, :), num_y(loop_order, :), num_z(loop_order, :))
      end do

      if (.not. bsat) then
        if (shell < kmesh_input%search_shells .and. print_output%iprint >= 3) then
          if (print_output%iprint > 0) write (stdout, '(1x,a,24x,a1)') '| B1 condition is not satisfied: Adding another shell', '|'

        elseif (shell == kmesh_input%search_shells) then

          if (print_output%iprint > 0) then
            write (stdout, *) ' '
            write (stdout, '(1x,a,i3,a)') 'Unable to satisfy the higher-order version of B1 with any of the first ' &
              , kmesh_input%search_shells, ' shells'
            write (stdout, '(1x,a)') 'Check that you have specified your unit cell to a high precision'
            write (stdout, '(1x,a)') 'Low precision might cause a loss of symmetry.'
            write (stdout, '(1x,a)') ' '
            write (stdout, '(1x,a)') 'If your cell is very long, or you have an irregular MP grid'
            write (stdout, '(1x,a)') 'Try increasing the parameter search_shells in the win file (default=30)'
            write (stdout, *) ' '
            call set_error_fatal(error, 'kmesh_get_automatic', comm)
            return
          end if

        end if
      end if

200   continue

      deallocate (tmp0, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating amat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (tmp1, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating amat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (tmp2, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating amat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (tmp3, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating amat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (amat, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating amat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (umat, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating umat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (vmat, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating vmat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (smat, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating smat in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (singv, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating singv in kmesh_shell_automatic', comm)
        return
      endif

      deallocate (num_x, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating num_x in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (num_y, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating num_y in kmesh_shell_automatic', comm)
        return
      endif
      deallocate (num_z, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating num_z in kmesh_shell_automatic', comm)
        return
      endif

      if (bsat) exit

    end do

    if (.not. bsat) then
      if (print_output%iprint > 0) then
        write (stdout, *) ' '
        write (stdout, '(1x,a,i3,a)') 'Unable to satisfy B1 with any of the first ', &
          kmesh_input%search_shells, ' shells'
        write (stdout, '(1x,a)') 'Your cell might be very long, or you may have an irregular MP grid'
        write (stdout, '(1x,a)') 'Try increasing the parameter search_shells in the win file (default=36)'
        write (stdout, *) ' '
      end if
      call set_error_fatal(error, 'kmesh_get_automatic', comm)
      return
    end if

    if (print_output%timing_level > 1) call io_stopwatch_stop('kmesh: shell_automatic', timer)

    return

  end subroutine kmesh_shell_automatic

  subroutine kmesh_shell_reconstruct(kmesh_input, num_kpts, multi, dnn, nnshell, bweight)
    !================================================
    !
    !!  Include more shells to calculate higher-order finite difference: 2b, 3b, ... Nb shells
    !!  Note: some shells are overwritten
    !================================================
    use w90_types, only: kmesh_input_type, print_output_type

    implicit none

    ! arguments
    type(kmesh_input_type), intent(inout) :: kmesh_input
    integer, intent(in) :: num_kpts
    integer, intent(inout) :: multi(max(kmesh_input%search_shells, 6*kmesh_input%higher_order_n))   ! the number of bvectors in the shell
    real(kind=dp), intent(inout) :: dnn(max(kmesh_input%search_shells, 6*kmesh_input%higher_order_n))
    integer, intent(inout) :: nnshell(num_kpts, max(kmesh_input%search_shells, 6*kmesh_input%higher_order_n))
    real(kind=dp), intent(inout) :: bweight(kmesh_input%max_shells_h)

    ! local variables
    real(kind=dp) :: bweight_temp, fact
    integer :: shell, order, loop_j, temp_multi(kmesh_input%num_shells), temp_nnshell(num_kpts, kmesh_input%num_shells)
    real(kind=dp) :: temp_dnn(kmesh_input%num_shells)

    ! update new shells (after simplify the first-order shell list, e.g. 1,4,6, ... -> 1,2,3,...)
    do shell = 1, kmesh_input%num_shells
      temp_multi(shell) = multi(kmesh_input%shell_list(shell))
      temp_nnshell(:, shell) = nnshell(:, kmesh_input%shell_list(shell))
      temp_dnn(shell) = dnn(kmesh_input%shell_list(shell))
    enddo

    do shell = 1, kmesh_input%num_shells
      kmesh_input%shell_list(shell) = shell
      multi(shell) = temp_multi(shell)
      nnshell(:, shell) = temp_nnshell(:, shell)
      dnn(shell) = temp_dnn(shell)
    enddo

    do order = 2, kmesh_input%higher_order_n
      do shell = 1, kmesh_input%num_shells
        kmesh_input%shell_list((order - 1)*kmesh_input%num_shells + shell) = &
          (order - 1)*kmesh_input%num_shells + shell
        multi((order - 1)*kmesh_input%num_shells + shell) = multi(shell)
        nnshell(:, (order - 1)*kmesh_input%num_shells + shell) = nnshell(:, shell)
        dnn((order - 1)*kmesh_input%num_shells + shell) = dnn(shell)*order
      enddo
    enddo

    ! calculate new bweights w_b, w_2b, ..., w_Nb
    do shell = 1, kmesh_input%num_shells
      bweight_temp = bweight(shell)
      do order = 1, kmesh_input%higher_order_n
        fact = 1.0_dp/REAL(order**2, DP)
        do loop_j = 1, kmesh_input%higher_order_n
          if (loop_j == order) cycle
          fact = (fact*REAL(loop_j**2, DP))/REAL(loop_j**2 - order**2, DP)
        enddo
        bweight(kmesh_input%num_shells*(order - 1) + shell) = bweight_temp*fact
      enddo
    enddo

    kmesh_input%num_shells = kmesh_input%num_shells*kmesh_input%higher_order_n

    return

  end subroutine kmesh_shell_reconstruct

  !================================================
  subroutine kmesh_get_amat(kmesh_input, amat, bvector, multi, loop_order, num_x, num_y, num_z)
    !================================================
    !
    !!  Find amat(coefficients to find bweight) and the numbers of x, y, z components for a given order
    !
    !================================================

    use w90_types, only: kmesh_input_type

    implicit none

    ! arguments
    type(kmesh_input_type), intent(inout) :: kmesh_input
    real(kind=dp), intent(inout) :: amat(:, :)
    integer, intent(inout) :: num_x(:)
    integer, intent(inout) :: num_y(:)
    integer, intent(inout) :: num_z(:)
    integer, intent(in) :: multi(kmesh_input%search_shells)   ! the number of kpoints in the shell
    integer, intent(in) :: loop_order
    real(kind=dp), intent(in) :: bvector(3, maxval(multi), kmesh_input%max_shells_h)

    ! local variables
    integer :: num_of_eqs, num_of_eqs_prev
    integer :: loop_i, loop_j, loop_s, loop_b

    num_of_eqs = (1 + loop_order)*(1 + 2*loop_order) !((3, 2n)) (combi. with repetition)
    num_of_eqs_prev = loop_order*(2*loop_order - 1)
    if (loop_order .eq. 1) num_of_eqs_prev = 0

    do loop_s = 1, kmesh_input%num_shells
      do loop_i = 1, num_of_eqs
        ! equation index, e.g. If loop_order == 3, (1,2,...,15) -> (xxx, xxy, xyy, yyy, xxz, xyz, ..., zzz)
        ! find the number of z components corresponding to the current loop_i
        do loop_j = 0, 2*loop_order
          if ((2*loop_order + 1)*loop_j - loop_j*(loop_j - 1)/2 <= loop_i - 1 &
              .and. (2*loop_order + 1)*(loop_j + 1) - (loop_j + 1)*loop_j/2 > loop_i - 1) then
            num_z(loop_i) = loop_j
            exit
          endif
        enddo
        ! find the number of x and y components
        num_y(loop_i) = loop_i - 1 - ((2*loop_order + 1)*num_z(loop_i) - num_z(loop_i) &
                                      *(num_z(loop_i) - 1)/2)
        num_x(loop_i) = 2*loop_order - num_y(loop_i) - num_z(loop_i)
        ! calculate sum_b bb...bbbb
        do loop_b = 1, multi(kmesh_input%shell_list(loop_s))
          amat(num_of_eqs_prev + loop_i, loop_s) = amat(num_of_eqs_prev + loop_i, loop_s) &
                                                   + (bvector(1, loop_b, loop_s)**num_x(loop_i)) &
                                                   *(bvector(2, loop_b, loop_s)**num_y(loop_i)) &
                                                   *(bvector(3, loop_b, loop_s)**num_z(loop_i))
        enddo
      enddo
    enddo
  end subroutine kmesh_get_amat

  !================================================
  subroutine kmesh_check_condition(kmesh_input, bsat, bvector, bweight, multi, loop_order, num_x, num_y, num_z)
    !================================================
    !
    !!  Check if the obtained bweight satisfy the conditions for finite-difference, including (B1).
    !
    !================================================

    use w90_types, only: kmesh_input_type

    implicit none

    ! arguments
    type(kmesh_input_type), intent(inout) :: kmesh_input
    logical, intent(inout) :: bsat
    integer, intent(inout) :: num_x(:)
    integer, intent(inout) :: num_y(:)
    integer, intent(inout) :: num_z(:)
    integer, intent(in) :: multi(kmesh_input%search_shells)   ! the number of kpoints in the shell
    integer, intent(in) :: loop_order
    real(kind=dp), intent(in) :: bvector(3, maxval(multi), kmesh_input%max_shells_h)
    real(kind=dp), intent(in) :: bweight(kmesh_input%max_shells_h)

    ! local variables
    integer :: num_of_eqs
    integer :: loop_i, loop_s, loop_b
    real(kind=dp) :: delta

    num_of_eqs = (1 + loop_order)*(1 + 2*loop_order)
    do loop_i = 1, num_of_eqs
      delta = 0.0_dp
      do loop_s = 1, kmesh_input%num_shells
        do loop_b = 1, multi(kmesh_input%shell_list(loop_s))
          delta = delta + bweight(loop_s)*(bvector(1, loop_b, loop_s)**num_x(loop_i)) &
                  *(bvector(2, loop_b, loop_s)**num_y(loop_i)) &
                  *(bvector(3, loop_b, loop_s)**num_z(loop_i))
        end do
      end do
      if (loop_order .eq. 1 .and. (loop_i .eq. 1 .or. loop_i .eq. 3 .or. loop_i .eq. 6)) then
        if (abs(delta - 1.0_dp) > kmesh_input%tol) bsat = .false.
      else
        if (abs(delta) > kmesh_input%tol) bsat = .false.
      endif
    enddo

  end subroutine kmesh_check_condition

  !================================================
  subroutine kmesh_shell_fixed(kmesh_input, print_output, bweight, dnn, kpt_cart, recip_lattice, &
                               lmn, multi, num_kpts, stdout, timer, error, comm)
    !================================================
    !
    !!  Find the B1 weights for a set of shells specified by the user
    !
    !================================================

    use w90_constants, only: eps7
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_input_type, print_output_type, timer_list_type

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(kmesh_input_type), intent(in) :: kmesh_input
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_kpts
    integer, intent(in) :: stdout
    integer, intent(in) :: lmn(:, :)
    integer, intent(in) :: multi(kmesh_input%search_shells)   ! the number of kpoints in the shell

    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) ::kpt_cart(:, :)
    real(kind=dp), intent(in) :: dnn(kmesh_input%search_shells) ! the bvectors
    real(kind=dp), intent(out) :: bweight(max_shells)

    ! local variables
    real(kind=dp), allocatable     :: bvector(:, :, :)

    real(kind=dp) :: singv(kmesh_input%num_shells)
    real(kind=dp) :: amat(max_shells, kmesh_input%num_shells)
    real(kind=dp) :: umat(max_shells, max_shells)
    real(kind=dp) :: vmat(kmesh_input%num_shells, kmesh_input%num_shells)
    real(kind=dp) :: smat(kmesh_input%num_shells, max_shells)
    integer, parameter :: lwork = max_shells*10
    real(kind=dp) :: work(lwork)
    real(kind=dp), parameter :: target(6) = (/1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
    logical :: b1sat
    integer :: ierr, loop_i, loop_j, loop_b, loop_s, info
    real(kind=dp) :: delta

    integer :: loop, shell

    if (print_output%timing_level > 1) call io_stopwatch_start('kmesh: shell_fixed', timer)

    allocate (bvector(3, maxval(multi), kmesh_input%num_shells), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating bvector in kmesh_shell_fixed', comm)
      return
    endif
    bvector = 0.0_dp; bweight = 0.0_dp
    amat = 0.0_dp; umat = 0.0_dp; vmat = 0.0_dp; smat = 0.0_dp; singv = 0.0_dp

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a)') '| The b-vectors are set in the win file                                      |'
    endif

    do shell = 1, kmesh_input%num_shells
      ! get the b vectors for this shell
      call kmesh_get_bvectors(kmesh_input, print_output, &
                              bvector(:, 1:multi(kmesh_input%shell_list(shell)), shell), kpt_cart, &
                              recip_lattice, dnn(kmesh_input%shell_list(shell)), lmn, 1, &
                              multi(kmesh_input%shell_list(shell)), num_kpts, timer, error, comm)
      if (allocated(error)) return
    end do

    if (print_output%iprint >= 3) then
      do shell = 1, kmesh_input%num_shells
        write (stdout, '(1x,a8,1x,I2,a14,1x,I2,49x,a)') '| Shell:', shell, ' Multiplicity:', &
          multi(kmesh_input%shell_list(shell)), '|'
        do loop = 1, multi(kmesh_input%shell_list(shell))
          write (stdout, '(1x,a10,I2,1x,a1,4x,3f12.6,5x,a9,9x,a)') '| b-vector ', loop, ':', &
            bvector(:, loop, shell)/print_output%lenconfac, '('//trim(print_output%length_unit)//'^-1)', '|'
        end do
      end do
    end if

    do loop_s = 1, kmesh_input%num_shells
      do loop_b = 1, multi(kmesh_input%shell_list(loop_s))
        amat(1, loop_s) = amat(1, loop_s) + bvector(1, loop_b, loop_s)*bvector(1, loop_b, loop_s)
        amat(2, loop_s) = amat(2, loop_s) + bvector(2, loop_b, loop_s)*bvector(2, loop_b, loop_s)
        amat(3, loop_s) = amat(3, loop_s) + bvector(3, loop_b, loop_s)*bvector(3, loop_b, loop_s)
        amat(4, loop_s) = amat(4, loop_s) + bvector(1, loop_b, loop_s)*bvector(2, loop_b, loop_s)
        amat(5, loop_s) = amat(5, loop_s) + bvector(2, loop_b, loop_s)*bvector(3, loop_b, loop_s)
        amat(6, loop_s) = amat(6, loop_s) + bvector(3, loop_b, loop_s)*bvector(1, loop_b, loop_s)
      end do
    end do

    info = 0
    call dgesvd('A', 'A', max_shells, kmesh_input%num_shells, amat, max_shells, singv, umat, &
                max_shells, vmat, kmesh_input%num_shells, work, lwork, info)
    if (info < 0) then
      if (print_output%iprint > 0) then
        write (stdout, '(1x,a,1x,I1,1x,a)') 'kmesh_shell_fixed: Argument', abs(info), &
          'of dgesvd is incorrect'
      endif
      call set_error_fatal(error, 'kmesh_shell_fixed: Problem with Singular Value Decomposition', comm)
      return
    else if (info > 0) then
      call set_error_fatal(error, 'kmesh_shell_fixed: Singular Value Decomposition did not converge', comm)
      return
    end if

    if (any(abs(singv) < eps7)) then
      call set_error_fatal(error, 'kmesh_shell_fixed: Singular Value Decomposition has found a very small singular value', comm)
      return
    endif

    smat = 0.0_dp
    do loop_s = 1, kmesh_input%num_shells
      smat(loop_s, loop_s) = 1/singv(loop_s)
    end do

    bweight(1:kmesh_input%num_shells) = matmul(transpose(vmat), matmul(smat, matmul(transpose(umat), target)))
    if (print_output%iprint >= 2) then
      do loop_s = 1, kmesh_input%num_shells
        write (stdout, '(1x,a,I2,a,f12.7,5x,a8,36x,a)') '| Shell: ', loop_s, &
          ' w_b ', bweight(loop_s)*print_output%lenconfac**2, '('//trim(print_output%length_unit)//'^2)', '|'
      end do
    end if

    !check b1

    b1sat = .true.
    if (.not. kmesh_input%skip_B1_tests) then
      do loop_i = 1, 3
        do loop_j = loop_i, 3
          delta = 0.0_dp
          do loop_s = 1, kmesh_input%num_shells
            do loop_b = 1, multi(kmesh_input%shell_list(loop_s))
              delta = delta + bweight(loop_s)*bvector(loop_i, loop_b, loop_s)*bvector(loop_j, loop_b, loop_s)
            end do
          end do
          if (loop_i == loop_j) then
            if (abs(delta - 1.0_dp) > kmesh_input%tol) b1sat = .false.
          end if
          if (loop_i /= loop_j) then
            if (abs(delta) > kmesh_input%tol) b1sat = .false.
          end if
        end do
      end do
    end if

    if (.not. b1sat) then
      call set_error_fatal(error, 'kmesh_shell_fixed: B1 condition not satisfied', comm)
      return
    endif

    if (print_output%timing_level > 1) call io_stopwatch_stop('kmesh: shell_fixed', timer)

    return

  end subroutine kmesh_shell_fixed

  !================================================
  subroutine kmesh_shell_from_file(kmesh_input, print_output, bvec_inp, bweight, dnn, kpt_cart, &
                                   recip_lattice, lmn, multi, num_kpts, seedname, stdout, timer, &
                                   error, comm)
    !================================================
    !!  Find the B1 weights for a set of b-vectors given in a file.
    !!  This routine is only activated via a devel_flag and is not
    !!  intended for regular use.
    !
    !================================================

    use w90_constants, only: eps7, maxlen
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_input_type, print_output_type, timer_list_type

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(kmesh_input_type), intent(inout) :: kmesh_input
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_kpts, stdout
    integer, intent(in) :: lmn(:, :)
    integer, intent(inout) :: multi(kmesh_input%search_shells)   ! the number of kpoints in the shell

    real(kind=dp), intent(in) :: recip_lattice(3, 3)
    real(kind=dp), intent(in) ::kpt_cart(:, :)
    real(kind=dp), intent(inout) :: bvec_inp(:, :, :)
    real(kind=dp), intent(in) :: dnn(kmesh_input%search_shells)  ! the bvectors
    real(kind=dp), intent(out) :: bweight(max_shells)

    character(len=50), intent(in)  :: seedname

    ! local variables
    real(kind=dp), allocatable     :: bvector(:, :)

    real(kind=dp), dimension(:), allocatable :: singv
    real(kind=dp), dimension(:, :), allocatable :: amat, umat, vmat, smat

    integer, parameter :: lwork = max_shells*10
    real(kind=dp) :: work(lwork)
    integer       :: bvec_list(num_nnmax, max_shells)
    real(kind=dp), parameter :: target(6) = (/1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
    logical :: b1sat
    integer :: ierr, loop_i, loop_j, loop_b, loop_s, info
    real(kind=dp) :: delta

    integer :: loop, shell, pos, kshell_in, counter, length, i, loop2, num_lines, tot_num_lines
    character(len=maxlen) :: dummy, dummy2

    if (print_output%timing_level > 1) call io_stopwatch_start('kmesh: shell_fixed', timer)

    allocate (bvector(3, sum(multi)), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating bvector in kmesh_shell_fixed', comm)
      return
    endif
    bvector = 0.0_dp; bweight = 0.0_dp

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a)') '| The b-vectors are defined in the kshell file                               |'
    endif

    counter = 1
    do shell = 1, kmesh_input%search_shells
      ! get the b vectors
      call kmesh_get_bvectors(kmesh_input, print_output, bvector(:, counter:counter + multi(shell) - 1), &
                              kpt_cart, recip_lattice, dnn(shell), lmn, 1, multi(shell), num_kpts, &
                              timer, error, comm)
      if (allocated(error)) return

      counter = counter + multi(shell)
    end do

    open (newunit=kshell_in, file=trim(seedname)//'.kshell', form='formatted', status='old', &
          action='read', iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error: Problem (1) opening input file '//trim(seedname)//'.kshell', comm)
      return
    endif

    num_lines = 0; tot_num_lines = 0
    do
      read (kshell_in, '(a)', iostat=ierr) dummy
      if (ierr == 0) then !read ok, proceed
        dummy = adjustl(dummy)
        tot_num_lines = tot_num_lines + 1
        if (.not. dummy(1:1) == '!' .and. .not. dummy(1:1) == '#') then
          if (len(trim(dummy)) > 0) num_lines = num_lines + 1
        endif
      else if (ierr > 0) then !error case
        call set_error_input(error, 'Error: Problem (2) reading input file '//trim(seedname)//'.kshell', comm)
        return
      else if (ierr < 0) then !end of record or end of file
        exit
      end if
    end do

    rewind (kshell_in)
    kmesh_input%num_shells = num_lines

    multi(:) = 0
    bvec_list = 1
    counter = 0
    do loop = 1, tot_num_lines
      read (kshell_in, '(a)', err=103, end=103) dummy2
      dummy2 = adjustl(dummy2)
      if (dummy2(1:1) == '!' .or. dummy2(1:1) == '#' .or. (len(trim(dummy2)) == 0)) cycle
      counter = counter + 1
      kmesh_input%shell_list(counter) = counter
      dummy = dummy2
      length = 1
      dummy = adjustl(dummy)
      do
        pos = index(dummy, ' ')
        dummy = dummy(pos + 1:)
        dummy = adjustl(dummy)
        if (len_trim(dummy) > 0) then
          length = length + 1
        else
          exit
        endif

      end do
      multi(counter) = length
      read (dummy2, *, err=230, end=230) (bvec_list(i, loop), i=1, length)
    end do

    bvec_inp = 0.0_dp
    do loop = 1, kmesh_input%num_shells
      do loop2 = 1, multi(loop)
        bvec_inp(:, loop2, loop) = bvector(:, bvec_list(loop2, loop))
      end do
    end do

    if (print_output%iprint >= 3) then
      do shell = 1, kmesh_input%num_shells
        write (stdout, '(1x,a8,1x,I2,a14,1x,I2,49x,a)') '| Shell:', shell, ' Multiplicity:', multi(shell), '|'
        do loop = 1, multi(shell)
          write (stdout, '(1x,a10,I2,1x,a1,4x,3f12.6,5x,a9,9x,a)') '| b-vector ', loop, ':', &
            bvec_inp(:, loop, shell)/print_output%lenconfac, '('//trim(print_output%length_unit)//'^-1)', '|'
        end do
      end do
    end if

    allocate (amat(max_shells, kmesh_input%num_shells), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating amat in kmesh_shell_from_file', comm)
      return
    endif
    allocate (umat(max_shells, max_shells), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating umat in kmesh_shell_from_file', comm)
      return
    endif
    allocate (vmat(kmesh_input%num_shells, kmesh_input%num_shells), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating vmat in kmesh_shell_from_file', comm)
      return
    endif
    allocate (smat(kmesh_input%num_shells, max_shells), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating smat in kmesh_shell_from_file', comm)
      return
    endif
    allocate (singv(kmesh_input%num_shells), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating singv in kmesh_shell_from_file', comm)
      return
    endif
    amat = 0.0_dp; umat = 0.0_dp; vmat = 0.0_dp; smat = 0.0_dp; singv = 0.0_dp

    do loop_s = 1, kmesh_input%num_shells
      do loop_b = 1, multi(loop_s)
        amat(1, loop_s) = amat(1, loop_s) + bvec_inp(1, loop_b, loop_s)*bvec_inp(1, loop_b, loop_s)
        amat(2, loop_s) = amat(2, loop_s) + bvec_inp(2, loop_b, loop_s)*bvec_inp(2, loop_b, loop_s)
        amat(3, loop_s) = amat(3, loop_s) + bvec_inp(3, loop_b, loop_s)*bvec_inp(3, loop_b, loop_s)
        amat(4, loop_s) = amat(4, loop_s) + bvec_inp(1, loop_b, loop_s)*bvec_inp(2, loop_b, loop_s)
        amat(5, loop_s) = amat(5, loop_s) + bvec_inp(2, loop_b, loop_s)*bvec_inp(3, loop_b, loop_s)
        amat(6, loop_s) = amat(6, loop_s) + bvec_inp(3, loop_b, loop_s)*bvec_inp(1, loop_b, loop_s)
      end do
    end do

    info = 0
    call dgesvd('A', 'A', max_shells, kmesh_input%num_shells, amat, max_shells, singv, umat, &
                max_shells, vmat, kmesh_input%num_shells, work, lwork, info)
    if (info < 0) then
      if (print_output%iprint > 0) then
        write (stdout, '(1x,a,1x,I1,1x,a)') 'kmesh_shell_fixed: Argument', abs(info), &
          'of dgesvd is incorrect'
      endif
      call set_error_fatal(error, 'kmesh_shell_fixed: Problem with Singular Value Decomposition', comm)
      return
    else if (info > 0) then
      call set_error_fatal(error, 'kmesh_shell_fixed: Singular Value Decomposition did not converge', comm)
      return
    end if

    if (any(abs(singv) < eps7)) then
      call set_error_fatal(error, 'kmesh_shell_fixed: Singular Value Decomposition has found a very small singular value', comm)
      return
    endif

    smat = 0.0_dp
    do loop_s = 1, kmesh_input%num_shells
      smat(loop_s, loop_s) = 1/singv(loop_s)
    end do

    bweight(1:kmesh_input%num_shells) = matmul(transpose(vmat), matmul(smat, matmul(transpose(umat), target)))
    if (print_output%iprint >= 2) then
      do loop_s = 1, kmesh_input%num_shells
        write (stdout, '(1x,a,I2,a,f12.7,5x,a8,36x,a)') '| Shell: ', loop_s, &
          ' w_b ', bweight(loop_s)*print_output%lenconfac**2, '('//trim(print_output%length_unit)//'^2)', '|'
      end do
    end if

    !check b1
    b1sat = .true.
    if (.not. kmesh_input%skip_B1_tests) then
      do loop_i = 1, 3
        do loop_j = loop_i, 3
          delta = 0.0_dp
          do loop_s = 1, kmesh_input%num_shells
            do loop_b = 1, multi(loop_s)
              delta = delta + bweight(loop_s)*bvec_inp(loop_i, loop_b, loop_s)*bvec_inp(loop_j, loop_b, loop_s)
            end do
          end do
          if (loop_i == loop_j) then
            if (abs(delta - 1.0_dp) > kmesh_input%tol) b1sat = .false.
          end if
          if (loop_i /= loop_j) then
            if (abs(delta) > kmesh_input%tol) b1sat = .false.
          end if
        end do
      end do
    end if

    if (.not. b1sat) then
      call set_error_fatal(error, 'kmesh_shell_fixed: B1 condition not satisfied', comm)
      return
    endif

    if (print_output%timing_level > 1) call io_stopwatch_stop('kmesh: shell_fixed', timer)

    return

103 call set_error_input(error, 'Error: Problem (3) reading input file '//trim(seedname)//'.kshell', comm)
    return
230 call set_error_input(error, 'Error: Problem reading in w90_readwrite_get_keyword_vector', comm)
    return

  end subroutine kmesh_shell_from_file

  !================================================
  function internal_maxloc(dist, nsupcell)
    !================================================
    !!  A reproducible maxloc function
    !!  so b-vectors come in the same
    !!  order each time
    !================================================

    use w90_constants, only: eps8

    implicit none
    integer, intent(in) :: nsupcell
    real(kind=dp), intent(in)  :: dist((2*nsupcell + 1)**3)
    !! Distances from the origin of the unit cells in the supercell.
    integer :: internal_maxloc

    integer :: guess(1), loop, counter
    integer :: list((2*nsupcell + 1)**3)

    list = 0
    counter = 1

    guess = maxloc(dist)
    list(1) = guess(1)
    ! look for any degenerate values
    do loop = 1, (2*nsupcell + 1)**3
      if (loop == guess(1)) cycle
      if (abs(dist(loop) - dist(guess(1))) < eps8) then
        counter = counter + 1
        list(counter) = loop
      endif
    end do
    ! and always return the lowest index
    internal_maxloc = minval(list(1:counter))

  end function internal_maxloc

  !================================================
  subroutine kmesh_bvectors_perm(bvec, bref, num_kpt, num_bvec, perm, invperm, revind, error, comm)
    !================================================
    !
    !!  Obtain possible permutation in ordering of b-vectors at different kpoints
    !
    !================================================
    implicit none

    ! arguments
    real(kind=dp), intent(in) :: bvec(:, :, :) ! set of bvecs for each k, possibly permuted, size (3,num_bvec,num_kpt)
    real(kind=dp), intent(in) :: bref(:, :) ! reference vector ordering, size (3,num_bvec)
    integer, intent(in) :: num_kpt, num_bvec
    integer, intent(inout) :: perm(:, :) ! assumed allocated
    integer, intent(inout) :: invperm(:, :) ! assumed allocated
    integer, intent(inout) :: revind(:, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    real(kind=dp), parameter :: tol = 1d-7 ! this should not be smaller than the k-point precision in the .win file
    integer :: ik, n, m
    logical :: found, found2

    do ik = 1, num_kpt
      do n = 1, num_bvec
        found = .false.
        found2 = .false.
        do m = 1, num_bvec
          if (all(abs(bvec(:, m, ik) - bref(:, n)) < tol)) then
            found = .true.
            perm(n, ik) = m
            invperm(m, ik) = n ! inverse mapping, used in postw90
            cycle
          endif
          if (all(abs(bvec(:, m, ik) + bref(:, n)) < tol)) then
            found2 = .true.
            revind(n, ik) = m
            cycle
          endif
        enddo
        if (.not. found .or. .not. found2) then
      call set_error_fatal(error, 'Unable to identify bk-vector permutation (kmesh_bvector_perm); consider k-point precision', comm)
          return
        endif
      enddo
    enddo

  end subroutine kmesh_bvectors_perm

end module w90_kmesh
