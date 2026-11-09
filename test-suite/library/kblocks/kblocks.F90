! Exercise comms_allgather_kblocks on a spread of k-point distributions.
!
! Each k-point is owned by exactly one rank, so assembling the global array out of the
! per-rank blocks is a copy. That is why this test demands bit-exact equality rather than a
! tolerance, and why it compares on EVERY rank rather than only on the root: the routine is
! expected to leave every rank holding the same assembled array.
!
! The global array is pre-filled with a sentinel, so a block the routine fails to write is
! a failure too, not a silently plausible value.
!
! Which of the routine's internal paths a given distribution takes is deliberately not
! asserted here -- it is not observable from outside and it is not part of the contract.
! What is asserted is that every distribution assembles the same correct array, and that a
! malformed one raises an error instead of corrupting it.
!
! Run with no arguments to sweep every case at several mesh sizes, which is what ctest
! does. Arguments <case> <nkpts> run a single combination, for debugging.
program kblocks

  use mpi
  use w90_constants, only: dp
  use w90_comms, only: comms_allgather_kblocks, w90_comm_type
  use w90_error_base, only: w90_error_type

  implicit none

  integer, parameter :: ncase = 9      ! 1-7 must succeed, 8-9 must be rejected
  integer, parameter :: nd = 4         ! leading dimensions of each block

  type(w90_comm_type) :: comm

  character(len=24) :: arg, case_name
  integer :: icase, ierr, imesh, jcase, nkpts, my_rank, num_nodes
  integer :: nfail, nfail_all, nmesh, nrun
  integer :: nk_list(3)
  integer, allocatable :: owner(:)
  logical :: single

  call mpi_init(ierr)
  comm%comm = MPI_COMM_WORLD
  call mpi_comm_rank(comm%comm, my_rank, ierr)
  call mpi_comm_size(comm%comm, num_nodes, ierr)

  single = (command_argument_count() >= 2)
  if (single) then
    call get_command_argument(1, arg); read (arg, *) icase
    call get_command_argument(2, arg); read (arg, *) nkpts
    nk_list = nkpts
    nmesh = 1
  else
    ! One k-point per rank is the degenerate case worth covering: several distributions
    ! collapse onto the identity there. The third size does not divide the rank count, so
    ! the blocks are unequal.
    nk_list = [num_nodes, 3*num_nodes, 3*num_nodes + 1]
    nmesh = 3
  end if

  nfail = 0
  nrun = 0
  if (my_rank == 0) then
    write (*, '(a,i0,a)') 'comms_allgather_kblocks: ', num_nodes, ' rank(s)'
    write (*, '(2x,a4,2x,a16,2x,a6,2x,a6,2x,a7,2x,a)') &
      'case', 'distribution', 'nkpts', 'bad', 'errors', 'verdict'
  end if

  if (single) then
    call one(icase, nkpts)
  else
    do imesh = 1, nmesh
      do jcase = 1, ncase
        call one(jcase, nk_list(imesh))
      end do
    end do
  end if

  call mpi_allreduce(nfail, nfail_all, 1, MPI_INTEGER, MPI_SUM, comm%comm, ierr)
  if (my_rank == 0) then
    if (nfail_all == 0) then
      write (*, '(a,i0,a)') 'PASS: ', nrun, ' combination(s), every one exact'
    else
      write (*, '(a,i0,a,i0)') 'FAIL: ', nfail_all, ' of ', nrun
    end if
  end if

  call mpi_finalize(ierr)
  ! Every rank holds the same reduced verdict, so they all stop the same way and ctest
  ! sees a non-zero status rather than needing the output parsed.
  if (nfail_all /= 0) stop 1

contains

  !> Build one distribution, call the routine, and check the assembled array.
  subroutine one(kase, nk)
    integer, intent(in) :: kase, nk

    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: array(:, :, :), array_loc(:, :, :)
    integer, allocatable :: global_k(:)
    integer :: i, ik, ikl, j, nbad, nbad_all, nd_loc, nerr, nerr_all, nloc, nloc_decl
    logical :: want_error, good
    character(len=7) :: verdict

    if (allocated(owner)) deallocate (owner)
    allocate (owner(nk))
    call distribute(kase, nk, case_name)

    nloc = count(owner == my_rank)
    allocate (global_k(nloc))
    ikl = 0
    do ik = 1, nk
      if (owner(ik) == my_rank) then
        ikl = ikl + 1
        global_k(ikl) = ik
      end if
    end do

    ! Cases 8 and 9 malform the call on purpose: a leading dimension that disagrees with
    ! the global array, and more local blocks than global_k accounts for. Both must be
    ! rejected by the routine's own checks.
    nd_loc = nd
    nloc_decl = nloc
    if (kase == 8) nd_loc = nd + 1
    if (kase == 9) nloc_decl = nloc + 1
    want_error = (kase == 8 .or. kase == 9)

    ! A rank that owns no k-point declares a zero-sized block: that is what the shape
    ! the routine is handed has to say, and padding it to one would be a shape mismatch.
    allocate (array_loc(nd_loc, nd, nloc_decl))
    array_loc = cmplx(0.0_dp, 0.0_dp, dp)
    do ikl = 1, nloc
      do j = 1, nd
        do i = 1, nd
          array_loc(i, j, ikl) = value_at(i, j, global_k(ikl))
        end do
      end do
    end do

    allocate (array(nd, nd, nk))
    array = cmplx(-999.0_dp, 999.0_dp, dp)     ! sentinel

    call comms_allgather_kblocks(array_loc, array, global_k, error, comm)

    nbad = 0
    nerr = 0
    if (allocated(error)) nerr = 1
    if (.not. allocated(error)) then
      do ik = 1, nk
        do j = 1, nd
          do i = 1, nd
            if (array(i, j, ik) /= reference(nk, i, j, ik)) nbad = nbad + 1
          end do
        end do
      end do
    end if
    call mpi_allreduce(nbad, nbad_all, 1, MPI_INTEGER, MPI_SUM, comm%comm, ierr)
    ! The error is raised per rank, and only rank 0 prints and counts. Reducing it keeps a
    ! failure on some other rank from going unnoticed -- and a routine that raised it on
    ! only some of the ranks would leave the rest waiting in the collective.
    call mpi_allreduce(nerr, nerr_all, 1, MPI_INTEGER, MPI_SUM, comm%comm, ierr)

    if (want_error) then
      good = (nerr_all == num_nodes)
    else
      good = (nerr_all == 0) .and. (nbad_all == 0)
    end if
    verdict = 'ok'
    if (.not. good) verdict = 'FAILED'
    if (.not. good .and. my_rank == 0) nfail = nfail + 1
    nrun = nrun + 1

    if (my_rank == 0) then
      write (*, '(2x,i4,2x,a16,2x,i6,2x,i6,2x,i7,2x,a)') &
        kase, case_name, nk, nbad_all, nerr_all, trim(verdict)
    end if

    deallocate (array, array_loc, global_k)
    if (allocated(error)) deallocate (error)
  end subroutine one

  !> The value k-point ik is supposed to carry: distinct for every (i, j, ik) and exactly
  !> representable, so equality is a fair demand.
  complex(kind=dp) function value_at(i1, i2, ik1)
    integer, intent(in) :: i1, i2, ik1
    value_at = cmplx(real(1000*ik1 + 10*i1 + i2, dp), real(-(3*ik1 + i1 - i2), dp), dp)
  end function value_at

  !> What the assembled array must hold. A k-point owned by nobody carries no data, so the
  !> fallback path has to leave a zero there rather than a sentinel.
  complex(kind=dp) function reference(nk, i1, i2, ik1)
    integer, intent(in) :: nk, i1, i2, ik1
    if (owner(ik1) < 0) then
      reference = cmplx(0.0_dp, 0.0_dp, dp)
    else
      reference = value_at(i1, i2, ik1)
    end if
  end function reference

  !> Fill owner(1:nk) with the rank that owns each k-point, -1 for nobody.
  subroutine distribute(kase, nk, name)
    integer, intent(in) :: kase, nk
    character(len=*), intent(out) :: name
    integer :: base, cc, ic, ip, kk, no, rest
    integer :: order(num_nodes)

    owner = -1
    do ip = 1, num_nodes
      order(ip) = ip - 1
    end do

    select case (kase)
    case (1)                              ! what comms_array_split produces
      name = 'contiguous'
      call runs(order(1:num_nodes), nk)
    case (2)                              ! contiguous, unequal counts
      name = 'unequal'
      call runs(order(1:num_nodes), nk)
      ! nk may already divide the rank count, in which case runs() hands out equal
      ! blocks. Move the first k-point of rank 1's run to rank 0: the counts then differ
      ! while every rank's run stays consecutive and in rank order.
      if (num_nodes > 1 .and. nk > num_nodes) then
        cc = count(owner == 0)
        owner(cc + 1) = 0
      end if
    case (3)                              ! contiguous, the last rank owns nothing
      name = 'empty-rank'
      if (num_nodes > 1) then
        call runs(order(1:num_nodes - 1), nk)
      else
        call runs(order(1:1), nk)
      end if
    case (4)                              ! blocks handed out in reverse rank order
      name = 'reversed'
      do ip = 1, num_nodes
        order(ip) = num_nodes - ip
      end do
      call runs(order(1:num_nodes), nk)
    case (5)                              ! strided: nothing consecutive anywhere
      name = 'round-robin'
      do kk = 1, nk
        owner(kk) = mod(kk - 1, num_nodes)
      end do
    case (6)                              ! contiguous but for a single swap
      name = 'one-swap'
      call runs(order(1:num_nodes), nk)
      if (nk > 1) then
        ic = owner(1); owner(1) = owner(nk); owner(nk) = ic
      end if
    case (7)                              ! each rank consecutive, one k-point unowned
      name = 'gap'
      if (num_nodes > 1 .and. nk > num_nodes) then
        cc = nk/num_nodes
        owner(1:cc) = 0                   ! k = cc+1 is left to nobody
        rest = nk - cc - 1
        base = rest/(num_nodes - 1)
        ic = mod(rest, num_nodes - 1)
        kk = cc + 2
        do ip = 1, num_nodes - 1
          no = base
          if (ip <= ic) no = no + 1
          if (no > 0) owner(kk:kk + no - 1) = ip
          kk = kk + no
        end do
      else
        call runs(order(1:num_nodes), nk)
        owner(nk) = -1
      end if
    case (8, 9)                           ! the guards, not the collective
      name = 'malformed'
      call runs(order(1:num_nodes), nk)
    end select
  end subroutine distribute

  !> Hand out consecutive runs of k-points to the ranks listed in order, remainder first.
  subroutine runs(order, nk)
    integer, intent(in) :: order(:), nk
    integer :: base, cc, ip, kk, np, rest
    np = size(order)
    base = nk/np
    rest = mod(nk, np)
    kk = 1
    do ip = 1, np
      cc = base
      if (ip <= rest) cc = cc + 1
      if (cc > 0) owner(kk:kk + cc - 1) = order(ip)
      kk = kk + cc
    end do
  end subroutine runs

end program kblocks
