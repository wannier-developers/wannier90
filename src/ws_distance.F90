!-*- mode: F90 -*-!
!------------------------------------------------------------!
! Copyright (C) 2026 Wannier Developer Group                 !
!                                                            !
! This library is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU Lesser General !
! Public License as published by the Free Software           !
! Foundation; either version 2.1 of the License, or (at your !
! option) any later version.                                 !
!                                                            !
! This library is distributed in the hope that it will be    !
! useful,but WITHOUT ANY WARRANTY; without even the implied  !
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR    !
! PURPOSE.  See the GNU Lesser General Public License for    !
! more details.                                              !
!                                                            !
! You should have received a copy of the GNU Lesser General  !
! Public License along with this library; if not, see        !
! <https://www.gnu.org/licenses/>.                           !
!                                                            !
! The webpage of the Wannier90 code is                       !
! <https://www.wannier.org>.                                 !
!                                                            !
! The Wannier90 code is hosted on GitHub                     !
! <https://github.com/wannier-developers/wannier90>          !
!------------------------------------------------------------!
!                                                            !
! ws_distance:                                               !
! Original implementation by Lorenzo Paulatto, with later    !
! modifications by Marco Gibertini, Dominik Gresch           !
! and Giovanni Pizzi                                         !
!                                                            !
!------------------------------------------------------------!

module w90_ws_distance

  !! This module computes the optimal Wigner-Seitz cell around each Wannier
  !! function to use for interpolation.

  ! Short documentation follows, for a longer explanation see the documentation
  ! of the use_ws_distance variable in the user guide.
  !
  ! Some comments:
  ! 1. This computation is done independently on all processors (when run in
  !    parallel). I think this shouldn't do a problem as the math is fairly simple
  !    and uses data already broadcasted (integer values, and the
  !    wannier_centres), but if there is the risk of having different
  !    degeneracies or similar things on different MPI processors, we should
  !    probably think to do the math on node 0, and then broadcast results.

  use w90_constants, only: dp
  use w90_error

  implicit none

  private

  public :: clean_ws_translate
  public :: ws_apply_ndegen
  public :: ws_expand_rvec
  public :: ws_translate_dist
  public :: ws_write_vec

  integer, parameter :: ndegenx = 8
  !! max number of unit cells that can touch
  !! in a single point (i.e.  vertex of cube)

contains

  !================================================!

  subroutine ws_translate_dist(ws_distance, ws_region, num_wann, wannier_centres, real_lattice, &
                               mp_grid, nrpts, irvec, error, comm, force_recompute)
    !================================================!
    !! Find the supercell translation (i.e. the translation by a integer number of
    !! supercell vectors, the supercell being defined by the mp_grid) that
    !! minimizes the distance between two given Wannier functions, i and j,
    !! the first in unit cell 0, the other in unit cell R.
    !! I.e., we find the translation to put WF j in the Wigner-Seitz of WF i.
    !! We also look for the number of equivalent translation, that happen when w_j,R
    !! is on the edge of the WS of w_i,0. The results are stored in global
    !! arrays wdist_ndeg, irdist_ws, crdist_ws.
    !================================================!

    use w90_utility, only: utility_cart_to_frac, utility_frac_to_cart, utility_inverse_mat
    use w90_types, only: ws_region_type, ws_distance_type

    implicit none

    type(ws_distance_type), intent(inout) :: ws_distance
    type(ws_region_type), intent(in) :: ws_region
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann
    integer, intent(in) :: nrpts
    integer, intent(in) :: irvec(:, :)

    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: wannier_centres(:, :)

    logical, optional, intent(in):: force_recompute ! set to true to force recomputing everything

    ! local variables
    real(kind=dp) :: inv_lattice(3, 3)
    integer  :: iw, jw, ideg, ir, ierr
    integer :: shifts(3, ndegenx)
    real(DP) :: irvec_cart(3), tmp(3), tmp_frac(3), R_out(3, ndegenx)

    ! The subroutine does nothing if called more than once, which may
    ! not be the best thing if you invoke it while the WFs are moving
    if (present(force_recompute)) then
      if (force_recompute) then
        call clean_ws_translate(ws_distance, error, comm)
        if (allocated(error)) return
      end if
    end if
    if (ws_distance%done) return
    ws_distance%done = .true.

    if (ndegenx*num_wann*nrpts <= 0) then
      call set_error_fatal(error, "unexpected dimensions in ws_translate_dist", comm)
      return
    end if

    allocate (ws_distance%irdist(3, ndegenx, num_wann, num_wann, nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating irdist_ws in ws_translate_dist', comm)
      return
    end if
    allocate (ws_distance%crdist(3, ndegenx, num_wann, num_wann, nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating crdist_ws in ws_translate_dist', comm)
      return
    end if
    allocate (ws_distance%ndeg(num_wann, num_wann, nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating wcenter_ndeg in ws_translate_dist', comm)
      return
    end if

    !translation_centre_frac = 0._dp
    ws_distance%ndeg = 0
    ws_distance%irdist = 0
    ws_distance%crdist = 0

    call utility_inverse_mat(real_lattice, inv_lattice)
    do ir = 1, nrpts
      do jw = 1, num_wann
        do iw = 1, num_wann
          call utility_frac_to_cart(REAL(irvec(:, ir), kind=dp), irvec_cart, real_lattice)
          ! function JW translated in the Wigner-Seitz around function IW
          ! and also find its degeneracy, and the integer shifts needed
          ! to identify it
          ! Note: the routine outputs R_out, but we don't really need it
          ! This is kept in case in the future we might want to use it
          ! R_out contains the actual vector between the two WFs. We
          ! calculate instead crdist_ws, that is the Bravais lattice vector
          ! between two supercell lattices, that is the only one we need
          ! later for interpolation etc.
          call r_wz_sc(-wannier_centres(:, iw) &
                       + (irvec_cart + wannier_centres(:, jw)), (/0._dp, 0._dp, 0._dp/), &
                       ws_distance%ndeg(iw, jw, ir), R_out, shifts, mp_grid, real_lattice, &
                       inv_lattice, ws_region%ws_search_size, ws_region%ws_distance_tol, &
                       error, comm)
          if (allocated(error)) return

          do ideg = 1, ws_distance%ndeg(iw, jw, ir)
            ws_distance%irdist(:, ideg, iw, jw, ir) = irvec(:, ir) + shifts(:, ideg)
            tmp_frac = REAL(ws_distance%irdist(:, ideg, iw, jw, ir), kind=dp)
            CALL utility_frac_to_cart(tmp_frac, tmp, real_lattice)
            ws_distance%crdist(:, ideg, iw, jw, ir) = tmp
          end do
        end do
      end do
    end do
  end subroutine ws_translate_dist

  !================================================!
  subroutine R_wz_sc(R_in, R0, ndeg, R_out, shifts, mp_grid, real_lattice, inv_lattice, &
                     ws_search_size, ws_distance_tol, error, comm)
    !================================================!
    !! Put R_in in the Wigner-Seitz cell centered around R0,
    !! and find all equivalent vectors to this (i.e., with same distance).
    !! Return their coordinates and the degeneracy, as well as the integer
    !! shifts needed to get the vector (these are always multiples of
    !! the mp_grid, i.e. they are supercell displacements in the large supercell)
    !================================================!

    use w90_utility, only: utility_cart_to_frac, utility_frac_to_cart

    implicit none

    ! arguments
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: ws_search_size(3)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: inv_lattice(3, 3)
    real(kind=dp), intent(in) :: ws_distance_tol
    real(DP), intent(in) :: R_in(3)
    real(DP), intent(in) :: R0(3)
    integer, intent(out) :: ndeg
    real(DP), intent(out) :: R_out(3, ndegenx)
    integer, intent(out) :: shifts(3, ndegenx)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    real(DP) :: R(3), R_f(3), R_in_f(3), R_bz(3), mod2_R_bz
    integer :: i, j, k

    ! init
    ndeg = 0
    R_out = 0._dp
    shifts = 0
    R_bz = R_in
    mod2_R_bz = SUM((R_bz - R0)**2)
    !
    ! take R_bz to cryst(frac) coord for translating
    call utility_cart_to_frac(R_bz, R_in_f, inv_lattice)

    ! In this first loop, I just look for the shortest vector that I obtain
    ! by trying to displace the second Wannier function by all
    ! 'large-supercell' vectors
    ! The size of the supercell, controlled by ws_search_size,
    ! is incremented by one unit in order to account for WFs whose centre
    ! wanders away from the original reference unit cell
    do i = -ws_search_size(1) - 1, ws_search_size(1) + 1
      do j = -ws_search_size(2) - 1, ws_search_size(2) + 1
        do k = -ws_search_size(3) - 1, ws_search_size(3) + 1

          R_f = R_in_f + REAL((/i*mp_grid(1), j*mp_grid(2), k*mp_grid(3)/), &
                              kind=DP)
          call utility_frac_to_cart(R_f, R, real_lattice)

          if (SUM((R - R0)**2) < mod2_R_bz) then
            R_bz = R
            mod2_R_bz = SUM((R_bz - R0)**2)
            ! I start to set a first shift that is applied to get R_bz.
            ! Note: I reset these every time I find a smaller vector.
            !
            ! At this stage, this is the same for all potentially degenerate
            ! points (hence the use of : in shifts(1,:), for instance)
            ! In the second loop below, this shift will be added to the
            ! additional shift that differs for each degenerate but
            ! equivalent point
            shifts(1, :) = i*mp_grid(1)
            shifts(2, :) = j*mp_grid(2)
            shifts(3, :) = k*mp_grid(3)
          end if
        end do
      end do
    end do

    ! Now, second loop to find the list of R_out that differ from R_in
    ! by a large-supercell lattice vector and are equally distant from R0
    ! (i.e. that are on the edges of the WS cell centered on R0)
    ! As above, the size of the supercell, controlled by ws_search_size,
    ! is incremented by one unit in order to account for WFs whose centre
    ! wanders away from the original reference unit cell

    ! I start from the last R_bz found
    mod2_R_bz = SUM((R_bz - R0)**2)
    ! check if R0 and R_in are the same vector
    if (mod2_R_bz < ws_distance_tol**2) then
      ndeg = 1
      R_out(:, 1) = R0
      ! I can safely return as 'shifts' is already set
      return
    end if
    !
    ! take R_bz to cryst(frac) coord for translating
    call utility_cart_to_frac(R_bz, R_in_f, inv_lattice)

    do i = -ws_search_size(1) - 1, ws_search_size(1) + 1
      do j = -ws_search_size(2) - 1, ws_search_size(2) + 1
        do k = -ws_search_size(3) - 1, ws_search_size(3) + 1

          r_f = r_in_f + real((/i*mp_grid(1), j*mp_grid(2), k*mp_grid(3)/), &
                              kind=DP)
          call utility_frac_to_cart(R_f, R, real_lattice)

          if (abs(sqrt(sum((r - r0)**2)) - sqrt(mod2_r_bz)) < ws_distance_tol) then
            ndeg = ndeg + 1
            if (ndeg > ndegenx) then
              call set_error_fatal(error, "surprising ndeg, I wouldn't expect a degeneracy larger than 8...", comm)
              return
            end if
            R_out(:, ndeg) = R
            ! I return/update also the shifts. Note that I have to sum these
            ! to the previous value since in this second loop I am using
            ! R_bz (from the first loop) as the 'central' reference point,
            ! that is already shifted by shift(:,ndeg)
            shifts(1, ndeg) = shifts(1, ndeg) + i*mp_grid(1)
            shifts(2, ndeg) = shifts(2, ndeg) + j*mp_grid(2)
            shifts(3, ndeg) = shifts(3, ndeg) + k*mp_grid(3)
          end if

        end do
      end do
    end do
    !================================================!
  end subroutine R_wz_sc
  !================================================!

  !================================================!
  subroutine ws_write_vec(ws_distance, nrpts, irvec, num_wann, use_ws_distance, &
                          write_ndegen_applied, seedname, error, comm)
    !================================================!
    !! Write to file the lattice vectors of the superlattice
    !! to be added to R vector in seedname_hr.dat, seedname_rmn.dat, etc.
    !! in order to have the second Wannier function inside the WS cell
    !! of the first one.
    !!
    !! With write_ndegen_applied those shifts are already folded into the
    !! real-space output files, whose R list no longer matches the one written
    !! here. The file is then informational only and the header says so.
    !================================================!

    use w90_io, only: io_date
    use w90_types, only: ws_distance_type

    implicit none

    type(ws_distance_type), intent(in) :: ws_distance
    type(w90_error_type), allocatable, intent(out) :: error
    integer, intent(in) :: num_wann
    logical, intent(in) :: use_ws_distance
    logical, intent(in) :: write_ndegen_applied
    character(len=50), intent(in)  :: seedname
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: nrpts
    integer, intent(in) :: irvec(3, nrpts)
    integer:: irpt, iw, jw, ideg, file_unit, ierr
    character(len=100) :: header
    character(len=40) :: applied_token
    character(len=9)  :: cdate, ctime

    call io_date(cdate, ctime)

    applied_token = ''
    if (write_ndegen_applied) applied_token = '  write_ndegen_applied=.true.'

    open (newunit=file_unit, file=trim(seedname)//'_wsvec.dat', form='formatted', &
          status='unknown', iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error: ws_write_vec: problem opening file '//trim(seedname)//'_ws_vec.dat', comm)
      return
    end if

    if (use_ws_distance) then
      header = '## written on '//cdate//' at '//ctime//' with use_ws_distance=.true.'// &
               trim(applied_token)
      write (file_unit, '(A)') trim(header)

      do irpt = 1, nrpts
        do iw = 1, num_wann
          do jw = 1, num_wann
            write (file_unit, '(5I5)') irvec(:, irpt), iw, jw
            write (file_unit, '(I5)') ws_distance%ndeg(iw, jw, irpt)
            do ideg = 1, ws_distance%ndeg(iw, jw, irpt)
              write (file_unit, '(5I5,2F12.6,I5)') ws_distance%irdist(:, ideg, iw, jw, irpt) - &
                irvec(:, irpt)
            end do
          end do
        end do
      end do
    else
      header = '## written on '//cdate//' at '//ctime//' with use_ws_distance=.false.'// &
               trim(applied_token)
      write (file_unit, '(A)') trim(header)

      do irpt = 1, nrpts
        do iw = 1, num_wann
          do jw = 1, num_wann
            write (file_unit, '(5I5)') irvec(:, irpt), &
              iw, jw
            write (file_unit, '(I5)') 1
            write (file_unit, '(3I5)') 0, 0, 0
          end do
        end do
      end do
    end if

    close (file_unit)
    !================================================!
  end subroutine ws_write_vec

  !================================================!
  subroutine ws_expand_rvec(ws_distance, use_ws_distance, num_wann, nrpts, irvec, ndegen, &
                            irvec_full, nrpts_full, ir_map, ir_origin, error, comm)
    !================================================!
    !! Build the fully expanded list of lattice vectors, i.e. the set of all
    !! R + T that occur in the Wigner-Seitz mapping computed by ws_translate_dist,
    !! together with the index map ir_map(ideg, i, j, ir) that sends a degenerate
    !! image of the pair (i, j) at the folded vector irvec(:, ir) to its position
    !! in that list.
    !!
    !! The expanded list is ordered lexicographically, so that it does not depend
    !! on the order in which the vectors are discovered.
    !!
    !! If use_ws_distance is false there is nothing to expand: irvec_full is
    !! irvec and ir_map is unused.
    !================================================!

    use w90_types, only: ws_distance_type

    implicit none

    type(ws_distance_type), intent(in) :: ws_distance
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    logical, intent(in) :: use_ws_distance
    integer, intent(in) :: num_wann
    integer, intent(in) :: nrpts
    integer, intent(in) :: irvec(3, nrpts)
    integer, intent(in) :: ndegen(nrpts)

    integer, allocatable, intent(out) :: irvec_full(:, :)
    integer, intent(out) :: nrpts_full
    integer, allocatable, intent(out) :: ir_map(:, :, :, :)
    integer, intent(out) :: ir_origin
    !! index of R = 0 in the expanded list

    ! local variables
    integer :: i, j, ideg, ir, i1, i2, i3, ierr, max_ndeg, rpt_origin
    integer :: ivdum(3), ivmin(3), ivmax(3)
    integer, allocatable :: index_box(:, :, :)

    rpt_origin = 0
    do ir = 1, nrpts
      if (all(irvec(:, ir) == 0)) rpt_origin = ir
    end do
    if (rpt_origin == 0) then
      call set_error_fatal(error, 'R=0 is not in the list of lattice vectors.', comm)
      return
    end if
    if (ndegen(rpt_origin) /= 1) then
      call set_error_fatal(error, 'ndegen for R=0 is not 1.', comm)
      return
    end if

    if (.not. use_ws_distance) then
      nrpts_full = nrpts
      ir_origin = rpt_origin

      allocate (irvec_full(3, nrpts_full), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating irvec_full in ws_expand_rvec', comm)
        return
      end if
      ! ws_apply_ndegen does not consult ir_map when there is nothing to expand,
      ! so allocate it only to have something to pass
      allocate (ir_map(1, 1, 1, 1), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating ir_map in ws_expand_rvec', comm)
        return
      end if
      ir_map = -1

      irvec_full = irvec
      return
    end if

    ! Check degeneracy factor ws_distance%ndeg for a Wannier function with itself,
    ! i.e. R = 0 and i = j, is 1.
    do ir = 1, nrpts
      do i = 1, num_wann
        do ideg = 1, ws_distance%ndeg(i, i, ir)
          if (all(ws_distance%irdist(:, ideg, i, i, ir) == 0)) then
            if (ws_distance%ndeg(i, i, ir) /= 1) then
              call set_error_fatal(error, 'ws_distance%ndeg for R=0 and i=j is not 1.', comm)
              return
            end if
          end if
        end do
      end do
    end do

    max_ndeg = maxval(ws_distance%ndeg)

    ! Mark every vector that occurs in irdist on an integer box spanning them all,
    ! then walk the box in lexicographic order to number the vectors found.
    ! Unused slots of irdist are zero, which is a vector of the list anyway.
    ! R_wz_sc bounds the box by +-2*(ws_search_size + 1)*mp_grid, and in practice
    ! it is a small multiple of mp_grid: a few MB of integers at worst.
    do i = 1, 3
      ivmin(i) = minval(ws_distance%irdist(i, :, :, :, :))
      ivmax(i) = maxval(ws_distance%irdist(i, :, :, :, :))
    end do

    allocate (index_box(ivmin(1):ivmax(1), ivmin(2):ivmax(2), ivmin(3):ivmax(3)), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating index_box in ws_expand_rvec', comm)
      return
    end if
    index_box = 0

    nrpts_full = 0
    do ir = 1, nrpts
      do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            ivdum = ws_distance%irdist(:, ideg, i, j, ir)
            if (index_box(ivdum(1), ivdum(2), ivdum(3)) == 0) then
              index_box(ivdum(1), ivdum(2), ivdum(3)) = 1
              nrpts_full = nrpts_full + 1
            end if
          end do
        end do
      end do
    end do

    allocate (irvec_full(3, nrpts_full), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating irvec_full in ws_expand_rvec', comm)
      return
    end if
    allocate (ir_map(max_ndeg, num_wann, num_wann, nrpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating ir_map in ws_expand_rvec', comm)
      return
    end if

    ! each marked slot is visited once, so a slot still holding the mark 1 has
    ! not been numbered yet
    ir = 0
    do i1 = ivmin(1), ivmax(1)
      do i2 = ivmin(2), ivmax(2)
        do i3 = ivmin(3), ivmax(3)
          if (index_box(i1, i2, i3) == 1) then
            ir = ir + 1
            index_box(i1, i2, i3) = ir
            irvec_full(:, ir) = (/i1, i2, i3/)
          end if
        end do
      end do
    end do
    ir_origin = index_box(0, 0, 0)

    ir_map = -1
    do ir = 1, nrpts
      do j = 1, num_wann
        do i = 1, num_wann
          do ideg = 1, ws_distance%ndeg(i, j, ir)
            ivdum = ws_distance%irdist(:, ideg, i, j, ir)
            ir_map(ideg, i, j, ir) = index_box(ivdum(1), ivdum(2), ivdum(3))
          end do
        end do
      end do
    end do

    deallocate (index_box, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating index_box in ws_expand_rvec', comm)
      return
    end if
    !================================================!
  end subroutine ws_expand_rvec

  !================================================!
  subroutine ws_apply_ndegen(ws_distance, use_ws_distance, num_wann, nrpts, ndegen, &
                             nrpts_full, ir_map, op_R, op_R_full)
    !================================================!
    !! Divide a real-space operator by its degeneracy weights and spread it over
    !! the expanded lattice-vector list built by ws_expand_rvec, so that it can be
    !! Fourier transformed with a plain sum over exp(i k.R), irrespective of
    !! use_ws_distance.
    !================================================!

    use w90_constants, only: cmplx_0
    use w90_types, only: ws_distance_type

    implicit none

    type(ws_distance_type), intent(in) :: ws_distance

    logical, intent(in) :: use_ws_distance
    integer, intent(in) :: num_wann
    integer, intent(in) :: nrpts
    integer, intent(in) :: ndegen(nrpts)
    integer, intent(in) :: nrpts_full
    integer, intent(in) :: ir_map(:, :, :, :)

    complex(kind=dp), intent(in) :: op_R(num_wann, num_wann, nrpts)
    !! operator on the folded grid, before applying the degeneracy weights
    complex(kind=dp), intent(out) :: op_R_full(num_wann, num_wann, nrpts_full)
    !! operator on the expanded grid, after applying the degeneracy weights

    integer :: ir, jr, i, j, ideg

    if (use_ws_distance) then
      op_R_full = cmplx_0
      do ir = 1, nrpts
        do j = 1, num_wann
          do i = 1, num_wann
            do ideg = 1, ws_distance%ndeg(i, j, ir)
              jr = ir_map(ideg, i, j, ir)
              op_R_full(i, j, jr) = op_R_full(i, j, jr) &
                                    + op_R(i, j, ir)/real(ndegen(ir)*ws_distance%ndeg(i, j, ir), dp)
            end do
          end do
        end do
      end do
    else
      ! nrpts_full == nrpts in this case
      do ir = 1, nrpts
        op_R_full(:, :, ir) = op_R(:, :, ir)/real(ndegen(ir), dp)
      end do
    end if
    !================================================!
  end subroutine ws_apply_ndegen

  !================================================!
  subroutine clean_ws_translate(ws_distance, error, comm)
    !================================================!
    use w90_types, only: ws_distance_type
    use w90_comms, only: w90_comm_type
    use w90_error, only: w90_error_type, set_error_dealloc

    implicit none

    type(ws_distance_type), intent(inout) :: ws_distance
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer :: ierr

    ws_distance%done = .false.
    if (allocated(ws_distance%irdist)) then
      deallocate (ws_distance%irdist, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ws_distance%irdist in clean_ws_translate', comm)
        return
      end if
    end if
    if (allocated(ws_distance%ndeg)) then
      deallocate (ws_distance%ndeg, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ws_distance%ndeg in clean_ws_translate', comm)
        return
      end if
    end if
    if (allocated(ws_distance%crdist)) then
      deallocate (ws_distance%crdist, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating ws_distance%crdist in clean_ws_translate', comm)
        return
      end if
    end if

    !================================================!
  end subroutine clean_ws_translate

end module w90_ws_distance
