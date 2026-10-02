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
! w90_sym_ibz:                                               !
! Reads the overlap, projection and eigenvalue matrices in   !
! the irreducible Brillouin zone (the ".immn", ".iamn" and   !
! ".ieig" files written by pw2wannier90.x with irr_bz=.true.,!
! together with the symmetry data in ".isym") and expands    !
! them onto the full BZ, so that no ".mmn", ".amn" or ".eig" !
! file is needed.                                            !
! Reference:                                                 !
!    T. Koretsune, Comput. Phys. Commun. 285, 108645 (2023)  !
!                                                            !
! Note on the symmetry operations:                           !
!    s(:,:,isym) . k          acts on reciprocal space       !
!    r . s(:,:,isym) - ft     acts on real space             !
! both in crystal coordinates.                               !
!                                                            !
!------------------------------------------------------------!

module w90_sym_ibz

  use w90_comms, only: w90_comm_type, comms_bcast, mpirank
  use w90_constants, only: dp, cmplx_0, cmplx_1, cmplx_i, twopi, eps6, eps8
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
                       set_error_file, set_error_input
  use w90_io, only: io_date
  use w90_types, only: kmesh_info_type, print_output_type, proj_type
  use w90_wannier90_types, only: sym_ibz_type

  implicit none

  private

  public :: sym_ibz_dealloc
  public :: sym_ibz_read
  public :: sym_ibz_read_eigvals
  public :: sym_ibz_read_overlaps

  !! tolerance used when matching k-points, b-vectors and lattice translations
  real(kind=dp), parameter :: tol_k = 1.0e-5_dp
  !! tolerance used when matching spinor rotation matrices
  real(kind=dp), parameter :: tol_u = 1.0e-5_dp
  !! tolerance used when matching the Wannier-centre shifts to a lattice vector
  real(kind=dp), parameter :: tol_shift = 1.0e-3_dp

contains

  !================================================!
  subroutine sym_ibz_read(sym, kpt_latt, mp_grid, num_bands, num_kpts, num_wann, print_output, &
                          seedname, stdout, error, comm)
    !================================================!
    !! Read the ".isym" file and work out the relation between the irreducible
    !! and the full k-point meshes.  Returns immediately if this was done before.
    !================================================!

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(sym_ibz_type), intent(inout) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt_latt(:, :)

    character(len=*), intent(in) :: seedname

    if (sym%ready) return ! already set up

    ! a previous, failed attempt may have left partially filled arrays behind
    call sym_ibz_reset(sym)

    call read_sym_file(sym, seedname, print_output, stdout, error, comm)

    if (.not. allocated(error)) then
      if (sym%nbnd /= num_bands) then
        call set_error_file(error, 'Error: '//trim(seedname)//'.isym has not the right number of bands', comm)
      else if (sym%num_wann /= num_wann) then
        call set_error_file(error, 'Error: '//trim(seedname)//'.isym has not the right number of projections', comm)
      end if
    end if

    if (.not. allocated(error)) call kpoint_map(sym, kpt_latt, mp_grid, num_kpts, seedname, error, comm)

    if (allocated(error)) then
      ! leave nothing half-initialised behind: a later call must not find
      ! (and use) an incomplete k-point map, and a retry must be able to
      ! allocate the arrays again
      call sym_ibz_reset(sym)
      return
    end if
    sym%ready = .true.

    if (print_output%iprint > 0) then
      write (stdout, '(1x,a,i0,a,i0,a)') 'Symmetry from '//trim(seedname)//'.isym: ', sym%nsym, &
        ' operations, ', sym%nks, ' irreducible k-points'
    end if
  end subroutine sym_ibz_read

  !================================================!
  subroutine sym_ibz_read_eigvals(sym, eigval, num_bands, num_kpts, print_output, seedname, &
                                  stdout, write_expanded, error, comm)
    !================================================!
    !! Read the eigenvalues of the irreducible k-points from ".ieig" and copy
    !! them onto the full k-point mesh; with write_expanded, also write the
    !! result to ".eig".
    !! sym_ibz_read() must have been called first.
    !================================================!

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: stdout

    logical, intent(in) :: write_expanded

    real(kind=dp), intent(inout) :: eigval(:, :)

    character(len=*), intent(in) :: seedname

    ! local variables
    real(kind=dp), allocatable :: eig_irr(:, :)
    integer :: eig_in, ierr, ik, iks, i, j, n

    if (.not. sym%ready) then
      call set_error_fatal(error, 'Error: sym_ibz_read must be called before sym_ibz_read_eigvals', comm)
      return
    end if

    allocate (eig_irr(num_bands, sym%nks), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating eig_irr in sym_ibz_read_eigvals', comm)
      return
    end if

    open (newunit=eig_in, file=trim(seedname)//'.ieig', form='formatted', status='old', &
          action='read', iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error: Problem opening input file '//trim(seedname)//'.ieig', comm)
      return
    end if

    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      'Reading eigenvalues from '//trim(seedname)//'.ieig'

    do iks = 1, sym%nks
      do n = 1, num_bands
        read (eig_in, *, iostat=ierr) i, j, eig_irr(n, iks)
        if (ierr /= 0) then
          close (eig_in)
          call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.ieig', comm)
          return
        end if
        if ((i /= n) .or. (j /= iks)) then
          close (eig_in)
          call set_error_file(error, 'Error: Found a mismatch in '//trim(seedname)//'.ieig', comm)
          return
        end if
      end do
    end do
    close (eig_in)

    do ik = 1, num_kpts
      eigval(:, ik) = eig_irr(:, sym%equiv(ik))
    end do

    if (write_expanded) then
      call write_full_eig(eigval, num_bands, num_kpts, print_output, seedname, stdout, error, comm)
      if (allocated(error)) return
    end if

    deallocate (eig_irr, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating eig_irr in sym_ibz_read_eigvals', comm)
      return
    end if
  end subroutine sym_ibz_read_eigvals

  !================================================!
  subroutine sym_ibz_read_overlaps(sym, kmesh_info, kpt_latt, proj_input, au_matrix, &
                                   m_matrix_local, num_bands, num_kpts, num_proj, num_wann, &
                                   print_output, use_bloch_phases, seedname, stdout, dist_k, &
                                   write_expanded, error, comm)
    !================================================!
    !! Read the overlap and projection matrices of the irreducible k-points
    !! from ".immn" and ".iamn" and expand them onto the full k-point mesh.
    !! With write_expanded, the expanded matrices are also written to ".mmn"
    !! and ".amn" (the latter not with use_bloch_phases, where no projections
    !! are read).
    !! sym_ibz_read() must have been called first.
    !================================================!

    implicit none

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(proj_type), allocatable, intent(in) :: proj_input(:)
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: dist_k(:)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_proj
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    complex(kind=dp), intent(inout) :: au_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :)

    logical, intent(in) :: use_bloch_phases
    logical, intent(in) :: write_expanded

    real(kind=dp), intent(in) :: kpt_latt(:, :)

    character(len=*), intent(in) :: seedname

    ! local variables
    complex(kind=dp), allocatable :: amn_irr(:, :, :), mmn_irr(:, :, :, :)
    real(kind=dp), allocatable :: bvec(:, :), pos(:, :), rshift(:, :, :)
    integer, allocatable :: bequiv(:, :), ib_of_nn(:, :), nn_of_ib(:, :), map_kpts(:)
    integer :: ierr, ik, iw, mmn_out, nb, nkp_loc, rank
    character(len=9) :: cdate, ctime

    if (.not. sym%ready) then
      call set_error_fatal(error, 'Error: sym_ibz_read must be called before sym_ibz_read_overlaps', comm)
      return
    end if

    if (num_proj /= num_wann) then
      call set_error_input(error, &
                           'Error: read_ibz requires as many projections as Wannier functions', comm)
      return
    end if

    nb = kmesh_info%nntot
    rank = mpirank(comm)

    ! index of each k-point within the rank-local part of m_matrix_local
    allocate (map_kpts(num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating map_kpts in sym_ibz_read_overlaps', comm)
      return
    end if
    map_kpts = 0
    nkp_loc = 1
    do ik = 1, num_kpts
      if (dist_k(ik) == rank) then
        map_kpts(ik) = nkp_loc
        nkp_loc = nkp_loc + 1
      end if
    end do

    call bvector_setup(kmesh_info, kpt_latt, num_kpts, bvec, ib_of_nn, nn_of_ib, error, comm)
    if (allocated(error)) return

    call bvector_equiv(sym, bvec, nb, bequiv, error, comm)
    if (allocated(error)) return

    ! Mmn
    allocate (mmn_irr(num_bands, num_bands, nb, sym%nks), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating mmn_irr in sym_ibz_read_overlaps', comm)
      return
    end if

    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      'Reading overlaps from '//trim(seedname)//'.immn'

    call read_immn(sym, kmesh_info, kpt_latt, nn_of_ib, mmn_irr, num_bands, nb, seedname, error, comm)
    if (allocated(error)) return

    mmn_out = 0
    if (write_expanded) then
      call open_expanded_file(mmn_out, trim(seedname)//'.mmn', error, comm)
      if (allocated(error)) return
      if (mmn_out /= 0) then
        call io_date(cdate, ctime)
        write (mmn_out, '(a)') 'Expanded from '//trim(seedname)//'.immn by wannier90.x (read_ibz) on '// &
          cdate//' at '//ctime
        write (mmn_out, '(3i12)') num_bands, num_kpts, nb
      end if
      if (print_output%iprint > 0) write (stdout, '(1x,a)') &
        'Writing the expanded overlaps to '//trim(seedname)//'.mmn'
    end if

    call expand_mmn(sym, kmesh_info, mmn_irr, m_matrix_local, bvec, bequiv, ib_of_nn, nn_of_ib, &
                    map_kpts, num_bands, num_kpts, nb, dist_k, rank, mmn_out, error, comm)
    if (mmn_out /= 0) close (mmn_out)
    if (allocated(error)) return

    deallocate (mmn_irr, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating mmn_irr in sym_ibz_read_overlaps', comm)
      return
    end if

    ! Amn
    if (use_bloch_phases) then
      au_matrix = cmplx_0
      do ik = 1, num_kpts
        do iw = 1, num_wann
          au_matrix(iw, iw, ik) = cmplx_1
        end do
      end do
    else
      allocate (pos(3, num_wann), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating pos in sym_ibz_read_overlaps', comm)
        return
      end if
      ! the centres of the projections are needed to fix the phases of the
      ! rotated projection functions; they come from the projections block
      if (.not. allocated(proj_input)) then
        call set_error_input(error, 'Error: read_ibz needs the projections block of the input file', comm)
        return
      end if
      if (size(proj_input) < num_wann) then
        call set_error_input(error, 'Error: read_ibz: too few projections defined in the input file', comm)
        return
      end if
      do iw = 1, num_wann
        pos(:, iw) = proj_input(iw)%site(:)
      end do

      call projection_shifts(sym, pos, num_wann, rshift, error, comm)
      if (allocated(error)) return

      allocate (amn_irr(num_bands, num_wann, sym%nks), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating amn_irr in sym_ibz_read_overlaps', comm)
        return
      end if

      if (print_output%iprint > 0) write (stdout, '(1x,a)') &
        'Reading projections from '//trim(seedname)//'.iamn'

      call read_iamn(sym, amn_irr, num_bands, num_wann, seedname, error, comm)
      if (allocated(error)) return

      call symmetrise_amn(sym, amn_irr, rshift, num_bands, num_wann, error, comm)
      if (allocated(error)) return

      call expand_amn(sym, amn_irr, au_matrix, rshift, num_bands, num_kpts, num_wann)

      if (write_expanded) then
        call write_full_amn(au_matrix, num_bands, num_kpts, num_wann, print_output, seedname, &
                            stdout, error, comm)
        if (allocated(error)) return
      end if

      deallocate (amn_irr, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error in deallocating amn_irr in sym_ibz_read_overlaps', comm)
        return
      end if
    end if
  end subroutine sym_ibz_read_overlaps

  !================================================!
  subroutine sym_ibz_dealloc(sym, error, comm)
    !================================================!
    !! Release the symmetry data read from the ".isym" file
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(inout) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    ! local variables
    integer :: ierr

    sym%ready = .false.
    ierr = 0
    if (allocated(sym%s)) deallocate (sym%s, stat=ierr)
    if (allocated(sym%t_rev)) deallocate (sym%t_rev, stat=ierr)
    if (allocated(sym%invs)) deallocate (sym%invs, stat=ierr)
    if (allocated(sym%ft)) deallocate (sym%ft, stat=ierr)
    if (allocated(sym%irr_kpt)) deallocate (sym%irr_kpt, stat=ierr)
    if (allocated(sym%u_spin)) deallocate (sym%u_spin, stat=ierr)
    if (allocated(sym%repmat)) deallocate (sym%repmat, stat=ierr)
    if (allocated(sym%rotmat)) deallocate (sym%rotmat, stat=ierr)
    if (allocated(sym%equiv)) deallocate (sym%equiv, stat=ierr)
    if (allocated(sym%equiv_sym)) deallocate (sym%equiv_sym, stat=ierr)
    if (allocated(sym%iks2ik)) deallocate (sym%iks2ik, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating sym in sym_ibz_dealloc', comm)
      return
    end if
  end subroutine sym_ibz_dealloc

  !================================================!
  !                private procedures              !
  !================================================!

  !================================================!
  subroutine sym_ibz_reset(sym)
    !================================================!
    !! Discard whatever sym_ibz_read() has stored so far, without reporting
    !! deallocation errors.  Used to clean up after a failed read, when an
    !! error has already been set.
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(inout) :: sym

    ! local variables
    integer :: ierr

    sym%ready = .false.
    sym%nsym = 0
    sym%nks = 0
    sym%nbnd = 0
    sym%num_wann = 0
    if (allocated(sym%s)) deallocate (sym%s, stat=ierr)
    if (allocated(sym%t_rev)) deallocate (sym%t_rev, stat=ierr)
    if (allocated(sym%invs)) deallocate (sym%invs, stat=ierr)
    if (allocated(sym%ft)) deallocate (sym%ft, stat=ierr)
    if (allocated(sym%irr_kpt)) deallocate (sym%irr_kpt, stat=ierr)
    if (allocated(sym%u_spin)) deallocate (sym%u_spin, stat=ierr)
    if (allocated(sym%repmat)) deallocate (sym%repmat, stat=ierr)
    if (allocated(sym%rotmat)) deallocate (sym%rotmat, stat=ierr)
    if (allocated(sym%equiv)) deallocate (sym%equiv, stat=ierr)
    if (allocated(sym%equiv_sym)) deallocate (sym%equiv_sym, stat=ierr)
    if (allocated(sym%iks2ik)) deallocate (sym%iks2ik, stat=ierr)
  end subroutine sym_ibz_reset

  !================================================!
  subroutine read_sym_file(sym, seedname, print_output, stdout, error, comm)
    !================================================!
    !! Read the symmetry operations, the irreducible k-points, the
    !! representation matrices of the little groups and the rotation matrices
    !! of the projections from the ".isym" file
    !================================================!

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(sym_ibz_type), intent(inout) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: stdout

    character(len=*), intent(in) :: seedname

    ! local variables
    real(kind=dp) :: rr, ri
    integer :: a, b, i, ib, ierr, ik, iks, isym, isym_in, m, n, nblocks, nlines, sym_in
    character(len=256) :: line

    open (newunit=sym_in, file=trim(seedname)//'.isym', form='formatted', status='old', &
          action='read', iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error: Problem opening input file '//trim(seedname)//'.isym', comm)
      return
    end if

    ! header line
    read (sym_in, '(a)', iostat=ierr) line
    if (ierr /= 0) goto 100
    if (print_output%iprint > 0) write (stdout, '(1x,a)') trim(adjustl(line))

    read (sym_in, *, iostat=ierr) sym%nsym, i
    if (ierr /= 0) goto 100
    sym%spinors = (i == 1)
    if (sym%nsym < 1) then
      close (sym_in)
      call set_error_file(error, 'Error: no symmetry operation in '//trim(seedname)//'.isym', comm)
      return
    end if

    allocate (sym%s(3, 3, sym%nsym), sym%ft(3, sym%nsym), sym%t_rev(sym%nsym), &
              sym%invs(sym%nsym), sym%u_spin(2, 2, sym%nsym), stat=ierr)
    if (ierr /= 0) then
      close (sym_in)
      call set_error_alloc(error, 'Error in allocating symmetry operations in read_sym_file', comm)
      return
    end if
    sym%u_spin = cmplx_0

    do isym = 1, sym%nsym
      read (sym_in, '(a)', iostat=ierr) line ! name of the operation
      if (ierr /= 0) goto 100
      do a = 1, 3
        read (sym_in, *, iostat=ierr) (sym%s(a, b, isym), b=1, 3)
        if (ierr /= 0) goto 100
      end do
      read (sym_in, *, iostat=ierr) (sym%ft(a, isym), a=1, 3)
      if (ierr /= 0) goto 100
      read (sym_in, *, iostat=ierr) sym%t_rev(isym)
      if (ierr /= 0) goto 100
      if (sym%spinors) then
        do a = 1, 2
          do b = 1, 2
            read (sym_in, *, iostat=ierr) rr, ri
            if (ierr /= 0) goto 100
            sym%u_spin(a, b, isym) = cmplx(rr, ri, kind=dp)
          end do
        end do
      else
        sym%u_spin(1, 1, isym) = cmplx_1
        sym%u_spin(2, 2, isym) = cmplx_1
      end if
      read (sym_in, *, iostat=ierr) sym%invs(isym)
      if (ierr /= 0) goto 100
      if (sym%invs(isym) < 1 .or. sym%invs(isym) > sym%nsym) then
        close (sym_in)
        call set_error_file(error, 'Error: bad inverse operation index in '//trim(seedname)//'.isym', comm)
        return
      end if
    end do

    ! irreducible k-points
    call next_nonblank(sym_in, line, ierr) ! section title
    if (ierr /= 0) goto 100
    read (sym_in, *, iostat=ierr) sym%nks
    if (ierr /= 0) goto 100
    if (sym%nks < 1) then
      close (sym_in)
      call set_error_file(error, 'Error: no irreducible k-point in '//trim(seedname)//'.isym', comm)
      return
    end if
    allocate (sym%irr_kpt(3, sym%nks), stat=ierr)
    if (ierr /= 0) then
      close (sym_in)
      call set_error_alloc(error, 'Error in allocating irr_kpt in read_sym_file', comm)
      return
    end if
    do iks = 1, sym%nks
      read (sym_in, *, iostat=ierr) (sym%irr_kpt(a, iks), a=1, 3)
      if (ierr /= 0) goto 100
    end do

    ! representation matrices of the little group of each irreducible k-point
    call next_nonblank(sym_in, line, ierr) ! section title
    if (ierr /= 0) goto 100
    read (sym_in, *, iostat=ierr) sym%nbnd, nblocks
    if (ierr /= 0) goto 100
    allocate (sym%repmat(sym%nbnd, sym%nbnd, sym%nsym, sym%nks), stat=ierr)
    if (ierr /= 0) then
      close (sym_in)
      call set_error_alloc(error, 'Error in allocating repmat in read_sym_file', comm)
      return
    end if
    sym%repmat = cmplx_0
    do ib = 1, nblocks
      read (sym_in, *, iostat=ierr) ik, isym, nlines
      if (ierr /= 0) goto 100
      if (ik < 1 .or. ik > sym%nks .or. isym < 1 .or. isym > sym%nsym) then
        close (sym_in)
        call set_error_file(error, 'Error: bad representation matrix index in '//trim(seedname)//'.isym', comm)
        return
      end if
      do i = 1, nlines
        read (sym_in, *, iostat=ierr) m, n, rr, ri
        if (ierr /= 0) goto 100
        if (m < 1 .or. m > sym%nbnd .or. n < 1 .or. n > sym%nbnd) then
          close (sym_in)
          call set_error_file(error, 'Error: bad representation matrix index in '//trim(seedname)//'.isym', comm)
          return
        end if
        sym%repmat(m, n, isym, ik) = cmplx(rr, ri, kind=dp)
      end do
    end do

    ! rotation matrices of the projection functions
    call next_nonblank(sym_in, line, ierr) ! section title
    if (ierr /= 0) goto 100
    read (sym_in, *, iostat=ierr) sym%num_wann
    if (ierr /= 0) goto 100
    allocate (sym%rotmat(sym%num_wann, sym%num_wann, sym%nsym), stat=ierr)
    if (ierr /= 0) then
      close (sym_in)
      call set_error_alloc(error, 'Error in allocating rotmat in read_sym_file', comm)
      return
    end if
    sym%rotmat = cmplx_0
    do isym = 1, sym%nsym
      read (sym_in, *, iostat=ierr) isym_in, nlines
      if (ierr /= 0) goto 100
      if (isym_in /= isym) then
        close (sym_in)
        call set_error_file(error, 'Error: bad rotation matrix index in '//trim(seedname)//'.isym', comm)
        return
      end if
      do i = 1, nlines
        read (sym_in, *, iostat=ierr) m, n, rr, ri
        if (ierr /= 0) goto 100
        if (m < 1 .or. m > sym%num_wann .or. n < 1 .or. n > sym%num_wann) then
          close (sym_in)
          call set_error_file(error, 'Error: bad rotation matrix index in '//trim(seedname)//'.isym', comm)
          return
        end if
        sym%rotmat(m, n, isym) = cmplx(rr, ri, kind=dp)
      end do
    end do

    close (sym_in)

    call repmat_rescale(sym)

    return

100 close (sym_in)
    call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.isym', comm)
    return
  end subroutine read_sym_file

  !================================================!
  subroutine next_nonblank(unit, line, ierr)
    !================================================!
    !! Read forward until a line which is not blank; used to step over the
    !! blank separator lines and pick up the section titles of the ".isym" file
    !================================================!

    implicit none

    ! arguments
    integer, intent(in) :: unit
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: line

    do
      read (unit, '(a)', iostat=ierr) line
      if (ierr /= 0) return
      if (len_trim(line) > 0) return
    end do
  end subroutine next_nonblank

  !================================================!
  subroutine repmat_rescale(sym)
    !================================================!
    !! Rescale the diagonal entries of the representation matrices which are
    !! diagonal in a given band but whose eigenvalue is not of unit modulus
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(inout) :: sym

    ! local variables
    real(kind=dp) :: val, vall
    integer :: iks, isym, n

    do iks = 1, sym%nks
      do isym = 1, sym%nsym
        do n = 1, sym%nbnd
          vall = sum(abs(sym%repmat(n, :, isym, iks)))
          val = abs(sym%repmat(n, n, isym, iks))
          if (abs(vall - val) < eps8 .and. nint(val) /= 0) then
            sym%repmat(n, n, isym, iks) = sym%repmat(n, n, isym, iks)*real(nint(val), dp)/val
          end if
        end do
      end do
    end do
  end subroutine repmat_rescale

  !================================================!
  subroutine kpoint_map(sym, kpt_latt, mp_grid, num_kpts, seedname, error, comm)
    !================================================!
    !! Relate the full k-point mesh to the irreducible one:
    !!   s(:,:,equiv_sym(ik)) . irr_kpt(:,equiv(ik)) = kpt_latt(:,ik)  (mod G)
    !! and iks2ik(iks) gives the full-mesh index of irr_kpt(:,iks)
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(inout) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_kpts

    real(kind=dp), intent(in) :: kpt_latt(:, :)

    character(len=*), intent(in) :: seedname

    ! local variables
    integer, allocatable :: kmap(:, :, :)
    integer :: idx(3), ierr, ik, iks, isym
    real(kind=dp) :: sk(3)

    if (num_kpts /= mp_grid(1)*mp_grid(2)*mp_grid(3)) then
      call set_error_input(error, 'Error: read_ibz requires a complete Monkhorst-Pack k-point mesh', comm)
      return
    end if

    allocate (kmap(0:mp_grid(1) - 1, 0:mp_grid(2) - 1, 0:mp_grid(3) - 1), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating kmap in kpoint_map', comm)
      return
    end if
    kmap = 0

    do ik = 1, num_kpts
      call grid_index(kpt_latt(:, ik), mp_grid, idx, ierr)
      if (ierr /= 0) then
        call set_error_input(error, &
                             'Error: read_ibz requires the k-points to lie on the Gamma-centred mp_grid mesh', comm)
        return
      end if
      if (kmap(idx(1), idx(2), idx(3)) /= 0) then
        call set_error_input(error, 'Error: read_ibz: two k-points of the mesh coincide', comm)
        return
      end if
      kmap(idx(1), idx(2), idx(3)) = ik
    end do

    allocate (sym%equiv(num_kpts), sym%equiv_sym(num_kpts), sym%iks2ik(sym%nks), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating equiv in kpoint_map', comm)
      return
    end if
    sym%equiv = 0
    sym%equiv_sym = 0
    sym%iks2ik = 0

    do iks = 1, sym%nks
      call grid_index(sym%irr_kpt(:, iks), mp_grid, idx, ierr)
      if (ierr /= 0) then
        call set_error_file(error, 'Error: an irreducible k-point of '//trim(seedname)// &
                            '.isym is not on the mp_grid mesh', comm)
        return
      end if
      sym%iks2ik(iks) = kmap(idx(1), idx(2), idx(3))

      do isym = 1, sym%nsym
        sk = matmul(real(sym%s(:, :, isym), dp), sym%irr_kpt(:, iks))
        if (sym%t_rev(isym) == 1) sk = -sk
        call grid_index(sk, mp_grid, idx, ierr)
        if (ierr /= 0) then
          call set_error_file(error, 'Error: a symmetry operation of '//trim(seedname)// &
                              '.isym does not map the k-point mesh onto itself', comm)
          return
        end if
        ik = kmap(idx(1), idx(2), idx(3))
        ! Keep the FIRST operation (in the order of the .isym file) which maps
        ! irr_kpt onto ik.  This is not a free gauge choice: pw2wannier90.x
        ! generates the Bloch states at the non-irreducible neighbours k_irr+b
        ! of the .immn file by the same rule, so equiv_sym must reproduce it
        ! (see expand_mmn).  Choosing a different operation would silently give
        ! wrong overlaps.
        if (sym%equiv(ik) == 0) then
          sym%equiv(ik) = iks
          sym%equiv_sym(ik) = isym
        else if (sym%equiv(ik) /= iks) then
          call set_error_file(error, &
                              'Error: read_ibz: a k-point is equivalent to two different irreducible k-points', comm)
          return
        end if
      end do
    end do

    if (any(sym%equiv == 0)) then
      call set_error_file(error, &
                          'Error: read_ibz: some k-points have no equivalent irreducible k-point', comm)
      return
    end if
  end subroutine kpoint_map

  !================================================!
  subroutine grid_index(kpt, mp_grid, idx, ierr)
    !================================================!
    !! Index of a k-point (crystal coordinates) on the Gamma-centred mp_grid
    !! mesh; ierr is non-zero if the k-point is not a mesh point
    !================================================!

    implicit none

    ! arguments
    integer, intent(in) :: mp_grid(3)
    integer, intent(out) :: idx(3)
    integer, intent(out) :: ierr

    real(kind=dp), intent(in) :: kpt(3)

    ! local variables
    integer :: i
    real(kind=dp) :: x

    ierr = 0
    do i = 1, 3
      x = kpt(i)*real(mp_grid(i), dp)
      idx(i) = nint(x)
      if (abs(x - real(idx(i), dp)) > 1.0e-4_dp) then
        ierr = 1
        return
      end if
      idx(i) = modulo(idx(i), mp_grid(i))
    end do
  end subroutine grid_index

  !================================================!
  subroutine bvector_setup(kmesh_info, kpt_latt, num_kpts, bvec, ib_of_nn, nn_of_ib, error, comm)
    !================================================!
    !! Set up the global list of b-vectors (crystal coordinates) in the order
    !! in which they appear at the first k-point.  This is the order in which
    !! the ".nnkp" file lists them, and hence the order in which pw2wannier90.x
    !! writes the blocks of the ".immn" file.  ib_of_nn and nn_of_ib relate this
    !! global order to the neighbour order of each k-point.
    !================================================!

    implicit none

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_kpts
    integer, allocatable, intent(out) :: ib_of_nn(:, :), nn_of_ib(:, :)

    real(kind=dp), intent(in) :: kpt_latt(:, :)
    real(kind=dp), allocatable, intent(out) :: bvec(:, :)

    ! local variables
    real(kind=dp) :: b(3)
    integer :: ib, ierr, ik, nb, nn

    nb = kmesh_info%nntot

    allocate (bvec(3, nb), ib_of_nn(nb, num_kpts), nn_of_ib(nb, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating bvec in bvector_setup', comm)
      return
    end if

    do nn = 1, nb
      bvec(:, nn) = kpt_latt(:, kmesh_info%nnlist(1, nn)) - kpt_latt(:, 1) &
                    + real(kmesh_info%nncell(:, 1, nn), dp)
    end do

    ib_of_nn = 0
    nn_of_ib = 0
    do ik = 1, num_kpts
      do nn = 1, nb
        b(:) = kpt_latt(:, kmesh_info%nnlist(ik, nn)) - kpt_latt(:, ik) &
               + real(kmesh_info%nncell(:, ik, nn), dp)
        ib = find_bvec(bvec, nb, b)
        if (ib == 0) then
          call set_error_fatal(error, 'Error: read_ibz: b-vector shells are not the same at every k-point', comm)
          return
        end if
        if (nn_of_ib(ib, ik) /= 0) then
          call set_error_fatal(error, 'Error: read_ibz: duplicate b-vector at a k-point', comm)
          return
        end if
        ib_of_nn(nn, ik) = ib
        nn_of_ib(ib, ik) = nn
      end do
    end do
  end subroutine bvector_setup

  !================================================!
  function find_bvec(bvec, nb, b) result(ib)
    !================================================!
    !! Index of b in the list bvec, or zero if it is not in the list
    !================================================!

    implicit none

    ! arguments
    integer, intent(in) :: nb
    real(kind=dp), intent(in) :: b(3)
    real(kind=dp), intent(in) :: bvec(:, :)
    integer :: ib

    ! local variables
    integer :: i

    ib = 0
    do i = 1, nb
      if (all(abs(bvec(:, i) - b(:)) < tol_k)) then
        ib = i
        return
      end if
    end do
  end function find_bvec

  !================================================!
  subroutine bvector_equiv(sym, bvec, nb, bequiv, error, comm)
    !================================================!
    !! bequiv(ib2,isym) is the index ib1 for which
    !!   s(:,:,isym) . bvec(:,ib1) = bvec(:,ib2)
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: nb
    integer, allocatable, intent(out) :: bequiv(:, :)

    real(kind=dp), intent(in) :: bvec(:, :)

    ! local variables
    real(kind=dp) :: sb(3)
    integer :: ib1, ib2, ierr, isym

    allocate (bequiv(nb, sym%nsym), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating bequiv in bvector_equiv', comm)
      return
    end if
    bequiv = 0

    do ib1 = 1, nb
      do isym = 1, sym%nsym
        sb = matmul(real(sym%s(:, :, isym), dp), bvec(:, ib1))
        if (sym%t_rev(isym) == 1) sb = -sb
        ib2 = find_bvec(bvec, nb, sb)
        if (ib2 > 0) bequiv(ib2, isym) = ib1
      end do
    end do

    if (any(bequiv == 0)) then
      call set_error_fatal(error, &
                           'Error: read_ibz: the set of b-vectors is not invariant under the crystal symmetry', comm)
      return
    end if
  end subroutine bvector_equiv

  !================================================!
  subroutine search_symop(sym, list_sym, list_inv, nop, isym_out, factor, tdiff, error, comm)
    !================================================!
    !! Multiply the operations listed in list_sym (each of them inverted when
    !! the corresponding entry of list_inv is negative) and return the index of
    !! the equivalent symmetry operation of the crystal.  factor is the sign
    !! relating the spinor rotation matrices, and tdiff is the residual lattice
    !! translation.
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: nop
    integer, intent(in) :: list_sym(nop), list_inv(nop)
    integer, intent(out) :: isym_out

    real(kind=dp), intent(out) :: factor
    real(kind=dp), intent(out) :: tdiff(3)

    ! local variables
    complex(kind=dp), parameter :: uspin_t(2, 2) = &
                                   reshape((/(0.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp), &
                                             (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)/), (/2, 2/))
    complex(kind=dp), parameter :: uspin_t_inv(2, 2) = &
                                   reshape((/(0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), &
                                             (-1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)/), (/2, 2/))

    complex(kind=dp) :: u0(2, 2), u1(2, 2), us(2, 2)
    real(kind=dp) :: t0(3), t1(3), td(3)
    integer :: s0(3, 3), s1(3, 3)
    integer :: i, iop, isym, trev

    isym_out = 0
    factor = 1.0_dp
    tdiff = 0.0_dp

    s0 = 0
    u0 = cmplx_0
    do i = 1, 3
      s0(i, i) = 1
    end do
    u0(1, 1) = cmplx_1
    u0(2, 2) = cmplx_1
    t0 = 0.0_dp
    trev = 0

    do iop = 1, nop
      isym = list_sym(iop)
      if (list_inv(iop) > 0) then
        ! r . s - ft
        s1 = sym%s(:, :, isym)
        t1 = sym%ft(:, isym)
        u1 = sym%u_spin(:, :, isym)
        if (sym%t_rev(isym) == 1) u1 = matmul(uspin_t, conjg(u1))
      else
        ! the inverse: (r + ft) . s^-1
        s1 = sym%s(:, :, sym%invs(isym))
        t1 = -matmul(sym%ft(:, isym), real(s1, dp))
        u1 = conjg(transpose(sym%u_spin(:, :, isym)))
        if (sym%t_rev(isym) == 1) u1 = matmul(u1, uspin_t_inv)
      end if

      s0 = matmul(s0, s1)
      t0 = matmul(t0, real(s1, dp)) + t1
      if (trev == 1) u1 = conjg(u1)
      u0 = matmul(u0, u1)
      trev = modulo(trev + sym%t_rev(isym), 2)
    end do

    do isym = 1, sym%nsym
      if (any(s0 /= sym%s(:, :, isym))) cycle
      if (trev /= sym%t_rev(isym)) cycle
      td = t0 - sym%ft(:, isym)
      if (any(abs(td - real(nint(td), dp)) > eps6)) cycle
      factor = 1.0_dp
      if (sym%spinors) then
        if (trev == 1) then
          us = matmul(uspin_t, conjg(sym%u_spin(:, :, isym)))
        else
          us = sym%u_spin(:, :, isym)
        end if
        if (all(abs(u0 - us) < tol_u)) then
          factor = 1.0_dp
        else if (all(abs(u0 + us) < tol_u)) then
          factor = -1.0_dp
        else
          cycle
        end if
      end if
      isym_out = isym
      tdiff = td
      return
    end do

    call set_error_fatal(error, &
                         'Error: read_ibz: the product of symmetry operations is not an operation of the crystal', comm)
  end subroutine search_symop

  !================================================!
  subroutine read_immn(sym, kmesh_info, kpt_latt, nn_of_ib, mmn_irr, num_bands, nb, seedname, &
                       error, comm)
    !================================================!
    !! Read the overlap matrices of the irreducible k-points from ".immn".
    !!
    !! pw2wannier90.x writes the blocks of each irreducible k-point in the
    !! order of the b-vectors of the ".nnkp" file, i.e. in the order of the
    !! global b-vector list of bvector_setup(), and heads each block with
    !!   iks  ikp  G(1:3)
    !! where ikp is the irreducible k-point equivalent to k_iks+b and G is
    !! the reciprocal lattice vector with  k_iks + b = (+/-) S ikp + G.  The
    !! header does not identify b uniquely (at high-symmetry k-points several
    !! b give the same ikp and G), so the block order has to be trusted; but
    !! every header is checked against the k-point mesh of this run, which
    !! catches a ".nnkp" file generated with a different neighbour set or a
    !! different symmetry convention.
    !================================================!

    implicit none

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: nb
    integer, intent(in) :: nn_of_ib(:, :)
    integer, intent(in) :: num_bands

    complex(kind=dp), intent(out) :: mmn_irr(:, :, :, :)
    real(kind=dp), intent(in) :: kpt_latt(:, :)

    character(len=*), intent(in) :: seedname

    ! local variables
    real(kind=dp) :: kb(3), kdiff(3), m_real, m_imag
    integer :: g(3), ib, idum(5), ierr, ik0, ikb, ikp, iks, isym, m, mmn_in, n, nb_tmp, nks_tmp, &
               nn, nn_tmp
    character(len=50) :: dummy
    character(len=200) :: msg

    open (newunit=mmn_in, file=trim(seedname)//'.immn', form='formatted', status='old', &
          action='read', iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error: Problem opening input file '//trim(seedname)//'.immn', comm)
      return
    end if

    read (mmn_in, '(a)', iostat=ierr) dummy
    if (ierr /= 0) goto 100
    read (mmn_in, *, iostat=ierr) nb_tmp, nks_tmp, nn_tmp
    if (ierr /= 0) goto 100

    if (nb_tmp /= num_bands) then
      close (mmn_in)
      call set_error_file(error, trim(seedname)//'.immn has not the right number of bands', comm)
      return
    end if
    if (nks_tmp /= sym%nks) then
      close (mmn_in)
      call set_error_file(error, trim(seedname)//'.immn has not the right number of k-points', comm)
      return
    end if
    if (nn_tmp /= nb) then
      close (mmn_in)
      call set_error_file(error, trim(seedname)//'.immn has not the right number of nearest neighbours', comm)
      return
    end if

    do iks = 1, sym%nks
      ik0 = sym%iks2ik(iks)
      do ib = 1, nb
        read (mmn_in, *, iostat=ierr) idum(:)
        if (ierr /= 0) goto 100

        ! the k-point k_iks + b of this block, as pw2wannier90.x must have
        ! found it: its irreducible equivalent and the lattice vector G
        nn = nn_of_ib(ib, ik0)
        ikb = kmesh_info%nnlist(ik0, nn)
        ikp = sym%equiv(ikb)
        isym = sym%equiv_sym(ikb)
        kb(:) = sym%irr_kpt(:, iks) + kpt_latt(:, ikb) - kpt_latt(:, ik0) &
                + real(kmesh_info%nncell(:, ik0, nn), dp)
        kdiff(:) = matmul(real(sym%s(:, :, isym), dp), sym%irr_kpt(:, ikp))
        if (sym%t_rev(isym) == 1) kdiff = -kdiff
        kdiff = kb - kdiff
        g = nint(kdiff)
        if (any(abs(kdiff - real(g, dp)) > tol_k)) then
          close (mmn_in)
          call set_error_fatal(error, 'Error: read_ibz: k+b is not related to an irreducible k-point '// &
                               'by a lattice vector (internal error)', comm)
          return
        end if
        if (idum(1) /= iks .or. idum(2) /= ikp .or. any(idum(3:5) /= g)) then
          close (mmn_in)
          ! (error messages are truncated at 128 characters)
          write (msg, '(a,i0,a,i0,a,5(1x,i0),a,5(1x,i0),a)') 'Error: '//trim(seedname)//'.immn k-point ', &
            iks, ' block ', ib, ': header', idum(:), ' but mesh gives', iks, ikp, g, '; .nnkp mismatch?'
          call set_error_file(error, trim(msg), comm)
          return
        end if

        do n = 1, num_bands
          do m = 1, num_bands
            read (mmn_in, *, iostat=ierr) m_real, m_imag
            if (ierr /= 0) goto 100
            mmn_irr(m, n, ib, iks) = cmplx(m_real, m_imag, kind=dp)
          end do
        end do
      end do
    end do
    close (mmn_in)

    return

100 close (mmn_in)
    call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.immn', comm)
    return
  end subroutine read_immn

  !================================================!
  subroutine read_iamn(sym, amn_irr, num_bands, num_wann, seedname, error, comm)
    !================================================!
    !! Read the projection matrices of the irreducible k-points from ".iamn"
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann

    complex(kind=dp), intent(out) :: amn_irr(:, :, :)

    character(len=*), intent(in) :: seedname

    ! local variables
    real(kind=dp) :: a_real, a_imag
    integer :: amn_in, icount, ierr, iks, m, n, nb_tmp, nks_tmp, np_tmp
    character(len=50) :: dummy

    open (newunit=amn_in, file=trim(seedname)//'.iamn', form='formatted', status='old', &
          action='read', iostat=ierr)
    if (ierr /= 0) then
      call set_error_file(error, 'Error: Problem opening input file '//trim(seedname)//'.iamn', comm)
      return
    end if

    read (amn_in, '(a)', iostat=ierr) dummy
    if (ierr /= 0) goto 100
    read (amn_in, *, iostat=ierr) nb_tmp, nks_tmp, np_tmp
    if (ierr /= 0) goto 100

    if (nb_tmp /= num_bands) then
      close (amn_in)
      call set_error_file(error, trim(seedname)//'.iamn has not the right number of bands', comm)
      return
    end if
    if (nks_tmp /= sym%nks) then
      close (amn_in)
      call set_error_file(error, trim(seedname)//'.iamn has not the right number of k-points', comm)
      return
    end if
    if (np_tmp /= num_wann) then
      close (amn_in)
      call set_error_file(error, trim(seedname)//'.iamn has not the right number of projections', comm)
      return
    end if

    amn_irr = cmplx_0
    do icount = 1, num_bands*num_wann*sym%nks
      read (amn_in, *, iostat=ierr) m, n, iks, a_real, a_imag
      if (ierr /= 0) goto 100
      if (m < 1 .or. m > num_bands .or. n < 1 .or. n > num_wann .or. iks < 1 .or. iks > sym%nks) then
        close (amn_in)
        call set_error_file(error, 'Error: Found a mismatch in '//trim(seedname)//'.iamn', comm)
        return
      end if
      amn_irr(m, n, iks) = cmplx(a_real, a_imag, kind=dp)
    end do
    close (amn_in)

    return

100 close (amn_in)
    call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.iamn', comm)
    return
  end subroutine read_iamn

  !================================================!
  subroutine projection_shifts(sym, pos, num_wann, rshift, error, comm)
    !================================================!
    !! rshift(:,iw,isym) is the lattice vector by which the symmetry operation
    !! isym displaces the centre of the projection iw away from the centre of
    !! the projections it is mapped onto
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann

    real(kind=dp), intent(in) :: pos(:, :)
    real(kind=dp), allocatable, intent(out) :: rshift(:, :, :)

    ! local variables
    real(kind=dp) :: r0(3)
    integer :: ierr, isym, iw, m
    logical :: found

    allocate (rshift(3, num_wann, sym%nsym), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating rshift in projection_shifts', comm)
      return
    end if

    do isym = 1, sym%nsym
      do iw = 1, num_wann
        found = .false.
        r0 = 0.0_dp
        do m = 1, num_wann
          if (abs(sym%rotmat(m, iw, isym)) < eps8) cycle
          if (.not. found) then
            r0(:) = pos(:, m)
            found = .true.
          else if (any(abs(pos(:, m) - r0(:)) > tol_k)) then
            call set_error_fatal(error, 'Error: read_ibz: a symmetry operation mixes projections '// &
                                 'sitting on different centres', comm)
            return
          end if
        end do
        if (.not. found) then
          call set_error_fatal(error, &
                               'Error: read_ibz: the rotation matrix of the projections has an empty column', comm)
          return
        end if
        rshift(:, iw, isym) = matmul(pos(:, iw), real(sym%s(:, :, isym), dp)) - sym%ft(:, isym) - r0(:)
        if (any(abs(rshift(:, iw, isym) - real(nint(rshift(:, iw, isym)), dp)) > tol_shift)) then
          call set_error_fatal(error, &
                               'Error: read_ibz: a rotated projection centre is not a lattice vector away '// &
                               'from its image', comm)
          return
        end if
        ! the shift is a lattice vector by construction; rounding it makes the
        ! phases below insensitive to the precision of the centres given in the
        ! input file
        rshift(:, iw, isym) = real(nint(rshift(:, iw, isym)), dp)
      end do
    end do
  end subroutine projection_shifts

  !================================================!
  subroutine symmetrise_amn(sym, amn_irr, rshift, num_bands, num_wann, error, comm)
    !================================================!
    !! Symmetrise the projections of each irreducible k-point over the little
    !! group of that k-point:
    !!   A = 1/N(h) sum_h <psi_m k| h |psi_l> <psi_l| h^-1 |g_n>
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_wann

    complex(kind=dp), intent(inout) :: amn_irr(:, :, :)
    real(kind=dp), intent(in) :: rshift(:, :, :)

    ! local variables
    complex(kind=dp), allocatable :: amn_sym(:, :, :), atmp(:, :)
    complex(kind=dp) :: phase
    real(kind=dp) :: kdiff(3), sk(3)
    integer :: ierr, iks, isym, iw, nh

    allocate (amn_sym(num_bands, num_wann, sym%nks), atmp(num_bands, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating amn_sym in symmetrise_amn', comm)
      return
    end if
    amn_sym = cmplx_0

    do iks = 1, sym%nks
      nh = 0
      do isym = 1, sym%nsym
        ! keep only the operations of the little group of this k-point
        sk = matmul(real(sym%s(:, :, isym), dp), sym%irr_kpt(:, iks))
        if (sym%t_rev(isym) == 1) sk = -sk
        kdiff = sym%irr_kpt(:, iks) - sk
        if (any(abs(kdiff - real(nint(kdiff), dp)) > tol_k)) cycle
        nh = nh + 1

        atmp = matmul(amn_irr(:, :, iks), sym%rotmat(:, :, isym))
        do iw = 1, num_wann
          phase = exp(-cmplx_i*twopi*dot_product(sym%irr_kpt(:, iks), rshift(:, iw, isym)))
          atmp(:, iw) = atmp(:, iw)*phase
        end do
        if (sym%t_rev(isym) == 1) atmp = conjg(atmp)
        amn_sym(:, :, iks) = amn_sym(:, :, iks) + matmul(sym%repmat(:, :, isym, iks), atmp)
      end do

      if (nh == 0) then
        call set_error_fatal(error, 'Error: read_ibz: empty little group for an irreducible k-point', comm)
        return
      end if
      amn_sym(:, :, iks) = amn_sym(:, :, iks)/real(nh, dp)
    end do

    amn_irr(:, :, :) = amn_sym(:, :, :)

    deallocate (amn_sym, atmp, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating amn_sym in symmetrise_amn', comm)
      return
    end if
  end subroutine symmetrise_amn

  !================================================!
  subroutine expand_amn(sym, amn_irr, au_matrix, rshift, num_bands, num_kpts, num_wann)
    !================================================!
    !! Rotate the symmetrised projections of the irreducible k-points onto the
    !! full k-point mesh
    !================================================!

    implicit none

    ! arguments
    type(sym_ibz_type), intent(in) :: sym

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann

    complex(kind=dp), intent(in) :: amn_irr(:, :, :)
    complex(kind=dp), intent(inout) :: au_matrix(:, :, :)
    real(kind=dp), intent(in) :: rshift(:, :, :)

    ! local variables
    complex(kind=dp) :: atmp(num_bands, num_wann)
    complex(kind=dp) :: phase
    integer :: ik, iks, isym, iw

    do ik = 1, num_kpts
      iks = sym%equiv(ik)
      isym = sym%equiv_sym(ik)
      atmp = matmul(amn_irr(:, :, iks), sym%rotmat(:, :, isym))
      do iw = 1, num_wann
        phase = exp(-cmplx_i*twopi*dot_product(sym%irr_kpt(:, iks), rshift(:, iw, isym)))
        atmp(:, iw) = atmp(:, iw)*phase
      end do
      if (sym%t_rev(isym) == 1) atmp = conjg(atmp)
      au_matrix(:, :, ik) = atmp(:, :)
    end do
  end subroutine expand_amn

  !================================================!
  subroutine expand_mmn(sym, kmesh_info, mmn_irr, m_matrix_local, bvec, bequiv, ib_of_nn, &
                        nn_of_ib, map_kpts, num_bands, num_kpts, nb, dist_k, rank, mmn_out, &
                        error, comm)
    !================================================!
    !! Rotate the overlap matrices of the irreducible k-points onto the full
    !! k-point mesh.  Only the k-points held by this MPI rank are stored, but
    !! the symmetry bookkeeping is done for all of them so that any error is
    !! raised on every rank.  If mmn_out is a non-zero unit (root rank only),
    !! every block is also computed and written to it in ".mmn" format.
    !!
    !! The Bloch state at a full-mesh k-point is defined as g(equiv_sym(k))
    !! applied to the state at its irreducible k-point.  For the bra (isym1)
    !! and for the ket at the full-mesh neighbour (isym3) this is a gauge
    !! choice and any operation of the coset would do; for the ket of the
    !! irreducible data, M(k_irr, b_irr) = <psi_k_irr| e^{-i b r} |psi_k_irr+b_irr>,
    !! isym2 must be the operation pw2wannier90.x used to build the state at
    !! k_irr+b_irr, i.e. the first matching operation in the order of the .isym
    !! file, which is what kpoint_map() stores.
    !================================================!

    implicit none

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(sym_ibz_type), intent(in) :: sym
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: bequiv(:, :), ib_of_nn(:, :), nn_of_ib(:, :), map_kpts(:)
    integer, intent(in) :: dist_k(:)
    integer, intent(in) :: mmn_out
    integer, intent(in) :: nb
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: rank

    complex(kind=dp), intent(in) :: mmn_irr(:, :, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :)
    real(kind=dp), intent(in) :: bvec(:, :)

    ! local variables
    complex(kind=dp) :: mtmp(num_bands, num_bands), rep(num_bands, num_bands)
    real(kind=dp) :: arg1, arg2, factor, tdiff(3)
    integer :: ib_full, ib_irr, ik0, ikb_full, ikb_irr, ikb_irr_eq, ik_irr
    integer :: ik, isym, isym1, isym2, isym3, list_inv(3), list_sym(3), m, n, nn, nn_irr

    do ik = 1, num_kpts
      ik_irr = sym%equiv(ik)
      isym1 = sym%equiv_sym(ik)
      ik0 = sym%iks2ik(ik_irr)

      do nn = 1, nb
        ! b-vector of this neighbour, in the global ordering, and the b-vector
        ! at the irreducible k-point which is mapped onto it by isym1
        ib_full = ib_of_nn(nn, ik)
        ib_irr = bequiv(ib_full, isym1)
        nn_irr = nn_of_ib(ib_irr, ik0)

        ikb_irr = kmesh_info%nnlist(ik0, nn_irr)
        ikb_full = kmesh_info%nnlist(ik, nn)
        isym2 = sym%equiv_sym(ikb_irr) ! fixed by the pw2wannier90 convention, see above
        isym3 = sym%equiv_sym(ikb_full)

        ! g(isym2)^-1 g(isym1)^-1 g(isym3)
        list_sym = (/isym2, isym1, isym3/)
        list_inv = (/-1, -1, 1/)
        call search_symop(sym, list_sym, list_inv, 3, isym, factor, tdiff, error, comm)
        if (allocated(error)) return

        if (dist_k(ik) /= rank .and. mmn_out == 0) cycle

        ikb_irr_eq = sym%equiv(ikb_irr)
        rep = sym%repmat(:, :, isym, ikb_irr_eq)
        if (sym%t_rev(isym2) == 1) rep = conjg(rep)
        mtmp = matmul(mmn_irr(:, :, ib_irr, ik_irr), rep)*factor
        if (sym%t_rev(isym1) == 1) mtmp = conjg(mtmp)

        ! exp(-i b_i . ft) comes from undoing the symmetry operation on
        ! exp(-i b_f . r); exp(-i (k_i+b_i) . tdiff) from the residual
        ! translation of the composed operation
        arg1 = -dot_product(bvec(:, ib_irr), sym%ft(:, isym1))
        arg2 = -dot_product(sym%irr_kpt(:, ikb_irr_eq), tdiff)
        if (sym%t_rev(isym2) == 1) arg2 = -arg2
        if (sym%t_rev(isym1) == 1) then
          arg1 = -arg1
          arg2 = -arg2
        end if
        if (sym%t_rev(isym) == 1) arg2 = -arg2
        mtmp = mtmp*exp(cmplx_i*twopi*(arg1 + arg2))

        if (mmn_out /= 0) then
          write (mmn_out, '(5i8)') ik, ikb_full, kmesh_info%nncell(:, ik, nn)
          do n = 1, num_bands
            do m = 1, num_bands
              write (mmn_out, '(2f18.12)') real(mtmp(m, n), dp), aimag(mtmp(m, n))
            end do
          end do
        end if

        if (dist_k(ik) == rank) m_matrix_local(:, :, nn, map_kpts(ik)) = mtmp(:, :)
      end do
    end do
  end subroutine expand_mmn

  !================================================!
  subroutine open_expanded_file(unit, filename, error, comm)
    !================================================!
    !! Open, on the root rank only, one of the full-BZ files written with
    !! write_ibz_expanded; unit is returned as zero on the other ranks.  The
    !! outcome is broadcast so that a failure is reported on every rank.
    !================================================!

    implicit none

    ! arguments
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(out) :: unit

    character(len=*), intent(in) :: filename

    ! local variables
    integer :: ierr

    unit = 0
    ierr = 0
    if (mpirank(comm) == 0) then
      open (newunit=unit, file=filename, form='formatted', status='replace', action='write', &
            iostat=ierr)
    end if
    call comms_bcast(ierr, 1, error, comm)
    if (allocated(error)) return
    if (ierr /= 0) then
      call set_error_file(error, 'Error: Problem opening output file '//trim(filename), comm)
      return
    end if
  end subroutine open_expanded_file

  !================================================!
  subroutine write_full_amn(au_matrix, num_bands, num_kpts, num_wann, print_output, seedname, &
                            stdout, error, comm)
    !================================================!
    !! Write the expanded projection matrices to ".amn" (root rank only)
    !================================================!

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout

    complex(kind=dp), intent(in) :: au_matrix(:, :, :)

    character(len=*), intent(in) :: seedname

    ! local variables
    integer :: amn_out, ik, m, n
    character(len=9) :: cdate, ctime

    call open_expanded_file(amn_out, trim(seedname)//'.amn', error, comm)
    if (allocated(error)) return
    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      'Writing the expanded projections to '//trim(seedname)//'.amn'
    if (amn_out == 0) return

    call io_date(cdate, ctime)
    write (amn_out, '(a)') 'Expanded from '//trim(seedname)//'.iamn by wannier90.x (read_ibz) on '// &
      cdate//' at '//ctime
    write (amn_out, '(3i12)') num_bands, num_kpts, num_wann
    do ik = 1, num_kpts
      do n = 1, num_wann
        do m = 1, num_bands
          write (amn_out, '(3i8,2f18.12)') m, n, ik, real(au_matrix(m, n, ik), dp), &
            aimag(au_matrix(m, n, ik))
        end do
      end do
    end do
    close (amn_out)
  end subroutine write_full_amn

  !================================================!
  subroutine write_full_eig(eigval, num_bands, num_kpts, print_output, seedname, stdout, error, comm)
    !================================================!
    !! Write the expanded eigenvalues to ".eig" (root rank only)
    !================================================!

    implicit none

    ! arguments
    type(print_output_type), intent(in) :: print_output
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: eigval(:, :)

    character(len=*), intent(in) :: seedname

    ! local variables
    integer :: eig_out, ik, n

    call open_expanded_file(eig_out, trim(seedname)//'.eig', error, comm)
    if (allocated(error)) return
    if (print_output%iprint > 0) write (stdout, '(1x,a)') &
      'Writing the expanded eigenvalues to '//trim(seedname)//'.eig'
    if (eig_out == 0) return

    do ik = 1, num_kpts
      do n = 1, num_bands
        write (eig_out, '(2i8,f18.12)') n, ik, eigval(n, ik)
      end do
    end do
    close (eig_out)
  end subroutine write_full_eig

end module w90_sym_ibz
