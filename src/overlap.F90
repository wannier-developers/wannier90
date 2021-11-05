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
!  w90_overlap: setup overlap and projection matrices        !
!                                                            !
!------------------------------------------------------------!

module w90_overlap

  !! This module reads in the overlap (Mmn) and Projections (Amn)
  !! and performs simple operations on them.

  use w90_constants, only: dp, cmplx_0, cmplx_1
  use w90_comms

  implicit none

  private

  public :: overlap_allocate
  public :: overlap_dealloc
  public :: overlap_project
  public :: overlap_project_gamma
  public :: overlap_read

contains

  !==================================================================!

  subroutine overlap_allocate(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, &
                              m_matrix_orig_local, u_matrix, u_matrix_opt, nntot, num_bands, &
                              num_kpts, num_wann, timing_level, seedname, stdout, comm)
    !==================================================================!
    !! Allocate memory to read Mmn and Amn from files
    !! This must be called before calling overlap_read
    !                                                                  !
    !==================================================================!

    use w90_io, only: io_error, io_stopwatch

    ! arguments
    integer, intent(in) :: nntot
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level

    complex(kind=dp), allocatable :: a_matrix(:, :, :)
    complex(kind=dp), allocatable :: m_matrix(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_local(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), allocatable :: u_matrix(:, :, :)
    complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)

    type(w90comm_type), intent(in) :: comm

    character(len=50), intent(in) :: seedname

    ! local variables
    integer, allocatable :: counts(:)
    integer, allocatable :: displs(:)
    integer :: ierr
    integer :: num_nodes, my_node_id
    logical :: disentanglement
    logical :: on_root = .false.

    disentanglement = (num_bands > num_wann)

    num_nodes = mpisize(comm)
    my_node_id = mpirank(comm)

    if (my_node_id == 0) on_root = .true.
    allocate (counts(0:num_nodes - 1))
    allocate (displs(0:num_nodes - 1))

    if (timing_level > 0) call io_stopwatch('overlap: allocate', 1, stdout, seedname)

    call comms_array_split(num_kpts, counts, displs, comm)

    allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating u_matrix in overlap_read', stdout, seedname)
    u_matrix = cmplx_0

    if (disentanglement) then
      if (on_root) then
        allocate (m_matrix_orig(num_bands, num_bands, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating m_matrix_orig in overlap_read', stdout, seedname)
        allocate (m_matrix(1, 1, 1, 1))
      else
        allocate (m_matrix_orig(1, 1, 1, 1))
        allocate (m_matrix(1, 1, 1, 1))
      endif

      allocate (m_matrix_orig_local(num_bands, num_bands, nntot, counts(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating m_matrix_orig_local in overlap_read', stdout, seedname)
      allocate (a_matrix(num_bands, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating a_matrix in overlap_read', stdout, seedname)
      allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating u_matrix_opt in overlap_read', stdout, seedname)

      allocate (m_matrix_local(1, 1, 1, 1))

    else
      if (on_root) then
        allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating m_matrix in overlap_read', stdout, seedname)
        m_matrix = cmplx_0
        allocate (m_matrix_orig(1, 1, 1, 1))
      else
        allocate (m_matrix(1, 1, 1, 1))
        allocate (m_matrix_orig(1, 1, 1, 1))
      endif

      allocate (m_matrix_local(num_wann, num_wann, nntot, counts(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating m_matrix_local in overlap_read', stdout, seedname)
      m_matrix_local = cmplx_0

      allocate (m_matrix_orig_local(1, 1, 1, 1))
      allocate (a_matrix(1, 1, 1))
      allocate (u_matrix_opt(1, 1, 1))

    endif

    if (timing_level > 0) call io_stopwatch('overlap: allocate', 2, stdout, seedname)

  end subroutine overlap_allocate

  !==================================================================!
  subroutine overlap_read(kmesh_info, select_projection, sitesym, a_matrix, m_matrix, &
                          m_matrix_local, m_matrix_orig, m_matrix_orig_local, u_matrix, &
                          u_matrix_opt, num_bands, num_kpts, num_proj, num_wann, timing_level, &
                          cp_pp, gamma_only, lsitesymmetry, use_bloch_phases, seedname, stdout, &
                          comm)
    !==================================================================!
    !! Read the Mmn and Amn from files
    !! Note: one needs to call overlap_allocate first!
    !                                                                  !
    !==================================================================!

    use w90_io, only: io_file_unit, io_error, io_stopwatch
    use w90_types, only: kmesh_info_type
    use w90_wannier90_types, only: select_projection_type, sitesym_type

    implicit none

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(select_projection_type), intent(in) :: select_projection
    type(sitesym_type), intent(in) :: sitesym
    type(w90comm_type), intent(in) :: comm

    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_proj
    integer, intent(in) :: num_wann
    integer, intent(in) :: timing_level
    integer, intent(in) :: stdout

    complex(kind=dp), intent(inout) :: a_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: u_matrix_opt(:, :, :)

    logical, intent(in) :: gamma_only
    logical, intent(in) :: lsitesymmetry
    logical, intent(in) :: cp_pp, use_bloch_phases

    character(len=50), intent(in)  :: seedname

    ! local variables
    integer :: mmn_in, amn_in, num_mmn, num_amn
    integer :: nb_tmp, nkp_tmp, nntot_tmp, np_tmp, ierr
    integer :: nkp, nkp2, inn, nn, n, m
    integer :: nnl, nnm, nnn, ncount
    logical :: nn_found
    real(kind=dp) :: m_real, m_imag, a_real, a_imag
    complex(kind=dp), allocatable :: mmn_tmp(:, :)
    character(len=50) :: dummy

    logical :: disentanglement
    integer :: num_nodes, my_node_id
    logical :: on_root = .false.
    integer, allocatable :: counts(:)
    integer, allocatable :: displs(:)

    disentanglement = (num_bands > num_wann)

    num_nodes = mpisize(comm)
    my_node_id = mpirank(comm)

    if (my_node_id == 0) on_root = .true.
    allocate (counts(0:num_nodes - 1))
    allocate (displs(0:num_nodes - 1))

    if (timing_level > 0) call io_stopwatch('overlap: read', 1, stdout, seedname)

    call comms_array_split(num_kpts, counts, displs, comm)

    if (disentanglement) then
      if (on_root) then
        m_matrix_orig = cmplx_0
      endif
      m_matrix_orig_local = cmplx_0
      a_matrix = cmplx_0
      u_matrix_opt = cmplx_0
    endif

    if (on_root) then

      ! Read M_matrix_orig from file
      mmn_in = io_file_unit()
      open (unit=mmn_in, file=trim(seedname)//'.mmn', &
            form='formatted', status='old', action='read', err=101)

      if (on_root) write (stdout, '(/a)', advance='no') ' Reading overlaps from '//trim(seedname)//'.mmn    : '

      ! Read the comment line
      read (mmn_in, '(a)', err=103, end=103) dummy
      if (on_root) write (stdout, '(a)') trim(dummy)

      ! Read the number of bands, k-points and nearest neighbours
      read (mmn_in, *, err=103, end=103) nb_tmp, nkp_tmp, nntot_tmp

      ! Checks
      if (nb_tmp .ne. num_bands) &
        call io_error(trim(seedname)//'.mmn has not the right number of bands', stdout, seedname)
      if (nkp_tmp .ne. num_kpts) &
        call io_error(trim(seedname)//'.mmn has not the right number of k-points', stdout, seedname)
      if (nntot_tmp .ne. kmesh_info%nntot) &
        call io_error(trim(seedname)//'.mmn has not the right number of nearest neighbours', stdout, seedname)

      ! Read the overlaps
      num_mmn = num_kpts*kmesh_info%nntot
      allocate (mmn_tmp(num_bands, num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating mmn_tmp in overlap_read', stdout, seedname)
      do ncount = 1, num_mmn
        read (mmn_in, *, err=103, end=103) nkp, nkp2, nnl, nnm, nnn
        do n = 1, num_bands
          do m = 1, num_bands
            read (mmn_in, *, err=103, end=103) m_real, m_imag
            mmn_tmp(m, n) = cmplx(m_real, m_imag, kind=dp)
          enddo
        enddo
        nn = 0
        nn_found = .false.
        do inn = 1, kmesh_info%nntot
          if ((nkp2 .eq. kmesh_info%nnlist(nkp, inn)) .and. &
              (nnl .eq. kmesh_info%nncell(1, nkp, inn)) .and. &
              (nnm .eq. kmesh_info%nncell(2, nkp, inn)) .and. &
              (nnn .eq. kmesh_info%nncell(3, nkp, inn))) then
            if (.not. nn_found) then
              nn_found = .true.
              nn = inn
            else
              call io_error('Error reading '//trim(seedname)// &
                            '.mmn. More than one matching nearest neighbour found', stdout, seedname)
            endif
          endif
        end do
        if (nn .eq. 0) then
          if (on_root) write (stdout, '(/a,i8,2i5,i4,2x,3i3)') &
            ' Error reading '//trim(seedname)//'.mmn:', ncount, nkp, nkp2, nn, nnl, nnm, nnn
          call io_error('Neighbour not found', stdout, seedname)
        end if
        if (disentanglement) then
          m_matrix_orig(:, :, nn, nkp) = mmn_tmp(:, :)
        else
          ! disentanglement=.false. means numbands=numwann, so no the dimensions are the same
          m_matrix(:, :, nn, nkp) = mmn_tmp(:, :)
        end if
      end do
      deallocate (mmn_tmp, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating mmn_tmp in overlap_read', stdout, seedname)
      close (mmn_in)
    endif

    if (disentanglement) then
      call comms_scatterv(m_matrix_orig_local, num_bands*num_bands*kmesh_info%nntot*counts(my_node_id), &
                          m_matrix_orig, num_bands*num_bands*kmesh_info%nntot*counts, &
                          num_bands*num_bands*kmesh_info%nntot*displs, stdout, seedname, comm)
    else
      call comms_scatterv(m_matrix_local, num_wann*num_wann*kmesh_info%nntot*counts(my_node_id), &
                          m_matrix, num_wann*num_wann*kmesh_info%nntot*counts, &
                          num_wann*num_wann*kmesh_info%nntot*displs, stdout, seedname, comm)
    endif

    if (.not. use_bloch_phases) then
      if (on_root) then

        ! Read A_matrix from file wannier.amn
        amn_in = io_file_unit()
        open (unit=amn_in, file=trim(seedname)//'.amn', form='formatted', status='old', err=102)

        if (on_root) write (stdout, '(/a)', advance='no') ' Reading projections from '//trim(seedname)//'.amn : '

        ! Read the comment line
        read (amn_in, '(a)', err=104, end=104) dummy
        if (on_root) write (stdout, '(a)') trim(dummy)

        ! Read the number of bands, k-points and wannier functions
        read (amn_in, *, err=104, end=104) nb_tmp, nkp_tmp, np_tmp

        ! Checks
        if (nb_tmp .ne. num_bands) &
          call io_error(trim(seedname)//'.amn has not the right number of bands', stdout, seedname)
        if (nkp_tmp .ne. num_kpts) &
          call io_error(trim(seedname)//'.amn has not the right number of k-points', stdout, seedname)
        if (np_tmp .ne. num_proj) &
          call io_error(trim(seedname)//'.amn has not the right number of projections', stdout, seedname)

        if (num_proj > num_wann .and. .not. select_projection%lselproj) &
          call io_error(trim(seedname)//'.amn has too many projections to be used without selecting a subset', stdout, seedname)

        ! Read the projections
        num_amn = num_bands*num_proj*num_kpts
        if (disentanglement) then
          do ncount = 1, num_amn
            read (amn_in, *, err=104, end=104) m, n, nkp, a_real, a_imag
            if (select_projection%proj2wann_map(n) < 0) cycle
            a_matrix(m, select_projection%proj2wann_map(n), nkp) = cmplx(a_real, a_imag, kind=dp)
          end do
        else
          do ncount = 1, num_amn
            read (amn_in, *, err=104, end=104) m, n, nkp, a_real, a_imag
            if (select_projection%proj2wann_map(n) < 0) cycle
            u_matrix(m, select_projection%proj2wann_map(n), nkp) = cmplx(a_real, a_imag, kind=dp)
          end do
        end if
        close (amn_in)
      endif

      if (disentanglement) then
        call comms_bcast(a_matrix(1, 1, 1), num_bands*num_wann*num_kpts, stdout, seedname, comm)
      else
        call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, stdout, seedname, comm)
      endif

    else

      do n = 1, num_kpts
        do m = 1, num_wann
          u_matrix(m, m, n) = cmplx_1
        end do
      end do

    end if

    ! If post-processing a Car-Parinello calculation (gamma only)
    ! then rotate M and A to the basis of Kohn-Sham eigenstates
    if (cp_pp) call overlap_rotate(a_matrix, m_matrix_orig, kmesh_info%nntot, num_bands, &
                                   timing_level, seedname, stdout)

    ! Check Mmn(k,b) is symmetric in m and n for gamma_only case
!~      if (gamma_only) call overlap_check_m_symmetry()

    ! If we don't need to disentangle we can now convert from A to U
    ! And rotate M accordingly
![ysl-b]
!       if(.not.disentanglement .and. (.not.cp_pp) .and. (.not. use_bloch_phases )) &
!            call overlap_project
!~       if((.not.cp_pp) .and. (.not. use_bloch_phases )) then
!~         if (.not.disentanglement) then
!~            if ( .not. gamma_only ) then
!~               call overlap_project
!~            else
!~               call overlap_project_gamma()
!~            endif
!~         else
!~            if (gamma_only) call overlap_symmetrize()
!~         endif
!~       endif
!
!~[aam]
    if ((.not. disentanglement) .and. (.not. cp_pp) .and. (.not. use_bloch_phases)) then
      if (.not. gamma_only) then
        call overlap_project(sitesym, m_matrix, m_matrix_local, u_matrix, kmesh_info%nnlist, &
                             kmesh_info%nntot, num_bands, num_kpts, num_wann, &
                             timing_level, lsitesymmetry, seedname, stdout, comm)
      else
        call overlap_project_gamma(m_matrix, u_matrix, kmesh_info%nntot, num_wann, &
                                   timing_level, seedname, stdout)
      endif
    endif
!~[aam]

    !~      if( gamma_only .and. use_bloch_phases ) then
    !~        write(stdout,'(1x,"+",76("-"),"+")')
    !~        write(stdout,'(3x,a)') 'WARNING: gamma_only and use_bloch_phases                 '
    !~        write(stdout,'(3x,a)') '         M must be calculated from *real* Bloch functions'
    !~        write(stdout,'(1x,"+",76("-"),"+")')
    !~      end if
![ysl-e]

    if (timing_level > 0) call io_stopwatch('overlap: read', 2, stdout, seedname)

    return
101 call io_error('Error: Problem opening input file '//trim(seedname)//'.mmn', stdout, seedname)
102 call io_error('Error: Problem opening input file '//trim(seedname)//'.amn', stdout, seedname)
103 call io_error('Error: Problem reading input file '//trim(seedname)//'.mmn', stdout, seedname)
104 call io_error('Error: Problem reading input file '//trim(seedname)//'.amn', stdout, seedname)

  end subroutine overlap_read

  !==================================================================!
  subroutine overlap_rotate(a_matrix, m_matrix_orig, nntot, num_bands, timing_level, seedname, &
                            stdout)
    !==================================================================!
    !
    !! Only used when interfaced to the CP code
    !! Not sure why this is done here and not in CP
    !                                                                  !
    !==================================================================!

    use w90_io, only: io_file_unit, io_error, io_stopwatch

    implicit none

    integer, intent(in) :: nntot
    integer, intent(in) :: stdout
    integer, intent(in) :: num_bands
    integer, intent(in) :: timing_level

    complex(kind=dp), intent(inout) :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), intent(inout) :: a_matrix(:, :, :)

    character(len=50), intent(in)  :: seedname

    integer       :: lam_unit, info, inn, i, j
    real(kind=DP) :: lambda(num_bands, num_bands)
    real(kind=DP) :: AP(num_bands*(num_bands + 1)/2)
    real(kind=DP) :: eig(num_bands), work(3*num_bands)

    if (timing_level > 1) call io_stopwatch('overlap: rotate', 1, stdout, seedname)

    lam_unit = io_file_unit()
    open (unit=lam_unit, file='lambda.dat', &
          form='unformatted', status='old', action='read')
!~    write(stdout,*) ' Reading lambda.dat...'
    read (lam_unit) lambda
!~    write(stdout,*) ' done'
    close (lam_unit)

    do j = 1, num_bands
      do i = 1, j
        AP(i + (j - 1)*j/2) = 0.5_dp*(lambda(i, j) + lambda(j, i))
      end do
    end do

    CALL DSPEV('V', 'U', num_bands, AP, eig, lambda, num_bands, work, info)
    if (info .ne. 0) &
      call io_error('Diagonalization of lambda in overlap_rotate failed', stdout, seedname)

    ! For debugging
!~    write(stdout,*) 'EIGENVALUES - CHECK WITH CP OUTPUT'
!~    do i=1,num_bands
!~       write(stdout,*) 13.6058*eig(i)
!~    end do

    ! Rotate M_mn
    do inn = 1, nntot
      m_matrix_orig(:, :, inn, 1) = &
        matmul(transpose(lambda), matmul(m_matrix_orig(:, :, inn, 1), lambda))
    end do

    ! Rotate A_mn
    a_matrix(:, :, 1) = matmul(transpose(lambda), a_matrix(:, :, 1))

    ! For debugging
!~    ! Write rotated A and M
!~    do i=1,num_bands
!~       do j=1,num_wann
!~          write(12,'(2i5,a,2f18.12)') i,j,'   1',a_matrix(i,j,1)
!~       enddo
!~    enddo
!~    do inn=1,nntot
!~       do i=1,num_bands
!~          do j=1,num_wann
!~             write(11,'(2i5,2a)') i,j,'    1','    1'
!~             write(11,'(2f18.12)') m_matrix_orig(i,j,inn,1)
!~          enddo
!~       enddo
!~    enddo
!~    stop

    if (timing_level > 1) call io_stopwatch('overlap: rotate', 2, stdout, seedname)

    return

  end subroutine overlap_rotate

  !==================================================================!
  subroutine overlap_dealloc(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, &
                             m_matrix_orig_local, u_matrix, u_matrix_opt, seedname, stdout, comm)
    !==================================================================!
    !                                                                  !
    !! Dellocate memory
    !                                                                  !
    !==================================================================!

    use w90_io, only: io_error

    implicit none

    ! arguments
    complex(kind=dp), allocatable, intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: a_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix_local(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix_orig_local(:, :, :, :)
    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname
    type(w90comm_type), intent(in) :: comm

    ! local variables
    integer :: ierr
    logical :: on_root = .false.

    if (mpirank(comm) == 0) on_root = .true.

    if (allocated(u_matrix_opt)) then
      deallocate (u_matrix_opt, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating u_matrix_opt in overlap_dealloc', stdout, seedname)
    end if
    if (allocated(a_matrix)) then
      deallocate (a_matrix, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating a_matrix in overlap_dealloc', stdout, seedname)
    end if
!    if (on_root) then
    if (allocated(m_matrix_orig)) then
      deallocate (m_matrix_orig, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_orig in overlap_dealloc', stdout, seedname)
    endif
!    endif
    if (allocated(m_matrix_orig_local)) then
      deallocate (m_matrix_orig_local, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_orig_local in overlap_dealloc', stdout, seedname)
    endif
    if (on_root) then
      if (allocated(m_matrix)) then
        deallocate (m_matrix, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating m_matrix in overlap_dealloc', stdout, seedname)
      endif
    endif
    if (allocated(m_matrix_local)) then
      deallocate (m_matrix_local, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_local in overlap_dealloc', stdout, seedname)
    endif
    if (allocated(u_matrix)) then
      deallocate (u_matrix, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating u_matrix in overlap_dealloc', stdout, seedname)
    endif

    return

  end subroutine overlap_dealloc

  !==================================================================!
  subroutine overlap_project(sitesym, m_matrix, m_matrix_local, u_matrix, nnlist, nntot, &
                             num_bands, num_kpts, num_wann, timing_level, lsitesymmetry, &
                             seedname, stdout, comm)
    !==================================================================!
    !!  Construct initial guess from the projection via a Lowdin transformation
    !!  See section 3 of the CPC 2008
    !!  Note that in this subroutine num_wann = num_bands
    !!  since, if we are here, then disentanglement = FALSE
    !                                                                  !
    !==================================================================!
    use w90_constants
    use w90_io, only: io_error, io_stopwatch
    use w90_utility, only: utility_zgemm
    use w90_sitesym, only: sitesym_symmetrize_u_matrix
    use w90_wannier90_types, only: sitesym_type

    implicit none

    type(sitesym_type), intent(in) :: sitesym
    type(w90comm_type), intent(in) :: comm

    integer, intent(in) :: nnlist(:, :)
    integer, intent(in) :: nntot
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level

    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :)

    logical, intent(in) :: lsitesymmetry

    character(len=50), intent(in)  :: seedname

    ! local variables
    integer :: i, j, m, nkp, info, ierr, nn, nkp2
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp)                 :: rwork(5*num_bands)
    complex(kind=dp)              :: ctmp2
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cvdag(:, :)

    ! pllel setup
    integer, allocatable :: counts(:)
    integer, allocatable :: displs(:)
    integer :: num_nodes, my_node_id
    logical :: on_root = .false.

    num_nodes = mpisize(comm)
    my_node_id = mpirank(comm)
    if (my_node_id == 0) on_root = .true.
    allocate (counts(0:num_nodes - 1))
    allocate (displs(0:num_nodes - 1))

    if (timing_level > 1) call io_stopwatch('overlap: project', 1, stdout, seedname)

    call comms_array_split(num_kpts, counts, displs, comm)

    allocate (svals(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating svals in overlap_project', stdout, seedname)
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cz in overlap_project', stdout, seedname)
    allocate (cvdag(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cvdag in overlap_project', stdout, seedname)
    allocate (cwork(4*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwork in overlap_project', stdout, seedname)

    ! Calculate the transformation matrix CU = CS^(-1/2).CA,
    ! where CS = CA.CA^\dagger.

    do nkp = 1, num_kpts
      !
      ! SINGULAR VALUE DECOMPOSITION
      !
      call ZGESVD('A', 'A', num_bands, num_bands, u_matrix(1, 1, nkp), &
                  num_bands, svals, cz, num_bands, cvdag, num_bands, cwork, &
                  4*num_bands, rwork, info)
      if (info .ne. 0) then
        write (stdout, *) ' ERROR: IN ZGESVD IN overlap_project'
        write (stdout, *) ' K-POINT NKP=', nkp, ' INFO=', info
        if (info .lt. 0) then
          write (stdout, *) ' THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
        endif
        call io_error('Error in ZGESVD in overlap_project', stdout, seedname)
      endif

!       u_matrix(:,:,nkp)=matmul(cz,cvdag)
      call utility_zgemm(u_matrix(:, :, nkp), cz, 'N', cvdag, 'N', num_wann)

      !
      ! CHECK UNITARITY
      !
      do i = 1, num_bands
        do j = 1, num_bands
          ctmp2 = cmplx_0
          do m = 1, num_bands
            ctmp2 = ctmp2 + u_matrix(m, j, nkp)*conjg(u_matrix(m, i, nkp))
          enddo
          if ((i .eq. j) .and. (abs(ctmp2 - cmplx_1) .gt. eps5)) then
            write (stdout, *) ' ERROR: unitarity of initial U'
            write (stdout, '(1x,a,i2)') 'nkp= ', nkp
            write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
            write (stdout, '(1x,a,f12.6,1x,f12.6)') &
              '[u_matrix.transpose(u_matrix)]_ij= ', &
              real(ctmp2, dp), aimag(ctmp2)
            call io_error('Error in unitarity of initial U in overlap_project', stdout, seedname)
          endif
          if ((i .ne. j) .and. (abs(ctmp2) .gt. eps5)) then
            write (stdout, *) ' ERROR: unitarity of initial U'
            write (stdout, '(1x,a,i2)') 'nkp= ', nkp
            write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
            write (stdout, '(1x,a,f12.6,1x,f12.6)') &
              '[u_matrix.transpose(u_matrix)]_ij= ', &
              real(ctmp2, dp), aimag(ctmp2)
            call io_error('Error in unitarity of initial U in overlap_project', stdout, seedname)
          endif
        enddo
      enddo
    enddo
    ! NKP

    if (lsitesymmetry) call sitesym_symmetrize_u_matrix(sitesym, u_matrix, num_bands, num_wann, &
                                                        num_kpts, num_wann, seedname, stdout) !RS: update U(Rk)

    ! so now we have the U's that rotate the wavefunctions at each k-point.
    ! the matrix elements M_ij have also to be updated
    do nkp = 1, counts(my_node_id)
      do nn = 1, nntot
        nkp2 = nnlist(nkp + displs(my_node_id), nn)
        ! cvdag = U^{dagger} . M   (use as workspace)
        call utility_zgemm(cvdag, u_matrix(:, :, nkp + displs(my_node_id)), 'C', &
                           m_matrix_local(:, :, nn, nkp), 'N', num_wann)
        ! cz = cvdag . U
        call utility_zgemm(cz, cvdag, 'N', u_matrix(:, :, nkp2), 'N', num_wann)
        m_matrix_local(:, :, nn, nkp) = cz(:, :)
      end do
    end do
    call comms_gatherv(m_matrix_local, num_wann*num_wann*nntot*counts(my_node_id), &
                       m_matrix, num_wann*num_wann*nntot*counts, num_wann*num_wann*nntot*displs, &
                       stdout, seedname, comm)

    deallocate (cwork, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cwork in overlap_project', stdout, seedname)
    deallocate (cvdag, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cvdag in overlap_project', stdout, seedname)
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cz in overlap_project', stdout, seedname)
    deallocate (svals, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating svals in overlap_project', stdout, seedname)

    if (timing_level > 1) call io_stopwatch('overlap: project', 2, stdout, seedname)

    return

  end subroutine overlap_project

![ysl-b]
  !==================================================================!
  subroutine overlap_project_gamma(m_matrix, u_matrix, nntot, num_wann, timing_level, seedname, &
                                   stdout)
    !==================================================================!
    !!  Construct initial guess from the projection via a Lowdin transformation
    !!  See section 3 of the CPC 2008
    !!  Note that in this subroutine num_wann = num_bands
    !!  since, if we are here, then disentanglement = FALSE
    !!  Gamma specific version
    !                                                                  !
    !==================================================================!
    use w90_constants
    use w90_io, only: io_error, io_stopwatch
    use w90_utility, only: utility_zgemm

    implicit none

    integer, intent(in) :: nntot
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level
    integer, intent(in) :: num_wann
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    character(len=50), intent(in)  :: seedname

    ! internal variables
    integer :: i, j, m, info, ierr, nn
    real(kind=dp)                 :: rtmp2
    real(kind=dp), allocatable :: u_matrix_r(:, :)
!~    real(kind=dp),    allocatable :: u_cmp(:)
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp), allocatable :: work(:)
    real(kind=dp), allocatable :: rz(:, :)
    real(kind=dp), allocatable :: rv(:, :)
!~    complex(kind=dp), allocatable :: ph(:)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cvdag(:, :)

!~ real(kind=dp),    allocatable :: u_cmp(:)
!~ integer :: n,mdev, ndev, nndev,p(1)
!~ real(kind=dp)                 :: dev, dev_tmp

    if (timing_level > 1) call io_stopwatch('overlap: project_gamma', 1, stdout, seedname)

!~    allocate(ph_g(num_wann),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating ph_g in overlap_project_gamma')
    ! internal variables
    allocate (u_matrix_r(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating u_matrix_r in overlap_project_gamma', stdout, seedname)
!~    allocate(u_cmp(num_wann),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating u_cmp in overlap_project_gamma')
    allocate (svals(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating svals in overlap_project_gamma', stdout, seedname)
    allocate (work(5*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating work in overlap_project_gamma', stdout, seedname)
    allocate (rz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rz in overlap_project_gamma', stdout, seedname)
    allocate (rv(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rv in overlap_project_gamma', stdout, seedname)
    allocate (cz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cz in overlap_project_gamma', stdout, seedname)
    allocate (cvdag(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cvdag in overlap_project_gamma', stdout, seedname)

    !
!~    ! If a wavefunction is real except for a phase factor e^(i*phi_m) = ph_g(m)
!~    ! U_mn = ph_g(m)^(-1)*<m_R|g_n>   (m_R implies a real wave function)
!~    ! At m_th row, find five elements with largest complex component
!~    ! -> ph_g(m) is calculated from those elements
!~    ! U_mn (new) = ph_g(m) * U_mn (old)
!~    !
!~    ph_g=cmplx_1
!~    do m=1,num_wann
!~       u_cmp(:)=abs(u_matrix(m,:,1))
!~       p=maxloc(u_cmp)
!~       ph_g(m)=conjg(u_matrix(m,p(1),1)/abs(u_matrix(m,p(1),1)))
!~       u_matrix_r(m,:)=real(ph_g(m)*u_matrix(m,:,1),dp)
!~    end do

    u_matrix_r(:, :) = real(u_matrix(:, :, 1), dp)

!~    ! M_mn (new) = ph(m) * M_mn (old) * conjg(ph(n))
!~
!~    do nn=1,nntot
!~       do n=1,num_wann
!~          do m=1,num_wann
!~             m_matrix(m,n,nn,1)=ph_g(m)*m_matrix(m,n,nn,1)*conjg(ph_g(n))
!~          end do
!~       end do
!~    end do
!~    !
!~    ! check whether M_mn is now symmetric
!~    !
!~    dev = 0.0_dp
!~    do nn=1,nntot
!~       do n=1,num_wann
!~          do m=1,n
!~             dev_tmp=abs(m_matrix(m,n,nn,1)-m_matrix(n,m,nn,1))
!~             if ( dev_tmp .gt. dev ) then
!~                dev = dev_tmp
!~                mdev  = m ; ndev  = n ;  nndev  = nn
!~             end if
!~          end do
!~       end do
!~    end do
!~    if ( dev .gt. eps ) then
!~       write(stdout,'(1x,"+",76("-"),"+")')
!~       write(stdout,'(3x,a)') 'WARNING: M is not strictly symmetric in overlap_project_gamma'
!~       write(stdout,'(3x,a,f12.8)') &
!~            'Largest deviation |M_mn-M_nm| at k : ',dev
!~       write(stdout,'(3(a5,i4))') &
!~            '   m=',mdev,',  n=',ndev,',  k=',nndev
!~       write(stdout,'(1x,"+",76("-"),"+")')
!~    end if
    !
    ! Calculate the transformation matrix RU = RS^(-1/2).RA,
    ! where RS = RA.RA^\dagger.
    !
    ! SINGULAR VALUE DECOMPOSITION
    !
    call DGESVD('A', 'A', num_wann, num_wann, u_matrix_r, num_wann, &
                svals, rz, num_wann, rv, num_wann, work, 5*num_wann, info)
    if (info .ne. 0) then
      write (stdout, *) ' ERROR: IN DGESVD IN overlap_project_gamma'
      if (info .lt. 0) then
        write (stdout, *) 'THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
      endif
      call io_error('overlap_project_gamma: problem in DGESVD 1', stdout, seedname)
    endif

    call dgemm('N', 'N', num_wann, num_wann, num_wann, 1.0_dp, &
               rz, num_wann, rv, num_wann, 0.0_dp, u_matrix_r, num_wann)
    !
    ! CHECK UNITARITY
    !
    do i = 1, num_wann
      do j = 1, num_wann
        rtmp2 = 0.0_dp
        do m = 1, num_wann
          rtmp2 = rtmp2 + u_matrix_r(m, j)*u_matrix_r(m, i)
        enddo
        if ((i .eq. j) .and. (abs(rtmp2 - 1.0_dp) .gt. eps5)) then
          write (stdout, *) ' ERROR: unitarity of initial U'
          write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
          write (stdout, '(1x,a,f12.6)') &
            '[u_matrix.transpose(u_matrix)]_ij= ', &
            rtmp2
          call io_error('Error in unitarity of initial U in overlap_project_gamma', stdout, seedname)
        endif
        if ((i .ne. j) .and. (abs(rtmp2) .gt. eps5)) then
          write (stdout, *) ' ERROR: unitarity of initial U'
          write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
          write (stdout, '(1x,a,f12.6,1x,f12.6)') &
            '[u_matrix.transpose(u_matrix)]_ij= ', &
            rtmp2
          call io_error('Error in unitarity of initial U in overlap_project_gamma', stdout, seedname)
        endif
      enddo
    enddo

    u_matrix(:, :, 1) = cmplx(u_matrix_r(:, :), 0.0_dp, dp)

    ! so now we have the U's that rotate the wavefunctions at each k-point.
    ! the matrix elements M_ij have also to be updated

    do nn = 1, nntot
      ! cvdag = U^{dagger} . M   (use as workspace)
      call utility_zgemm(cvdag, u_matrix(:, :, 1), 'C', m_matrix(:, :, nn, 1), 'N', num_wann)
      ! cz = cvdag . U
      call utility_zgemm(cz, cvdag, 'N', u_matrix(:, :, 1), 'N', num_wann)
      m_matrix(:, :, nn, 1) = cz(:, :)
    end do

    deallocate (cvdag, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cvdag in overlap_project_gamma', stdout, seedname)
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cz in overlap_project_gamma', stdout, seedname)
    deallocate (rv, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rv in overlap_project_gamma', stdout, seedname)
    deallocate (rz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rz in overlap_project_gamma', stdout, seedname)
    deallocate (work, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating work in overlap_project_gamma', stdout, seedname)
    deallocate (svals, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating svals in overlap_project_gamma', stdout, seedname)
!~    deallocate(u_cmp,stat=ierr)
!~    if (ierr/=0) call io_error('Error in deallocating u_cmp in overlap_project_gamma')
    deallocate (u_matrix_r, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating u_matrix_r in overlap_project_gamma', stdout, seedname)

    if (timing_level > 1) call io_stopwatch('overlap: project_gamma', 2, stdout, seedname)

    return

  end subroutine overlap_project_gamma
![ysl-e]

end module w90_overlap
