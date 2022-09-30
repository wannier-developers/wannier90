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
  public :: overlap_read
  ! these are only public for old wannier_lib
  public :: overlap_project
  public :: overlap_project_gamma

contains

  !================================================!

  subroutine overlap_allocate(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, &
                              m_matrix_orig_local, u_matrix, u_matrix_opt, nntot, num_bands, &
                              num_kpts, num_wann, timing_level, timer, dist_k, error, comm)
    !================================================!
    !! Allocate memory to read Mmn and Amn from files
    !! This must be called before calling overlap_read
    !
    !================================================!

    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: timer_list_type
    use w90_error

    ! arguments
    integer, intent(in) :: nntot
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: timing_level
    integer, intent(in) :: dist_k(:)

    complex(kind=dp), allocatable :: a_matrix(:, :, :)
    complex(kind=dp), allocatable :: m_matrix(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_local(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), allocatable :: m_matrix_orig_local(:, :, :, :)
    complex(kind=dp), allocatable :: u_matrix(:, :, :)
    complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)

    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    integer :: ierr
    integer :: my_node_id, nkl
    logical :: disentanglement
    logical :: on_root = .false.

    disentanglement = (num_bands > num_wann)

    my_node_id = mpirank(comm)
    nkl = count(dist_k == my_node_id) ! number of k on this rank

    if (my_node_id == 0) on_root = .true.

    if (timing_level > 0) call io_stopwatch_start('overlap: allocate', timer)

    if (disentanglement) then
      if (on_root) then
        allocate (m_matrix_orig(num_bands, num_bands, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) then
          call set_error_alloc(error, 'Error in allocating m_matrix_orig in overlap_allocate', comm)
          return
        endif
      else
        allocate (m_matrix_orig(0, 0, 0, 0))
      endif
      allocate (m_matrix_orig_local(num_bands, num_bands, nntot, nkl), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating m_matrix_orig_local in overlap_allocate', comm)
        return
      endif
      m_matrix_orig = cmplx_0
      m_matrix_orig_local = cmplx_0
    else
      allocate (m_matrix_orig_local(0, 0, 0, 0))
      allocate (m_matrix_orig(0, 0, 0, 0))
    endif

    if (on_root) then
      allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating m_matrix in overlap_allocate', comm)
        return
      endif
      m_matrix = cmplx_0
    else
      allocate (m_matrix(0, 0, 0, 0))
    endif
    allocate (m_matrix_local(num_wann, num_wann, nntot, nkl), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating m_matrix_local in overlap_allocate', comm)
      return
    endif
    m_matrix_local = cmplx_0

    allocate (a_matrix(num_bands, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating a_matrix in overlap_allocate', comm)
      return
    endif
    allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating u_matrix in overlap_allocate', comm)
      return
    endif
    u_matrix = cmplx_0
    allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating u_matrix_opt in overlap_allocate', comm)
      return
    endif
    u_matrix_opt = cmplx_0

    if (timing_level > 0) call io_stopwatch_stop('overlap: allocate', timer)

  end subroutine overlap_allocate

  !================================================!
  subroutine overlap_read(kmesh_info, select_projection, au_matrix, m_matrix_local, num_bands, &
                          num_kpts, num_proj, num_wann, print_output, timing_level, cp_pp, &
                          use_bloch_phases, seedname, stdout, timer, dist_k, error, comm)
    !================================================!
    !! Read the Mmn and Amn from files
    !! Note: one needs to call overlap_allocate first!
    !
    !================================================!

    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: kmesh_info_type, print_output_type, timer_list_type
    use w90_wannier90_types, only: select_projection_type
    use w90_error

    implicit none

    ! arguments
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(print_output_type), intent(in) :: print_output
    type(select_projection_type), intent(in) :: select_projection
    type(timer_list_type), intent(inout) :: timer
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: dist_k(:)
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_proj
    integer, intent(in) :: num_wann
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level

    complex(kind=dp), intent(inout) :: au_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :)

    logical, intent(in) :: cp_pp, use_bloch_phases

    character(len=50), intent(in) :: seedname

    ! local variables
    complex(kind=dp), allocatable :: mmn_tmp(:, :)
    real(kind=dp) :: m_real, m_imag, a_real, a_imag

    integer, allocatable :: map_kpts(:)
    integer :: mmn_in, amn_in, num_mmn, num_amn
    integer :: my_node_id
    integer :: nb_tmp, nkp_tmp, nntot_tmp, np_tmp, ierr
    integer :: nkp, nkp2, nkp_loc, inn, nn, n, m
    integer :: nnl, nnm, nnn, ncount

    logical :: disentanglement
    logical :: nn_found
    logical :: on_root = .false.

    character(len=50) :: dummy

    disentanglement = (num_bands > num_wann)

    my_node_id = mpirank(comm)
    if (my_node_id == 0) then
      on_root = .true.
    endif
    allocate (map_kpts(num_kpts))
    nkp_loc = 1
    do nkp = 1, num_kpts
      if (dist_k(nkp) == my_node_id) then
        map_kpts(nkp) = nkp_loc
        nkp_loc = nkp_loc + 1
      endif
    enddo

    if (timing_level > 0) call io_stopwatch_start('overlap: read', timer)

    !if (on_root) then - read on local bits on all nodes

    ! Read M_matrix_orig from file
    open (newunit=mmn_in, file=trim(seedname)//'.mmn', &
          form='formatted', status='old', action='read', err=101)

    if (print_output%iprint > 0) write (stdout, '(/a)', advance='no') ' Reading overlaps from '//trim(seedname)//'.mmn    : '

    ! Read the comment line
    read (mmn_in, '(a)', err=103, end=103) dummy
    if (print_output%iprint > 0) write (stdout, '(a)') trim(dummy)

    ! Read the number of bands, k-points and nearest neighbours
    read (mmn_in, *, err=103, end=103) nb_tmp, nkp_tmp, nntot_tmp

    ! Checks
    if (nb_tmp .ne. num_bands) then
      call set_error_file(error, trim(seedname)//'.mmn has not the right number of bands', comm)
      return
    endif
    if (nkp_tmp .ne. num_kpts) then
      call set_error_file(error, trim(seedname)//'.mmn has not the right number of k-points', comm)
      return
    endif
    if (nntot_tmp .ne. kmesh_info%nntot) then
      write (*, *) my_node_id, nntot_tmp, kmesh_info%nntot
      call set_error_file(error, trim(seedname)//'.mmn has not the right number of nearest neighbours', comm)
      return
    endif

    ! Read the overlaps
    num_mmn = num_kpts*kmesh_info%nntot
    allocate (mmn_tmp(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating mmn_tmp in overlap_read', comm)
      return
    endif
    do ncount = 1, num_mmn
      read (mmn_in, *, err=103, end=103) nkp, nkp2, nnl, nnm, nnn
      do n = 1, num_bands
        do m = 1, num_bands
          read (mmn_in, *, err=103, end=103) m_real, m_imag
          mmn_tmp(m, n) = cmplx(m_real, m_imag, kind=dp)
        enddo
      enddo
      if (dist_k(nkp) == my_node_id) then
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
              call set_error_file(error, 'Error reading '//trim(seedname)// &
                                  '.mmn. More than one matching nearest neighbour found', comm)
              return
            endif
          endif
        end do
        if (nn .eq. 0) then
          if (on_root) write (stdout, '(/a,i8,2i5,i4,2x,3i3)') &
            ' Error reading '//trim(seedname)//'.mmn:', ncount, nkp, nkp2, nn, nnl, nnm, nnn
          call set_error_file(error, 'Neighbour not found', comm)
          return
        end if
        m_matrix_local(:, :, nn, map_kpts(nkp)) = mmn_tmp(:, :)
      end if
    end do
    deallocate (mmn_tmp, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating mmn_tmp in overlap_read', comm)
      return
    endif
    close (mmn_in)

    !if (disentanglement) then
    !  w = num_bands*num_bands*kmesh_info%nntot
    !  call comms_scatterv(m_matrix_orig_local, w*counts(my_node_id), m_matrix_orig, w*counts, w*displs, error, comm)
    !  if (allocated(error)) return
    !else
    !  w = num_wann*num_wann*kmesh_info%nntot
    !  call comms_scatterv(m_matrix_local, w*counts(my_node_id), m_matrix, w*counts, w*displs, error, comm)
    !  if (allocated(error)) return
    !endif

    if (.not. use_bloch_phases) then
      !if (on_root) then read on all nodes

      ! Read A_matrix from file wannier.amn
      open (newunit=amn_in, file=trim(seedname)//'.amn', form='formatted', status='old', err=102)

      if (print_output%iprint > 0) write (stdout, '(/a)', advance='no') ' Reading projections from '//trim(seedname)//'.amn : '

      ! Read the comment line
      read (amn_in, '(a)', err=104, end=104) dummy
      if (print_output%iprint > 0) write (stdout, '(a)') trim(dummy)

      ! Read the number of bands, k-points and wannier functions
      read (amn_in, *, err=104, end=104) nb_tmp, nkp_tmp, np_tmp

      ! Checks
      if (nb_tmp .ne. num_bands) then
        call set_error_file(error, trim(seedname)//'.amn has not the right number of bands', comm)
        return
      endif
      if (nkp_tmp .ne. num_kpts) then
        call set_error_file(error, trim(seedname)//'.amn has not the right number of k-points', comm)
        return
      endif
      if (np_tmp .ne. num_proj) then
        call set_error_file(error, trim(seedname)//'.amn has not the right number of projections', comm)
        return
      endif

      if (num_proj > num_wann .and. .not. select_projection%lselproj) then
        call set_error_file(error, trim(seedname)//'.amn has too many projections to be used without selecting a subset', comm)
        return
      endif

      if (.not. allocated(select_projection%proj2wann_map)) then
        call set_error_fatal(error, 'select_projection%proj2wann_map not allocated in overlap_read call', comm)
        return
      endif

      ! Read the projections
      num_amn = num_bands*num_proj*num_kpts
      do ncount = 1, num_amn
        read (amn_in, *, err=104, end=104) m, n, nkp, a_real, a_imag
        if (select_projection%proj2wann_map(n) < 0) cycle
        au_matrix(m, select_projection%proj2wann_map(n), nkp) = cmplx(a_real, a_imag, kind=dp)
      end do
      close (amn_in)
      !endif

      !call comms_bcast(au_matrix(1, 1, 1), num_bands*num_wann*num_kpts, error, comm)
      !if (allocated(error)) return

    else

      do n = 1, num_kpts
        do m = 1, num_wann
          au_matrix(m, m, n) = cmplx_1
        end do
      end do

    end if

    ! If post-processing a Car-Parinello calculation (gamma only)
    ! then rotate M and A to the basis of Kohn-Sham eigenstates
    if (cp_pp) call overlap_rotate(au_matrix, m_matrix_local, kmesh_info%nntot, num_bands, &
                                   timing_level, timer, error, comm)
    if (allocated(error)) return

    deallocate (map_kpts)
    if (timing_level > 0) call io_stopwatch_stop('overlap: read', timer)

    return
101 call set_error_file(error, 'Error: Problem opening input file '//trim(seedname)//'.mmn', comm)
    return
102 call set_error_file(error, 'Error: Problem opening input file '//trim(seedname)//'.amn', comm)
    return
103 call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.mmn', comm)
    return
104 call set_error_file(error, 'Error: Problem reading input file '//trim(seedname)//'.amn', comm)
    return

    !if (on_root) deallocate(m_matrix_orig)

  end subroutine overlap_read

!~[aam]
!~  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~ subroutine overlap_check_m_symmetry()
!~ !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~
!~   use w90_constants,  only : eps6
!~   use w90_parameters, only : num_bands,num_wann,a_matrix,m_matrix_orig,u_matrix,m_matrix, &
!~                              nntot,timing_level,disentanglement
!~   use w90_io,         only : io_error,io_stopwatch
!~
!~   implicit none
!~
!~   integer       :: nn,i,j,m,n, p(1), mdev, ndev, nndev,ierr, num_tmp
!~  real(kind=dp) :: dev,dev_tmp
!~
!~   if (timing_level>1) call io_stopwatch('overlap: check_m_sym',1)
!~
!~   if (disentanglement) then
!~      num_tmp=num_bands
!~   else
!~      num_tmp=num_wann
!~   endif
!~
!~   ! check whether M_mn is symmetric
!~   dev = 0.0_dp
!~   do nn=1,nntot
!~      do n=1,num_tmp
!~         do m=1,n
!~            if (disentanglement) then
!~               dev_tmp=abs(m_matrix_orig(m,n,nn,1)-m_matrix_orig(n,m,nn,1))
!~            else
!~               dev_tmp=abs(m_matrix(m,n,nn,1)-m_matrix(n,m,nn,1))
!~            endif
!~            if ( dev_tmp .gt. dev ) then
!~               dev = dev_tmp
!~               mdev  = m ; ndev  = n ;  nndev  = nn
!~            end if
!~         end do
!~      end do
!~   end do
!~
!~   if ( dev .gt. eps6 ) then
!~      write(stdout,'(1x,"+",76("-"),"+")')
!~      write(stdout,'(3x,a)') 'WARNING: M is not strictly symmetric in overlap_check_m_symmetry'
!~      write(stdout,'(3x,a,f12.8)') &
!~           'Largest deviation |M_mn-M_nm| at k : ',dev
!~      write(stdout,'(3(a5,i4))') &
!~           '   m=',mdev,',  n=',ndev,',  k=',nndev
!~      write(stdout,'(1x,"+",76("-"),"+")')
!~   end if
!~
!~   if (timing_level>1) call io_stopwatch('overlap: check_m_sym',2)
!~
!~   return
!~
!~ end subroutine overlap_check_m_symmetry
!~[aam]

!~![ysl-b]
!~  !%%%%%%%%%%%%%%%%%%%%%
!~  subroutine overlap_symmetrize
!~  !%%%%%%%%%%%%%%%%%%%%%
!~
!~    use w90_parameters, only : num_bands,num_wann,a_matrix,m_matrix_orig,u_matrix,m_matrix, &
!~                               nntot,timing_level,disentanglement
!~    use w90_io,         only : io_error,io_stopwatch
!~
!~    implicit none
!~
!~    integer       :: nn,i,j,m,n, p(1), mdev, ndev, nndev,ierr
!~    real(kind=dp) :: eps,dev,dev_tmp
!~    real(kind=dp),    allocatable :: a_cmp(:)
!~
!~
!~    if (timing_level>1) call io_stopwatch('overlap: symmetrize',1)
!~
!~    allocate(a_cmp(num_wann),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating a_cmp in overlap_symmetrize')
!~    allocate(ph_g(num_bands),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating ph_g in overlap_symmetrize')
!~
!~    eps = 1.0e-6_dp
!~    ph_g=cmplx_0
!~
!~    ! disentanglement - m_matrix_orig
!~    ! localzation only - m_matrix
!~
!~    ! If a wavefunction is real except for a phase factor e^(i*phi_m) = ph_g(m)
!~    ! A_mn = ph_g(m)^(-1)*<m_R|g_n>   (m_R implies a real wave function)
!~    ! At m_th row, find five elements with largest complex component
!~    ! -> ph_g(m) is calculated from those elements
!~    !
!~   do m=1,num_bands
!~       a_cmp(:)=abs(a_matrix(m,:,1))
!~       p=maxloc(a_cmp)
!~       ph_g(m)=conjg(a_matrix(m,p(1),1)/abs(a_matrix(m,p(1),1)))
!~   end do
!~
!~    ! M_mn (new) = ph_g(m) * M_mn (old) * conjg(ph_g(n))
!~    ! a_mn (new) = ph_g(m) * a_mn (old)
!~
!~    do nn=1,nntot
!~       do n=1,num_bands
!~          do m=1,num_bands
!~             m_matrix_orig(m,n,nn,1)=ph_g(m)*m_matrix_orig(m,n,nn,1)*conjg(ph_g(n))
!~          end do
!~       end do
!~    end do
!~    do n=1,num_wann
!~       do m=1,num_bands
!~          a_matrix(m,n,1)=ph_g(m)*a_matrix(m,n,1)
!~       end do
!~    end do
!~
!~    ! check whether M_mn is now symmetric
!~    dev = 0.0_dp
!~    do nn=1,nntot
!~       do n=1,num_bands
!~          do m=1,n
!~             dev_tmp=abs(m_matrix_orig(m,n,nn,1)-m_matrix_orig(n,m,nn,1))
!~             if ( dev_tmp .gt. dev ) then
!~                dev = dev_tmp
!~                mdev  = m ; ndev  = n ;  nndev  = nn
!~             end if
!~          end do
!~       end do
!~    end do
!~    if ( dev .gt. eps ) then
!~       write(stdout,'(1x,"+",76("-"),"+")')
!~       write(stdout,'(3x,a)') 'WARNING: M is not strictly symmetric in overlap_symmetrize'
!~       write(stdout,'(3x,a,f12.8)') &
!~            'Largest deviation |M_mn-M_nm| at k : ',dev
!~       write(stdout,'(3(a5,i4))') &
!~            '   m=',mdev,',  n=',ndev,',  k=',nndev
!~       write(stdout,'(1x,"+",76("-"),"+")')
!~    end if
!~
!~    deallocate(a_cmp,stat=ierr)
!~    if (ierr/=0) call io_error('Error in deallocating a_cmp in overlap_symmetrize')
!~
!~    if (timing_level>1) call io_stopwatch('overlap: symmetrize',2)
!~
!~    return
!~
!~  end subroutine overlap_symmetrize
!~![ysl-e]

  !================================================!
  subroutine overlap_rotate(a_matrix, m_matrix_orig, nntot, num_bands, timing_level, timer, error, &
                            comm)
    !================================================!
    !
    !! Only used when interfaced to the CP code
    !! Not sure why this is done here and not in CP
    !
    !================================================!

    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_error, only: w90_error_type, set_error_fatal
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: nntot
    integer, intent(in) :: num_bands
    integer, intent(in) :: timing_level

    complex(kind=dp), intent(inout) :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), intent(inout) :: a_matrix(:, :, :)

    ! local variables
    integer       :: lam_unit, info, inn, i, j
    real(kind=dp) :: lambda(num_bands, num_bands)
    real(kind=dp) :: AP(num_bands*(num_bands + 1)/2)
    real(kind=dp) :: eig(num_bands), work(3*num_bands)

    if (timing_level > 1) call io_stopwatch_start('overlap: rotate', timer)

    open (newunit=lam_unit, file='lambda.dat', &
          form='unformatted', status='old', action='read')
    read (lam_unit) lambda
    close (lam_unit)

    do j = 1, num_bands
      do i = 1, j
        AP(i + (j - 1)*j/2) = 0.5_dp*(lambda(i, j) + lambda(j, i))
      end do
    end do

    CALL DSPEV('V', 'U', num_bands, AP, eig, lambda, num_bands, work, info)
    if (info .ne. 0) then
      call set_error_fatal(error, 'Diagonalization of lambda in overlap_rotate failed', comm)
      return
    endif

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

    if (timing_level > 1) call io_stopwatch_stop('overlap: rotate', timer)

    return

  end subroutine overlap_rotate

  !================================================!
  subroutine overlap_dealloc(a_matrix, m_matrix, m_matrix_local, m_matrix_orig, &
                             m_matrix_orig_local, u_matrix, u_matrix_opt, error, comm)
    !================================================!
    !
    !! Dellocate memory
    !
    !================================================!

    !use w90_io, only: io_error
    use w90_error, only: set_error_dealloc, w90_error_type

    implicit none

    ! arguments
    complex(kind=dp), allocatable, intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix_orig(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: a_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: u_matrix_opt(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix_local(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: m_matrix_orig_local(:, :, :, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    integer :: ierr
    logical :: on_root = .false.

    if (mpirank(comm) == 0) on_root = .true.

    if (allocated(u_matrix_opt)) then
      deallocate (u_matrix_opt, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating u_matrix_opt in overlap_dealloc', comm)
        return
      endif
    end if
    if (allocated(a_matrix)) then
      deallocate (a_matrix, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating a_matrix in overlap_dealloc', comm)
        return
      endif
    end if
    if (allocated(m_matrix_orig)) then
      deallocate (m_matrix_orig, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating m_matrix_orig in overlap_dealloc', comm)
        return
      endif
    endif
    if (allocated(m_matrix_orig_local)) then
      deallocate (m_matrix_orig_local, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating m_matrix_orig_local in overlap_dealloc', comm)
        return
      endif
    endif
    if (on_root) then
      if (allocated(m_matrix)) then
        deallocate (m_matrix, stat=ierr)
        if (ierr /= 0) then
          call set_error_dealloc(error, 'Error deallocating m_matrix in overlap_dealloc', comm)
          return
        endif
      endif
    endif
    if (allocated(m_matrix_local)) then
      deallocate (m_matrix_local, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating m_matrix_local in overlap_dealloc', comm)
        return
      endif
    endif
    if (allocated(u_matrix)) then
      deallocate (u_matrix, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating u_matrix in overlap_dealloc', comm)
        return
      endif
    endif

    return

  end subroutine overlap_dealloc

  !================================================!
  subroutine overlap_project(sitesym, m_matrix_local, u_matrix, nnlist, nntot, &
                             num_bands, num_kpts, num_wann, timing_level, lsitesymmetry, stdout, &
                             timer, dist_k, error, comm)
    !================================================!
    !!  Construct initial guess from the projection via a Lowdin transformation
    !!  See section 3 of the CPC 2008
    !!  Note that in this subroutine num_wann = num_bands
    !!  since, if we are here, then disentanglement = FALSE
    !
    !================================================!
    use w90_constants
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_error, only: w90_error_type, set_error_alloc, set_error_fatal, set_error_dealloc, &
      set_error_fatal
    use w90_utility, only: utility_zgemm
    use w90_sitesym, only: sitesym_symmetrize_u_matrix
    use w90_wannier90_types, only: sitesym_type
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    type(sitesym_type), intent(in) :: sitesym
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: dist_k(:)
    integer, intent(in) :: nnlist(:, :)
    integer, intent(in) :: nntot
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: num_wann
    integer, intent(in) :: timing_level
    integer, intent(in) :: stdout

    !complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(inout) :: m_matrix_local(:, :, :, :)

    logical, intent(in) :: lsitesymmetry

    ! local variables
    integer :: i, j, m, nkp, nkp_loc, info, ierr, nn, nkp2
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp)                 :: rwork(5*num_bands)
    complex(kind=dp)              :: ctmp2
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cvdag(:, :)

    ! pllel setup
    integer :: my_node_id
    logical :: on_root = .false.

    my_node_id = mpirank(comm)
    if (my_node_id == 0) on_root = .true.

    if (timing_level > 1) call io_stopwatch_start('overlap: project', timer)

    allocate (svals(num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating svals in overlap_project', comm)
      return
    endif
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cz in overlap_project', comm)
      return
    endif
    allocate (cvdag(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cvdag in overlap_project', comm)
      return
    endif
    allocate (cwork(4*num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cwork in overlap_project', comm)
      return
    endif

    ! Calculate the transformation matrix CU = CS^(-1/2).CA,
    ! where CS = CA.CA^\dagger.

    do nkp = 1, num_kpts
      if (dist_k(nkp) == my_node_id) then
        !
        ! SINGULAR VALUE DECOMPOSITION
        !
        call zgesvd('A', 'A', num_bands, num_bands, u_matrix(1, 1, nkp), num_bands, svals, cz, &
                    num_bands, cvdag, num_bands, cwork, 4*num_bands, rwork, info)
        if (info .ne. 0) then
          write (stdout, *) ' ERROR: IN ZGESVD IN overlap_project'
          write (stdout, *) ' K-POINT NKP=', nkp, ' INFO=', info
          if (info .lt. 0) then
            write (stdout, *) ' THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
          endif
          call set_error_fatal(error, 'Error in ZGESVD in overlap_project', comm)
          return
        endif

        ! u_matrix(:,:,nkp)=matmul(cz,cvdag)
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
              call set_error_fatal(error, 'Error in unitarity of initial U in overlap_project', comm)
              return
            endif
            if ((i .ne. j) .and. (abs(ctmp2) .gt. eps5)) then
              write (stdout, *) ' ERROR: unitarity of initial U'
              write (stdout, '(1x,a,i2)') 'nkp= ', nkp
              write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
              write (stdout, '(1x,a,f12.6,1x,f12.6)') &
                '[u_matrix.transpose(u_matrix)]_ij= ', &
                real(ctmp2, dp), aimag(ctmp2)
              call set_error_fatal(error, 'Error in unitarity of initial U in overlap_project', comm)
              return
            endif
          enddo
        enddo
      else
        u_matrix(:, :, nkp) = 0.0_dp
      endif
    enddo
    ! NKP
    call comms_allreduce(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts, 'SUM', error, comm)
    if (allocated(error)) return

    if (lsitesymmetry) then
      call sitesym_symmetrize_u_matrix(sitesym, u_matrix, num_bands, num_wann, num_kpts, num_wann, &
                                       stdout, error, comm) !RS: update U(Rk)
      if (allocated(error)) return
    endif

    ! so now we have the U's that rotate the wavefunctions at each k-point.
    ! the matrix elements M_ij have also to be updated
    nkp_loc = 1
    do nkp = 1, num_kpts
      if (dist_k(nkp) == my_node_id) then
        do nn = 1, nntot
          nkp2 = nnlist(nkp, nn)
          ! cvdag = U^{dagger} . M   (use as workspace)
          call utility_zgemm(cvdag, u_matrix(:, :, nkp), 'C', &
                             m_matrix_local(:, :, nn, nkp_loc), 'N', num_wann)
          ! cz = cvdag . U
          call utility_zgemm(cz, cvdag, 'N', u_matrix(:, :, nkp2), 'N', num_wann)
          m_matrix_local(:, :, nn, nkp_loc) = cz(:, :)
        end do
        nkp_loc = nkp_loc + 1
      endif
    end do
    !call comms_reduce(m_matrix(1,1,1,1), num_wann*num_wann*nntot*num_kpts, 'SUM', error, comm)
    !if (allocated(error)) return

    deallocate (cwork, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cwork in overlap_project', comm)
      return
    endif
    deallocate (cvdag, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cvdag in overlap_project', comm)
      return
    endif
    deallocate (cz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cz in overlap_project', comm)
      return
    endif
    deallocate (svals, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating svals in overlap_project', comm)
      return
    endif

    if (timing_level > 1) call io_stopwatch_stop('overlap: project', timer)

    return

  end subroutine overlap_project

![ysl-b]
  !================================================!
  subroutine overlap_project_gamma(m_matrix, u_matrix, nntot, num_wann, timing_level, stdout, &
                                   timer, error, comm)
    !================================================!
    !!  Construct initial guess from the projection via a Lowdin transformation
    !!  See section 3 of the CPC 2008
    !!  Note that in this subroutine num_wann = num_bands
    !!  since, if we are here, then disentanglement = FALSE
    !!  Gamma specific version
    !
    !================================================!
    use w90_constants
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_error, only: w90_error_type, set_error_alloc, set_error_fatal, set_error_dealloc, &
      set_error_fatal
    use w90_utility, only: utility_zgemm
    use w90_types, only: timer_list_type

    implicit none

    ! arguments
    integer, intent(in) :: nntot
    integer, intent(in) :: stdout
    integer, intent(in) :: timing_level
    integer, intent(in) :: num_wann
    complex(kind=dp), intent(inout) :: m_matrix(:, :, :, :)
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    ! internal variables
    integer :: i, j, m, info, ierr, nn
    real(kind=dp) :: rtmp2
    real(kind=dp), allocatable :: u_matrix_r(:, :)
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp), allocatable :: work(:)
    real(kind=dp), allocatable :: rz(:, :)
    real(kind=dp), allocatable :: rv(:, :)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cvdag(:, :)

    if (timing_level > 1) call io_stopwatch_start('overlap: project_gamma', timer)

    allocate (u_matrix_r(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating u_matrix_r in overlap_project_gamma', comm)
      return
    endif
    allocate (svals(num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating svals in overlap_project_gamma', comm)
      return
    endif
    allocate (work(5*num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating work in overlap_project_gamma', comm)
      return
    endif
    allocate (rz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating rz in overlap_project_gamma', comm)
      return
    endif
    allocate (rv(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating rv in overlap_project_gamma', comm)
      return
    endif
    allocate (cz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cz in overlap_project_gamma', comm)
      return
    endif
    allocate (cvdag(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating cvdag in overlap_project_gamma', comm)
      return
    endif

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
    call dgesvd('A', 'A', num_wann, num_wann, u_matrix_r, num_wann, svals, rz, num_wann, rv, &
                num_wann, work, 5*num_wann, info)
    if (info .ne. 0) then
      write (stdout, *) ' ERROR: IN DGESVD IN overlap_project_gamma'
      if (info .lt. 0) then
        write (stdout, *) 'THE ', -info, '-TH ARGUMENT HAD ILLEGAL VALUE'
      endif
      call set_error_fatal(error, 'overlap_project_gamma: problem in DGESVD 1', comm)
      return
    endif

    call dgemm('N', 'N', num_wann, num_wann, num_wann, 1.0_dp, rz, num_wann, rv, num_wann, 0.0_dp, &
               u_matrix_r, num_wann)
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
          call set_error_fatal(error, 'Error in unitarity of initial U in overlap_project_gamma', comm)
          return
        endif
        if ((i .ne. j) .and. (abs(rtmp2) .gt. eps5)) then
          write (stdout, *) ' ERROR: unitarity of initial U'
          write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
          write (stdout, '(1x,a,f12.6,1x,f12.6)') &
            '[u_matrix.transpose(u_matrix)]_ij= ', &
            rtmp2
          call set_error_fatal(error, 'Error in unitarity of initial U in overlap_project_gamma', comm)
          return
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
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cvdag in overlap_project_gamma', comm)
      return
    endif
    deallocate (cz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating cz in overlap_project_gamma', comm)
      return
    endif
    deallocate (rv, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating rv in overlap_project_gamma', comm)
      return
    endif
    deallocate (rz, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating rz in overlap_project_gamma', comm)
      return
    endif
    deallocate (work, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating work in overlap_project_gamma', comm)
      return
    endif
    deallocate (svals, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating svals in overlap_project_gamma', comm)
      return
    endif
    deallocate (u_matrix_r, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating u_matrix_r in overlap_project_gamma', comm)
      return
    endif

    if (timing_level > 1) call io_stopwatch_stop('overlap: project_gamma', timer)

    return

  end subroutine overlap_project_gamma
![ysl-e]

end module w90_overlap
