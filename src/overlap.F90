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

module w90_overlap
  !! This module reads in the overlap (Mmn) and Projections (Amn)
  !! and performs simple operations on them.

  use w90_constants, only: dp, cmplx_0, cmplx_1
  use w90_parameters, only: disentanglement
  use w90_io, only: stdout
  use w90_comms, only: on_root, comms_bcast

  implicit none

  private

!~  public :: overlap_dis_read
  public :: overlap_allocate
  public :: overlap_read
  public :: overlap_dealloc
  public :: overlap_project
  public :: overlap_project_gamma  ![ysl]
!~  public :: overlap_check_m_symmetry

contains

  !%%%%%%%%%%%%%%%%%%%%%
  subroutine overlap_allocate()
    !%%%%%%%%%%%%%%%%%%%%%
    !! Allocate memory to read Mmn and Amn from files
    !! This must be called before calling overlap_read

    use w90_parameters, only: num_bands, num_wann, num_kpts, nntot, timing_level, &
      u_matrix, m_matrix_orig, m_matrix_orig_local, a_matrix, &
      u_matrix_opt, m_matrix, m_matrix_local
    use w90_io, only: io_error, io_stopwatch
    use w90_comms, only: my_node_id, num_nodes, comms_array_split

    integer :: ierr
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs

    if (timing_level > 0) call io_stopwatch('overlap: allocate', 1)

    call comms_array_split(num_kpts, counts, displs)

    allocate (u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating u_matrix in overlap_read')
    u_matrix = cmplx_0

    if (disentanglement) then
      if (on_root) then
        allocate (m_matrix_orig(num_bands, num_bands, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating m_matrix_orig in overlap_read')
      endif
      allocate (m_matrix_orig_local(num_bands, num_bands, nntot, counts(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating m_matrix_orig_local in overlap_read')
      allocate (a_matrix(num_bands, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating a_matrix in overlap_read')
      allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating u_matrix_opt in overlap_read')
    else
      if (on_root) then
        allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
        if (ierr /= 0) call io_error('Error in allocating m_matrix in overlap_read')
        m_matrix = cmplx_0
      endif
      allocate (m_matrix_local(num_wann, num_wann, nntot, counts(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating m_matrix_local in overlap_read')
      m_matrix_local = cmplx_0
    endif

    if (timing_level > 0) call io_stopwatch('overlap: allocate', 2)

  end subroutine overlap_allocate

  !%%%%%%%%%%%%%%%%%%%%%
  subroutine overlap_read()
    !%%%%%%%%%%%%%%%%%%%%%
    !! Read the Mmn and Amn from files
    !! Note: one needs to call overlap_allocate first!

    use w90_parameters, only: num_bands, num_wann, num_kpts, nntot, nncell, nnlist, &
      num_proj, lselproj, proj2wann_map, &
      devel_flag, u_matrix, m_matrix, a_matrix, timing_level, &
      m_matrix_orig, u_matrix_opt, cp_pp, use_bloch_phases, gamma_only, & ![ysl]
      m_matrix_local, m_matrix_orig_local, lhasproj
    use w90_io, only: io_file_unit, io_error, seedname, io_stopwatch
    use w90_comms, only: my_node_id, num_nodes, &
      comms_array_split, comms_scatterv

    implicit none

    integer :: nkp, nkp2, inn, nn, n, m, i, j
    integer :: mmn_in, amn_in, num_mmn, num_amn
    integer :: nnl, nnm, nnn, ncount
    integer :: nb_tmp, nkp_tmp, nntot_tmp, np_tmp, ierr
    real(kind=dp) :: m_real, m_imag, a_real, a_imag, mu_tmp, sigma_tmp
    complex(kind=dp), allocatable :: mmn_tmp(:, :)
    character(len=50) :: dummy
    logical :: nn_found
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs

    if (timing_level > 0) call io_stopwatch('overlap: read', 1)

    call comms_array_split(num_kpts, counts, displs)

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
        call io_error(trim(seedname)//'.mmn has not the right number of bands')
      if (nkp_tmp .ne. num_kpts) &
        call io_error(trim(seedname)//'.mmn has not the right number of k-points')
      if (nntot_tmp .ne. nntot) &
        call io_error(trim(seedname)//'.mmn has not the right number of nearest neighbours')

      ! Read the overlaps
      num_mmn = num_kpts*nntot
      allocate (mmn_tmp(num_bands, num_bands), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating mmn_tmp in overlap_read')
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
        do inn = 1, nntot
          if ((nkp2 .eq. nnlist(nkp, inn)) .and. &
              (nnl .eq. nncell(1, nkp, inn)) .and. &
              (nnm .eq. nncell(2, nkp, inn)) .and. &
              (nnn .eq. nncell(3, nkp, inn))) then
            if (.not. nn_found) then
              nn_found = .true.
              nn = inn
            else
              call io_error('Error reading '//trim(seedname)// &
                            '.mmn. More than one matching nearest neighbour found')
            endif
          endif
        end do
        if (nn .eq. 0) then
          if (on_root) write (stdout, '(/a,i8,2i5,i4,2x,3i3)') &
            ' Error reading '//trim(seedname)//'.mmn:', ncount, nkp, nkp2, nn, nnl, nnm, nnn
          call io_error('Neighbour not found')
        end if
        if (disentanglement) then
          m_matrix_orig(:, :, nn, nkp) = mmn_tmp(:, :)
        else
          ! disentanglement=.false. means numbands=numwann, so no the dimensions are the same
          m_matrix(:, :, nn, nkp) = mmn_tmp(:, :)
        end if
      end do
      deallocate (mmn_tmp, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating mmn_tmp in overlap_read')
      close (mmn_in)
    endif

    if (disentanglement) then
!       call comms_bcast(m_matrix_orig(1,1,1,1),num_bands*num_bands*nntot*num_kpts)
      call comms_scatterv(m_matrix_orig_local, num_bands*num_bands*nntot*counts(my_node_id), &
                          m_matrix_orig, num_bands*num_bands*nntot*counts, num_bands*num_bands*nntot*displs)
    else
!       call comms_bcast(m_matrix(1,1,1,1),num_wann*num_wann*nntot*num_kpts)
      call comms_scatterv(m_matrix_local, num_wann*num_wann*nntot*counts(my_node_id), &
                          m_matrix, num_wann*num_wann*nntot*counts, num_wann*num_wann*nntot*displs)
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
          call io_error(trim(seedname)//'.amn has not the right number of bands')
        if (nkp_tmp .ne. num_kpts) &
          call io_error(trim(seedname)//'.amn has not the right number of k-points')
        if (np_tmp .ne. num_proj) &
          call io_error(trim(seedname)//'.amn has not the right number of projections')

        if (num_proj > num_wann .and. .not. lselproj) &
          call io_error(trim(seedname)//'.amn has too many projections to be used without selecting a subset')

        ! Read the projections
        num_amn = num_bands*num_proj*num_kpts
        if (disentanglement) then
          do ncount = 1, num_amn
            read (amn_in, *, err=104, end=104) m, n, nkp, a_real, a_imag
            if (proj2wann_map(n) < 0) cycle
            a_matrix(m, proj2wann_map(n), nkp) = cmplx(a_real, a_imag, kind=dp)
          end do
        else
          do ncount = 1, num_amn
            read (amn_in, *, err=104, end=104) m, n, nkp, a_real, a_imag
            if (proj2wann_map(n) < 0) cycle
            u_matrix(m, proj2wann_map(n), nkp) = cmplx(a_real, a_imag, kind=dp)
          end do
        end if
        close (amn_in)
      endif

      if (disentanglement) then
        call comms_bcast(a_matrix(1, 1, 1), num_bands*num_wann*num_kpts)
      else
        call comms_bcast(u_matrix(1, 1, 1), num_wann*num_wann*num_kpts)
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
    if (cp_pp) call overlap_rotate()

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
        call overlap_project()
      else
        call overlap_project_gamma()
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

    if (timing_level > 0) call io_stopwatch('overlap: read', 2)

    return
101 call io_error('Error: Problem opening input file '//trim(seedname)//'.mmn')
102 call io_error('Error: Problem opening input file '//trim(seedname)//'.amn')
103 call io_error('Error: Problem reading input file '//trim(seedname)//'.mmn')
104 call io_error('Error: Problem reading input file '//trim(seedname)//'.amn')

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

  !%%%%%%%%%%%%%%%%%%%%%
  subroutine overlap_rotate
    !%%%%%%%%%%%%%%%%%%%%%
    !! Only used when interfaced to the CP code
    !! Not sure why this is done here and not in CP

    use w90_parameters, only: num_bands, a_matrix, m_matrix_orig, nntot, timing_level
    use w90_io, only: io_file_unit, io_error, io_stopwatch

    implicit none

    integer       :: lam_unit, info, inn, i, j
    real(kind=DP) :: lambda(num_bands, num_bands)
    real(kind=DP) :: AP(num_bands*(num_bands + 1)/2)
    real(kind=DP) :: eig(num_bands), work(3*num_bands)

    if (timing_level > 1) call io_stopwatch('overlap: rotate', 1)

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
      call io_error('Diagonalization of lambda in overlap_rotate failed')

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

    if (timing_level > 1) call io_stopwatch('overlap: rotate', 2)

    return

  end subroutine overlap_rotate

  !%%%%%%%%%%%%%%%%%%%%%
  subroutine overlap_dealloc()
    !%%%%%%%%%%%%%%%%%%%%%
    !! Dellocate memory

    use w90_parameters, only: u_matrix, m_matrix, m_matrix_orig, &
      a_matrix, u_matrix_opt, &
      m_matrix_local, m_matrix_orig_local
    use w90_io, only: io_error

    implicit none

    integer :: ierr

    if (allocated(u_matrix_opt)) then
      deallocate (u_matrix_opt, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating u_matrix_opt in overlap_dealloc')
    end if
    if (allocated(a_matrix)) then
      deallocate (a_matrix, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating a_matrix in overlap_dealloc')
    end if
    if (on_root) then
    if (allocated(m_matrix_orig)) then
      deallocate (m_matrix_orig, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_orig in overlap_dealloc')
    endif
    endif
    if (allocated(m_matrix_orig_local)) then
      deallocate (m_matrix_orig_local, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_orig_local in overlap_dealloc')
    endif
!~![ysl-b]
!~    if (allocated( ph_g)) then
!~       deallocate( ph_g, stat=ierr )
!~       if (ierr/=0) call io_error('Error deallocating ph_g in overlap_dealloc')
!~    endif
!~![ysl-e]

!    if (on_root) then
!    deallocate ( m_matrix, stat=ierr )
!    if (ierr/=0) call io_error('Error deallocating m_matrix in overlap_dealloc')
!    endif
!    deallocate ( m_matrix_local, stat=ierr )
!    if (ierr/=0) call io_error('Error deallocating m_matrix_local in overlap_dealloc')
!    deallocate ( u_matrix, stat=ierr )
!    if (ierr/=0) call io_error('Error deallocating u_matrix in overlap_dealloc')
    if (on_root) then
      if (allocated(m_matrix)) then
        deallocate (m_matrix, stat=ierr)
        if (ierr /= 0) call io_error('Error deallocating m_matrix in overlap_dealloc')
      endif
    endif
    if (allocated(m_matrix_local)) then
      deallocate (m_matrix_local, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating m_matrix_local in overlap_dealloc')
    endif
    if (allocated(u_matrix)) then
      deallocate (u_matrix, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating u_matrix in overlap_dealloc')
    endif

    return

  end subroutine overlap_dealloc

  !==================================================================!
  subroutine overlap_project()
    !==================================================================!
    !!  Construct initial guess from the projection via a Lowdin transformation
    !!  See section 3 of the CPC 2008
    !!  Note that in this subroutine num_wann = num_bands
    !!  since, if we are here, then disentanglement = FALSE
    !                                                                  !
    !                                                                  !
    !==================================================================!
    use w90_constants
    use w90_io, only: io_error, io_stopwatch
    use w90_parameters, only: num_bands, num_wann, num_kpts, timing_level, &
      u_matrix, m_matrix, nntot, nnlist, &
      m_matrix_local
    use w90_utility, only: utility_zgemm
    use w90_parameters, only: lsitesymmetry !RS:
    use w90_sitesym, only: sitesym_symmetrize_u_matrix !RS:
    use w90_comms, only: my_node_id, num_nodes, &
      comms_array_split, comms_scatterv, comms_gatherv

    implicit none

    ! internal variables
    integer :: i, j, m, nkp, info, ierr, nn, nkp2
    real(kind=dp), allocatable :: svals(:)
    real(kind=dp)                 :: rwork(5*num_bands)
    complex(kind=dp)              :: ctmp2
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:, :)
    complex(kind=dp), allocatable :: cvdag(:, :)
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs

    if (timing_level > 1) call io_stopwatch('overlap: project', 1)

    call comms_array_split(num_kpts, counts, displs)

    allocate (svals(num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating svals in overlap_project')
    allocate (cz(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cz in overlap_project')
    allocate (cvdag(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cvdag in overlap_project')
    allocate (cwork(4*num_bands), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cwork in overlap_project')

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
        call io_error('Error in ZGESVD in overlap_project')
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
            call io_error('Error in unitarity of initial U in overlap_project')
          endif
          if ((i .ne. j) .and. (abs(ctmp2) .gt. eps5)) then
            write (stdout, *) ' ERROR: unitarity of initial U'
            write (stdout, '(1x,a,i2)') 'nkp= ', nkp
            write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
            write (stdout, '(1x,a,f12.6,1x,f12.6)') &
              '[u_matrix.transpose(u_matrix)]_ij= ', &
              real(ctmp2, dp), aimag(ctmp2)
            call io_error('Error in unitarity of initial U in overlap_project')
          endif
        enddo
      enddo
    enddo
    ! NKP

    if (lsitesymmetry) call sitesym_symmetrize_u_matrix(num_wann, u_matrix) !RS: update U(Rk)

    ! so now we have the U's that rotate the wavefunctions at each k-point.
    ! the matrix elements M_ij have also to be updated
    do nkp = 1, counts(my_node_id)
      do nn = 1, nntot
        nkp2 = nnlist(nkp + displs(my_node_id), nn)
        ! cvdag = U^{dagger} . M   (use as workspace)
        call utility_zgemm(cvdag, u_matrix(:, :, nkp + displs(my_node_id)), 'C', m_matrix_local(:, :, nn, nkp), 'N', num_wann)
        ! cz = cvdag . U
        call utility_zgemm(cz, cvdag, 'N', u_matrix(:, :, nkp2), 'N', num_wann)
        m_matrix_local(:, :, nn, nkp) = cz(:, :)
      end do
    end do
    call comms_gatherv(m_matrix_local, num_wann*num_wann*nntot*counts(my_node_id), &
                       m_matrix, num_wann*num_wann*nntot*counts, num_wann*num_wann*nntot*displs)

    deallocate (cwork, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cwork in overlap_project')
    deallocate (cvdag, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cvdag in overlap_project')
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cz in overlap_project')
    deallocate (svals, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating svals in overlap_project')

    if (timing_level > 1) call io_stopwatch('overlap: project', 2)

    return

  end subroutine overlap_project

![ysl-b]
  !==================================================================!
  subroutine overlap_project_gamma()
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
    use w90_parameters, only: num_wann, timing_level, &
      u_matrix, m_matrix, nntot!,num_kpts,nnlist
    use w90_utility, only: utility_zgemm

    implicit none

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

    if (timing_level > 1) call io_stopwatch('overlap: project_gamma', 1)

!~    allocate(ph_g(num_wann),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating ph_g in overlap_project_gamma')
    ! internal variables
    allocate (u_matrix_r(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating u_matrix_r in overlap_project_gamma')
!~    allocate(u_cmp(num_wann),stat=ierr)
!~    if (ierr/=0) call io_error('Error in allocating u_cmp in overlap_project_gamma')
    allocate (svals(num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating svals in overlap_project_gamma')
    allocate (work(5*num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating work in overlap_project_gamma')
    allocate (rz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rz in overlap_project_gamma')
    allocate (rv(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating rv in overlap_project_gamma')
    allocate (cz(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cz in overlap_project_gamma')
    allocate (cvdag(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating cvdag in overlap_project_gamma')

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
      call io_error('overlap_project_gamma: problem in DGESVD 1')
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
          call io_error('Error in unitarity of initial U in overlap_project_gamma')
        endif
        if ((i .ne. j) .and. (abs(rtmp2) .gt. eps5)) then
          write (stdout, *) ' ERROR: unitarity of initial U'
          write (stdout, '(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j
          write (stdout, '(1x,a,f12.6,1x,f12.6)') &
            '[u_matrix.transpose(u_matrix)]_ij= ', &
            rtmp2
          call io_error('Error in unitarity of initial U in overlap_project_gamma')
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
    if (ierr /= 0) call io_error('Error in deallocating cvdag in overlap_project_gamma')
    deallocate (cz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating cz in overlap_project_gamma')
    deallocate (rv, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rv in overlap_project_gamma')
    deallocate (rz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating rz in overlap_project_gamma')
    deallocate (work, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating work in overlap_project_gamma')
    deallocate (svals, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating svals in overlap_project_gamma')
!~    deallocate(u_cmp,stat=ierr)
!~    if (ierr/=0) call io_error('Error in deallocating u_cmp in overlap_project_gamma')
    deallocate (u_matrix_r, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating u_matrix_r in overlap_project_gamma')

    if (timing_level > 1) call io_stopwatch('overlap: project_gamma', 2)

    return

  end subroutine overlap_project_gamma
![ysl-e]

end module w90_overlap
