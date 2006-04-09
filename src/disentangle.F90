!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
! Copyright (C) 2004,2006 Jonathan Yates, Arash Mostofi,
!            Nicola Marzari, Ivo Souza, David Vanderbilt
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------

module disentangle

  use constants, only: dp,twopi,cmplx_i,cmplx_0,cmplx_1
  use io, only: io_error,stdout
  use parameters

  implicit none

  private

  ! Variables local to this module
  complex(kind=dp), allocatable :: clamp(:,:,:)
  real(kind=dp),    allocatable :: eigval_opt(:,:)

  ! Determined in dis_windows
  logical                :: linner                               
  logical, allocatable   :: lfrozen(:,:) 
  integer, allocatable   :: nfirstwin(:)
  integer, allocatable   :: ndimfroz(:)
  integer, allocatable   :: indxfroz(:,:)
  integer, allocatable   :: indxnfroz(:,:)
  


  public :: dis_main


contains


  !==================================================================!
  subroutine dis_main()
  !==================================================================!
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !==================================================================!  


    ! internal variables
    integer :: nkp,nn,i,j,m,n,l
    integer :: nkp2,info,ierr

    complex(kind=dp) :: ctmp
    complex(kind=dp), allocatable :: cmtmp(:,:)
    complex(kind=dp), allocatable :: caa(:,:,:)

    ! For ZGESVD
    real(kind=dp),    allocatable :: svals(:)
    real(kind=dp),    allocatable :: rwork(:)
    complex(kind=dp), allocatable :: cv(:,:)
    complex(kind=dp), allocatable :: cz(:,:)
    complex(kind=dp), allocatable :: cwork(:)

    write(stdout,'(/1x,a)') &
         '*------------------------------- DISENTANGLE --------------------------------*'

    allocate(eigval_opt(num_bands,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating eigval_opt in dis_main')
    eigval_opt=eigval


    ! Set up energy windows
    call dis_windows()


    ! Allocate arrays
    allocate(clamp(num_bands,num_bands,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating clamp in dis_main')
    allocate(cmtmp(num_bands,num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cmtmp in dis_main')
    allocate(caa(num_bands,num_bands,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating caa in dis_main')
    
    ! here we construct the U (unitarized projection) as described above
    write(stdout,'(/1x,a)') &
         '                  Unitarised projection of Wannier functions                  '
    write(stdout,'(1x,a)') &
         '                  ------------------------------------------                  '
    write(stdout,'(3x,a)') 'A_mn = <psi_m|g_n> --> S = A.A^+ --> U = S^-1/2.A'
    write(stdout,'(3x,a)',advance='no') 'In dis_project...' 
    call dis_project()
    write(stdout,'(a)') ' done'
 

    do nkp = 1, num_kpts  
       do j = 1, num_wann  
          do i = 1, ndimwin(nkp)  
             clamp(i,j,nkp) = u_matrix_opt(i,j,nkp)  
          enddo
       enddo
    enddo


    ! If there is an inner window, need to modify projection procedure
    ! (Sec. III.G SMV)
    if (linner) then
       write(stdout,'(3x,a)') 'Using an inner window (linner = T)'  
       write(stdout,'(3x,a)',advance='no') 'In dis_proj_froz...' 
       call dis_proj_froz()
       write(stdout,'(a)') ' done'
    else
       write(stdout,'(3x,a)') 'No inner window (linner = F)'         
    endif


    ! DEBUG
    ! check that the states in the columns of the final matrix clamp are ort
    ! at every k-point (i.e., that the matrix is unitary in the sense that
    ! conjg(lamp) . lamp = 1 - but not lamp . conjg(lamp) = 1 - )
    ! In particular, this checks whether the projected gaussians are indeed
    ! orthogonal to the frozen states, at those k-points where both are pres
    ! the trial subspace.
    do nkp = 1, num_kpts  
       do l = 1, num_wann  
          do m = 1, l  
             ctmp = cmplx_0
             do j = 1, ndimwin(nkp)  
                ctmp = ctmp + conjg(clamp(j,m,nkp)) * clamp(j,l,nkp)
             enddo
             if (l.eq.m) then  
                if (abs(ctmp - cmplx_1).gt.1.0e-8_dp) then  
                   write(stdout,*) ' *** ERROR ***'  
                   write(stdout,*) ' projected gaussians in clamp not orthonormal'
                   write(stdout,'(a11,i4)') ' at k-point ', nkp  
                   write(stdout,'(i2,2x,i2,f16.12,1x,f16.12)') l, m, ctmp  
                   call io_error('dis_main: orthonormal error 1') 
                endif
             else  
                if (abs(ctmp).gt.1.0e-8_dp) then  
                   write(stdout,*) ' *** ERROR ***'  
                   write(stdout,*) ' projected gaussians in clamp not orthonormal'
                   write(stdout,'(a11,i4)') ' at k-point ', nkp  
                   write(stdout,'(i2,2x,i2,f16.12,1x,f16.12)') l, m, ctmp  
                   call io_error('dis_main: orthonormal error 2') 
                endif
             endif
          enddo
       enddo
    enddo
    ! ENDDEBUG


!!$    write(stdout,'(1x,a)') &
!!$         '+----------------------------------------------------------------------------+'


    ! Slim down the original overlap matrices, removing rows and columns
    ! corresponding to u_nks that fall outside the outer energy window
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          nkp2 = nnlist(nkp,nn)  
          do i = 1, ndimwin(nkp)  
             do j = 1, ndimwin(nkp2)  
                m = nfirstwin(nkp) + i - 1  
                n = nfirstwin(nkp2) + j - 1  
                cmtmp(i,j) = m_matrix_orig(m,n,nn,nkp)  
             enddo
          enddo
          m_matrix_orig(:,:,nn,nkp) = cmplx_0
!          do i = 1, num_bands  
!             do j = 1, num_bands  
!                m_matrix_orig(i,j,nkp,nn) = cmplx_0
!             enddo
!          enddo
          do i = 1, ndimwin(nkp)  
             do j = 1, ndimwin(nkp2)  
                m_matrix_orig(i,j,nn,nkp) = cmtmp(i,j)  
             enddo
          enddo
       enddo
    enddo


    lwindow=.false.
    do nkp=1,num_kpts
       do j=nfirstwin(nkp),nfirstwin(nkp)+ndimwin(nkp)-1
          lwindow(j,nkp)=.true.
       end do
    end do



    ! Extract the optimally-connected num_wann-dimensional subspaces
    write(stdout,'(/1x,a)') &
         '                  Extraction of optimally-connected subspace                  '
    write(stdout,'(1x,a)') &
         '                  ------------------------------------------                  '
    call dis_extract()
    write(stdout,'(1x,a/)') &
         '+----------------------------------------------------------------------------+'


    ! Find the num_wann x num_wann overlap matrices between the basis states of
    ! the optimal subspaces
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          nkp2 = nnlist(nkp,nn)
          do i = 1, num_wann  
             do j = 1, num_wann  
                cmtmp(i,j) = cmplx_0
                do m = 1, ndimwin(nkp)  
                   do n = 1, ndimwin(nkp2)  
                      cmtmp(i,j) = cmtmp(i,j) + conjg(clamp(m,i,nkp)) &
                           * m_matrix_orig(m,n,nn,nkp) * clamp(n,j,nkp2)
                   enddo
                enddo
             enddo
          enddo
          m_matrix_orig(:,:,nn,nkp) = cmplx_0
!          do i = 1, num_bands  
!             do j = 1, num_bands  
!                m_matrix_orig(i,j,nkp,nn) = cmplx_0
!             enddo
!          enddo
          do i = 1, num_wann  
             do j = 1, num_wann  
                m_matrix_orig(j,i,nn,nkp) = cmtmp(j,i)  
             enddo
          enddo
       enddo
    enddo
    
    ! Allocate arrays needed for ZGESVD
    allocate(svals(num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating svals in dis_main')
    allocate(rwork(5*num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating rwork in dis_main')
    allocate(cv(num_bands,num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cv in dis_main')
    allocate(cz(num_bands,num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cz in dis_main')
    allocate(cwork(4*num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cwork in dis_main')
    
    ! Find the initial guess for the square unitary rotation matrix cu. The
    ! is similar to Sec. III.D of SMV, but with square instead of rectangula
    ! matrices:
    do nkp = 1, num_kpts  
       ! First find caa, the square overlap matrix <psitilde_nk|g_m>, where psi
       ! eigenstate of the optimal subspace. Note that, contrary to what is imp
       ! in Sec. III.E of SMV, this does *not* need to be computed by brute for
       ! instead we take advantage of the previous computation of overlaps with
       ! the same gaussians that was used to initiate the minimization of Omega
       ! note: |psi> clamp = |psitilde> and obviously <psitilde| = clamp^dagger
       do j = 1, num_wann  
          do i = 1, num_wann  
             caa(i,j,nkp) = cmplx_0
             do l = 1, ndimwin(nkp)  
                caa(i,j,nkp) = caa(i,j,nkp) &
                     + conjg(clamp(l,i,nkp)) * a_matrix(l,j,nkp)
             enddo
          enddo
       enddo
       
       ! Singular-value decomposition
       call ZGESVD ('A', 'A', num_wann, num_wann, caa(1,1,nkp), num_bands, &
            svals, cz, num_bands, cv, num_bands, cwork, 4*num_bands, rwork, info)
       if (info.ne.0) then  
          write(stdout,*) ' ERROR: IN ZGESVD IN dis_main'  
          write(stdout,*) 'K-POINT NKP=', nkp, ' INFO=', info  
          if (info.lt.0) then  
             write(stdout,*) 'THE ',  -info, '-TH ARGUMENT HAD ILLEGAL VALUE'  
          endif
          call io_error('dis_main: problem in ZGESVD 1')  
       endif
       
       ! cu is the initial guess for the unitary rotation of the basis states g
       ! out by the subroutine extract
       do j = 1, num_wann  
          do i = 1, num_wann  
             u_matrix(i,j,nkp) = cmplx_0
             do l = 1, num_wann  
                u_matrix(i,j,nkp) = u_matrix(i,j,nkp) + cz(i,l) * cv(l,j) 
             enddo
          enddo
       enddo
    enddo
    ! loop over k
    
    ! so now we have the U's that rotate the wavefunctions at each k-point.
    ! the matrix elements M_ij have also to be updated
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          nkp2 = nnlist(nkp,nn)  
          do i = 1, num_wann  
             do j = 1, num_wann  
                cmtmp(i,j) = cmplx_0
                do m = 1, num_wann  
                   do n = 1, num_wann  
                      cmtmp(i,j) = cmtmp(i,j) &
                           + conjg(u_matrix(m,i,nkp)) * u_matrix(n,j,nkp2) &
                           * m_matrix_orig(m,n,nn,nkp)
                   enddo
                enddo
             enddo
          enddo
          do i = 1, num_wann  
             do j = 1, num_wann  
                m_matrix(i,j,nn,nkp) = cmtmp(i,j)  
             enddo
          enddo
       enddo
    enddo

!    eigval_opt=0.0_dp
!    ! Copy the eigenvalues in the optimal space to eigval_opt (for plotting etc)
!    do nkp=1,num_kpts
!       do i=1,num_wann
!          eigval_opt(i,nkp)=eigval(i,nkp)
!       enddo
!    enddo

    u_matrix_opt=cmplx_0

    do nkp = 1, num_kpts  
       do j = 1, num_wann  
          do i = 1, ndimwin(nkp)  
             u_matrix_opt(i,j,nkp)  = clamp(i,j,nkp) 
          enddo
       enddo
    enddo




    ! Deallocate arrays for ZGESVD
    deallocate(cwork,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating cwork in dis_main')
    deallocate(cz,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating cz in dis_main')
    deallocate(cv,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating cv in dis_main')
    deallocate(rwork,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating rwork in dis_main')
    deallocate(svals,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating svals in dis_main')

    ! Deallocate module arrays allocated in dis_windows
    deallocate(lfrozen,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating lfrozen in dis_main')
    deallocate(indxnfroz,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating indxnfroz in dis_main')
    deallocate(indxfroz,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating indxfroz in dis_main')
    deallocate(ndimfroz,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating ndimfroz in dis_main')
    deallocate(nfirstwin,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating nfirstwin in dis_main')

    ! Deallocate arrays allocated in dis_main
    deallocate(caa,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating caa in dis_main')
    deallocate(cmtmp,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating cmtmp in dis_main')
    deallocate(clamp,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating clamp in dis_main')

    deallocate(eigval_opt,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating eigval_opt in dis_main')


    return

  end subroutine dis_main



  !==================================================================!
  subroutine dis_windows()
  !==================================================================!
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !==================================================================!  
    ! here, we select the states that are inside the outer window (i.e.
    ! the energy window out of which we fish out the optimally-connected
    ! subspace) and those that are inside the inner window (that make up
    ! the frozen manifold, and are straightfowardly included as they are).
    ! this, in practice, amounts to slimming down the original num_wann x num_wann
    ! overlap matrices, removing rows and columns that belong to u_nks that
    ! have been excluded forever, and marking (indexing) the rows and column
    ! that correspond to frozen states.

    ! note - in windows eigval_opt are shifted, so the lowest ones go
    ! from nfirstwin(nkp) to nfirstwin(nkp)+ndimwin(nkp)-1, and above
    ! they are set to zero

    implicit none

    ! internal variables
    integer :: i,j,nkp,ierr
    integer :: imin,imax,kifroz_min,kifroz_max

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


    ! Allocate module arrays
    allocate(nfirstwin(num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating nfirstwin in dis_windows')
    allocate(ndimfroz(num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ndimfroz in dis_windows')
    allocate(indxfroz(num_bands,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating indxfroz in dis_windows')
    allocate(indxnfroz(num_bands,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating indxnfroz in dis_windows')
    allocate(lfrozen(num_bands,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating lfrozen in dis_windows')

    linner = .false.  

    write(stdout,'(1x,a)') &
         '+----------------------------------------------------------------------------+'
    write(stdout,'(1x,a)') &
         '|                              Energy  Windows                               |'
    write(stdout,'(1x,a)') &
         '|                              ---------------                               |'
    write(stdout,'(1x,a,f10.5,a,f10.5,a)') &
         '|                   Outer: ',dis_win_min,'  to ',dis_win_max,&
         '  (eV)                   |'
    if (frozen_states) then
       write(stdout,'(1x,a,f10.5,a,f10.5,a)') &
            '|                   Inner: ',dis_froz_min,'  to ',dis_froz_max,&
            '  (eV)                   |'
    else
       write(stdout,'(1x,a)') &
            '|                   No frozen states were specified                          |'
    endif
    write(stdout,'(1x,a)') &
         '+----------------------------------------------------------------------------+'

    do nkp = 1, num_kpts  
       ! Check which eigenvalues fall within the outer window
       if ( (eigval_opt(1,nkp).gt.dis_win_max).or.&
            (eigval_opt(num_bands,nkp).lt.dis_win_min) ) then
          write(stdout,*) ' ERROR AT K-POINT: ', nkp  
          write(stdout,*) ' ENERGY WINDOW (eV):    [',dis_win_min,  ',', dis_win_max,     ']'
          write(stdout,*) ' EIGENVALUE RANGE (eV): [',eigval_opt(1,nkp),',',eigval_opt(num_bands,nkp),']'
          call io_error('DISENTANGLE: The outer energy window contains no eigenvalues')
       endif


       ! Note: we assume that eigvals are ordered from the bottom up
       imin = 0  
       do i = 1, num_bands
          if (imin.eq.0) then  
             if ( (eigval_opt(i,nkp).ge.dis_win_min).and.&
                  (eigval_opt(i,nkp).le.dis_win_max) ) imin = i
             imax = i  
          endif
          if (eigval_opt(i,nkp).le.dis_win_max) imax = i  
       enddo

       ndimwin(nkp) = imax - imin + 1  

       nfirstwin(nkp) = imin  

       if (ndimwin(nkp).lt.num_wann) then  
          write(stdout,483) 'Error at k-point ',nkp,' ndimwin=',ndimwin(nkp),' num_wann=',num_wann
483       format(1x,a17,i4,a8,i3,a9,i3)  
          call io_error('DISENTANGLE: Energy window contains fewer states than number of target WFs') 
       endif


       do i = 1, ndimwin(nkp)  
          lfrozen(i,nkp) = .false.  
       enddo

       ! Check which eigenvalues fall within the inner window
       kifroz_min = 0  
       kifroz_max = -1

       ! (note that the above obeys kifroz_max-kifroz_min+1=kdimfroz=0, as we w
       if (frozen_states) then
          do i = imin, imax  
             if (kifroz_min.eq.0) then 
                if ( (eigval_opt(i,nkp).ge.dis_froz_min).and.&
                     (eigval_opt(i,nkp).le.dis_froz_max) ) then
                   ! Relative to bottom of outer window
                   kifroz_min = i - imin + 1  
                   kifroz_max = i - imin + 1  
                endif
             elseif (eigval_opt(i,nkp).le.dis_froz_max) then  
                kifroz_max = kifroz_max + 1  
                ! DEBUG
                ! if(kifroz_max.ne.i-imin+1) stop 'something wrong...'
                ! ENDDEBUG
             endif
          enddo
       endif

       ndimfroz(nkp) = kifroz_max - kifroz_min + 1  

       if (ndimfroz(nkp).gt.num_wann) then  
          write(stdout,401) nkp, ndimfroz(nkp),num_wann  
401       format(' ERROR AT K-POINT ',i4,' THERE ARE ',i2, &
               ' BANDS INSIDE THE INNER WINDOW AND ONLY',i2, &
               ' TARGET BANDS')
          write(stdout,402) (eigval_opt(i,nkp),i = imin, imax)  
402       format('BANDS: (eV)',10(F10.5,1X))  
          call io_error('DISENTANGLE: More states in the frozen window than target WFs')
       endif

       if (ndimfroz(nkp).gt.0) linner = .true.  
       ! DEBUG
       !         write(*,'(a,i4,a,i2,a,i2)') 'k point ',nkp,
       !     &    ' lowest band in outer win is # ',imin,
       !     &    '   # frozen states is ',ndimfroz(nkp)
       ! ENDDEBUG
       ! Generate index array for frozen states (those inside inner window)
       if (ndimfroz(nkp).gt.0) then  
          do i = 1, ndimfroz(nkp)  
             indxfroz(i,nkp) = kifroz_min + i - 1  
             lfrozen(indxfroz(i,nkp),nkp) = .true.  
          enddo
          if (indxfroz(ndimfroz(nkp),nkp).ne.kifroz_max) then  
             write(stdout,*) ' Error at k-point ', nkp, ' frozen band #', i  
             write(stdout,*) ' ndimfroz=', ndimfroz(nkp)  
             write(stdout,*) ' kifroz_min=', kifroz_min  
             write(stdout,*) ' kifroz_max=', kifroz_max  
             write(stdout,*) ' indxfroz(i,nkp)=', indxfroz(i,nkp)  
             call io_error('DISENTANGLE: Something fishy...')
          endif
       endif

       ! Generate index array for non-frozen states
       i = 0  
       do j = 1, ndimwin(nkp)  
          !           if (lfrozen(j,nkp).eqv..false.) then  
          if (.not.lfrozen(j,nkp)) then  
             i = i + 1  
             indxnfroz(i,nkp) = j  
          endif
       enddo

       if ( i.ne.ndimwin(nkp) - ndimfroz(nkp) ) then  
          write(stdout,*) ' Error at k-point: ',nkp
          write(stdout,'(3(a,i5))') ' i: ',i,' ndimwin: ',ndimwin(nkp),&
               ' ndimfroz: ',ndimfroz(nkp)
          call io_error('DISENTANGLE: i .ne. (ndimwin-ndimfroz) at k-point')
       endif

       ! Slim down eigval vector at present k
       do i = 1, ndimwin(nkp)  
          j = nfirstwin(nkp) + i - 1  
          eigval_opt(i,nkp) = eigval_opt(j,nkp)  
       enddo

       do i = ndimwin(nkp) + 1, num_bands
          eigval_opt(i,nkp) = 0.0_dp
       enddo

    enddo
    ! [k-point loop (nkp)]

    if (iprint>1) then
       write(stdout,'(1x,a)') &
            '|                        K-points with Frozen States                         |'
       write(stdout,'(1x,a)') &
            '|                        ---------------------------                         |'
       i=0
       do nkp=1,num_kpts
          if (ndimfroz(nkp).gt.0) then
             i=i+1
             if (i.eq.1) then
                write(stdout,'(1x,a,i6)',advance='no') '|',nkp
             else if ((i.gt.1) .and. (i.lt.12)) then
                write(stdout,'(i6)',advance='no') nkp
             else if (i.eq.12) then 
                write(stdout,'(i6,a)') nkp,'    |'
                i=0
             endif
          endif
       enddo
       if (i.ne.0) then
          do j=1,12-i
             write(stdout,'(6x)',advance='no')
          enddo
          write(stdout,'(a)') '    |'
       endif
       write(stdout,'(1x,a)') &
            '+----------------------------------------------------------------------------+'
    endif

    write(stdout,'(3x,a,i4)') 'Number of target bands to extract: ',num_wann
!!$    write(stdout,'(1x,a,76("-"),a)')'+','+'
    if (iprint>1) then
       write(stdout,'(1x,a)') &
            '+----------------------------------------------------------------------------+'
       write(stdout,'(1x,a)') &
            '|                                  Windows                                   |'
       write(stdout,'(1x,a)') &
            '|                                  -------                                   |'
       write(stdout,'(1x,a)') &
            '|               K-point      Ndimwin     Ndimfroz    Nfirstwin               |'
       write(stdout,'(1x,a)') &
            '|               ----------------------------------------------               |'
       do nkp=1,num_kpts
          write(stdout,403) nkp,ndimwin(nkp),ndimfroz(nkp),nfirstwin(nkp)
       enddo
403    format(1x,'|',14x,i6,7x,i6,7x,i6,6x,i6,18x,'|')
       write(stdout,'(1x,a)') &
            '+----------------------------------------------------------------------------+'
!!$    write(stdout,'(1x,a,76("-"),a)')'+','+'
    endif

!!$    do nkp=1,num_kpts
!!$       write(stdout,'(a,i5,a,i3,a,i3,a,i3)') '  K-point ',nkp,&
!!$            ' ndimwin: ',ndimwin(nkp),' ndimfroz: ',ndimfroz(nkp),&
!!$            ' nfirstwin: ',nfirstwin(nkp)
!!$    enddo
!!$    write(stdout,*) ' '  

    return  

  end subroutine dis_windows


  !==================================================================!
  subroutine dis_project()
  !==================================================================!
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !==================================================================!  


    implicit none


    ! internal variables
    integer :: i,j,l,m,nkp,info,ierr
    real(kind=dp),    allocatable :: svals(:)
    real(kind=dp),    allocatable :: rwork(:) 
    complex(kind=dp)              :: ctmp2
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:,:)
    complex(kind=dp), allocatable :: cvdag(:,:)
    complex(kind=dp), allocatable :: catmpmat(:,:,:)


    allocate(catmpmat(num_bands,num_bands,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating catmpmat in dis_project')
    allocate(svals(num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating svals in dis_project')
    allocate(rwork(5*num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating rwork in dis_project')
    allocate(cvdag(num_bands,num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cvdag in dis_project')
    allocate(cz(num_bands,num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cz in dis_project')
    allocate(cwork(4*num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cwork in dis_project')


    ! here we slim down the ca matrix
    ! up to here num_bands(=num_bands) X num_wann(=num_wann)
    do nkp = 1, num_kpts
       do i = 1, ndimwin(nkp)
          do j = 1, num_wann
             catmpmat(i,j,nkp) = a_matrix(nfirstwin(nkp)+i-1,j,nkp)  
          enddo
       enddo
       do i = 1, ndimwin(nkp)  
          do j = 1, num_wann  
             a_matrix(i,j,nkp) = catmpmat(i,j,nkp)  
          enddo
       enddo
       do i = ndimwin(nkp) + 1, num_bands  
          do j = 1, num_wann  
             a_matrix(i,j,nkp) = cmplx_0  
          enddo
       enddo
    enddo

    ! original notes from Nicola (refer only to the square case) -----------
    !
    ! now I calculate the transformation matrix CU = CS^(-1/2).CA,
    ! where CS = CA.CA^\dagger. CS is diagonalized with a Schur
    ! factorization, to be on the safe side of numerical stability.
    ! From LAPACK:
    ! ZGEES computes for an N-by-N complex nonsymmetric matrix Y, the eigen-
    ! values, the Schur form T, and, optionally, the matrix of Schur vectors
    ! This gives the Schur factorization Y = Z*T*(Z**H).
    !
    ! Optionally, it also orders the eigenvalues on the diagonal of the Schu
    ! form so that selected eigenvalues are at the top left.  The leading co
    ! of Z then form an orthonormal basis for the invariant subspace corresp
    ! ing to the selected eigenvalues.
    !
    ! A complex matrix is in Schur form if it is upper triangular.
    ! NOTES FROM IVO DISENTANGLING (e.g. NON-SQUARE) PROJECTION
    ! (See Sec. III.D of SMV paper)
    ! COMPUTE THE ndimwin(K) x num_wann MATRIX CU THAT YIELDS, FROM THE ndimwin
    ! ORIGINAL BLOCH STATES, THE num_wann BLOCH-LIKE STATES WITH MAXIMAL PROJEC
    ! ONTO THE num_wann GAUSSIANS:
    !
    ! CU = CA.CS^{-1/2}, CS = transpose(CA).CA
    !
    ! USE THE SINGULAR-VALUE DECOMPOSITION OF THE MATRIX CA:
    !
    ! CA = CZ.CD.CVdagger (note: zgesvd spits out CVdagger)
    !
    ! WHICH YIELDS
    !
    ! CU = CZ.CD.CD^{-1}.CVdag
    !
    ! WHERE CZ IS ndimwin(NKP) x ndimwin(NKP) AND UNITARY, CD IS
    ! ndimwin(NKP) x num_wann AND DIAGONAL, CD^{-1} IS num_wann x num_wann AND
    ! DIAGONAL, AND CVdag IS num_wann x num_wann AND UNITARY.

    do nkp = 1, num_kpts  
       !
       ! SINGULAR VALUE DECOMPOSITION
       !
       call ZGESVD('A', 'A', ndimwin(nkp), num_wann, a_matrix(1,1,nkp), &
            num_bands, svals, cz, num_bands, cvdag, num_bands, cwork, &
            4*num_bands, rwork, info)
       if (info.ne.0) then  
          write(stdout,*) ' ERROR: IN ZGESVD IN dis_project'  
          write(stdout,*) ' K-POINT NKP=', nkp, ' INFO=', info  
          if (info.lt.0) then  
             write(stdout,*) ' THE ',  -info, '-TH ARGUMENT HAD ILLEGAL VALUE'  
          endif
          call io_error('DISENTANGLE: Error in ZGESVD in dis_project')
          call io_error('dis_project: problem in ZGESVD 1')   
       endif

       !
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
       ! it seemed to be identical to the input ca (as it should be)

       do j = 1, num_wann  
          do i = 1, ndimwin(nkp)  
             u_matrix_opt(i,j,nkp) = cmplx_0
             a_matrix(i,j,nkp) = cmplx_0
             do l = 1, num_wann  
                u_matrix_opt(i,j,nkp) = u_matrix_opt(i,j,nkp) + cz(i,l) * cvdag(l,j)  
                a_matrix(i,j,nkp) = a_matrix(i,j,nkp) + cz(i,l) * svals(l) * cvdag(l,j)  
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
                ctmp2 = ctmp2 + u_matrix_opt(m,j,nkp) * conjg(u_matrix_opt(m,i,nkp))  
             enddo
             if ( (i.eq.j).and.(abs(ctmp2-cmplx_1).gt.0.00001_dp) ) then
                write(stdout,*) ' ERROR: unitarity of initial U'  
                write(stdout,'(1x,a,i2)') 'nkp= ', nkp  
                write(stdout,'(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j  
                write(stdout,'(1x,a,f12.6,1x,f12.6)') &
                     '[u_matrix_opt.transpose(u_matrix_opt)]_ij= ',&
                     real(ctmp2,dp),aimag(ctmp2)
                call io_error('DISENTANGLE: Error in unitarity of initial U in dis_project')
             endif
             if ( (i.ne.j) .and. (abs(ctmp2).gt.0.00001_dp) ) then  
                write(stdout,*) ' ERROR: unitarity of initial U'  
                write(stdout,'(1x,a,i2)') 'nkp= ', nkp  
                write(stdout,'(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j  
                write(stdout,'(1x,a,f12.6,1x,f12.6)') &
                     '[u_matrix_opt.transpose(u_matrix_opt)]_ij= ', &
                     real(ctmp2,dp),aimag(ctmp2)
                call io_error('DISENTANGLE: Error in unitarity of initial U in dis_project')
             endif
          enddo
       enddo
    enddo
    ! NKP

    deallocate(cwork,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cwork in dis_project')
    deallocate(cz,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cz in dis_project')
    deallocate(cvdag,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cvdag in dis_project')
    deallocate(rwork,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rwork in dis_project')
    deallocate(svals,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating svals in dis_project')
    deallocate(catmpmat,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating catmpmat in dis_project')

    return  

  end subroutine dis_project



  !==================================================================!
  subroutine dis_proj_froz()
  !==================================================================!
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !                                                                  !
  !==================================================================!  

      implicit none

      ! COMPUTES THE LEADING EIGENVECTORS OF Q_froz . P_s . Q_froz, WHERE P_s
      ! PROJECTOR OPERATOR ONTO THE SUBSPACE S OF THE PROJECTED GAUSSIANS, P_f
      ! THE PROJECTOR ONTO THE FROZEN STATES, AND Q_froz = 1 - P_froz, ALL EXP
      ! IN THE BASIS OF THE BLOCH EIGENSTATES INSIDE THE OUTER ENERGY WINDOW
      ! (See Eq. (27) in Sec. III.G of SMV)
      ! INPUT: num_wann,ndimwin,ndimfroz,indxfroz,lfrozen
      ! MODIFIED: clamp (At input it contains the gaussians projected onto the
      !                 window states in the routine project.f. At output the
      !                 entries with the second index from 1 to ndimfroz(nkp) contain
      !                 the frozen (inner window) states, while those from
      !                 ndimfroz(nkp)+1 to num_wann have been replaced by the new trial
      !                 states outside the inner window.)

      ! *********************************************************
      ! VARIABLES USED BY LAPACK'S ZHPEVX DIAGONALIZATION ROUTINE
      ! *********************************************************
      
      integer,          allocatable :: iwork(:)
      integer,          allocatable :: ifail(:)
      real(kind=dp),    allocatable :: w(:)
      real(kind=dp),    allocatable :: rwork(:)
      complex(kind=dp), allocatable :: cap(:)
      complex(kind=dp), allocatable :: cwork(:)
      complex(kind=dp), allocatable :: cz(:,:)

      ! *********
      ! INTERNAL:
      ! *********
      !
      ! CP_S(M,N)          PROJECTION OPERATOR ONTO THE SUBSPACE OF THE PROJEC
      !                   GAUSSIANS, EXPRESSED IN THE BASIS OF THE ORIGINAL BL
      !                   EIGENSTATES INSIDE THE OUTER WINDOW (FOR THE PRESENT
      !                   K-POINT)
      ! CQ_FROZ(M,N)       PROJECTION OPERATOR ONTO THE SUBSPACE OF THE STATES
      !                   THE SPACE OF FROZEN STATES (BUT INSIDE THE OUTER WIN
      !                   EXPRESSED IN THE BASIS OF THE ORIGINAL BLOCH EIGENST
      !                   INSIDE THE OUTER WINDOW (FOR THE PRESENT K-POINT)
      ! CPQ(M,N)           THE MATRIX cp_s . cq_froz FOR THE PRESENT K-POINT
      ! CQPQ(M,N)          THE MATRIX cq_froz . cp_s . cq_froz FOR THE PRESENT
      !

      integer :: il,iu,nkp,l,j,n,m,info,ierr
      complex(kind=dp), allocatable :: cp_s(:,:)
      complex(kind=dp), allocatable :: cq_froz(:,:)
      complex(kind=dp), allocatable :: cpq(:,:)
      complex(kind=dp), allocatable :: cqpq(:,:)


      allocate(iwork(5*num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating iwork in dis_proj_froz')
      allocate(ifail(num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating ifail in dis_proj_froz')
      allocate(w(num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating w in dis_proj_froz')
      allocate(rwork(7*num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating rwork in dis_proj_froz')
      allocate(cap((num_bands*(num_bands+1))/2),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cap in dis_proj_froz')
      allocate(cwork(2*num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cwork in dis_proj_froz')
      allocate(cz(num_bands,num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cz in dis_proj_froz')
      
      allocate(cp_s(num_bands,num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cp_s in dis_proj_froz')
      allocate(cq_froz(num_bands,num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cq_froz in dis_proj_froz')
      allocate(cpq(num_bands,num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cpq in dis_proj_froz')
      allocate(cqpq(num_bands,num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cqpq in dis_proj_froz')

      do nkp = 1, num_kpts  
         ! Put the frozen states in the lowest columns of clamp
         if (ndimfroz(nkp).gt.0) then  
            do l = 1, ndimfroz(nkp)  
               do j = 1, ndimwin(nkp)  
                  clamp(j,l,nkp) = cmplx_0  
               enddo
               clamp(indxfroz(l,nkp),l,nkp) = cmplx_1
            enddo
         endif

         ! If there are non-frozen states, compute the num_wann-ndimfroz(nkp) leadin
         ! eigenvectors of cqpq
         if (num_wann.gt.ndimfroz(nkp)) then  
            do n = 1, ndimwin(nkp)  
               do m = 1, ndimwin(nkp)  
                  cq_froz(m,n) = cmplx_0  
                  cp_s(m,n) = cmplx_0  
                  do l = 1, num_wann  
                     cp_s(m,n) = cp_s(m,n) + clamp(m,l,nkp) * conjg(clamp(n,l,nkp))
                  enddo
               enddo
               if (.not.lfrozen(n,nkp)) cq_froz(n,n) = cmplx_1
            enddo
            do n = 1, ndimwin(nkp)  
               do m = 1, ndimwin(nkp)  
                  cpq(m,n) = cmplx_0  
                  do l = 1, ndimwin(nkp)  
                     cpq(m,n) = cpq(m,n) + cp_s(m,l) * cq_froz(l,n)  
                  enddo
               enddo
            enddo
            do n = 1, ndimwin(nkp)  
               do m = 1, ndimwin(nkp)  
                  cqpq(m,n) = cmplx_0  
                  do l = 1, ndimwin(nkp)  
                     cqpq(m,n) = cqpq(m,n) + cq_froz(m,l) * cpq(l,n)  
                  enddo
               enddo
            enddo

            ! DEBUG
            ! check hermiticity of cqpq
            do n = 1, ndimwin(nkp)  
               do m = 1, n  
                  if (abs(cqpq(m,n) - conjg(cqpq(n,m))).gt.1.0e-8_dp) then
                     write(stdout,*) ' matrix CQPQ is not hermitian'  
                     write(stdout,*) ' k-point ', nkp  
                     call io_error('dis_proj_froz: error')  
                  endif
               enddo
            enddo

            ! ENDDEBUG
            do n = 1, ndimwin(nkp)  
               do m = 1, n  
                  cap(m + (n - 1) * n / 2) = cqpq(m,n)  
               enddo
            enddo
            il = ndimwin(nkp) - (num_wann - ndimfroz(nkp) ) + 1  
            iu = ndimwin(nkp)  
            call ZHPEVX ('V', 'A', 'U', ndimwin(nkp), cap, 0.0_dp, 0.0_dp, il, &
                 iu, -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)
            ! DEBUG
            if (info.lt.0) then  
               write(stdout,*) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING CQPQ MATRIX'
               write(stdout,*) ' THE ',  -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
               call io_error('dis_proj_frozen: error')  
            elseif (info.gt.0) then  
               write(stdout,*) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING CQPQ MATRIX'
               write(stdout,*) info, 'EIGENVECTORS FAILED TO CONVERGE'  
               call io_error('dis_proj_frozen: error') 
            endif
            ! ENDDEBUG
            ! DEBUG
            if (m.ne.ndimwin(nkp)) then  
               write(stdout,*) ' *** ERROR *** in dis_proj_froz'  
               write(stdout,*) ' Number of eigenvalues/vectors obtained is', &
                    m, ' not equal to the number asked,', ndimwin(nkp)
               call io_error('dis_proj_frozen: error')  
            endif
            ! ENDDEBUG
            ! DEBUG
            ! check that the eigenvalues are between 0 and 1
            if (iprint>2) then
               write(stdout,'(/a,i3,a,i3,a,i3,a)') ' K-point ', nkp, ' ndimwin: ', &
                    ndimwin(nkp),' we want the',num_wann - ndimfroz(nkp),&
                    ' leading eigenvector(s) of QPQ'
            endif
            do j = 1, ndimwin(nkp)  
               if (iprint>2) write(stdout,'(a,i3,a,f16.12)') '  lambda(', j, ')=', w(j)  
               if ( (w(j).lt.-1.0e-8_dp).or.(w(j).gt.1.0d0 + 1.0e-8_dp) ) then
                  call io_error('dis_proj_frozen: error - Eigenvalues not between 0 and 1') 
               endif
            enddo
            ! ENDDEBUG

            ! PICK THE num_wann-nDIMFROZ(NKP) LEADING EIGENVECTORS AS TRIAL STATES
            ! and PUT THEM RIGHT AFTER THE FROZEN STATES IN cLAMP
            do l = ndimfroz(nkp) + 1, num_wann  
               do j = 1, ndimwin(nkp)  
                  clamp(j,l,nkp) = cz(j,il)  
               enddo
               il = il + 1  
            enddo

            ! DEBUG
            if (il - 1.ne.iu) then  
               call io_error('dis_proj_frozen: error -  il-1.ne.iu')
            endif
            ! ENDDEBUG

            ! num_wann>nDIMFROZ(NKP)

         endif
         ! NKP

      enddo

      deallocate(cqpq,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cqpq in dis_proj_froz')
      deallocate(cpq,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cpq in dis_proj_froz')
      deallocate(cq_froz,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cq_froz in dis_proj_froz')
      deallocate(cp_s,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cp_s in dis_proj_froz')

      deallocate(cz,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cz in dis_proj_froz')
      deallocate(cwork,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cwork in dis_proj_froz')
      deallocate(cap,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cap in dis_proj_froz')
      deallocate(rwork,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating rwork in dis_proj_froz')
      deallocate(w,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating w in dis_proj_froz')
      deallocate(ifail,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating ifail in dis_proj_froz')
      deallocate(iwork,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating iwork in dis_proj_froz')

      return  

    end subroutine dis_proj_froz



    !==================================================================!
    subroutine dis_extract()
    !==================================================================!
    !                                                                  !
    !                                                                  !
    !                                                                  !
    !                                                                  !
    !==================================================================!  
      
      use io, only: io_time

      implicit none

      ! Extracts an num_wann-dimensional subspace at each k by minimizing Omega_I
      ! INPUT: num_wann,cm,ndimwin,ndimfroz,indxnfroz,nnlist,nshells,nntot,
      ! nwhich,wb,alphafixe,nitere,eigval
      ! MODIFIED:
      !
      !           clamp (At input it contains the initial guess for the optima
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


      ! Internal variables
      integer :: i,j,l,m,nkp,info,ierr
      integer :: icompflag,iter,nnx,ndnnx,ndnn,nnsh,ndiff,k_pls_b
      real(kind=dp) :: womegai,wkomegai,womegai1,rlambda
      real(kind=dp), allocatable :: wkomegai1(:)
      complex(kind=dp), allocatable :: ceamp(:,:,:)
      complex(kind=dp), allocatable :: camp(:,:,:)
      complex(kind=dp), allocatable :: cmk(:,:,:)
      complex(kind=dp), allocatable :: czmat_in(:,:,:)
      complex(kind=dp), allocatable :: czmat_out(:,:,:)
      complex(kind=dp), allocatable :: cham(:,:,:)

      integer,          allocatable :: iwork(:)
      integer,          allocatable :: ifail(:)
      real(kind=dp),    allocatable :: w(:)
      real(kind=dp),    allocatable :: rwork(:)
      complex(kind=dp), allocatable :: cap(:)
      complex(kind=dp), allocatable :: cwork(:)
      complex(kind=dp), allocatable :: cz(:,:)

      allocate(iwork(5*num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating iwork in dis_extract')
      allocate(ifail(num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating ifail in dis_extract')
      allocate(w(num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating w in dis_extract')
      allocate(rwork(7*num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating rwork in dis_extract')
      allocate(cap((num_bands*(num_bands+1))/2),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cap in dis_extract')
      allocate(cwork(2*num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cwork in dis_extract')
      allocate(cz(num_bands,num_bands),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cz in dis_extract')

      allocate(wkomegai1(num_kpts),stat=ierr)
      if (ierr/=0) call io_error('Error allocating wkomegai1 in dis_extract')
      allocate(ceamp(num_bands,num_bands,num_kpts),stat=ierr)
      if (ierr/=0) call io_error('Error allocating ceamp in dis_extract')
      allocate(camp(num_bands,num_bands,num_kpts),stat=ierr)
      if (ierr/=0) call io_error('Error allocating camp in dis_extract')
      allocate(cmk(num_bands,num_bands,nntot),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cmk in dis_extract')
      allocate(czmat_in(num_bands,num_bands,num_kpts),stat=ierr)
      if (ierr/=0) call io_error('Error allocating czmat_in in dis_extract')
      allocate(czmat_out(num_bands,num_bands,num_kpts),stat=ierr)
      if (ierr/=0) call io_error('Error allocating czmat_out in dis_extract')
      allocate(cham(num_bands,num_bands,num_kpts),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cham in dis_extract')


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
      ! CLAMP(J,L,NKP)    AMPLITUDE OF THE J-TH ENERGY EIGENVECTOR INSIDE THE
      !                   ENERGY WINDOW AT THE NKP-TH K-POINT IN THE EXPANSION
      !                   THE L-TH LEADING RLAMBDA EIGENVECTOR AT THE SAME K-P
      !                   If there are M_k frozen states, they occupy the lowe
      !                   entries of the second index of clamp, and the leadin
      !                   nabnds-M_k eigenvectors of the Z matrix occupy the
      !                   remaining slots
      ! CAMP(J,L,NKP)     SAME AS CLAMP, BUT FOR THE COMPLEMENT SUBSPACE INSID
      !                   ENERGY WINDOW (I.E., THE NON-LEADING RLAMBDA EIGENVE
      ! CEAMP(J,L,NKPTS)  SAME AS CLAMP, BUT INSTEAD OF RLAMBDA EIGENVECTOR, I
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
      if (iprint>2) then
         write(stdout,'(a,/)') '  Original eigenvalues inside outer window:'  
         do nkp = 1, num_kpts  
            write(stdout,'(a,i3,3x,20(f9.5,1x))') '  K-point ', nkp,&
                 ( eigval_opt(i, nkp), i = 1, ndimwin (nkp) )
         enddo
      endif
      ! ENDDEBUG

!!$      write(stdout,121) num_wann  
!!$121   format(/,'  DIMENSION OF SUBSPACE: num_wann=',I6)  
!!$      write(stdout,122) dis_mix_ratio  
!!$122   format(/,'  MIXING PARAMETER: DIS_MIX_RATIO =',F6.3,/)  

      ! TO DO: Check if this is the best place to initialize icompflag
      icompflag = 0  

      write(stdout,'(1x,a)') &
           '+---------------------------------------------------------------------+<-- DIS'
      write(stdout,'(1x,a)') &
           '|  Iter     Omega_I(i-1)      Omega_I(i)      Delta (%)        Time   |<-- DIS'
      write(stdout,'(1x,a)') &
           '+---------------------------------------------------------------------+<-- DIS'

      ! ------------------
      ! BIG ITERATION LOOP
      ! ------------------
      do iter = 1, dis_num_iter  

         if (iter.eq.1) then  
            ! Initialize Z matrix at k points w/ non-frozen states
            do nkp = 1, num_kpts  
               if ( num_wann.gt.ndimfroz(nkp) ) then  
                  do nnx = 1, nntot  
                     do j = 1, ndimwin(nnlist(nkp,nnx))  
                        do i = 1, ndimwin(nkp)  
                          cmk(i,j,nnx) = m_matrix_orig(i,j,nnx,nkp)  
                        enddo
                     enddo
                  enddo
                  call dis_zmatrix(nkp,cmk,czmat_in(1,1,nkp))
               endif
            enddo
         else  
            ! [iter.ne.1]
            ! Update Z matrix at k points with non-frozen states, using a mixing sch
            do nkp = 1, num_kpts  
               if (num_wann.gt.ndimfroz(nkp) ) then  
                  do i = 1, ndimwin(nkp) - ndimfroz(nkp)  
                     do j = 1, i  
                        czmat_in(j,i,nkp) = dis_mix_ratio * czmat_out(j,i,nkp) &
                             + (1.0_dp - dis_mix_ratio) * czmat_in(j,i,nkp)
                        ! hermiticity
                        czmat_in(i,j,nkp) = conjg(czmat_in(j,i,nkp))  
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

         do nkp = 1, num_kpts  
            wkomegai1(nkp) = real(num_wann,dp) * wbtot  
            if ( ndimfroz(nkp).gt.0 ) then  
               do nnx = 1, nntot  
                     do j = 1, ndimwin(nnlist(nkp,nnx))  
                        do i = 1, ndimwin(nkp)  
                           cmk(i,j,nnx) = m_matrix_orig(i,j,nnx,nkp)  
                        enddo
                     enddo
               enddo
               do m = 1, ndimfroz(nkp)  
                  rlambda = dis_zeig(nkp,m,cmk) 
                  wkomegai1(nkp) = wkomegai1(nkp) - rlambda  
               enddo
            endif
         enddo


         ! Refine optimal subspace at k points w/ non-frozen states
         do nkp = 1, num_kpts  
            if ( num_wann.gt.ndimfroz(nkp) ) then  
               ! Diagonalize Z matrix
               do j = 1, ndimwin(nkp) - ndimfroz(nkp)  
                  do i = 1, j  
                     cap(i + ( (j - 1) * j) / 2) = czmat_in(i,j,nkp)  
                  enddo
               enddo
               ndiff = ndimwin(nkp) - ndimfroz(nkp)  
               call ZHPEVX ('V', 'A', 'U', ndiff, cap, 0.0_dp, 0.0_dp, 0, 0, &
                    -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)
               if (info.lt.0) then  
                  write(stdout,*) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING Z MATRIX'
                  write(stdout,*) ' THE ',  -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
                  call io_error(' dis_extract: error')  
               endif
               if (info.gt.0) then  
                  write(stdout,*) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING Z MATRIX'
                  write(stdout,*) info, ' EIGENVECTORS FAILED TO CONVERGE'  
                  call io_error(' dis_extract: error')  
               endif

               ! Update the optimal subspace by incorporating the num_wann-ndimfroz(nkp) l
               ! eigenvectors of the Z matrix into clamp. Also, add contribution from
               ! non-frozen states to wkomegai1(nkp) (minus the corresponding eigenvalu
               m = ndimfroz(nkp)  
               do j = ndimwin(nkp) - num_wann + 1, ndimwin(nkp) - ndimfroz(nkp)
                  m = m + 1  
                  wkomegai1(nkp) = wkomegai1(nkp) - w(j)  
                  do i = 1, ndimwin(nkp)  
                     clamp(i,m,nkp) = cmplx_0  
                  enddo
                  do i = 1, ndimwin(nkp) - ndimfroz(nkp)  
                     clamp(indxnfroz(i,nkp),m,nkp) = cz(i,j)  
                  enddo
               enddo
            endif
            ! [if num_wann>ndimfroz(nkp)]

            ! Now that we have contribs. from both frozen and non-frozen states to
            ! wkomegai1(nkp), add it to womegai1
            womegai1 = womegai1 + wkomegai1(nkp)  

            ! AT THE LAST ITERATION FIND A BASIS FOR THE (NDIMWIN(NKP)-num_wann)-DIMENS
            ! COMPLEMENT SPACE
            if (iter.eq.dis_num_iter) then  
               if (ndimwin(nkp).gt.num_wann) then  
                  do j = 1, ndimwin(nkp) - num_wann  
                     if ( num_wann.gt.ndimfroz(nkp) ) then  
                        ! USE THE NON-LEADING EIGENVECTORS OF THE Z-MATRIX
                        do i = 1, ndimwin(nkp)  
                           camp(i,j,nkp) = cz(i,j)  
                        enddo
                     else  
                        ! Then num_wann=NDIMFROZ(NKP)
                        ! USE THE ORIGINAL NON-FROZEN BLOCH EIGENSTATES
                        do i = 1,ndimwin(nkp)  
                           camp(i,j,nkp) = cmplx_0  
                           if (i.eq.indxnfroz(j,nkp)) camp(i,j,nkp) = cmplx_1
                        enddo
                     endif
                  enddo
               else  
                  icompflag = 1
               endif
            endif

         enddo
         ! [Loop over k points (nkp)]


         womegai1 = womegai1 / real(num_kpts,dp)  

         ! DEBUG
         ! Orthonormality check
         !         do nkp=1,nkpts
         !           write(*,*) ' '
         !           write(*,'(a8,i4)') 'k-point ',nkp
         !           do l=1,num_wann
         !           do m=1,l
         !             ctmp=czero
         !             do j=1,ndimwin(nkp)
         !               ctmp=ctmp+conjg(clamp(j,m,nkp))*clamp(j,l,nkp)
         !             enddo
         !             write(*,'(i2,2x,i2,f16.12,1x,f16.12)') l,m,ctmp
         !             if(l.eq.m) then
         !               if(abs(ctmp-cmplx(1.0d0,0.0d0)).gt.1.0e-8) then
         !                 write(*,'(a49,i4)')
         !     1           '*** ERROR *** with iterative subspace at k-point ',
         !     2           nkp
         !                 write(*,*) 'vectors in clamp not orthonormal'
         !                 stop
         !               endif
         !             else
         !               if(abs(ctmp).gt.1.0e-8) then
         !                 write(*,'(a49,i4)')
         !     1           '*** ERROR *** with iterative subspace at k-point ',
         !     2           nkp
         !                 write(*,*) 'vectors in clamp not orthonormal'
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
            wkomegai = real(num_wann,dp) * wbtot  
            do nnx = 1, nntot  
                  do j = 1, ndimwin(nnlist(nkp,nnx))  
                     do i = 1, ndimwin(nkp)  
                        cmk(i,j,nnx) = m_matrix_orig(i,j,nnx,nkp)  
                     enddo
                  enddo
            enddo
            do m = 1, num_wann  
               rlambda = dis_zeig(nkp,m,cmk) 
               wkomegai = wkomegai - rlambda  
            enddo

            womegai = womegai + wkomegai  

         enddo
         ! [Loop over k (nkp)]

         womegai = womegai / real(num_kpts,dp)  

         write(stdout,124) iter,womegai1*lenconfac**2,womegai*lenconfac**2,&
              100.0_dp*(womegai1-womegai)/womegai,io_time()


124      format(2x,i6,3x,f14.8,3x,f14.8,3x,es13.6,2x,f8.2,4x,'<-- DIS')

         ! Construct the updated Z matrix, CZMAT_OUT, at k points w/ non-frozen s
         do nkp = 1, num_kpts  
            if (num_wann.gt.ndimfroz(nkp) ) then  
               do nnx = 1, nntot  
                     k_pls_b = nnlist(nkp,nnx)  
                     do j = 1, ndimwin(k_pls_b)  
                        do i = 1, ndimwin(nkp)  
                           cmk(i,j,nnx) = m_matrix_orig(i,j,nnx,nkp)  
                        enddo
                     enddo
               enddo
               call dis_zmatrix(nkp,cmk,czmat_out(1,1,nkp))
            endif
         enddo

      enddo
      ! [BIG ITERATION LOOP (iter)]

      if (icompflag.eq.1) then
         if (iprint>2) then
            write(stdout,('(/4x,a)')) &
                 'WARNING: Complement subspace has zero dimensions at the following k-points:'
            i=0
            write(stdout,'(4x)',advance='no')
            do nkp=1,num_kpts
               if (ndimwin(nkp).eq.num_wann) then  
                  i=i+1
                  if (i.le.12) then
                     write(stdout,'(i6)',advance='no') nkp
                  else
                     i=1
                     write(stdout,'(/4x)',advance='no')
                     write(stdout,'(i6)',advance='no') nkp
                  endif
               endif
            enddo
         endif
      endif

      ! Write the final womegai. This should remain unchanged during the
      ! subsequent minimization of Omega_tilde in wannierise.f90
      ! We store it in the checkpoint file as a sanity check
      write(stdout,*) ' '
      write(stdout,'(/8x,a,f14.8,a)') 'Final Omega_I ',womegai*lenconfac**2,' ('//trim(length_unit)//'^2)'
      omega_invariant=womegai

      ! DIAGONALIZE THE HAMILTONIAN WITHIN THE OPTIMIZED SUBSPACES
      do nkp = 1, num_kpts  
         do j = 1, num_wann  
            do i = 1, num_wann  
               cham(i,j,nkp) = cmplx_0  
               do l = 1, ndimwin(nkp)  
                  cham(i,j,nkp) = cham(i,j,nkp) + conjg(clamp(l,i,nkp)) &
                       * clamp(l,j,nkp) * eigval_opt(l,nkp)
               enddo
            enddo
         enddo

         do j = 1, num_wann  
            do i = 1, j  
               cap(i + ( (j - 1) * j) / 2) = cham(i,j,nkp)  
            enddo
         enddo

         call ZHPEVX ('V', 'A', 'U', num_wann, cap, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
              m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)

         if (info.lt.0) then  
            write(stdout,*) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
            write(stdout,*) ' THE ',  -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
            call io_error(' dis_extract: error')   
         endif
         if (info.gt.0) then  
            write(stdout,*) ' *** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
            write(stdout,*) info, 'EIGENVECTORS FAILED TO CONVERGE'  
            call io_error(' dis_extract: error')   
         endif

         ! Store the energy eigenvalues of the optimal subspace (used in wann_ban
         do i = 1, num_wann  
            eigval_opt(i,nkp) = w(i)  
         enddo

         ! CALCULATE AMPLITUDES OF THE CORRESPONDING ENERGY EIGENVECTORS IN TERMS
         ! THE ORIGINAL ("WINDOW SPACE") ENERGY EIGENVECTORS
         do j = 1, num_wann  
            do i = 1, ndimwin(nkp)  
               ceamp(i,j,nkp) = cmplx_0  
               do l = 1, num_wann  
                  ceamp(i,j,nkp) = ceamp(i,j,nkp) + cz(l,j) * clamp(i,l,nkp)
               enddo
            enddo
         enddo
         ! NKP
      enddo

      ! DEBUG
      if (iprint>2) then
         write(stdout,'(/,a,/)') '  Eigenvalues inside optimal subspace:'  
         do nkp = 1, num_kpts  
            write(stdout,'(a,i3,2x,20(f9.5,1x))') '  K-point ', &
                 nkp, (eigval_opt(i,nkp), i = 1, num_wann)
         enddo
      endif
      ! ENDDEBUG

      ! Replace clamp by ceamp, which is given out to wannier.f (Both span the
      ! same space, but the latter is more convenient for the purpose of obtai
      ! an optimal Fourier-interpolated band structure: see Sec. III.E of SMV.
      do nkp = 1, num_kpts  
         do j = 1, num_wann  
            do i = 1, ndimwin(nkp)  
               clamp(i,j,nkp) = ceamp(i,j,nkp)  
            enddo
         enddo
      enddo

      if (icompflag.eq.1) then  
         if (iprint>2) then
            write(stdout,*) 'AT SOME K-POINT(S) COMPLEMENT SUBSPACE HAS ZERO DIMENSIONALITY'
            write(stdout,*) '=> DID NOT CREATE FILE COMPSPACE.DAT'  
         endif
      else  
         ! DIAGONALIZE THE HAMILTONIAN IN THE COMPLEMENT SUBSPACE, WRITE THE
         ! CORRESPONDING EIGENFUNCTIONS AND ENERGY EIGENVALUES
         do nkp = 1, num_kpts  
            do j = 1, ndimwin(nkp) - num_wann  
               do i = 1, ndimwin(nkp) - num_wann  
                  cham(i,j,nkp) = cmplx_0  
                  do l = 1, ndimwin(nkp)  
                     cham(i,j,nkp) = cham(i,j,nkp) + conjg(camp(l,i,nkp)) &
                          * camp(l,j,nkp) * eigval_opt(l,nkp)
                  enddo
               enddo
            enddo
            do j = 1, ndimwin(nkp) - num_wann  
               do i = 1, j  
                  cap(i + ( (j - 1) * j) / 2) = cham(i,j,nkp)  
               enddo
            enddo
            ndiff = ndimwin(nkp) - num_wann  
            call ZHPEVX ('V', 'A', 'U', ndiff, cap, 0.0_dp, 0.0_dp, 0, 0, &
                 -1.0_dp, m, w, cz, num_bands, cwork, rwork, iwork, ifail, info)
            if (info.lt.0) then  
               write(stdout,*) '*** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
               write(stdout,*) 'THE ',  -info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
               call io_error(' dis_extract: error')   
            endif
            if (info.gt.0) then  
               write(stdout,*) '*** ERROR *** ZHPEVX WHILE DIAGONALIZING HAMILTONIAN'
               write(stdout,*) info, 'EIGENVECTORS FAILED TO CONVERGE'  
               call io_error(' dis_extract: error')   
            endif

            ! CALCULATE AMPLITUDES OF THE ENERGY EIGENVECTORS IN THE COMPLEMENT SUBS
            ! TERMS OF THE ORIGINAL ENERGY EIGENVECTORS
            do j = 1, ndimwin(nkp) - num_wann  
               do i = 1, ndimwin(nkp)  
                  camp(i,j,nkp) = cmplx_0  
                  do l = 1, ndimwin(nkp) - num_wann  
                     camp(i,j,nkp) = camp(i,j,nkp) + cz(l,j) * clamp(i,l,nkp)
                  enddo
               enddo
            enddo
         enddo
         ! [loop over k points (nkp)]

      endif
      ! [if icompflag=1]


      deallocate(cham,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cham in dis_extract')
      deallocate(czmat_out,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating czmat_out in dis_extract')
      deallocate(czmat_in,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating czmat_in in dis_extract')
      deallocate(cmk,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cmk in dis_extract')
      deallocate(camp,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating camp in dis_extract')
      deallocate(ceamp,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating ceamp in dis_extract')
      deallocate(wkomegai1,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating wkomegai1 in dis_extract')

      deallocate(cz,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cz in dis_extract')
      deallocate(cwork,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cwork in dis_extract')
      deallocate(cap,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cap in dis_extract')
      deallocate(rwork,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating rwork in dis_extract')
      deallocate(w,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating w in dis_extract')
      deallocate(ifail,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating ifail in dis_extract')
      deallocate(iwork,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating iwork in dis_extract')


      return  


    contains


      !==================================================================!
      function dis_zeig(nkp,m,cmk)
      !==================================================================!
      !                                                                  !
      !                                                                  !
      !                                                                  !
      !                                                                  !
      !==================================================================!  
        
        ! Computes <lambda>_mk = sum_{n=1}^N sum_b w_b |<u_{mk}|u_{n,k+b}>|^2
        ! [See Eqs. (12) and (17) of SMV]

        implicit none

        integer, intent(in) :: nkp
        integer, intent(in) :: m
        complex(kind=dp), intent(in) :: cmk(num_bands,num_bands,nntot)

        ! Internal variables
        real(kind=dp) :: dis_zeig
        complex(kind=dp) :: cdot_bloch
        integer :: n,nnx,ndnnx,ndnn,nnsh,nkp2,l,j

        dis_zeig=0.0_dp

        do n = 1, num_wann  
           ! Loop over b-vectors
           do nnx = 1, nntot  
                 nkp2 = nnlist(nkp,nnx)  
                 ! Dotproduct
                 cdot_bloch = cmplx_0          
                 do l = 1, ndimwin(nkp)  
                    do j = 1, ndimwin(nkp2)  
                       cdot_bloch = cdot_bloch + &
                            conjg(clamp(l,m,nkp)) * clamp(j,n,nkp2) * cmk(l,j,nnx)
                    enddo
                 enddo
                 dis_zeig = dis_zeig + wb(nnx) * abs(cdot_bloch)**2  
              enddo
        enddo


        return  

      end function dis_zeig



      !==================================================================!
      subroutine dis_zmatrix(nkp,cmk,cmtrx)
      !==================================================================!
      !                                                                  !
      !                                                                  !
      !                                                                  !
      !                                                                  !
      !==================================================================!  

        implicit none

        integer, intent(in) :: nkp
        complex(kind=dp), intent(in) :: cmk(num_bands,num_bands,nntot)
        complex(kind=dp), intent(inout) :: cmtrx(num_bands,num_bands)

        ! OUTPUT:
        !
        ! CMTRX(M,N)        (M,N)-TH ENTRY IN THE
        !                   (NDIMWIN(NKP)-NDIMFROZ(NKP)) x (NDIMWIN(NKP)-NDIMFRO
        !                   HERMITIAN MATRIX AT THE NKP-TH K-POINT
        

        ! Internal variables
        integer :: l,m,n,j,nnx,ndnn,ndnnx,nnsh,k_pls_b
        complex(kind=dp) :: cdot_bloch1,cdot_bloch2

        ! LOOP OVER INDEPENDENT MATRIX ENTRIES
        do m = 1, ndimwin(nkp) - ndimfroz(nkp)  
           do n = 1, m  
              cmtrx(m,n) = cmplx_0  
              do l = 1, num_wann  
                 ! LOOP OVER B-VECTORS
                 do nnx = 1, nntot  
                       k_pls_b = nnlist(nkp,nnx)  
                       ! CALCULATE THE DOTPRODUCTS
                       cdot_bloch1 = cmplx_0  
                       cdot_bloch2 = cmplx_0  
                       do j = 1, ndimwin(k_pls_b)  
                          cdot_bloch1 = cdot_bloch1 + &
                               clamp(j,l,k_pls_b) * cmk(indxnfroz(m,nkp),j,nnx)
                          cdot_bloch2 = cdot_bloch2 + &
                               conjg(clamp(j,l,k_pls_b) * cmk(indxnfroz(n,nkp),j,nnx))
                       enddo
                       ! ADD CONTRIBUTION TO CMTRX(M,N)
                       cmtrx(m,n) = cmtrx(m,n) + wb(nnx) * cdot_bloch1 * cdot_bloch2
                       ! NNSH
                 enddo
                 ! L
              enddo
              ! hermiticity
              cmtrx(n,m) = conjg(cmtrx(m,n))  
              ! N
           enddo
           ! M
        enddo


        return  

      end subroutine dis_zmatrix


    end subroutine dis_extract



end module disentangle
