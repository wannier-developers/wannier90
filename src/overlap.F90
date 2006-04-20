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

module overlap
 
  use constants, only : dp,cmplx_0
  use parameters, only : disentanglement
  use io, only : stdout

  implicit none
 
  private

!!$  public :: overlap_dis_read
  public :: overlap_read
  public :: overlap_dealloc


contains

    

  !%%%%%%%%%%%%%%%%%%%%%
  subroutine overlap_read( )
  !%%%%%%%%%%%%%%%%%%%%%
    
    use parameters, only : num_bands, num_wann, num_kpts, nntot, nncell, nnlist,&
                           devel_flag, u_matrix, m_matrix, a_matrix, &
                           m_matrix_orig, u_matrix_opt, cp_pp
    use io,         only : io_file_unit, io_error, seedname

    implicit none

    integer :: nkp, nkp2, inn, nn, n, m, i, j
    integer :: mmn_in, amn_in, num_mmn, num_amn
    integer :: nnl, nnm, nnn, ncount
    integer :: nb_tmp, nkp_tmp, nntot_tmp, nw_tmp, ierr
    real(kind=dp) :: m_real, m_imag, a_real, a_imag
    complex(kind=dp), allocatable :: mmn_tmp(:,:)
    character(len=50) :: dummy
    logical :: nn_found


    allocate ( u_matrix( num_wann,num_wann,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating u_matrix in overlap_read')
    allocate ( m_matrix( num_wann,num_wann,nntot,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating m_matrix in overlap_read')

    if (disentanglement) then
       allocate(m_matrix_orig(num_bands,num_bands,nntot,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating m_matrix_orig in overlap_read')
       allocate(a_matrix(num_bands,num_wann,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating a_matrix in overlap_read')
       allocate(u_matrix_opt(num_bands,num_wann,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating u_matrix_opt in overlap_read')
    endif

    u_matrix = cmplx_0
    m_matrix = cmplx_0
    
    if (disentanglement) then
       m_matrix_orig = cmplx_0
       a_matrix      = cmplx_0
       u_matrix_opt  = cmplx_0
    endif


    if (index(devel_flag,'f77input')>0) then
       ! This block left for the short term as a mean
       ! to quickly benchmark against the old f77 code
       ! Read U_matrix and M_matrix from file 
       open(20,file='wannier0.dat',form='formatted',status='unknown')
       do i=1,num_wann
          do j=1,num_wann
             do nkp=1,num_kpts
                read(20,*) u_matrix(i,j,nkp)
                do nn=1,nntot
                   read(20,*) m_matrix(i,j,nn,nkp)
                end do
             end do
          end do
       end do
       close(20)

    else

       ! Read M_matrix_orig from file
       mmn_in=io_file_unit()
       open(unit=mmn_in,file=trim(seedname)//'.mmn',&
            form='formatted',status='old',action='read',err=101)
              
       write(stdout,'(/a)',advance='no') ' Reading overlaps from '//trim(seedname)//'.mmn    : '

       ! Read the comment line
       read(mmn_in,'(a)',err=103,end=103) dummy
       write(stdout,'(a)') trim(dummy)

       ! Read the number of bands, k-points and nearest neighbours
       read(mmn_in,*,err=103,end=103) nb_tmp,nkp_tmp,nntot_tmp

       ! Checks
       if (nb_tmp.ne.num_bands) &
            call io_error(trim(seedname)//'.mmn has not the right number of bands')
       if (nkp_tmp.ne.num_kpts) &
            call io_error(trim(seedname)//'.mmn has not the right number of k-points')
       if (nntot_tmp.ne.nntot) &
            call io_error(trim(seedname)//'.mmn has not the right number of nearest neighbours')

       ! Read the overlaps
       num_mmn=num_kpts*nntot
       allocate(mmn_tmp(num_bands,num_bands),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating mmn_tmp in overlap_read')
       do ncount = 1, num_mmn
          read(mmn_in,*,err=103,end=103) nkp,nkp2,nnl,nnm,nnn
          do n=1,num_bands
             do m=1,num_bands
                read(mmn_in,*,err=103,end=103) m_real, m_imag
                mmn_tmp(m,n) = cmplx(m_real,m_imag,kind=dp)
             enddo
          enddo
          nn=0
          nn_found=.false.
          do inn = 1, nntot
             if ((nkp2.eq.nnlist(nkp,inn)).and. &
                  (nnl.eq.nncell(1,nkp,inn)).and. &
                  (nnm.eq.nncell(2,nkp,inn)).and. &
                  (nnn.eq.nncell(3,nkp,inn)) ) then
                if (.not.nn_found) then
                   nn_found=.true.
                   nn=inn
                else
                   call io_error('Error reading '//trim(seedname)//'.mmn.&
                        & More than one matching nearest neighbour found')
                endif
             endif
          end do
          if (nn.eq.0) then
             write(stdout,'(/a,i8,2i5,i4,2x,3i3)') &
                  ' Error reading '//trim(seedname)//'.mmn:',ncount,nkp,nkp2,nn,nnl,nnm,nnn
             call io_error('Neighbour not found')
          end if
          if (disentanglement) then
             m_matrix_orig(:,:,nn,nkp) = mmn_tmp(:,:)
          else
             m_matrix(:,:,nn,nkp) = mmn_tmp(:,:)
          end if
       end do
       deallocate(mmn_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating mmn_tmp in overlap_read')
 
       close(mmn_in)

       ! Read A_matrix from file wannier.amn
       amn_in=io_file_unit()
       open(unit=amn_in,file=trim(seedname)//'.amn',form='formatted',status='old',err=102)
       
       write(stdout,'(/a)',advance='no') ' Reading projections from '//trim(seedname)//'.amn : '

       ! Read the comment line
       read(amn_in,'(a)',err=104,end=104) dummy
       write(stdout,'(a)') trim(dummy)

       ! Read the number of bands, k-points and wannier functions
       read(amn_in,*,err=104,end=104) nb_tmp, nkp_tmp, nw_tmp
       
       ! Checks
       if (nb_tmp.ne.num_bands) &
            call io_error(trim(seedname)//'.amn has not the right number of bands')
       if (nkp_tmp.ne.num_kpts) &
            call io_error(trim(seedname)//'.amn has not the right number of k-points')
       if (nw_tmp.ne.num_wann) &
            call io_error(trim(seedname)//'.amn has not the right number of Wannier functions')
       
       ! Read the projections
       num_amn = num_bands*num_wann*num_kpts
       if (disentanglement) then
          do ncount = 1, num_amn
             read(amn_in,*,err=104,end=104) m,n,nkp,a_real,a_imag
             a_matrix(m,n,nkp) = cmplx(a_real,a_imag,kind=dp)
          end do
       else
          do ncount = 1, num_amn
             read(amn_in,*,err=104,end=104) m,n,nkp,a_real,a_imag
             u_matrix(m,n,nkp) = cmplx(a_real,a_imag,kind=dp)
          end do
       end if

       close(amn_in)
       
       ! If post-processing a Car-Parinello calculation (gamma only)
       ! then rotate M and A to the basis of Kohn-Sham eigenstates
       if (cp_pp) call overlap_rotate()

       
       ! If we don't need to disentangle we can now convert from A to U
       ! And rotate M accordingly
       if(.not.disentanglement .and. (.not.cp_pp) ) call overlap_project

    endif


return 
101    call io_error('Error: Problem opening input file '//trim(seedname)//'.mmn')
102    call io_error('Error: Problem opening input file '//trim(seedname)//'.amn')
103    call io_error('Error: Problem reading input file '//trim(seedname)//'.mmn')
104    call io_error('Error: Problem reading input file '//trim(seedname)//'.amn')


  end subroutine overlap_read

  
  !%%%%%%%%%%%%%%%%%%%%%
  subroutine overlap_rotate
  !%%%%%%%%%%%%%%%%%%%%%

    use parameters, only : num_bands,a_matrix,m_matrix_orig,nntot
    use io,         only : io_file_unit,io_error

    implicit none

    integer       :: lam_unit,info,inn,i,j
    real(kind=DP) :: lambda(num_bands,num_bands)
    real(kind=DP) :: AP(num_bands*(num_bands+1)/2)
    real(kind=DP) :: eig(num_bands),work(3*num_bands)

    lam_unit=io_file_unit()
    open(unit=lam_unit,file='lambda.dat',&
         form='unformatted',status='old',action='read')
!!$    write(stdout,*) ' Reading lambda.dat...' 
    read(lam_unit) lambda
!!$    write(stdout,*) ' done'
    close(lam_unit)

    do j=1,num_bands
       do i=1,j
          AP(i+(j-1)*j/2)=0.5_dp*(lambda(i,j)+lambda(j,i))
       end do
    end do

    CALL DSPEV('V','U', num_bands, AP, eig, lambda, num_bands, work, info)
    if(info.ne.0) &
         call io_error('Diagonalization of lambda in overlap_rotate failed')

    ! For debugging
!!$    write(stdout,*) 'EIGENVALUES - CHECK WITH CP OUTPUT'
!!$    do i=1,num_bands
!!$       write(stdout,*) 13.6058*eig(i)
!!$    end do

    ! Rotate M_mn
    do inn=1,nntot
       m_matrix_orig(:,:,inn,1) = &
            matmul(transpose(lambda),matmul(m_matrix_orig(:,:,inn,1),lambda))
    end do

    ! Rotate A_mn
    a_matrix(:,:,1) = matmul(transpose(lambda),a_matrix(:,:,1))

    ! For debugging
!!$    ! Write rotated A and M
!!$    do i=1,num_bands
!!$       do j=1,num_wann
!!$          write(12,'(2i5,a,2f18.12)') i,j,'   1',a_matrix(i,j,1)
!!$       enddo
!!$    enddo
!!$    do inn=1,nntot
!!$       do i=1,num_bands
!!$          do j=1,num_wann
!!$             write(11,'(2i5,2a)') i,j,'    1','    1'
!!$             write(11,'(2f18.12)') m_matrix_orig(i,j,inn,1)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    stop

    return

  end subroutine overlap_rotate



  !%%%%%%%%%%%%%%%%%%%%%
  subroutine overlap_dealloc( )
  !%%%%%%%%%%%%%%%%%%%%%

    use parameters, only : u_matrix,m_matrix,m_matrix_orig,&
                       a_matrix,u_matrix_opt
    use io,     only : io_error

    implicit none

    integer :: ierr

    if (allocated(u_matrix_opt)) then
       deallocate( u_matrix_opt, stat=ierr )
       if (ierr/=0) call io_error('Error deallocating u_matrix_opt in overlap_dealloc')
    end if
    if (allocated( a_matrix) ) then
       deallocate( a_matrix, stat=ierr )
       if (ierr/=0) call io_error('Error deallocating a_matrix in overlap_dealloc')
    end if
    if (allocated( m_matrix_orig)) then
       deallocate( m_matrix_orig, stat=ierr )
       if (ierr/=0) call io_error('Error deallocating m_matrix_orig in overlap_dealloc')
    endif

    deallocate ( m_matrix, stat=ierr )
    if (ierr/=0) call io_error('Error deallocating m_matrix in overlap_dealloc')
    deallocate ( u_matrix, stat=ierr )
    if (ierr/=0) call io_error('Error deallocating u_matrix in overlap_dealloc')

    return

  end subroutine overlap_dealloc



  !==================================================================!
  subroutine overlap_project()
  !==================================================================!
  !                                                                  !
  !  Note that in this subroutine num_wann = num_bands               ! 
  !  since, if we are here, then disentanglement = FALSE             !
  !                                                                  !
  !                                                                  !
  !==================================================================!  
    use constants
    use io,         only : io_error
    use parameters, only : num_bands,num_wann,num_kpts,&
                           u_matrix,m_matrix,nntot,nnlist
    use utility,    only : utility_zgemm

    implicit none


    ! internal variables
    integer :: i,j,m,nkp,info,ierr,nn,nkp2
    real(kind=dp),    allocatable :: svals(:)
    real(kind=dp)                 :: rwork(5*num_bands) 
    complex(kind=dp)              :: ctmp2
    complex(kind=dp), allocatable :: cwork(:)
    complex(kind=dp), allocatable :: cz(:,:)
    complex(kind=dp), allocatable :: cvdag(:,:)


    allocate(svals(num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating svals in overlap_project')
    allocate(cz(num_bands,num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cz in overlap_project')
    allocate(cvdag(num_bands,num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cvdag in overlap_project')
    allocate(cwork(4*num_bands),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cwork in overlap_project')


    ! Calculate the transformation matrix CU = CS^(-1/2).CA,
    ! where CS = CA.CA^\dagger.

    do nkp = 1, num_kpts  
       !
       ! SINGULAR VALUE DECOMPOSITION
       !
       call ZGESVD('A', 'A', num_bands, num_bands, u_matrix(1,1,nkp), &
            num_bands, svals, cz, num_bands, cvdag, num_bands, cwork, &
            4*num_bands, rwork, info)
       if (info.ne.0) then  
          write(stdout,*) ' ERROR: IN ZGESVD IN overlap_project'  
          write(stdout,*) ' K-POINT NKP=', nkp, ' INFO=', info  
          if (info.lt.0) then  
             write(stdout,*) ' THE ',  -info, '-TH ARGUMENT HAD ILLEGAL VALUE'  
          endif
          call io_error('DISENTANGLE: Error in ZGESVD in overlap_project')
       endif

!       u_matrix(:,:,nkp)=matmul(cz,cvdag)
       call utility_zgemm(u_matrix(:,:,nkp),cz,'N',cvdag,'N',num_wann)

       !
       ! CHECK UNITARITY
       !
       do i = 1, num_bands
          do j = 1, num_bands 
             ctmp2 = cmplx_0  
             do m = 1, num_bands  
                ctmp2 = ctmp2 + u_matrix(m,j,nkp) * conjg(u_matrix(m,i,nkp))  
             enddo
             if ( (i.eq.j).and.(abs(ctmp2-cmplx_1).gt.0.00001_dp) ) then
                write(stdout,*) ' ERROR: unitarity of initial U'  
                write(stdout,'(1x,a,i2)') 'nkp= ', nkp  
                write(stdout,'(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j  
                write(stdout,'(1x,a,f12.6,1x,f12.6)') &
                     '[u_matrix.transpose(u_matrix)]_ij= ',&
                     real(ctmp2,dp),aimag(ctmp2)
                call io_error('Error in unitarity of initial U in overlap_project')
             endif
             if ( (i.ne.j) .and. (abs(ctmp2).gt.0.00001_dp) ) then  
                write(stdout,*) ' ERROR: unitarity of initial U'  
                write(stdout,'(1x,a,i2)') 'nkp= ', nkp  
                write(stdout,'(1x,a,i2,2x,a,i2)') 'i= ', i, 'j= ', j  
                write(stdout,'(1x,a,f12.6,1x,f12.6)') &
                     '[u_matrix.transpose(u_matrix)]_ij= ', &
                     real(ctmp2,dp),aimag(ctmp2)
                call io_error('Error in unitarity of initial U in overlap_project')
             endif
          enddo
       enddo
    enddo
    ! NKP


    ! so now we have the U's that rotate the wavefunctions at each k-point.
    ! the matrix elements M_ij have also to be updated 
    do nkp=1, num_kpts 
       do nn=1,nntot
          nkp2=nnlist(nkp,nn)
          ! cvdag = U^{dagger} . M   (use as workspace)
          call utility_zgemm(cvdag,u_matrix(:,:,nkp),'C',m_matrix(:,:,nn,nkp),'N',num_wann)
          ! cz = cvdag . U
          call utility_zgemm(cz,cvdag,'N',u_matrix(:,:,nkp2),'N',num_wann)
          m_matrix(:,:,nn,nkp) = cz(:,:)
       end do
    end do
    
    deallocate(cwork,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cwork in overlap_project')
    deallocate(cvdag,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cvdag in overlap_project')
    deallocate(cz,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cz in overlap_project')
    deallocate(svals,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating svals in overlap_project')


    return  

  end subroutine overlap_project



end module overlap
