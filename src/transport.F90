!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2004-2007 Jonathan Yates, Arash Mostofi,     !
! Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt  !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!
!
!-----------------------------------------------------------------------!
!  Based on                                                             !
!  < dosqc_1.0 >                                                        !
!  Density Of States and Quantum Conductance - Version 1.0              !
!  Marco Buongiorno Nardelli, January 2000.                             !
!                                                                       !
!  Reference:                                                           !
!  - M. Buongiorno Nardelli, "Electronic transport in extended systems: !
!  application to carbon nanotubes", Phys. Rev. B, vol. 60(11), 7828    !
!  (1999)                                                               !
!-----------------------------------------------------------------------!

!=======================================================================!
! Definition of parameters used in w90_transport                        !
!=======================================================================!
!                                                                       !
! transport_mode       = 'bulk' or 'lcr'                                !
! tran_win_min         = minimum E                                      !
! tran_win_max         = maximum E                                      !
! tran_energy_step     = delta E                                        !
! tran_num_bb          = # of WFs in a principal layer of a perfectly   !
!                        periodic bulk system                           !
! tran_num_ll          = # of WFs in a principal layer of a left lead   !            
! tran_num_rr          = # of WFs in a principal layer of a right lead  !             
! tran_num_cc          = # of WFs in a disordered conducter cell        !       
! tran_num_lc          = # of WFs in a disordered conducter cell that   !
!                        are used to calculate interaction with a       !
!                        left lead                                      !
! tran_num_cr          = # of WFs in a disordered conducter cell that   !
!                        are used to calculate interaction with a       !
!                        right lead                                     !
! tran_num_bandc       = width of band-diagonal hC matrix               !
! tran_read_ht         = .true. => read H matrix from h*.dat files      !
! tran_write_ht        = .true. => write H matrix from h*.dat files     !   
! tran_use_same_lead   = .true. => in L-C-R construction, left and      !
!                                  right lead are the same kind         !
! tran_num_cell_ll     = # of unit cells in a PL of left lead           !
! tran_num_cell_rr     = # of unit cells in a PL of right lead          !
!                        (equal to tran_num_cell_ll for now)            !
! tran_group_threshold = distance defining the grouping of WFs          !
! tran_wf_threshold    = threshold of integrals in parity subroutine    !
!                        to be considered zero                          !
!=======================================================================!

module w90_transport

  use w90_constants, only : dp

  implicit none

  private

! small complex number 
  complex(kind=dp), parameter :: eta=(0.d0,0.0005d0)

! nterx  = # of maximum iteration to calculate transfer matrix
  integer, parameter :: nterx=50
  ! cartesian axis to which real_lattice(:,one_dim_vec) is parallel
  integer :: one_dim_vec
  integer :: nrpts_one_dim
  ! num_pl : number of unit cell in a principal layer
  integer :: num_pl

  real(kind=dp), allocatable :: hr_one_dim(:,:,:)
  real(kind=dp), allocatable :: hB0(:,:)
  real(kind=dp), allocatable :: hB1(:,:)
  real(kind=dp), allocatable :: hL0(:,:)
  real(kind=dp), allocatable :: hL1(:,:)
  real(kind=dp), allocatable :: hR0(:,:)
  real(kind=dp), allocatable :: hR1(:,:)
  real(kind=dp), allocatable :: hC(:,:)
  real(kind=dp), allocatable :: hLC(:,:)
  real(kind=dp), allocatable :: hCR(:,:)


  public :: tran_main

contains
  !==================================================================!
  subroutine tran_main()
    !==================================================================!

    use w90_io,         only : stdout,io_stopwatch
    use w90_parameters, only : transport_mode,tran_read_ht,timing_level,hr_plot,num_wann
    use w90_hamiltonian,only : hamiltonian_get_hr,hamiltonian_write_hr,hamiltonian_setup

    implicit none

    integer,dimension(num_wann) :: tran_sorted_idx
 
    if (timing_level>0) call io_stopwatch('tran: main',1) 

    write(stdout,'(/1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,'(1x,a)') '|                              TRANSPORT                                    |'
    write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,*)

    if (index(transport_mode,'bulk')>0 ) then
      write(stdout,'(/1x,a/)') 'Calculation of Quantum Conductance and DoS: bulk mode'
       if (.not.tran_read_ht) then
          call hamiltonian_setup()
          call hamiltonian_get_hr()
          if (hr_plot) call hamiltonian_write_hr()
          call tran_reduce_hr()
          call tran_cut_hr_one_dim()
          call tran_get_ht()
       end if
       call tran_bulk()
    end if

    if (index(transport_mode,'lcr')>0 ) then
       write(stdout,'(/1x,a/)') 'Calculation of Quantum Conductance and DoS: lead-conductor-lead mode'
       if (.not.tran_read_ht) then
          call hamiltonian_setup()
          call hamiltonian_get_hr()
          call tran_lcr_sort(tran_sorted_idx)
          call tran_hr_parity_shift(tran_sorted_idx)
          call tran_lcr_build_ham(tran_sorted_idx) 
       endif 
       call tran_lcr()
    end if
  
    if (timing_level>0) call io_stopwatch('tran: main',2)

  end subroutine tran_main

  !==================================================================!
  subroutine tran_reduce_hr()
    !==================================================================!
    !
    ! reduce ham_r from 3-d to 1-d
    !
    use w90_constants,   only : dp, eps8
    use w90_io,          only : io_error, io_stopwatch, stdout
    use w90_parameters,  only : one_dim_dir,real_lattice,num_wann, &
                                mp_grid,timing_level
    use w90_hamiltonian, only : irvec,nrpts,ham_r

    integer :: ierr
    integer :: irvec_max, irvec_tmp(3), two_dim_vec(2)
    integer :: i, j
    integer :: i1, i2, i3, n1, nrpts_tmp, loop_rpt

    if (timing_level>1) call io_stopwatch('tran: reduce_hr',1)

    ! Find one_dim_vec which is parallel to one_dim_dir
    ! two_dim_vec - the other two lattice vectors
    j = 0
    do i=1,3
       if ( abs(abs(real_lattice(one_dim_dir,i)) &
            - sqrt(dot_product(real_lattice(:,i),real_lattice(:,i)))) .lt. eps8 ) then
          one_dim_vec = i
          j = j +1 
       end if
    end do
    if ( j .ne. 1 ) then
       write(stdout,'(i3,a)') j,' : 1-D LATTICE VECTOR NOT DEFINED'
       call io_error('Error: 1-d lattice vector not defined in tran_reduce_hr')
    end if
       
    j=0
    do i=1,3
       if ( i .ne. one_dim_vec ) then
          j = j +1
          two_dim_vec(j)=i
       end if
    end do

    ! starting H matrix should include all W-S supercell where 
    ! the center of the cell spans the full space of the home cell
    ! adding one more buffer layer when mp_grid(one_dim_vec) is an odd number

    !irvec_max = (mp_grid(one_dim_vec)+1)/2
    irvec_tmp = maxval(irvec,DIM=2)+1
    irvec_max = irvec_tmp(one_dim_vec)
    nrpts_one_dim = 2*irvec_max+1 
    allocate(hr_one_dim(num_wann,num_wann,-irvec_max:irvec_max),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hr_one_dim in tran_reduce_hr') 
    hr_one_dim = 0.0_dp
       
    ! check imaginary part
    write(stdout,'(1x,a,F12.6)') 'Maximum imaginary part of the real-space Hamiltonian: ',maxval(abs(aimag(ham_r)))

    ! select a subset of ham_r, where irvec is 0 along the two other lattice vectors

    nrpts_tmp = 0
loop_n1: do n1 = -irvec_max, irvec_max
       do loop_rpt = 1,nrpts
          i1 = mod(n1 - irvec(one_dim_vec,loop_rpt),mp_grid(one_dim_vec))
          i2 = irvec(two_dim_vec(1),loop_rpt)
          i3 = irvec(two_dim_vec(2),loop_rpt)
          if (i1.eq.0 .and. i2.eq.0 .and. i3.eq.0 ) then
             nrpts_tmp = nrpts_tmp+1
             hr_one_dim(:,:,n1) = real(ham_r(:,:,loop_rpt),dp)
             cycle loop_n1
          end if
       end do
    end do loop_n1   

    if (nrpts_tmp .ne. nrpts_one_dim ) then   
       write(stdout,'(a)') 'FAILED TO EXTRACT 1-D HAMILTONIAN'
       call io_error('Error: cannot extract 1d hamiltonian in tran_reduce_hr')
    end if

    if (timing_level>1) call io_stopwatch('tran: reduce_hr',2)

    return

  end subroutine tran_reduce_hr

  !==================================================================!
  subroutine tran_cut_hr_one_dim()
    !==================================================================!
    !
    use w90_constants,   only : dp
    use w90_io,          only : io_stopwatch,stdout
    use w90_parameters,  only : num_wann,mp_grid,timing_level,real_lattice,&
                                hr_cutoff,dist_cutoff,dist_cutoff_mode, &
                                one_dim_dir,length_unit
    use w90_hamiltonian, only : wannier_centres_translated

    implicit none
    !
    integer :: irvec_max
    integer :: i, j, n1
    real(kind=dp) :: hr_max
    real(kind=dp) :: dist
    real(kind=dp) :: dist_vec(3)
    real(kind=dp) :: dist_ij_vec(3)
    real(kind=dp) :: shift_vec(3,-nrpts_one_dim/2:nrpts_one_dim/2)
    real(kind=dp) :: hr_tmp(num_wann,num_wann)
    !
    if (timing_level>1) call io_stopwatch('tran: cut_hr_one_dim',1)
    !
    irvec_max = nrpts_one_dim/2
    ! maximum possible dist_cutoff
    dist = real(mp_grid(one_dim_vec),dp)*abs(real_lattice(one_dim_dir,one_dim_vec))/2.0_dp

    if ( dist_cutoff .gt. dist ) then
       write(stdout,'(1x,a,1x,F10.5,1x,a)') 'dist_cutoff',dist_cutoff,trim(length_unit),'is too large'
       dist_cutoff = dist
       write(stdout,'(4x,a,1x,F10.5,1x,a)') 'reset to',dist_cutoff,trim(length_unit)
    end if

    do n1 = -irvec_max, irvec_max
       shift_vec(:,n1) = real(n1,dp)*(real_lattice(:,one_dim_vec))
    !       write(stdout,'(a,3f10.6)') 'shift_vec', shift_vec(:,n1)
    end do

    ! apply dist_cutoff first
    if ( index(dist_cutoff_mode,'one_dim')>0 ) then
       do i=1,num_wann
          do j=1,num_wann
             dist_ij_vec(one_dim_dir)=wannier_centres_translated(one_dim_dir,i)-wannier_centres_translated(one_dim_dir,j)
             do n1 = -irvec_max, irvec_max
                dist_vec(one_dim_dir) = dist_ij_vec(one_dim_dir)+ shift_vec(one_dim_dir,n1)
                dist = abs(dist_vec(one_dim_dir))
                if ( dist .gt. dist_cutoff ) hr_one_dim(j,i,n1)=0.0_dp
             end do
          end do
        end do
    else
        do i=1,num_wann
           do j=1,num_wann
              dist_ij_vec(:)=wannier_centres_translated(:,i) - wannier_centres_translated(:,j)
              do n1 = -irvec_max,irvec_max
                 dist_vec(:) =  dist_ij_vec(:)+ shift_vec(:,n1) 
                 dist = sqrt(dot_product(dist_vec,dist_vec))
                 if ( dist .gt. dist_cutoff ) hr_one_dim(j,i,n1)=0.0_dp
              end do
           end do
        end do
    end if

    ! output maximum to check a decay of H as a function of lattice vector R
    write(stdout,'(/1x,a78)') repeat('-',78)
    write(stdout,'(1x,4x,a)') &
                'Maximum real part of the real-space Hamiltonian at each lattice point'
    write(stdout,'(1x,8x,a62)') repeat('-',62)
    write(stdout,'(1x,11x,a,11x,a)') 'Lattice point R', 'Max |H_ij(R)|'
    ! calculate number of units inside a principal layer
    num_pl = 0
    do n1=-irvec_max,irvec_max
       hr_tmp(:,:)=abs(hr_one_dim(:,:,n1))       
       hr_max = maxval(hr_tmp) 
       if ( hr_max .gt. hr_cutoff ) then
          if (abs(n1) .gt. num_pl) num_pl = abs(n1)
       else
          hr_one_dim(:,:,n1)=0.0_dp
       end if
       write(stdout,'(1x,9x,5x,I5,5x,12x,F12.6)') n1, hr_max 
    end do
    write(stdout,'(1x,8x,a62)') repeat('-',62)

    write(stdout,'(/1x,a,I6)') 'Number of unit cells inside the principal layer:',num_pl 
    write(stdout,'(1x,a,I6)')  'Number of Wannier Functions inside the principal layer:',num_pl*num_wann 
    ! apply hr_cutoff to each element inside the principal layer
    do n1 = -num_pl , num_pl
       do i=1,num_wann
          do j=1,num_wann
             if ( abs(hr_one_dim(j,i,n1)) .lt. hr_cutoff ) hr_one_dim(j,i,n1)=0.0_dp
          end do
       end do
    end do

    if (timing_level>1) call io_stopwatch('tran: cut_hr_one_dim',2)

    return
  
  end subroutine tran_cut_hr_one_dim

  !==================================================================!
    subroutine tran_get_ht()
    !==================================================================!
    !  construct h00 and h01
    !==================================================================!
    !
    use w90_constants,  only : dp
    use w90_io,         only : io_error, io_stopwatch, seedname, io_date, &
                               io_file_unit
    use w90_parameters, only : num_wann, tran_num_bb, tran_write_ht, &
                               fermi_energy, timing_level
    !
    implicit none
    !
    integer :: ierr, file_unit
    integer :: i, j, n1, im, jm
    character(len=9)   :: cdate, ctime
    !
    if (timing_level>1) call io_stopwatch('tran: get_ht',1)
    !
    tran_num_bb=num_pl*num_wann
    !
    allocate(hB0(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hB0 in tran_get_ht')
    allocate(hB1(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hB1 in tran_get_ht')
    !
    hB0 = 0.0_dp 
    hB1 = 0.0_dp 
    !
    ! h00
    do j=0, num_pl-1
       do i=0, num_pl-1
          n1 = i-j
          im=i*num_wann
          jm=j*num_wann
          hB0(jm+1:jm+num_wann,im+1:im+num_wann)=hr_one_dim(:,:,n1)
       end do
    end do
  
    ! h01
    do j=1, num_pl
       do i=0, j-1
          n1 = i-(j-1)+num_pl
          im=i*num_wann
          jm=(j-1)*num_wann
          hB1(jm+1:jm+num_wann,im+1:im+num_wann)=hr_one_dim(:,:,n1)
       end do
    end do

    ! shift by fermi_energy
    do i=1,tran_num_bb
       hB0(i,i)=hB0(i,i)-fermi_energy
    end do

    if ( tran_write_ht ) then
  
       file_unit = io_file_unit()
       open(file_unit,file=trim(seedname)//'_htB.dat',status='unknown',form='formatted',action='write')

       call io_date(cdate,ctime)
       write(file_unit,*) 'written on '//cdate//' at '//ctime ! Date and time
       write(file_unit,'(I6)') tran_num_bb
       write(file_unit,'(6F12.6)') ((hB0(j,i),j=1,tran_num_bb),i=1,tran_num_bb)
       write(file_unit,'(I6)') tran_num_bb
       write(file_unit,'(6F12.6)') ((hB1(j,i),j=1,tran_num_bb),i=1,tran_num_bb)
     
       close(file_unit)
 
    end if
 
    if (timing_level>1) call io_stopwatch('tran: get_ht',2)

    return
  
  end subroutine tran_get_ht

  !==================================================================!
  subroutine tran_bulk()
    !==================================================================!

    use w90_constants,  only : dp, cmplx_0, cmplx_1, cmplx_i, pi
    use w90_io,         only : io_error, io_stopwatch, seedname, io_date, &
                               io_file_unit, stdout
    use w90_parameters, only : tran_num_bb, tran_read_ht,  &
                               tran_win_min, tran_win_max, tran_energy_step, &
                               timing_level

    implicit none

    integer :: qc_unit, dos_unit
    integer :: ierr
    integer :: n_e, n, i
    real(kind=dp) ::  qc, dos 
    real(kind=dp) ::  e_scan
    complex(kind=dp) :: e_scan_cmp
    complex(kind=dp), allocatable, dimension(:,:) :: tot, tott
    complex(kind=dp), allocatable, dimension(:,:) :: g_B, gR, gL
    complex(kind=dp), allocatable, dimension(:,:) :: sLr, sRr
    complex(kind=dp), allocatable, dimension(:,:) :: s1, s2, c1
    character(len=50) :: filename
    character(len=9)  :: cdate, ctime
 
    if (timing_level>1) call io_stopwatch('tran: bulk',1)

    allocate (tot(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tot in tran_bulk')
    allocate (tott(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tott in tran_bulk')
    allocate (g_B(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating g_B in tran_bulk')
    allocate (gL(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating gL in tran_bulk')
    allocate (gR(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating gR in tran_bulk')
    allocate (sLr(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating sLr in tran_bulk')
    allocate (sRr(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating sRr in tran_bulk')
    allocate (s1(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s1 in tran_bulk')
    allocate (s2(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s2 in tran_bulk')
    allocate (c1(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating c1 in tran_bulk')

    call io_date(cdate,ctime)

    qc_unit = io_file_unit()
    open(qc_unit, file=trim(seedname)//'_qc.dat',status='unknown', &
         form='formatted',action='write')
    write(qc_unit,*) '## written on '//cdate//' at '//ctime ! Date and time

    dos_unit = io_file_unit()
    open(dos_unit, file=trim(seedname)//'_dos.dat',status='unknown', &
         form='formatted',action='write') 
    write(dos_unit,*) '## written on '//cdate//' at '//ctime ! Date and time

    !   set up the layer hamiltonians

    if (tran_read_ht)  then
       allocate(hB0(tran_num_bb,tran_num_bb),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating hB0 in tran_bulk')
       allocate(hB1(tran_num_bb,tran_num_bb),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating hB1 in tran_bulk')
       filename = trim(seedname)//'_htB.dat'
       call tran_read_htX(tran_num_bb,hB0,hB1,filename)                     
    end if

    !   loop over the energies

    n_e = floor((tran_win_max-tran_win_min)/tran_energy_step)+1

    write(stdout,'(/1x,a)',advance='no') 'Calculating quantum conductance and density of states...'

    do n=1,n_e
       e_scan = tran_win_min + real(n-1,dp)*tran_energy_step
 
       ! compute conductance according to Fisher and Lee
       ! retarded Green

       e_scan_cmp = e_scan+eta 
       call tran_transfer(tot,tott,hB0,hB1,e_scan_cmp,tran_num_bb)
       call tran_green(tot,tott,hB0,hB1,e_scan,g_B,0,1,tran_num_bb) 
 
       ! compute S_Lr and S_Rr

       c1(:,:) = cmplx(hB1(:,:),kind=dp)

       ! Self-energy (Sigma_L^r) : sLr = (hB1)^+ * tott
       ! Self-energy (Sigma_R^r) : sRr = (hB1)   * tot
       sLr = cmplx_0
       sRr = cmplx_0
       call ZGEMM('C','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,c1,tran_num_bb,tott,tran_num_bb,cmplx_0,sLr,tran_num_bb)
       call ZGEMM('N','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,c1,tran_num_bb,tot, tran_num_bb,cmplx_0,sRr,tran_num_bb)

       ! Gamma_L = i(Sigma_L^r-Sigma_L^a)
       gL = cmplx_i*(sLr - conjg(transpose(sLr)))
       ! Gamma_R = i(Sigma_R^r-Sigma_R^a)
       gR = cmplx_i*(sRr - conjg(transpose(sRr)))

       s1 = cmplx_0
       s2 = cmplx_0
       c1 = cmplx_0
       ! s1 = Gamma_L * g_B^r
       call ZGEMM('N','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,gL,tran_num_bb,g_B,tran_num_bb,cmplx_0,s1,tran_num_bb)
       ! s2 = Gamma_L * g_B^r * Gamma_R 
       call ZGEMM('N','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,s1,tran_num_bb,gR, tran_num_bb,cmplx_0,s2,tran_num_bb)
       ! c1 = Gamma_L * g_B^r * Gamma_R * g_B^a
       call ZGEMM('N','C',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,s2,tran_num_bb,g_B,tran_num_bb,cmplx_0,c1,tran_num_bb)
           
       qc = 0.0_dp
       do i=1,tran_num_bb
          qc = qc + real(c1(i,i),dp)
       end do
       write(qc_unit,'(f12.6,f15.6)') e_scan, qc

       dos = 0.0_dp
       do i=1,tran_num_bb
          dos = dos - aimag(g_B(i,i))
       end do
       dos = dos / pi
       write(dos_unit,'(f12.6,f15.6)') e_scan, dos

    end do

    write(stdout,'(a)') ' done'

    close(qc_unit)
    close(dos_unit)

    deallocate (c1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating c1 in tran_bulk')
    deallocate (s2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s2 in tran_bulk')
    deallocate (s1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s1 in tran_bulk')
    deallocate (sRr,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sRr in tran_bulk')
    deallocate (sLr,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sLr in tran_bulk')
    deallocate (gR,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating gR in tran_bulk')
    deallocate (gL,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating gL in tran_bulk')
    deallocate (g_B,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating g_B in tran_bulk')
    deallocate (tott,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tott in tran_bulk')
    deallocate (tot,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tot in tran_bulk')


    if (timing_level>1) call io_stopwatch('tran: bulk',2)

    return

  end subroutine tran_bulk 

  !==================================================================!
  subroutine tran_lcr()
    !==================================================================!

    use w90_constants,  only : dp, cmplx_0, cmplx_1, cmplx_i, pi
    use w90_io,         only : io_error, io_stopwatch, seedname, io_date, &
                               stdout, io_file_unit
    use w90_parameters, only : tran_num_ll, tran_num_rr, tran_num_cc, tran_num_lc, &
                               tran_num_cr, tran_num_bandc, &
                               tran_win_min, tran_win_max, tran_energy_step,      &
                               tran_use_same_lead, timing_level, tran_read_ht

    implicit none

    integer :: qc_unit, dos_unit
    integer :: ierr                
    integer :: KL, KU, KC
    integer :: n_e, n, i, j, k, info
    integer, allocatable :: ipiv(:)
    real(kind=dp) ::  qc, dos
    real(kind=dp) ::  e_scan
    real(kind=dp), allocatable, dimension(:,:) :: hCband
    complex(kind=dp) :: e_scan_cmp
    complex(kind=dp), allocatable, dimension(:,:) :: hLC_cmp, hCR_cmp,     &
                              totL, tottL, totR, tottR,                    &
                              g_surf_L, g_surf_R, g_C, g_C_inv,            &
                              gR, gL, sLr, sRr, s1, s2, c1, c2
    character(len=50) :: filename 
    character(len=9)  :: cdate, ctime

    if (timing_level>1) call io_stopwatch('tran: lcr',1)

    call io_date(cdate,ctime)

    qc_unit = io_file_unit()
    open(qc_unit, file=trim(seedname)//'_qc.dat',status='unknown', &
         form='formatted',action='write')
    write(qc_unit,*) '## written on '//cdate//' at '//ctime ! Date and time

    dos_unit = io_file_unit()
    open(dos_unit, file=trim(seedname)//'_dos.dat',status='unknown', &
         form='formatted',action='write') 
    write(dos_unit,*) '## written on '//cdate//' at '//ctime ! Date and time

    KL = max(tran_num_lc, tran_num_cr, tran_num_bandc) - 1
    KU = KL
    KC = max(tran_num_lc, tran_num_cr)

    allocate (hCband(2*KL+KU+1,tran_num_cc),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hCband in tran_lcr')
    allocate (hLC_cmp(tran_num_ll,tran_num_lc),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hLC_cmp in tran_lcr')
    allocate (hCR_cmp(tran_num_cr,tran_num_rr),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hCR_cmp in tran_lcr')

    !If construct used only when reading matrices from file
    if (tran_read_ht) then
        allocate (hL0(tran_num_ll,tran_num_ll),stat=ierr)
        if (ierr/=0) call io_error('Error in allocating hL0 in tran_lcr')
        allocate (hL1(tran_num_ll,tran_num_ll),stat=ierr)
        if (ierr/=0) call io_error('Error in allocating hL1 in tran_lcr')
        allocate (hC(tran_num_cc,tran_num_cc),stat=ierr)
        if (ierr/=0) call io_error('Error in allocating hC in tran_lcr')
        allocate (hLC(tran_num_ll,tran_num_lc),stat=ierr)
        if (ierr/=0) call io_error('Error in allocating hLC in tran_lcr')
        allocate (hCR(tran_num_cr,tran_num_rr),stat=ierr)
        if (ierr/=0) call io_error('Error in allocating hCR in tran_lcr')

        filename = trim(seedname)//'_htL.dat'
        call tran_read_htX(tran_num_ll,hL0,hL1,filename)
    
        if (.not. tran_use_same_lead ) then
           allocate (hR0(tran_num_rr,tran_num_rr),stat=ierr)
           if (ierr/=0) call io_error('Error in allocating hR0 in tran_lcr')
           allocate (hR1(tran_num_rr,tran_num_rr),stat=ierr)
           if (ierr/=0) call io_error('Error in allocating hR1 in tran_lcr')
           filename = trim(seedname)//'_htR.dat'
           call tran_read_htX(tran_num_rr,hR0,hR1,filename)
        end if

        filename = trim(seedname)//'_htC.dat'
        call tran_read_htC(tran_num_cc,hC,filename)
        filename = trim(seedname)//'_htLC.dat'
        call tran_read_htXY(tran_num_ll,tran_num_lc,hLC,filename)
        filename = trim(seedname)//'_htCR.dat'
        call tran_read_htXY(tran_num_cr,tran_num_rr,hCR,filename)
    endif

    !  Banded matrix H_C  :  save memory !
    do j=1, tran_num_cc
       do i=max(1,j-KU), min(tran_num_cc, j+KL)
          hCband(KL+KU+1+i-j,j)=hC(i,j)
       end do
    end do
    deallocate(hC,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating hC in tran_lcr')

    !  H_LC : to a complex matrix
    hLC_cmp(:,:) = cmplx(hLC(:,:),kind=dp)
    deallocate(hLC,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating hLC in tran_lcr')

    !  H_CR : to a complex matrix
    hCR_cmp(:,:) = cmplx(hCR(:,:),kind=dp)
    deallocate(hCR,stat=ierr)      
    if (ierr/=0) call io_error('Error in deallocating hCR in tran_lcr')

    allocate (totL(tran_num_ll,tran_num_ll),stat=ierr) 
    if (ierr/=0) call io_error('Error in allocating totL in tran_lcr')
    allocate (tottL(tran_num_ll,tran_num_ll),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tottL in tran_lcr')
    if (.not. tran_use_same_lead) then
       allocate (totR(tran_num_rr,tran_num_rr),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating totR in tran_lcr')
       allocate (tottR(tran_num_rr,tran_num_rr),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating tottR in tran_lcr')
    end if
    allocate(g_surf_L(tran_num_ll,tran_num_ll),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating g_surf_L in tran_lcr')
    allocate(g_surf_R(tran_num_rr,tran_num_rr),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating g_surf_R in tran_lcr')
    allocate(g_C_inv(2*KL+KU+1,tran_num_cc),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating g_C_inv in tran_lcr')
    allocate(g_C(tran_num_cc,tran_num_cc),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating g_C in tran_lcr')
    allocate(sLr(tran_num_lc,tran_num_lc),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating sLr in tran_lcr')
    allocate(sRr(tran_num_cr,tran_num_cr),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating sRr in tran_lcr')
    allocate(gL(tran_num_lc,tran_num_lc),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating gL in tran_lcr')
    allocate(gR(tran_num_cr,tran_num_cr),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating gR in tran_lcr')
    allocate(c1(tran_num_lc,tran_num_ll),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating c1 in tran_lcr')
    allocate(c2(tran_num_cr,tran_num_rr),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating c2 in tran_lcr')
    allocate(s1(KC,KC),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s1 in tran_lcr')
    allocate(s2(KC,KC),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s2 in tran_lcr')
    allocate(ipiv(tran_num_cc),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ipiv in tran_lcr')

    !  Loop over the energies
    n_e = floor((tran_win_max-tran_win_min)/tran_energy_step)+1

    write(stdout,'(/1x,a)',advance='no') 'Calculating quantum conductance and density of states...'

    do n=1,n_e

       e_scan = tran_win_min + real(n-1,dp)*tran_energy_step
 
       !    compute conductance according to Fisher and Lee
       !    compute self-energies following Datta

       e_scan_cmp = e_scan+eta

       ! Surface green function for the left lead : g_surf_L
       call tran_transfer(totL,tottL,hL0,hL1,e_scan_cmp,tran_num_ll)
       call tran_green(totL,tottL,hL0,hL1,e_scan,g_surf_L,-1,1,tran_num_ll)

       ! Self-energy (Sigma_L) : sLr = (hLC_cmp)^+ * g_surf_L * hLC_cmp
       c1 = cmplx_0
       sLr = cmplx_0
       call ZGEMM('C','N',tran_num_lc,tran_num_ll,tran_num_ll,cmplx_1,&
          hLC_cmp,tran_num_ll,g_surf_L,tran_num_ll,cmplx_0,c1,tran_num_lc)
       call ZGEMM('N','N',tran_num_lc,tran_num_lc,tran_num_ll,cmplx_1,&
          c1,tran_num_lc,hLC_cmp,tran_num_ll,cmplx_0,sLr,tran_num_lc)

       ! Surface green function for the right lead : g_surf_R
       if (tran_use_same_lead) then
          call tran_green(totL,tottL,hL0,hL1,e_scan,g_surf_R,1,1,tran_num_rr)
       else
          call tran_transfer(totR,tottR,hR0,hR1,e_scan_cmp,tran_num_rr)
          call tran_green(totR,tottR,hR0,hR1,e_scan,g_surf_R,1,1,tran_num_rr)
       end if 

       ! Self-energy (Sigma_R) : sRr = hCR_cmp * g_surf_R * (hCR_cmp)^+
       c2 = cmplx_0
       sRr = cmplx_0
       call ZGEMM('N','N',tran_num_cr,tran_num_rr,tran_num_rr,cmplx_1,&
         hCR_cmp,tran_num_cr,g_surf_R,tran_num_rr,cmplx_0,c2,tran_num_cr)
       call ZGEMM('N','C',tran_num_cr,tran_num_cr,tran_num_rr,cmplx_1,&
         c2,tran_num_cr,hCR_cmp,tran_num_cr,cmplx_0,sRr,tran_num_cr)

       ! g_C^-1 = -H
       g_C_inv(:,:) = cmplx(-hCband(:,:),kind=dp)

       ! g_C^-1 = -H - Sigma_L^r        
       do j=1,tran_num_lc
          do i=max(1,j-KU),min(tran_num_lc,j+KL)
             g_C_inv(KL+KU+1+i-j,j)=g_C_inv(KL+KU+1+i-j,j)-sLr(i,j)
          end do
       end do

       ! g_C^-1 = -H - Sigma_L^r - Sigma_R^r
       do j=(tran_num_cc-tran_num_cr)+1,tran_num_cc
          do i=max((tran_num_cc-tran_num_cr)+1,j-(tran_num_cr-1)),min(tran_num_cc,j+(tran_num_cr-1))
             g_C_inv(KL+KU+1+i-j,j)=g_C_inv(KL+KU+1+i-j,j)-sRr(i-(tran_num_cc-tran_num_cr),j-(tran_num_cc-tran_num_cr))
          end do
       end do

       ! g_C^-1 = eI - H - Sigma_L^r - Sigma_R^r
       do i=1,tran_num_cc
          g_C_inv(KL+KU+1,i) = e_scan + g_C_inv(KL+KU+1,i)
       end do

       ! invert g_C^-1 => g_C
       g_C = cmplx_0
       do i=1,tran_num_cc
          g_C(i,i)= cmplx_1
       end do

       call ZGBSV(tran_num_cc,KL,KU,tran_num_cc,g_C_inv,2*KL+KU+1,ipiv,g_C,tran_num_cc,info)
       if (info.ne.0) then
          write(stdout,*) 'ERROR: IN ZGBSV IN tran_lcr, INFO=', info 
          call io_error('tran_lcr: problem in ZGBSV')
       end if

       ! Gamma_L = i(Sigma_L^r-Sigma_L^a)
       gL = cmplx_i*(sLr - conjg(transpose(sLr)))

       ! s1 = Gamma_L * g_C^r
       s1 = cmplx_0
       do j=1,KC
          do i=1,tran_num_lc
             do k=1,tran_num_lc
                s1(i,j)=s1(i,j)+gL(i,k)*g_C(k,j+(tran_num_cc-KC))
             end do
          end do
       end do

       ! Gamma_R = i(Sigma_R^r-Sigma_R^a)
       gR = cmplx_i*(sRr - conjg(transpose(sRr)))

       ! s2 = Gamma_R * g_C^a
       s2=cmplx_0
       do j=1,KC
          do i=1,tran_num_cr
             do k=1,tran_num_cr
                s2(i+(KC-tran_num_cr),j)=s2(i+(KC-tran_num_cr),j)               &
                   +gR(i,k)*conjg(g_C(j,k+(tran_num_cc-tran_num_cr)))
             end do
          end do
       end do

       qc = 0.0_dp
       do i=1,KC
          do j=1,KC
             qc = qc + real(s1(i,j)*s2(j,i),dp)
          end do
       end do
       write(qc_unit,'(f12.6,f15.6)') e_scan, qc

       ! compute density of states for the conductor layer

       dos = 0.0_dp
       do i=1,tran_num_cc
          dos = dos - aimag(g_C(i,i))
       end do
       dos = dos / pi
       write(dos_unit,'(f12.6,f15.6)') e_scan, dos

    end do

    write(stdout,'(a)') ' done'

    close(qc_unit)
    close(dos_unit)

    deallocate(ipiv,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating ipiv in tran_lcr')
    deallocate(s2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s2 in tran_lcr')
    deallocate(s1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s1 in tran_lcr')
    deallocate(c2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating c2 in tran_lcr')
    deallocate(c1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating c1 in tran_lcr')
    deallocate(gR,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating gR in tran_lcr')
    deallocate(gL,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating gL in tran_lcr')
    deallocate(sRr,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sRr in tran_lcr')
    deallocate(sLr,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sLr in tran_lcr')
    deallocate(g_C,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating g_C in tran_lcr')
    deallocate(g_C_inv,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating g_C_inv in tran_lcr')
    deallocate(g_surf_R,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating g_surf_R in tran_lcr')
    deallocate(g_surf_L,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating g_surf_L in tran_lcr')
    if(allocated(tottR)) deallocate (tottR,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tottR in tran_lcr')
    if(allocated(totR)) deallocate (totR,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating totR in tran_lcr')
    deallocate (tottL,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tottL in tran_lcr')
    deallocate (totL,stat=ierr) 
    if (ierr/=0) call io_error('Error in deallocating totL in tran_lcr')
    deallocate (hCR_cmp,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating hCR_cmp in tran_lcr')
    deallocate (hLC_cmp,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating hLC_cmp in tran_lcr')
    deallocate (hCband,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating hCband in tran_lcr')

    if (timing_level>1) call io_stopwatch('tran: lcr',2)

    return

  end subroutine tran_lcr 

  !==================================================================!
  subroutine tran_transfer(tot,tott,h_00,h_01,e_scan_cmp,nxx) 
    !==================================================================!
    !                                                                  !
    ! iterative construction of the transfer matrix                    !
    ! as Lopez-Sancho^2&Rubio, J.Phys.F:Met.Phys., v.14, 1205 (1984)   !
    ! and ibid. v.15, 851 (1985)                                       !
    !                                                                  !
    !===================================================================

    use w90_constants, only : dp, cmplx_0, cmplx_1, eps7
    use w90_io, only : stdout, io_error

    implicit none

    integer, intent(in) :: nxx
    complex(kind=dp), intent(in) ::  e_scan_cmp
    complex(kind=dp), intent(out) ::  tot(nxx,nxx)
    complex(kind=dp), intent(out) ::  tott(nxx,nxx)
    real(kind=dp), intent(in) :: h_00(nxx,nxx)
    real(kind=dp), intent(in) :: h_01(nxx,nxx)
    !
    integer  :: ierr, info
    integer  :: i, j, n, nxx2
    integer, allocatable :: ipiv(:)
    real(kind=dp) :: conver,conver2
    complex(kind=dp), allocatable, dimension(:,:) :: tsum, tsumt
    complex(kind=dp), allocatable, dimension(:,:) :: t11, t12
    complex(kind=dp), allocatable, dimension(:,:) :: s1, s2
    complex(kind=dp), allocatable, dimension(:,:,:) :: tau, taut

    allocate(ipiv(nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ipiv in tran_transfer')
    allocate(tsum(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tsum in tran_transfer')
    allocate(tsumt(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tsumt in tran_transfer')
    allocate(t11(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating t11 in tran_transfer')
    allocate(t12(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating t12 in tran_transfer')
    allocate(s1(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s1 in tran_transfer')
    allocate(s2(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s2 in tran_transfer')
    allocate(tau(nxx,nxx,2),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tau in tran_transfer')
    allocate(taut(nxx,nxx,2),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating taut in tran_transfer')

    nxx2 = nxx*nxx
    
    tot = cmplx_0
    tott = cmplx_0

    ! construction of the transfer matrix
    ! t12 = e - h_00
    t12(:,:) = cmplx(-h_00(:,:),kind=dp)
    do i=1,nxx
       t12(i,i) = e_scan_cmp + t12(i,i)
    end do

    ! compute (e - h_00)^-1 and store it in t11 
    t11 = cmplx_0
    do i=1,nxx
       t11(i,i) = cmplx_1
    end do

    ! inverse of t12 -> t11
    call ZGESV(nxx,nxx,t12,nxx,ipiv,t11,nxx,info)
    if (info.ne.0) then
       write(stdout,*) 'ERROR:  IN ZGESV IN tran_transfer, INFO=', info 
       call io_error('tran_transfer: problem in ZGESV 1') 
    end if

    ! compute intermediate t-matrices (defined as tau(nxx,nxx,niter)
    ! and taut(...)):
    tau = cmplx_0
    taut = cmplx_0

    ! t_0:
    t12(:,:) = cmplx(h_01(:,:),kind=dp)

    !  tau  = ( e - H_00 )^-1 * H_01^+
    call ZGEMM('N','C',nxx,nxx,nxx,cmplx_1,t11,nxx,t12,nxx,cmplx_0,tau(1,1,1),nxx)
    !  taut = ( e - H_00 )^-1 * H_01
    call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,t11,nxx,t12,nxx,cmplx_0,taut(1,1,1),nxx)

    !   initialize T:
    tot(:,:) =tau(:,:,1)
    tsum(:,:)=taut(:,:,1)

    !   initialize T^bar:
    tott(:,:) =taut(:,:,1)
    tsumt(:,:)=tau(:,:,1)

    !   main loop:
    do n=1,nterx
                
       t11 = cmplx_0
       t12 = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tau(1,1,1),nxx,taut(1,1,1),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,taut(1,1,1),nxx,tau(1,1,1),nxx,cmplx_0,t12,nxx)

       s1(:,:) = -t11(:,:)-t12(:,:)
       do i=1,nxx
          s1(i,i) = cmplx_1+s1(i,i)
       end do

       s2 = cmplx_0
       do i=1,nxx
          s2(i,i)=cmplx_1 
       end do

       call ZGESV(nxx,nxx,s1,nxx,ipiv,s2,nxx,info)
       if (info.ne.0) then
          write(stdout,*) 'ERROR:  IN ZGESV IN tran_transfer, INFO=', info 
          call io_error('tran_transfer: problem in ZGESV 2') 
       end if

       t11 = cmplx_0
       t12 = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tau(1,1,1),nxx,tau(1,1,1),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,taut(1,1,1),nxx,taut(1,1,1),nxx,cmplx_0,t12,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,s2,nxx,t11,nxx,cmplx_0,tau(1,1,2),nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,s2,nxx,t12,nxx,cmplx_0,taut(1,1,2),nxx)

       !   put the transfer matrices together

       t11 = cmplx_0
       s1  = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsum,nxx,tau(1,1,2),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsum,nxx,taut(1,1,2),nxx,cmplx_0,s1,nxx)
       call ZCOPY(nxx2,t11,1,s2,1)
       call ZAXPY(nxx2,cmplx_1,tot,1,s2,1)

       tot(:,:) = s2(:,:)
       tsum(:,:)= s1(:,:)

       t11 = cmplx_0
       s1  = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsumt,nxx,taut(1,1,2),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsumt,nxx,tau(1,1,2),nxx,cmplx_0,s1,nxx) 
       call ZCOPY(nxx2,t11,1,s2,1)
       call ZAXPY(nxx2,cmplx_1,tott,1,s2,1)

       tott(:,:) = s2(:,:)
       tsumt(:,:)= s1(:,:)

       tau(:,:,1) = tau(:,:,2)
       taut(:,:,1)= taut(:,:,2)

       ! convergence check on the t-matrices

       conver = 0.0_dp
       conver2 = 0.0_dp

       do j=1,nxx
          do i=1,nxx
              conver=conver+sqrt(real(tau(i,j,2),dp)**2+aimag(tau(i,j,2))**2)
              conver2=conver2+sqrt(real(taut(i,j,2),dp)**2+aimag(taut(i,j,2))**2)
          end do
       end do

       if (conver.lt.eps7 .and. conver2.lt.eps7) return
    end do 

    if (conver.gt.eps7 .or. conver2.gt.eps7) &
       call io_error('Error in converging transfer matrix in tran_transfer') 

    deallocate(ipiv,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating ipiv in tran_transfer')
    deallocate(tsum,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tsum in tran_transfer')
    deallocate(tsumt,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tsumt in tran_transfer')
    deallocate(t11,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating t11 in tran_transfer')
    deallocate(t12,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating t12 in tran_transfer')
    deallocate(s1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s1 in tran_transfer')
    deallocate(s2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s2 in tran_transfer')
    deallocate(tau,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tau in tran_transfer')
    deallocate(taut,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating taut in tran_transfer')

    return 

  end subroutine tran_transfer

  !==================================================================!
  subroutine tran_green(tot,tott,h_00,h_01,e_scan, g,igreen,invert,nxx)
    !==================================================================!
    !   construct green's functions
    !   
    !   igreen = -1  left surface
    !   igreen =  1  right surface
    !   igreen =  0  bulk
   
    !   invert = 0 computes g^-1
    !   invert = 1 computes g^-1 and g
    !==================================================================!

    use w90_constants, only : dp, cmplx_0, cmplx_1
    use w90_io, only : stdout, io_error

    implicit none

    integer, intent(in) :: nxx
    integer, intent(in) :: igreen
    integer, intent(in) :: invert
    real(kind=dp),     intent(in) :: e_scan 
    real(kind=dp),     intent(in) :: h_00(nxx,nxx), h_01(nxx,nxx)
    complex(kind=dp),  intent(in) :: tot(nxx,nxx),tott(nxx,nxx)
    complex(kind=dp), intent(out) :: g(nxx,nxx)

    integer :: ierr, info
    integer :: i
    integer, allocatable :: ipiv(:)
    complex(kind=dp), allocatable, dimension(:,:) :: g_inv, eh_00
    complex(kind=dp), allocatable, dimension(:,:) :: s1, s2, c1           

    allocate(ipiv(nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ipiv in tran_green')
    allocate(g_inv(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating g_inv in tran_green')
    allocate(eh_00(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating eh_00 in tran_green')
    allocate(c1(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating c1 in tran_green')
    allocate(s1(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating s1 in tran_green')
    allocate(s2(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating s2 in tran_green')

    c1(:,:)=cmplx(h_01(:,:),kind=dp)

    select case(igreen)

      case(1) 

       ! construct the surface green's function g00 

       s1 = cmplx_0
       ! s1 = H_01 * T
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tot,nxx,cmplx_0,s1,nxx)

       ! eh_00 =  -H_00 - H_01*T
       eh_00(:,:) = cmplx(-h_00(:,:),kind=dp)-s1(:,:)
       ! eh_00 = e_scan -H_00 - H_01*T
       do i=1,nxx
          eh_00(i,i) = cmplx(e_scan,kind=dp) + eh_00(i,i)
       end do

       g_inv(:,:) = eh_00(:,:)

       ! identity
       g = cmplx_0
       do i=1,nxx
          g(i,i)=cmplx_1 
       end do
    
       if (invert.eq.1) then
          call ZGESV(nxx,nxx,eh_00,nxx,ipiv,g,nxx,info)
          if (info.ne.0) then
             write(stdout,*) 'ERROR:  IN ZGESV IN tran_green, INFO=', info 
             call io_error('tran_green: problem in ZGESV 1') 
          end if
       end if

      case(-1)

       !  construct the dual surface green's function gbar00 

       s1 = cmplx_0
       ! s1 = H_01^+ * T^bar
       call ZGEMM('C','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tott,nxx,cmplx_0,s1,nxx)

       ! s1 = -H_00 - H_01^+ * T^bar
       eh_00(:,:) = cmplx(-h_00(:,:),kind=dp)-s1(:,:)
       ! s1 = e_scan - H_00 - H_01^+ * T^bar
       do i=1,nxx
          eh_00(i,i) = cmplx(e_scan,kind=dp) + eh_00(i,i)
       end do
     
       g_inv(:,:) = eh_00(:,:)

       ! identity
       g = cmplx_0
       do i=1,nxx
          g(i,i)=cmplx_1
       end do

       if (invert.eq.1) then
          call ZGESV(nxx,nxx,eh_00,nxx,ipiv,g,nxx,info)
          if (info.ne.0) then
             write(stdout,*) 'ERROR:  IN ZGESV IN tran_green, INFO=', info 
             call io_error('tran_green: problem in ZGESV 2') 
          end if
       end if
    
      case(0)

      !  construct the bulk green's function gnn or (if surface=.true.) the
      !  sub-surface green's function

       s1 = cmplx_0
       s2 = cmplx_0
       ! s1 = H_01 * T
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tot,nxx,cmplx_0,s1,nxx)
       ! s2 = H_01^+ * T^bar
       call ZGEMM('C','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tott,nxx,cmplx_0,s2,nxx)

       eh_00(:,:) = cmplx(-h_00(:,:),kind=dp)-s1(:,:)-s2(:,:)
       do i=1,nxx
          eh_00(i,i) = cmplx(e_scan,kind=dp) + eh_00(i,i)
       end do

       g_inv(:,:) = eh_00(:,:)

       ! identity
       g = cmplx_0
       do i=1,nxx
          g(i,i)=cmplx_1 
       end do
     
       if (invert.eq.1) then
          call ZGESV(nxx,nxx,eh_00,nxx,ipiv,g,nxx,info)
          if (info.ne.0) then
             write(stdout,*) 'ERROR:  IN ZGESV IN tran_green, INFO=', info 
             call io_error('tran_green: problem in ZGESV 3') 
          end if
       end if

    end select

    deallocate(s2)
    if (ierr/=0) call io_error('Error in deallocating s2 in tran_green')
    deallocate(s1)
    if (ierr/=0) call io_error('Error in deallocating s1 in tran_green')
    deallocate(c1)
    if (ierr/=0) call io_error('Error in deallocating c1 in tran_green')
    deallocate(eh_00)
    if (ierr/=0) call io_error('Error in deallocating eh_00 in tran_green')
    deallocate(g_inv)
    if (ierr/=0) call io_error('Error in deallocating g_inv in tran_green')
    deallocate(ipiv)
    if (ierr/=0) call io_error('Error in deallocating ipiv in tran_green')

    return

  end subroutine tran_green 

  !============================================!
   subroutine tran_read_htX(nxx,h_00,h_01,h_file)
    !============================================!

    use w90_constants, only : dp
    use w90_io, only : stdout, io_file_unit, io_error, maxlen
 
    implicit none

    integer, intent(in) ::  nxx
    real(kind=dp), intent(out) :: h_00(nxx,nxx), h_01(nxx,nxx)
    character(len=50), intent(in) :: h_file
    !
    integer :: i, j, nw, file_unit
    character(len=maxlen) :: dummy

    file_unit = io_file_unit()

    open(unit=file_unit, file=h_file, form='formatted', &
         status='old', action='read',err=101)

    write(stdout,'(/a)',advance='no') ' Reading H matrix from '//h_file//'  : '

    read(file_unit,'(a)',err=102,end=102) dummy
    write(stdout, '(a)') trim(dummy)

    read(file_unit,*,err=102,end=102) nw
    if (nw.ne.nxx) call io_error('wrong matrix size in transport: read_htX')            
    read(file_unit,*) ((h_00(i,j),i=1,nxx),j=1,nxx)
    read(file_unit,*,err=102,end=102) nw
    if (nw.ne.nxx) call io_error('wrong matrix size in transport: read_htX')
    read(file_unit,*,err=102,end=102) ((h_01(i,j),i=1,nxx),j=1,nxx)

    close(unit=file_unit)

    return

101    call io_error('Error: Problem opening input file '//h_file)
102    call io_error('Error: Problem reading input file '//h_file)

  end subroutine tran_read_htX

  !============================================!
  subroutine tran_read_htC(nxx,h_00,h_file)
    !============================================!

    use w90_constants, only : dp
    use w90_io, only : stdout, io_file_unit, io_error, maxlen
 
    implicit none

    integer, intent(in) ::  nxx
    real(kind=dp), intent(out) :: h_00(nxx,nxx)
    character(len=50), intent(in) :: h_file
    ! 
    integer :: i, j, nw, file_unit
    character(len=maxlen) :: dummy

    file_unit = io_file_unit()

    open(unit=file_unit, file=h_file, form='formatted', &
         status='old', action='read',err=101)

    write(stdout,'(/a)',advance='no') ' Reading H matrix from '//h_file//'  : '

    read(file_unit,'(a)',err=102,end=102) dummy
    write(stdout, '(a)') trim(dummy)

    read(file_unit,*,err=102,end=102) nw
    if (nw.ne.nxx) call io_error('wrong matrix size in transport: read_htC')
    read(file_unit,*,err=102,end=102) ((h_00(i,j),i=1,nxx),j=1,nxx)

    close(unit=file_unit)

    return

101    call io_error('Error: Problem opening input file '//h_file)
102    call io_error('Error: Problem reading input file '//h_file)

  end subroutine tran_read_htC

  !============================================!
  subroutine tran_read_htXY(nxx1,nxx2,h_01,h_file)
    !============================================!

    use w90_constants, only : dp
    use w90_io, only : stdout, io_file_unit, io_error, maxlen
 
    implicit none

    integer, intent(in) ::  nxx1,nxx2
    real(kind=dp), intent(out) :: h_01(nxx1,nxx2)
    character(len=50), intent(in) :: h_file
    ! 
    integer :: i, j, nw1, nw2, file_unit
    character(len=maxlen) :: dummy

    file_unit = io_file_unit()

    open(unit=file_unit, file=h_file, form='formatted', &
         status='old', action='read',err=101)

    write(stdout,'(/a)',advance='no') ' Reading H matrix from '//h_file//'  : '

    read(file_unit,'(a)',err=102,end=102) dummy
    write(stdout, '(a)') trim(dummy)

    read(file_unit,*,err=102,end=102) nw1, nw2

    if (nw1.ne.nxx1 .or. nw2.ne.nxx2) call io_error('wrong matrix size in transport: read_htXY')
    
    read(file_unit,*,err=102,end=102) ((h_01(i,j),i=1,nxx1),j=1,nxx2)

    close(unit=file_unit)  

    return

101    call io_error('Error: Problem opening input file '//h_file)
102    call io_error('Error: Problem reading input file '//h_file)
 
  end subroutine tran_read_htXY

  !========================================!
  subroutine tran_lcr_sort(tran_sorted_idx)
    !========================================!

    use w90_constants,          only : dp
    use w90_io,                 only : io_error,stdout,io_stopwatch
    use w90_parameters,         only : one_dim_dir,tran_num_ll,tran_num_rr,num_wann,&
                                       real_lattice,tran_group_threshold,iprint,timing_level
    use w90_hamiltonian,        only : wannier_centres_translated

    implicit none

    integer,intent(out),dimension(num_wann)           :: tran_sorted_idx
 
    real(dp),dimension(2,num_wann)                    :: centres_non_sorted,centres_initial_sorted
    real(dp),dimension(2,tran_num_ll)                 :: PL1,PL2,PL3,PL4,PL
    real(dp),dimension(2,num_wann-(4*tran_num_ll))    :: central_region
    real(dp)                                          :: reference_position,cell_length,distance

    integer                                           :: i,j,k,PL_selector,max_i,iterator
    integer,allocatable,dimension(:)                  :: PL_groups,PL1_groups,PL2_groups,PL3_groups,PL4_groups,central_region_groups
    integer,dimension(3)                              :: coord

    character(30)                                     :: fmt_1

    if (timing_level>1) call io_stopwatch('tran: lcr_sort',1)

    !
    !Check Plotting regime has been entered and translated centres have been found
    !
    if (size(wannier_centres_translated) .eq. 0) then
        call io_error('Plotting regime must find translated centres in order to perform lcr transport')
    endif

    !read one_dim_dir and creates an array (coord) that correspond to the 
    !conduction direction (coord(1)) and the two perpendicular directions 
    !(coord(2),coord(3)), such that a right-handed set is formed
    !    
    if (one_dim_dir .eq. 1) then
        coord(1)=1
        coord(2)=2
        coord(3)=3
    elseif (one_dim_dir .eq. 2) then
        coord(1)=2
        coord(2)=3
        coord(3)=1
    elseif (one_dim_dir .eq. 3) then
        coord(1)=3
        coord(2)=1
        coord(3)=2
    endif
    !
    !Check
    !
    if (((real_lattice(coord(1),coord(2)) .ne. 0) .or. (real_lattice(coord(1),coord(3)) .ne. 0)) .or. &
        ((real_lattice(coord(2),coord(1)) .ne. 0) .or. (real_lattice(coord(3),coord(1)) .ne. 0))) then
        call io_error(&
        'Lattice vector in conduction direction must point along x,y or z direction and be orthogonal to the remaining lattice vectors.')
    endif
    !
    !Check
    !
    if (num_wann < 4*tran_num_ll) then 
        call io_error('Principle layers are too big.')
    endif

    write(stdout,*)'------------------------- 2c2 Calculation Type: ------------------------------'
    write(stdout,*)' '

100 continue
    !
    !Extract a 2d array of the wannier_indices and their coord(1) from wannier_centers_translated
    !
    do i=1,num_wann
        centres_non_sorted(1,i)=i
        centres_non_sorted(2,i)=wannier_centres_translated(coord(1),i)
    enddo
    write(stdout,*)' Sorting WFs into principal layers'
    !
    !Initial sorting according to coord(1).
    !
    call sort(centres_non_sorted,centres_initial_sorted)
    !
    !Extract principal layers. WARNING: This extraction implies the structure of the supercell is
    !2 principal layers of lead on the left and on the right of a central conductor.
    !
    PL1=centres_initial_sorted(:,1:tran_num_ll)
    PL2=centres_initial_sorted(:,tran_num_ll+1:2*tran_num_ll)
    PL3=centres_initial_sorted(:,num_wann-(2*tran_num_ll-1):num_wann-(tran_num_ll))
    PL4=centres_initial_sorted(:,num_wann-(tran_num_ll-1):)

    if (iprint .ge. 4) then
        write(stdout,*)' Group Breakdown of each principal layer'
    endif

    do i=1,4 
        !
        !Creating a variable PL_selector which choose the appropriate PL
        !
        PL_selector=i                                 
        select case(PL_selector)
        case(1)
            PL=PL1
        case(2)
            PL=PL2
        case(3)
            PL=PL3
        case(4)
            PL=PL4
        endselect
        !
        !Grouping wannier functions with similar coord(1)
        !
        call group(PL,PL_groups)

        if (iprint .ge. 4) then
            !
            !Print group breakdown
            !
            write(fmt_1,'(i5)')size(PL_groups)
            fmt_1=adjustl(fmt_1)
            fmt_1='(a3,i1,a1,i5,a2,'//trim(fmt_1)//'i4,a1)'
            write(stdout,fmt_1)' PL',i,' ',size(PL_groups),' (',(PL_groups(j),j=1,size(PL_groups)),')'
        endif
        !
        !Returns the sorted PL
        !
        call master_sort_and_group(PL,PL_groups,coord,tran_num_ll) 
        
        select case(PL_selector)
        case(1)
            PL1=PL
            allocate(PL1_groups(size(PL_groups)))
            PL1_groups=PL_groups
        case(2)
            PL2=PL
            allocate(PL2_groups(size(PL_groups)))
            PL2_groups=PL_groups
        case(3)
            PL3=PL
            allocate(PL3_groups(size(PL_groups)))
            PL3_groups=PL_groups
        case(4)
            PL4=PL
            allocate(PL4_groups(size(PL_groups)))
            PL4_groups=PL_groups
        endselect  

        deallocate(PL_groups)        

    enddo
    !
    !Checks
    !
    if ((size(PL1_groups) .ne. size(PL2_groups)) .or.&
        (size(PL2_groups) .ne. size(PL3_groups)) .or.&
        (size(PL3_groups) .ne. size(PL4_groups))) then 
        call io_error('Inconsistent number of groups among principal layers')
    endif

    do i=1,size(PL1_groups)
        if ((PL1_groups(i) .ne. PL2_groups(i)) .or.&
            (PL2_groups(i) .ne. PL3_groups(i)) .or.&
            (PL3_groups(i) .ne. PL4_groups(i))) then 
            call io_error('Inconsitent number of wannier function among similar groups within principal layers')
        endif
    enddo
    !
    !Grouping and sorting of central conductor region
    !
    !Define central region
    !
    central_region=centres_initial_sorted(:,2*tran_num_ll+1:num_wann-(2*tran_num_ll))
    !
    !Group central region
    !
    call group(central_region,central_region_groups)
    !
    !Print central region group breakdown
    !
    if (iprint .ge. 4) then
        write(stdout,*)' Group Breakdown of central region'
        write(fmt_1,'(i5)')size(central_region_groups)
        fmt_1=adjustl(fmt_1)
        fmt_1='(a5,i5,a2,'//trim(fmt_1)//'i4,a1)'
        write(stdout,fmt_1)'     ',size(central_region_groups),' (',(central_region_groups(j),j=1,size(central_region_groups)),')'
    endif
    !
    !Returns sorted central group region
    !
    call master_sort_and_group(central_region,central_region_groups,coord,num_wann-(4*tran_num_ll))

    !
    !Build the sorted index array
    !
    tran_sorted_idx                                                     =centres_initial_sorted(1,:)
    tran_sorted_idx(1:tran_num_ll)                                      =PL1(1,:)
    tran_sorted_idx(tran_num_ll+1:2*tran_num_ll)                        =PL2(1,:)
    tran_sorted_idx(2*tran_num_ll+1:num_wann-(2*tran_num_ll))           =central_region(1,:)
    tran_sorted_idx(num_wann-(2*tran_num_ll-1):num_wann-(tran_num_ll))  =PL3(1,:)
    tran_sorted_idx(num_wann-(tran_num_ll-1):)                          =PL4(1,:)
    !
    !Now we check that the leftmost group and the rightmost group aren't
    !supposed to be the same group
    !
    reference_position=wannier_centres_translated(coord(1),tran_sorted_idx(1))
    cell_length=real_lattice(coord(1),coord(1))
    iterator=1
    do i=1,tran_num_ll
        distance=abs(abs(reference_position-wannier_centres_translated(coord(1),tran_sorted_idx(num_wann-i+1)))-cell_length)
        if (distance .lt. tran_group_threshold) then
            wannier_centres_translated(coord(1),tran_sorted_idx(num_wann-i+1))= &
            wannier_centres_translated(coord(1),tran_sorted_idx(num_wann-i+1))-cell_length
            iterator=iterator+1
        endif
    enddo

    if (iterator .gt. 1) then
        if (iprint .ge. 4) then
            write(stdout,*)' Grouping inconsistency: Restarting sorting'
            write(stdout,*)' '
        endif
        deallocate(PL1_groups)
        deallocate(PL2_groups)
        deallocate(PL3_groups)
        deallocate(PL4_groups)
        goto 100
    endif

    write(stdout,*)' '
    write(stdout,*)'------------------------- Sorted Wannier Centres -----------------------------'
    do i=1,num_wann
        write(stdout,FMT='(A12,I4,3F12.6)')' WF centre  ',tran_sorted_idx(i),wannier_centres_translated(1,tran_sorted_idx(i)), &
                                                          wannier_centres_translated(2,tran_sorted_idx(i)), &
                                                          wannier_centres_translated(3,tran_sorted_idx(i))
    enddo
    write(stdout,*)'------------------------------------------------------------------------------'
    write(stdout,*)' '

    if (timing_level>1) call io_stopwatch('tran: lcr_sort',2)

    return

  end subroutine tran_lcr_sort


  !========================================!
  subroutine master_sort_and_group(Array,Array_groups,coord,Array_size)
    !========================================!

    use w90_constants,          only : dp
    use w90_io,                 only : io_error,stdout,io_stopwatch
    use w90_parameters,         only : one_dim_dir,tran_num_ll,iprint,timing_level
    use w90_hamiltonian,        only : wannier_centres_translated

    implicit none

    integer,intent(in),dimension(:)                 :: Array_groups
    integer,intent(in),dimension(3)                 :: coord
    integer,intent(in)                              :: Array_size 

    real(dp),intent(inout),dimension(2,Array_size)  :: Array

    integer                                         :: i,j,k,Array_num_groups,increment, &
                                                       subgroup_increment,group_num_subgroups
    integer,allocatable,dimension(:)                :: group_subgroups
    
    real(dp),allocatable,dimension(:,:)             :: group_array,sorted_group_array, &
                                                       subgroup_array,sorted_subgroup_array
    character(30)                                   :: fmt_2

    if (timing_level>2) call io_stopwatch('tran: lcr_sort: master_sort',1)
    !
    !Number of groups inside the principal layer
    !
    Array_num_groups=size(Array_groups)                  
    !
    !Convenient variable which will be amended later. Used to appropriately extract the group array from the Array
    !
    increment=1
    !
    !Loop over groups inside Array
    !
    do j=1,Array_num_groups
        allocate(group_array(2,Array_groups(j)))
        allocate(sorted_group_array(2,Array_groups(j)))
        !
        !Extract the group from the Array
        !
        group_array=Array(:,increment:increment+Array_groups(j)-1)
        !
        !Updating group_array to contain coord(2)
        !
        do k=1,Array_groups(j)
            group_array(2,k)=wannier_centres_translated(coord(2),int(group_array(1,k)))
        enddo                                

        call sort(group_array,sorted_group_array)
        call group(sorted_group_array,group_subgroups)

        group_num_subgroups=size(group_subgroups)
        
        if (iprint .ge. 4) then
            !
            !Printing subgroup breakdown    
            !
            write(fmt_2,'(i5)')group_num_subgroups
            fmt_2=adjustl(fmt_2)
            fmt_2='(a7,i3,a1,i5,a2,'//trim(fmt_2)//'i4,a1)'
            write(stdout,fmt_2)' Group ',j,' ',group_num_subgroups,' (',(group_subgroups(i),i=1,group_num_subgroups),')'
        endif
        !
        !Convenient variable which will be amended later. Used to appropriately extract the subgroup array from the group_array
        !
        subgroup_increment=1
        !
        !Loop over subgroups inside group
        !
        do k=1,group_num_subgroups
            allocate(subgroup_array(2,group_subgroups(k)))
            allocate(sorted_subgroup_array(2,group_subgroups(k)))
            !
            !Extract the subgroup from the group
            !
            subgroup_array=sorted_group_array(:,subgroup_increment:subgroup_increment+group_subgroups(k)-1)
            !
            !Updating subgroup_array to contain coord(3)
            !
            do i=1,group_subgroups(k)
                subgroup_array(2,i)=wannier_centres_translated(coord(3),int(subgroup_array(1,i)))
            enddo

            call sort(subgroup_array,sorted_subgroup_array)
            !
            !Update sorted_group array with the sorted subgroup array
            !
            sorted_group_array(:,subgroup_increment:subgroup_increment+group_subgroups(k)-1)=sorted_subgroup_array
            !
            !Update the subgroup_increment
            !
            subgroup_increment=subgroup_increment+group_subgroups(k)
            deallocate(subgroup_array)
            deallocate(sorted_subgroup_array)
        enddo
        !
        !Update Array with the sorted group array
        !
        Array(:,increment:increment+Array_groups(j)-1)=sorted_group_array    
        !
        !Update the group increment
        !
        increment=increment+Array_groups(j)
        deallocate(group_array)
        deallocate(sorted_group_array)
        deallocate(group_subgroups)
    enddo
            
    if (timing_level>2) call io_stopwatch('tran: lcr_sort: master_sort',2)

    return

  end subroutine master_sort_and_group


  !========================================!
  subroutine sort(non_sorted,sorted)
    !========================================!

    use w90_constants,          only : dp
    use w90_io,                 only : io_error

    use w90_hamiltonian,        only : wannier_centres_translated

    implicit none

    real(dp),intent(inout),dimension(:,:)       :: non_sorted
    real(dp),intent(out),dimension(:,:)         :: sorted    
  
    integer,dimension(1)                        :: min_loc
    integer                                     :: num_col,i

    num_col=size(non_sorted,2)

    do i=1,num_col
        !
        !look for the location of the minimum value of the coordinates in non_sorted
        !
        min_loc=minloc(non_sorted(2,:))
        !
        !now the index in the first row of sorted is the index non_sorted(1,min_loc)
        !
        sorted(1,i)=non_sorted(1,min_loc(1))
        !
        !here is the corresponding coordinate 
        !
        sorted(2,i)=non_sorted(2,min_loc(1))
        !
        !here one replaces the minimum coordinate with 10**10 such that this value
        !will not be picked-up again by minloc
        !
        non_sorted(2,min_loc(1))=10**10
    enddo

    return    

    endsubroutine sort

  !========================================!
  subroutine group(array,array_groups)
    !========================================!

    use w90_constants,          only : dp
    use w90_io,                 only : io_error

    use w90_parameters,         only : tran_group_threshold 

    implicit none

    real(dp),intent(in),dimension(:,:)           :: array
    integer,intent(out),allocatable,dimension(:) :: array_groups

    integer,allocatable,dimension(:)             :: dummy_array
    logical,allocatable,dimension(:)             :: logic
    integer                                      :: array_idx,i,j,group_number,array_size


    array_size=size(array,2)

    allocate(dummy_array(array_size))
    allocate(logic(array_size))
    !
    !Initialise dummy array
    !
    dummy_array=0  
    !
    !Initialise logic to false
    !
    logic=.false.
    !
    !Define counter of number of groups  
    !
    array_idx=1
    !
    !Loop over columns of array (ie array_size)
    !
    do i=1,array_size
        !
        !If an element of logic is true then it means the wannier function has already been grouped
        !
        if (logic(i) .eq. .false.) then 
            !
            !Create a group for the wannier function
            !
            logic(i)=.true.
            !
            !Initialise the number of wannier functions in this group to be 1
            !
            group_number=1
            !
            !Loop over the rest of wannier functions in array
            !
            do j=min(i+1,array_size),array_size
                !
                !Special termination cases
                !
                if ((j .eq. 1) .or. (i .eq. array_size)) then
                    dummy_array(array_idx)=group_number
                    exit
                endif
                if (j .eq. array_size .and. (abs(array(2,j)-array(2,i)) .le. tran_group_threshold)) then
                    group_number=group_number+1
                    dummy_array(array_idx)=group_number
                    logic(j)=.true.
                    exit
                endif
                !
                !Check distance between wannier function_i and wannier function_j
                !
                if (abs(array(2,j)-array(2,i)) .le. tran_group_threshold) then
                    !
                    !Increment number of wannier functions in group
                    !
                    group_number=group_number+1
                    !
                    !Assigns wannier function to the group
                    !
                    logic(j)=.true.
                else
                    !
                    !Group is finished and store number of wanniers in the group to dummy_array
                    !
                    dummy_array(array_idx)=group_number
                    !
                    !Increment number of groups
                    !
                    array_idx=array_idx+1
                    exit
                endif
            enddo
        endif
    enddo
    !
    !Copy elements of dummy_array to array_groups
    !
    allocate(array_groups(array_idx))
    array_groups=dummy_array(:array_idx)

    deallocate(dummy_array)
    deallocate(logic)

    return

    endsubroutine group

  !========================================!
  subroutine tran_lcr_build_ham(tran_sorted_idx)
    !========================================!

    use w90_constants,          only : dp
    use w90_io,                 only : io_error,stdout,seedname,io_file_unit,io_date,io_stopwatch
    use w90_parameters,         only : tran_num_cell_ll,num_wann,tran_num_ll,kpt_cart,fermi_energy,&
                                       tran_write_ht,tran_num_rr,tran_num_lc,tran_num_cr,tran_num_cc,&
                                       timing_level
    use w90_hamiltonian,        only : ham_r,wannier_centres_translated
    
    implicit none

    integer,intent(in),dimension(num_wann) :: tran_sorted_idx
    integer                                :: i,j,k,num_wann_cell_ll,file_unit
    
    real(dp),allocatable,dimension(:,:)    :: sub_block

    character(len=9)                       :: cdate, ctime

    if (timing_level>1) call io_stopwatch('tran: lcr_build_ham',1)

    allocate(hL0(tran_num_ll,tran_num_ll))
    allocate(hL1(tran_num_ll,tran_num_ll))
    allocate(hR0(tran_num_ll,tran_num_ll))
    allocate(hR1(tran_num_ll,tran_num_ll))
    allocate(hC(num_wann-(2*tran_num_ll),num_wann-(2*tran_num_ll)))
    allocate(hLC(tran_num_ll,tran_num_ll))
    allocate(hCR(tran_num_ll,tran_num_ll))    
    !
    !This checks that only the gamma point is used in wannierisation
    !This is necessary since this calculation only makes sense if we
    !have periodicity over the supercell.
    !
    if ((size(kpt_cart,2) .ne. 1 ) .and. (kpt_cart(1,1) .eq. 0.d0) &
                                   .and. (kpt_cart(2,1) .eq. 0.d0) &
                                   .and. (kpt_cart(3,1) .eq. 0.d0) ) then 
        call io_error('Calculation must be performed at gamma only')
    endif
    

    num_wann_cell_ll=tran_num_ll/tran_num_cell_ll

    allocate(sub_block(num_wann_cell_ll,num_wann_cell_ll))
    !
    !Build hL0 & hL1
    !
    hL0=0.d0
    hL1=0.d0
    !    
    !Loop over the sub_blocks corresponding to distinct unit cells inside the principal layer
    !
    do i=1,tran_num_cell_ll
        !
        !Each sub_block will be duplicated along the corresponding diagonal. This ensures the correct symmetry for the leads.
        !
        sub_block=0.d0
        !
        !Extract matrix elements from ham_r needed for hL0 (and all but the diagonal sub_block of hL1)
        !
        do j=1,num_wann_cell_ll
            do k=1,num_wann_cell_ll
                sub_block(j,k)=real(ham_r(tran_sorted_idx(j),tran_sorted_idx((i-1)*num_wann_cell_ll+k),1),dp)
            enddo
        enddo
        !
        !Filling up hL0 sub_block by sub_block
        !
        do j=1,tran_num_cell_ll-i+1
            !
            !Fill diagonal and upper diagonal sub_blocks
            !
            hL0((j-1)*num_wann_cell_ll+1                        :j*num_wann_cell_ll,&
                (j-1)*num_wann_cell_ll+1+(i-1)*num_wann_cell_ll :j*num_wann_cell_ll+(i-1)*num_wann_cell_ll)=sub_block
            !
            !Fill lower diagonal sub_blocks
            !
            if (i .gt. 1) then
                hL0((j-1)*num_wann_cell_ll+1+(i-1)*num_wann_cell_ll :j*num_wann_cell_ll+(i-1)*num_wann_cell_ll,&
                    (j-1)*num_wann_cell_ll+1                        :j*num_wann_cell_ll)=transpose(sub_block)
            endif
        enddo
        !
        !Filling up non-diagonal hL1 sub_blocks (nothing need be done for i=1)
        !
        if (i .gt. 1) then
            do j=1,i-1
                hL1((tran_num_cell_ll-(i-j))*num_wann_cell_ll+1 : (tran_num_cell_ll-(i-1-j))*num_wann_cell_ll,&
                    (j-1)*num_wann_cell_ll+1                    : j*num_wann_cell_ll)=sub_block                
            enddo
        endif
    enddo
    !
    !Special case tran_num_cell_ll=1, the diagonal sub-block of hL1 is hL1, so cannot be left as zero
    !
    if (tran_num_cell_ll .eq. 1) then
        do j=1,num_wann_cell_ll
            do k=num_wann_cell_ll+1,2*num_wann_cell_ll
                hL1(j,k)=real(ham_r(tran_sorted_idx(j),tran_sorted_idx(k),1),dp)
            enddo
        enddo
    endif

    !
    !Build hR0 & hR1
    !
    hR0=0.d0
    hR1=0.d0
    !
    !Loop over the sub_blocks corresponding to distinct unit cells inside the principal layer
    !
    do i=1,tran_num_cell_ll
        !
        !Each sub_block will be duplicated along the corresponding diagonal. This ensures the correct symmetry for the leads.
        !
        sub_block=0.d0
        !
        !Extract matrix elements from ham_r needed for hR0 (and all but the diagonal sub_block of hR1)
        !
        do j=1,num_wann_cell_ll
            do k=1,num_wann_cell_ll
                sub_block(j,k)=real(ham_r(tran_sorted_idx(num_wann-i*(num_wann_cell_ll)+j),tran_sorted_idx(num_wann-num_wann_cell_ll+k),1),dp)
            enddo
        enddo
        !
        !Filling up hR0 sub_block by sub_block
        !
        do j=1,tran_num_cell_ll-i+1
            !
            !Fill diagonal and upper diagonal sub_blocks
            !
            hR0((j-1)*num_wann_cell_ll+1                        :j*num_wann_cell_ll,&
                (j-1)*num_wann_cell_ll+1+(i-1)*num_wann_cell_ll :j*num_wann_cell_ll+(i-1)*num_wann_cell_ll)=sub_block
            !
            !Fill lower diagonal sub_blocks
            !
            if (i .gt. 1) then
                hR0((j-1)*num_wann_cell_ll+1+(i-1)*num_wann_cell_ll :j*num_wann_cell_ll+(i-1)*num_wann_cell_ll,&
                    (j-1)*num_wann_cell_ll+1                        :j*num_wann_cell_ll)=transpose(sub_block)
            endif
        enddo
        !
        !Filling up non-diagonal hR1 sub_blocks (nothing need be done for i=1)
        !
        if (i .gt. 1) then
            do j=1,i-1
                hR1((tran_num_cell_ll-(i-j))*num_wann_cell_ll+1 : (tran_num_cell_ll-(i-1-j))*num_wann_cell_ll,&
                    (j-1)*num_wann_cell_ll+1                    : j*num_wann_cell_ll)=sub_block
            enddo
        endif
    enddo
    !
    !Special case tran_num_cell_ll=1, the diagonal sub-block of hR1 is hR1, so cannot be left as zero
    !
    if (tran_num_cell_ll .eq. 1) then
        do j=1,num_wann_cell_ll
            do k=num_wann_cell_ll+1,2*num_wann_cell_ll
                hR1(j,k)=real(ham_r(tran_sorted_idx(j),tran_sorted_idx(k),1),dp)
            enddo
        enddo
    endif

    !
    !Building hLC
    !
    hLC=0.d0
    do i=1,tran_num_ll
        do j=tran_num_ll+1,2*tran_num_ll
            hLC(i,j-tran_num_ll)=real(ham_r(tran_sorted_idx(i),tran_sorted_idx(j),1),dp)
        enddo
    enddo
    if (tran_num_cell_ll .gt. 1) then
        do j=1,tran_num_cell_ll
            do k=1,tran_num_cell_ll
                if (k .ge. j) then
                    hLC((j-1)*num_wann_cell_ll+1:j*num_wann_cell_ll,(k-1)*num_wann_cell_ll+1:k*num_wann_cell_ll)=0.d0
                endif
            enddo
        enddo
    endif
    !
    !Building hC
    !
    hC=0.d0
    do i=tran_num_ll+1,num_wann-tran_num_ll
        do j=tran_num_ll+1,num_wann-tran_num_ll
            hC(i-tran_num_ll,j-tran_num_ll)=real(ham_r(tran_sorted_idx(i),tran_sorted_idx(j),1),dp)
        enddo
    enddo
    !
    !Building hCR
    !
    hCR=0.d0
    do i=num_wann-2*tran_num_ll+1,num_wann-tran_num_ll
        do j=num_wann-tran_num_ll+1,num_wann
            hCR(i-(num_wann-2*tran_num_ll),j-(num_wann-tran_num_ll))=real(ham_r(tran_sorted_idx(i),tran_sorted_idx(j),1),dp)
        enddo
    enddo
    if (tran_num_cell_ll .gt. 1) then
        do j=1,tran_num_cell_ll
            do k=1,tran_num_cell_ll
                if (k .ge. j) then
                    hCR((j-1)*num_wann_cell_ll+1:j*num_wann_cell_ll,(k-1)*num_wann_cell_ll+1:k*num_wann_cell_ll)=0.d0
                endif
            enddo
        enddo
    endif
    !
    !Subtract the Fermi energy from the diagonal elements of hC,hL0,hR0
    !
    do i=1,tran_num_ll
        hL0(i,i)=hL0(i,i)-fermi_energy
        hR0(i,i)=hR0(i,i)-fermi_energy
    enddo
    do i=1,num_wann-(2*tran_num_ll)
        hC(i,i)=hC(i,i)-fermi_energy
    enddo
    !
    !Define tran_num_** parameters that are used later in tran_lcr
    !
    tran_num_rr=tran_num_ll
    tran_num_lc=tran_num_ll
    tran_num_cr=tran_num_ll
    tran_num_cc=num_wann-(2*tran_num_ll)

    deallocate(sub_block)
    !
    !Writing to file:
    !
    if ( tran_write_ht ) then
        write(stdout,*)'------------------------------- Writing ht files  ----------------------------'  
        !
        file_unit = io_file_unit()
        open(file_unit,file=trim(seedname)//'_htL.dat',status='unknown',form='formatted',action='write')

        call io_date(cdate,ctime)
        write(file_unit,*) 'written on '//cdate//' at '//ctime ! Date and time
        write(file_unit,'(I6)') tran_num_ll
        write(file_unit,'(6F12.6)') ((hL0(j,i),j=1,tran_num_ll),i=1,tran_num_ll)
        write(file_unit,'(I6)') tran_num_ll
        write(file_unit,'(6F12.6)') ((hL1(j,i),j=1,tran_num_ll),i=1,tran_num_ll)

        close(file_unit)
        write(stdout,*)' '//trim(seedname)//'_htL.dat  written'
        !
        !hR
        !
        file_unit = io_file_unit()
        open(file_unit,file=trim(seedname)//'_htR.dat',status='unknown',form='formatted',action='write')

        call io_date(cdate,ctime)
        write(file_unit,*) 'written on '//cdate//' at '//ctime ! Date and time
        write(file_unit,'(I6)') tran_num_rr
        write(file_unit,'(6F12.6)') ((hR0(j,i),j=1,tran_num_rr),i=1,tran_num_rr)
        write(file_unit,'(I6)') tran_num_rr
        write(file_unit,'(6F12.6)') ((hR1(j,i),j=1,tran_num_rr),i=1,tran_num_rr)

        close(file_unit)
        write(stdout,*)' '//trim(seedname)//'_htR.dat  written'
        !
        !hLC
        !
        file_unit = io_file_unit()
        open(file_unit,file=trim(seedname)//'_htLC.dat',status='unknown',form='formatted',action='write')

        call io_date(cdate,ctime)
        write(file_unit,*) 'written on '//cdate//' at '//ctime ! Date and time
        write(file_unit,'(2I6)') tran_num_ll,tran_num_lc
        write(file_unit,'(6F12.6)') ((hLC(j,i),j=1,tran_num_lc),i=1,tran_num_lc)

        close(file_unit)
        write(stdout,*)' '//trim(seedname)//'_htLC.dat written'
        !
        !hCR
        !
        file_unit = io_file_unit()
        open(file_unit,file=trim(seedname)//'_htCR.dat',status='unknown',form='formatted',action='write')

        call io_date(cdate,ctime)
        write(file_unit,*) 'written on '//cdate//' at '//ctime ! Date and time
        write(file_unit,'(2I6)') tran_num_cr,tran_num_rr
        write(file_unit,'(6F12.6)') ((hCR(j,i),j=1,tran_num_cr),i=1,tran_num_cr)

        close(file_unit)
        write(stdout,*)' '//trim(seedname)//'_htCR.dat written'
        !
        !hC
        !
        file_unit = io_file_unit()
        open(file_unit,file=trim(seedname)//'_htC.dat',status='unknown',form='formatted',action='write')

        call io_date(cdate,ctime)
        write(file_unit,*) 'written on '//cdate//' at '//ctime ! Date and time
        write(file_unit,'(I6)') tran_num_cc
        write(file_unit,'(6F12.6)') ((hC(j,i),j=1,tran_num_cc),i=1,tran_num_cc)

        close(file_unit)
        write(stdout,*)' '//trim(seedname)//'_htC.dat  written'

        write(stdout,*)'------------------------------------------------------------------------------'
    end if

    if (timing_level>1) call io_stopwatch('tran: lcr_build_ham',2)

    return 

  end subroutine tran_lcr_build_ham


!=============================================
 subroutine tran_hr_parity_shift(tran_sorted_idx)
  !==============================================
  ! Find parity of each MLWF and set all to be positive
  ! for consistency from unit cell to unit cell within 
  ! the PLs.
  ! The parity of the MLWFs are calculated by summing 
  ! their values over all space, assigning positive parity
  ! to those which have a positive sum. (Negative for
  ! negative sum)
  !
  ! Also further work is needed for mulitple MLWF with
  ! the same centre. These must be consitently ordered 
  ! using the sequence of signs of the higher order 
  ! integrals as a signiture.
  ! 
  ! M. Shelley 03/09


    use w90_constants,          only : dp,cmplx_0,twopi,pi
    use w90_io,                 only : io_error,stdout,seedname,io_file_unit,io_date,&
                                       io_stopwatch
    use w90_parameters,         only : num_wann,spin,wvfn_formatted,have_disentangled,&
                                       ndimwin,lwindow,num_bands,u_matrix,u_matrix_opt,&
                                       real_lattice,tran_wf_threshold,tran_num_ll,&
                                       tran_num_cell_ll,iprint,timing_level
    use w90_hamiltonian,        only : ham_r,wannier_centres_translated
    
    implicit none

    integer,intent(in),dimension(num_wann)            :: tran_sorted_idx
    complex(kind=dp), allocatable                     :: r_wvfn(:,:),uu(:),tran_u_matrix(:,:)
    real(kind=dp)                                     :: w_real,w_imag,cell_len(3)
    real(kind=dp),allocatable                         :: r(:,:)
    real(dp),dimension(tran_num_ll/tran_num_cell_ll)  :: parity_identifier
    character(len=11)                                 :: wfnname
    character(len=30)                                 :: i_char
    logical                                           :: have_file,inc_band(num_bands),&
                                                         ham_r_shift(num_wann),integral_found(num_wann)
    integer                                           :: file_unit,n,m,p,nx,ny,nz,counter,&
                                                         nk,nbnd,ngx,ngy,ngz,ierr,&
                                                         integral_counter,order,n1,n2,n3,npoint,&
                                                         i,j,num_wann_cell_ll,p_max

    if (timing_level>1) call io_stopwatch('tran: hr_parity_shift',1)

    write(wfnname,200 ) 1,spin
    inquire(file=wfnname,exist=have_file)
    if(.not.have_file) call io_error('tran_hr_parity_shift: file '//wfnname//' not found')

    file_unit=io_file_unit()
    if(wvfn_formatted) then
       open(unit=file_unit,file=wfnname,form='formatted')
       read(file_unit,*) ngx,ngy,ngz,nk,nbnd
    else
       open(unit=file_unit,file=wfnname,form='unformatted')
       read(file_unit) ngx,ngy,ngz,nk,nbnd
    end if

200 format ('UNK',i5.5,'.',i1) 

    if(have_disentangled) then
        allocate(r_wvfn(ngx*ngy*ngz,num_bands),stat=ierr )
        if (ierr/=0) call io_error('Error in allocating r_wvfn in tran_hr_parity_shift')
        allocate(tran_u_matrix(num_bands,num_wann), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating tran_u_matrix in tran_hr_parity_shift')
    else
        allocate(r_wvfn(ngx*ngy*ngz,num_wann),stat=ierr )
        if (ierr/=0) call io_error('Error in allocating r_wvfn in tran_hr_parity_shift')
        allocate(tran_u_matrix(num_wann,num_wann), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating tran_u_matrix in tran_hr_parity_shift')
    endif

    allocate(uu(num_wann),stat=ierr )
    if (ierr/=0) call io_error('Error in allocating uu in tran_hr_parity_shift')

    allocate(r(ngx*ngy*ngz,3),stat=ierr )
    if (ierr/=0) call io_error('Error in allocating r in tran_hr_parity_shift')

    !
    !Create logical array set to true if band is in the outer window at Gamma. 
    !
    inc_band=.true.
    if(have_disentangled) then
        inc_band(:)=lwindow(:,1)
    end if
    !
    !Read U_nk file
    !
    write(stdout,*)'    Reading '//wfnname//' file'
    if(have_disentangled) then
        counter=1
        do n=1,num_bands
            if(counter>ndimwin(1)) exit
            if(wvfn_formatted) then
                do nx=1,ngx*ngy*ngz
                    read(file_unit,*) w_real, w_imag
                    r_wvfn(nx,counter) = cmplx(w_real,w_imag,kind=dp)
                end do
            else
                read(file_unit) (r_wvfn(nx,counter),nx=1,ngx*ngy*ngz)
            end if
            if(inc_band(n)) counter=counter+1
        end do
    else
        do n=1,num_wann
           if(wvfn_formatted) then
                do nx=1,ngx*ngy*ngz
                    read(file_unit,*) w_real, w_imag
                    r_wvfn(nx,n) = cmplx(w_real,w_imag,kind=dp)
                end do
            else
                read(file_unit) (r_wvfn(nx,n),nx=1,ngx*ngy*ngz)
            end if
        end do
    end if
    close(file_unit)

    tran_u_matrix=cmplx_0
    uu=cmplx_0
    ham_r_shift=.false.
    integral_found=.false.

    do i=1,3
        cell_len(i)=dsqrt(real_lattice(i,1)**2+real_lattice(i,2)**2+real_lattice(i,3)**2)
    enddo


    write(stdout,*)' Calculating parities of MLWFs'
    write(stdout,*)' '

    !
    !Loop over wannier functions
    !
    write(stdout,*)'WF       integral_counter       integral value'
    do n=1,num_wann
        !
        !Disentanglement step
        !
        if(have_disentangled) then
            do p=1,num_bands
                do m=1,num_wann
                   tran_u_matrix(p,n)=tran_u_matrix(p,n)+u_matrix(m,n,1)*u_matrix_opt(p,m,1)
                enddo
            enddo
            p_max=num_bands
        else
            tran_u_matrix=u_matrix(:,:,1)
            p_max=num_wann
        endif
        !
        !Calculate uu for 0th order
        !      
        do p=1,p_max
            uu(n)=uu(n)+tran_u_matrix(p,n)*sum(r_wvfn(:,p))
        enddo
        !
        !Normalise
        !
        uu(n)=uu(n)/(ngx*ngy*ngz)
        if (iprint .ge. 4) write(stdout,*)n,'0',real(uu(n))
        !
        !Integral not significantly away from zero, move to higher orders
        !
        if (abs(uu(n)) .lt. tran_wf_threshold) then
            order=1
            integral_counter=1
            !
            !Define vectors from WF centre to each grid point
            !
            do nz=1,ngz
                do ny=1,ngy
                    do nx=1,ngx
                        npoint=nx+(ny-1)*ngx+(nz-1)*ngx*ngy
                        do i=1,3
                            r(npoint,i)=(nx-1)*real_lattice(1,i)/ngx &
                                       +(ny-1)*real_lattice(2,i)/ngy &
                                       +(nz-1)*real_lattice(3,i)/ngz-wannier_centres_translated(i,n)
                        enddo
                    enddo
                enddo
            enddo
            !
            !Loop over powers of integral pre-factor g(:)=(sin(r(:,1)*2*pi/cell_len(1))^n1)&
            !                                            *(sin(r(:,2)*2*pi/cell_len(2))^n2)&
            !                                            *(sin(r(:,3)*2*pi/cell_len(3))^n3)
            !The sin's enforce periodicity of the supercell
            !
444 continue
            do n1=order,0,-1
                do n2=order,0,-1
                    do n3=order,0,-1
                        if ( (n1+n2+n3 .eq. order) .and. .not. integral_found(n) ) then
                            do p=1,p_max
                                !
                                !Calculate double sum with prefactor (disentanglement step already performed at order=0)
                                !
                                uu(n)=uu(n)+tran_u_matrix(p,n)*sum((dsin(r(:,1)*2.d0*pi/cell_len(1))**n1)&
                                                                  *(dsin(r(:,2)*2.d0*pi/cell_len(2))**n2)&
                                                                  *(dsin(r(:,3)*2.d0*pi/cell_len(3))**n3)*r_wvfn(:,p))
                            enddo
                            !
                            !Normalise
                            !
                            uu(n)=uu(n)/(ngx*ngy*ngz)
                            if (iprint .ge. 4) write(stdout,*)n,integral_counter,real(uu(n))
                            !
                            !Check if integral is sufficiently far from zero
                            !
                            if (abs(uu(n)) .lt. tran_wf_threshold) then
                                integral_counter=integral_counter+1
                                !
                                !Define cases which increase order
                                !
                                select case(integral_counter)
                                case(4)
                                    order=2
                                    goto 444
                                case(10)
                                    order=3
                                    goto 444
                                case(20) 
                                    call io_error('tran_hr_parity_unkg: Parity integral not found')
                                endselect
                            else
                                integral_found(n)=.true.
                            endif
                        endif            
                    enddo
                enddo
            enddo
        endif !order > 0 if
        !
        !Save parity
        !
        if (real(uu(n)) .lt. 0) ham_r_shift(n)=.true.
    enddo !Wannier function loop

    !
    !Set all negative parity to be positive
    !
    do n=1,num_wann
        do m=1,num_wann
            if (ham_r_shift(n)) then
                ham_r(n,m,1)=-ham_r(n,m,1)
                ham_r(m,n,1)=-ham_r(m,n,1)
            endif
        enddo
    enddo
    write(stdout,*)' Consistent parities enforced throughout PL1, PL2, PL3 and PL4'

    !
    !Quick parity double check through PL1 & 2
    !
    num_wann_cell_ll=tran_num_ll/tran_num_cell_ll
!    do i=1,2*tran_num_cell_ll
!        parity_identifier=0
!        do j=1,num_wann_cell_ll
!            parity_identifier(j)=real(ham_r(tran_sorted_idx(1),tran_sorted_idx(j),1))*&
!                                 real(ham_r(tran_sorted_idx((i-1)*num_wann_cell_ll+1),&
!                                            tran_sorted_idx((i-1)*num_wann_cell_ll+j),1))
            !
            !if parity_identifier(j) < 0 parities must be incorrect
            !
!            if (parity_identifier(j) .lt. 0 ) then
!                write(i_char,*)i
!                call io_error('Parities found to be inconsistent in unit cell '//trim(i_char)//' from left hand side')
!            endif
!        enddo
!    enddo
    !
    !PL3 & 4
    !
!    do i=1,2*tran_num_cell_ll
!        parity_identifier=0
!        do j=1,num_wann_cell_ll
!            parity_identifier(j)=real(ham_r(tran_sorted_idx(1),tran_sorted_idx(j),1))*&
!                                 real(ham_r(tran_sorted_idx(num_wann-(2*tran_num_ll)+(i-1)*num_wann_cell_ll+1),&
!                                            tran_sorted_idx(num_wann-(2*tran_num_ll)+(i-1)*num_wann_cell_ll+j),1))
            !
            !if parity_identifier(j) < 0 parities must be incorrect
            !
!            if (parity_identifier(j) .lt. 0 ) then
!                write(i_char,*)2*tran_num_cell_ll-i+1
!                call io_error('Parities found to be inconsistent in unit cell '//trim(i_char)//' from right hand side')
!            endif
!        enddo
!    enddo


    deallocate(r_wvfn,stat=ierr )
    if (ierr/=0) call io_error('Error deallocating r_wvfn in tran_hr_parity_shift')
    deallocate(uu,stat=ierr )
    if (ierr/=0) call io_error('Error deallocating uu in tran_hr_parity_shift')
    deallocate(r,stat=ierr )
    if (ierr/=0) call io_error('Error deallocating r in tran_hr_parity_shift')
    deallocate(tran_u_matrix,stat=ierr )
    if (ierr/=0) call io_error('Error deallocating tran_u_matrix in tran_hr_parity_shift')


    if (timing_level>1) call io_stopwatch('tran: hr_parity_shift',2)

    return

 end subroutine tran_hr_parity_shift


end module w90_transport
