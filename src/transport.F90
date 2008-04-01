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

!=====================================================================!
! Definition of parameters used in w90_transport                      !
!=====================================================================!
!                                                                     !
! transport_mode     = 'bulk' or 'lcr'                                !
! tran_win_min   = minimum E                                          !
! tran_win_max   = maximum E                                          !
! tran_energy_step   = delta E                                        !
! tran_num_bb = # of WFs in a principal layer of a perfectly          !
!                 periodic bulk system                                !
! tran_num_ll = # of WFs in a principal layer of a left lead          !            
! tran_num_rr = # of WFs in a principal layer of a right lead         !             
! tran_num_cc = # of WFs in a disordered conducter cell               !       
! tran_num_lc = # of WFs in a disordered conducter cell that are used !
!                 to calculate interaction with a left lead           !
! tran_num_cr = # of WFs in a disordered conducter cell that are used !
!                 to calculate interaction with a right lead          !
! tran_num_bandc     = width of band-diagonal hC matrix               !
! tran_read_ht       = .true. => read H matrix from h*.dat files      !
! tran_write_ht      = .true. => write H matrix from h*.dat files     !
! tran_use_same_lead = .true. => in L-C-R construction, left and      !
!                                right lead are the same kind         !
!                                                                     !
!=====================================================================!

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
    use w90_parameters, only : transport_mode,tran_read_ht,timing_level,hr_plot
    use w90_hamiltonian,only : hamiltonian_get_hr,hamiltonian_write_hr,hamiltonian_setup

    implicit none
 
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
                               tran_use_same_lead, timing_level

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

end module w90_transport
