!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2007-13 Jonathan Yates, Arash Mostofi,       !
!                Giovanni Pizzi, Young-Su Lee,               !
!                Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

module w90_plot

  implicit none

  private

  public :: plot_main


contains

  !============================================!
  subroutine plot_main( )
    !============================================!

    use w90_constants,   only : eps6
    use w90_io,          only : stdout,io_stopwatch
    use w90_parameters,  only : num_kpts,bands_plot,dos_plot,hr_plot, &
                                kpt_latt,fermi_surface_plot, &
                                wannier_plot,timing_level
    use w90_hamiltonian, only : hamiltonian_get_hr,hamiltonian_write_hr, &
                                hamiltonian_setup

    implicit none

    integer :: nkp
    logical :: have_gamma

    if (timing_level>0) call io_stopwatch('plot: main',1)

    write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,'(1x,a)') '|                               PLOTTING                                    |'
    write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,*)

    if(bands_plot .or. dos_plot .or. fermi_surface_plot .or. hr_plot) then
       ! Check if the kmesh includes the gamma point
       have_gamma=.false.
       do nkp=1,num_kpts
           if (all(abs(kpt_latt(:,nkp))<eps6)) have_gamma=.true.       
       end do
       if(.not. have_gamma) &
            write(stdout,'(1x,a)') '!!!! Kpoint grid does not include Gamma. ' // &
            & ' Interpolation may be incorrect. !!!!'
       ! Transform Hamiltonian to WF basis
       !
       call hamiltonian_setup()
       !
       call hamiltonian_get_hr()
       !
       if(bands_plot) call plot_interpolate_bands
       !
       if(fermi_surface_plot) call plot_fermi_surface
       !
       if(hr_plot) call hamiltonian_write_hr()
       !
    end if

    if(wannier_plot) call plot_wannier

    if (timing_level>0) call io_stopwatch('plot: main',2)


  end subroutine plot_main



  !-----------------------------------!
  !----------------- Private Routines ------------------------------
  !-----------------------------------!



  !============================================!
  subroutine plot_interpolate_bands
    !============================================!
    !                                            !
    !     Plots the interpolated band structure  !
    !                                            !
    !============================================!

    use w90_constants,  only  : dp,cmplx_0,cmplx_i,twopi
    use w90_io,         only  : io_error,stdout,io_file_unit,seedname,&
                                io_time,io_stopwatch
    use w90_parameters, only  : num_wann,bands_num_points,recip_metric,&
                                bands_num_spec_points,timing_level, &
                                bands_spec_points,bands_label,bands_plot_format, &
                                bands_plot_mode,num_bands_project,bands_plot_project
    use w90_hamiltonian, only : irvec,nrpts,ndegen,ham_r

    implicit none

    complex(kind=dp),allocatable  :: ham_r_cut(:,:,:)
    complex(kind=dp),allocatable  :: ham_pack(:)
    complex(kind=dp)              :: fac
    complex(kind=dp),allocatable  :: ham_kprm(:,:)
    complex(kind=dp),allocatable  :: U_int(:,:)
    complex(kind=dp),allocatable  :: cwork(:)     
    real(kind=dp),allocatable     :: rwork(:)     
    real(kind=dp)      :: kpath_len(bands_num_spec_points/2)     
    integer            :: kpath_pts(bands_num_spec_points/2)     
    real(kind=dp), allocatable :: xval(:)
    real(kind=dp), allocatable :: eig_int(:,:),plot_kpoint(:,:)
    real(kind=dp), allocatable :: bands_proj(:,:)
    real(kind=dp)        :: rdotk,vec(3),emin,emax,time0
    integer, allocatable :: irvec_cut(:,:)
    integer              :: irvec_max(3)
    integer              :: nrpts_cut
    integer, allocatable :: iwork(:),ifail(:)
    integer              :: info,i,j
    integer              :: irpt,nfound,loop_kpt,counter
    integer              :: loop_spts,total_pts,loop_i,nkp
    integer              :: num_paths,num_spts,ierr
    integer              :: bndunit,gnuunit,loop_w,loop_p
    character(len=3),allocatable   :: glabel(:)
    character(len=10),allocatable  :: xlabel(:)
    character(len=10),allocatable  :: ctemp(:)
    !
    if (timing_level>1) call io_stopwatch('plot: interpolate_bands',1)
    !
    time0=io_time()
    write(stdout,*) 
    write(stdout,'(1x,a)') 'Calculating interpolated band-structure'
    write(stdout,*) 
    !
    allocate(ham_pack((num_wann*(num_wann+1))/2),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ham_pack in plot_interpolate_bands')
    allocate(ham_kprm(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ham_kprm in plot_interpolate_bands')
    allocate(U_int(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating U_int in plot_interpolate_bands')
    allocate(cwork(2*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cwork in plot_interpolate_bands')
    allocate(rwork(7*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating rwork in plot_interpolate_bands')
    allocate(iwork(5*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating iwork in plot_interpolate_bands')
    allocate(ifail(num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ifail in plot_interpolate_bands')
    !
    ! Work out how many points in the total path and the positions of the special points
    !
    num_paths=bands_num_spec_points/2
    num_spts=num_paths+1
    do loop_spts=1,num_paths
       vec=bands_spec_points(:,2*loop_spts)-bands_spec_points(:,2*loop_spts-1)
       kpath_len(loop_spts)=sqrt(dot_product(vec,(matmul(recip_metric,vec))))
       if(loop_spts==1) then
          kpath_pts(loop_spts)=bands_num_points
       else 
          kpath_pts(loop_spts)=nint(real(bands_num_points,dp)*kpath_len(loop_spts)/kpath_len(1))
       end if
    end do
    total_pts=sum(kpath_pts)+1
    allocate(plot_kpoint(3,total_pts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating plot_kpoint in plot_interpolate_bands')
    allocate(xval(total_pts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating xval in plot_interpolate_bands')
    allocate(eig_int(num_wann,total_pts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating eig_int in plot_interpolate_bands')
    allocate(bands_proj(num_wann,total_pts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating bands_proj in plot_interpolate_bands')
    allocate(glabel(num_spts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating num_spts in plot_interpolate_bands')
    allocate(xlabel(num_spts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating xlabel in plot_interpolate_bands')
    allocate(ctemp(bands_num_spec_points),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ctemp in plot_interpolate_bands')
    eig_int=0.0_dp;bands_proj=0.0_dp
    !
    ! Find the position of each kpoint in the path
    !
    counter=0
    do loop_spts=1,num_paths
       do loop_i=1,kpath_pts(loop_spts)
          counter=counter+1
          if(counter==1) then
             xval(counter)=0.0_dp
          else
             xval(counter)=xval(counter-1)+kpath_len(loop_spts)/real(kpath_pts(loop_spts),dp)
          endif
          plot_kpoint(:,counter)=bands_spec_points(:,2*loop_spts-1)+ &
               (bands_spec_points(:,2*loop_spts)-bands_spec_points(:,2*loop_spts-1))* &
               (real(loop_i-1,dp)/real(kpath_pts(loop_spts),dp))
       end do
    end do
    xval(total_pts)=sum(kpath_len)
    plot_kpoint(:,total_pts)=bands_spec_points(:,bands_num_spec_points)
    !
    ! Write out the kpoints in the path
    !
    bndunit=io_file_unit()
    open(bndunit,file=trim(seedname)//'_band.kpt',form='formatted')
    write(bndunit,*) total_pts
    do loop_spts=1,total_pts
       write(bndunit,'(3f12.6,3x,a)') (plot_kpoint(loop_i,loop_spts),loop_i=1,3),"1.0"
    end do
    close(bndunit)
    !
    ! Cut H matrix in real-space
    !
    if (index(bands_plot_mode,'cut').ne.0)  call plot_cut_hr()
    !
    ! Interpolate the Hamiltonian at each kpoint
    !
    do loop_kpt=1,total_pts
       ham_kprm=cmplx_0
       if (index(bands_plot_mode,'s-k').ne.0) then
          do irpt=1,nrpts
             rdotk=twopi*dot_product(plot_kpoint(:,loop_kpt),irvec(:,irpt))
             fac=exp(cmplx_i*rdotk)/real(ndegen(irpt),dp)
             ham_kprm=ham_kprm+fac*ham_r(:,:,irpt)
          end do
       elseif (index(bands_plot_mode,'cut').ne.0) then
          do irpt=1,nrpts_cut
             rdotk=twopi*dot_product(plot_kpoint(:,loop_kpt),irvec_cut(:,irpt))
!!$[aam] check divide by ndegen?
             fac=exp(cmplx_i*rdotk)
             ham_kprm=ham_kprm+fac*ham_r_cut(:,:,irpt)
          end do
       else
          call io_error('Error in plot_interpolate bands: value of bands_plot_mode not recognised')
       endif
       ! Diagonalise H_k (->basis of eigenstates)
       do j=1,num_wann
          do i=1,j
             ham_pack(i+((j-1)*j)/2)=ham_kprm(i,j)
          enddo
       enddo
       call ZHPEVX('V','A','U',num_wann,ham_pack,0.0_dp,0.0_dp,0,0,-1.0_dp, &
            nfound,eig_int(1,loop_kpt),U_int,num_wann,cwork,rwork,iwork,ifail,info)
       if(info < 0) then
          write(stdout,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
          call io_error('Error in plot_interpolate_bands')
       endif
       if(info > 0) then
          write(stdout,'(i3,a)') info,' EIGENVECTORS FAILED TO CONVERGE'
          call io_error('Error in plot_interpolate_bands')
       endif
       ! Compute projection onto WF if requested
       if(num_bands_project>0) then
       do loop_w=1,num_wann
         do loop_p=1,num_wann 
          if(any(bands_plot_project==loop_p)) then
            bands_proj(loop_w,loop_kpt)=bands_proj(loop_w,loop_kpt)+abs(U_int(loop_p,loop_w))**2
          end if 
         end do
       end do
       end if
       !
    end do
    !
    ! Interpolation Finished!
    ! Now we write plotting files
    !
    emin=minval(eig_int)-1.0_dp
    emax=maxval(eig_int)+1.0_dp


    if(index(bands_plot_format,'gnu')>0)  call plot_interpolate_gnuplot
    if(index(bands_plot_format,'xmgr')>0) call plot_interpolate_xmgrace


    write(stdout,'(1x,a,f11.3,a)')  &
         'Time to calculate interpolated band structure ',io_time()-time0,' (sec)'
    write(stdout,*)

    if (allocated(ham_r_cut)) deallocate(ham_r_cut,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating ham_r_cut in plot_interpolate_bands')
    if (allocated(irvec_cut)) deallocate(irvec_cut,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating irvec_cut in plot_interpolate_bands')
    !
    if (timing_level>1) call io_stopwatch('plot: interpolate_bands',2)
    !

    contains

  !============================================!
  subroutine plot_cut_hr() 
    !============================================!
    !                                                           
    !  In real-space picture, ham_r(j,i,k) is an interaction between                   
    !  j_th WF at 0 and i_th WF at the lattice point translated 
    !  by matmul(real_lattice(:,:),irvec(:,k))                          
    !  We truncate Hamiltonian matrix when                                                               
    !   1) |  r_i(0) - r_j (R) | > dist_cutoff               
    !   2) |  ham_r(i,j,k)     | < hr_cutoff              
    !  while the condition 1) is essential to get a meaningful band structure, 
    !    ( dist_cutoff must be smaller than the shortest distance from
    !      the center of W-S supercell to the points at the cell boundaries )
    !  the condition 2) is optional.
    !      
    !  limitation: when bands_plot_dim .ne. 3
    !      one_dim_vec must be parallel to one of the cartesian axis
    !      and perpendicular to the other two primitive lattice vectors                         
    !============================================!

    use w90_constants,   only : dp,cmplx_0, eps8
    use w90_io,          only : io_error,stdout
    use w90_parameters,  only : num_wann, mp_grid, real_lattice,   &
                                one_dim_dir, bands_plot_dim,       &
                                hr_cutoff, dist_cutoff, dist_cutoff_mode
    use w90_hamiltonian, only : wannier_centres_translated

    implicit none
    !
    integer :: nrpts_tmp
    integer :: one_dim_vec, two_dim_vec(2)
    integer :: i, j, n1, n2, n3, i1, i2, i3
    real(kind=dp), allocatable :: ham_r_tmp(:,:)
    real(kind=dp), allocatable :: shift_vec(:,:)
    real(kind=dp) :: dist_ij_vec(3)
    real(kind=dp) :: dist_vec(3)
    real(kind=dp) :: dist

    allocate(ham_r_tmp(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ham_r_tmp in plot_cut_hr')
   
    irvec_max = maxval(irvec,DIM=2)+1

    if (bands_plot_dim .ne. 3) then
       ! Find one_dim_vec which is parallel to one_dim_dir
       ! two_dim_vec - the other two lattice vectors
       ! Along the confined directions, take only irvec=0
       j = 0
       do i=1,3
          if ( abs(abs(real_lattice(one_dim_dir,i)) &
               - sqrt(dot_product(real_lattice(:,i),real_lattice(:,i)))) .lt. eps8 ) then
             one_dim_vec = i
             j = j +1
          end if
       end do
       if ( j .ne. 1 ) call io_error('Error: 1-d lattice vector not defined in plot_cut_hr')
       j=0
       do i=1,3
          if ( i .ne. one_dim_vec ) then
             j = j +1
             two_dim_vec(j)=i
          end if
       end do
       if (bands_plot_dim .eq. 1) then
          irvec_max(two_dim_vec(1))=0
          irvec_max(two_dim_vec(2))=0
       end if
       if (bands_plot_dim .eq. 2) irvec_max(one_dim_vec)=0
    end if  

    nrpts_cut = (2*irvec_max(1)+1)*(2*irvec_max(2)+1)*(2*irvec_max(3)+1)
    allocate(ham_r_cut(num_wann,num_wann,nrpts_cut),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ham_r_cut in plot_cut_hr')
    allocate(irvec_cut(3,nrpts_cut),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating irvec_cut in plot_cut_hr')
    allocate(shift_vec(3,nrpts_cut),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating shift_vec in plot_cut_hr')

    nrpts_tmp = 0
    do n1 = -irvec_max(1), irvec_max(1)    
       do n2 = -irvec_max(2), irvec_max(2)    
loop_n3:  do n3 = -irvec_max(3), irvec_max(3)
             do irpt = 1, nrpts
                i1 = mod(n1 - irvec(1,irpt),mp_grid(1))
                i2 = mod(n2 - irvec(2,irpt),mp_grid(2))
                i3 = mod(n3 - irvec(3,irpt),mp_grid(3))
                if (i1.eq.0 .and. i2.eq.0 .and. i3.eq.0 ) then
                   nrpts_tmp = nrpts_tmp+1 
                   ham_r_cut(:,:,nrpts_tmp) = ham_r(:,:,irpt)
                   irvec_cut(1,nrpts_tmp)=n1
                   irvec_cut(2,nrpts_tmp)=n2
                   irvec_cut(3,nrpts_tmp)=n3
                   cycle loop_n3
                end if
             end do
          end do loop_n3
       end do
    end do

    if ( nrpts_tmp .ne. nrpts_cut) then
       write(stdout,'(a)') 'FAILED TO EXPAND ham_r'
       call io_error('Error in plot_cut_hr')
    end if

    ! AAM: 29/10/2009 Bug fix thanks to Dr Shujun Hu, NIMS, Japan. 
    do irpt = 1, nrpts_cut
       ! line below is incorrect for non-orthorhombic cells
       !shift_vec(:,irpt) = matmul(real_lattice(:,:),real(irvec_cut(:,irpt),dp))
       ! line below is the same as calculating 
       ! matmul(transpose(real_lattice(:,:)),irvec_cut(:,irpt))
       shift_vec(:,irpt) = matmul(real(irvec_cut(:,irpt),dp),real_lattice(:,:))
    end do

    ! note: dist_cutoff_mode does not necessarily follow bands_plot_dim
    ! e.g. for 1-d system (bands_plot_dim=1) we can still apply 3-d dist_cutoff (dist_cutoff_mode=three_dim)
    if (index(dist_cutoff_mode, 'one_dim')>0) then
       do i=1,num_wann
           do j=1,num_wann
              dist_ij_vec(one_dim_dir)= &
                   wannier_centres_translated(one_dim_dir,i) - wannier_centres_translated(one_dim_dir,j)
              do irpt =1, nrpts_cut
                 dist_vec(one_dim_dir) = dist_ij_vec(one_dim_dir)+ shift_vec(one_dim_dir,irpt)
                 dist = abs(dist_vec(one_dim_dir))
                 if ( dist .gt. dist_cutoff ) &
                    ham_r_cut(j,i,irpt) = cmplx_0
              end do
           end do
       end do
    else if (index(dist_cutoff_mode, 'two_dim')>0) then
       do i=1,num_wann
           do j=1,num_wann
              dist_ij_vec(:)=wannier_centres_translated(:,i) - wannier_centres_translated(:,j)
              do irpt =1, nrpts_cut
                 dist_vec(:) = dist_ij_vec(:)+ shift_vec(:,irpt)
                 dist_vec(one_dim_dir) = 0.0_dp
                 dist = sqrt(dot_product(dist_vec,dist_vec))
                 if ( dist .gt. dist_cutoff ) &
                    ham_r_cut(j,i,irpt) = cmplx_0
              end do
           end do
       end do
    else
       do i=1,num_wann
           do j=1,num_wann
              dist_ij_vec(:)=wannier_centres_translated(:,i) - wannier_centres_translated(:,j)
              do irpt =1, nrpts_cut
                 dist_vec(:) = dist_ij_vec(:)+ shift_vec(:,irpt)
                 dist = sqrt(dot_product(dist_vec,dist_vec))
                 if ( dist .gt. dist_cutoff ) &
                    ham_r_cut(j,i,irpt) = cmplx_0
              end do
           end do
       end do
    end if

    do irpt = 1,nrpts_cut
       do i=1, num_wann
          do j=1,num_wann
             if (abs(ham_r_cut(j,i,irpt)) .lt. hr_cutoff) &
                ham_r_cut(j,i,irpt) = cmplx_0 
          end do
       end do
    end do

    write(stdout,'(/1x,a78)') repeat('-',78)
    write(stdout,'(1x,4x,a)') &
                 'Maximum absolute value of Real-space Hamiltonian at each lattice point'
    write(stdout,'(1x,8x,a62)') repeat('-',62)
    write(stdout,'(1x,11x,a,11x,a)') 'Lattice point R', 'Max |H_ij(R)|'
    !  output maximum ham_r_cut at each lattice point
    do irpt = 1, nrpts_cut
       ham_r_tmp(:,:)=abs(ham_r_cut(:,:,irpt))
       write(stdout,'(1x,9x,3I5,9x,F12.6)') irvec_cut(:,irpt),maxval(ham_r_tmp)
    end do 
    !
    return

  end subroutine plot_cut_hr           

  !============================================!
  subroutine plot_interpolate_gnuplot
    !============================================!
    !                                            !
    !     Plots the interpolated band structure  !
    !           in gnuplot format                !
    !============================================!

    use w90_constants,  only : dp
    use w90_io,         only : io_file_unit,seedname
    use w90_parameters, only : num_wann,bands_num_spec_points, &
                               bands_label,num_bands_project

    implicit none

    !
    bndunit=io_file_unit()
    open(bndunit,file=trim(seedname)//'_band.dat',form='formatted')
    gnuunit=io_file_unit()
    open(gnuunit,file=trim(seedname)//'_band.gnu',form='formatted')
    !
    ! Gnuplot format
    !
    do i=1,num_wann
       do nkp=1,total_pts
          if(num_bands_project>0) then
             write(bndunit,'(3E16.8)') xval(nkp),eig_int(i,nkp),bands_proj(i,nkp)
          else
             write(bndunit,'(2E16.8)') xval(nkp),eig_int(i,nkp)
          end if
       enddo
       write(bndunit,*) ' '
    enddo
    close(bndunit)
    ! Axis labels
    glabel(1)=' '//bands_label(1)//' '
    do i=2,num_paths
       if(bands_label(2*(i-1))/=bands_label(2*(i-1)+1)) then
          glabel(i)=bands_label(2*(i-1))//'/'//bands_label(2*(i-1)+1)
       else
          glabel(i)=' '//bands_label(2*(i-1))//' '
       end if
    end do
    glabel(num_spts)=' '//bands_label(bands_num_spec_points)//' '
    ! gnu file
    write(gnuunit,701) xval(total_pts),emin,emax
    do i = 1, num_paths-1
       write(gnuunit,705) sum(kpath_len(1:i)),emin,sum(kpath_len(1:i)),emax
    enddo
    write(gnuunit,702, advance="no") glabel(1),0.0_dp,(glabel(i+1),sum(kpath_len(1:i)),i=1,bands_num_spec_points/2-1)
    write(gnuunit,703) glabel(1+bands_num_spec_points/2),sum(kpath_len(:))
    write(gnuunit,*) 'plot ','"'//trim(seedname)//'_band.dat','"' 
    close(gnuunit)

    if(num_bands_project>0) then
     gnuunit=io_file_unit()
     open(gnuunit,file=trim(seedname)//'_band_proj.gnu',form='formatted')
     write(gnuunit,'(a)') '#File to plot a colour-mapped Bandstructure' 
     write(gnuunit,'(a)') 'set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )'
     write(gnuunit,'(a)') 'unset ztics'
     write(gnuunit,'(a)') 'unset key'
     write(gnuunit,'(a)') '# can make pointsize smaller (~0.5). Too small and nothing is printed'
     write(gnuunit,'(a)') 'set pointsize 0.8'
     write(gnuunit,'(a)') 'set grid xtics'
     write(gnuunit,'(a)') 'set view 0,0'
     write(gnuunit,'(a,f9.5,a)') 'set xrange [0:',xval(total_pts),']'
     write(gnuunit,'(a,f9.5,a,f9.5,a)') 'set yrange [',emin,':',emax,']'
     write(gnuunit,702, advance="no") glabel(1),0.0_dp,(glabel(i+1),sum(kpath_len(1:i)),i=1,bands_num_spec_points/2-1)
     write(gnuunit,703) glabel(1+bands_num_spec_points/2),sum(kpath_len(:))

     write(gnuunit,'(a,a,a,a)') 'splot ','"'//trim(seedname)//'_band.dat','"',' u 1:2:3 w p pt 13 palette'
     write(gnuunit,'(a)') '#use the next lines to make a nice figure for a paper'
     write(gnuunit,'(a)') '#set term postscript enhanced eps color lw 0.5 dl 0.5'
     write(gnuunit,'(a)') '#set pointsize 0.275'
    end if
    !
701 format('set style data dots',/,'set nokey',/, 'set xrange [0:',F8.5,']',/,'set yrange [',F9.5,' :',F9.5,']')
702 format('set xtics (',:20('"',A3,'" ',F8.5,','))
703 format(A3,'" ',F8.5,')')
705 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead')

  end subroutine plot_interpolate_gnuplot

  subroutine plot_interpolate_xmgrace
    !============================================!
    !                                            !
    !     Plots the interpolated band structure  !
    !         in Xmgrace format                  !
    !============================================!

    use w90_io,         only : io_file_unit,seedname,io_date
    use w90_parameters, only : num_wann,bands_num_spec_points

    implicit none

    character (len=9) :: cdate, ctime

    call io_date(cdate, ctime)

    ! Axis labels

    ! Switch any G to Gamma

    do i=1,bands_num_spec_points
       if(bands_label(i)=='G') then
          ctemp(i)='\xG\0'
       else
          ctemp(i)=bands_label(i)
       end if
    end do


    xlabel(1)=' '//trim(ctemp(1))//' '
    do i=2,num_paths
       if(ctemp(2*(i-1))/=ctemp(2*(i-1)+1)) then
          xlabel(i)=trim(ctemp(2*(i-1)))//'/'//trim(ctemp(2*(i-1)+1))
       else
          xlabel(i)=ctemp(2*(i-1))
       end if
    end do
    xlabel(num_spts)=ctemp(bands_num_spec_points)

    gnuunit=io_file_unit()
    open(gnuunit,file=trim(seedname)//'_band.agr',form='formatted')
    !
    ! Xmgrace format
    !
    write(gnuunit,'(a)') '# Grace project file                      '
    write(gnuunit,'(a)') '# written using Wannier90 www.wannier.org '
    write(gnuunit,'(a)') '@version 50113                            '
    write(gnuunit,'(a)') '@page size 792, 612                       '
    write(gnuunit,'(a)') '@page scroll 5%                           '
    write(gnuunit,'(a)') '@page inout 5%                            '
    write(gnuunit,'(a)') '@link page off                            '
    write(gnuunit,'(a)') '@timestamp def "'//cdate//' at '//ctime//'" ' 
    write(gnuunit,'(a)') '@with g0'                                  
    write(gnuunit,'(a)') '@    world xmin 0.00'
    write(gnuunit,'(a,f10.5)') '@    world xmax ',xval(total_pts)
    write(gnuunit,'(a,f10.5)') '@    world ymin ',emin
    write(gnuunit,'(a,f10.5)') '@    world ymax ',emax
    write(gnuunit,'(a)') '@default linewidth 1.5'
    write(gnuunit,'(a)') '@    xaxis  tick on'
    write(gnuunit,'(a)') '@    xaxis  tick major 1'
    write(gnuunit,'(a)') '@    xaxis  tick major color 1'
    write(gnuunit,'(a)') '@    xaxis  tick major linestyle 3'
    write(gnuunit,'(a)') '@    xaxis  tick major grid on'
    write(gnuunit,'(a)') '@    xaxis  tick spec type both'
    write(gnuunit,'(a,i0)') '@    xaxis  tick spec ',1+bands_num_spec_points/2
    write(gnuunit,'(a)') '@    xaxis  tick major 0, 0'
    do i=1,bands_num_spec_points/2
       write(gnuunit,'(a,i0,a,a)') '@    xaxis  ticklabel ',i-1,',', '"'//trim(adjustl(xlabel(i)))//'"'
       write(gnuunit,'(a,i0,a,f10.5)') '@    xaxis  tick major ',i,' , ',sum(kpath_len(1:i))
    end do
    write(gnuunit,'(a,i0,a)') '@    xaxis  ticklabel ',bands_num_spec_points/2 &
         ,',"'//trim(adjustl(xlabel(1+bands_num_spec_points/2)))//'"'
    write(gnuunit,'(a)') '@    xaxis  ticklabel char size 1.500000'
    write(gnuunit,'(a)') '@    yaxis  tick major 10'
    write(gnuunit,'(a)') '@    yaxis  label "Band Energy (eV)"'
    write(gnuunit,'(a)') '@    yaxis  label char size 1.500000'
    write(gnuunit,'(a)') '@    yaxis  ticklabel char size 1.500000'
    do i=1,num_wann
       write(gnuunit,'(a,i0,a)') '@    s',i-1,' line color 1'
    end do
    do i=1,num_wann
       write(gnuunit,'(a,i0)') '@target G0.S',i-1
       write(gnuunit,'(a)') '@type xy'
       do nkp=1,total_pts
          write(gnuunit,'(2E16.8)') xval(nkp),eig_int(i,nkp)
       end do
       write(gnuunit,'(a,i0)') '&'
    end do


  end subroutine plot_interpolate_xmgrace

end subroutine plot_interpolate_bands


  !===========================================================!
  subroutine plot_fermi_surface
    !===========================================================!
    !                                                           !
    !  Prepares a Xcrysden bxsf file to view the fermi surface  !
    !                                                           !
    !===========================================================!

    use w90_constants,   only : dp,cmplx_0,cmplx_i,twopi
    use w90_io,          only : io_error,stdout,io_file_unit,seedname,&
                                io_date,io_time,io_stopwatch
    use w90_parameters,  only : num_wann,fermi_surface_num_points,timing_level,&
                                recip_lattice,nfermi,fermi_energy_list
    use w90_hamiltonian, only : irvec,nrpts,ndegen,ham_r

    implicit none

    complex(kind=dp) ,allocatable :: ham_pack(:)
    complex(kind=dp)   :: fac
    complex(kind=dp) ,allocatable :: ham_kprm(:,:)
    complex(kind=dp) ,allocatable :: U_int(:,:)
    complex(kind=dp) ,allocatable :: cwork(:)
    real(kind=dp) ,allocatable    :: rwork(:)
    real(kind=dp),allocatable  :: eig_int(:,:)
    real(kind=dp)      :: rdotk,time0
    integer, allocatable :: iwork(:),ifail(:)
    integer              :: loop_x,loop_y,loop_z,INFO,ikp,i,j,ierr
    integer              :: irpt,nfound,npts_plot,loop_kpt,bxsf_unit
    character(len=9)     :: cdate, ctime
    !
    if (timing_level>1) call io_stopwatch('plot: fermi_surface',1)
    time0=io_time()
    write(stdout,*)
    write(stdout,'(1x,a)') 'Calculating Fermi surface'
    write(stdout,*)
    !
    if(nfermi>1) call io_error("Error in plot: nfermi>1. Set the fermi level "&
         //"using the input parameter 'fermi_level'") 
    !
    allocate(ham_pack((num_wann*(num_wann+1))/2),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ham_pack plot_fermi_surface') 
    allocate(ham_kprm(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ham_kprm plot_fermi_surface')
    allocate(U_int(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating U_int in plot_fermi_surface')
    allocate(cwork(2*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cwork in plot_fermi_surface')
    allocate(rwork(7*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating rwork in plot_fermi_surface')
    allocate(iwork(5*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating iwork in plot_fermi_surface')
    allocate(ifail(num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ifail in plot_fermi_surface')
    !
    npts_plot=(fermi_surface_num_points+1)**3
    allocate(eig_int(num_wann,npts_plot),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating eig_int in plot_fermi_surface')
    eig_int=0.0_dp
    U_int=(0.0_dp,0.0_dp)
    !
    ikp=0
    do loop_x=1,fermi_surface_num_points+1
       do loop_y=1,fermi_surface_num_points+1
          do loop_z=1,fermi_surface_num_points+1
             ikp=ikp+1

             ham_kprm=cmplx_0
             do irpt=1,nrpts
                rdotk=twopi*real( (loop_x-1)*irvec(1,irpt)+ &
                     (loop_y-1)*irvec(2,irpt) + (loop_z-1)* &
                     irvec(3,irpt) ,dp)/real(fermi_surface_num_points,dp)
                fac=exp(cmplx_i*rdotk)/real(ndegen(irpt),dp)
                ham_kprm=ham_kprm+fac*ham_r(:,:,irpt)
             end do
             ! Diagonalise H_k (->basis of eigenstates)
             do j=1,num_wann
                do i=1,j
                   ham_pack(i+((j-1)*j)/2)=ham_kprm(i,j)
                enddo
             enddo
             call ZHPEVX('N','A','U',num_wann,ham_pack,0.0_dp,0.0_dp,0,0,-1.0_dp, &
                  nfound,eig_int(1,ikp),U_int,num_wann,cwork,rwork,iwork,ifail,info)
             if(info < 0) then
                write(stdout,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
                call io_error('Error in plot_fermi_surface')
             endif
             if(info > 0) then
                write(stdout,'(i3,a)') info,' EIGENVECTORS FAILED TO CONVERGE'
                call io_error('Error in plot_fermi_surface')
             endif
          end do
       end do
    end do

    call io_date(cdate, ctime)
    bxsf_unit=io_file_unit()
    open(bxsf_unit,FILE=trim(seedname)//'.bxsf',STATUS='UNKNOWN',FORM='FORMATTED')
    write(bxsf_unit,*) ' BEGIN_INFO'
    write(bxsf_unit,*) '      #'
    write(bxsf_unit,*) '      # this is a Band-XCRYSDEN-Structure-File'
    write(bxsf_unit,*) '      # for Fermi Surface Visualisation'
    write(bxsf_unit,*) '      #'
    write(bxsf_unit,*) '      # Generated by the Wannier90 code http://www.wannier.org'
    write(bxsf_unit,*) '      # On ',cdate,' at ',ctime
    write(bxsf_unit,*) '      #'
    write(bxsf_unit,*) '      Fermi Energy:', fermi_energy_list(1)
    write(bxsf_unit,*) ' END_INFO'
    write(bxsf_unit,*) 
    write(bxsf_unit,*) ' BEGIN_BLOCK_BANDGRID_3D'
    write(bxsf_unit,*) 'from_wannier_code'
    write(bxsf_unit,*) ' BEGIN_BANDGRID_3D_fermi'
    write(bxsf_unit,*) num_wann
    write(bxsf_unit,*) fermi_surface_num_points+1,fermi_surface_num_points+1,fermi_surface_num_points+1
    write(bxsf_unit,*) '0.0 0.0 0.0'
    write(bxsf_unit,*) (recip_lattice(1,i), i=1,3)
    write(bxsf_unit,*) (recip_lattice(2,i), i=1,3)
    write(bxsf_unit,*) (recip_lattice(3,i), i=1,3)
    do i=1,num_wann
       write(bxsf_unit,*) 'BAND: ',i
       do loop_kpt=1,npts_plot
          write(bxsf_unit,'(2E16.8)') eig_int(i,loop_kpt)
       enddo
    enddo
    write(bxsf_unit,*) 'END_BANDGRID_3D'
    write(bxsf_unit,*) ' END_BLOCK_BANDGRID_3D'
    close(bxsf_unit)

    write(stdout,'(1x,a,f11.3,a)') 'Time to calculate interpolated Fermi surface ',io_time()-time0,' (sec)'
    write(stdout,*)
    !    
    if (timing_level>1) call io_stopwatch('plot: fermi_surface',2)
    !
    return

  end subroutine plot_fermi_surface


  !============================================!
  subroutine plot_wannier
    !============================================!
    !                                            !
    ! Plot the WF in Xcrysden format             !
    !  based on code written by Michel Posternak !
    !                                            !
    !============================================!

    use w90_constants,  only : dp,cmplx_0,cmplx_i,twopi,cmplx_1
    use w90_io,         only : io_error,stdout,io_file_unit,seedname, &
                               io_date,io_stopwatch
    use w90_parameters, only : num_wann,num_bands,num_kpts,u_matrix,spin, &
         ngs=>wannier_plot_supercell,kpt_latt,num_species,atoms_species_num, &
         atoms_symbol,atoms_pos_cart,num_atoms,real_lattice,have_disentangled, &
         ndimwin,lwindow,u_matrix_opt,num_wannier_plot,wannier_plot_list, &
         wannier_plot_mode,wvfn_formatted,timing_level,wannier_plot_format

    implicit none

    real(kind=dp) :: scalfac,tmax,tmaxx,x_0ang,y_0ang,z_0ang
    real(kind=dp) :: fxcry(3),dirl(3,3),w_real,w_imag,ratmax,ratio
    complex(kind=dp), allocatable :: wann_func(:,:,:,:)
    complex(kind=dp), allocatable :: r_wvfn(:,:)
    complex(kind=dp), allocatable :: r_wvfn_tmp(:,:)
    complex(kind=dp) :: catmp,wmod

    logical :: have_file
    integer :: i,j,nsp,nat,nbnd,counter,ierr
    integer :: loop_kpt,ik,ix,iy,iz,nk,ngx,ngy,ngz,nxx,nyy,nzz
    integer :: loop_b,nx,ny,nz,npoint,file_unit,loop_w,num_inc
    character(len=11) :: wfnname
    character(len=60) :: wanxsf,wancube
    character(len=9)  :: cdate, ctime
    logical           :: inc_band(num_bands)  
    !
    if (timing_level>1) call io_stopwatch('plot: wannier',1)
    !
    write(wfnname,200 ) 1,spin
    inquire(file=wfnname,exist=have_file)
    if(.not.have_file) call io_error('plot_wannier: file '//wfnname//' not found') 

    file_unit=io_file_unit()
    if(wvfn_formatted) then
       open(unit=file_unit,file=wfnname,form='formatted')
       read(file_unit,*) ngx,ngy,ngz,nk,nbnd
    else
       open(unit=file_unit,file=wfnname,form='unformatted')
       read(file_unit) ngx,ngy,ngz,nk,nbnd
    end if
    close(file_unit)

200 format ('UNK',i5.5,'.',i1)

    allocate(wann_func(-((ngs(1))/2)*ngx:((ngs(1)+1)/2)*ngx-1,&
         -((ngs(2))/2)*ngy:((ngs(2)+1)/2)*ngy-1,&
         -((ngs(3))/2)*ngz:((ngs(3)+1)/2)*ngz-1,num_wannier_plot),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating wann_func in plot_wannier')
    if(have_disentangled) then
       allocate(r_wvfn_tmp(ngx*ngy*ngz,maxval(ndimwin)),stat=ierr )
       if (ierr/=0) call io_error('Error in allocating r_wvfn_tmp in plot_wannier')
    end if
    allocate(r_wvfn(ngx*ngy*ngz,num_wann),stat=ierr )
    if (ierr/=0) call io_error('Error in allocating r_wvfn in plot_wannier')
    wann_func=cmplx_0


    call io_date(cdate, ctime)
    do loop_kpt=1,num_kpts

       inc_band=.true.
       num_inc=num_wann
       if(have_disentangled) then
          inc_band(:)=lwindow(:,loop_kpt)
          num_inc=ndimwin(loop_kpt)
       end if

       write(wfnname,200 ) loop_kpt,spin
       file_unit=io_file_unit()
       if(wvfn_formatted) then
          open(unit=file_unit,file=wfnname,form='formatted')
          read(file_unit,*) ix,iy,iz,ik,nbnd
       else
          open(unit=file_unit,file=wfnname,form='unformatted')
          read(file_unit) ix,iy,iz,ik,nbnd
       end if

       if ( (ix/=ngx) .or. (iy/=ngy) .or. (iz/=ngz) .or. (ik/=loop_kpt) ) then
          write(stdout,'(1x,a,a)') 'WARNING: mismatch in file',trim(wfnname)
          write(stdout,'(1x,5(a6,I5))')  '   ix=',ix ,'   iy=',iy ,'   iz=',iz ,'   ik=',ik      ,' nbnd=',nbnd
          write(stdout,'(1x,5(a6,I5))')  '  ngx=',ngx,'  ngy=',ngy,'  ngz=',ngz,'  kpt=',loop_kpt,'bands=',num_bands
          call io_error('plot_wannier')
       end if

       if(have_disentangled) then
          counter=1
          do loop_b=1,num_bands
             if(counter>num_inc) exit
             if(wvfn_formatted) then
                do nx=1,ngx*ngy*ngz
                   read(file_unit,*) w_real, w_imag
                   r_wvfn_tmp(nx,counter) = cmplx(w_real,w_imag,kind=dp)
                end do
             else
                read(file_unit) (r_wvfn_tmp(nx,counter),nx=1,ngx*ngy*ngz)
             end if
             if(inc_band(loop_b)) counter=counter+1
          end do
       else
          do loop_b=1,num_bands 
             if(wvfn_formatted) then
                do nx=1,ngx*ngy*ngz
                   read(file_unit,*) w_real, w_imag
                   r_wvfn(nx,loop_b) = cmplx(w_real,w_imag,kind=dp)
                end do
             else
                read(file_unit) (r_wvfn(nx,loop_b),nx=1,ngx*ngy*ngz)
             end if
          end do
       end if

       close(file_unit)

       if(have_disentangled) then
          r_wvfn=cmplx_0
          do loop_w=1,num_wann
             do loop_b=1,num_inc
                r_wvfn(:,loop_w)=r_wvfn(:,loop_w)+ &
                     u_matrix_opt(loop_b,loop_w,loop_kpt)*r_wvfn_tmp(:,loop_b)
             end do
          end do
       end if


       ! nxx, nyy, nzz span a parallelogram in the real space mesh, of side
       ! 2*nphir, and centered around the maximum of phi_i, nphimx(i, 1 2 3)
       !
       ! nx ny nz are the nxx nyy nzz brought back to the unit cell in
       ! which u_nk(r)=cptwrb(r,n)  is represented
       !
       ! There is a big performance improvement in looping over num_wann
       ! in the inner loop. This is poor memory access for wann_func and
       ! but the reduced number of operations wins out. 

       do nzz =-((ngs(3))/2)*ngz,((ngs(3)+1)/2)*ngz-1
          nz=mod(nzz,ngz)
          if(nz.lt.1) nz=nz+ngz
          do nyy=-((ngs(2))/2)*ngy,((ngs(2)+1)/2)*ngy-1
             ny=mod(nyy,ngy)
             if(ny.lt.1) ny=ny+ngy
             do nxx=-((ngs(1))/2)*ngx,((ngs(1)+1)/2)*ngx-1
                nx=mod(nxx,ngx)
                if(nx.lt.1) nx=nx+ngx

                scalfac=kpt_latt(1,loop_kpt)*real(nxx-1,dp)/real(ngx,dp)+ &
                     kpt_latt(2,loop_kpt)*real(nyy-1,dp)/real(ngy,dp)+ &
                     kpt_latt(3,loop_kpt)*real(nzz-1,dp)/real(ngz,dp)
                npoint=nx+(ny-1)*ngx+(nz-1)*ngy*ngx
                catmp=exp(twopi*cmplx_i*scalfac)
                do loop_b=1,num_wann
                   do loop_w=1,num_wannier_plot            
                      wann_func(nxx,nyy,nzz,loop_w)=wann_func(nxx,nyy,nzz,loop_w)+ &
                           u_matrix(loop_b,wannier_plot_list(loop_w),loop_kpt)*r_wvfn(npoint,loop_b)*catmp
                   end do
                end do
             end do
          end do

       end do

    end do !loop over kpoints

    ! fix the global phase by setting the wannier to
    ! be real at the point where it has max. modulus

    do loop_w=1,num_wannier_plot
       tmaxx=0.0
       wmod=cmplx_1
       do nzz=-((ngs(3))/2)*ngz,((ngs(3)+1)/2)*ngz-1
          do nyy=-((ngs(2))/2)*ngy,((ngs(2)+1)/2)*ngy-1
             do nxx=-((ngs(1))/2)*ngx,((ngs(1)+1)/2)*ngx-1
                wann_func(nxx,nyy,nzz,loop_w)= wann_func(nxx,nyy,nzz,loop_w)/ real(num_kpts,dp)
                tmax=real(wann_func(nxx,nyy,nzz,loop_w)* & 
                     conjg(wann_func(nxx,nyy,nzz,loop_w)),dp)
                if (tmax>tmaxx) then
                   tmaxx=tmax
                   wmod=wann_func(nxx,nyy,nzz,loop_w)
                end if
             end do
          end do
       end do
       wmod=wmod/sqrt(real(wmod)**2+aimag(wmod)**2)
       wann_func(:,:,:,loop_w)=wann_func(:,:,:,loop_w)/wmod
    end do
    !
    ! Check the 'reality' of the WF
    !
    do loop_w=1,num_wannier_plot
       ratmax=0.0_dp
       do nzz=-((ngs(3))/2)*ngz,((ngs(3)+1)/2)*ngz-1
          do nyy=-((ngs(2))/2)*ngy,((ngs(2)+1)/2)*ngy-1
             do nxx=-((ngs(1))/2)*ngx,((ngs(1)+1)/2)*ngx-1
                if (abs(real(wann_func(nxx,nyy,nzz,loop_w),dp))>=0.01_dp) then
                   ratio=abs(aimag(wann_func(nxx,nyy,nzz,loop_w)))/ &
                        abs(real(wann_func(nxx,nyy,nzz,loop_w),dp))
                   ratmax=max(ratmax,ratio)
                end if
             end do
          end do
       end do
       write(stdout,'(6x,a,i4,7x,a,f11.6)') 'Wannier Function Num: ',wannier_plot_list(loop_w),&
            'Maximum Im/Re Ratio = ',ratmax
    end do
    write(stdout,*) ' '

    if (wannier_plot_format.eq.'xcrysden') then
       call internal_xsf_format()
    elseif (wannier_plot_format.eq.'cube') then
       call internal_cube_format()
    else
       call io_error('wannier_plot_format not recognised in wannier_plot')
    endif

    if (timing_level>1) call io_stopwatch('plot: wannier',2)

    return

  contains


    !============================================!
    subroutine internal_cube_format()
    !============================================!
    !                                            !
    ! Write WFs in Gaussian cube format.         !
    !                                            !
    !============================================!


      use w90_constants,  only: bohr
      use w90_parameters, only: recip_lattice,iprint,&
           wannier_plot_radius,wannier_centres,atoms_symbol, &
           translate_home_cell
      use w90_utility,    only: utility_translate_home

      implicit none

      real(kind=dp), allocatable :: wann_cube(:,:,:)
      real(kind=dp) :: rstart(3),rend(3),rlength(3),orig(3),dgrid(3)
      real(kind=dp) :: moda(3),modb(3)
      real(kind=dp) :: radius,val_Q
      real(kind=dp) :: wc(3,num_wann) 
      integer :: ierr,iname,max_elements,iw
      integer :: isp,iat,nzz,nyy,nxx,loop_w,qxx,qyy,qzz,wann_index
      integer :: istart(3),iend(3),ilength(3)
      integer :: ixx,iyy,izz
      integer, allocatable :: atomic_Z(:)
      character(len=2), dimension(109) :: periodic_table= (/ &
           & 'H ',                                                                                'He', &
           & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
           & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
           & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
           & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
           & 'Cs','Ba', &
           & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
           & 'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
           & 'Fr','Ra', &
           & 'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr', &
           & 'Rf','Db','Sg','Bh','Hs','Mt' /)      

      allocate(atomic_Z(num_species),stat=ierr)
      if (ierr.ne.0) call io_error('Error: allocating atomic_Z in wannier_plot')

      radius = wannier_plot_radius
      val_Q = 1.0_dp ! dummy value for cube file

      ! Assign atomic numbers to species
      max_elements=size(periodic_table)
      do isp=1,num_species
        do iname=1,max_elements
           if ( atoms_symbol(isp).eq.periodic_table(iname) ) then
              atomic_Z(isp)=iname
              exit
           endif
        enddo
     end do
        
202   format (a,'_',i5.5,'.cube')

      ! Lengths of real and reciprocal lattice vectors 
      do i=1,3
         moda(i) = sqrt( real_lattice(i,1)*real_lattice(i,1) &
              + real_lattice(i,2)*real_lattice(i,2) &
              + real_lattice(i,3)*real_lattice(i,3) )
         modb(i) = sqrt( recip_lattice(i,1)*recip_lattice(i,1) &
              + recip_lattice(i,2)*recip_lattice(i,2) &
              + recip_lattice(i,3)*recip_lattice(i,3) )
      enddo
      
      ! Grid spacing in each lattice direction
      dgrid(1) = moda(1)/ngx; dgrid(2) = moda(2)/ngy; dgrid(3)=moda(3)/ngz

      ! Translate centres to home unit cell
      wc = wannier_centres
      if (translate_home_cell) then
         do iw=1,num_wann
            call utility_translate_home(wc(:,iw),real_lattice,recip_lattice)
         enddo
      endif

      ! Loop over WFs
      do loop_w=1,num_wannier_plot
         
         wann_index = wannier_plot_list(loop_w)
         write(wancube,202) trim(seedname),wann_index

         ! Find start and end of cube wrt simulation cell origin
         do i=1,3
            ! ... in terms of distance along each lattice vector direction i
            rstart(i) = ( wc(1,wann_index)*recip_lattice(i,1) &
                 + wc(2,wann_index)*recip_lattice(i,2) &
                 + wc(3,wann_index)*recip_lattice(i,3) - radius*modb(i) ) * moda(i) / twopi
            rend(i) = ( wc(1,wann_index)*recip_lattice(i,1) &
                 + wc(2,wann_index)*recip_lattice(i,2) &
                 + wc(3,wann_index)*recip_lattice(i,3) + radius*modb(i) ) * moda(i) / twopi
         enddo

         rlength(:) = rend(:) - rstart(:)
         ilength(:) = ceiling(rlength(:)/dgrid(:))
         
         ! ... in terms of integer gridpoints along each lattice vector direction i
         istart(:)  = floor(rstart(:)/dgrid(:)) + 1
         iend(:)    = istart(:) + ilength(:) - 1 

         ! Origin of cube wrt simulation cell in Cartesian co-ordinates
         do i=1,3
            orig(i) = real(istart(1)-1,dp)*dgrid(1)*real_lattice(1,i)/moda(1) &
                 + real(istart(2)-1,dp)*dgrid(2)*real_lattice(2,i)/moda(2) &
                 + real(istart(3)-1,dp)*dgrid(3)*real_lattice(3,i)/moda(3)
         enddo

         ! Debugging
         if (iprint>3) then
            write(stdout,'(a,3i12)')     'ngi     =', ngx,ngy,ngz
            write(stdout,'(a,3f12.6)')   'dgrid   =',(dgrid(i),i=1,3)
            write(stdout,'(a,3f12.6)')   'rstart  =',(rstart(i),i=1,3)
            write(stdout,'(a,3f12.6)')   'rend    =',(rend(i),i=1,3)
            write(stdout,'(a,3f12.6)')   'rlength =',(rlength(i),i=1,3)
            write(stdout,'(a,3i12)')     'istart  =',(istart(i),i=1,3)
            write(stdout,'(a,3i12)')     'iend    =',(iend(i),i=1,3)
            write(stdout,'(a,3i12)')     'ilength =',(ilength(i),i=1,3)
            write(stdout,'(a,3f12.6)')   'orig    =',(orig(i),i=1,3)
            write(stdout,'(a,3f12.6,/)') 'wann_cen=',(wannier_centres(i,wann_index),i=1,3)
         endif

         allocate(wann_cube(1:ilength(1),1:ilength(2),1:ilength(3)),stat=ierr)
         if (ierr.ne.0) call io_error('Error: allocating wann_cube in wannier_plot')

         ! initialise
         wann_cube = 0.0_dp

         do nzz=1,ilength(3)
            qzz=nzz+istart(3)-1
            izz=int((abs(qzz)-1)/ngz)
!            if (qzz.lt.-ngz) qzz=qzz+izz*ngz
!            if (qzz.gt.(ngs(3)-1)*ngz-1) then
            if (qzz.lt.(-((ngs(3))/2)*ngz)) qzz=qzz+izz*ngz
            if (qzz.gt.((ngs(3)+1)/2)*ngz-1) then
               write(stdout,*) 'Error plotting WF cube. Try one of the following:'
               write(stdout,*) '   (1) increase wannier_plot_supercell;'
               write(stdout,*) '   (2) decrease wannier_plot_radius;'
               write(stdout,*) '   (3) set wannier_plot_format=xcrysden'
               call io_error('Error plotting WF cube.')
            endif
            do nyy=1,ilength(2)
               qyy=nyy+istart(2)-1
               iyy=int((abs(qyy)-1)/ngy)
!               if (qyy.lt.-ngy) qyy=qyy+iyy*ngy
!               if (qyy.gt.(ngs(2)-1)*ngy-1) then
               if (qyy.lt.(-((ngs(2))/2)*ngy)) qyy=qyy+iyy*ngy
               if (qyy.gt.((ngs(2)+1)/2)*ngy-1) then
                  write(stdout,*) 'Error plotting WF cube. Try one of the following:'
                  write(stdout,*) '   (1) increase wannier_plot_supercell;'
                  write(stdout,*) '   (2) decrease wannier_plot_radius;'
                  write(stdout,*) '   (3) set wannier_plot_format=xcrysden'
                  call io_error('Error plotting WF cube.')
               endif
               do nxx=1,ilength(1)
                  qxx=nxx+istart(1)-1
                  ixx=int((abs(qxx)-1)/ngx)
!                  if (qxx.lt.-ngx) qxx=qxx+ixx*ngx
!                  if (qxx.gt.(ngs(1)-1)*ngx-1) then
                  if (qxx.lt.(-((ngs(1))/2)*ngx)) qxx=qxx+ixx*ngx
                  if (qxx.gt.((ngs(1)+1)/2)*ngx-1) then
                     write(stdout,*) 'Error plotting WF cube. Try one of the following:'
                     write(stdout,*) '   (1) increase wannier_plot_supercell;'
                     write(stdout,*) '   (2) decrease wannier_plot_radius;'
                     write(stdout,*) '   (3) set wannier_plot_format=xcrysden'
                     call io_error('Error plotting WF cube.')                     
                  endif
                  wann_cube(nxx,nyy,nzz) = real(wann_func(qxx,qyy,qzz,loop_w),dp)
               enddo
            enddo
         enddo

         ! Write cube file (everything in Bohr)
         file_unit=io_file_unit()
         open(unit=file_unit,file=trim(wancube),form='formatted',status='unknown')
         ! First two lines are comments
         write(file_unit,*) '     Generated by Wannier90 code http://www.wannier.org'
         write(file_unit,*) '     On ',cdate,' at ',ctime
         ! Number of atoms, origin of cube (Cartesians) wrt simulation cell
         write(file_unit,'(i4,3f13.5)') num_atoms, orig(1)/bohr, orig(2)/bohr, orig(3)/bohr
         ! Number of grid points in each direction, lattice vector
         write(file_unit,'(i4,3f13.5)') ilength(1), real_lattice(1,1)/(real(ngx,dp)*bohr), &
              real_lattice(1,2)/(real(ngy,dp)*bohr), real_lattice(1,3)/(real(ngz,dp)*bohr)
         write(file_unit,'(i4,3f13.5)') ilength(2), real_lattice(2,1)/(real(ngx,dp)*bohr), &
              real_lattice(2,2)/(real(ngy,dp)*bohr), real_lattice(2,3)/(real(ngz,dp)*bohr)
         write(file_unit,'(i4,3f13.5)') ilength(3), real_lattice(3,1)/(real(ngx,dp)*bohr), &
              real_lattice(3,2)/(real(ngy,dp)*bohr), real_lattice(3,3)/(real(ngz,dp)*bohr)

         ! Atomic number, valence charge, position of atom
         do isp=1,num_species
            do iat=1,atoms_species_num(isp)
               write(file_unit,'(i4,4f13.5)') atomic_Z(isp), val_Q, (atoms_pos_cart(i,iat,isp)/bohr,i=1,3)
            end do
         end do

         ! Volumetric data in batches of 6 values per line, 'z'-direction first.
         do nxx=1,ilength(1)
            do nyy=1,ilength(2)
               do nzz=1,ilength(3)
                  write(file_unit,'(E13.5)',advance='no') wann_cube(nxx,nyy,nzz)
                  if ((mod(nzz,6).eq.0) .or. (nzz.eq.ilength(3))) write(file_unit,'(a)') ''
               enddo
            enddo
         enddo

         deallocate(wann_cube,stat=ierr)
         if (ierr.ne.0) call io_error('Error: deallocating wann_cube in wannier_plot')
        
      end do

      deallocate(atomic_Z,stat=ierr)
      if (ierr.ne.0) call io_error('Error: deallocating atomic_Z in wannier_plot')

      return

    end subroutine internal_cube_format




    subroutine internal_xsf_format()

      implicit none

201   format (a,'_',i5.5,'.xsf')

      ! this is to create the WF...xsf output, to be read by XCrySDen
      ! (coordinates + isosurfaces)

      x_0ang=-real(((ngs(1))/2)*ngx+1,dp)/real(ngx,dp)*real_lattice(1,1)- &
           real(((ngs(2))/2)*ngy+1,dp)/real(ngy,dp)*real_lattice(2,1)-  &
           real(((ngs(3))/2)*ngz+1,dp)/real(ngz,dp)*real_lattice(3,1)
      y_0ang=-real(((ngs(1))/2)*ngx+1,dp)/real(ngx,dp)*real_lattice(1,2)- &
           real(((ngs(2))/2)*ngy+1,dp)/real(ngy,dp)*real_lattice(2,2)-  &
           real(((ngs(3))/2)*ngz+1,dp)/real(ngz,dp)*real_lattice(3,2)
      z_0ang=-real(((ngs(1))/2)*ngx+1,dp)/real(ngx,dp)*real_lattice(1,3)- &
           real(((ngs(2))/2)*ngy+1,dp)/real(ngy,dp)*real_lattice(2,3)-  &
           real(((ngs(3))/2)*ngz+1,dp)/real(ngz,dp)*real_lattice(3,3)

      fxcry(1)=real(ngs(1)*ngx-1,dp)/real(ngx,dp)
      fxcry(2)=real(ngs(2)*ngy-1,dp)/real(ngy,dp)
      fxcry(3)=real(ngs(3)*ngz-1,dp)/real(ngz,dp)
      do j=1,3
         dirl(:,j)=fxcry(:)*real_lattice(:,j)
      end do

      do loop_b=1,num_wannier_plot

         write(wanxsf,201) trim(seedname),wannier_plot_list(loop_b)

         file_unit=io_file_unit()
         open(unit=file_unit,file=trim(wanxsf),form='formatted',status='unknown')
         write(file_unit,*) '      #'
         write(file_unit,*) '      # Generated by the Wannier90 code http://www.wannier.org'
         write(file_unit,*) '      # On ',cdate,' at ',ctime
         write(file_unit,*) '      #'
         ! should pass this into the code
         if (index(wannier_plot_mode,'mol')>0 ) then
            write(file_unit,'("ATOMS")')
         else
            write(file_unit,'("CRYSTAL")')
            write(file_unit,'("PRIMVEC")')
            write(file_unit,'(3f12.7)') real_lattice(1,1),real_lattice(1,2),real_lattice(1,3)
            write(file_unit,'(3f12.7)') real_lattice(2,1),real_lattice(2,2),real_lattice(2,3)
            write(file_unit,'(3f12.7)') real_lattice(3,1),real_lattice(3,2),real_lattice(3,3)
            write(file_unit,'("CONVVEC")')
            write(file_unit,'(3f12.7)') real_lattice(1,1),real_lattice(1,2),real_lattice(1,3)
            write(file_unit,'(3f12.7)') real_lattice(2,1),real_lattice(2,2),real_lattice(2,3)
            write(file_unit,'(3f12.7)') real_lattice(3,1),real_lattice(3,2),real_lattice(3,3)
            write(file_unit,'("PRIMCOORD")')
            write(file_unit,'(i6,"  1")') num_atoms
         endif
         do nsp=1,num_species
            do nat=1,atoms_species_num(nsp)
               write(file_unit,'(a2,3x,3f12.7)') atoms_symbol(nsp),(atoms_pos_cart(i,nat,nsp),i=1,3)
            end do
         end do

         write(file_unit,'(/)')
         write(file_unit,'("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
         write(file_unit,'(3i6)') ngs(1)*ngx,ngs(2)*ngy,ngs(3)*ngz
         write(file_unit,'(3f12.6)') x_0ang,y_0ang,z_0ang
         write(file_unit,'(3f12.7)') dirl(1,1),dirl(1,2),dirl(1,3)
         write(file_unit,'(3f12.7)') dirl(2,1),dirl(2,2),dirl(2,3)
         write(file_unit,'(3f12.7)') dirl(3,1),dirl(3,2),dirl(3,3)
         write(file_unit,'(6e13.5)') &
              (((real(wann_func(nx,ny,nz,loop_b)),nx=-((ngs(1))/2)*ngx,((ngs(1)+1)/2)*ngx-1), &
              ny=-((ngs(2))/2)*ngy,((ngs(2)+1)/2)*ngy-1),nz=-((ngs(3))/2)*ngz,((ngs(3)+1)/2)*ngz-1)
         write(file_unit,'("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
         close(file_unit)

      end do

      return

    end subroutine internal_xsf_format

  end subroutine plot_wannier

end module w90_plot
 
