!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2004,2006 Jonathan Yates, Arash Mostofi,     !
!            Nicola Marzari, Ivo Souza, David Vanderbilt     !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

module w90_plot

  use w90_constants, only : dp
  implicit none

  private

  complex(kind=dp),allocatable  :: ham_r(:,:,:)
  integer, allocatable          :: irvec(:,:)
  integer, allocatable          :: ndegen(:)
  integer                       :: nrpts

  public :: plot_main


contains

  !============================================!
  subroutine plot_main( )
    !============================================!

    use w90_constants, only : cmplx_0
    use w90_io, only        : io_error,stdout,io_stopwatch
    use w90_parameters, only    : num_kpts,bands_plot,dos_plot,h_plot, &
         mp_grid,kpt_latt,fermi_surface_plot,num_wann,wannier_plot,timing_level

    implicit none

    integer :: ierr,nkp
    logical :: have_gamma

    if (timing_level>0) call io_stopwatch('plot: main',1)

    write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,'(1x,a)') '|                               PLOTTING                                    |'
    write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,*)

    if(bands_plot .or. dos_plot .or. fermi_surface_plot .or. h_plot) then
       ! Check if the kmesh includes the gamma point
       have_gamma=.false.
       do nkp=1,num_kpts
           if (all(kpt_latt(:,nkp)<0.000001_dp)) have_gamma=.true.       
       end do
       if(.not. have_gamma) &
    write(stdout,'(1x,a)') '!!!! Kpoint grid does not include Gamma. Interpolation may be incorrect. !!!!'
       ! Find the number of points in the Wigner-Seitz cell
       call wigner_seitz(count_pts=.true.)
       allocate(irvec(3,3*num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating irvec in plot_main')
       irvec=0
       allocate(ndegen(3*num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating ndegen in plot_main')
       ndegen=0
       allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating ham_r in plot_main')
       ham_r=cmplx_0
       ! Set up the wigner_seitz vectors and transform Hamiltonian to WF basis
       call wigner_seitz(count_pts=.false.)
       call plot_get_hr
       !
       if(bands_plot) call plot_interpolate_bands
       !
       if(fermi_surface_plot) call plot_fermi_surface
       !
       if(h_plot) call plot_ham_r
       !
       deallocate(ham_r,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating ham_r in plot_main')
       deallocate(ndegen,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating ndegen in plot_main')
       deallocate(irvec,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating irvec in plot_main')
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

    use w90_constants,  only : dp,cmplx_0,cmplx_i,twopi
    use w90_io,         only : io_error,stdout,io_file_unit,seedname,&
                               io_time,io_stopwatch
    use w90_parameters, only : num_wann,bands_num_points,recip_metric,&
                               bands_num_spec_points,timing_level, &
                               bands_spec_points,bands_label

    implicit none

    complex(kind=dp)   :: ham_pack((num_wann*(num_wann+1))/2)
    complex(kind=dp)   :: fac
    complex(kind=dp)   :: ham_kprm(num_wann,num_wann)
    complex(kind=dp)   :: U_int(num_wann,num_wann)
    complex(kind=dp)   :: cwork(2*num_wann)     
    real(kind=dp)      :: rwork(7*num_wann)     
    real(kind=dp)      :: kpath_len(bands_num_spec_points/2)     
    integer            :: kpath_pts(bands_num_spec_points/2)     
    real(kind=dp), allocatable :: xval(:)
    real(kind=dp), allocatable :: eig_int(:,:),plot_kpoint(:,:)
    real(kind=dp)        :: rdotk,vec(3),emin,emax,time0
    integer              :: iwork(5*num_wann)
    integer              :: info,ifail,i,j
    integer              :: loop_rpt,nfound,loop_kpt,counter
    integer              :: loop_spts,total_pts,loop_i,nkp
    integer              :: num_paths,num_spts,ierr
    integer              :: bndunit,gnuunit
    character(len=3),allocatable  :: glabel(:)
    !
    if (timing_level>1) call io_stopwatch('plot: interpolate_bands',1)
    !
    time0=io_time()
    write(stdout,*) 
    write(stdout,'(1x,a)') 'Calculating interpolated band-structure'
    write(stdout,*) 
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
    allocate(glabel(num_spts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating num_spts in plot_interpolate_bands')
    eig_int=0.0_dp
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
       write(bndunit,'(3f12.6,3x,f3.1)') (plot_kpoint(loop_i,loop_spts),loop_i=1,3),1.0_dp
    end do
    close(bndunit)
    !
    ! Interpolate the Hamiltonian at each kpoint
    !
    do loop_kpt=1,total_pts
       ham_kprm=cmplx_0
       do loop_rpt=1,nrpts
          rdotk=twopi*dot_product(plot_kpoint(:,loop_kpt),irvec(:,loop_rpt))
          fac=exp(cmplx_i*rdotk)/real(ndegen(loop_rpt),dp)
          ham_kprm=ham_kprm+fac*ham_r(:,:,loop_rpt)
       end do
       ! Diagonalise H_k (->basis of eigenstates)
       do j=1,num_wann
          do i=1,j
             ham_pack(i+((j-1)*j)/2)=ham_kprm(i,j)
          enddo
       enddo
       call ZHPEVX('N','A','U',num_wann,ham_pack,0.0_dp,0.0_dp,0,0,-1.0_dp, &
            nfound,eig_int(1,loop_kpt),U_int,num_wann,cwork,rwork,iwork,ifail,info)
       if(info < 0) then
          write(stdout,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
          call io_error('Error in plot_interpolate_bands')
       endif
       if(info > 0) then
          write(stdout,'(i3,a)') info,' EIGENVECTORS FAILED TO CONVERGE'
          call io_error('Error in plot_interpolate_bands')
       endif
    end do
    !
    ! Interpolation Finished!
    ! Now we write plotting files
    !
    emin=minval(eig_int)-1.0_dp
    emax=maxval(eig_int)+1.0_dp
    bndunit=io_file_unit()
    open(bndunit,file=trim(seedname)//'_band.dat',form='formatted')
    gnuunit=io_file_unit()
    open(gnuunit,file=trim(seedname)//'_band.gnu',form='formatted')
    !
    ! Gnuplot format
    !
    do i=1,num_wann
       do nkp=1,total_pts
          write(bndunit,'(2E16.8)') xval(nkp),eig_int(i,nkp)
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

    write(stdout,'(1x,a,f11.3,a)')  'Time to calculate interpolated band structure ',io_time()-time0,' (sec)'
    write(stdout,*)
    !
    if (timing_level>1) call io_stopwatch('plot: interpolate_bands',2)
    !
701 format('set data style dots',/,'set nokey',/, 'set xrange [0:',F8.5,']',/,'set yrange [',F9.5,' :',F9.5,']')
702 format('set xtics (',:20('"',A3,'" ',F8.5,','))
703 format(A3,'" ',F8.5,')')
705 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead')

  end subroutine plot_interpolate_bands

  !===========================================================!
  subroutine plot_fermi_surface
    !===========================================================!
    !                                                           !
    !  Prepares a Xcrysden bxsf file to view the fermi surface  !
    !                                                           !
    !===========================================================!

    use w90_constants,  only : dp,cmplx_0,cmplx_i,twopi
    use w90_io,         only : io_error,stdout,io_file_unit,seedname,&
         io_date,io_time,io_stopwatch
    use w90_parameters, only : num_wann,fermi_surface_num_points,timing_level,&
         recip_lattice,fermi_energy

    implicit none

    complex(kind=dp)   :: ham_pack((num_wann*(num_wann+1))/2)
    complex(kind=dp)   :: fac
    complex(kind=dp)   :: ham_kprm(num_wann,num_wann)
    complex(kind=dp)   :: U_int(num_wann,num_wann)
    complex(kind=dp)   :: cwork(2*num_wann)
    real(kind=dp)      :: rwork(7*num_wann)     
    real(kind=dp),allocatable  :: eig_int(:,:)
    real(kind=dp)      :: rdotk,time0
    integer            :: iwork(5*num_wann)     
    integer            :: loop_x,loop_y,loop_z,INFO,ikp,ifail,i,j,ierr
    integer            :: loop_rpt,nfound,npts_plot,loop_kpt,bxsf_unit
    character(len=9)   :: cdate, ctime
    !
    if (timing_level>1) call io_stopwatch('plot: fermi_surface',1)
    !
    time0=io_time()
    write(stdout,*) 
    write(stdout,'(1x,a)') 'Calculating Fermi surface'
    write(stdout,*) 
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
             do loop_rpt=1,nrpts
                rdotk=twopi*real( (loop_x-1)*irvec(1,loop_rpt)+ &
                     (loop_y-1)*irvec(2,loop_rpt) + (loop_z-1)* &
                     irvec(3,loop_rpt) ,dp)/real(fermi_surface_num_points,dp)
                fac=exp(cmplx_i*rdotk)/real(ndegen(loop_rpt),dp)
                ham_kprm=ham_kprm+fac*ham_r(:,:,loop_rpt)
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
    write(bxsf_unit,*) '      Fermi Energy:', fermi_energy 
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

    use w90_constants, only : dp,cmplx_0,cmplx_i,twopi,cmplx_1
    use w90_io, only        : io_error,stdout,io_file_unit,seedname,io_time, &
         io_date,io_stopwatch
    use w90_parameters, only    : num_wann,num_bands,num_kpts,u_matrix,spin, &
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

    allocate(wann_func(-ngx:(ngs-1)*ngx-1,-ngy:(ngs-1)*ngy-1,-ngz:(ngs-1)*ngz-1,num_wannier_plot),stat=ierr ) 
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

       do nzz=-ngz,(ngs-1)*ngz-1
          nz=mod(nzz,ngz)
          if(nz.lt.1) nz=nz+ngz
          do nyy=-ngy,(ngs-1)*ngy-1
             ny=mod(nyy,ngy)
             if(ny.lt.1) ny=ny+ngy
             do nxx=-ngx,(ngs-1)*ngx-1
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
       do nzz=-ngz,(ngs-1)*ngz-1
          do nyy=-ngy,(ngs-1)*ngy-1
             do nxx=-ngx,(ngs-1)*ngx-1
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
       do nzz=-ngz,(ngs-1)*ngz-1
          do nyy=-ngy,(ngs-1)*ngy-1
             do nxx=-ngx,(ngs-1)*ngx-1
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


      use w90_constants,  only: bohr!,periodic_table
      use w90_parameters, only: recip_lattice,iprint,&
           wannier_plot_radius,wannier_centres,atoms_symbol
 
      implicit none

      real(kind=dp), allocatable :: wann_cube(:,:,:)
      real(kind=dp) :: rstart(3),rend(3),rlength(3),orig(3),dgrid(3)
      real(kind=dp) :: moda(3),modb(3)
      real(kind=dp) :: radius,val_Q
      integer :: ierr,iname,max_elements
      integer :: isp,iat,nzz,nyy,nxx,loop_w,qxx,qyy,qzz,wann_index
      integer :: istart(3),iend(3),ilength(3)
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

      ! Loop over WFs
      do loop_w=1,num_wannier_plot
         
         wann_index = wannier_plot_list(loop_w)
         write(wancube,202) trim(seedname),wann_index

         ! Find start and end of cube wrt simulation cell origin
         do i=1,3
            ! ... in terms of distance along each lattice vector direction i
            rstart(i) = ( wannier_centres(1,wann_index)*recip_lattice(i,1) &
                 + wannier_centres(2,wann_index)*recip_lattice(i,2) &
                 + wannier_centres(3,wann_index)*recip_lattice(i,3) - radius*modb(i) ) * moda(i) / twopi
            rend(i) = ( wannier_centres(1,wann_index)*recip_lattice(i,1) &
                 + wannier_centres(2,wann_index)*recip_lattice(i,2) &
                 + wannier_centres(3,wann_index)*recip_lattice(i,3) + radius*modb(i) ) * moda(i) / twopi
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
            write(stdout,'(a,3i12)')    'ngi     =', ngx,ngy,ngz
            write(stdout,'(a,3f12.6)') 'dgrid   =',(dgrid(i),i=1,3)
            write(stdout,'(a,3f12.6)') 'rstart  =',(rstart(i),i=1,3)
            write(stdout,'(a,3f12.6)') 'rend    =',(rend(i),i=1,3)
            write(stdout,'(a,3f12.6)') 'rlength =',(rlength(i),i=1,3)
            write(stdout,'(a,3i12)')    'istart  =',(istart(i),i=1,3)
            write(stdout,'(a,3i12)')    'iend    =',(iend(i),i=1,3)
            write(stdout,'(a,3i12)')    'ilength =',(ilength(i),i=1,3)
            write(stdout,'(a,3f12.6)') 'orig    =',(orig(i),i=1,3)
            write(stdout,'(a,3f12.6,/)') 'wann_cen=',(wannier_centres(i,wann_index),i=1,3)
         endif

         allocate(wann_cube(1:ilength(1),1:ilength(2),1:ilength(3)),stat=ierr)
         if (ierr.ne.0) call io_error('Error: allocating wann_cube in wannier_plot')

         ! Copy WF to cube
         do nzz=1,ilength(3)
            do nyy=1,ilength(2)
               do nxx=1,ilength(1)
                  qxx=nxx+istart(1)-ngx-1
                  if (qxx.lt.-ngx) qxx=qxx+ngx
                  if (qxx.gt.-1)   qxx=qxx-ngx
                  qyy=nyy+istart(2)-ngy-1
                  if (qyy.lt.-ngy) qyy=qyy+ngy
                  if (qyy.gt.-1)   qyy=qyy-ngy
                  qzz=nzz+istart(3)-ngz-1
                  if (qzz.lt.-ngz) qzz=qzz+ngz
                  if (qzz.gt.-1)   qzz=qzz-ngz
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

      x_0ang=-real(ngx+1,dp)/real(ngx,dp)*real_lattice(1,1)- &
           real(ngy+1,dp)/real(ngy,dp)*real_lattice(2,1)-  &
           real(ngz+1,dp)/real(ngz,dp)*real_lattice(3,1)
      y_0ang=-real(ngx+1,dp)/real(ngx,dp)*real_lattice(1,2)- &
           real(ngy+1,dp)/real(ngy,dp)*real_lattice(2,2)-  &
           real(ngz+1,dp)/real(ngz,dp)*real_lattice(3,2)
      z_0ang=-real(ngx+1,dp)/real(ngx,dp)*real_lattice(1,3)- &
           real(ngy+1,dp)/real(ngy,dp)*real_lattice(2,3)-  &
           real(ngz+1,dp)/real(ngz,dp)*real_lattice(3,3)

      fxcry(1)=real(ngs*ngx-1,dp)/real(ngx,dp)
      fxcry(2)=real(ngs*ngy-1,dp)/real(ngy,dp)
      fxcry(3)=real(ngs*ngz-1,dp)/real(ngz,dp)
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
         write(file_unit,'(3i6)') ngs*ngx,ngs*ngy,ngs*ngz
         write(file_unit,'(3f12.6)') x_0ang,y_0ang,z_0ang
         write(file_unit,'(3f12.7)') dirl(1,1),dirl(1,2),dirl(1,3)
         write(file_unit,'(3f12.7)') dirl(2,1),dirl(2,2),dirl(2,3)
         write(file_unit,'(3f12.7)') dirl(3,1),dirl(3,2),dirl(3,3)
         write(file_unit,'(6e13.5)') &
              (((real(wann_func(nx,ny,nz,loop_b)),nx=-ngx,(ngs-1)*ngx-1), &
              ny=-ngy,(ngs-1)*ngy-1),nz=-ngz,(ngs-1)*ngz-1)
         write(file_unit,'("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
         close(file_unit)

      end do

      return

    end subroutine internal_xsf_format

  end subroutine plot_wannier

 


 !============================================!
  subroutine plot_get_hr
  !============================================!
  !                                            !
  !  Calculate the Hamiltonian in the WF basis !
  !                                            !
  !============================================!

    use w90_constants, only : dp,cmplx_0,cmplx_i,twopi
    use w90_io, only        : io_error,io_stopwatch
    use w90_parameters, only    : num_kpts,u_matrix,eigval,num_wann,kpt_latt, &
         u_matrix_opt,num_bands,lwindow,have_disentangled,ndimwin,timing_level

    implicit none
  
    complex(kind=dp)   :: ham_k(num_wann,num_wann,num_kpts)
    complex(kind=dp)   :: fac
    real(kind=dp)      :: rdotk
    real(kind=dp)      :: eigval_opt(num_bands,num_kpts)
    real(kind=dp)      :: eigval2(num_wann,num_kpts)
    integer            :: loop_kpt,i,j,m,loop_rpt,counter

    if (timing_level>1) call io_stopwatch('plot: get_hr',1)

    ham_k=cmplx_0
    eigval_opt=0.0_dp
    eigval2=0.0_dp


    if(have_disentangled) then

       ! slim down eigval to contain states within the outer window

       do loop_kpt=1,num_kpts
          counter=0
          do j=1,num_bands
             if(lwindow(j,loop_kpt)) then
                counter=counter+1
                eigval_opt(counter,loop_kpt)=eigval(j,loop_kpt)
             end if
          end do
       end do
       
       ! rotate eigval into the optimal subspace
       ! in general eigval would be a matrix at each kpoints
       ! but we choose u_matrix_opt such that the Hamiltonian is
       ! diagonal at each kpoint. (I guess we should check it here)
       
       do loop_kpt=1,num_kpts
          do j=1,num_wann
             do m=1,ndimwin(loop_kpt)
                eigval2(j,loop_kpt)=eigval2(j,loop_kpt)+eigval_opt(m,loop_kpt)* &
                     real(conjg(u_matrix_opt(m,j,loop_kpt))*u_matrix_opt(m,j,loop_kpt),dp)
             enddo
          enddo
       enddo

    else
       eigval2(1:num_wann,:)=eigval(1:num_wann,:)
    end if


    ! At this point eigval2 contains num_wann values which belong to the wannier subspace.


    ! Rotate Hamiltonian into the basis of smooth bloch states 
    !          H(k)=U^{dagger}(k).H_0(k).U(k)
    ! Note: we enforce hermiticity here

    do loop_kpt=1,num_kpts
       do j=1,num_wann
          do i=1,j
             do m=1,num_wann
                ham_k(i,j,loop_kpt)=ham_k(i,j,loop_kpt)+eigval2(m,loop_kpt)* &
                     conjg(u_matrix(m,i,loop_kpt))*u_matrix(m,j,loop_kpt)
             enddo
             if(i.lt.j) ham_k(j,i,loop_kpt)=conjg(ham_k(i,j,loop_kpt))
          enddo
       enddo
    enddo


    ! Fourier transform rotated hamiltonian into WF basis
    ! H_ij(k) --> H_ij(R) = (1/N_kpts) sum_k e^{-ikR} H_ij(k)

    do loop_rpt=1,nrpts
       do loop_kpt=1,num_kpts
          rdotk=twopi*dot_product(kpt_latt(:,loop_kpt),real(irvec(:,loop_rpt),dp))
          fac=exp(-cmplx_i*rdotk)/real(num_kpts,dp)
          ham_r(:,:,loop_rpt)=ham_r(:,:,loop_rpt)+fac*ham_k(:,:,loop_kpt)
       enddo
    enddo
    
    if (timing_level>1) call io_stopwatch('plot: get_hr',2)

  end subroutine plot_get_hr



 !============================================!
  subroutine plot_ham_r
  !============================================!
  !                                            !
  !  Write the Hamiltonian in the WF basis:    !
  !  input for the conductance calculation     !
  !                                            !
  !  Only for 1-dim at this point              ! 
  !  mp_grid should be 1 except one direction  !
  !   - will be checked                        !
  !  *Gamma-point : simply write ham_r         !
  !                                            !
  !  ex) 5 k-pts, in terms of irvec            !
  !      z : zero                              !
  !                                            !
  !     |  0  1  2 |     |  z  z  z |          !
  ! H00=| -1  0  1 | H01=|  2  z  z |          !
  !     | -2 -1  0 |     |  1  2  z |          !      
  !                                            !
  !============================================!

    use w90_constants, only : dp
    use w90_io,        only : io_error,io_stopwatch,io_file_unit, &
                              stdout,seedname
    use w90_parameters, only : mp_grid,num_wann,ham_r_max,timing_level

    implicit none
  
    integer            :: kdir,nrx,file_unit,ierr
    integer            :: i,j,im,jm,m,loop_rpt,counter
    real(kind=dp)      :: ham_rr(num_wann,num_wann), max_hamr(nrpts)
    real(kind=dp), allocatable  :: hr(:,:)

    if (timing_level>1) call io_stopwatch('plot: ham_r',1)

    ! check mp_grid - one dimension?

    counter=0
    do i=1,3
       if ( mp_grid(i) .eq. 1 ) then
          counter=counter+1
       else
          kdir = i
       end if
    end do
   
    if( counter .lt. 2 ) then 
      write(stdout,'(i3,a)') counter,' : NOT A ONE-DIMENSIONAL SYSTEM'
      call io_error('Error in plot_ham_r: incorrect mp_grid') 
    end if

    ! write the  whole matrix with all indices for a crosscheck
    file_unit=io_file_unit()
    open(unit=file_unit,file=trim(seedname)//'.h_all.dat',form='formatted',status='unknown')
    if ( counter .eq. 3 ) then
       do i=1,num_wann
          do j=1,num_wann
             write( file_unit,'(2I5,2F12.6)') j,i,dreal(ham_r(j,i,1)),dimag(ham_r(j,i,1))
          end do
       end do
    else
       do loop_rpt=1,nrpts
          do i=1,num_wann
             do j=1,num_wann
                write( file_unit,'(4I5,2F12.6)') loop_rpt, irvec(kdir,loop_rpt), j,i, &
                                     dreal(ham_r(j,i,loop_rpt)),dimag(ham_r(j,i,loop_rpt))
             end do
          end do
       end do
    end if   
    close(file_unit)

    ! stop here for Gamma sampling

    if (counter .eq. 3 ) return

    ! degeneracy - simply divide by ndegen ( should be fine if ham_r elements are very small at the boundary )

    do loop_rpt=1,nrpts
       if ( ndegen(loop_rpt) .ne. 1 ) ham_r(:,:,loop_rpt) = ham_r(:,:,loop_rpt) / real(ndegen(loop_rpt),dp)
    end do

    ! output largest elements in each lattice point inside W-S cell  
    write(stdout,'(a)') ' Output ham_r for conductance/dos calculation'
    write(stdout,'(a)') ' --check the biggest component in ham_r for each lattice in Wigner-Seitz cell'
    do loop_rpt=1,nrpts
       ham_rr=abs(real(ham_r(:,:,loop_rpt)))
       max_hamr(loop_rpt) = maxval( ham_rr )
       write(stdout,'(I5,"  irvec=",I5,F12.6)') loop_rpt, irvec(kdir,loop_rpt), max_hamr(loop_rpt)   
    end do
    write(stdout,'(a)') ' '
   
    ! H00 will be (nrx+1)*num_wann by (nrx+1)*num_wann

    nrx = nrpts/2

    ! check irvec
    ! irvec should be [-nrx,nrx] in sequence
    m=-nrx
    do loop_rpt=1,nrpts
       if (irvec(kdir,loop_rpt) .ne. m) then
          write(stdout,'(a8,i5,a)') 'nrpts=',nrpts,'YOUR GRID IS NOT THE ONE YOU THINK'
          call io_error('Error in plot_ham_r: wrong guess for nrpts') 
       end if
       m=m+1
    end do 

    ! apply cutoff ham_r_max -> new nrx
   
    nrx = 0
    do loop_rpt=1,nrpts 
       if (max_hamr(loop_rpt) .gt. ham_r_max ) then
          if ( abs(irvec(kdir,loop_rpt)) .gt. nrx )  nrx = abs(irvec(kdir,loop_rpt))
       end if
    end do
    write(stdout,'(a,f10.6)') '  ham_r_max :',ham_r_max
    write(stdout,'(a,i3,a)') '  ham_r truncated at',nrx,'th repeated cell'

    allocate(hr((nrx+1)*num_wann,(nrx+1)*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hr in plot_ham_r')

    file_unit=io_file_unit()
    open(unit=file_unit,file=trim(seedname)//'.h.dat',form='formatted',status='unknown')
  
    ! H00  
    write(file_unit,'(I6)') (nrx+1)*num_wann 
    hr = 0.d0
    do j=0,nrx
       do i=0,nrx  
          m = i-j+nrpts/2+1
          im=i*num_wann
          jm=j*num_wann
          hr(jm+1:jm+num_wann,im+1:im+num_wann)=real(ham_r(:,:,m))
       end do
    end do
    write(file_unit,'(6F15.10)') hr
    
    ! H01  
    write(file_unit,'(I6)') (nrx+1)*num_wann 
    hr = 0.d0
    do j=0,nrx
       do i=0,j-1
          m = i-j+(nrpts/2+1)+(nrx+1) ! starting point shifted by nrx+1  
          im=i*num_wann
          jm=j*num_wann
          hr(jm+1:jm+num_wann,im+1:im+num_wann)=real(ham_r(:,:,m))
       end do
    end do
    write(file_unit,'(6F15.10)') hr

    close(file_unit)

    deallocate(hr,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating hr in plot_ham_r')
    
    if (timing_level>1) call io_stopwatch('plot: ham_r',2)

  end subroutine plot_ham_r



  !================================================================================!
  subroutine wigner_seitz(count_pts)
  !================================================================================!
  ! Calculates a grid of points that fall inside of (and eventually on the         !
  ! surface of) the Wigner-Seitz supercell centered on the origin of the B         !
  ! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3  !
  !================================================================================!

    use w90_constants, only : dp
    use w90_io, only        : stdout,io_error,io_stopwatch
    use w90_parameters, only    : mp_grid,real_metric,iprint,timing_level

  ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
  !                   irvec(1:3,irpt) in the basis of the lattice vectors
  ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
  ! nrpts             number of Wigner-Seitz grid points

  implicit none

  logical, intent(in) :: count_pts 

  integer       :: ndiff (3)
  real(kind=dp) :: dist(125),tot,dist_min
  integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j

  if (timing_level>1) call io_stopwatch('plot: wigner_seitz',1)

  ! The Wannier functions live in a supercell of the real space unit cell
  ! this supercell is mp_grid unit cells long in each direction
  !
  ! We loop over grid points r on a unit cell that is 8 times larger than this
  ! primitive supercell. 
  !
  ! One of these points is in the W-S cell if it is closer to R=0 than any of the
  ! other points, R (where R are the translation vectors of the supercell)

  ! In the end nrpts contains the total number of grid
  ! points that have been found in the Wigner-Seitz cell

  nrpts = 0  
  do n1 = -mp_grid(1) , mp_grid(1)  
     do n2 = -mp_grid(2), mp_grid(2)  
        do n3 = -mp_grid(3),  mp_grid(3)  
           ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=1, or icnt=14
           icnt = 0  
           do i1 = -2, 2  
              do i2 = -2, 2  
                 do i3 = -2, 2  
                    icnt = icnt + 1  
                    ! Calculate distance squared |r-R|^2
                    ndiff(1) = n1 - i1 * mp_grid(1)  
                    ndiff(2) = n2 - i2 * mp_grid(2)  
                    ndiff(3) = n3 - i3 * mp_grid(3)  
                    dist(icnt) = 0.0_dp  
                    do i = 1, 3  
                       do j = 1, 3  
                          dist(icnt) = dist(icnt) + real(ndiff(i),dp) * real_metric(i,j) &
                               * real(ndiff(j),dp)
                       enddo
                    enddo
                 enddo
              enddo


           enddo
           dist_min=minval(dist)
           if (abs(dist(63) - dist_min ) .lt.1.e-7_dp) then
              nrpts = nrpts + 1  
              if(.not. count_pts) then
                 ndegen(nrpts)=0
                do i=1,125
                   if (abs (dist (i) - dist_min) .lt.1.e-7_dp) ndegen(nrpts)=ndegen(nrpts)+1
                end do
                irvec(1, nrpts) = n1  
                irvec(2, nrpts) = n2   
                irvec(3, nrpts) = n3   
              endif
           end if

           !n3
        enddo
        !n2
     enddo
     !n1
  enddo
  !
  if(count_pts) return


  if(iprint>=3) then
     write(stdout,'(1x,i4,a,/)') nrpts,  ' lattice points in Wigner-Seitz supercell:'
     do i=1,nrpts
        write(stdout,'(4x,a,3(i3,1x),a,i2)') '  vector ', irvec(1,i),irvec(2,i),&
             irvec(3,i),'  degeneracy: ', ndegen(i)
     enddo
  endif
  ! Check the "sum rule"
  tot = 0.0_dp  
  do i = 1, nrpts  
     tot = tot + 1.0_dp/real(ndegen(i),dp)  
  enddo
  if (abs (tot - real(mp_grid(1) * mp_grid(2) * mp_grid(3),dp) ) > 1.e-8_dp) then
     call io_error('ERROR in plot_wigner_seitz: error in finding Wigner-Seitz points')
  endif

  if (timing_level>1) call io_stopwatch('plot: wigner_seitz',2)

  return  
end subroutine wigner_seitz




end module w90_plot
 
