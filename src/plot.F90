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


module plot
 
  use constants, only : dp
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

    use constants, only : cmplx_0
    use io, only        : io_error,stdout
    use parameters, only    : num_kpts,bands_plot,dos_plot,&
         fermi_surface_plot,num_wann

    implicit none

    integer :: ierr

    write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,'(1x,a)') '|                               PLOTTING                                    |'
    write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(stdout,*)

    if(bands_plot .or. dos_plot .or. fermi_surface_plot) then
       allocate(irvec(3,3*num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating irvec in plot_main')
       irvec=0
       allocate(ndegen(3*num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating ndegen in plot_main')
       ndegen=0
       call wigner_seitz
       allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating ham_r in plot_main')
       ham_r=cmplx_0
       call plot_get_hr
       !
       if(bands_plot) call plot_interpolate_bands
       !
       if(fermi_surface_plot) call plot_fermi_surface
       !
       deallocate(ham_r,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating ham_r in plot_main')
       deallocate(ndegen,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating ndegen in plot_main')
       deallocate(irvec,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating irvec in plot_main')
       !
    end if

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

    use constants, only : dp,cmplx_0,cmplx_i,twopi
    use io, only        : io_error,stdout,io_file_unit,seedname,io_time
    use parameters, only    : num_wann,bands_num_points,&
         recip_metric,bands_num_spec_points, &
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
       call ZHPEVX('N','A','U',num_wann,ham_pack,0.d0,0.d0,0,0,-1.d0, &
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
          write(bndunit,'(2E16.8)') XVAL(NKP),eig_int(i,nkp)
       enddo
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
    do i = 1, num_paths
       write(gnuunit,705) sum(kpath_len(1:i)),emin,sum(kpath_len(1:i)),emax
    enddo
    write(gnuunit,702, advance="no") glabel(1),0.d0,(glabel(i+1),sum(kpath_len(1:i)),i=1,bands_num_spec_points/2-1)
    write(gnuunit,703) glabel(1+bands_num_spec_points/2),sum(kpath_len(:))
    write(gnuunit,*) 'plot ','"'//trim(seedname)//'_band.dat','"' 
    close(gnuunit)

    write(stdout,'(1x,a,f11.3,a)')  'Time to calculate interpolated band structure ',io_time()-time0,' (sec)'
    write(stdout,*)

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

    use constants, only : dp,cmplx_0,cmplx_i,twopi
    use io, only        : io_error,stdout,io_file_unit,seedname,io_date,io_time
    use parameters, only    : num_wann,fermi_surface_num_points,&
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
             call ZHPEVX('N','A','U',num_wann,ham_pack,0.d0,0.d0,0,0,-1.d0, &
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
    write(bxsf_unit,*) '      # Generated by the Wannier code http://www.wannier.org'
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


  end subroutine plot_fermi_surface

 
 !============================================!
  subroutine plot_get_hr
  !============================================!
  !                                            !
  !  Calculate the Hamiltonian in the WF basis !
  !                                            !
  !============================================!

    use constants, only : dp,cmplx_0,cmplx_i,twopi
    use io, only        : io_error
    use parameters, only    : num_kpts,u_matrix,eigval_opt,num_wann,kpt_latt

    implicit none
  
    complex(kind=dp)   :: ham_k(num_wann,num_wann,num_kpts)
    complex(kind=dp)   :: fac
    real(kind=dp)      :: rdotk
    integer            :: loop_kpt,i,j,m,loop_rpt


    ham_k=cmplx_0

    ! Rotate Hamiltonian into the basis of smooth bloch states 
    !          H(k)=U^{dagger}(k).H_0(k).U(k)
    ! Note: we enforce hermiticity here
 
    do loop_kpt=1,num_kpts
       do j=1,num_wann
          do i=1,j
             do m=1,num_wann
                ham_k(i,j,loop_kpt)=ham_k(i,j,loop_kpt)+eigval_opt(m,loop_kpt)* &
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
    

  end subroutine plot_get_hr


  !================================================================================!
  subroutine wigner_seitz( )
  !================================================================================!
  ! Calculates a grid of points that fall inside of (and eventually on the         !
  ! surface of) the Wigner-Seitz supercell centered on the origin of the B         !
  ! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3  !
  !================================================================================!

    use constants, only : dp
    use io, only        : stdout,io_error
    use parameters, only    : mp_grid,real_metric,num_kpts,iprint

  ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
  !                   irvec(1:3,irpt) in the basis of the lattice vectors
  ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
  ! nrpts             number of Wigner-Seitz grid points

  implicit none

  integer       :: ndiff (3)
  real(kind=dp) :: dist(27),tot,dist_min
  integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j


  ! Loop over grid points r on a unit cell that is 8 times larger than a
  ! primitive supercell. In the end nrpts contains the total number of grid
  ! points that have been found in the Wigner-Seitz cell

  nrpts = 0  
  do n1 = 0, 2 * mp_grid(1)  
     do n2 = 0, 2 * mp_grid(2)  
        do n3 = 0, 2 * mp_grid(3)  
           ! Loop over the 27 points R. R=0 corresponds to i1=i2=i3=1, or icnt=14
           icnt = 0  
           do i1 = 0, 2  
              do i2 = 0, 2  
                 do i3 = 0, 2  
                    icnt = icnt + 1  
                    ! Calculate distance squared |r-R|^2
                    ndiff(1) = n1 - i1 * mp_grid(1)  
                    ndiff(2) = n2 - i2 * mp_grid(2)  
                    ndiff(3) = n3 - i3 * mp_grid(3)  
                    dist(icnt) = 0.d0  
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
           if (abs(dist(14) - dist_min ) .lt.1.d-7) then
              nrpts = nrpts + 1  
              ! Hopefully this will never happen, i.e., I think 3*nkpts is
              ! an upper bound to the number of lattice points in (or on
              ! the surface of) the Wigner-Seitz supercell
              if (nrpts.gt.3 * num_kpts) then  
                 call io_error('Error: nrpts too small in plot_wigner_seitz &
                      & edit plot.f90 and contact code developers')
              endif
              ndegen(nrpts)=0
              do i=1,27
                 if (abs (dist (i) - dist_min) .lt.1.d-7) ndegen(nrpts)=ndegen(nrpts)+1
              end do

              irvec(1, nrpts) = n1 - mp_grid(1)  
              irvec(2, nrpts) = n2 - mp_grid(2)  
              irvec(3, nrpts) = n3 - mp_grid(3)  

           end if

           !n3
        enddo
        !n2
     enddo
     !n1
  enddo
  !
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
  if (abs (tot - real(mp_grid(1) * mp_grid(2) * mp_grid(3),dp) ) .gt.1.d-8) then
     call io_error('ERROR in plot_wigner_seitz: error in finding Wigner-Seitz points')
  endif


  return  
end subroutine wigner_seitz




  end module plot
 
