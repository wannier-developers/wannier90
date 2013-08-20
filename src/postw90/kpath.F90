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

module w90_kpath

  ! Calculates one of the following along a specified k-path:
  ! 
  !  - Energy bands (eventually colored by the spin) 
  !
  !  - (Berry curvature)x(-1) summed over occupied bands
  !
  !  - Integrand of orbital magnetization Morb=LCtil+ICtil

  use w90_constants, only : dp

  implicit none

  public

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine k_path

    use w90_comms
    use w90_constants,  only     : dp,cmplx_0,cmplx_i,twopi,eps8
    use w90_io,         only     : io_error,io_file_unit,seedname,&
                                   io_time,io_stopwatch,stdout
    use w90_utility, only        : utility_diagonalize
    use w90_postw90_common, only : fourier_R_to_k
    use w90_parameters, only     : num_wann,recip_metric,kpath_task,&
                                   kpath_num_points,bands_num_spec_points,&
                                   bands_spec_points,bands_label,&
                                   kpath_bands_colour,nfermi,fermi_energy_list,&
                                   berry_curv_unit
    use w90_get_oper, only       : get_HH_R,HH_R,get_AA_R,get_BB_R,get_CC_R,&
                                   get_FF_R,get_SS_R
    use w90_spin, only           : get_spin_nk
    use w90_berry, only          : get_imf_k_list,get_imfgh_k_list
    use w90_constants, only      : bohr

    integer           :: i,j,n,num_paths,num_spts,loop_path,loop_kpt,&
                         total_pts,counter,loop_i,dataunit,gnuunit,pyunit,&
                         kpath_pts(bands_num_spec_points/2)
    real(kind=dp)     :: ymin,ymax,vec(3),kpt(3),spn_k(num_wann),&
                         imf_k_list(3,3,nfermi),img_k_list(3,3,nfermi),&
                         imh_k_list(3,3,nfermi),Morb_k(3,3),&
                         kpath_len(bands_num_spec_points/2),range
    logical           :: plot_bands,plot_curv,plot_morb
    character(len=20) :: file_name

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    real(kind=dp), allocatable    :: xval(:),eig(:,:),curv(:,:),& 
                                     morb(:,:),color(:,:),&
                                     plot_kpoint(:,:) 
    character(len=3),allocatable  :: glabel(:)

    ! Everything is done on the root node (not worthwhile parallelizing) 
    ! However, we still have to read and distribute the data if we 
    ! are in parallel. So calls to get_oper are done on all nodes at the moment
    !
    plot_bands=.false.
    if(index(kpath_task,'bands')>0) then
       plot_bands=.true.
    end if
    plot_curv=.false.
    if(index(kpath_task,'curv')>0) then
       plot_curv=.true.
    end if
    plot_morb=.false.
    if(index(kpath_task,'morb')>0) then
       plot_morb=.true.
    end if
    ! Set up the needed Wannier matrix elements
    call get_HH_R
    if(plot_curv.or.plot_morb) then
       call get_AA_R
    endif
    if(plot_morb) then
       call get_BB_R
       call get_CC_R
    endif
    if(plot_bands .and. kpath_bands_colour=='spin') call get_SS_R

    if(on_root) then

       write(stdout,'(/,/,1x,a)')&
            'Properties calculated in module  k p a t h'
       write(stdout,'(1x,a)')&
            '------------------------------------------'

       if(plot_bands) then
          select case(kpath_bands_colour)
          case("none")
             write(stdout,'(/,3x,a)') '* Energy bands in eV'
          case("spin")
             write(stdout,'(/,3x,a)') '* Energy bands in eV, coloured by spin'
          end select
       end if
       if(plot_curv) then
          if(berry_curv_unit=='ang2') then
             write(stdout,'(/,3x,a)') '* Negative Berry curvature in Ang^2'
          elseif(berry_curv_unit=='bohr2') then
             write(stdout,'(/,3x,a)') '* Negative Berry curvature in Bohr^2'
          endif
          if(nfermi/=1) call io_error('Need to specify one value of '&
               //'the fermi energy when kpath_task=curv')
       end if
       if(plot_morb) then
          write(stdout,'(/,3x,a)')& 
               '* Orbital magnetization k-space integrand in eV.Ang^2'
          if(nfermi/=1) call io_error('Need to specify one value of '&
               //'the fermi energy when kpath_task=morb')
       end if

       ! Work out how many points there are in the total path, and the 
       ! positions of the special points
       !
       num_paths=bands_num_spec_points/2 ! number of straigh line segments 
       ! (paths)
       num_spts=num_paths+1 ! number of path endpoints (special pts)
       do loop_path=1,num_paths
          vec=bands_spec_points(:,2*loop_path)&
               -bands_spec_points(:,2*loop_path-1)
          kpath_len(loop_path)=sqrt(&
               dot_product(vec,(matmul(recip_metric,vec)))&
               )
          !
          ! kpath_pts(loop_path) is the number of points in path number 
          ! loop_path (all segments have the same density of points)
          !
          if(loop_path==1) then
             kpath_pts(loop_path)=kpath_num_points
          else 
             kpath_pts(loop_path)=nint(real(kpath_num_points,dp)&
                  *kpath_len(loop_path)/kpath_len(1))
          end if
       end do
       total_pts=sum(kpath_pts)+1

       ! Reciprocal-lattice coordinates of the k-points along the path
       !
       allocate(plot_kpoint(3,total_pts))

       ! Value of the horizontal coordinate in the actual plots (units of 
       ! distance in k-space)
       !
       allocate(xval(total_pts))

       ! Value of the vertical coordinate in the actual plots: energy bands 
       !
       if(plot_bands) then
          allocate(HH(num_wann,num_wann))
          allocate(UU(num_wann,num_wann))
          allocate(eig(num_wann,total_pts))
          if(kpath_bands_colour/='none') allocate(color(num_wann,total_pts))
       end if

       ! Value of the vertical coordinate in the actual plots
       !
       if(plot_curv) allocate(curv(total_pts,3)) 
       if(plot_morb) allocate(morb(total_pts,3))

       allocate(glabel(num_spts))

       ! Find the position of each kpoint along the path
       !
       counter=0
       do loop_path=1,num_paths
          do loop_i=1,kpath_pts(loop_path)
             counter=counter+1
             if(counter==1) then
                xval(counter)=0.0_dp
             else
                xval(counter)=xval(counter-1)&
                     +kpath_len(loop_path)/real(kpath_pts(loop_path),dp)
             endif
             plot_kpoint(:,counter)=bands_spec_points(:,2*loop_path-1)&
                  +( bands_spec_points(:,2*loop_path)&
                  -bands_spec_points(:,2*loop_path-1)&
                  )&
                  *(real(loop_i-1,dp)/real(kpath_pts(loop_path),dp))
          end do
       end do
       !
       ! Last point
       !
       xval(total_pts)=sum(kpath_len)
       plot_kpoint(:,total_pts)=bands_spec_points(:,bands_num_spec_points)

       if(plot_bands) then
          !
          ! Write out the kpoints in the path in a format that can be inserted
          ! directly in the pwscf input file (the '1.0_dp' in the second column is 
          ! a k-point weight, expected by pwscf)
          !
          dataunit=io_file_unit()
          open(dataunit,file=trim(seedname)//'-path.kpt',form='formatted')
          write(dataunit,*) total_pts
          do loop_kpt=1,total_pts
             write(dataunit,'(3f12.6,3x,f4.1)')&
                  (plot_kpoint(loop_i,loop_kpt),loop_i=1,3),1.0_dp
          end do
          close(dataunit)
       endif

       ! Loop over k-points on the path and evaluate the requested quantities
       !
       do loop_kpt=1,total_pts       
          kpt(:)=plot_kpoint(:,loop_kpt)

          if(plot_bands) then
             call fourier_R_to_k(kpt,HH_R,HH,0)
             call utility_diagonalize(HH,num_wann,eig(:,loop_kpt),UU)
             !
             ! Color-code energy bands with the spin projection along the
             ! chosen spin quantization axis
             !
             if(kpath_bands_colour=='spin') then
                call get_spin_nk(kpt,spn_k)
                color(:,loop_kpt)=spn_k(:)
                !
                ! The following is needed to prevent bands from disappearing 
                ! when the magnitude of the Wannier interpolated spn_k (very 
                ! slightly) exceeds 1.0 (e.g. in bcc Fe along N--G--H)
                !
                do n=1,num_wann
                   if(color(n,loop_kpt)>1.0_dp-eps8) then
                      color(n,loop_kpt)=1.0_dp-eps8
                   elseif(color(n,loop_kpt)<-1.0_dp+eps8) then
                      color(n,loop_kpt)=-1.0_dp+eps8
                   endif
                enddo
             end if
          endif

          if(plot_curv) then
             call get_imf_k_list(kpt,imf_k_list)
             curv(loop_kpt,1)=sum(imf_k_list(:,1,1))
             curv(loop_kpt,2)=sum(imf_k_list(:,2,1))
             curv(loop_kpt,3)=sum(imf_k_list(:,3,1))
          end if

          if(plot_morb) then
             call get_imfgh_k_list(kpt,imf_k_list,img_k_list,imh_k_list)
             Morb_k=img_k_list(:,:,1)+imh_k_list(:,:,1)&
                   -2.0_dp*fermi_energy_list(1)*imf_k_list(:,:,1)
             Morb_k=-Morb_k/2.0_dp ! differs by -1/2 from Eq.97 LVTS12
             morb(loop_kpt,1)=sum(Morb_k(:,1))
             morb(loop_kpt,2)=sum(Morb_k(:,2))
             morb(loop_kpt,3)=sum(Morb_k(:,3))
          end if

       end do !loop_kpt

       if(plot_curv .and. berry_curv_unit=='bohr2') curv=curv/bohr**2

       ! Axis labels
       !
       glabel(1)=' '//bands_label(1)//' '
       do i=2,num_paths
          if(bands_label(2*(i-1))/=bands_label(2*(i-1)+1)) then
             glabel(i)=bands_label(2*(i-1))//'/'//bands_label(2*(i-1)+1)
          else
             glabel(i)=' '//bands_label(2*(i-1))//' '
          end if
       end do
       glabel(num_spts)=' '//bands_label(bands_num_spec_points)//' '

       ! Now write the plotting files

       write(stdout,'(/,1x,a)')     'Output files:'

       if(plot_bands) then
          file_name=trim(seedname)//'-path.kpt'
          write(stdout,'(/,3x,a)') file_name
          !
          ! Data file
          !
          dataunit=io_file_unit()
          file_name=trim(seedname)//'-bands.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do i=1,num_wann
             do loop_kpt=1,total_pts
                if(kpath_bands_colour=='none') then
                   write(dataunit,'(2E16.8)') xval(loop_kpt),eig(i,loop_kpt)
                else
                   write(dataunit,'(3E16.8)') xval(loop_kpt),&
                        eig(i,loop_kpt),color(i,loop_kpt)  
                end if
             enddo
             write(dataunit,*) ' '
          enddo
          close(dataunit)
       endif

       if(plot_bands .and. .not.plot_curv .and. .not.plot_morb) then
          !
          ! Gnuplot script
          !
          ymin=minval(eig)-1.0_dp
          ymax=maxval(eig)+1.0_dp
          gnuunit=io_file_unit()
          file_name=trim(seedname)//'-bands.gnu'
          write(stdout,'(/,3x,a)') file_name
          open(gnuunit,file=file_name,form='formatted')
          do i = 1,num_paths-1
             write(gnuunit,705) sum(kpath_len(1:i)),ymin,&
                  sum(kpath_len(1:i)),ymax
          enddo
          if(kpath_bands_colour=='none') then
             write(gnuunit,701) xval(total_pts),ymin,ymax
             write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
                  (glabel(i+1),sum(kpath_len(1:i)),i=1,num_paths-1)
             write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
             write(gnuunit,*) 'plot ','"'//trim(seedname)//'-bands.dat','"' 
          elseif(kpath_bands_colour=='spin') then
             !
             ! Only works with gnuplot v4.2 and higher
             !
             write(gnuunit,706) xval(total_pts),ymin,ymax
             write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
                  (glabel(i+1),sum(kpath_len(1:i)),i=1,num_paths-1)
             write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
             write(gnuunit,*)&
                  'set palette defined (-1 "blue", 0 "green", 1 "red")'
             write(gnuunit,*) 'set pm3d map'
             write(gnuunit,*) 'set zrange [-1:1]'
             write(gnuunit,*) 'splot ','"'//trim(seedname)//'-bands.dat',& 
                  '" with dots palette' 
          end if
          close(gnuunit)
          !
          ! python script
          !
          pyunit=io_file_unit()
          file_name=trim(seedname)//'-bands.py'
          write(stdout,'(/,3x,a)') file_name
          open(pyunit,file=file_name,form='formatted')             
          write(pyunit,'(a)') 'import pylab as pl'
          write(pyunit,'(a)') 'import numpy as np'
          write(pyunit,'(a)') "data = np.loadtxt('"//trim(seedname)//&
               "-bands.dat')"
          write(pyunit,'(a)') "x=data[:,0]"
          write(pyunit,'(a)') "y=data[:,1]"
          if(kpath_bands_colour=='spin') write(pyunit,'(a)') "z=data[:,2]"
          write(pyunit,'(a)') "tick_labels=[]"
          write(pyunit,'(a)') "tick_locs=[]"
          do j=1,num_spts
             if(trim(glabel(j))==' G') then
                write(pyunit,'(a)') "tick_labels.append('$\Gamma$')"
             else
                write(pyunit,'(a)') "tick_labels.append('"//trim(glabel(j))&
                     //"'.strip())"
             endif
             if(j==1) then
                write(pyunit,'(a,F12.6,a)') "tick_locs.append(0)"
             else
                write(pyunit,'(a,F12.6,a)') "tick_locs.append(",&
                     sum(kpath_len(1:j-1)),")"
             endif
          enddo
          if(kpath_bands_colour=='none') then
             write(pyunit,'(a)') "pl.scatter(x,y,color='k',marker='+',s=0.1)"
          elseif(kpath_bands_colour=='spin') then
             write(pyunit,'(a)')&
                  "pl.scatter(x,y,c=z,marker='+',s=1,cmap=pl.cm.jet)"
          endif
          write(pyunit,'(a)') "pl.xlim([0,max(x)])"
          write(pyunit,'(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y)),"&
               //"max(y)+0.025*(max(y)-min(y))])"
          write(pyunit,'(a)') "pl.xticks(tick_locs,tick_labels)"
          write(pyunit,'(a)') "for n in range(1,len(tick_locs)):"
          write(pyunit,'(a)') "   pl.plot([tick_locs[n],tick_locs[n]],"&
               //"[pl.ylim()[0],pl.ylim()[1]],color='gray',"&
               //"linestyle='-',linewidth=0.5)"
          write(pyunit,'(a)') "pl.ylabel('Energy [eV]')"
          if(kpath_bands_colour=='spin') then
             write(pyunit,'(a)')&
                  "pl.axes().set_aspect(aspect=0.65*max(x)/(max(y)-min(y)))"
             write(pyunit,'(a)') "pl.colorbar(shrink=0.7)"
          endif
          write(pyunit,'(a)') "outfile = '"//trim(seedname)//"-bands.pdf'"
          write(pyunit,'(a)') "pl.savefig(outfile)"
          write(pyunit,'(a)') "pl.show()"
          
       endif ! plot_bands .and. .not.plot_curv .and. .not.plot_morb

       if(plot_curv) then
          ! It is conventional to plot the negative curvature
          curv=-curv
          dataunit=io_file_unit()
          file_name=trim(seedname)//'-curv.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do loop_kpt=1,total_pts
             write(dataunit,'(4E16.8)') xval(loop_kpt),&
                  curv(loop_kpt,:) 
          enddo
          write(dataunit,*) ' '
          close(dataunit)
       endif

       if(plot_curv .and. .not.plot_bands) then

          do i=1,3
             !   
             ! gnuplot script
             !
             gnuunit=io_file_unit()
             file_name=trim(seedname)//'-curv_'//achar(119+i)//'.gnu'
             write(stdout,'(/,3x,a)') file_name
             open(gnuunit,file=file_name,form='formatted')
             ymin=minval(curv(:,i))
             ymax=maxval(curv(:,i))
             range=ymax-ymin
             ymin=ymin-0.02_dp*range
             ymax=ymax+0.02_dp*range
             write(gnuunit,707) xval(total_pts),ymin,ymax
             do j=1,num_paths-1
                write(gnuunit,705) sum(kpath_len(1:j)),ymin,&
                                   sum(kpath_len(1:j)),ymax
             enddo
             write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
                  (glabel(j+1),sum(kpath_len(1:j)),j=1,num_paths-1)
             write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
             write(gnuunit,*)&
                  'plot ','"'//trim(seedname)//'-curv.dat','" u 1:'//achar(49+i)
             close(gnuunit)
             !   
             ! python script
             !
             pyunit=io_file_unit()
             file_name=trim(seedname)//'-curv_'//achar(119+i)//'.py'
             write(stdout,'(/,3x,a)') file_name
             open(pyunit,file=file_name,form='formatted')             
             write(pyunit,'(a)') 'import pylab as pl'
             write(pyunit,'(a)') 'import numpy as np'
             write(pyunit,'(a)') "data = np.loadtxt('"//trim(seedname)//&
                  "-curv.dat')"
             write(pyunit,'(a)') "x=data[:,0]"
             write(pyunit,'(a)') "y=data[:,"//achar(48+i)//"]"
             write(pyunit,'(a)') "tick_labels=[]"
             write(pyunit,'(a)') "tick_locs=[]"
             do j=1,num_spts
                if(trim(glabel(j))==' G') then
                   write(pyunit,'(a)') "tick_labels.append('$\Gamma$')"
                else
                   write(pyunit,'(a)') "tick_labels.append('"&
                        //trim(glabel(j))//"'.strip())"
                endif
                if(j==1) then
                   write(pyunit,'(a,F12.6,a)') "tick_locs.append(0)"
                else
                   write(pyunit,'(a,F12.6,a)') "tick_locs.append(",&
                        sum(kpath_len(1:j-1)),")"
                endif
             enddo
             write(pyunit,'(a)') "pl.plot(x,y,color='k')"
             write(pyunit,'(a)') "pl.xlim([0,max(x)])"
             write(pyunit,'(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y)),"&
                  //"max(y)+0.025*(max(y)-min(y))])"
             write(pyunit,'(a)') "pl.xticks(tick_locs,tick_labels)"
             write(pyunit,'(a)') "for n in range(1,len(tick_locs)):"
             write(pyunit,'(a)') "   pl.plot([tick_locs[n],tick_locs[n]],"&
                  //"[pl.ylim()[0],pl.ylim()[1]],color='gray',"&
                  //"linestyle='-',linewidth=0.5)"
             if(berry_curv_unit=='ang2') then
                write(pyunit,'(a)') "pl.ylabel('$-\Omega_"//achar(119+i)&
                     //"(\mathbf{k})$  [ $\AA^2$ ]')"
             elseif(berry_curv_unit=='bohr2') then
                write(pyunit,'(a)') "pl.ylabel('$-\Omega_"//achar(119+i)&
                     //"(\mathbf{k})$  [ bohr$^2$ ]')"
             endif
             write(pyunit,'(a)') "outfile = '"//trim(seedname)//&
                  "-curv_"//achar(119+i)//".pdf'"
             write(pyunit,'(a)') "pl.savefig(outfile)"
             write(pyunit,'(a)') "pl.show()"
          enddo
          
       end if ! plot_curv .and. .not.plot_bands

       if(plot_morb) then
          dataunit=io_file_unit()
          file_name=trim(seedname)//'-morb.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do loop_kpt=1,total_pts
             write(dataunit,'(4E16.8)') xval(loop_kpt),morb(loop_kpt,:)
          enddo
          write(dataunit,*) ' '
          close(dataunit)
       endif

       if(plot_morb .and. .not.plot_bands) then
          do i=1,3
             !
             ! gnuplot script
             !
             gnuunit=io_file_unit()
             file_name=trim(seedname)//'-morb_'//achar(119+i)//'.gnu'
             write(stdout,'(/,3x,a)') file_name
             open(gnuunit,file=file_name,form='formatted')
             ymin=minval(morb(:,i))
             ymax=maxval(morb(:,i))
             range=ymax-ymin
             ymin=ymin-0.02_dp*range
             ymax=ymax+0.02_dp*range
             write(gnuunit,707) xval(total_pts),ymin,ymax
             do j=1,num_paths-1
                write(gnuunit,705) sum(kpath_len(1:j)),ymin,&
                                   sum(kpath_len(1:j)),ymax
             enddo
             write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
                  (glabel(j+1),sum(kpath_len(1:j)),j=1,num_paths-1)
             write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
             write(gnuunit,*)&
                  'plot ','"'//trim(seedname)//'-morb.dat','" u 1:'//achar(49+i)
             close(gnuunit)
             !   
             ! python script
             !
             pyunit=io_file_unit()
             file_name=trim(seedname)//'-morb_'//achar(119+i)//'.py'
             write(stdout,'(/,3x,a)') file_name
             open(pyunit,file=file_name,form='formatted')             
             write(pyunit,'(a)') 'import pylab as pl'
             write(pyunit,'(a)') 'import numpy as np'
             write(pyunit,'(a)') "data = np.loadtxt('"//trim(seedname)//&
                  "-morb.dat')"
              write(pyunit,'(a)') "x=data[:,0]"
              write(pyunit,'(a)') "y=data[:,"//achar(48+i)//"]"
              write(pyunit,'(a)') "tick_labels=[]"
              write(pyunit,'(a)') "tick_locs=[]"
              do j=1,num_spts
                 if(trim(glabel(j))==' G') then
                    write(pyunit,'(a)') "tick_labels.append('$\Gamma$')"
                 else
                    write(pyunit,'(a)')&
                         "tick_labels.append('"//trim(glabel(j))//"'.strip())"
                 endif
                 if(j==1) then
                    write(pyunit,'(a,F12.6,a)') "tick_locs.append(0)"
                 else
                    write(pyunit,'(a,F12.6,a)') "tick_locs.append(",&
                         sum(kpath_len(1:j-1)),")"
                 endif
              enddo
              write(pyunit,'(a)') "pl.plot(x,y,color='k')"
              write(pyunit,'(a)') "pl.xlim([0,max(x)])"
              write(pyunit,'(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y)),"&
                                         //"max(y)+0.025*(max(y)-min(y))])"
              write(pyunit,'(a)') "pl.xticks(tick_locs,tick_labels)"
              write(pyunit,'(a)') "for n in range(1,len(tick_locs)):"
              write(pyunit,'(a)') "   pl.plot([tick_locs[n],tick_locs[n]],"&
                   //"[pl.ylim()[0],pl.ylim()[1]],color='gray',"&
                   //"linestyle='-',linewidth=0.5)"
              write(pyunit,'(a)') "pl.ylabel(r'$M^{\rm{orb}}_z(\mathbf{k})$"&
                   //"  [ Ry$\cdot\AA^2$ ]')"
              write(pyunit,'(a)') "outfile = '"//trim(seedname)//&
                   "-morb_"//achar(119+i)//".pdf'"
              write(pyunit,'(a)') "pl.savefig(outfile)"
              write(pyunit,'(a)') "pl.show()"
          enddo

       end if ! plot_morb .and. .not.plot_bands

       if(plot_bands .and. (plot_curv .or. plot_morb)) then
          !
          ! python script
          !
          do i=1,3
             pyunit=io_file_unit()
             if(plot_curv) then
                file_name=trim(seedname)//'-bands+curv_'//achar(119+i)//'.py'
             elseif(plot_morb) then
                file_name=trim(seedname)//'-bands+morb_'//achar(119+i)//'.py'
             endif
             write(stdout,'(/,3x,a)') file_name
             open(pyunit,file=file_name,form='formatted')             
             write(pyunit,'(a)') 'import pylab as pl'
             write(pyunit,'(a)') 'import numpy as np'
             write(pyunit,'(a)') 'from matplotlib.gridspec import GridSpec'
             write(pyunit,'(a)') "tick_labels=[]"
             write(pyunit,'(a)') "tick_locs=[]"
             do j=1,num_spts
                if(trim(glabel(j))==' G') then
                   write(pyunit,'(a)') "tick_labels.append('$\Gamma$')"
                else
                   write(pyunit,'(a)') "tick_labels.append('"//trim(glabel(j))&
                        //"'.strip())"
                endif
                if(j==1) then
                   write(pyunit,'(a,F12.6,a)') "tick_locs.append(0)"
                else
                   write(pyunit,'(a,F12.6,a)') "tick_locs.append(",&
                        sum(kpath_len(1:j-1)),")"
                endif
             enddo
             write(pyunit,'(a)') "fig = pl.figure()"
             write(pyunit,'(a)') "gs = GridSpec(2, 1,hspace=0.00)"
             !
             ! upper panel (energy bands)
             !
             write(pyunit,'(a)') "axes1 = pl.subplot(gs[0, 0:])"
             write(pyunit,'(a)') "data = np.loadtxt('"//trim(seedname)//&
                  "-bands.dat')"
             write(pyunit,'(a)') "x=data[:,0]"
             write(pyunit,'(a,F12.6)') "y=data[:,1]-",fermi_energy_list(1)
             if(kpath_bands_colour=='spin') write(pyunit,'(a)') "z=data[:,2]"
             if(kpath_bands_colour=='none') then
                write(pyunit,'(a)') "pl.scatter(x,y,color='k',marker='+',s=0.1)"
             elseif(kpath_bands_colour=='spin') then
                write(pyunit,'(a)')&
                     "pl.scatter(x,y,c=z,marker='+',s=1,cmap=pl.cm.jet)"
             endif
             write(pyunit,'(a)') "pl.xlim([0,max(x)])"
             write(pyunit,'(a)') "pl.ylim([-0.65,0.65]) # Adjust this range as needed"
             write(pyunit,'(a)') "pl.plot([tick_locs[0],tick_locs[-1]],[0,0],"&
                  //"color='black',linestyle='--',linewidth=0.5)"
             write(pyunit,'(a)') "pl.xticks(tick_locs,tick_labels)"
             write(pyunit,'(a)') "for n in range(1,len(tick_locs)):"
             write(pyunit,'(a)') "   pl.plot([tick_locs[n],tick_locs[n]],"&
                  //"[pl.ylim()[0],pl.ylim()[1]],color='gray',"&
                  //"linestyle='-',linewidth=0.5)"
             write(pyunit,'(a)') "pl.ylabel('Energy$-$E$_F$ [eV]')"
             write(pyunit,'(a)') "pl.tick_params(axis='x',"&
                  //"which='both',bottom='off',top='off',labelbottom='off')"
             !
             ! lower panel (curvature or orbital magnetization)
             !
             write(pyunit,'(a)') "axes2 = pl.subplot(gs[1, 0:])"
             if(plot_curv) then
                write(pyunit,'(a)') "data = np.loadtxt('"//trim(seedname)//&
                     "-curv.dat')"
             elseif(plot_morb) then
                write(pyunit,'(a)') "data = np.loadtxt('"//trim(seedname)//&
                     "-morb.dat')"
             endif
             write(pyunit,'(a)') "x=data[:,0]"
             write(pyunit,'(a)') "y=data[:,"//achar(48+i)//"]"
             write(pyunit,'(a)') "pl.plot(x,y,color='k')"
             write(pyunit,'(a)') "pl.xlim([0,max(x)])"
             write(pyunit,'(a)') "pl.ylim([min(y)-0.025*(max(y)-min(y)),"&
                  //"max(y)+0.025*(max(y)-min(y))])"
             write(pyunit,'(a)') "pl.xticks(tick_locs,tick_labels)"
             write(pyunit,'(a)') "for n in range(1,len(tick_locs)):"
             write(pyunit,'(a)') "   pl.plot([tick_locs[n],tick_locs[n]],"&
                  //"[pl.ylim()[0],pl.ylim()[1]],color='gray',"&
                  //"linestyle='-',linewidth=0.5)"
             if(plot_curv) then
                if(berry_curv_unit=='ang2') then
                   write(pyunit,'(a)') "pl.ylabel('$-\Omega_"//achar(119+i)&
                        //"(\mathbf{k})$  [ $\AA^2$ ]')"
                elseif(berry_curv_unit=='bohr2') then
                   write(pyunit,'(a)') "pl.ylabel('$-\Omega_"//achar(119+i)&
                        //"(\mathbf{k})$  [ bohr$^2$ ]')"
                endif
                write(pyunit,'(a)') "outfile = '"//trim(seedname)//&
                     "-bands+curv_"//achar(119+i)//".pdf'"
             elseif(plot_morb) then
                write(pyunit,'(a)') "pl.ylabel(r'$M^{\rm{orb}}_z(\mathbf{k})$"&
                     //"  [ Ry$\cdot\AA^2$ ]')"
                write(pyunit,'(a)') "outfile = '"//trim(seedname)//&
                     "-morb_"//achar(119+i)//".pdf'"
             endif
             write(pyunit,'(a)') "pl.savefig(outfile)"
             write(pyunit,'(a)') "pl.show()"
          enddo

       endif ! plot_bands .and. plot_curv


    end if ! on_root

701 format('set style data dots',/,'unset key',/,&
         'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
702 format('set xtics (',:20('"',A3,'" ',F8.5,','))
703 format(A3,'" ',F8.5,')')
704 format('set palette defined (',F8.5,' "red", 0 "green", ',F8.5,' "blue")')
705 format('set arrow from ',F16.8,',',F16.8,' to ',F16.8,',',F16.8, ' nohead')
706 format('unset key',/,&
         'set xrange [0:',F9.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
707 format('set style data lines',/,'set nokey',/,&
         'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')

  end subroutine k_path

end module w90_kpath
