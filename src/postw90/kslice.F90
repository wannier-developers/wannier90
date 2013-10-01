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
!  !                                                          !
!------------------------------------------------------------!

module w90_kslice

  ! Makes a heatmap plot on a slice in k-space of:
  ! 
  !  - The negative Berry curvature summed over occupied bands 
  !
  !  - The integrand of the k-space orbital magnetization formula
  !
  ! Plots the intersections of constant-energy isosurfaces with the slice
  !
  ! The slice is defined by three input variables, all in reciprocal
  ! lattice coordinates:
  !
  !    kslice_corner(1:3) is the lower left corner 
  !    kslice_b1(1:3) and kslice_b2(1:3) are the vectors subtending the slice
  !
  !---------------------------------
  ! TO DO: Parallelize over k-points
  !---------------------------------

  implicit none

  public

  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine k_slice

    use w90_comms
    use w90_constants,  only     : dp,twopi,eps8
    use w90_io,         only     : io_error,io_file_unit,seedname,&
                                   io_time,io_stopwatch,stdout
    use w90_utility, only        : utility_diagonalize
    use w90_postw90_common, only : fourier_R_to_k
    use w90_parameters, only     : num_wann,kslice,kslice_task,&
                                   kslice_2dkmesh,kslice_corner,kslice_b1,&
                                   kslice_b2,kslice_fermi_level,&
                                   found_kslice_fermi_level,&
                                   kslice_fermi_lines_colour,recip_lattice,&
                                   nfermi,fermi_energy_list,berry_curv_unit
    use w90_get_oper, only       : get_HH_R,HH_R,get_AA_R,get_BB_R,get_CC_R,&
                                   get_SS_R
    use w90_wan_ham, only        : get_eig_deleig
    use w90_spin, only           : get_spin_nk
    use w90_berry, only          : get_imf_k_list,get_imfgh_k_list
    use w90_constants, only      : bohr

    integer           :: loop_xy,loop_x,loop_y,n,n1,n2,n3,i
    integer           :: zdataunit,coorddataunit,& 
                         bandsunit,scriptunit,dataunit
    real(kind=dp)     :: bvec(3,3),yvec(3),zvec(3),b1mod,b2mod,ymod,cosb1b2,&
                         areab1b2,cosyb2,kpt(3),kpt_x,kpt_y,k1,k2,&
                         imf_k_list(3,3,nfermi),img_k_list(3,3,nfermi),&
                         imh_k_list(3,3,nfermi),Morb_k(3,3),curv(3),morb(3),&
                         spn_k(num_wann),del_eig(num_wann,3),Delta_k,Delta_E,&
                         zhat(3),vdum(3),db1,db2
    logical           :: plot_fermi_lines,plot_curv,plot_morb,fermi_lines_color,&
                         heatmap
    character(len=40) :: filename,square

    integer,          allocatable :: bnddataunit(:)
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    real(kind=dp),    allocatable :: eig(:)

    ! Everything is done on the root node.  However, we still have to
    ! read and distribute the data if we are in parallel, so calls to
    ! get_oper are done on all nodes
   
    plot_fermi_lines=.false.
    if(index(kslice_task,'fermi_lines')>0) plot_fermi_lines=.true.
    plot_curv=.false.
    if(index(kslice_task,'curv')>0) plot_curv=.true.
    plot_morb=.false.
    if(index(kslice_task,'morb')>0) plot_morb=.true.
    fermi_lines_color=.false.
    if(kslice_fermi_lines_colour/='none') fermi_lines_color=.true.
    heatmap=.false.
    if(plot_curv .or. plot_morb) heatmap=.true.

    if(plot_fermi_lines .and. fermi_lines_color .and. heatmap)&
         call io_error('Error: spin-colored Fermi lines not allowed in '&
         //'curv/morb heatmap plots') 

    call get_HH_R
    if(plot_curv.or.plot_morb) call get_AA_R
    if(plot_morb) then
       call get_BB_R
       call get_CC_R
    endif
    if(plot_fermi_lines .and. kslice_fermi_lines_colour=='spin') call get_SS_R

    if(on_root) then

       ! Set Cartesian components of the vectors (b1,b2) spanning the slice, 
       !
       bvec(1,:)=matmul(kslice_b1(:),recip_lattice(:,:))
       bvec(2,:)=matmul(kslice_b2(:),recip_lattice(:,:))
       ! z_vec (orthogonal to the slice)
       zvec(1)=bvec(1,2)*bvec(2,3)-bvec(1,3)*bvec(2,2)
       zvec(2)=bvec(1,3)*bvec(2,1)-bvec(1,1)*bvec(2,3)
       zvec(3)=bvec(1,1)*bvec(2,2)-bvec(1,2)*bvec(2,1)
       ! y_vec (orthogonal to b1=x_vec)
       yvec(1)=zvec(2)*bvec(1,3)-zvec(3)*bvec(1,2)
       yvec(2)=zvec(3)*bvec(1,1)-zvec(1)*bvec(1,3)
       yvec(3)=zvec(1)*bvec(1,2)-zvec(2)*bvec(1,1)
       ! Area (modulus b1 x b2 = z_vec)
       areab1b2=sqrt(zvec(1)**2+zvec(2)**2+zvec(3)**2)
       ! Moduli b1,b2,y_vec
       b1mod=sqrt(bvec(1,1)**2+bvec(1,2)**2+bvec(1,3)**2)
       b2mod=sqrt(bvec(2,1)**2+bvec(2,2)**2+bvec(2,3)**2)
       ymod=sqrt(yvec(1)**2+yvec(2)**2+yvec(3)**2)
       ! Cosine of the angle between y_vec and b2
       cosyb2=yvec(1)*bvec(2,1)+yvec(2)*bvec(2,2)+yvec(3)*bvec(2,3)
       cosyb2=cosyb2/(ymod*b2mod)
       ! Cosine of the angle between b1=x_vec and b2
       cosb1b2=bvec(1,1)*bvec(2,1)+bvec(1,2)*bvec(2,2)+bvec(1,3)*bvec(2,3)
       cosb1b2=cosb1b2/(b1mod*b2mod)       
!       if (abs(cosb1b2)<eps8 .and. b1mod==b2mod) then
       if (abs(cosb1b2)<eps8 .and. abs(b1mod-b2mod)<eps8) then
         square='True'
       else
         square='False'
       end if  

       write(stdout,'(/,/,1x,a)')&
            'Properties calculated in module  k s l i c e'
       write(stdout,'(1x,a)')&
            '--------------------------------------------'

       if(plot_fermi_lines) then
          if(.not.found_kslice_fermi_level) call io_error&
               ('Error: must specify either fermi_energy or'&
               //' kslice_fermi_level when kslice_task = fermi_lines')
          select case(kslice_fermi_lines_colour)
          case("none")
             write(stdout,'(/,3x,a)') '* Fermi lines'
          case("spin")
             write(stdout,'(/,3x,a)') '* Fermi lines coloured by spin'
          end select
          write(stdout,'(/,7x,a,f10.4,1x,a)')&
               '(Fermi level: ',kslice_fermi_level,'eV)'
       endif
       if(plot_curv) then
          if(berry_curv_unit=='ang2') then
             write(stdout,'(/,3x,a)') '* Negative Berry curvature in Ang^2'
          elseif(berry_curv_unit=='bohr2') then
             write(stdout,'(/,3x,a)') '* Negative Berry curvature in Bohr^2'
          endif
          if(nfermi/=1) call io_error('Need to specify one value of '&
               //'the fermi energy when kslice_task=curv')
       elseif(plot_morb) then
          write(stdout,'(/,3x,a)')&
               '* Orbital magnetization k-space integrand in eV.Ang^2'
          if(nfermi/=1) call io_error('Need to specify one value of '&
               //'the fermi energy when kslice_task=morb')
       endif

       write(stdout,'(/,/,1x,a)') 'Output files:' 

       if(.not.fermi_lines_color) then
          coorddataunit=io_file_unit() 
          filename=trim(seedname)//'-kslice-coord.dat'
          write(stdout,'(/,3x,a)') filename
          open(coorddataunit,file=filename,form='formatted')
       endif

       if(plot_fermi_lines) then
          allocate(HH(num_wann,num_wann))
          allocate(UU(num_wann,num_wann))
          allocate(eig(num_wann))
          if(fermi_lines_color) then
             allocate(delHH(num_wann,num_wann,3))
             dataunit=io_file_unit()
             filename=trim(seedname)//'-kslice-fermi-spn.dat'
             write(stdout,'(/,3x,a)') filename
             open(dataunit,file=filename,form='formatted')
          else
             bandsunit=io_file_unit()
             filename=trim(seedname)//'-kslice-bands.dat'
             write(stdout,'(/,3x,a)') filename
             open(bandsunit,file=filename,form='formatted')
             if(.not.heatmap) then
                allocate(bnddataunit(num_wann))
                do n=1,num_wann
                   n1=n/100
                   n2=(n-n1*100)/10
                   n3=n-n1*100-n2*10
                   bnddataunit(n)=io_file_unit()
                   filename=trim(seedname)//'-bnd_'&
                        //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
                   write(stdout,'(/,3x,a)') filename
                   open(bnddataunit(n),file=filename,form='formatted')
                enddo
             endif
          endif
       endif

       if(plot_curv) then
          zdataunit=io_file_unit()
          filename=trim(seedname)//'-kslice-curv.dat'
          write(stdout,'(/,3x,a)') filename
          open(zdataunit,file=filename,form='formatted')
       elseif(plot_morb) then
          zdataunit=io_file_unit()
          filename=trim(seedname)//'-kslice-morb.dat'
          write(stdout,'(/,3x,a)') filename
          open(zdataunit,file=filename,form='formatted')
       end if
     
       db1=1.0_dp/real(kslice_2dkmesh(1),dp)
       db2=1.0_dp/real(kslice_2dkmesh(2),dp)

       ! Loop over uniform mesh of k-points on the slice
       !
       do loop_xy=0,product(kslice_2dkmesh)-1
          loop_x=loop_xy/kslice_2dkmesh(2)
          loop_y=loop_xy-loop_x*kslice_2dkmesh(2)          
          ! k1 and k2 are the coefficients of the k-point in the basis
          ! (kslice_b1,kslice_b2)
          k1=loop_x*db1
          k2=loop_y*db2             
          kpt=kslice_corner+k1*kslice_b1+k2*kslice_b2
          ! Convert to (kpt_x,kpt_y), the 2D Cartesian coordinates
          ! with x along x_vec=b1 and y along y_vec
          kpt_x=k1*b1mod+k2*b2mod*cosb1b2
          kpt_y=k2*b2mod*cosyb2
          if(.not.fermi_lines_color) write(coorddataunit,'(2E16.8)') kpt_x,kpt_y

          if(plot_fermi_lines) then
             if(fermi_lines_color) then
                call get_spin_nk(kpt,spn_k)
                do n=1,num_wann
                   if(spn_k(n)>1.0_dp-eps8) then
                      spn_k(n)=1.0_dp-eps8
                   elseif(spn_k(n)<-1.0_dp+eps8) then
                      spn_k(n)=-1.0_dp+eps8
                   endif
                enddo
                call get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
                Delta_k=max(b1mod*db1,b2mod*db2)
             else
                call fourier_R_to_k(kpt,HH_R,HH,0)
                call utility_diagonalize(HH,num_wann,eig,UU)
             endif
             do n=1,num_wann
                if(.not.fermi_lines_color) then
                   ! For python
                   write(bandsunit,'(E16.8)') eig(n)
                   ! For gnuplot, using 'grid data' format
                    if(.not.heatmap) then
                       write(bnddataunit(n),'(3E16.8)') kpt_x,kpt_y,eig(n)
                       if(loop_y==kslice_2dkmesh(2)-1 .and. &
                            loop_x/=kslice_2dkmesh(1)-1) write (bnddataunit(n),*) ' '
                    endif
                elseif(kslice_fermi_lines_colour=='spin') then
                   ! vdum = dE/dk projected on the k-slice
                   zhat=zvec/sqrt(dot_product(zvec,zvec))
                   vdum(:)=del_eig(n,:)-dot_product(del_eig(n,:),zhat)*zhat(:)
                   Delta_E=sqrt(dot_product(vdum,vdum))*Delta_k
!                   Delta_E=Delta_E*sqrt(2.0_dp) ! optimize this factor
                   if(abs(eig(n)-kslice_fermi_level)<Delta_E)&
                        write(dataunit,'(3E16.8)') kpt_x,kpt_y,spn_k(n)
                endif
             enddo
          endif

          if(plot_curv) then
             call get_imf_k_list(kpt,imf_k_list)
             curv(1)=sum(imf_k_list(:,1,1))
             curv(2)=sum(imf_k_list(:,2,1))
             curv(3)=sum(imf_k_list(:,3,1))
             if(berry_curv_unit=='bohr2') curv=curv/bohr**2   
             ! Print the negative Berry curvature 
             write(zdataunit,'(3E16.8)') -curv(:)
          end if

          if(plot_morb) then
             call get_imfgh_k_list(kpt,imf_k_list,img_k_list,imh_k_list)
             Morb_k=img_k_list(:,:,1)+imh_k_list(:,:,1)&
                   -2.0_dp*fermi_energy_list(1)*imf_k_list(:,:,1)
             Morb_k=-Morb_k/2.0_dp ! differs by -1/2 from Eq.97 LVTS12
             morb(1)=sum(Morb_k(:,1))
             morb(2)=sum(Morb_k(:,2))
             morb(3)=sum(Morb_k(:,3))
             write(zdataunit,'(3E16.8)') morb(:)
          end if

       end do !loop_xy
       
       if(.not.fermi_lines_color) then
          write(coorddataunit,*) ' '
          close(coorddataunit)
       endif
       
       if(heatmap) then
          write(zdataunit,*) ' '
          close(zdataunit)
       endif
       if(plot_fermi_lines) then
          if(fermi_lines_color) then
             close(dataunit)
          else
             write(bandsunit,*) ' '
             close(bandsunit)
             if(.not.heatmap) then
                do n=1,num_wann
                   write(bnddataunit(n),*) ' '
                   close(bnddataunit(n))
                enddo
             endif
          endif
       endif

       if(plot_fermi_lines .and. .not.fermi_lines_color .and. .not.heatmap) then
          !
          ! gnuplot script for black Fermi lines
          !
          scriptunit=io_file_unit()
          filename=trim(seedname)//'-kslice-fermi_lines.gnu'
          write(stdout,'(/,3x,a)') filename
          open(scriptunit,file=filename,form='formatted')
          write(scriptunit,'(a)') 'unset surface'
          write(scriptunit,'(a)') 'set contour'
          write(scriptunit,'(a)') 'set view map'
          write(scriptunit,'(a,f9.5)') 'set cntrparam levels discrete ',&
               kslice_fermi_level
          write(scriptunit,'(a)') 'set cntrparam bspline'
          do n=1,num_wann
             n1=n/100
             n2=(n-n1*100)/10
             n3=n-n1*100-n2*10
             write(scriptunit,'(a)') 'set table "bnd_'&
                  //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat"'
             write(scriptunit,'(a)') 'splot "'//trim(seedname)//'-bnd_'&
                  //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat"'
             write(scriptunit,'(a)') 'unset table'
          enddo
          write(scriptunit,'(a)')&
               '#Uncomment next two lines to create postscript'
          write(scriptunit,'(a)') '#set term post eps enh'
          write(scriptunit,'(a)')&
               '#set output "'//trim(seedname)//'-kslice-fermi_lines.eps"'
          write(scriptunit,'(a)') 'set size ratio -1'
          write(scriptunit,'(a)') 'unset tics'
          write(scriptunit,'(a)') 'unset key'
          write(scriptunit,'(a)')&
               '#For postscript try changing lw 1 --> lw 2 in the next line'
          write(scriptunit,'(a)') 'set style line 1 lt 1 lw 1'
          if(num_wann==1) then
             write(scriptunit,'(a)')&
                  'plot "bnd_001.dat" using 1:2 w lines ls 1'
          else
             write(scriptunit,'(a)')&
                  'plot "bnd_001.dat" using 1:2 w lines ls 1,'&
                  //achar(92)
          endif
          do n=2,num_wann-1
             n1=n/100
             n2=(n-n1*100)/10
             n3=n-n1*100-n2*10
             write(scriptunit,'(a)') '     "bnd_'&
                  //achar(48+n1)//achar(48+n2)//achar(48+n3)&
                  //'.dat" using 1:2 w lines ls 1,'//achar(92)
          enddo
          n=num_wann
          n1=n/100
          n2=(n-n1*100)/10
          n3=n-n1*100-n2*10
          write(scriptunit,'(a)') '     "bnd_'&
               //achar(48+n1)//achar(48+n2)//achar(48+n3)&
               //'.dat" using 1:2 w lines ls 1'
          close(scriptunit)
          !
          ! Python script for black Fermi lines
          !  
          scriptunit=io_file_unit()
          filename=trim(seedname)//'-kslice-fermi_lines.py'
          write(stdout,'(/,3x,a)') filename
          open(scriptunit,file=filename,form='formatted')      
          write(scriptunit,'(a)') 'import pylab as pl' 
          write(scriptunit,'(a)') 'import numpy as np'
          write(scriptunit,'(a)') 'import matplotlib.mlab as ml'
          write(scriptunit,'(a)') 'from collections import OrderedDict'
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)') "points = np.loadtxt('"//trim(seedname)//&
                                       "-kslice-coord.dat')"
          write(scriptunit,'(a)') 'points_x=points[:,0]'
          write(scriptunit,'(a)') 'points_y=points[:,1]'
          write(scriptunit,'(a)') 'num_pt=len(points)'             
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a,f12.6)') 'area=', areab1b2
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)') 'square= '//square
          write(scriptunit,'(a)') ' '

          write(scriptunit,'(a)') 'if square:'
          write(scriptunit,'(a)')&
               '  x_coord=list(OrderedDict.fromkeys(points_x))'
          write(scriptunit,'(a)')&
               '  y_coord=list(OrderedDict.fromkeys(points_y))'
          write(scriptunit,'(a)') '  dimx=len(x_coord)'
          write(scriptunit,'(a)') '  dimy=len(y_coord)'
          write(scriptunit,'(a)') 'else:'
          write(scriptunit,'(a)') '  xmin=np.min(points_x)'
          write(scriptunit,'(a)') '  ymin=np.min(points_y)'
          write(scriptunit,'(a)') '  xmax=np.max(points_x)'
          write(scriptunit,'(a)') '  ymax=np.max(points_y)'  
          write(scriptunit,'(a)')&
               '  a=np.max(np.array([xmax-xmin,ymax-ymin]))'
          write(scriptunit,'(a)')&
               '  num_int=int(round(np.sqrt(num_pt*a**2/area)))'
          write(scriptunit,'(a)') '  xint = np.linspace(xmin,xmin+a,num_int)'
          write(scriptunit,'(a)') '  yint = np.linspace(ymin,ymin+a,num_int)'
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)')&
              '# Energy level for isocontours (typically the Fermi level)'
          write(scriptunit,'(a,f12.6)') 'ef=',kslice_fermi_level
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)')&
               "bands=np.loadtxt('"//trim(seedname)//"-kslice-bands.dat')"
          write(scriptunit,'(a)') 'numbands=bands.size/num_pt'
          write(scriptunit,'(a)') 'if square:'
          write(scriptunit,'(a)')&
               '  bbands=bands.reshape((dimx,dimy,numbands))'
          write(scriptunit,'(a)') '  for i in range(numbands):'
          write(scriptunit,'(a)') '    pl.contour(x_coord,'&
               //'y_coord,bbands[:,:,i],[ef],colors="black")'
          write(scriptunit,'(a)') 'else:'
          write(scriptunit,'(a)') '  bbands=bands.reshape((num_pt,numbands))'
          write(scriptunit,'(a)') '  bandint=[]'
          write(scriptunit,'(a)') '  for i in range(numbands):'
          write(scriptunit,'(a)') '    bandint.append(ml.griddata'&
               //'(points_x,points_y, bbands[:,i], xint, yint))'
          write(scriptunit,'(a)') '    pl.contour(xint,yint,bandint[i],'&
               //'[ef],colors="black")'                             
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)') '# Remove the axes'
          write(scriptunit,'(a)') 'ax = pl.gca()'
          write(scriptunit,'(a)') 'ax.xaxis.set_visible(False)'
          write(scriptunit,'(a)') 'ax.yaxis.set_visible(False)'
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)') "pl.axes().set_aspect('equal')"
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)') "outfile = '"//trim(seedname)//&
               "-fermi_lines.pdf'"
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)') ' '
          write(scriptunit,'(a)') 'pl.savefig(outfile)'
          write(scriptunit,'(a)') 'pl.show()'
          close(scriptunit)
       endif !plot_fermi_lines .and. kslice_fermi_lines_colour=='none'
             !.and. .not.heatmap

       if(plot_fermi_lines .and. fermi_lines_color .and. .not.heatmap) then
          !
          ! gnuplot script for spin-colored Fermi lines
          !
          scriptunit=io_file_unit()
          filename=trim(seedname)//'-kslice-fermi_lines.gnu'
          write(stdout,'(/,3x,a)') filename
          open(scriptunit,file=filename,form='formatted')
          write(scriptunit,'(a)') 'unset key'
          write(scriptunit,'(a)') 'unset tics'
          write(scriptunit,'(a)') 'set cbtics'
          write(scriptunit,'(a)')&
               'set palette defined (-1 "blue", 0 "green", 1 "red")'
          write(scriptunit,'(a)') 'set pm3d map'
          write(scriptunit,'(a)') 'set zrange [-1:1]'
          write(scriptunit,'(a)') 'set size ratio -1'
          write(scriptunit,'(a)')&
               '#Uncomment next two lines to create postscript'
           write(scriptunit,'(a)') '#set term post eps enh'
          write(scriptunit,'(a)') '#set output "'&
               //trim(seedname)//'-kslice-fermi_lines.eps"'
          write(scriptunit,'(a)') 'splot "'&
               //trim(seedname)//'-kslice-fermi-spn.dat" with dots palette'
          !
          ! python script for spin-colored Fermi lines
          !
          scriptunit=io_file_unit()
          filename=trim(seedname)//'-kslice-fermi_lines.py'
          write(stdout,'(/,3x,a)') filename
          open(scriptunit,file=filename,form='formatted')
          write(scriptunit,'(a)') 'import pylab as pl' 
          write(scriptunit,'(a)') 'import numpy as np'
          write(scriptunit,'(a)') "data = np.loadtxt('"//trim(seedname)//&
               "-kslice-fermi-spn.dat')"
          write(scriptunit,'(a)') 'x=data[:,0]'
          write(scriptunit,'(a)') 'y=data[:,1]'
          write(scriptunit,'(a)') 'z=data[:,2]'
          write(scriptunit,'(a)')&
               "pl.scatter(x,y,c=z,marker='+',s=2,cmap=pl.cm.jet)"
          write(scriptunit,'(a,F12.6,a)')&
               "pl.plot([0,",kpt_x,"],[0,0],color='black',linestyle='-',"&
               //"linewidth=0.5)"
          write(scriptunit,'(a,F12.6,a,F12.6,a,F12.6,a)')&
               "pl.plot([",kpt_x,",",kpt_x,"],[0,",kpt_y,"],color='black',"&
               //"linestyle='-',linewidth=0.5)"
          write(scriptunit,'(a,F12.6,a,F12.6,a,F12.6,a)')&
               "pl.plot([0,",kpt_x,"],[",kpt_y,",",kpt_y,&
               "],color='black',linestyle='-',linewidth=0.5)"
          write(scriptunit,'(a,F12.6,a)') "pl.plot([0,0],[0,",kpt_y,&
               "],color='black',linestyle='-',linewidth=0.5)"
          write(scriptunit,'(a,F12.6,a)') 'pl.xlim([0,',kpt_x,'])'
          write(scriptunit,'(a,F12.6,a)') 'pl.ylim([0,',kpt_y,'])'
          write(scriptunit,'(a)') 'cbar=pl.colorbar()'
          write(scriptunit,'(a)') 'ax = pl.gca()'
          write(scriptunit,'(a)') 'ax.xaxis.set_visible(False)'
          write(scriptunit,'(a)') 'ax.yaxis.set_visible(False)'
          write(scriptunit,'(a)') "pl.savefig('"//trim(seedname)//&
               "-kslice-fermi_lines.pdf')"
          write(scriptunit,'(a)') 'pl.show()'
          close(scriptunit)
       endif ! plot_fermi_lines .and. fermi_lines_color .and. .not.heatmap

       if(heatmap) then
          !
          ! python script for curvature/Morb heatmaps [+ black Fermi lines]
          !
          do i=1,3

             scriptunit=io_file_unit()
             if(plot_curv .and. .not.plot_fermi_lines) then
                filename=trim(seedname)//'-kslice-curv_'//achar(119+i)//'.py'
                write(stdout,'(/,3x,a)') filename
                open(scriptunit,file=filename,form='formatted')
             elseif(plot_curv .and. plot_fermi_lines) then
                filename=trim(seedname)//'-kslice-curv_'//achar(119+i)//&
                     '+fermi_lines.py'
                write(stdout,'(/,3x,a)') filename
                open(scriptunit,file=filename,form='formatted')
             elseif(plot_morb .and. .not.plot_fermi_lines) then
                filename=trim(seedname)//'-kslice-morb_'//achar(119+i)//'.py'
                write(stdout,'(/,3x,a)') filename
                open(scriptunit,file=filename,form='formatted')
             elseif(plot_morb .and. plot_fermi_lines) then
                filename=trim(seedname)//'-kslice-morb_'//achar(119+i)//&
                     '+fermi_lines.py'
                write(stdout,'(/,3x,a)') filename
                open(scriptunit,file=filename,form='formatted')
             endif
             write(scriptunit,'(a)') 'import pylab as pl'
             write(scriptunit,'(a)') 'import numpy as np'
             write(scriptunit,'(a)') 'import matplotlib.mlab as ml'
             write(scriptunit,'(a)') 'from collections import OrderedDict'
             write(scriptunit,'(a)') ' '
             write(scriptunit,'(a)') "points = np.loadtxt('"//trim(seedname)//&
                                          "-kslice-coord.dat')"
             write(scriptunit,'(a)') 'points_x=points[:,0]'
             write(scriptunit,'(a)') 'points_y=points[:,1]'
             write(scriptunit,'(a)') 'num_pt=len(points)'             
             write(scriptunit,'(a)') ' '
             write(scriptunit,'(a,f12.6)') 'area=', areab1b2
             write(scriptunit,'(a)') ' '
             write(scriptunit,'(a)') 'square= '//square
             write(scriptunit,'(a)') ' '
             
             write(scriptunit,'(a)') 'if square:'
             write(scriptunit,'(a)')&
                  '  x_coord=list(OrderedDict.fromkeys(points_x))'
             write(scriptunit,'(a)')&
                  '  y_coord=list(OrderedDict.fromkeys(points_y))'
             write(scriptunit,'(a)') '  dimx=len(x_coord)'
             write(scriptunit,'(a)') '  dimy=len(y_coord)'
             write(scriptunit,'(a)') 'else:'
             write(scriptunit,'(a)') '  xmin=np.min(points_x)'
             write(scriptunit,'(a)') '  ymin=np.min(points_y)'
             write(scriptunit,'(a)') '  xmax=np.max(points_x)'
             write(scriptunit,'(a)') '  ymax=np.max(points_y)'  
             write(scriptunit,'(a)')&
                  '  a=np.max(np.array([xmax-xmin,ymax-ymin]))'
             write(scriptunit,'(a)')&
                  '  num_int=int(round(np.sqrt(num_pt*a**2/area)))'
             write(scriptunit,'(a)') '  xint = np.linspace(xmin,xmin+a,num_int)'
             write(scriptunit,'(a)')&
                  '  yint = np.linspace(ymin,ymin+a,num_int)'
             write(scriptunit,'(a)') ' '
             
             if(plot_fermi_lines) then
                write(scriptunit,'(a)')&
                    '# Energy level for isocontours (typically the Fermi level)'
                write(scriptunit,'(a,f12.6)') 'ef=',kslice_fermi_level
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)')&
                     "bands=np.loadtxt('"//trim(seedname)//"-kslice-bands.dat')"
                write(scriptunit,'(a)') 'numbands=bands.size/num_pt'
                write(scriptunit,'(a)') 'if square:'
                write(scriptunit,'(a)')&
                     '  bbands=bands.reshape((dimx,dimy,numbands))'
                write(scriptunit,'(a)') '  for i in range(numbands):'
                write(scriptunit,'(a)') '    pl.contour(x_coord,y_coord,'&
                     //'bbands[:,:,i],[ef],colors="black")'
                write(scriptunit,'(a)') 'else:'
                write(scriptunit,'(a)') '  bbands=bands.reshape((num_pt,'&
                     //'numbands))'
                write(scriptunit,'(a)') '  bandint=[]'
                write(scriptunit,'(a)') '  for i in range(numbands):'
                write(scriptunit,'(a)') '    bandint.append(ml.griddata'&
                     //'(points_x,points_y, bbands[:,i], xint, yint))'
                write(scriptunit,'(a)') '    pl.contour(xint,yint,'&
                     //'bandint[i],[ef],colors="black")'     
             endif
             
             if(plot_curv) then
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)') "outfile = '"//trim(seedname)//&
               "-kslice-curv_"//achar(119+i)//".pdf'"
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)')&
                     "val = np.loadtxt('"//trim(seedname)//&
                     "-kslice-curv.dat', usecols=("//achar(47+i)//",))"
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)')&
                     'val_log=np.array([np.log10(abs(elem))*np.sign(elem) &
                 &if abs(elem)>10 else elem/10.0 for elem in val])'
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)') 'if square: '
                write(scriptunit,'(a)')&
                     '  vval=val_log.reshape(dimx,dimy).transpose()'
                write(scriptunit,'(a)') '  mn=int(np.floor(vval.min()))'
                write(scriptunit,'(a)') '  mx=int(np.ceil(vval.max()))' 
                write(scriptunit,'(a)') '  ticks=range(mn,mx+1)'
                write(scriptunit,'(a)') "  pl.contourf(x_coord,y_coord,"&
                     //"vval,ticks,origin='lower')"
                write(scriptunit,'(a)') '  #pl.imshow(vval,origin="lower",'&
                     //'extent=(min(x_coord),max(x_coord),min(y_coord),'&
                     //'max(y_coord)))'
                write(scriptunit,'(a)') 'else: '
                write(scriptunit,'(a)') '  valint = ml.griddata(points_x,'&
                     //'points_y, val_log, xint, yint)'  
                write(scriptunit,'(a)') '  mn=int(np.floor(valint.min()))'
                write(scriptunit,'(a)') '  mx=int(np.ceil(valint.max()))' 
                write(scriptunit,'(a)') '  ticks=range(mn,mx+1)'
                write(scriptunit,'(a)') '  pl.contourf(xint,yint,valint,ticks)'
                write(scriptunit,'(a)') '  #pl.imshow(valint,origin="lower",'&
                     //'extent=(min(xint),max(xint),min(yint),max(yint)))'
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)') 'ticklabels=[]'
                write(scriptunit,'(a)') 'for n in ticks:'
                write(scriptunit,'(a)') ' if n<0: '
                write(scriptunit,'(a)')&
                     "  ticklabels.append('-$10^{%d}$' % abs(n))"
                write(scriptunit,'(a)') ' elif n==0:'
                write(scriptunit,'(a)') "  ticklabels.append(' $%d$' %  n)" 
                write(scriptunit,'(a)') ' else:'
                write(scriptunit,'(a)') "  ticklabels.append(' $10^{%d}$' % n)" 
                write(scriptunit,'(a)') ' '           
                write(scriptunit,'(a)') 'cbar=pl.colorbar()'              
                write(scriptunit,'(a)') 'cbar.set_ticks(ticks)'
                write(scriptunit,'(a)') 'cbar.set_ticklabels(ticklabels)'
         
             elseif(plot_morb) then
               
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)') "outfile = '"//trim(seedname)//&
               "-kslice-morb_"//achar(119+i)//".pdf'"
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)')&
                     "val = np.loadtxt('"//trim(seedname)//&
                     "-kslice-morb.dat', usecols=("//achar(47+i)//",))"
               write(scriptunit,'(a)') ' '
               write(scriptunit,'(a)') 'if square: '
               write(scriptunit,'(a)')&
                    '  vval=val.reshape(dimx,dimy).transpose()'
               write(scriptunit,'(a)') '  pl.imshow(vval,origin="lower",'&
                    //'extent=(min(x_coord),max(x_coord),min(y_coord),'&
                    //'max(y_coord)))'
               write(scriptunit,'(a)') 'else: '
               write(scriptunit,'(a)') '  valint = ml.griddata(points_x,'&
                    //'points_y, val, xint, yint)' 
               write(scriptunit,'(a)') '  pl.imshow(valint,origin="lower",'&
                    //'extent=(min(xint),max(xint),min(yint),max(yint)))'
               write(scriptunit,'(a)') 'pl.colorbar()'

            endif
                       
            write(scriptunit,'(a)') ' '
            write(scriptunit,'(a)') 'ax = pl.gca()'
            write(scriptunit,'(a)') 'ax.xaxis.set_visible(False)'
            write(scriptunit,'(a)') 'ax.yaxis.set_visible(False)'
            write(scriptunit,'(a)') ' '
            write(scriptunit,'(a)') 'pl.savefig(outfile)'
            write(scriptunit,'(a)') 'pl.show()'
            
            close(scriptunit)
          enddo
          !
       endif !heatmap

       write(stdout,*) ' '

    end if ! on_root
 
end subroutine k_slice

end module w90_kslice
