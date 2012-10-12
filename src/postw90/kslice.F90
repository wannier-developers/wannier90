  !-*- mode: F90; mode: font-lock -*-!

module w90_kslice

  ! Makes a heatmap plot on a slice in k-space of:
  ! 
  !  - The negative Berry curvature summed over occupied bands 
  !
  !  - The integrand of the k-space orbital magnetization formula
  !
  ! Calculates the intersections of constant-energy isosurfaces with the slice
  !
  ! The slice is defined by three input variables, all in reciprocal
  ! lattice coordinates:
  !
  !    slice_corner(1:3) is the lower left corner 
  !    slice_b1(1:3) and slice_b2(1:3) are the vectors subtending the slice

  !---------------------------------------------------------------------
  ! TO DO:  
  !
  !      * Add the possibility to color the energy contours by the spin
  !
  !      * Add python and octave script for energy_cntr (in addition to
  !        gnuplot), allow to choose between them (uncomment 
  !        kslice_plot_format in parameters.F90)
  !---------------------------------------------------------------------

  implicit none

  public

  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine k_slice

    use w90_comms
    use w90_constants,  only     : dp,twopi
    use w90_io,         only     : io_error,io_file_unit,seedname,&
                                   io_time,io_stopwatch,stdout
    use w90_utility, only        : utility_diagonalize
    use w90_postw90_common, only : fourier_R_to_k
    use w90_parameters, only     : num_wann,kslice,kslice_task,&
                                   kslice_interp_mesh,kslice_corner,kslice_b1,&
                                   kslice_b2,kslice_cntr_energy,&
                                   found_kslice_cntr_energy,recip_lattice,&
                                   found_fermi_energy,fermi_energy
    use w90_get_oper, only       : get_HH_R,HH_R,get_AA_R,get_BB_R,get_CC_R
    use w90_berry, only          : get_imf_k,get_img_k,get_imh_k
    use w90_utility, only        : utility_recip_lattice
    use w90_constants, only      : bohr,ev_au

    integer           :: loop_tot,loop_x,loop_y,n,n1,n2,n3,i
    integer           :: xdataunit,ydataunit,zdataunit,bandsunit,scriptunit
    real(kind=dp)     :: avec(3,3),bvec(3,3),recip_vol,b1mod,b2mod,cosb1b2,&
                         a2mod,kpt(3),kpt_x,kpt_y,k1,k2,k_cart(3),&
                         imf_k(3,3),img_k(3,3),imh_k(3,3),&
                         curv_au(3),morb_au(3),Morb_k(3,3)
    logical           :: plot_curv_heatmap,plot_morb_heatmap,plot_energy_cntr
    character(len=20) :: filename

    integer, allocatable          :: xyzdataunit(:)
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    real(kind=dp),    allocatable :: eig(:)

    ! Everything is done on the root node (not worthwhile parallelizing). 
    ! However, we still have to read and distribute the data if we 
    ! are in parallel. So calls to get_oper are done on all nodes at the moment

    plot_curv_heatmap=.false.
    if(index(kslice_task,'curv_heatmap')>0) then
       plot_curv_heatmap=.true.
    end if
    plot_morb_heatmap=.false.
    if(index(kslice_task,'morb_heatmap')>0) then
       plot_morb_heatmap=.true.
    end if
    plot_energy_cntr=.false.
    if(index(kslice_task,'energy_cntr')>0) then
       plot_energy_cntr=.true.
    end if
    ! Set up the needed Wannier matrix elements
    call get_HH_R
    if(plot_curv_heatmap.or.plot_morb_heatmap) call get_AA_R
    if(plot_morb_heatmap) then
       call get_BB_R
       call get_CC_R
    endif

    if(on_root) then

       ! Set Cartesian components of the vectors (b_1,b_2) spanning the slice, 
       ! and their reciprocals. Store as *rows* in bvec and avec
       !
       bvec(1,:)=matmul(kslice_b1(:),recip_lattice(:,:))
       bvec(2,:)=matmul(kslice_b2(:),recip_lattice(:,:))
       ! Need also b_3 = b_1 \times b_2
       bvec(3,1)=bvec(1,2)*bvec(2,3)-bvec(1,3)*bvec(2,2)
       bvec(3,2)=bvec(1,3)*bvec(2,1)-bvec(1,1)*bvec(2,3)
       bvec(3,3)=bvec(1,1)*bvec(2,2)-bvec(1,2)*bvec(2,1)
       ! Cosine of the angle between b_1 and b_2
       b1mod=sqrt(bvec(1,1)**2+bvec(1,2)**2+bvec(1,3)**2)
       b2mod=sqrt(bvec(2,1)**2+bvec(2,2)**2+bvec(2,3)**2)
       cosb1b2=bvec(1,1)*bvec(2,1)+bvec(1,2)*bvec(2,2)+bvec(1,3)*bvec(2,3)
       cosb1b2=cosb1b2/(b1mod*b2mod)
      ! Find the reciprocal of bvec, store in avec
       call utility_recip_lattice(bvec,avec,recip_vol)
       a2mod=sqrt(avec(2,1)**2+avec(2,2)**2+avec(2,3)**2)

       write(stdout,'(/,/,1x,a)')&
            'Properties calculated in module  k s l i c e'
       write(stdout,'(1x,a)')&
            '--------------------------------------------'

       if(plot_energy_cntr) then
          write(stdout,'(/,3x,a)') '* Constant energy contours'
          if(.not.found_kslice_cntr_energy) call io_error&
               ('Error: must specify either fermi_energy or kslice_cntr_energy when kslice_task = energy_cntr')
          write(stdout,'(/,7x,a,f10.4,1x,a)')&
               '(Energy isocontour level: ',kslice_cntr_energy,'eV)'
       endif
       if(plot_curv_heatmap) then
          write(stdout,'(/,3x,a)') '* Negative Berry curvature in a.u. (bohr^2)'
          if(.not.found_fermi_energy) call io_error&
               (&
 'Need to specify either "fermi_energy" or "num_elec_cell" when plot_curv_heatmap=T'&
               )
       elseif(plot_morb_heatmap) then
          write(stdout,'(/,3x,a)') '* Orbital magnetization in a.u. (Ha.bohr^2)'
          if(.not.found_fermi_energy) call io_error&
               (&
'Need to set either "fermi_energy" or "num_elec_cell" when plot_morb_heatmap=T'&
               )
       endif

       write(stdout,'(/,/,1x,a,/)') 'Output files:' 
       !
       ! octave and python+matplotlib: write the x- and y-axis mesh values
       xdataunit=io_file_unit() 
       filename=trim(seedname)//'_slice_x.dat'
       write(stdout,'(/,3x,a)') filename
       open(xdataunit,file=filename,form='formatted')
       ydataunit=io_file_unit() 
       filename=trim(seedname)//'_slice_y.dat'
       write(stdout,'(/,3x,a)') filename
       open(ydataunit,file=filename,form='formatted')

       if(plot_energy_cntr) then
          allocate(HH(num_wann,num_wann))
          allocate(UU(num_wann,num_wann))
          allocate(eig(num_wann))
          bandsunit=io_file_unit()
          filename=trim(seedname)//'-slice_bands.dat'
          write(stdout,'(/,3x,a)') filename
          open(bandsunit,file=filename,form='formatted')
          allocate(xyzdataunit(num_wann))
          do n=1,num_wann
             n1=n/100
             n2=(n-n1*100)/10
             n3=n-n1*100-n2*10
             xyzdataunit(n)=io_file_unit()
             filename=trim(seedname)//'_xyz_'&
                  //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
             write(stdout,'(/,3x,a)') filename
             open(xyzdataunit(n),file=filename,form='formatted')
          enddo
       endif

       if(plot_curv_heatmap) then
          zdataunit=io_file_unit()
          filename=trim(seedname)//'-slice_curv.dat'
          write(stdout,'(/,3x,a)') filename
          open(zdataunit,file=filename,form='formatted')
       elseif(plot_morb_heatmap) then
          zdataunit=io_file_unit()
          filename=trim(seedname)//'-slice_morb.dat'
          write(stdout,'(/,3x,a)') filename
          open(zdataunit,file=filename,form='formatted')
       end if
     
       ! Loop over uniform mesh of k-points on the slice
       !
       do loop_tot=0,product(kslice_interp_mesh)-1
          loop_x=loop_tot/kslice_interp_mesh(2)
          loop_y=loop_tot-loop_x*kslice_interp_mesh(2)
          kpt(1)=kslice_corner(1)&
               +kslice_b1(1)*real(loop_x,dp)/real(kslice_interp_mesh(1),dp)&
               +kslice_b2(1)*real(loop_y,dp)/real(kslice_interp_mesh(2),dp)
          kpt(2)=kslice_corner(2)&
               +kslice_b1(2)*real(loop_x,dp)/real(kslice_interp_mesh(1),dp)&
               +kslice_b2(2)*real(loop_y,dp)/real(kslice_interp_mesh(2),dp)
          kpt(3)=kslice_corner(3)&
               +kslice_b1(3)*real(loop_x,dp)/real(kslice_interp_mesh(1),dp)&
               +kslice_b2(3)*real(loop_y,dp)/real(kslice_interp_mesh(2),dp)

          ! Convert to (x,y) Cartesian coordinates, with slice_b1 along x and
          ! (slice_b1,slice_b2) in the xy-plane
          !
          k_cart(:)=matmul(kpt(:),recip_lattice(:,:))
          ! k1 and k2 are the coefficients of the k-point in the basis
          ! (slice_b1,slice_b2)
          k1=(k_cart(1)*avec(1,1)+k_cart(2)*avec(1,2)+k_cart(3)*avec(1,3))/twopi
          k2=(k_cart(1)*avec(2,1)+k_cart(2)*avec(2,2)+k_cart(3)*avec(2,3))/twopi
          kpt_x=k1*b1mod+k2*b2mod*cosb1b2
          kpt_y=twopi*k2/a2mod

          if(loop_x==0) write(ydataunit,'(e16.8)') kpt_y
          if(loop_y==0) write(xdataunit,'(e16.8)') kpt_x

          if(plot_curv_heatmap) then
             call get_imf_k(kpt,imf_k)
             curv_au(1)=sum(imf_k(:,1))/bohr**2
             curv_au(2)=sum(imf_k(:,2))/bohr**2
             curv_au(3)=sum(imf_k(:,3))/bohr**2
   
             ! Heatmap plot of the negative Berry curvature using a 
             ! "log" scale as in Yao's 2004 PRL
             !
             curv_au=-curv_au
             do i=1,3
                if(curv_au(i)>0) then
                   if(curv_au(i)>10) then
                      curv_au(i)=log10(curv_au(i))
                   else
                      curv_au(i)=0.5_dp ! 1.0 would work as well
                   endif
                else
                   if(curv_au(i)<-10) then
                      curv_au(i)=-log10(-curv_au(i))
                   else
                      curv_au(i)=-0.5_dp ! -0.99 would work, but not -1.0!
                   endif
                endif
             enddo
             write(zdataunit,'(3E16.8)') curv_au(:)
          end if

          if(plot_morb_heatmap) then
             call get_imf_k(kpt,imf_k)
             call get_img_k(kpt,img_k)
             call get_imh_k(kpt,imh_k)
             Morb_k=img_k+imh_k-2.0_dp*fermi_energy*imf_k
             morb_au(1)=sum(Morb_k(:,1))
             morb_au(2)=sum(Morb_k(:,2))
             morb_au(3)=sum(Morb_k(:,3))
             morb_au=morb_au*eV_au/bohr**2
             write(zdataunit,'(3E16.8)') morb_au(:)
          end if

          if(plot_energy_cntr) then
             call fourier_R_to_k(kpt,HH_R,HH,0)
             call utility_diagonalize(HH,num_wann,eig,UU)
             do n=1,num_wann
                ! For python/octave
                write(bandsunit,'(E16.8)') eig(n)
                ! For gnuplot, using 'grid data' format
                write(xyzdataunit(n),'(3E16.8)') kpt_x,kpt_y,eig(n)
                if(loop_y==kslice_interp_mesh(2)-1 .and. &
                   loop_x/=kslice_interp_mesh(1)-1) write (xyzdataunit(n),*) ' '
             enddo
          endif

       end do !loop_tot
       
       write(xdataunit,*) ' '
       close(xdataunit)
       write(ydataunit,*) ' '
       close(ydataunit)
       if(plot_curv_heatmap.or.plot_morb_heatmap) then
          write(zdataunit,*) ' '
          close(zdataunit)
       endif
       if(plot_energy_cntr) then
          write(bandsunit,*) ' '
          close(bandsunit)
          do n=1,num_wann
             write(xyzdataunit(n),*) ' '
             close(xyzdataunit(n))
          enddo
       endif

       if(plot_energy_cntr) then
          scriptunit=io_file_unit()
          open(scriptunit,file=trim(seedname)//'-energy_cntr.gnu',&
               form='formatted')
          write(scriptunit,'(a)') 'unset surface'
          write(scriptunit,'(a)') 'set contour'
          write(scriptunit,'(a)') 'set view map'
          write(scriptunit,'(a,f9.5)') 'set cntrparam levels discrete ',&
               kslice_cntr_energy
          write(scriptunit,'(a)') 'set cntrparam bspline'
          do n=1,num_wann
             n1=n/100
             n2=(n-n1*100)/10
             n3=n-n1*100-n2*10
             write(scriptunit,'(a)') 'set table "xyz_'&
                  //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat"'
             write(scriptunit,'(a)') 'splot "'//trim(seedname)//'_xyz_'&
                  //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat"'
             write(scriptunit,'(a)') 'unset table'
          enddo
          write(scriptunit,'(a)')&
               '#Uncomment next two lines to create postscript'
          write(scriptunit,'(a)') '#set term post eps enh'
          write(scriptunit,'(a)')&
               '#set output "'//trim(seedname)//'_econtour.eps"'
          write(scriptunit,'(a)') 'set size ratio -1'
          write(scriptunit,'(a)') 'unset tics'
          write(scriptunit,'(a)') 'set nokey'
          write(scriptunit,'(a)')&
               '#For postscript try changing lw 1 --> lw 2 in the next line'
          write(scriptunit,'(a)') 'set style line 1 lt 1 lw 1'
          if(num_wann==1) then
             write(scriptunit,'(a)')&
                  'plot "xyz_001.dat" using 1:2 w lines ls 1'
          else
             write(scriptunit,'(a)')&
                  'plot "xyz_001.dat" using 1:2 w lines ls 1,'&
                  //achar(92)
          endif
          do n=2,num_wann-1
             n1=n/100
             n2=(n-n1*100)/10
             n3=n-n1*100-n2*10
             write(scriptunit,'(a)') '     "xyz_'&
                  //achar(48+n1)//achar(48+n2)//achar(48+n3)&
                  //'.dat" using 1:2 w lines ls 1,'//achar(92)
          enddo
          n=num_wann
          n1=n/100
          n2=(n-n1*100)/10
          n3=n-n1*100-n2*10
          write(scriptunit,'(a)') '     "xyz_'&
               //achar(48+n1)//achar(48+n2)//achar(48+n3)&
               //'.dat" using 1:2 w lines ls 1'
       endif !plot_energy_cntr


       if(plot_curv_heatmap .or. plot_morb_heatmap) then
          
          do i=1,3

             scriptunit=io_file_unit()
             if(plot_curv_heatmap) then
                open(scriptunit,file=trim(seedname)//&
                     '-curv_'//achar(119+i)//'-heatmap.py',form='formatted')
             elseif(plot_morb_heatmap) then
                open(scriptunit,file=trim(seedname)//&
                     '-morb_'//achar(119+i)//'-heatmap.py',form='formatted')
             endif
             write(scriptunit,'(a)') 'import pylab as pl'
             write(scriptunit,'(a)') 'import numpy as np'
             write(scriptunit,'(a)') 'import mpl_toolkits'
             write(scriptunit,'(a)') ' '
             write(scriptunit,'(a)') "x = np.loadtxt('"//trim(seedname)//&
                  "_slice_x.dat')"
             write(scriptunit,'(a)') 'dimx=x.size'
             write(scriptunit,'(a)') ' '
             write(scriptunit,'(a)') "y = np.loadtxt('"//trim(seedname)//&
                  "_slice_y.dat')"
             write(scriptunit,'(a)') 'dimy=y.size'
             write(scriptunit,'(a)') ' '
             
             if(plot_energy_cntr) then
                write(scriptunit,'(a)')&
                    '# Energy level for isocontours (typically the Fermi level)'
                write(scriptunit,'(a,f12.6)') 'ef=',kslice_cntr_energy
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)')&
                     "bands=np.loadtxt('"//trim(seedname)//"-slice_bands.dat')"
                write(scriptunit,'(a)') 'numbands=bands.size/dimx/dimy'
                write(scriptunit,'(a)')&
                     'bands2=bands.reshape((dimx,dimy,numbands))'
                write(scriptunit,'(a)') 'for i in range(numbands):'
                write(scriptunit,'(a)')&
               "    pl.contour(np.transpose(bands2[:,:,i]),[ef],colors='black')"
             endif

             if(plot_curv_heatmap) then
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)') "outfile = '"//trim(seedname)//&
               "-curv_"//achar(119+i)//"-heatmap.pdf'"
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)')&
                     "z = np.loadtxt('"//trim(seedname)//&
                     "-slice_curv.dat', usecols=("//achar(47+i)//",))"
             elseif(plot_morb_heatmap) then
                write(scriptunit,'(a)') ' '
                write(scriptunit,'(a)') "outfile = '"//trim(seedname)//&
               "-morb_"//achar(119+i)//"-heatmap.pdf'"
!                write(scriptunit,'(a)')&
!                     "z = np.loadtxt('"//trim(seedname)//"-slice_morb.dat')"
                write(scriptunit,'(a)')&
                     "z = np.loadtxt('"//trim(seedname)//&
                     "-slice_morb.dat', usecols("//achar(47+i)//"))"
             endif
             write(scriptunit,'(a)') 'zz=z.reshape((dimx,dimy)).transpose()'
             write(scriptunit,'(a)') ' '

             if(plot_curv_heatmap) then
                ! Stepped color scale
                write(scriptunit,'(a)') 'pl.contourf(zz)'
             elseif(plot_morb_heatmap) then
                ! Gradual color scale
                write(scriptunit,'(a)')&
                     "imshow(zz,interpolation='bicubic',origin='lower')"
             endif

             write(scriptunit,'(a)') 'pl.colorbar()'
             write(scriptunit,'(a)') ' '
             write(scriptunit,'(a)') '# Remove the axes'
             write(scriptunit,'(a)') 'ax = pl.gca()'
             write(scriptunit,'(a)') 'ax.xaxis.set_visible(False)'
             write(scriptunit,'(a)') 'ax.yaxis.set_visible(False)'
             write(scriptunit,'(a)') ' '
             write(scriptunit,'(a)') 'pl.savefig(outfile)'
             write(scriptunit,'(a)') 'pl.show()'

          enddo

       endif !plot_curv_heatmap .or. plot_morb_heatmap
              
    end if ! on_root
 
end subroutine k_slice

end module w90_kslice
