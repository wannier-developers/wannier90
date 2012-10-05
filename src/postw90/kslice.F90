  !-*- mode: F90; mode: font-lock -*-!

module w90_kslice

  ! Makes a heatmap plot on a slice in k-space of one of the following:
  ! 
  !  - Minus the Berry curvature summed over occupied bands 
  !
  !  - Integrand of orbital magnetization formula (LCtil, ICtil, and Omega)
  !
  ! Calculates the intersections of the constant-energy surfaces with the slice
  !
  ! ******
  ! TO DO: Add the possibility to color the isoenergy contours by the spin 
  ! ****** 
  !
  ! The slice is defined by three input variables, all in reciprocal
  ! lattice coordinates:
  !
  !    slice_corner(1:3) is the lower left corner 
  !    slice_b1(1:3) and slice_b2(1:3) are the vectors subtending the slice


  implicit none

  public

  contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine k_slice

    use w90_comms
    use w90_constants,  only   : dp,twopi
    use w90_io,         only   : io_error,io_file_unit,seedname,&
                                 io_time,io_stopwatch,stdout
    use w90_utility, only     : utility_diagonalize
    use w90_postw90_common, only : fourier_R_to_k
    use w90_parameters, only   : num_wann,kslice,kslice_task,&
                                 kslice_num_points,&
                                 kslice_corner,kslice_b1,kslice_b2,&
                                 kslice_cntr_energy,&
                                 found_kslice_cntr_energy,&
                                 recip_lattice,found_fermi_energy,fermi_energy
    use w90_get_oper, only     : get_HH_R,HH_R,get_AA_R,get_BB_R,get_CC_R
    use w90_berry, only : get_imf_ab_k,get_img_ab_k,get_imh_ab_k
    use w90_utility, only      : utility_recip_lattice
    use w90_constants, only    : bohr,ev_au

    integer           :: loop_kpt,loop_x,loop_y,n,n1,n2,n3
    integer           :: xdataunit,ydataunit,zdataunit,bandsunit,gnuunit
    real(kind=dp)     :: avec(3,3),bvec(3,3),recip_vol,b1mod,b2mod,cosb1b2,&
                         a2mod,kpt(3),kpt_x,kpt_y,k1,k2,k_cart(3),&
                         imf_ab_k(3),img_ab_k(3),imh_ab_k(3),&
                         curv_au,morb_au(3)
    logical           :: plot_curv_heatmap,plot_morb_heatmap,plot_energy_cntr,&
                         got_it
    character(len=20) :: filename

    integer, allocatable          :: xyzdataunit(:)
    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    real(kind=dp),    allocatable :: eig(:)    !,wan_c_2(:,:,:)

    ! Everything is done on the root node (not worthwhile parallelizing). 
    !
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
    endif

    got_it=.false.
    plot_curv_heatmap=.false.
    if(index(kslice_task,'curv_heatmap')>0) then
       plot_curv_heatmap=.true.
       got_it=.true.
    end if
    plot_morb_heatmap=.false.
    if(index(kslice_task,'morb_heatmap')>0) then
       plot_morb_heatmap=.true.
       got_it=.true.
    end if
    plot_energy_cntr=.false.
    if(index(kslice_task,'energy_cntr')>0) then
       plot_energy_cntr=.true.
       got_it=.true.
    end if
    if(plot_curv_heatmap .and. plot_morb_heatmap)  then
       call io_error(&
           '(kslice_task cannot include both "curv_heatmap" and "morb_heatmap"')
       stop
    endif
    if(.not.got_it) then
       call io_error(&
            '(kslice_task must include only one of the keywords "curv_heatmap" and "morb_heatmap," and/or "energy_cntr"')
       stop
    end if

    ! Set up the needed Wannier matrix elements
    call get_HH_R
    if(plot_curv_heatmap .or. plot_morb_heatmap) call get_AA_R
    if(plot_morb_heatmap) then
       call get_BB_R
       call get_CC_R
    endif
    
    if(on_root) then

       write(stdout,'(/,/,1x,a)') '============================'
       write(stdout,'(1x,a)')     'Plotting on a k-point slice:'
       write(stdout,'(1x,a)')     '============================'

       ! ***TO DO*** Write in stdout vectors defining slice, in latt and cart

       ! octave and python+matplotlib: write the x- and y-axis mesh values
       xdataunit=io_file_unit() 
       filename=trim(seedname)//'_slice_x.dat'
       write(stdout,'(/,3x,a)') filename
       open(xdataunit,file=filename,form='formatted')
       ydataunit=io_file_unit() 
       filename=trim(seedname)//'_slice_y.dat'
       write(stdout,'(/,3x,a)') filename
       open(ydataunit,file=filename,form='formatted')

       if(plot_curv_heatmap) then
          write(stdout,'(/,3x,a)') '* Berry curvature in a.u.'
          if(.not.found_fermi_energy) call io_error&
               (&
 'Need to set either "fermi_energy" or "num_elec_cell" when plot_curv_heatmap=T'&
               )
          zdataunit=io_file_unit()
          filename=trim(seedname)//'_slice_curv.dat'
          write(stdout,'(/,3x,a)') filename
          open(zdataunit,file=filename,form='formatted')
       elseif(plot_morb_heatmap) then
          write(stdout,'(/,3x,a)')&
  '* LC_tilde-mu.Omega, IC_tilde-mu.Omega, and their sum (total M_orb) in a.u.'
          if(.not.found_fermi_energy) call io_error&
               (&
'Need to set either "fermi_energy" or "num_elec_cell" when plot_morb_heatmap=T'&
               )
          zdataunit=io_file_unit()
          filename=trim(seedname)//'_slice_morb.dat'
          write(stdout,'(/,3x,a)') filename
          open(zdataunit,file=filename,form='formatted')
       end if

       if(plot_energy_cntr) then
          write(stdout,'(/,3x,a)') '* Fermi surface contours'
          if(.not.found_kslice_cntr_energy) call io_error&
               ('Error: must specify either fermi_energy or kslice_cntr_energy when kslice_task = energy_cntr')
          write(stdout,'(/,1x,a,f12.8,1x,a)')&
               'Energy isocontour level: ',kslice_cntr_energy,'eV'
          allocate(HH(num_wann,num_wann))
          allocate(UU(num_wann,num_wann))
          allocate(eig(num_wann))
          bandsunit=io_file_unit()
          filename=trim(seedname)//'_slice_bands.dat'
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
     
       ! Loop over uniform mesh of k-points on the slice
       !
       do loop_kpt=0,kslice_num_points**2-1
          loop_x=loop_kpt/kslice_num_points
          loop_y=loop_kpt-loop_x*kslice_num_points
          kpt(1)=kslice_corner(1)&
               +kslice_b1(1)*real(loop_x,dp)/kslice_num_points&
               +kslice_b2(1)*real(loop_y,dp)/kslice_num_points
          kpt(2)=kslice_corner(2)&
               +kslice_b1(2)*real(loop_x,dp)/kslice_num_points&
               +kslice_b2(2)*real(loop_y,dp)/kslice_num_points
          kpt(3)=kslice_corner(3)&
               +kslice_b1(3)*real(loop_x,dp)/kslice_num_points&
               +kslice_b2(3)*real(loop_y,dp)/kslice_num_points

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
             call get_imf_ab_k(kpt,imf_ab_k)
             curv_au=sum(imf_ab_k)/bohr**2
   
             ! The goal here is to make a plot as in Fig. 3 of 
             ! Yao et al, PRL 92, 037204 (2004)
             !
             ! First flip the sign, since they plot *minus* the Berry curvature
             !
             curv_au=-curv_au
             if(curv_au>0) then
                if(curv_au>10) then
                   curv_au=log10(curv_au)
                else
                   curv_au=1.0_dp
                endif
             else
                if(curv_au<-10) then
                   curv_au=-log10(-curv_au)
                else
!                   curv_au=-1.0_dp
                   curv_au=-0.99_dp ! -1.0 would not work!
                endif
             endif
             write(zdataunit,'(E16.8)') curv_au
          end if

          if(plot_morb_heatmap) then
             call get_imf_ab_k(kpt,imf_ab_k)
             call get_img_ab_k(kpt,img_ab_k)
             call get_imh_ab_k(kpt,imh_ab_k)
             morb_au(1)=sum(img_ab_k(:)-fermi_energy*imf_ab_k(:))
             morb_au(2)=sum(imh_ab_k(:)-fermi_energy*imf_ab_k(:))
             morb_au(3)=morb_au(1)+morb_au(2)
             morb_au=morb_au*eV_au/bohr**2
             write(zdataunit,'(3E16.8)') morb_au
          end if

          if(plot_energy_cntr) then
             call fourier_R_to_k(kpt,HH_R,HH,0)
             call utility_diagonalize(HH,num_wann,eig,UU)
!             write(bandsunit,'(99E16.8)') eig
             do n=1,num_wann
                ! For python/octave
                write(bandsunit,'(E16.8)') eig(n)
                ! For gnuplot, using 'grid data' format
                write(xyzdataunit(n),'(3E16.8)') kpt_x,kpt_y,eig(n)
                if(loop_y==kslice_num_points-1 .and. &
                   loop_x/=kslice_num_points-1) write (xyzdataunit(n),*) ' '
             enddo
          endif

       end do !loop_kpt
       
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

       gnuunit=io_file_unit()
       open(gnuunit,file=trim(seedname)//'_energy_cntr.gnu',form='formatted')
       write(gnuunit,'(a)') 'unset surface'
       write(gnuunit,'(a)') 'set contour'
       write(gnuunit,'(a)') 'set view map'
       write(gnuunit,'(a,f9.5)') 'set cntrparam levels discrete ',&
            kslice_cntr_energy
       write(gnuunit,'(a)') 'set cntrparam bspline'
       do n=1,num_wann
          n1=n/100
          n2=(n-n1*100)/10
          n3=n-n1*100-n2*10
          write(gnuunit,'(a)') 'set table "xyz_'&
               //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat"'
          write(gnuunit,'(a)') 'splot "'//trim(seedname)//'_xyz_'&
               //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat"'
          write(gnuunit,'(a)') 'unset table'
       enddo
       write(gnuunit,'(a)') '#Uncomment next two lines to create postscript'
       write(gnuunit,'(a)') '#set term post eps enh'
       write(gnuunit,'(a)') '#set output "'//trim(seedname)//'_econtour.eps"'
       write(gnuunit,'(a)') 'set size ratio -1'
       write(gnuunit,'(a)') 'unset tics'
       write(gnuunit,'(a)') 'set nokey'
       write(gnuunit,'(a)')&
            '#For postscript try changing lw 1 --> lw 2 in the next line'
       write(gnuunit,'(a)') 'set style line 1 lt 1 lw 1'
       if(num_wann==1) then
          write(gnuunit,'(a)') 'plot "xyz_001.dat" using 1:2 w lines ls 1'
       else
          write(gnuunit,'(a)') 'plot "xyz_001.dat" using 1:2 w lines ls 1,'&
               //achar(92)
       endif
       do n=2,num_wann-1
          n1=n/100
          n2=(n-n1*100)/10
          n3=n-n1*100-n2*10
          write(gnuunit,'(a)') '     "xyz_'&
               //achar(48+n1)//achar(48+n2)//achar(48+n3)&
               //'.dat" using 1:2 w lines ls 1,'//achar(92)
       enddo
       n=num_wann
       n1=n/100
       n2=(n-n1*100)/10
       n3=n-n1*100-n2*10
       write(gnuunit,'(a)') '     "xyz_'&
            //achar(48+n1)//achar(48+n2)//achar(48+n3)&
            //'.dat" using 1:2 w lines ls 1'
    
    end if ! on_root

!701 format('set style data dots',/,'set nokey',/,&
!         'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
!702 format('set xtics (',:20('"',A3,'" ',F8.5,','))
!703 format(A3,'" ',F8.5,')')
!704 format('set palette defined (',F8.5,' "red", 0 "green", ',F8.5,' "blue")')
!705 format('set arrow from ',F16.8,',',F16.8,' to ',F16.8,',',F16.8, ' nohead')
!706 format('set nokey',/,&
!         'set xrange [0:',F9.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
!707 format('set style data lines',/,'set nokey',/,&
!         'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
 
end subroutine k_slice

end module w90_kslice
