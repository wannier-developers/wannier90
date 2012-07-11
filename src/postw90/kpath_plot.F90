!-*- mode: F90; mode: font-lock -*-!

module w90_kpath_plot

  ! Calculates one of the following along a specified k-path:
  ! 
  !  - Energy bands (eventually colored by the spin) 
  !
  !  - Berry curvature summed over occupied bands
  !
  !  - Integrand of orbital magnetization Morb=LCtil+ICtil

  implicit none

  public

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       ! 
  !===========================================================!

  subroutine k_path

    use w90_comms
    use w90_constants,  only   : dp,cmplx_0,cmplx_i,twopi
    use w90_io,         only   : io_error,io_file_unit,seedname,&
         io_time,io_stopwatch,stdout
    use w90_utility, only      : utility_diagonalize
    use w90_wanint_common, only : fourier_R_to_k
    use w90_parameters, only   : num_wann,recip_metric,kpath_task,&
         kpath_num_points,bands_num_spec_points,&
         bands_spec_points,bands_label,bands_color,&
         band_by_band,found_fermi_energy,fermi_energy,&
         omega_from_FF
    use w90_get_oper, only     : get_HH_R,HH_R,get_AA_R,get_BB_R,get_CC_R,&
         get_FF_R,get_SS_R
    use w90_spin_wanint, only  : get_spn_nk
    use w90_berry_wanint, only : get_imf_ab_k,get_img_ab_k,get_imh_ab_k
    use w90_constants, only    : bohr,ev_au

    integer           :: i,num_paths,num_spts,loop_path,loop_kpt,&
         total_pts,counter,loop_i,dataunit,gnuunit,&
         kpath_pts(bands_num_spec_points/2)
    real(kind=dp)     :: ymin,ymax,vec(3),kpt(3),spn_nk(num_wann),&
         imf_ab_k(3),img_ab_k(3),imh_ab_k(3),&
         LCtil(3),ICtil(3),&
         kpath_len(bands_num_spec_points/2),range 
    logical           :: plot_bands,plot_curv,plot_orb,got_it
    character(len=20) :: file_name

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    real(kind=dp), allocatable    :: xval(:),curv_au(:),curv_decomp_au(:,:),&
         orb_au(:,:),orb_decomp_au(:,:),eig_n(:,:),&
         color_n(:,:),plot_kpoint(:,:)
    character(len=3),allocatable  :: glabel(:)

    ! Everything is done on the root node (not worthwhile parallelizing) 
    ! However, we still have to read and distribute the data if we 
    ! are in parallel. So calls to get_oper are done on all nodes at the moment

    got_it=.false.
    plot_bands=.false.
    if(index(kpath_task,'bands')>0) then
       plot_bands=.true.
       got_it=.true.
    end if
    plot_curv=.false.
    if(index(kpath_task,'curv')>0) then
       plot_curv=.true.
       got_it=.true.
    end if
    plot_orb=.false.
    if(index(kpath_task,'orb')>0) then
       plot_orb=.true.
       got_it=.true.
    end if
    if(.not.got_it) then
       call io_error(&
            '(kpath_task must include one or more of the keywords "bands" "curv" and "orb"'&
            )
       stop
    end if

    if(on_root) then
       write(stdout,'(/,/,1x,a)') '=============================='
       write(stdout,'(1x,a)')     'Plotting along a k-point path:'
       write(stdout,'(1x,a)')     '=============================='

       if(plot_bands) then
          select case(bands_color)
          case("none")
             write(stdout,'(/,3x,a)') '* Energy bands in eV'
          case("spin")
             write(stdout,'(/,3x,a)') '* Energy bands in eV colored by spin'
          end select
       end if
       if(plot_curv) then
          write(stdout,'(/,3x,a)') '* Total Berry curvature in a.u.'
          if(.not.found_fermi_energy) call io_error&
               (&
               'Need to set either "fermi_energy" or "num_elec_cell" when plot_curv=T'&
               )
       end if
       if(plot_orb) then
          write(stdout,'(/,3x,a)')& 
               '* LC_tilde, IC_tilde  and total magnetization in a.u.'
          if(.not.found_fermi_energy) call io_error&
               (&
               'Need to set either "fermi_energy" or "num_elec_cell" when plot_orb=T'&
               )
       end if
    endif

    ! Set up the needed Wannier matrix elements
    call get_HH_R
    if(plot_curv.or.plot_orb) then
       call get_AA_R
       if(omega_from_FF) call get_FF_R
    endif
    if(plot_orb) then
       call get_BB_R
       call get_CC_R
    endif
    if(plot_bands .and. bands_color=='spin') call get_SS_R

    if(on_root) then

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

       allocate(HH(num_wann,num_wann))
       allocate(UU(num_wann,num_wann))

       ! Value of the vertical coordinate in the actual plots: energy bands 
       !
       if(plot_bands) then
          allocate(eig_n(num_wann,total_pts))
          if(bands_color/='none') then
             allocate(color_n(num_wann,total_pts))
          end if
       end if

       ! Value of the vertical coordinate in the actual plots: total Berry
       ! curvature summed over occupied states
       !
       if(plot_curv) then
          allocate(curv_au(total_pts))
          allocate(curv_decomp_au(total_pts,3))
       end if
       if(plot_orb) then
          allocate(orb_au(total_pts,3))
          allocate(orb_decomp_au(total_pts,9))
       end if

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

       ! Write out the kpoints in the path in a format that can be inserted
       ! directly in the pwscf input file (the '1.0_dp' in the second column is 
       ! a k-point weight, expected by pwscf)
       !
       dataunit=io_file_unit()
       open(dataunit,file=trim(seedname)//'_band.kpt',form='formatted')
       write(dataunit,*) total_pts
       do loop_kpt=1,total_pts
          write(dataunit,'(3f12.6,3x,f4.1)')&
               (plot_kpoint(loop_i,loop_kpt),loop_i=1,3),1.0_dp
       end do
       close(dataunit)

       ! Loop over k-points on the path and evaluate the requested quantities
       !
       do loop_kpt=1,total_pts       
          kpt(:)=plot_kpoint(:,loop_kpt)

          ! Energy bands
          !
          if(plot_bands) then
             call fourier_R_to_k(kpt,HH_R,HH,0)
             call utility_diagonalize(HH,num_wann,eig_n(:,loop_kpt),UU)
             !
             ! Color-code energy bands with the spin projection along the
             ! chosen quantization axis
             !
             if(bands_color=='spin') then
                call get_spn_nk(kpt,spn_nk)
                color_n(:,loop_kpt)=spn_nk(:)
             end if
          endif

          ! Berry curvature summed over occupied states
          !
          if(plot_curv) then
             call get_imf_ab_k(kpt,imf_ab_k)
             ! 
             ! Use atomic units, to facilitate comparison with WYSV06
             !
             curv_au(loop_kpt)=sum(imf_ab_k)/bohr**2
             !
             ! Decompose into J0 (Omega_bar), J1 (DA) and J2 (DD) parts
             !
             curv_decomp_au(loop_kpt,1)=imf_ab_k(1)/bohr**2
             curv_decomp_au(loop_kpt,2)=imf_ab_k(2)/bohr**2
             curv_decomp_au(loop_kpt,3)=imf_ab_k(3)/bohr**2
          end if

          if(plot_orb) then

             call get_imf_ab_k(kpt,imf_ab_k)
             call get_img_ab_k(kpt,img_ab_k)
             call get_imh_ab_k(kpt,imh_ab_k)

             LCtil(:)=img_ab_k(:)-fermi_energy*imf_ab_k(:)
             orb_au(loop_kpt,1)=sum(LCtil)
             ICtil(:)=imh_ab_k(:)-fermi_energy*imf_ab_k(:)
             orb_au(loop_kpt,2)=sum(ICtil)
             orb_au(loop_kpt,3)=orb_au(loop_kpt,1)+orb_au(loop_kpt,2)

             ! Decompose each into J0, J1, and J2 parts
             !
             ! J0 part of LCtil
             orb_decomp_au(loop_kpt,1)=LCtil(1)
             ! J1 part of LCtil
             orb_decomp_au(loop_kpt,2)=LCtil(2)
             ! J2 part of LCtil
             orb_decomp_au(loop_kpt,3)=LCtil(3)
             !
             ! J0 part of ICtil
             orb_decomp_au(loop_kpt,4)=ICtil(1)
             ! J1 part of ICtil
             orb_decomp_au(loop_kpt,5)=ICtil(2)
             ! J2 part of ICtil
             orb_decomp_au(loop_kpt,6)=ICtil(3)
             !
             ! J0 part of M_tot
             orb_decomp_au(loop_kpt,7)=LCtil(1)+ICtil(1)
             ! J1 part of M_tot
             orb_decomp_au(loop_kpt,8)=LCtil(2)+ICtil(2)
             ! J2 part of M_tot
             orb_decomp_au(loop_kpt,9)=LCtil(3)+ICtil(3)

             ! Convert to atomic units (right now it is an energy in eV
             ! times a Berry curvature in Ang^2)
             !
             orb_au(loop_kpt,:)=orb_au(loop_kpt,:)*eV_au/bohr**2
             orb_decomp_au(loop_kpt,:)=&
                  orb_decomp_au(loop_kpt,:)*eV_au/bohr**2
          end if

       end do !loop_kpt

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

       write(stdout,'(/,/,1x,a)') '------------------'
       write(stdout,'(1x,a)')     'Output data files:'
       write(stdout,'(1x,a)')     '------------------'

       ! This file (created earlier) is written no matter which task
       !
       file_name=trim(seedname)//'_band.kpt'
       write(stdout,'(/,3x,a)') file_name

       if(plot_bands) then

          ! Data file, gnuplot format
          !
          dataunit=io_file_unit()
          file_name=trim(seedname)//'_band.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do i=1,num_wann
             do loop_kpt=1,total_pts
                if(bands_color=='none') then
                   write(dataunit,'(2E16.8)') xval(loop_kpt),eig_n(i,loop_kpt)
                else
                   write(dataunit,'(3E16.8)') xval(loop_kpt),&
                        eig_n(i,loop_kpt),color_n(i,loop_kpt)  
                end if
             enddo
             write(dataunit,*) ' '
          enddo
          close(dataunit)

          ymin=minval(eig_n)-1.0_dp
          ymax=maxval(eig_n)+1.0_dp

          ! Gnuplot script
          !
          gnuunit=io_file_unit()
          file_name=trim(seedname)//'_band.gnu'
          write(stdout,'(/,3x,a)') file_name
          open(gnuunit,file=file_name,form='formatted')
          do i = 1,num_paths-1
             write(gnuunit,705) sum(kpath_len(1:i)),ymin,sum(kpath_len(1:i)),&
                  ymax
          enddo
          if(bands_color=='none') then
             write(gnuunit,701) xval(total_pts),ymin,ymax
             write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
                  (glabel(i+1),sum(kpath_len(1:i)),i=1,num_paths-1)
             write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
             write(gnuunit,*) 'plot ','"'//trim(seedname)//'_band.dat','"' 
          elseif(bands_color=='spin') then
             !
             ! Only works with gnuplot v4 (4.2?) and higher
             !
             write(gnuunit,706) xval(total_pts),ymin,ymax
             write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
                  (glabel(i+1),sum(kpath_len(1:i)),i=1,num_paths-1)
             write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
             write(gnuunit,*)&
                  'set palette defined (-1 "red", 0 "green", 1 "blue")'
             write(gnuunit,*) 'set pm3d map'
             write(gnuunit,*) 'set zrange [-1:1]'
             write(gnuunit,*) 'splot ','"'//trim(seedname)//'_band.dat',& 
                  '" with dots palette' 
!          elseif(bands_color=='curv') then
!             write(gnuunit,706) xval(total_pts),ymin,ymax
!             write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
!                  (glabel(i+1),sum(kpath_len(1:i)),i=1,num_paths-1)
!             write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
!             write(gnuunit,704) minval(color_n),maxval(color_n)
!             write(gnuunit,*) 'set pm3d map'
!             write(gnuunit,*) 'splot ','"'//trim(seedname)//'_band.dat',& 
!                  '" with dots palette' 
          end if
          close(gnuunit)

       end if ! plot_bands

       if(plot_curv) then

          dataunit=io_file_unit()
          file_name=trim(seedname)//'_curv.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do loop_kpt=1,total_pts
             write(dataunit,'(2E16.8)') xval(loop_kpt),&
                  curv_au(loop_kpt)
          enddo
          write(dataunit,*) ' '
          close(dataunit)
          !
          dataunit=io_file_unit()
          file_name=trim(seedname)//'_curv_decomp.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do loop_kpt=1,total_pts
             write(dataunit,'(4E16.8)')&
                  xval(loop_kpt),curv_decomp_au(loop_kpt,:)
          enddo
          write(dataunit,*) ' '
          close(dataunit)

          ymin=minval(curv_au)
          ymax=maxval(curv_au)
          range=ymax-ymin
          ymin=ymin-0.02_dp*range
          ymax=ymax+0.02_dp*range

          gnuunit=io_file_unit()
          file_name=trim(seedname)//'_curv.gnu'
          write(stdout,'(/,3x,a)') file_name
          open(gnuunit,file=file_name,form='formatted')
          write(gnuunit,707) xval(total_pts),ymin,ymax
          do i = 1, num_paths-1
             write(gnuunit,705) sum(kpath_len(1:i)),ymin,sum(kpath_len(1:i)),&
                  ymax
          enddo
          write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
               (glabel(i+1),sum(kpath_len(1:i)),i=1,num_paths-1)
          write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
          write(gnuunit,*) 'plot ','"'//trim(seedname)//'_curv.dat','"' 
          close(gnuunit)

       end if ! plot_curv

       if(plot_orb) then

          dataunit=io_file_unit()
          file_name=trim(seedname)//'_orb.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do loop_kpt=1,total_pts
             write(dataunit,'(4E16.8)') xval(loop_kpt),orb_au(loop_kpt,:)
          enddo
          write(dataunit,*) ' '
          close(dataunit)
          !
          dataunit=io_file_unit()
          file_name=trim(seedname)//'_orb_decomp.dat'
          write(stdout,'(/,3x,a)') file_name
          open(dataunit,file=file_name,form='formatted')
          do loop_kpt=1,total_pts
             write(dataunit,'(10E16.8)')&
                  xval(loop_kpt),orb_decomp_au(loop_kpt,:)
          enddo
          write(dataunit,*) ' '
          close(dataunit)

          ymin=minval(orb_au(:,:))
          ymax=maxval(orb_au(:,:))
          range=ymax-ymin
          ymin=ymin-0.02_dp*range
          ymax=ymax+0.02_dp*range

          gnuunit=io_file_unit()
          file_name=trim(seedname)//'_orb.gnu'
          write(stdout,'(/,3x,a)') file_name
          open(gnuunit,file=file_name,form='formatted')
          write(gnuunit,707) xval(total_pts),ymin,ymax
          do i = 1, num_paths-1
             write(gnuunit,705) sum(kpath_len(1:i)),ymin,sum(kpath_len(1:i)),&
                  ymax
          enddo
          write(gnuunit,702, advance="no") glabel(1),0.0_dp,&
               (glabel(i+1),sum(kpath_len(1:i)),i=1,num_paths-1)
          write(gnuunit,703) glabel(1+num_paths),sum(kpath_len(:))
          ! u 1:4 plots the total M(k) along the path (LCtil+ICtil)
          write(gnuunit,*) 'plot ','"'//trim(seedname)//'_orb.dat','" u 1:4' 
          close(gnuunit)

       end if ! plot_orb

    end if ! on_root

701 format('set style data dots',/,'set nokey',/,&
         'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
702 format('set xtics (',:20('"',A3,'" ',F8.5,','))
703 format(A3,'" ',F8.5,')')
704 format('set palette defined (',F8.5,' "red", 0 "green", ',F8.5,' "blue")')
705 format('set arrow from ',F16.8,',',F16.8,' to ',F16.8,',',F16.8, ' nohead')
706 format('set nokey',/,&
         'set xrange [0:',F9.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
707 format('set style data lines',/,'set nokey',/,&
         'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')

  end subroutine k_path

end module w90_kpath_plot
