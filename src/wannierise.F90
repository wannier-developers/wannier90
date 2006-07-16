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

module w90_wannierise

  use w90_constants

  implicit none

  private

  public :: wann_main

  ! Data to avoid large allocation within iteration loop
  complex(kind=dp), allocatable  :: cr (:,:,:,:)   
  complex(kind=dp), allocatable  :: crt (:,:,:,:)  
  real(kind=dp),    allocatable  :: rnkb (:,:,:)   
  real(kind=dp),    allocatable  :: ln_tmp(:,:,:)

  type localisation_vars
     real(kind=dp) :: om_i   
     real(kind=dp) :: om_d   
     real(kind=dp) :: om_od  
     real(kind=dp) :: om_tot 
!!$     real(kind=dp) :: om_1   
!!$     real(kind=dp) :: om_2   
!!$     real(kind=dp) :: om_3   
  end type localisation_vars


contains

  !==================================================================!
  subroutine wann_main
    !==================================================================!
    !                                                                  !
    ! Calculate the Unitary Rotations to give                          !
    !            Maximally Localised Wannier Functions                 !
    !                                                                  !
    !===================================================================  
    use w90_constants,  only : dp,cmplx_1,cmplx_0
    use w90_io,         only : stdout,io_error,io_time,io_stopwatch,io_file_unit
    use w90_parameters, only : num_wann,num_cg_steps,num_iter,wb,nnlist, &
         nntot,wbtot,u_matrix,m_matrix,num_kpts,iprint, &
         num_print_cycles,num_dump_cycles,omega_invariant, &
         param_read_um,param_write_um,length_unit,lenconfac, &
         proj_site,real_lattice,write_r2mn,guiding_centres,&
         num_guide_cycles,num_no_guide_iter,&
         trial_step,fixed_step,lfixstep,write_proj,have_disentangled
    use w90_utility,    only : utility_frac_to_cart,utility_zgemm

    implicit none


    type(localisation_vars) :: old_spread
    type(localisation_vars) :: wann_spread
    type(localisation_vars) :: trial_spread


    ! guiding centres
    real(kind=dp), allocatable :: rguide (:,:)  
    integer :: irguide

    ! local arrays used and passed in subroutines
    complex(kind=dp), allocatable :: csheet(:,:,:)
    complex(kind=dp), allocatable :: cdodq(:,:,:)  
    complex(kind=dp), allocatable :: cdodq1(:,:,:)  
    complex(kind=dp), allocatable :: cdodq2(:,:,:)  
    complex(kind=dp), allocatable :: cdodq3(:,:,:)  
    real(kind=dp),    allocatable :: sheet (:,:,:)
    real(kind=dp),    allocatable :: rave(:,:),r2ave(:),rave2(:)  

    !local arrays not passed into subroutines
    complex(kind=dp), allocatable  :: cwschur1 (:), cwschur2 (:)  
    complex(kind=dp), allocatable  :: cwschur3 (:), cwschur4 (:)  
    complex(kind=dp), allocatable  :: cdq(:,:,:),cdqkeep(:,:,:)  
    complex(kind=dp), allocatable  :: cz (:,:)  
    complex(kind=dp), allocatable  :: cmtmp(:,:),tmp_cdq(:,:) 
    complex(kind=dp), allocatable  :: m0(:,:,:,:),u0(:,:,:)
    complex(kind=dp), allocatable  :: cwork(:)
    real(kind=dp),    allocatable  :: evals(:)
    real(kind=dp),    allocatable  :: rwork(:)

    real(kind=dp) :: doda0
    real(kind=dp) :: falphamin,alphamin
    real(kind=dp) :: gcfac,gcnorm1,gcnorm0
    integer :: i,n,iter,ind,ierr,iw,ncg,info
    logical :: lprint,ldump,lquad


    ! Allocate stuff

    ! module data
    allocate(  m0 (num_wann, num_wann, nntot, num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating m0 in wann_main')
    allocate(  u0 (num_wann, num_wann, num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating u0 in wann_main')
    allocate(  cr (num_wann, num_wann, nntot, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cr in wann_main')
    allocate(  crt (num_wann, num_wann, nntot, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating crt in wann_main')
    allocate( rnkb (num_wann, nntot, num_kpts),stat=ierr    )     
    if (ierr/=0) call io_error('Error in allocating rnkb in wann_main')
    allocate( ln_tmp (num_wann, nntot, num_kpts), stat=ierr    )
    if (ierr/=0) call io_error('Error in allocating ln_tmp in wann_main')

    cr=cmplx_0;  crt=cmplx_0;  rnkb=cmplx_0

    ! sub vars passed into other subs 
    allocate( csheet (num_wann, nntot, num_kpts), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating csheet in wann_main')
    allocate( cdodq (num_wann, num_wann, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cdodq in wann_main')
    allocate( cdodq1 (num_wann, num_wann, num_kpts),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cdodq1 in wann_main')
    allocate( cdodq2 (num_wann, num_wann, num_kpts),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cdodq2 in wann_main')
    allocate( cdodq3 (num_wann, num_wann, num_kpts),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cdodq3 in wann_main')
    allocate( sheet (num_wann, nntot, num_kpts), stat=ierr    )
    if (ierr/=0) call io_error('Error in allocating sheet in wann_main')
    allocate( rave (3, num_wann), r2ave (num_wann),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating rave in wann_main')
    allocate( rave2 (num_wann),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating rave2 in wann_main')
    allocate( rguide (3, num_wann)   )
    if (ierr/=0) call io_error('Error in allocating rguide in wann_main')

    csheet=cmplx_1;cdodq=cmplx_0;cdodq1=cmplx_0;cdodq2=cmplx_0;cdodq3=cmplx_0
    rave2=0.0_dp;sheet=0.0_dp;rave=0.0_dp;r2ave=0.0_dp;sheet=0.0_dp

    ! sub vars not passed into other subs
    allocate( cwschur1 (num_wann), cwschur2 (10 * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cwshur1 in wann_main')
    allocate( cwschur3 (num_wann), cwschur4 (num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cwshur3 in wann_main')
    allocate( cdq (num_wann, num_wann, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cdq in wann_main')
    allocate( cz (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cz in wann_main')
    allocate( cmtmp (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cmtmp in wann_main')
    allocate( cdqkeep (num_wann, num_wann, num_kpts),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cdqkeep in wann_main')
    allocate(tmp_cdq(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tmp_cdq in wann_main')
    allocate( evals (num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating evals in wann_main')
    allocate( cwork (4*num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating cwork in wann_main')
    allocate( rwork (3*num_wann-2),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating rwork in wann_main')


    cwschur1=cmplx_0; cwschur2=cmplx_0; cwschur3=cmplx_0; cwschur4=cmplx_0
    cdq=cmplx_0; cz=cmplx_0; cmtmp=cmplx_0; cdqkeep=cmplx_0
    
    gcnorm1=0.0_dp; gcnorm0=0.0_dp

    ! initialise rguide to projection centres (Cartesians in units of Ang)
    if( guiding_centres) then
       do n=1,num_wann
          call utility_frac_to_cart(proj_site(:,n),rguide(:,n),real_lattice)
        enddo
    end if

    write(stdout,*)
    write(stdout,'(1x,a)') '*------------------------------- WANNIERISE ---------------------------------*'
    write(stdout,'(1x,a)') '+--------------------------------------------------------------------+<-- CONV'
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(1x,a)') '| Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV'
    else
       write(stdout,'(1x,a)') '| Iter  Delta Spread     RMS Gradient      Spread (Bohr^2)     Time  |<-- CONV'
    endif
    write(stdout,'(1x,a)') '+--------------------------------------------------------------------+<-- CONV'
    write(stdout,*)


    irguide=0
    if (guiding_centres.and.(num_no_guide_iter.le.0)) then
       call phases(csheet,sheet,rguide,irguide)
       irguide=1
    endif

    ! calculate initial centers and spread
    call omega(csheet,sheet,rave,r2ave,rave2,wann_spread)
    omega_invariant = wann_spread%om_i

    if (lfixstep) lquad=.false.
    ncg  = 0
    iter = 0
    old_spread%om_tot = 0.0_dp

    ! print initial state
    write(stdout,'(1x,a78)') repeat('-',78) 
    write(stdout,'(1x,a)') 'Initial State'
    do iw=1,num_wann
       write(stdout,1000) iw,(rave(ind,iw)*lenconfac,ind=1,3),&
            (r2ave(iw) - rave2(iw))*lenconfac**2
    end do
    write(stdout,1001) (sum(rave(ind,:))*lenconfac,ind=1,3), (sum(r2ave)-sum(rave2))*lenconfac**2
    write(stdout,*)
    write(stdout,'(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
         iter,(wann_spread%om_tot-old_spread%om_tot)*lenconfac**2,sqrt(abs(gcnorm1))*lenconfac,&
         wann_spread%om_tot*lenconfac**2,io_time(),'<-- CONV'
    write(stdout,'(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
         'O_D=',wann_spread%om_d*lenconfac**2,' O_OD=',wann_spread%om_od*lenconfac**2,&
         ' O_TOT=',wann_spread%om_tot*lenconfac**2,' <-- SPRD'
    write(stdout,'(1x,a78)') repeat('-',78) 


    ! main iteration loop
    do iter=1,num_iter

       lprint=.false.
       if ( (mod(iter,num_print_cycles).eq.0) .or. (iter.eq.1) &
            .or. (iter.eq.num_iter)) lprint=.true.

       ldump=.false.
       if ( (num_dump_cycles.gt.0) .and. (mod(iter,num_dump_cycles).eq.0) ) ldump=.true.

       if(lprint) write(stdout,'(1x,a,i6)') 'Cycle: ',iter

       if ( guiding_centres.and.(iter.gt.num_no_guide_iter) & 
            .and.(mod(iter,num_guide_cycles).eq.0) ) then
          call phases(csheet,sheet,rguide,irguide)
          irguide=1
       endif

       ! calculate gradient of omega
       call domega(csheet,sheet,rave,cdodq1,cdodq2,cdodq3,cdodq)

       if ( lprint .and. iprint>2 ) &
            write(stdout,*) ' LINE --> Iteration                     :',iter

       ! calculate search direction (cdq)
       call internal_search_direction()

       ! save search direction 
       cdqkeep(:,:,:) = cdq(:,:,:)

       ! check whether we're doing fixed step lengths
       if (lfixstep) then

          alphamin=fixed_step

       ! or a parabolic line search
       else

          ! take trial step
          cdq(:,:,:)=cdqkeep(:,:,:)*( trial_step / (4.0_dp*wbtot) ) 
          
          ! store original U and M before rotating
          u0=u_matrix ; m0=m_matrix

          ! update U and M
          call internal_new_u_and_m()

          ! calculate spread at trial step
          call omega(csheet,sheet,rave,r2ave,rave2,trial_spread)

          ! Calculate optimal step (alphamin)
          call internal_optimal_step()

       endif

       ! print line search information
       if ( lprint .and. iprint>2 ) then
          write(stdout,*) ' LINE --> Spread at initial point       :',wann_spread%om_tot*lenconfac**2
          if (.not.lfixstep) &
               write(stdout,*) ' LINE --> Spread at trial step          :',trial_spread%om_tot*lenconfac**2
          write(stdout,*) ' LINE --> Slope along search direction  :',doda0*lenconfac**2
          write(stdout,*) ' LINE --> ||SD gradient||^2             :',gcnorm1*lenconfac**2
          if (.not.lfixstep) then
             write(stdout,*) ' LINE --> Trial step length             :',trial_step
             if (lquad) then
                write(stdout,*) ' LINE --> Optimal parabolic step length :',alphamin
                write(stdout,*) ' LINE --> Spread at predicted minimum   :',falphamin*lenconfac**2
             endif
          else
             write(stdout,*) ' LINE --> Fixed step length             :',fixed_step
          endif
             write(stdout,*) ' LINE --> CG coefficient                :',gcfac
       endif


       ! if taking a fixed step or if parabolic line search was successful
       if (lfixstep.or.lquad) then

          ! take optimal step
          cdq(:,:,:) = cdqkeep(:,:,:) * ( alphamin / (4.0_dp*wbtot) ) 
          
          ! if doing a line search then restore original U and M before rotating 
          if (.not.lfixstep) then 
             u_matrix=u0 ; m_matrix=m0
          endif

          ! update U and M
          call internal_new_u_and_m()
          
          call wann_spread_copy(wann_spread,old_spread)
          
          ! calculate the new centers and spread
          call omega(csheet,sheet,rave,r2ave,rave2,wann_spread)
        
       ! parabolic line search was unsuccessful, use trial step already taken
       else 

          call wann_spread_copy(trial_spread,wann_spread)
          call wann_spread_copy(wann_spread,old_spread)

       endif
 

       ! print the new centers and spreads
       if(lprint) then
          do iw=1,num_wann
             write(stdout,1000) iw,(rave(ind,iw)*lenconfac,ind=1,3),&
                  (r2ave(iw) - rave2(iw))*lenconfac**2
          end do
          write(stdout,1001) (sum(rave(ind,:))*lenconfac,ind=1,3), &
               (sum(r2ave)-sum(rave2))*lenconfac**2
          write(stdout,*)
          write(stdout,'(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
               iter,(wann_spread%om_tot-old_spread%om_tot)*lenconfac**2,&
               sqrt(abs(gcnorm1))*lenconfac,&
               wann_spread%om_tot*lenconfac**2,io_time(),'<-- CONV'
          write(stdout,'(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
               'O_D=',wann_spread%om_d*lenconfac**2,&
               ' O_OD=',wann_spread%om_od*lenconfac**2,&
               ' O_TOT=',wann_spread%om_tot*lenconfac**2,' <-- SPRD'
          write(stdout,'(1x,a,E15.7,a,E15.7,a,E15.7,a)') &
               'Delta: O_D=',(wann_spread%om_d-old_spread%om_d)*lenconfac**2,&
               ' O_OD=',(wann_spread%om_od-old_spread%om_od)*lenconfac**2,&
               ' O_TOT=',(wann_spread%om_tot-old_spread%om_tot)*lenconfac**2,' <-- DLTA'
          write(stdout,'(1x,a78)') repeat('-',78) 
       end if

       if (ldump) call param_write_um

    enddo
    ! end of the minimization loop


    write(stdout,'(1x,a)') 'Final State'
    do iw=1,num_wann
       write(stdout,1000) iw,(rave(ind,iw)*lenconfac,ind=1,3),&
            (r2ave(iw) - rave2(iw))*lenconfac**2
    end do
    write(stdout,1001) (sum(rave(ind,:))*lenconfac,ind=1,3),&
         (sum(r2ave)-sum(rave2))*lenconfac**2
    write(stdout,*)
    write(stdout,'(3x,a21,a,f15.9)') '     Spreads ('//trim(length_unit)//'^2)',&
         '       Omega I      = ',wann_spread%om_i*lenconfac**2
    write(stdout,'(3x,a,f15.9)') '     ================       Omega D      = ',&
         wann_spread%om_d*lenconfac**2
    write(stdout,'(3x,a,f15.9)') '                            Omega OD     = ',&
         wann_spread%om_od*lenconfac**2
    write(stdout,'(3x,a21,a,f15.9)') 'Final Spread ('//trim(length_unit)//'^2)',&
         '       Omega Total  = ',wann_spread%om_tot*lenconfac**2  
    write(stdout,'(1x,a78)') repeat('-',78) 


    if (guiding_centres) call phases(csheet,sheet,rguide,irguide)


    ! unitarity is checked
    call internal_check_unitarity()

    ! write extra info regarding omega_invariant
    if (iprint>2) call internal_svd_omega_i()

    ! write matrix elements <m|r^2|n> to file
    if (write_r2mn) call internal_write_r2mn()

    ! write U and M to file
    call param_write_um

    ! calculate and write projection of WFs on original bands in outer window
    if (have_disentangled .and. write_proj) call wann_calc_projection()

    ! deallocate sub vars not passed into other subs
    deallocate(rwork,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rwork in wann_main')
    deallocate(cwork,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cwork in wann_main')
    deallocate(evals,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating evals in wann_main')
    deallocate(tmp_cdq,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tmp_cdq in wann_main')
    deallocate(cdqkeep,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cdqkeep in wann_main')
    deallocate(cmtmp,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cmtmp in wann_main')
    deallocate(cz,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cz in wann_main')
    deallocate(cdq,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cdq in wann_main')
    deallocate(cwschur3,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cwschur3 in wann_main')
    deallocate(cwschur1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cwschur1 in wann_main')

    ! deallocate sub vars passed into other subs
    deallocate(rguide,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rguide in wann_main')
    deallocate(rave2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rave2 in wann_main')
    deallocate(rave,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rave in wann_main')
    deallocate(sheet,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sheet in wann_main')
    deallocate(cdodq3,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cdodq3 in wann_main')
    deallocate(cdodq2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cdodq2 in wann_main')
    deallocate(cdodq1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cdodq1 in wann_main')
    deallocate(cdodq,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cdodq in wann_main')
    deallocate(csheet,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating csheet in wann_main')
    
    ! deallocate module data
    deallocate( ln_tmp , stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating ln_tmp in wann_main')
    deallocate( rnkb,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating rnkb in wann_main')
    deallocate(  crt,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating crt in wann_main') 
    deallocate(  cr,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating cr in wann_main') 

    deallocate(u0, stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating u0 in wann_main')
    deallocate(m0, stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating m0 in wann_main')
    
    return

1000 format(2x,'WF centre and spread', &
         &       i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)

1001 format(2x,'Sum of centres and spreads', &
         &       1x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)


  contains



    !===============================================!
    subroutine internal_search_direction()
      !===============================================!
      !                                               !
      ! Calculate the conjugate gradients search      !
      ! direction using the Fletcher-Reeves formula:  !
      !                                               !
      !     cg_coeff = [g(i).g(i)]/[g(i-1).g(i-1)]    !
      !                                               !
      !===============================================!

      implicit none

      complex(kind=dp) :: zdotc

      ! gcnorm1 = Tr[gradient . gradient] -- NB gradient is anti-Hermitian
      gcnorm1 = real(zdotc(num_kpts*num_wann*num_wann,cdodq,1,cdodq,1),dp)

      ! calculate cg_coefficient
      if ( (iter.eq.1) .or. (ncg.ge.num_cg_steps) ) then
         gcfac = 0.0_dp                 ! Steepest descents
         ncg   = 0
      else
         if (gcnorm0.gt.epsilon(1.0_dp)) then
            gcfac = gcnorm1/gcnorm0     ! Fletcher-Reeves CG coefficient
            ! prevent CG coefficient from getting too large
            if (gcfac.gt.3.0_dp) then
               if ( lprint .and. iprint>2 ) &
                    write(stdout,*) ' LINE --> CG coeff too large. Resetting :',gcfac
               gcfac = 0.0_dp
               ncg = 0
            else
               ncg = ncg + 1
            endif
         else
            gcfac = 0.0_dp
            ncg   = 0
         endif
      endif

      ! save for next iteration
      gcnorm0 = gcnorm1

      ! calculate search direction
      cdq(:,:,:) = cdodq(:,:,:) + cdqkeep(:,:,:) * gcfac

      ! calculate gradient along search direction - Tr[gradient . search direction]
      ! NB gradient is anti-hermitian
      doda0 = -real(zdotc(num_kpts*num_wann*num_wann,cdodq,1,cdq,1),dp)

      doda0 = doda0 / (4.0_dp*wbtot)

      ! check search direction is not uphill
      if (doda0.gt.0.0_dp) then
         ! if doing a CG step then reset CG
         if (ncg.gt.0) then
            if ( lprint .and. iprint>2 ) &
                 write(stdout,*) ' LINE --> Search direction uphill: resetting CG'
            cdq(:,:,:) = cdodq(:,:,:)
            ncg = 0
            gcfac = 0.0_dp
            ! re-calculate gradient along search direction
            doda0 = -real(zdotc(num_kpts*num_wann*num_wann,cdodq,1,cdq,1),dp)
            doda0 = doda0 / (4.0_dp*wbtot)
            ! if search direction still uphill then reverse search direction
            if (doda0.gt.0.0_dp) then
               if ( lprint .and. iprint>2 ) &
                    write(stdout,*) ' LINE --> Search direction still uphill: reversing'
               cdq(:,:,:) = -cdq(:,:,:)
               doda0 = -doda0
            endif
            ! if doing a SD step then reverse search direction
         else
            if ( lprint .and. iprint>2 ) &
                 write(stdout,*) ' LINE --> Search direction uphill: reversing'
            cdq(:,:,:) = -cdq(:,:,:)
            doda0 = -doda0
         endif
      endif

      ! calculate search direction
      cdq(:,:,:) = cdodq(:,:,:) + cdqkeep(:,:,:) * gcfac

      return

    end subroutine internal_search_direction


    !===============================================!
    subroutine internal_optimal_step()
      !===============================================!
      !                                               !
      ! Calculate the optimal step length based on a  !
      ! parabolic line search                         !
      !                                               !
      !===============================================!

      implicit none

      real(kind=dp) :: fac,shift,eqa,eqb


      fac = trial_spread%om_tot - wann_spread%om_tot
      if ( abs(fac) .gt. tiny(1.0_dp) ) then
         fac   = 1.0_dp/fac
         shift = 1.0_dp
      else
         fac    = 1.0e6_dp
         shift = fac*trial_spread%om_tot - fac*wann_spread%om_tot
      endif
      eqb = fac*doda0  
      eqa = shift - eqb*trial_step
      if ( abs(eqa/(fac*wann_spread%om_tot)).gt.epsilon(1.0_dp) ) then
         lquad=.true.
         alphamin  = - 0.5_dp * eqb / eqa * (trial_step**2)
         falphamin = wann_spread%om_tot &
              - 0.25_dp * eqb * eqb / (fac * eqa) * (trial_step**2)
      else
         if ( lprint .and. iprint>2 ) write(stdout,*) &
              ' LINE --> Parabolic line search unstable: using trial step'
         lquad=.false.
         alphamin  = trial_step
         falphamin = trial_spread%om_tot
      endif

      if (doda0*alphamin.gt.0.0_dp) then
         if ( lprint .and. iprint>2 ) write(stdout,*) &
              ' LINE --> Line search unstable : using trial step'
         lquad=.false.
         alphamin=trial_step
         falphamin=trial_spread%om_tot
      endif


      return

    end subroutine internal_optimal_step


    !===============================================!
    subroutine internal_new_u_and_m()               
      !===============================================!
      !                                               !
      ! Update U and M matrices after a trial step    !
      !                                               !
      !===============================================!

      implicit none

      integer :: nkp,nn,nkp2,nsdim
      logical :: ltmp

      call io_stopwatch('wann_main: u and m',1)

      do nkp=1,num_kpts
         ! cdq(nkp) is anti-Hermitian; tmp_cdq = i*cdq  is Hermitian
         tmp_cdq(:,:) = cmplx_i * cdq(:,:,nkp)
         ! Hermitian matrix eigen-solver
         call zheev('V','U',num_wann,tmp_cdq,num_wann,evals,cwork,4*num_wann,rwork,info)
         if (info.ne.0) then  
            write(stdout,*) &
                 'wann_main: ZHEEV in internal_new_u_and_m failed, info= ',info
            write(stdout,*) '           trying Schur decomposition instead'
!!$            call io_error('wann_main: problem in ZHEEV in internal_new_u_and_m') 
            tmp_cdq(:,:) = cdq(:,:,nkp)
            call zgees ('V', 'N', ltmp, num_wann, tmp_cdq, num_wann, nsdim, &
                 cwschur1, cz, num_wann, cwschur2, 10 * num_wann, cwschur3, &
                 cwschur4, info)
            if (info.ne.0) then  
               write(stdout,*) 'wann_main: SCHUR failed, info= ', info  
               call io_error('wann_main: problem computing schur form 1') 
            endif
            do i=1,num_wann
               tmp_cdq(:,i) = cz(:,i) * exp(cwschur1(i))
            enddo
            ! cmtmp   = tmp_cdq . cz^{dagger}
            call utility_zgemm(cmtmp,tmp_cdq,'N',cz,'C',num_wann)
            cdq(:,:,nkp)=cmtmp(:,:)
         else
            do i=1,num_wann
               cmtmp(:,i) = tmp_cdq(:,i) * exp(-cmplx_i * evals(i))
            enddo
            ! cdq(nkp)   = cmtmp . tmp_cdq^{dagger}
            call utility_zgemm(cdq(:,:,nkp),cmtmp,'N',tmp_cdq,'C',num_wann)
         endif
      enddo

!!$      do nkp = 1, num_kpts  
!!$         tmp_cdq(:,:) = cdq(:,:,nkp)
!!$         call zgees ('V', 'N', ltmp, num_wann, tmp_cdq, num_wann, nsdim, &
!!$              cwschur1, cz, num_wann, cwschur2, 10 * num_wann, cwschur3, &
!!$              cwschur4, info)
!!$         if (info.ne.0) then  
!!$            write(stdout,*) 'SCHUR: ', info  
!!$            call io_error('wann_main: problem computing schur form 1') 
!!$         endif
!!$         do i=1,num_wann
!!$            tmp_cdq(:,i) = cz(:,i) * exp(cwschur1(i))
!!$         enddo
!!$         ! cmtmp   = tmp_cdq . cz^{dagger}
!!$         call utility_zgemm(cmtmp,tmp_cdq,'N',cz,'C',num_wann)
!!$         cdq(:,:,nkp)=cmtmp(:,:)
!!$      enddo

      ! the orbitals are rotated
      do nkp=1,num_kpts
         ! cmtmp = U(k) . cdq(k)
         call utility_zgemm(cmtmp,u_matrix(:,:,nkp),'N',cdq(:,:,nkp),'N',num_wann)
         u_matrix(:,:,nkp)=cmtmp(:,:)
      enddo

      ! and the M_ij are updated
      do nkp = 1, num_kpts  
         do nn = 1, nntot  
            nkp2 = nnlist (nkp, nn)  
            ! tmp_cdq = cdq^{dagger} . M
            call utility_zgemm(tmp_cdq,cdq(:,:,nkp),'C',m_matrix(:,:,nn,nkp),'N',num_wann)
            ! cmtmp = tmp_cdq . cdq
            call utility_zgemm(cmtmp,tmp_cdq,'N',cdq(:,:,nkp2),'N',num_wann)
            m_matrix(:,:,nn,nkp) = cmtmp(:,:)
         enddo
      enddo

      call io_stopwatch('wann_main: u and m',2)

      return

    end subroutine internal_new_u_and_m


    !========================================!
    subroutine internal_check_unitarity()
    !========================================!

      implicit none

      integer :: nkp,i,j,m
      complex(kind=dp) :: ctmp1,ctmp2

      do nkp = 1, num_kpts  
         do i = 1, num_wann  
            do j = 1, num_wann  
               ctmp1 = cmplx_0  
               ctmp2 = cmplx_0  
               do m = 1, num_wann  
                  ctmp1 = ctmp1 + u_matrix (i, m, nkp) * conjg (u_matrix (j, m, nkp) )  
                  ctmp2 = ctmp2 + u_matrix (m, j, nkp) * conjg (u_matrix (m, i, nkp) )  
               enddo
               if ( (i.eq.j) .and. (abs (ctmp1 - cmplx_1 ) .gt.0.00001_dp) ) &
                    then
                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                       ctmp1
                  call io_error('wann_main: unitariety error 1')  
               endif
               if ( (i.eq.j) .and. (abs (ctmp2 - cmplx_1 ) .gt.0.00001_dp) ) &
                    then
                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                       ctmp2
                  call io_error('wann_main: unitariety error 2')  
               endif
               if ( (i.ne.j) .and. (abs (ctmp1) .gt.0.00001_dp) ) then  
                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                       ctmp1
                  call io_error('wann_main: unitariety error 3')  
               endif
               if ( (i.ne.j) .and. (abs (ctmp2) .gt.0.00001_dp) ) then  
                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                       ctmp2
                  call io_error('wann_main: unitariety error 4')  
               endif
            enddo
         enddo
      enddo


      return

    end subroutine internal_check_unitarity


    !========================================!
    subroutine internal_write_r2mn()
    !========================================!
    !                                        !
    ! Write seedname.r2mn file               !
    !                                        !
    !========================================!
      use w90_io, only: seedname,io_file_unit,io_error
      
      implicit none

      integer :: r2mnunit,nw1,nw2,nkp,nn
      real(kind=dp) :: r2ave_mn,delta

      ! note that here I use formulas analogue to Eq. 23, and not to the
      ! shift-invariant Eq. 32 .
      r2mnunit=io_file_unit()
      open(r2mnunit,file=trim(seedname)//'.r2mn',form='formatted',err=158)
      do nw1 = 1, num_wann  
         do nw2 = 1, num_wann  
            r2ave_mn = 0.0_dp  
            delta = 0.0_dp  
            if (nw1.eq.nw2) delta = 1.0_dp  
            do nkp = 1, num_kpts  
               do nn = 1, nntot
                  r2ave_mn = r2ave_mn + wb(nn) * &
                       ( 2.0_dp * delta - real(m_matrix(nw1,nw2,nn,nkp) - &
                       conjg(m_matrix(nw2,nw1,nn,nkp)),kind=dp) )
               enddo
            enddo
            r2ave_mn = r2ave_mn / real(num_kpts,dp)  
            write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn  
         enddo
      enddo
      close(r2mnunit)  
      
      return
      
158   call io_error('Error opening file '//trim(seedname)//'.r2mn in wann_main')
      
    end subroutine internal_write_r2mn



    !========================================!
    subroutine internal_svd_omega_i()
    !========================================!

      implicit none

      complex(kind=dp), allocatable  :: cv1(:,:),cv2(:,:)
      complex(kind=dp), allocatable  :: cw1(:),cw2(:)  
      complex(kind=dp), allocatable  :: cpad1 (:)  
      real(kind=dp),    allocatable  :: singvd (:)  

      integer :: nkp,nn,nb,na,ind
      real(kind=dp) :: omt1,omt2,omt3

      allocate( cw1 (10 * num_wann),stat=ierr  )
      if (ierr/=0) call io_error('Error in allocating cw1 in wann_main')
      allocate( cw2 (10 * num_wann),stat=ierr  )
      if (ierr/=0) call io_error('Error in allocating cw2 in wann_main')
      allocate( cv1 (num_wann, num_wann),stat=ierr  )
      if (ierr/=0) call io_error('Error in allocating cv1 in wann_main')
      allocate( cv2 (num_wann, num_wann),stat=ierr  )
      if (ierr/=0) call io_error('Error in allocating cv2 in wann_main')
      allocate( singvd (num_wann),stat=ierr  )
      if (ierr/=0) call io_error('Error in allocating singvd in wann_main')
      allocate( cpad1 (num_wann * num_wann),stat=ierr  )
      if (ierr/=0) call io_error('Error in allocating cpad1 in wann_main')

      cw1=cmplx_0; cw2=cmplx_0; cv1=cmplx_0; cv2=cmplx_0; cpad1=cmplx_0 
      singvd=0.0_dp

      ! singular value decomposition
      omt1 = 0.0_dp ; omt2 = 0.0_dp ; omt3 = 0.0_dp
      do nkp = 1, num_kpts  
         do nn = 1, nntot  
            ind = 1  
            do nb = 1, num_wann  
               do na = 1, num_wann  
                  cpad1 (ind) = m_matrix (na, nb, nn, nkp)  
                  ind = ind+1  
               enddo
            enddo
            call zgesvd ('A', 'A', num_wann, num_wann, cpad1, num_wann, singvd, cv1, &
                 num_wann, cv2, num_wann, cw1, 10 * num_wann, cw2, info)
            if (info.ne.0) then  
               call io_error('ERROR: Singular value decomp. zgesvd failed')  
            endif

            do nb = 1, num_wann  
               omt1 = omt1 + wb(nn) * (1.0_dp - singvd (nb) **2)  
               omt2 = omt2 - wb(nn) * (2.0_dp * log (singvd (nb) ) )  
               omt3 = omt3 + wb(nn) * (acos (singvd (nb) ) **2)  
            enddo
         enddo
      enddo
      omt1 = omt1 / real(num_kpts,dp)  
      omt2 = omt2 / real(num_kpts,dp)  
      omt3 = omt3 / real(num_kpts,dp)  
      write ( stdout , * ) ' '  
      write(stdout,'(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ',&
           omt1*lenconfac**2,'('//trim(length_unit)//'^2)'
      write(stdout,'(2x,a,f15.9,1x,a)') '                 -2log s = ',&
           omt2*lenconfac**2,'('//trim(length_unit)//'^2)'
      write(stdout,'(2x,a,f15.9,1x,a)') '                  acos^2 = ',&
           omt3*lenconfac**2,'('//trim(length_unit)//'^2)'

      deallocate(cpad1,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating cpad1 in wann_main')
      deallocate(singvd,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating singvd in wann_main')
      deallocate(cv2,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating cv2 in wann_main')
      deallocate(cv1,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating cv1 in wann_main')
      deallocate(cw2,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating cw2 in wann_main')
      deallocate(cw1,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating cw1 in wann_main')

      return
      
    end subroutine internal_svd_omega_i


  end subroutine wann_main


  !==================================================================!
  subroutine phases (csheet, sheet, rguide, irguide)
    !==================================================================!
    !                                                                  !
    !                                                                  !
    !===================================================================  
    use w90_parameters,     only : num_wann,m_matrix,nntot,neigh, &
         nnh,bk,bka,num_kpts
    use w90_io,         only : io_error,io_stopwatch
    use w90_utility,    only : utility_inv3

    implicit none

    complex(kind=dp), intent(out)   :: csheet (:,:,:)
    real(kind=dp)   , intent(out)   :: sheet (:,:,:)
    real(kind=dp)   , intent(inout) :: rguide (:,:)
    integer         , intent(in)    :: irguide

    !local
    complex(kind=dp) :: csum (nntot/2)  
    real(kind=dp)    ::  xx(nntot/2)
    real(kind=dp)    :: smat(3,3),svec(3),sinv(3,3)  
    real(kind=dp)    :: xx0,det,brn
    complex(kind=dp) :: csumt
    integer :: loop_wann,na,nkp,i,j,nn,ind,m


    call io_stopwatch('phases',1)


    csum=cmplx_0; xx=0.0_dp

    ! report problem to solve
    ! for each band, csum is determined and then its appropriate
    ! guiding center rguide(3,nwann)

    do loop_wann = 1,num_wann  
       ! get average phase for each unique bk direction
       do na = 1, nnh  
          csum (na) = cmplx_0
          do nkp = 1, num_kpts  
             nn = neigh (nkp, na)  
             csum (na) = csum (na) + m_matrix (loop_wann, loop_wann, nn, nkp)  
          enddo
       enddo

       ! now analyze that information to get good guess at
       ! wannier center
       !      write(*,*)
       !      do na=1,nnh
       !       write(*,'a,3f10.5,a,2f10.5)')
       !    &    ' bka=',(bka(j,na),j=1,3),'  csum=',csum(na)
       !      end do
       ! problem is to find a real-space 3-vector rguide such that
       !   phase of csum(nn) ~= phase of exp[ -i bka(nn) dot rguide ]
       ! or, letting
       !   xx(nn) = - Im ln csum(nn)  (modulo 2*pi)
       ! then
       !   bka(nn) dot rguide ~= xx(nn)
       !
       ! we take an arbitrary branch cut for first three xx(nn)
       ! and determine rguide from these; then for each additional bka
       ! vector, we first determine the most consistent branch cut,
       ! and then update rguide
       !
       ! rguide is obtained by minimizing
       !   sum_nn [ bka(nn) dot rguide - xx(nn) ] ^2
       ! or, setting the derivative with respect to rcenter to zero,
       !   sum_i smat(j,i) * rguide(i,nwann) = svec(j)
       ! where
       !   smat(j,i) = sum_nn bka(j,nn) * bka(i,nn)
       !   svec(j)   = sum_nn bka(j,nn) * xx(nn)
       ! initialize smat and svec

       smat=0.0_dp
       svec=0.0_dp

       do nn = 1, nnh  
          if (nn.le.3) then  
             !         obtain xx with arbitrary branch cut choice
             xx (nn) = - aimag (log (csum (nn) ) )  
          else  
             !         obtain xx with branch cut choice guided by rguide
             xx0 = 0.0_dp  
             do j = 1, 3  
                xx0 = xx0 + bka (j, nn) * rguide (j, loop_wann)  
             enddo
             !         xx0 is expected value for xx
!             csumt = exp (ci * xx0)  
             csumt = exp (cmplx_i * xx0)  
             !         csumt has opposite of expected phase of csum(nn)
             xx (nn) = xx0 - aimag (log (csum (nn) * csumt) )  
          endif

          !       write(*,'(a,i5,3f7.3,2f10.5)') 'nn, bka, xx, mag =',
          !    1    nn,(bka(j,nn),j=1,3),xx(nn),abs(csum(nn))/float(num_kpts)
          !       update smat and svec
          do j = 1, 3  
             do i = 1, 3  
                smat (j, i) = smat (j, i) + bka (j, nn) * bka (i, nn)  
             enddo
             svec (j) = svec (j) + bka (j, nn) * xx (nn)  
          enddo

          if (nn.ge.3) then  
             !         determine rguide
             call utility_inv3 (smat, sinv, det)  
             !         the inverse of smat is sinv/det
             if (abs (det) .gt.1.e-06_dp) then  
                !          to check that the first nn bka vectors are not
                !          linearly dependent - this is a change from original code
                if (irguide.ne.0) then  
                   do j = 1, 3  
                      rguide (j, loop_wann) = 0.0_dp  
                      do i = 1, 3  
                         rguide (j, loop_wann) = rguide (j, loop_wann) + sinv (j, i) &
                              * svec (i) / det
                      enddo
                   enddo
                endif
             endif
          endif

       enddo

    enddo

    !     obtain branch cut choice guided by rguid
    sheet = 0.0_dp
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          do loop_wann = 1, num_wann  
             ! sheet (loop_wann, nn, nkp) = 0.d0
             do j = 1, 3  
                sheet(loop_wann,nn,nkp) = sheet(loop_wann,nn,nkp) &
                     + bk(j,nn,nkp) * rguide(j,loop_wann)
             enddo
             ! csheet (loop_wann, nn, nkp) = exp (ci * sheet (loop_wann, nn, nkp) )  
          enddo
       enddo
    enddo
    csheet = exp(cmplx_i*sheet)

    ! now check that we picked the proper sheet for the log
    ! of m_matrix. criterion: q_n^{k,b}=Im(ln(M_nn^{k,b})) + b \cdot r_n are
    ! circa 0 for a good solution, circa multiples of 2 pi  for a bad one.
    ! I use the guiding center, instead of r_n, to understand which could be
    ! right sheet

    rnkb = 0.0_dp
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          do m = 1, num_wann  
             !           rnkb (m, nn, nkp) = 0.0_dp  
             brn = 0.0_dp
             do ind = 1, 3  
                brn = brn + bk(ind,nn,nkp) * rguide(ind,m)  
             enddo
             rnkb (m, nn, nkp) = rnkb (m, nn, nkp) + brn  
          enddo
       enddo
    enddo
!    write ( stdout , * ) ' '  
!    write ( stdout , * ) ' PHASES ARE SET USING THE GUIDING CENTERS'  
!    write ( stdout , * ) ' '  
!    do nkp = 1, num_kpts  
!       do n = 1, num_wann  
!          do nn = 1, nntot  
!             pherr = aimag(log(csheet(n,nn,nkp)*m_matrix(n,n,nn,nkp))) &
!                  - sheet(n,nn,nkp)+rnkb(n,nn,nkp)-aimag(log(m_matrix(n,n,nn,nkp)))
!          enddo
!       enddo
!    enddo

    call io_stopwatch('phases',2)

    return  

  end subroutine phases



  !==================================================================!
  subroutine omega(csheet,sheet,rave,r2ave,rave2,wann_spread)
    !==================================================================!
    !                                                                  !
    !   Calculate the Wannier Function spread                          !
    !                                                                  !
    !===================================================================  
    use w90_parameters, only : num_wann,m_matrix,nntot,wb,bk,num_kpts,&
                           omega_invariant
    use w90_io,         only : io_error,io_stopwatch

    implicit none

    complex(kind=dp), intent(in)  :: csheet (:,:,:)
    real(kind=dp)   , intent(in)  :: sheet (:,:,:)
    real(kind=dp)   , intent(out) :: rave (:,:)
    real(kind=dp)   , intent(out) :: r2ave (:)
    real(kind=dp)   , intent(out) :: rave2 (:)
    type(localisation_vars)    , intent(out)  :: wann_spread

    !local variables
    real(kind=dp) :: summ,mnn2
    real(kind=dp) :: brn
    integer :: ind,nkp,nn,m,n,iw
    logical, save :: first_pass=.true.

    call io_stopwatch('omega',1)


    do nkp = 1, num_kpts
       do nn = 1, nntot
          do n = 1, num_wann
             ln_tmp(n,nn,nkp)=( aimag(log(csheet(n,nn,nkp) &
                     * m_matrix(n,n,nn,nkp))) - sheet(n,nn,nkp) )
          end do
      end do
    end do


    rave  = 0.0_dp
    do iw = 1, num_wann  
       do ind = 1, 3  
          do nkp = 1, num_kpts  
             do nn = 1, nntot  
                rave(ind,iw) = rave(ind,iw) + wb(nn) * bk(ind,nn,nkp) &
                      *ln_tmp(iw,nn,nkp)
             enddo
          enddo
       enddo
    enddo
    rave = -rave/real(num_kpts,dp)

    rave2 = 0.0_dp
    do iw = 1, num_wann  
       rave2(iw) = sum(rave(:,iw)*rave(:,iw))
    enddo

    ! aam: is this useful?
!!$    rtot=0.0_dp
!!$    do ind = 1, 3  
!!$       do loop_wann = 1, num_wann  
!!$          rtot (ind) = rtot (ind) + rave (ind, loop_wann)  
!!$       enddo
!!$    enddo

    r2ave = 0.0_dp
    do iw = 1, num_wann  
       do nkp = 1, num_kpts  
          do nn = 1, nntot  
             mnn2 = real(m_matrix(iw,iw,nn,nkp)*conjg(m_matrix(iw,iw,nn,nkp)),kind=dp)
             r2ave(iw) = r2ave(iw) + wb(nn) * ( 1.0_dp - mnn2 + ln_tmp(iw,nn,nkp)**2 )
          enddo
       enddo
    enddo
    r2ave = r2ave/real(num_kpts,dp)

!!$    wann_spread%om_1 = 0.0_dp  
!!$    do nkp = 1, num_kpts  
!!$       do nn = 1, nntot  
!!$          do loop_wann = 1, num_wann  
!!$             wann_spread%om_1 = wann_spread%om_1 + wb(nn) * &
!!$                  ( 1.0_dp - m_matrix(loop_wann,loop_wann,nn,nkp) * &
!!$                  conjg(m_matrix(loop_wann,loop_wann,nn,nkp)) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$    wann_spread%om_1 = wann_spread%om_1 / real(num_kpts,dp)  
!!$
!!$    wann_spread%om_2 = 0.0_dp  
!!$    do loop_wann = 1, num_wann  
!!$       sqim = 0.0_dp  
!!$       do nkp = 1, num_kpts  
!!$          do nn = 1, nntot  
!!$             sqim = sqim + wb(nn) * &
!!$                  ( (aimag(log(csheet(loop_wann,nn,nkp) * &
!!$                  m_matrix(loop_wann,loop_wann,nn,nkp))) - &
!!$                  sheet(loop_wann,nn,nkp))**2 )
!!$          enddo
!!$       enddo
!!$       sqim = sqim / real(num_kpts,dp)  
!!$       wann_spread%om_2 = wann_spread%om_2 + sqim  
!!$    enddo
!!$
!!$    wann_spread%om_3 = 0.0_dp  
!!$    do loop_wann = 1, num_wann  
!!$       bim = 0.0_dp
!!$       do ind = 1, 3  
!!$          do nkp = 1, num_kpts  
!!$             do nn = 1, nntot  
!!$                bim(ind) = bim(ind) &
!!$                     + wb(nn) * bk(ind,nn,nkp) &
!!$                     * ( aimag(log(csheet(loop_wann,nn,nkp) &
!!$                     * m_matrix(loop_wann,loop_wann,nn,nkp))) &
!!$                     - sheet(loop_wann,nn,nkp) )
!!$             enddo
!!$          enddo
!!$       enddo
!!$       bim = bim/real(num_kpts,dp)
!!$       bim2 = 0.0_dp  
!!$       do ind = 1, 3  
!!$          bim2 = bim2 + bim (ind) * bim (ind)  
!!$       enddo
!!$       wann_spread%om_3 = wann_spread%om_3 - bim2  
!!$    enddo

    !jry: Either the above (om1,2,3) or the following is redundant
    !     keep it in the code base for testing

    ! wann_spread%om_i only needs to be calculated on the first pass
    ! on subsequent passes it may be set to omega_invariant
    if (first_pass) then
       wann_spread%om_i = 0.0_dp  
       do nkp = 1, num_kpts  
          do nn = 1, nntot  
             summ = 0.0_dp  
             do m = 1, num_wann  
                do n = 1, num_wann  
                   summ = summ &
                        + real(m_matrix(n,m,nn,nkp)*conjg(m_matrix(n,m,nn,nkp)),kind=dp)
                enddo
             enddo
             wann_spread%om_i = wann_spread%om_i &
                  + wb(nn) * (real(num_wann,dp) - summ)
          enddo
       enddo
       wann_spread%om_i = wann_spread%om_i / real(num_kpts,dp)
       first_pass=.false.
    else 
       wann_spread%om_i=omega_invariant
    endif

    wann_spread%om_od = 0.0_dp  
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          do m = 1, num_wann  
             do n = 1, num_wann  
                if (m.ne.n) wann_spread%om_od = wann_spread%om_od &
                     + wb(nn) * real( m_matrix(n,m,nn,nkp) &
                     * conjg(m_matrix(n,m,nn,nkp)), kind=dp )
             enddo
          enddo
       enddo
    enddo
    wann_spread%om_od = wann_spread%om_od / real(num_kpts,dp)  


    wann_spread%om_d = 0.0_dp  
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          do n = 1, num_wann  
             brn = sum(bk(:,nn,nkp)*rave(:,n))
             wann_spread%om_d = wann_spread%om_d + wb(nn) &
                  * ( ln_tmp(n,nn,nkp) + brn)**2
          enddo
       enddo
    enddo
    wann_spread%om_d = wann_spread%om_d / real(num_kpts,dp)  

    wann_spread%om_tot = wann_spread%om_i + wann_spread%om_d + wann_spread%om_od

    call io_stopwatch('omega',2)

    return  


  end subroutine omega


  !==================================================================!
  subroutine domega(csheet,sheet,rave,cdodq1,cdodq2,cdodq3,cdodq)
    !==================================================================!
    !                                                                  !
    !   Calculate the Gradient of the Wannier Function spread          !
    !                                                                  !
    !===================================================================  
    use w90_parameters, only : num_wann,wb,bk,nntot,m_matrix,num_kpts
    use w90_io,         only : io_error,io_stopwatch

    implicit none

    complex(kind=dp), intent(in)  :: csheet (:,:,:)    
    complex(kind=dp), intent(out) :: cdodq (:,:,:)     
    complex(kind=dp), intent(out) :: cdodq1 (:,:,:)    
    complex(kind=dp), intent(out) :: cdodq2 (:,:,:)    
    complex(kind=dp), intent(out) :: cdodq3 (:,:,:)    
    real(kind=dp),    intent(in)  :: sheet (:,:,:)     
    real(kind=dp),    intent(out) :: rave (:,:)        

    ! local
    integer :: iw,ind,nkp,nn,m,n
    complex(kind=dp) :: mnn

    call io_stopwatch('domega',1)

    do nkp = 1, num_kpts
       do nn = 1, nntot
          do n = 1, num_wann
             ln_tmp(n,nn,nkp)=wb(nn)*( aimag(log(csheet(n,nn,nkp) &
                     * m_matrix(n,n,nn,nkp))) - sheet(n,nn,nkp) )
          end do
      end do
    end do

    ! recalculate rave
    rave = 0.0_dp
    do iw = 1, num_wann  
       do ind = 1, 3  
          do nkp = 1, num_kpts  
             do nn = 1, nntot  
                rave(ind,iw) = rave(ind,iw) +  bk(ind,nn,nkp) &
                     * ln_tmp(iw,nn,nkp)
             enddo
          enddo
       enddo
    enddo
    rave = -rave/real(num_kpts,dp)

    ! R_mn=M_mn/M_nn and q_m^{k,b} = Im phi_m^{k,b} + b.r_n are calculated
    rnkb = 0.0_dp
    do nkp=1,num_kpts
       do nn=1,nntot
          do n=1,num_wann
             mnn = m_matrix(n,n,nn,nkp)
             crt(:,n,nn,nkp) = m_matrix(:,n,nn,nkp) / mnn
             cr(:,n,nn,nkp)  = m_matrix(:,n,nn,nkp) * conjg(mnn)
             rnkb(n,nn,nkp) = sum(bk(:,nn,nkp)*rave(:,n))
          enddo
       enddo
    enddo

    ! cd0dq(m,n,nkp) is calculated
    cdodq1 = cmplx_0 ; cdodq2 = cmplx_0 ; cdodq3 = cmplx_0
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          do n = 1, num_wann  
             do m = 1, num_wann  
                ! A[R^{k,b}]=(R-Rdag)/2
                cdodq1(m,n,nkp) = cdodq1(m,n,nkp) &
                     + wb(nn) * cmplx(0.5_dp,0.0_dp,kind=dp) &
                     *( cr(m,n,nn,nkp) - conjg(cr(n,m,nn,nkp)) )
                ! -S[T^{k,b}]=-(T+Tdag)/2i ; T_mn = Rt_mn q_n
                cdodq2(m,n,nkp) = cdodq2(m,n,nkp) -  &
                      ( crt(m,n,nn,nkp) * ln_tmp(n,nn,nkp)  &
                     + conjg( crt(n,m,nn,nkp) * ln_tmp(m,nn,nkp) ) ) &
                     * cmplx(0.0_dp,-0.5_dp,kind=dp)
                cdodq3(m,n,nkp) = cdodq3(m,n,nkp) - wb(nn) &
                     * ( crt(m,n,nn,nkp) * rnkb(n,nn,nkp) + conjg(crt(n,m,nn,nkp) &
                     * rnkb(m,nn,nkp)) ) * cmplx(0.0_dp,-0.5_dp,kind=dp)
             enddo
          enddo
       enddo
    enddo
    cdodq1 = cdodq1 / cmplx(num_kpts,0.0_dp,kind=dp) * cmplx(4.0_dp,0.0_dp,kind=dp)
    cdodq2 = cdodq2 / cmplx(num_kpts,0.0_dp,kind=dp) * cmplx(4.0_dp,0.0_dp,kind=dp)
    cdodq3 = cdodq3 / cmplx(num_kpts,0.0_dp,kind=dp) * cmplx(4.0_dp,0.0_dp,kind=dp)
    cdodq  = cdodq1 + cdodq2 + cdodq3

    call io_stopwatch('domega',2)

    return  


  end subroutine domega


  !==================================================================!
  subroutine wann_spread_copy(orig,copy)
  !==================================================================!
  !                                                                  !
  !==================================================================!

    implicit none

    type(localisation_vars), intent(in)  :: orig
    type(localisation_vars), intent(out) :: copy

    copy%om_i   =  orig%om_i  
    copy%om_d   =  orig%om_d  
    copy%om_od  =  orig%om_od 
    copy%om_tot =  orig%om_tot
!!$    copy%om_1   =  orig%om_1  
!!$    copy%om_2   =  orig%om_2  
!!$    copy%om_3   =  orig%om_3  

    return

  end subroutine wann_spread_copy


  !==================================================================!
  subroutine wann_calc_projection()
  !==================================================================!
  !                                                                  !
  ! Calculates and writes the projection of each Wannier function    !
  ! on the original bands within the outer window.                   !
  !                                                                  !
  !==================================================================!

    use w90_parameters, only : num_bands,num_wann,num_kpts,&
                           u_matrix_opt,eigval,lwindow
    use w90_io,         only : stdout

    implicit none

    integer :: nw,nb,nkp,counter
    real(kind=dp) :: summ

    write(stdout,'(/1x,a78)') repeat('-',78)
    write(stdout,'(1x,9x,a)') &
         'Projection of Bands in Outer Window on all Wannier Functions'
    write(stdout,'(1x,8x,62a)') repeat('-',62)
    write(stdout,'(1x,16x,a)') '   Kpt  Band      Eigval      |Projection|^2'
    write(stdout,'(1x,16x,a47)') repeat('-',47) 

    do nkp=1,num_kpts
       counter=0
       do nb=1,num_bands
          if (lwindow(nb,nkp)) then
             counter=counter+1
             summ=0.0_dp
             do nw=1,num_wann
                summ=summ+abs(u_matrix_opt(counter,nw,nkp))**2
             enddo
             write(stdout,'(1x,16x,i5,1x,i5,1x,f14.6,2x,f14.8)') &
                  nkp,nb,eigval(nb,nkp),summ
          endif
       enddo
    enddo
    write(stdout,'(1x,a78/)') repeat('-',78)

    return

  end subroutine wann_calc_projection

end module w90_wannierise
