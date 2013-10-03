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

module w90_wannierise

  use w90_constants

  implicit none

  private

  public :: wann_main
  public :: wann_main_gamma  ![ysl]

  ! Data to avoid large allocation within iteration loop
  real(kind=dp),    allocatable  :: rnkb (:,:,:)   
  real(kind=dp),    allocatable  :: ln_tmp(:,:,:)

  ! The next variable is used to trigger the calculation of the invarient spread
  ! we only need to do this on entering wann_main (_gamma)
  logical :: first_pass


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
    use w90_io,         only : stdout,io_error,io_time,io_stopwatch &
         ,io_file_unit
    use w90_parameters, only : num_wann,num_cg_steps,num_iter,nnlist, &
         nntot,wbtot,u_matrix,m_matrix,num_kpts,iprint,num_print_cycles, &
         num_dump_cycles,omega_invariant,param_write_chkpt,length_unit, &
         lenconfac,proj_site,real_lattice,write_r2mn,guiding_centres, &
         num_guide_cycles,num_no_guide_iter,timing_level,trial_step,spinors, &
         fixed_step,lfixstep,write_proj,have_disentangled,conv_tol,num_proj, &
         conv_window,conv_noise_amp,conv_noise_num,wannier_centres,write_xyz, &
         wannier_spreads,omega_total,omega_tilde,optimisation,write_vdw_data,&
         write_hr_diag
    use w90_utility,    only : utility_frac_to_cart,utility_zgemm

    !ivo
    use w90_hamiltonian, only : hamiltonian_setup,hamiltonian_get_hr,ham_r,&
                                rpt_origin

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
    integer       :: i,n,iter,ind,ierr,iw,ncg,info
    logical       :: lprint,ldump,lquad
    real(kind=dp), allocatable :: history(:)
    real(kind=dp)              :: save_spread
    logical                    :: lconverged,lrandom,lfirst
    integer                    :: conv_count,noise_count,page_unit

    if (timing_level>0) call io_stopwatch('wann: main',1)

    first_pass=.true.

    ! Allocate stuff

    allocate(history(conv_window),stat=ierr)
    if (ierr/=0) call io_error('Error allocating history in wann_main')

    ! module data
    if(optimisation>0) then
       allocate(  m0 (num_wann, num_wann, nntot, num_kpts),stat=ierr)
    end if
    if (ierr/=0) call io_error('Error in allocating m0 in wann_main')
    allocate(  u0 (num_wann, num_wann, num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating u0 in wann_main')
    allocate( rnkb (num_wann, nntot, num_kpts),stat=ierr    )     
    if (ierr/=0) call io_error('Error in allocating rnkb in wann_main')
    allocate( ln_tmp (num_wann, nntot, num_kpts), stat=ierr    )
    if (ierr/=0) call io_error('Error in allocating ln_tmp in wann_main')

    rnkb=0.0_dp

    ! sub vars passed into other subs 
    allocate( csheet (num_wann, nntot, num_kpts), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating csheet in wann_main')
    allocate( cdodq (num_wann, num_wann, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cdodq in wann_main')
    allocate( sheet (num_wann, nntot, num_kpts), stat=ierr    )
    if (ierr/=0) call io_error('Error in allocating sheet in wann_main')
    allocate( rave (3, num_wann), stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating rave in wann_main')
    allocate( r2ave (num_wann),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating r2ave in wann_main')
    allocate( rave2 (num_wann),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating rave2 in wann_main')
    allocate( rguide (3, num_wann)   )
    if (ierr/=0) call io_error('Error in allocating rguide in wann_main')

    csheet=cmplx_1;cdodq=cmplx_0
    sheet=0.0_dp;rave=0.0_dp;r2ave=0.0_dp;rave2=0.0_dp;rguide=0.0_dp

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
       do n=1,num_proj
          call utility_frac_to_cart(proj_site(:,n),rguide(:,n),real_lattice)
       enddo
!       if(spinors) then ! not needed with new changes to spinor proj 2013 JRY
!          do n=1,num_proj
!             call utility_frac_to_cart(proj_site(:,n),rguide(:,n+num_proj),real_lattice)
!          enddo
!       end if
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
       call wann_phases(csheet,sheet,rguide,irguide)
       irguide=1
    endif

    ! calculate initial centers and spread
    call wann_omega(csheet,sheet,rave,r2ave,rave2,wann_spread)

    ! public variables
    omega_total = wann_spread%om_tot
    omega_invariant = wann_spread%om_i
    omega_tilde = wann_spread%om_d + wann_spread%om_od
  
    ! public arrays of Wannier centres and spreads
    wannier_centres = rave
    wannier_spreads = r2ave - rave2

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

    lconverged=.false. ; lfirst=.true. ; lrandom=.false.
    conv_count=0 ; noise_count=0

    if(.not.lfixstep .and.optimisation<=0) then
       page_unit=io_file_unit()
       open(unit=page_unit,status='scratch',form='unformatted')
    endif

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
          call wann_phases(csheet,sheet,rguide,irguide)
          irguide=1
       endif

       ! calculate gradient of omega
       call wann_domega(csheet,sheet,rave,cdodq)

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
          u0=u_matrix 

          if(optimisation<=0) then
             write(page_unit)   m_matrix
             rewind(page_unit)
          else
             m0=m_matrix
          endif

          ! update U and M
          call internal_new_u_and_m()

          ! calculate spread at trial step
          call wann_omega(csheet,sheet,rave,r2ave,rave2,trial_spread)

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
             u_matrix=u0
             if(optimisation<=0) then
                read(page_unit)  m_matrix
                rewind(page_unit)
             else
                m_matrix=m0
             endif
          endif

          ! update U and M
          call internal_new_u_and_m()
          
          call wann_spread_copy(wann_spread,old_spread)
          
          ! calculate the new centers and spread
          call wann_omega(csheet,sheet,rave,r2ave,rave2,wann_spread)
        
       ! parabolic line search was unsuccessful, use trial step already taken
       else 

          call wann_spread_copy(wann_spread,old_spread)
          call wann_spread_copy(trial_spread,wann_spread)

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

       ! Public array of Wannier centres and spreads
       wannier_centres = rave
       wannier_spreads = r2ave - rave2

       ! Public variables
       omega_total = wann_spread%om_tot
       omega_tilde = wann_spread%om_d + wann_spread%om_od

       if (ldump) call param_write_chkpt('postdis')

       if (conv_window.gt.1) call internal_test_convergence()

       if (lconverged) then 
          write(stdout,'(/13x,a,es10.3,a,i2,a)') &
               '<<<     Delta <',conv_tol,&
               '  over ',conv_window,' iterations     >>>'
          write(stdout,'(13x,a/)')  '<<< Wannierisation convergence criteria satisfied >>>'
          exit
       endif

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

    if (write_xyz) call wann_write_xyz()

    if(write_hr_diag) then
       call hamiltonian_setup()
       call hamiltonian_get_hr()
       write(stdout,*)
       write(stdout,'(1x,a)') 'On-site Hamiltonian matrix elements'
       write(stdout,'(3x,a)') '  n        <0n|H|0n> (eV)'
       write(stdout,'(3x,a)') '-------------------------'
       do i=1,num_wann
          write(stdout,'(3x,i3,5x,f12.6)') i,real(ham_r(i,i,rpt_origin),kind=dp)
       enddo
       write(stdout,*)
    endif

    if (guiding_centres) call wann_phases(csheet,sheet,rguide,irguide)

    ! unitarity is checked
!!$    call internal_check_unitarity()
    call wann_check_unitarity()

    ! write extra info regarding omega_invariant
!!$    if (iprint>2) call internal_svd_omega_i()
    if (iprint>2) call wann_svd_omega_i()

    ! write matrix elements <m|r^2|n> to file
!!$    if (write_r2mn) call internal_write_r2mn()
    if (write_r2mn) call wann_write_r2mn()

    ! calculate and write projection of WFs on original bands in outer window
    if (have_disentangled .and. write_proj) call wann_calc_projection()

    ! aam: write data required for vdW utility
    if (write_vdw_data) call wann_write_vdw_data()

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
    deallocate(cdodq,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cdodq in wann_main')
    deallocate(csheet,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating csheet in wann_main')
    
    ! deallocate module data
    deallocate( ln_tmp , stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating ln_tmp in wann_main')
    deallocate( rnkb,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating rnkb in wann_main')

    deallocate(u0, stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating u0 in wann_main')
    if(optimisation>0) then
       deallocate(m0, stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating m0 in wann_main')
    end if

    deallocate(history,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating history in wann_main')

    if (timing_level>0) call io_stopwatch('wann: main',2)

    return

1000 format(2x,'WF centre and spread', &
         &       i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)

1001 format(2x,'Sum of centres and spreads', &
         &       1x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)


  contains

    
    !===============================================!
    subroutine internal_test_convergence()
      !===============================================!
      !                                               !
      ! Determine whether minimisation of non-gauge   !
      ! invariant spread is converged                 !
      !                                               !
      !===============================================!

      implicit none

      real(kind=dp) :: delta_omega
      integer :: j,ierr
      real(kind=dp), allocatable :: temp_hist(:)
        
      allocate(temp_hist(conv_window),stat=ierr)
      if (ierr/=0) call io_error('Error allocating temp_hist in wann_main')

      delta_omega=wann_spread%om_tot-old_spread%om_tot
       
      if (iter.le.conv_window) then
         history(iter) = delta_omega
      else
         temp_hist = eoshift(history,1,delta_omega)
         history = temp_hist
      endif

      conv_count=conv_count+1

      if (conv_count.lt.conv_window) then
         return
      else
!!$         write(stdout,*) (history(j),j=1,conv_window)
         do j=1,conv_window
            if ( abs(history(j)).gt.conv_tol ) return
         enddo
      endif

      if ((conv_noise_amp.gt.0.0_dp) .and. (noise_count.lt.conv_noise_num)) then
         if (lfirst) then
            lfirst=.false.
            save_spread=wann_spread%om_tot
            lrandom=.true.
            conv_count=0
         else
            if (abs(save_spread-wann_spread%om_tot).lt.conv_tol) then
               lconverged=.true.
               return
            else
               save_spread=wann_spread%om_tot
               lrandom=.true.
               conv_count=0
            endif
         endif
      else
         lconverged=.true.
      endif

      if (lrandom) noise_count = noise_count + 1

      deallocate(temp_hist,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating temp_hist in wann_main')

      return

    end subroutine internal_test_convergence



    !===============================================!
    subroutine internal_random_noise()
      !===============================================!
      !                                               !
      ! Add some random noise to the search direction !
      ! to help escape from local minima              !
      !                                               !
      !===============================================!

      implicit none

      integer :: ikp,iw,jw
      real(kind=dp), allocatable :: noise_real(:,:), noise_imag(:,:)
      complex(kind=dp), allocatable :: cnoise(:,:)

      ! Allocate
      allocate(noise_real(num_wann,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error allocating noise_real in wann_main')
      allocate(noise_imag(num_wann,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error allocating noise_imag in wann_main')
      allocate(cnoise(num_wann,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error allocating cnoise in wann_main')
      
      ! Initialise
      cnoise=cmplx_0; noise_real=0.0_dp; noise_imag=0.0_dp

      ! cdq is a num_wann x num_wann x num_kpts anti-hermitian array
      ! to which we add a random anti-hermitian matrix 

      do ikp=1,num_kpts
         do iw=1,num_wann
            call random_seed()
            call random_number(noise_real(:,iw))
            call random_seed()
            call random_number(noise_imag(:,iw))
         enddo
         do jw=1,num_wann
            do iw=1,jw
               if (iw.eq.jw) then 
                  cnoise(iw,jw) = cmplx(0.0_dp,noise_imag(iw,jw),dp)
               else
                  cnoise(iw,jw) = cmplx(noise_real(iw,jw),noise_imag(iw,jw),dp)
               endif
               cnoise(jw,iw) = -conjg(cnoise(iw,jw))
            enddo
         enddo
         ! Add noise to search direction
         cdq(:,:,ikp) = cdq(:,:,ikp) + conv_noise_amp * cnoise(:,:)
      enddo

      ! Deallocate
      deallocate(cnoise,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating cnoise in wann_main')
      deallocate(noise_imag,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating noise_imag in wann_main')
      deallocate(noise_real,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating noise_real in wann_main')

      return

    end subroutine internal_random_noise



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

      if (timing_level>1) call io_stopwatch('wann: main: search_direction',1)

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

      ! add some random noise to search direction, if required
      if (lrandom) then
         write(stdout,'(a,i3,a,i3,a)') &
              ' [ Adding random noise to search direction. Time ',noise_count,' / ',conv_noise_num,' ]'
         call internal_random_noise()
      endif
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
            if (lrandom) call internal_random_noise()
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

 !!$     ! calculate search direction
 !!$     cdq(:,:,:) = cdodq(:,:,:) + cdqkeep(:,:,:) * gcfac

      if (timing_level>1) call io_stopwatch('wann: main: search_direction',2)

      lrandom=.false.

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

      if (timing_level>1) call io_stopwatch('wann: main: optimal_step',1)

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

      if (timing_level>1) call io_stopwatch('wann: main: optimal_step',2)

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

      if (timing_level>1) call io_stopwatch('wann: main: u_and_m',1)

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

      if (timing_level>1) call io_stopwatch('wann: main: u_and_m',2)

      return

    end subroutine internal_new_u_and_m


!!$    !========================================!
!!$    subroutine internal_check_unitarity()
!!$    !========================================!
!!$
!!$      implicit none
!!$
!!$      integer :: nkp,i,j,m
!!$      complex(kind=dp) :: ctmp1,ctmp2
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',1)
!!$
!!$      do nkp = 1, num_kpts  
!!$         do i = 1, num_wann  
!!$            do j = 1, num_wann  
!!$               ctmp1 = cmplx_0  
!!$               ctmp2 = cmplx_0  
!!$               do m = 1, num_wann  
!!$                  ctmp1 = ctmp1 + u_matrix (i, m, nkp) * conjg (u_matrix (j, m, nkp) )  
!!$                  ctmp2 = ctmp2 + u_matrix (m, j, nkp) * conjg (u_matrix (m, i, nkp) )  
!!$               enddo
!!$               if ( (i.eq.j) .and. (abs (ctmp1 - cmplx_1 ) .gt. eps5) ) &
!!$                    then
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp1
!!$                  call io_error('wann_main: unitariety error 1')  
!!$               endif
!!$               if ( (i.eq.j) .and. (abs (ctmp2 - cmplx_1 ) .gt. eps5) ) &
!!$                    then
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp2
!!$                  call io_error('wann_main: unitariety error 2')  
!!$               endif
!!$               if ( (i.ne.j) .and. (abs (ctmp1) .gt. eps5) ) then  
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp1
!!$                  call io_error('wann_main: unitariety error 3')  
!!$               endif
!!$               if ( (i.ne.j) .and. (abs (ctmp2) .gt. eps5) ) then  
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp2
!!$                  call io_error('wann_main: unitariety error 4')  
!!$               endif
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',2)
!!$
!!$      return
!!$
!!$    end subroutine internal_check_unitarity


!!$    !========================================!
!!$    subroutine internal_write_r2mn()
!!$    !========================================!
!!$    !                                        !
!!$    ! Write seedname.r2mn file               !
!!$    !                                        !
!!$    !========================================!
!!$      use w90_io, only: seedname,io_file_unit,io_error
!!$      
!!$      implicit none
!!$
!!$      integer :: r2mnunit,nw1,nw2,nkp,nn
!!$      real(kind=dp) :: r2ave_mn,delta
!!$
!!$      ! note that here I use formulas analogue to Eq. 23, and not to the
!!$      ! shift-invariant Eq. 32 .
!!$      r2mnunit=io_file_unit()
!!$      open(r2mnunit,file=trim(seedname)//'.r2mn',form='formatted',err=158)
!!$      do nw1 = 1, num_wann  
!!$         do nw2 = 1, num_wann  
!!$            r2ave_mn = 0.0_dp  
!!$            delta = 0.0_dp  
!!$            if (nw1.eq.nw2) delta = 1.0_dp  
!!$            do nkp = 1, num_kpts  
!!$               do nn = 1, nntot
!!$                  r2ave_mn = r2ave_mn + wb(nn) * &
!!$                       ! [GP-begin, Apr13, 2012: corrected sign inside "real"]
!!$                       ( 2.0_dp * delta - real(m_matrix(nw1,nw2,nn,nkp) + &
!!$                       conjg(m_matrix(nw2,nw1,nn,nkp)),kind=dp) )
!!$                       ! [GP-end]
!!$               enddo
!!$            enddo
!!$            r2ave_mn = r2ave_mn / real(num_kpts,dp)  
!!$            write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn  
!!$         enddo
!!$      enddo
!!$      close(r2mnunit)  
!!$      
!!$      return
!!$      
!!$158   call io_error('Error opening file '//trim(seedname)//'.r2mn in wann_main')
!!$      
!!$    end subroutine internal_write_r2mn



!!$    !========================================!
!!$    subroutine internal_svd_omega_i()
!!$    !========================================!
!!$
!!$      implicit none
!!$
!!$      complex(kind=dp), allocatable  :: cv1(:,:),cv2(:,:)
!!$      complex(kind=dp), allocatable  :: cw1(:),cw2(:)  
!!$      complex(kind=dp), allocatable  :: cpad1 (:)  
!!$      real(kind=dp),    allocatable  :: singvd (:)  
!!$
!!$      integer :: nkp,nn,nb,na,ind
!!$      real(kind=dp) :: omt1,omt2,omt3
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',1)
!!$
!!$      allocate( cw1 (10 * num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cw1 in wann_main')
!!$      allocate( cw2 (10 * num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cw2 in wann_main')
!!$      allocate( cv1 (num_wann, num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cv1 in wann_main')
!!$      allocate( cv2 (num_wann, num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cv2 in wann_main')
!!$      allocate( singvd (num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating singvd in wann_main')
!!$      allocate( cpad1 (num_wann * num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cpad1 in wann_main')
!!$
!!$      cw1=cmplx_0; cw2=cmplx_0; cv1=cmplx_0; cv2=cmplx_0; cpad1=cmplx_0 
!!$      singvd=0.0_dp
!!$
!!$      ! singular value decomposition
!!$      omt1 = 0.0_dp ; omt2 = 0.0_dp ; omt3 = 0.0_dp
!!$      do nkp = 1, num_kpts  
!!$         do nn = 1, nntot  
!!$            ind = 1  
!!$            do nb = 1, num_wann  
!!$               do na = 1, num_wann  
!!$                  cpad1 (ind) = m_matrix (na, nb, nn, nkp)  
!!$                  ind = ind+1  
!!$               enddo
!!$            enddo
!!$            call zgesvd ('A', 'A', num_wann, num_wann, cpad1, num_wann, singvd, cv1, &
!!$                 num_wann, cv2, num_wann, cw1, 10 * num_wann, cw2, info)
!!$            if (info.ne.0) then  
!!$               call io_error('ERROR: Singular value decomp. zgesvd failed')  
!!$            endif
!!$
!!$            do nb = 1, num_wann  
!!$               omt1 = omt1 + wb(nn) * (1.0_dp - singvd (nb) **2)  
!!$               omt2 = omt2 - wb(nn) * (2.0_dp * log (singvd (nb) ) )  
!!$               omt3 = omt3 + wb(nn) * (acos (singvd (nb) ) **2)  
!!$            enddo
!!$         enddo
!!$      enddo
!!$      omt1 = omt1 / real(num_kpts,dp)  
!!$      omt2 = omt2 / real(num_kpts,dp)  
!!$      omt3 = omt3 / real(num_kpts,dp)  
!!$      write ( stdout , * ) ' '  
!!$      write(stdout,'(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ',&
!!$           omt1*lenconfac**2,'('//trim(length_unit)//'^2)'
!!$      write(stdout,'(2x,a,f15.9,1x,a)') '                 -2log s = ',&
!!$           omt2*lenconfac**2,'('//trim(length_unit)//'^2)'
!!$      write(stdout,'(2x,a,f15.9,1x,a)') '                  acos^2 = ',&
!!$           omt3*lenconfac**2,'('//trim(length_unit)//'^2)'
!!$
!!$      deallocate(cpad1,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cpad1 in wann_main')
!!$      deallocate(singvd,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating singvd in wann_main')
!!$      deallocate(cv2,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cv2 in wann_main')
!!$      deallocate(cv1,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cv1 in wann_main')
!!$      deallocate(cw2,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cw2 in wann_main')
!!$      deallocate(cw1,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cw1 in wann_main')
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',2)
!!$
!!$      return
!!$      
!!$    end subroutine internal_svd_omega_i


  end subroutine wann_main


  !==================================================================!
  subroutine wann_phases (csheet, sheet, rguide, irguide, m_w)
    !==================================================================!
    !                                                                  !
    !                                                                  !
    !===================================================================  
    use w90_constants,  only : eps6
    use w90_parameters, only : num_wann,m_matrix,nntot,neigh, &
         nnh,bk,bka,num_kpts,timing_level
    use w90_io,         only : io_stopwatch
    use w90_utility,    only : utility_inv3

    implicit none

    complex(kind=dp), intent(out)   :: csheet (:,:,:)
    real(kind=dp)   , intent(out)   :: sheet (:,:,:)
    real(kind=dp)   , intent(inout) :: rguide (:,:)
    integer         , intent(in)    :: irguide
    real(kind=dp), intent(in), optional :: m_w(:,:,:)

    !local
    complex(kind=dp) :: csum (nnh)  
    real(kind=dp)    ::  xx(nnh)
    real(kind=dp)    :: smat(3,3),svec(3),sinv(3,3)  
    real(kind=dp)    :: xx0,det,brn
    complex(kind=dp) :: csumt
    integer :: loop_wann,na,nkp,i,j,nn,ind,m


    if (timing_level>1) call io_stopwatch('wann: phases',1)


    csum=cmplx_0; xx=0.0_dp

    ! report problem to solve
    ! for each band, csum is determined and then its appropriate
    ! guiding center rguide(3,nwann)

    do loop_wann = 1,num_wann  

       if (.not. present(m_w) ) then
          ! get average phase for each unique bk direction
          do na = 1, nnh  
             csum (na) = cmplx_0
             do nkp = 1, num_kpts  
                nn = neigh (nkp, na)  
                csum (na) = csum (na) + m_matrix (loop_wann, loop_wann, nn, nkp)  
             enddo
          enddo

       else
       
          do na = 1, nnh  
             csum (na) = cmplx_0
             do nkp = 1, num_kpts  
                nn = neigh (nkp, na)  
                csum (na) = csum (na) &
                + cmplx(m_w(loop_wann,loop_wann,2*nn-1),m_w(loop_wann,loop_wann,2*nn),dp)  
             enddo
          enddo

       end if
          
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
             if (abs (det) .gt.eps6) then  
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

    if (timing_level>1) call io_stopwatch('wann: phases',2)

    return  

  end subroutine wann_phases



  !==================================================================!
  subroutine wann_omega(csheet,sheet,rave,r2ave,rave2,wann_spread)
    !==================================================================!
    !                                                                  !
    !   Calculate the Wannier Function spread                          !
    !                                                                  !
    !===================================================================  
    use w90_parameters, only : num_wann,m_matrix,nntot,wb,bk,num_kpts,&
                           omega_invariant,timing_level
    use w90_io,         only : io_stopwatch

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

    if (timing_level>1) call io_stopwatch('wann: omega',1)


    do nkp = 1, num_kpts
       do nn = 1, nntot
          do n = 1, num_wann
             ! Note that this ln_tmp is defined differently wrt the one in wann_domega
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

    if (timing_level>1) call io_stopwatch('wann: omega',2)

    return  


  end subroutine wann_omega


  !==================================================================!
  subroutine wann_domega(csheet,sheet,rave,cdodq)
    !==================================================================!
    !                                                                  !
    !   Calculate the Gradient of the Wannier Function spread          !
    !                                                                  !
    !===================================================================  
    use w90_parameters, only : num_wann,wb,bk,nntot,m_matrix,num_kpts,timing_level
    use w90_io,         only : io_stopwatch,io_error

    implicit none

    complex(kind=dp), intent(in)  :: csheet (:,:,:)    
    real(kind=dp),    intent(in)  :: sheet (:,:,:)     
    real(kind=dp),    intent(out) :: rave (:,:)        
    complex(kind=dp), intent(out) :: cdodq (:,:,:)     

    !local
    complex(kind=dp), allocatable  :: cr (:,:)   
    complex(kind=dp), allocatable  :: crt (:,:)  

    ! local
    integer :: iw,ind,nkp,nn,m,n,ierr
    complex(kind=dp) :: mnn

    if (timing_level>1) call io_stopwatch('wann: domega',1)

    allocate(  cr (num_wann, num_wann),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cr in wann_domega')
    allocate(  crt (num_wann, num_wann),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating crt in wann_domega')


    do nkp = 1, num_kpts
       do nn = 1, nntot
          do n = 1, num_wann
             ! Note that this ln_tmp is defined differently wrt the one in wann_omega
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
             rnkb(n,nn,nkp) = sum(bk(:,nn,nkp)*rave(:,n))
          enddo
       enddo
    enddo

    ! cd0dq(m,n,nkp) is calculated
    cdodq=cmplx_0
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          do n=1,num_wann
             mnn = m_matrix(n,n,nn,nkp)
             crt(:,n) = m_matrix(:,n,nn,nkp) / mnn
             cr(:,n)  = m_matrix(:,n,nn,nkp) * conjg(mnn)
          enddo

          do n = 1, num_wann  
             do m = 1, num_wann  
                ! A[R^{k,b}]=(R-Rdag)/2
                cdodq(m,n,nkp) = cdodq(m,n,nkp) &
                     + wb(nn) * 0.5_dp &
                     *( cr(m,n) - conjg(cr(n,m)) )
                ! -S[T^{k,b}]=-(T+Tdag)/2i ; T_mn = Rt_mn q_n
                cdodq(m,n,nkp) = cdodq(m,n,nkp) -  &
                      ( crt(m,n) * ln_tmp(n,nn,nkp)  &
                     + conjg( crt(n,m) * ln_tmp(m,nn,nkp) ) ) &
                     * cmplx(0.0_dp,-0.5_dp,kind=dp)
                cdodq(m,n,nkp) = cdodq(m,n,nkp) - wb(nn) &
                     * ( crt(m,n) * rnkb(n,nn,nkp) + conjg(crt(n,m) &
                     * rnkb(m,nn,nkp)) ) * cmplx(0.0_dp,-0.5_dp,kind=dp)
             enddo
          enddo
       enddo
    enddo
    cdodq = cdodq / real(num_kpts,dp) * 4.0_dp


    if (timing_level>1) call io_stopwatch('wann: domega',2)

    return  


  end subroutine wann_domega


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
                           u_matrix_opt,eigval,lwindow,timing_level
    use w90_io,         only : stdout,io_stopwatch

    implicit none

    integer :: nw,nb,nkp,counter
    real(kind=dp) :: summ

    if (timing_level>1) call io_stopwatch('wann: calc_projection',1)

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

    if (timing_level>1) call io_stopwatch('wann: calc_projection',2)

    return

  end subroutine wann_calc_projection

  !=====================================!
  subroutine wann_write_xyz()
    !=====================================!
    !                                     !
    ! Write xyz file with Wannier centres !
    !                                     !
    !=====================================!

    use w90_io,         only: seedname,io_file_unit,io_date,stdout
    use w90_parameters, only: translate_home_cell,num_wann,wannier_centres, &
                              lenconfac,real_lattice,recip_lattice,iprint, &
                              num_atoms,atoms_symbol,atoms_pos_cart,&
                              num_species,atoms_species_num
    use w90_utility,    only : utility_translate_home

    implicit none

    integer          :: iw,ind,xyz_unit,nsp,nat
    character(len=9) :: cdate, ctime
    real(kind=dp)    :: wc(3,num_wann) 

    wc = wannier_centres

    if (translate_home_cell) then
       do iw=1,num_wann
          call utility_translate_home(wc(:,iw),real_lattice,recip_lattice)
       enddo
    endif

    if (iprint>2) then
       write(stdout,'(1x,a)') 'Final centres (translated to home cell for writing xyz file)'
       do iw=1,num_wann
          write(stdout,888) iw,(wc(ind,iw)*lenconfac,ind=1,3)
       end do
       write(stdout,'(1x,a78)') repeat('-',78)
       write(stdout,*)
    endif

    xyz_unit=io_file_unit()
    open(xyz_unit,file=trim(seedname)//'_centres.xyz',form='formatted')
    write(xyz_unit,'(i6)') num_wann + num_atoms
    call io_date(cdate,ctime)
    write(xyz_unit,*) 'Wannier centres, written by Wannier90 on'//cdate//' at '//ctime
    do iw=1,num_wann
       write(xyz_unit,'("X",6x,3(f14.8,3x))') (wc(ind,iw),ind=1,3)
    end do
    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          write(xyz_unit,'(a2,5x,3(f14.8,3x))') atoms_symbol(nsp),atoms_pos_cart(:,nat,nsp)
       end do
    end do
    close(xyz_unit)

    write(stdout,'(/a)') ' Wannier centres written to file '//trim(seedname)//'_centres.xyz'
    
    return

888 format(2x,'WF centre',i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )')
    
  end subroutine wann_write_xyz


    !=================================================================!
     subroutine wann_write_vdw_data()
    !=================================================================!
    !                                                                 !
    ! Write a file with Wannier centres, spreads and occupations for  !
    ! post-processing computation of vdW C6 coeffients.               !
    !                                                                 !
    ! Based on code written by Lampros Andrinopoulos.                 !
    !=================================================================!

    use w90_io,         only: seedname,io_file_unit,io_date,stdout,io_error
    use w90_parameters, only: translate_home_cell,num_wann,wannier_centres, &
                              lenconfac,real_lattice,recip_lattice,iprint, &
                              atoms_symbol,atoms_pos_cart,num_species, &
                              atoms_species_num,wannier_spreads,u_matrix, &
                              u_matrix_opt,num_kpts,have_disentangled, &
                              num_valence_bands,num_elec_per_state,write_vdw_data
    use w90_utility,    only :  utility_translate_home
    use w90_constants,  only :  cmplx_0,eps6
!!$    use w90_disentangle, only : ndimfroz

    implicit none
    
    integer          :: iw,vdw_unit,r,s,k,m,ierr,ndim
    real(kind=dp)    :: wc(3,num_wann)
    real(kind=dp)    :: ws(num_wann)
    complex(kind=dp), allocatable :: f_w(:,:),v_matrix(:,:) !f_w2(:,:)
 
    wc = wannier_centres
    ws = wannier_spreads 

    ! translate Wannier centres to the home unit cell
    do iw=1,num_wann
       call utility_translate_home(wc(:,iw),real_lattice,recip_lattice)
    enddo

    allocate(f_w(num_wann, num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating f_w in wann_write_vdw_data')
    
!!$    ! aam: remove f_w2 at end
!!$    allocate(f_w2(num_wann, num_wann),stat=ierr)
!!$    if (ierr/=0) call io_error('Error in allocating f_w2 in wann_write_vdw_data')

    if (have_disentangled) then

       ! dimension of occupied subspace
       if (num_valence_bands.le.0) call io_error('Please set num_valence_bands in seedname.win')
       ndim=num_valence_bands

       allocate(v_matrix(ndim, num_wann),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating V_matrix in wann_write_vdw_data')
       
       ! aam: initialise
       f_w(:,:) = cmplx_0
       v_matrix(:,:) = cmplx_0
!!$       f_w2(:,:) = cmplx_0

       ! aam: IN THE END ONLY NEED DIAGONAL PART, SO COULD SIMPLIFY...
       ! aam: calculate V = U_opt . U    
       do s=1,num_wann
          do k=1,ndim
             do m=1,num_wann
                v_matrix(k,s) = v_matrix(k,s) + u_matrix_opt(k,m,1) * u_matrix(m,s,1)
             enddo
          enddo
       enddo
       
       ! aam: calculate f = V^dagger . V
       do r=1,num_wann
          do s=1,num_wann
             do k=1,ndim
                f_w(r,s) = f_w(r,s) + v_matrix(k,s) * conjg(v_matrix(k,r))
             enddo
          enddo
       enddo
       
!!$       ! original formulation
!!$       do r=1,num_wann
!!$          do s=1,num_wann
!!$             do nkp=1,num_kpts
!!$                do k=1,ndimfroz(nkp)
!!$                   do m=1,num_wann
!!$                      do l=1,num_wann
!!$                         f_w2(r,s) = f_w2(r,s) + &
!!$                              u_matrix_opt(k,m,nkp) * u_matrix(m,s,nkp) * &
!!$                              conjg(u_matrix_opt(k,l,nkp)) * conjg(u_matrix(l,r,nkp))
!!$                      end do
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$       end do

!!$       ! test equivalence
!!$       do r=1,num_wann
!!$          do s=1,num_wann
!!$             if (abs(real(f_w(r,s),dp)-real(f_w2(r,s),dp)).gt.eps6) then
!!$                write(*,'(i6,i6,f16.10,f16.10)') r,s,real(f_w(r,s),dp),real(f_w2(r,s),dp)
!!$             endif
!!$             if (abs(aimag(f_w(r,s))-aimag(f_w2(r,s))).gt.eps6) then
!!$                write(*,'(a,i6,i6,f16.10,f16.10)') 'Im: ',r,s,aimag(f_w(r,s)),aimag(f_w2(r,s))
!!$             endif
!!$          enddo
!!$       enddo
!!$       write(*,*) ' done vdw ' 

    else
       ! for valence only, all occupancies are unity
       f_w(:,:) = 1.0_dp
    endif

    ! aam: write the seedname.vdw file directly here
    vdw_unit=io_file_unit()
    open(unit=vdw_unit,file=trim(seedname)//'.vdw',action='write')
    if (have_disentangled) then
       write(vdw_unit,'(a)') 'disentangle T'
    else
       write(vdw_unit,'(a)') 'disentangle F'
    endif
    write(vdw_unit,'(a)') 'amalgamate F'
    write(vdw_unit,'(a,i3)') 'degeneracy', num_elec_per_state
    write(vdw_unit,'(a)') 'num_frag 2'
    write(vdw_unit,'(a)') 'num_wann'
    write(vdw_unit,'(i3,1x,i3)') num_wann/2, num_wann/2
    write(vdw_unit,'(a)') 'tol_occ 0.9'
    write(vdw_unit,'(a)') 'pxyz'
    write(vdw_unit,'(a)') 'F F F'
    write(vdw_unit,'(a)') 'F F F'
    write(vdw_unit,'(a)') 'tol_dist 0.05'
    write(vdw_unit,'(a)') 'centres_spreads_occ'
    write(vdw_unit,'(a)') 'ang'
    do iw=1,num_wann
       write(vdw_unit,'(4(f13.10,1x),1x,f11.8)') wc(1:3,iw),ws(iw),real(f_w(iw,iw))
    end do
    close(vdw_unit)

    write(stdout,'(/a/)') ' vdW data written to file '//trim(seedname)//'.vdw'

    if (have_disentangled) then
       deallocate(v_matrix,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating v_matrix in wann_write_vdw_data')
    endif

!!$    deallocate(f_w2,stat=ierr)
!!$    if (ierr/=0) call io_error('Error in deallocating f_w2 in wann_write_vdw_data')
    deallocate(f_w,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating f_w in wann_write_vdw_data')
    
    return
    
  end subroutine wann_write_vdw_data


  !========================================!
  subroutine wann_check_unitarity()
    !========================================!

    use w90_constants,  only : dp,cmplx_1,cmplx_0,eps5
    use w90_io,         only : io_stopwatch,io_error,stdout
    use w90_parameters, only : num_kpts, num_wann,  &
                               u_matrix, timing_level 

    implicit none

    integer :: nkp,i,j,m
    complex(kind=dp) :: ctmp1,ctmp2

    if (timing_level>1) call io_stopwatch('wann: check_unitarity',1)

    do nkp = 1, num_kpts  
       do i = 1, num_wann  
          do j = 1, num_wann  
             ctmp1 = cmplx_0  
             ctmp2 = cmplx_0  
             do m = 1, num_wann  
                ctmp1 = ctmp1 + u_matrix (i, m, nkp) * conjg (u_matrix (j, m, nkp) )  
                ctmp2 = ctmp2 + u_matrix (m, j, nkp) * conjg (u_matrix (m, i, nkp) )  
             enddo
             if ( (i.eq.j) .and. (abs (ctmp1 - cmplx_1 ) .gt. eps5) ) &
                  then
                write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                     ctmp1
                call io_error('wann_check_unitarity: error 1')  
             endif
             if ( (i.eq.j) .and. (abs (ctmp2 - cmplx_1 ) .gt. eps5) ) &
                  then
                write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                     ctmp2
                call io_error('wann_check_unitarity: error 2')  
             endif
             if ( (i.ne.j) .and. (abs (ctmp1) .gt. eps5) ) then  
                write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                     ctmp1
                call io_error('wann_check_unitarity: error 3')  
             endif
             if ( (i.ne.j) .and. (abs (ctmp2) .gt. eps5) ) then  
                write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
                     ctmp2
                call io_error('wann_check_unitarity: error 4')  
             endif
          enddo
       enddo
    enddo

    if (timing_level>1) call io_stopwatch('wann: check_unitarity',2)

    return

  end subroutine wann_check_unitarity


  !========================================!
  subroutine wann_write_r2mn()
    !========================================!
    !                                        !
    ! Write seedname.r2mn file               !
    !                                        !
    !========================================!

    use w90_constants,  only : dp
    use w90_io,         only : seedname,io_file_unit,io_error
    use w90_parameters, only : num_kpts, num_wann, nntot, wb,  &
                               m_matrix
      
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
                     ! [GP-begin, Apr13, 2012: corrected sign inside "real"]
                     ( 2.0_dp * delta - real(m_matrix(nw1,nw2,nn,nkp) + &
                     conjg(m_matrix(nw2,nw1,nn,nkp)),kind=dp) )
                     ! [GP-end]
             enddo
          enddo
          r2ave_mn = r2ave_mn / real(num_kpts,dp)  
          write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn  
       enddo
    enddo
    close(r2mnunit)  
    
    return
    
158   call io_error('Error opening file '//trim(seedname)//'.r2mn in wann_write_r2mn')
      
  end subroutine wann_write_r2mn


  !========================================!
  subroutine wann_svd_omega_i()
  !========================================!
  
    use w90_constants,  only : dp,cmplx_0
    use w90_io,         only : io_stopwatch,io_error,stdout
    use w90_parameters, only: num_wann, num_kpts, nntot, wb, &
                              m_matrix, lenconfac, length_unit, &
                              timing_level

    implicit none

    complex(kind=dp), allocatable  :: cv1(:,:),cv2(:,:)
    complex(kind=dp), allocatable  :: cw1(:),cw2(:)  
    complex(kind=dp), allocatable  :: cpad1 (:)  
    real(kind=dp),    allocatable  :: singvd (:)  

    integer :: ierr, info
    integer :: nkp,nn,nb,na,ind
    real(kind=dp) :: omt1,omt2,omt3

    if (timing_level>1) call io_stopwatch('wann: svd_omega_i',1)

    allocate( cw1 (10 * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cw1 in wann_svd_omega_i')
    allocate( cw2 (10 * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cw2 in wann_svd_omega_i')
    allocate( cv1 (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cv1 in wann_svd_omega_i')
    allocate( cv2 (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cv2 in wann_svd_omega_i')
    allocate( singvd (num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating singvd in wann_svd_omega_i')
    allocate( cpad1 (num_wann * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cpad1 in wann_svd_omega_i')

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
    if (ierr/=0) call io_error('Error in deallocating cpad1 in wann_svd_omega_i')
    deallocate(singvd,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating singvd in wann_svd_omega_i')
    deallocate(cv2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cv2 in wann_svd_omega_i')
    deallocate(cv1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cv1 in wann_svd_omega_i')
    deallocate(cw2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cw2 in wann_svd_omega_i')
    deallocate(cw1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cw1 in wann_svd_omega_i')

    if (timing_level>1) call io_stopwatch('wann: svd_omega_i',2)

    return
    
  end subroutine wann_svd_omega_i

  !==================================================================!
  subroutine wann_main_gamma
    !==================================================================!
    !                                                                  !
    ! Calculate the Unitary Rotations to give                          !
    !            Maximally Localised Wannier Functions                 !
    !                      Gamma version                               !
    !===================================================================  
    use w90_constants,  only : dp,cmplx_1,cmplx_0
    use w90_io,         only : stdout,io_error,io_time,io_stopwatch
    use w90_parameters, only : num_wann,num_iter,wb, &
         nntot,u_matrix,m_matrix,num_kpts,iprint, &
         num_print_cycles,num_dump_cycles,omega_invariant, &
         param_write_chkpt,length_unit,lenconfac, &
         proj_site,real_lattice,write_r2mn,guiding_centres, &
         num_guide_cycles,num_no_guide_iter,timing_level, &
         write_proj,have_disentangled,conv_tol,conv_window, &
         wannier_centres,write_xyz,wannier_spreads,omega_total, &
         omega_tilde,write_vdw_data
    use w90_utility,    only : utility_frac_to_cart,utility_zgemm

    implicit none

    type(localisation_vars) :: old_spread
    type(localisation_vars) :: wann_spread

    ! guiding centres
    real(kind=dp), allocatable :: rguide (:,:)  
    integer :: irguide

    ! local arrays used and passed in subroutines
    real(kind=dp),    allocatable :: m_w(:,:,:)  
    complex(kind=dp), allocatable :: csheet(:,:,:)
    real(kind=dp),    allocatable :: sheet (:,:,:)
    real(kind=dp),    allocatable :: rave(:,:),r2ave(:),rave2(:)  

    !local arrays not passed into subroutines
    complex(kind=dp), allocatable  :: u0(:,:,:)
    complex(kind=dp), allocatable  :: uc_rot(:,:)
    real(kind=dp),    allocatable  :: ur_rot(:,:)
    complex(kind=dp), allocatable  :: cz (:,:)  

    real(kind=dp) :: sqwb
    integer       :: i,n,nn,iter,ind,ierr,iw
    integer       :: tnntot
    logical       :: lprint,ldump
    real(kind=dp), allocatable :: history(:)
    logical                    :: lconverged

    if (timing_level>0) call io_stopwatch('wann: main_gamma',1)

    first_pass=.true.

    ! Allocate stuff

    allocate(history(conv_window),stat=ierr)
    if (ierr/=0) call io_error('Error allocating history in wann_main_gamma')

!!$    if (.not.allocated(ph_g)) then
!!$       allocate(  ph_g(num_wann),stat=ierr )
!!$       if (ierr/=0) call io_error('Error in allocating ph_g in wann_main_gamma')
!!$       ph_g = cmplx_1
!!$    endif

    ! module data
    allocate( rnkb (num_wann, nntot, num_kpts),stat=ierr    )
    if (ierr/=0) call io_error('Error in allocating rnkb in wann_main_gamma')
    allocate( ln_tmp (num_wann, nntot, num_kpts), stat=ierr    )
    if (ierr/=0) call io_error('Error in allocating ln_tmp in wann_main_gamma')
    
    rnkb=0.0_dp
    tnntot=2*nntot

    ! sub vars passed into other subs 
    allocate( m_w(num_wann,num_wann,tnntot),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating m_w in wann_main_gamma')
    allocate( csheet (num_wann, nntot, num_kpts), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating csheet in wann_main_gamma')
    allocate( sheet (num_wann, nntot, num_kpts), stat=ierr    )
    if (ierr/=0) call io_error('Error in allocating sheet in wann_main_gamma')
    allocate( rave (3, num_wann), stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating rave in wann_main_gamma')
    allocate( r2ave (num_wann), stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating r2ave in wann_main_gamma')
    allocate( rave2 (num_wann), stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating rave2 in wann_main_gamma')
    allocate( rguide (3, num_wann)   )
    if (ierr/=0) call io_error('Error in allocating rguide in wann_main_gamma')

    csheet=cmplx_1
    sheet=0.0_dp;rave=0.0_dp;r2ave=0.0_dp;rave2=0.0_dp;rguide=0.0_dp

    ! sub vars not passed into other subs
    allocate(u0 (num_wann, num_wann, num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating u0 in wann_main_gamma')
    allocate(uc_rot(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating uc_rot in wann_main_gamma')
    allocate(ur_rot(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ur_rot in wann_main_gamma')
    allocate( cz (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cz in wann_main_gamma')

    cz=cmplx_0
    
    ! store original U before rotating
!!$    ! phase factor ph_g is applied to u_matrix
!!$    ! NB: ph_g is applied to u_matrix_opt if (have_disentangled)
!!$    if (have_disentangled) then
!!$       u0=u_matrix 
!!$    else
!!$       do iw=1,num_wann
!!$          u0(iw,:,:)= conjg(ph_g(iw))*u_matrix(iw,:,:)
!!$       end do
!!$    endif
    u0=u_matrix

!!$    lguide = .false.
    ! guiding centres are not neede for orthorhombic systems
    if (nntot .eq. 3) guiding_centres = .false.

    if (guiding_centres ) then
    ! initialise rguide to projection centres (Cartesians in units of Ang)
!!$       if ( use_bloch_phases) then
!!$          lguide = .true.
!!$       else
       do n=1,num_wann
          call utility_frac_to_cart(proj_site(:,n),rguide(:,n),real_lattice)
       enddo
!!$       endif
    endif

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
!!$    if (guiding_centres.and.(num_no_guide_iter.le.0)) then
!!$       if (nntot.gt.3) call wann_phases(csheet,sheet,rguide,irguide)
!!$       irguide=1
!!$    endif
    if (guiding_centres.and.(num_no_guide_iter.le.0)) then
       call wann_phases(csheet,sheet,rguide,irguide)
       irguide=1
    endif

    !  weight m_matrix first to reduce number of operations
    !  m_w : weighted real matrix
    do nn = 1, nntot  
       sqwb=sqrt(wb(nn))
       m_w(:,:,2*nn-1)=sqwb*real(m_matrix(:,:,nn,1),dp)
       m_w(:,:,2*nn)=sqwb*aimag(m_matrix(:,:,nn,1))
    end do

    ! calculate initial centers and spread
    call wann_omega_gamma(m_w,csheet,sheet,rave,r2ave,rave2,wann_spread)

    ! public variables
    omega_total = wann_spread%om_tot
    omega_invariant = wann_spread%om_i
    omega_tilde = wann_spread%om_d + wann_spread%om_od

    ! Public array of Wannier centres and spreads
    wannier_centres = rave
    wannier_spreads = r2ave - rave2

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
    write(stdout,'(1x,i6,2x,E12.3,19x,F18.10,3x,F8.2,2x,a)') &
         iter,(wann_spread%om_tot-old_spread%om_tot)*lenconfac**2,&
         wann_spread%om_tot*lenconfac**2,io_time(),'<-- CONV'
    write(stdout,'(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
         'O_D=',wann_spread%om_d*lenconfac**2,' O_OD=',wann_spread%om_od*lenconfac**2,&
         ' O_TOT=',wann_spread%om_tot*lenconfac**2,' <-- SPRD'
    write(stdout,'(1x,a78)') repeat('-',78) 

    lconverged=.false. 

    ! initialize ur_rot
    ur_rot=0.0_dp
    do i=1,num_wann
       ur_rot(i,i)=1.0_dp
    end do

    ! main iteration loop

    do iter=1,num_iter

       lprint=.false.
       if ( (mod(iter,num_print_cycles).eq.0) .or. (iter.eq.1) &
            .or. (iter.eq.num_iter)) lprint=.true.

       ldump=.false.
       if ( (num_dump_cycles.gt.0) .and. (mod(iter,num_dump_cycles).eq.0) ) ldump=.true.

       if(lprint) write(stdout,'(1x,a,i6)') 'Cycle: ',iter

!!$       ! initialize rguide as rave for use_bloch_phases
!!$       if ( (iter.gt.num_no_guide_iter) .and. lguide ) then 
!!$          rguide(:,:) = rave(:,:)
!!$          lguide = .false.
!!$       endif
!!$       if ( guiding_centres.and.(iter.gt.num_no_guide_iter) &
!!$            .and.(mod(iter,num_guide_cycles).eq.0) ) then
!!$          if(nntot.gt.3) call wann_phases(csheet,sheet,rguide,irguide)    
!!$          irguide=1
!!$       endif

       if ( guiding_centres.and.(iter.gt.num_no_guide_iter) &
            .and.(mod(iter,num_guide_cycles).eq.0)) then
          call wann_phases(csheet,sheet,rguide,irguide,m_w)    
          irguide=1
       endif

       call internal_new_u_and_m_gamma()

       call wann_spread_copy(wann_spread,old_spread)
       
       ! calculate the new centers and spread
       call wann_omega_gamma(m_w,csheet,sheet,rave,r2ave,rave2,wann_spread)

       ! print the new centers and spreads
       if ( lprint ) then
          do iw=1,num_wann
             write(stdout,1000) iw,(rave(ind,iw)*lenconfac,ind=1,3),&
                  (r2ave(iw) - rave2(iw))*lenconfac**2
          end do
          write(stdout,1001) (sum(rave(ind,:))*lenconfac,ind=1,3), &
               (sum(r2ave)-sum(rave2))*lenconfac**2
          write(stdout,*)
          write(stdout,'(1x,i6,2x,E12.3,19x,F18.10,3x,F8.2,2x,a)') &
               iter,(wann_spread%om_tot-old_spread%om_tot)*lenconfac**2,&
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

       ! Public array of Wannier centres and spreads
       wannier_centres = rave
       wannier_spreads = r2ave - rave2

       ! Public variables
       omega_total = wann_spread%om_tot
       omega_tilde = wann_spread%om_d + wann_spread%om_od


       if (ldump) then
          uc_rot(:,:)=cmplx(ur_rot(:,:),0.0_dp,dp)
          call  utility_zgemm(u_matrix,u0,'N',uc_rot,'N',num_wann)
          call param_write_chkpt('postdis')
       endif

       if (conv_window.gt.1) call internal_test_convergence_gamma()

       if (lconverged) then
          write(stdout,'(/13x,a,es10.3,a,i2,a)') &
               '<<<     Delta <',conv_tol,&
               '  over ',conv_window,' iterations     >>>'
          write(stdout,'(13x,a/)')  '<<< Wannierisation convergence criteria satisfied >>>'
          exit
       endif

    enddo
    ! end of the minimization loop

    ! update M
    do nn = 1, nntot
       sqwb=1.0_dp/sqrt(wb(nn))  
       m_matrix(:,:,nn,1)=sqwb*cmplx(m_w(:,:,2*nn-1),m_w(:,:,2*nn),dp)
    end do
    ! update U
    uc_rot(:,:)=cmplx(ur_rot(:,:),0.0_dp,dp)
    call  utility_zgemm(u_matrix,u0,'N',uc_rot,'N',num_wann)

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

    if (write_xyz) call wann_write_xyz()

    if (guiding_centres) call wann_phases(csheet,sheet,rguide,irguide)

    ! unitarity is checked
!!$    call internal_check_unitarity()
    call wann_check_unitarity()

    ! write extra info regarding omega_invariant
!!$    if (iprint>2) call internal_svd_omega_i()
    if (iprint>2) call wann_svd_omega_i()

    ! write matrix elements <m|r^2|n> to file
!!$    if (write_r2mn) call internal_write_r2mn()
    if (write_r2mn) call wann_write_r2mn()

    ! calculate and write projection of WFs on original bands in outer window
    if (have_disentangled .and. write_proj) call wann_calc_projection()

    ! aam: write data required for vdW utility
    if (write_vdw_data) call wann_write_vdw_data()

    ! deallocate sub vars not passed into other subs
    deallocate(cz,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating cz in wann_main_gamma')
    deallocate(ur_rot,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating ur_rot in wann_main_gamma')
    deallocate(uc_rot,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating uc_rot in wann_main_gamma')
    deallocate(u0, stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating u0 in wann_main_gamma')

    ! deallocate sub vars passed into other subs
    deallocate(rguide,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rguide in wann_main_gamma')
    deallocate(rave2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rave2 in wann_main_gamma')
    deallocate(rave,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating rave in wann_main_gamma')
    deallocate(sheet,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sheet in wann_main_gamma')
    deallocate(csheet,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating csheet in wann_main_gamma')
    deallocate(m_w,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating m_w in wann_main_gamma')
    
    ! deallocate module data
    deallocate( ln_tmp , stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating ln_tmp in wann_main_gamma')
    deallocate( rnkb,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating rnkb in wann_main_gamma')

    deallocate(history,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating history in wann_main_gamma')

    if (timing_level>0) call io_stopwatch('wann: main_gamma',2)

    return

1000 format(2x,'WF centre and spread', &
         &       i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)

1001 format(2x,'Sum of centres and spreads', &
         &       1x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)

  contains
    
    !===============================================!
    subroutine internal_new_u_and_m_gamma()
      !===============================================!

      use w90_constants, only : pi, eps10

      implicit none

      real(kind=dp) :: theta, twotheta
      real(kind=dp) :: a11, a12, a21, a22
      real(kind=dp) :: cc, ss, rtmp1, rtmp2
      real(kind=dp), parameter :: pifour=0.25_dp*pi
      integer       :: nn, nw1, nw2, nw3

      if (timing_level>1) call io_stopwatch('wann: main_gamma: new_u_and_m_gamma',1)

loop_nw1: do nw1=1,num_wann
loop_nw2: do nw2=nw1+1,num_wann
   
           a11=0.0_dp ; a12=0.0_dp ; a22=0.0_dp 
           do nn=1,tnntot
              a11=a11+(m_w(nw1,nw1,nn)-m_w(nw2,nw2,nn))**2
              a12=a12+m_w(nw1,nw2,nn)*(m_w(nw1,nw1,nn)-m_w(nw2,nw2,nn))
              a22=a22+m_w(nw1,nw2,nn)**2
           end do
           a12 = 2.0_dp * a12
           a22 = 4.0_dp * a22
           a21 = a22-a11
           if(abs(a12).gt.eps10) then
             twotheta=0.5_dp*(a21+sqrt(a21**2+4.0_dp*a12**2))/a12
             theta=0.5_dp*atan(twotheta)
           elseif (a21.lt.eps10 ) then
             theta=0.0_dp
           else
             theta=pifour
           endif
           cc=cos(theta)
           ss=sin(theta)

           ! update M               
           do nn=1,tnntot
              ! MR
              do nw3=1,num_wann
                 rtmp1=m_w(nw3,nw1,nn)*cc+m_w(nw3,nw2,nn)*ss
                 rtmp2=-m_w(nw3,nw1,nn)*ss+m_w(nw3,nw2,nn)*cc
                 m_w(nw3,nw1,nn)=rtmp1
                 m_w(nw3,nw2,nn)=rtmp2
              end do
              ! R^+ M R
              do nw3=1,num_wann
                 rtmp1=cc*m_w(nw1,nw3,nn)+ss*m_w(nw2,nw3,nn)
                 rtmp2=-ss*m_w(nw1,nw3,nn)+cc*m_w(nw2,nw3,nn)
                 m_w(nw1,nw3,nn)=rtmp1
                 m_w(nw2,nw3,nn)=rtmp2
              end do
           end do
           ! update U : U=UR
           do nw3=1,num_wann
              rtmp1=ur_rot(nw3,nw1)*cc+ur_rot(nw3,nw2)*ss
              rtmp2=-ur_rot(nw3,nw1)*ss+ur_rot(nw3,nw2)*cc
              ur_rot(nw3,nw1)=rtmp1
              ur_rot(nw3,nw2)=rtmp2
           end do
        end do loop_nw2
      end do loop_nw1

      if (timing_level>1) call io_stopwatch('wann: main_gamma: new_u_and_m_gamma',2)

      return

    end subroutine internal_new_u_and_m_gamma

    !===============================================!
    subroutine internal_test_convergence_gamma()
      !===============================================!
      !                                               !
      ! Determine whether minimisation of non-gauge-  !
      ! invariant spread is converged                 !
      !                                               !
      !===============================================!

      implicit none

      real(kind=dp) :: delta_omega
      integer :: j,ierr
      real(kind=dp), allocatable :: temp_hist(:)
        
      allocate(temp_hist(conv_window),stat=ierr)
      if (ierr/=0) call io_error('Error allocating temp_hist in wann_main')

      delta_omega=wann_spread%om_tot-old_spread%om_tot
       
      if (iter.le.conv_window) then
         history(iter) = delta_omega
      else
         temp_hist = eoshift(history,1,delta_omega)
         history = temp_hist
      endif

      lconverged = .false.

      if (iter.ge.conv_window) then
!!$         write(stdout,*) (history(j),j=1,conv_window)
         do j=1,conv_window
            if ( abs(history(j)).gt.conv_tol ) exit
            lconverged = .true.
         enddo
      endif

      deallocate(temp_hist,stat=ierr)
      if (ierr/=0) call io_error('Error deallocating temp_hist in wann_main_gamma')

      return

    end subroutine internal_test_convergence_gamma



!!$    !========================================!
!!$    subroutine internal_check_unitarity()
!!$    !========================================!
!!$
!!$      implicit none
!!$
!!$      integer :: nkp,i,j,m
!!$      complex(kind=dp) :: ctmp1,ctmp2
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',1)
!!$
!!$      do nkp = 1, num_kpts  
!!$         do i = 1, num_wann  
!!$            do j = 1, num_wann  
!!$               ctmp1 = cmplx_0  
!!$               ctmp2 = cmplx_0  
!!$               do m = 1, num_wann  
!!$                  ctmp1 = ctmp1 + u_matrix (i, m, nkp) * conjg (u_matrix (j, m, nkp) )  
!!$                  ctmp2 = ctmp2 + u_matrix (m, j, nkp) * conjg (u_matrix (m, i, nkp) )  
!!$               enddo
!!$               if ( (i.eq.j) .and. (abs (ctmp1 - cmplx_1 ) .gt. eps5) ) &
!!$                    then
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp1
!!$                  call io_error('wann_main: unitariety error 1')  
!!$               endif
!!$               if ( (i.eq.j) .and. (abs (ctmp2 - cmplx_1 ) .gt. eps5) ) &
!!$                    then
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp2
!!$                  call io_error('wann_main: unitariety error 2')  
!!$               endif
!!$               if ( (i.ne.j) .and. (abs (ctmp1) .gt. eps5) ) then  
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp1
!!$                  call io_error('wann_main: unitariety error 3')  
!!$               endif
!!$               if ( (i.ne.j) .and. (abs (ctmp2) .gt. eps5) ) then  
!!$                  write ( stdout , * ) ' ERROR: unitariety of final U', nkp, i, j, &
!!$                       ctmp2
!!$                  call io_error('wann_main: unitariety error 4')  
!!$               endif
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: check_unitarity',2)
!!$
!!$      return
!!$
!!$    end subroutine internal_check_unitarity


!!$    !========================================!
!!$    subroutine internal_write_r2mn()
!!$    !========================================!
!!$    !                                        !
!!$    ! Write seedname.r2mn file               !
!!$    !                                        !
!!$    !========================================!
!!$      use w90_io, only: seedname,io_file_unit,io_error
!!$      
!!$      implicit none
!!$
!!$      integer :: r2mnunit,nw1,nw2,nkp,nn
!!$      real(kind=dp) :: r2ave_mn,delta
!!$
!!$      ! note that here I use formulas analogue to Eq. 23, and not to the
!!$      ! shift-invariant Eq. 32 .
!!$      r2mnunit=io_file_unit()
!!$      open(r2mnunit,file=trim(seedname)//'.r2mn',form='formatted',err=158)
!!$      do nw1 = 1, num_wann  
!!$         do nw2 = 1, num_wann  
!!$            r2ave_mn = 0.0_dp  
!!$            delta = 0.0_dp  
!!$            if (nw1.eq.nw2) delta = 1.0_dp  
!!$            do nkp = 1, num_kpts  
!!$               do nn = 1, nntot
!!$                  r2ave_mn = r2ave_mn + wb(nn) * &
!!$                       ! [GP-begin, Apr13, 2012: corrected sign inside "real"]
!!$                       ( 2.0_dp * delta - real(m_matrix(nw1,nw2,nn,nkp) + &
!!$                       conjg(m_matrix(nw2,nw1,nn,nkp)),kind=dp) )
!!$               enddo
!!$            enddo
!!$            r2ave_mn = r2ave_mn / real(num_kpts,dp)  
!!$            write (r2mnunit, '(2i6,f20.12)') nw1, nw2, r2ave_mn  
!!$         enddo
!!$      enddo
!!$      close(r2mnunit)  
!!$      
!!$      return
!!$      
!!$158   call io_error('Error opening file '//trim(seedname)//'.r2mn in wann_main')
!!$      
!!$    end subroutine internal_write_r2mn



!!$    !========================================!
!!$    subroutine internal_svd_omega_i()
!!$    !========================================!
!!$
!!$      implicit none
!!$
!!$      complex(kind=dp), allocatable  :: cv1(:,:),cv2(:,:)
!!$      complex(kind=dp), allocatable  :: cw1(:),cw2(:)  
!!$      complex(kind=dp), allocatable  :: cpad1 (:)  
!!$      real(kind=dp),    allocatable  :: singvd (:)  
!!$
!!$      integer :: nkp,nn,nb,na,ind
!!$      real(kind=dp) :: omt1,omt2,omt3
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',1)
!!$
!!$      allocate( cw1 (10 * num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cw1 in wann_main')
!!$      allocate( cw2 (10 * num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cw2 in wann_main')
!!$      allocate( cv1 (num_wann, num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cv1 in wann_main')
!!$      allocate( cv2 (num_wann, num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cv2 in wann_main')
!!$      allocate( singvd (num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating singvd in wann_main')
!!$      allocate( cpad1 (num_wann * num_wann),stat=ierr  )
!!$      if (ierr/=0) call io_error('Error in allocating cpad1 in wann_main')
!!$
!!$      cw1=cmplx_0; cw2=cmplx_0; cv1=cmplx_0; cv2=cmplx_0; cpad1=cmplx_0 
!!$      singvd=0.0_dp
!!$
!!$      ! singular value decomposition
!!$      omt1 = 0.0_dp ; omt2 = 0.0_dp ; omt3 = 0.0_dp
!!$      do nkp = 1, num_kpts  
!!$         do nn = 1, nntot  
!!$            ind = 1  
!!$            do nb = 1, num_wann  
!!$               do na = 1, num_wann  
!!$                  cpad1 (ind) = m_matrix (na, nb, nn, nkp)  
!!$                  ind = ind+1  
!!$               enddo
!!$            enddo
!!$            call zgesvd ('A', 'A', num_wann, num_wann, cpad1, num_wann, singvd, cv1, &
!!$                 num_wann, cv2, num_wann, cw1, 10 * num_wann, cw2, info)
!!$            if (info.ne.0) then  
!!$               call io_error('ERROR: Singular value decomp. zgesvd failed')  
!!$            endif
!!$
!!$            do nb = 1, num_wann  
!!$               omt1 = omt1 + wb(nn) * (1.0_dp - singvd (nb) **2)  
!!$               omt2 = omt2 - wb(nn) * (2.0_dp * log (singvd (nb) ) )  
!!$               omt3 = omt3 + wb(nn) * (acos (singvd (nb) ) **2)  
!!$            enddo
!!$         enddo
!!$      enddo
!!$      omt1 = omt1 / real(num_kpts,dp)  
!!$      omt2 = omt2 / real(num_kpts,dp)  
!!$      omt3 = omt3 / real(num_kpts,dp)  
!!$      write ( stdout , * ) ' '  
!!$      write(stdout,'(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ',&
!!$           omt1*lenconfac**2,'('//trim(length_unit)//'^2)'
!!$      write(stdout,'(2x,a,f15.9,1x,a)') '                 -2log s = ',&
!!$           omt2*lenconfac**2,'('//trim(length_unit)//'^2)'
!!$      write(stdout,'(2x,a,f15.9,1x,a)') '                  acos^2 = ',&
!!$           omt3*lenconfac**2,'('//trim(length_unit)//'^2)'
!!$
!!$      deallocate(cpad1,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cpad1 in wann_main')
!!$      deallocate(singvd,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating singvd in wann_main')
!!$      deallocate(cv2,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cv2 in wann_main')
!!$      deallocate(cv1,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cv1 in wann_main')
!!$      deallocate(cw2,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cw2 in wann_main')
!!$      deallocate(cw1,stat=ierr)
!!$      if (ierr/=0) call io_error('Error in deallocating cw1 in wann_main')
!!$
!!$      if (timing_level>1) call io_stopwatch('wann: main: svd_omega_i',2)
!!$
!!$      return
!!$      
!!$    end subroutine internal_svd_omega_i


  end subroutine wann_main_gamma

  !==================================================================!
  subroutine wann_omega_gamma(m_w,csheet,sheet,rave,r2ave,rave2,wann_spread)
    !==================================================================!
    !                                                                  !
    !   Calculate the Wannier Function spread                          !
    !                                                                  !
    !===================================================================  
    use w90_parameters, only : num_wann,nntot,wbtot,wb,bk,&
                           omega_invariant,timing_level
    use w90_io,         only : io_error,io_stopwatch

    implicit none

    real(kind=dp)   , intent(in)  :: m_w (:,:,:)
    complex(kind=dp), intent(in)  :: csheet (:,:,:)
    real(kind=dp)   , intent(in)  :: sheet (:,:,:)
    real(kind=dp)   , intent(out) :: rave (:,:)
    real(kind=dp)   , intent(out) :: r2ave (:)
    real(kind=dp)   , intent(out) :: rave2 (:)
    type(localisation_vars)    , intent(out)  :: wann_spread

    !local variables
    real(kind=dp) :: summ, brn
    real(kind=dp), allocatable :: m_w_nn2(:)
    integer :: ind, nn,m,n,iw, rn, cn, ierr

    if (timing_level>1) call io_stopwatch('wann: omega_gamma',1)

    allocate( m_w_nn2(num_wann),stat=ierr )
    if (ierr/=0) call io_error('Error in allocating m_w_nn2 in wann_omega_gamma')

    if (nntot .eq. 3 ) then
       do nn = 1, nntot
          rn = 2*nn-1
          cn = 2*nn
          do n = 1, num_wann
             ln_tmp(n,nn,1)= atan2(m_w(n,n,cn),m_w(n,n,rn)) 
          end do
       end do
    else
       do nn = 1, nntot
          rn = 2*nn-1
          cn = 2*nn
          do n = 1, num_wann
             ln_tmp(n,nn,1)= aimag(log(csheet(n,nn,1)*cmplx(m_w(n,n,rn),m_w(n,n,cn),dp))) &
                             -sheet(n,nn,1) 
           end do
       end do
    endif

    rave  = 0.0_dp
    do iw = 1, num_wann  
       do ind = 1, 3  
          do nn = 1, nntot  
             rave(ind,iw) = rave(ind,iw) - wb(nn) * bk(ind,nn,1) &
                   *ln_tmp(iw,nn,1)
          enddo
       enddo
    enddo

    rave2 = 0.0_dp
    do iw = 1, num_wann  
       rave2(iw) = sum(rave(:,iw)*rave(:,iw))
    enddo

    m_w_nn2 = 0.0_dp
    r2ave = wbtot
    do iw = 1, num_wann  
       do nn = 1, nntot  
          rn = 2*nn -1 
          cn = 2*nn
          m_w_nn2(iw) = m_w_nn2(iw) + m_w(iw,iw,rn)**2 + m_w(iw,iw,cn)**2
          r2ave(iw) = r2ave(iw) + wb(nn)*ln_tmp(iw,nn,1)**2 
       enddo
       r2ave(iw) = r2ave(iw) - m_w_nn2(iw) 
    enddo

    if (first_pass) then
       summ = 0.0_dp  
       do nn = 1, nntot  
          rn = 2*nn -1
          cn = 2*nn 
          do m = 1, num_wann  
             do n = 1, num_wann  
                summ = summ + m_w(n,m,rn)**2+m_w(n,m,cn)**2
             enddo
          enddo
       enddo
       wann_spread%om_i = wbtot*real(num_wann,dp) - summ
       first_pass=.false.
    else 
       wann_spread%om_i=omega_invariant
    endif

    wann_spread%om_od = wbtot*real(num_wann,dp) - sum(m_w_nn2(:)) - wann_spread%om_i 

    if (nntot.eq.3) then 
       wann_spread%om_d = 0.0_dp  
    else
       wann_spread%om_d = 0.0_dp  
       do nn = 1, nntot  
          do n = 1, num_wann  
             brn = sum(bk(:,nn,1)*rave(:,n))
             wann_spread%om_d = wann_spread%om_d + wb(nn) &
                  * ( ln_tmp(n,nn,1) + brn)**2
          enddo
       enddo
    end if

    wann_spread%om_tot = wann_spread%om_i + wann_spread%om_d + wann_spread%om_od

    deallocate( m_w_nn2,stat=ierr )
    if (ierr/=0) call io_error('Error in deallocating m_w_nn2 in wann_omega_gamma')

    if (timing_level>1) call io_stopwatch('wann: omega_gamma',2)

    return  

  end subroutine wann_omega_gamma

end module w90_wannierise
