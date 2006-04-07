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

module wannierise

  use constants

  implicit none

  private

  public :: wann_main

  ! Data to avoid large allocation within itteration loop
  complex(kind=dp), allocatable  :: cr (:,:,:,:)   
  complex(kind=dp), allocatable  :: crt (:,:,:,:)  
  real(kind=dp),    allocatable  :: rnkb (:,:,:)   

  type localisation_vars
     real(kind=dp) :: om_i   
     real(kind=dp) :: om_d   
     real(kind=dp) :: om_od  
     real(kind=dp) :: om_tot 
     real(kind=dp) :: om_1   
     real(kind=dp) :: om_2   
     real(kind=dp) :: om_3   
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
    use constants,  only : dp,cmplx_1,cmplx_0
    use io,         only : stdout,io_error,io_time,io_stopwatch
    use parameters, only : num_wann,num_cg_steps,num_iter,wb,nnlist, &
                           nntot,wbtot,u_matrix,m_matrix,num_kpts,iprint, &
                           num_print_cycles,num_dump_cycles,omega_invariant, &
                           param_read_um,param_write_um,length_unit,lenconfac

    implicit none


    type(localisation_vars) :: old_spread
    type(localisation_vars) :: wann_spread
    type(localisation_vars) :: trial_spread


    ! parameters from input file

    logical :: lcg, lrguide


    ! local arrays used and passed in subroutines
    complex(kind=dp), allocatable :: csheet (:,:,:)
    complex(kind=dp), allocatable :: cdodq1 (:,:,:)  
    complex(kind=dp), allocatable :: cdodq2 (:,:,:)  
    complex(kind=dp), allocatable :: cdodq3 (:,:,:)  
    real(kind=dp), allocatable    :: sheet (:,:,:)
    real(kind=dp), allocatable    :: rave (:,:), r2ave (:)  , rave2 (:)  
    real(kind=dp), allocatable    :: rguide (:,:)  


    !local arrays not passed into subroutines
    complex(kind=dp), allocatable  :: cwschur1 (:), cwschur2 (:)  
    complex(kind=dp), allocatable  :: cwschur3 (:), cwschur4 (:)  
    complex(kind=dp), allocatable  :: cpad1 (:)  
    complex(kind=dp), allocatable  :: cpad2 (:)  
    complex(kind=dp), allocatable  :: cpad3 (:)  
    complex(kind=dp), allocatable  :: cw1 (:)  
    complex(kind=dp), allocatable  :: cw2 (:)  
    complex(kind=dp), allocatable  :: cu0 (:,:,:)  ! copy of unitary
    complex(kind=dp), allocatable  :: cm0 (:,:,:,:)  ! copy of overlaps
    complex(kind=dp), allocatable  :: cdq (:,:,:)  
    complex(kind=dp), allocatable  :: cz (:,:)  
    complex(kind=dp), allocatable  :: cv1 (:,:)  
    complex(kind=dp), allocatable  :: cv2 (:,:)  
    complex(kind=dp), allocatable  :: cdodq (:,:,:)  
    complex(kind=dp), allocatable  :: cmtmp (:,:)  
    complex(kind=dp), allocatable  :: cdqkeep (:,:,:)  
    real(kind=dp),    allocatable  :: singvd (:)  

    complex(kind=dp), allocatable  :: tmp_cdq(:,:)

    complex(kind=dp) :: cfunc_exp1,cfunc_exp2,ctmp1,ctmp2
    complex(kind=dp) :: cfunc_exp,cfunc_exp3
    real(kind=dp) :: omt1,omt2,omt3,doda0
    real(kind=dp) :: alpha,fac,shift,trial_step
    real(kind=dp) :: falphamin,eqa,eqb,eqc,delta,alphamin,r2ave_mn
    real(kind=dp) :: gcfac,gcnorm1,gcnorm0
    integer :: nkp,i,j,nn,nb,na,m,n,iter,ind,ierr,loop_wann,ncg,bis_loop
    integer :: nw1,nw2,nkp2,ncgfix,nsdim,nrguide,irguide,info
    logical :: linput,select,lprint,ldump
    real(kind=dp) :: alphamin_quad,falphamin_quad,om_trial

    ! == bisection parameters -- slightly experimental ==!
    integer, parameter       :: max_bis_iter = 3
    real(kind=dp), parameter :: bis_factor   = 0.1_dp
    real(kind=dp), parameter :: mono_thresh_omega = 1.0e-11_dp
    real(kind=dp), parameter :: mono_thresh_grad  = 1.0e-6_dp
    ! ===================================================!

    linput = .false.
    nrguide = 0
    lrguide = .false.
    lcg = .true.

    ! Allocate stuff

    ! module data
    allocate(  cr (num_wann, num_wann, nntot, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cr in wann_main')
    allocate(  crt (num_wann, num_wann, nntot, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating crt in wann_main')
    allocate( rnkb (num_wann, nntot, num_kpts),stat=ierr    )     
    if (ierr/=0) call io_error('Error in allocating rnkb in wann_main')
    cr=(0.d0,0.d0);  crt=(0.d0,0.d0);  rnkb=(0.d0,0.d0)


    ! sub vars passed into other subs 
    allocate( csheet (num_wann, nntot, num_kpts), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating csheet in wann_main')
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
    csheet=(0.d0,0.d0);cdodq1=(0.d0,0.d0);cdodq2=(0.d0,0.d0)
    cdodq3=(0.d0,0.d0);rave2=0.d0
    sheet=0.d0;rave=0.d0; r2ave=0.d0; rguide=0.d0

    ! sub vars not passed into other subs
    allocate( cwschur1 (num_wann), cwschur2 (10 * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cwshur1 in wann_main')
    allocate( cwschur3 (num_wann), cwschur4 (num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cwshur3 in wann_main')
    allocate( cpad1 (num_wann * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cpad1 in wann_main')
    allocate( cpad2 (num_wann * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cpad2 in wann_main')
    allocate( cpad3 (num_wann * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cpad3 in wann_main')
    allocate( cw1 (10 * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cw1 in wann_main')
    allocate( cw2 (10 * num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cw2 in wann_main')
    allocate( cu0 (num_wann, num_wann, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cu0 in wann_main')
    allocate( cm0 (num_wann, num_wann, nntot, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cm0 in wann_main')
    allocate( cdq (num_wann, num_wann, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cdq in wann_main')
    allocate( cz (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cz in wann_main')
    allocate( cv1 (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cv1 in wann_main')
    allocate( cv2 (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cv2 in wann_main')
    allocate( cdodq (num_wann, num_wann, num_kpts),stat=ierr ) 
    if (ierr/=0) call io_error('Error in allocating cdodq in wann_main')
    allocate( cmtmp (num_wann, num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cmtmp in wann_main')
    allocate( cdqkeep (num_wann, num_wann, num_kpts),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating cdqkeep in wann_main')
    allocate( singvd (num_wann),stat=ierr  )
    if (ierr/=0) call io_error('Error in allocating singvd in wann_main')
    cwschur1=(0.d0,0.d0); cwschur2=(0.d0,0.d0); cwschur3=(0.d0,0.d0)
    cwschur4=(0.d0,0.d0); cpad1 =(0.d0,0.d0);cpad2=(0.d0,0.d0)
    cpad3=(0.d0,0.d0);  cw1=(0.d0,0.d0); cw2=(0.d0,0.d0); cu0=(0.d0,0.d0)
    cm0 =(0.d0,0.d0);cdq =(0.d0,0.d0);cz=(0.d0,0.d0);  cv1=(0.d0,0.d0)
    cv2=(0.d0,0.d0); cdodq=(0.d0,0.d0); cmtmp=(0.d0,0.d0); cdqkeep=(0.d0,0.d0)
    singvd=0.d0
    sheet=0.d0
    !  csheet=exp((0.d0,1.d0)*sheet)
    csheet=cmplx_1
    gcnorm1=0.d0
    gcnorm0=0.d0

    
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

    allocate(tmp_cdq(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error allocating tmp_cdq in wann_main')

    ! parameter for line minimisation
    trial_step = 2.0_dp

    ! calculate initial centers and spread
    call omega(csheet,sheet,rave,r2ave,rave2,wann_spread)

    ncg  = 0
    iter = 0
    old_spread%om_tot = 0.0_dp
    gcnorm1=0

    ! print initial state
    write(stdout,'(1x,a78)') repeat('-',78) 
    write(stdout,'(1x,a)') 'Initial State'
    do loop_wann=1,num_wann
       write(stdout,1000) loop_wann,(rave(ind,loop_wann)*lenconfac,ind=1,3),&
            (r2ave(loop_wann) - rave2(loop_wann))*lenconfac**2
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

    omega_invariant = wann_spread%om_i

    ! main iteration loop
    do iter=1,num_iter

       lprint=.false.
       if ( (mod(iter,num_print_cycles).eq.0) .or. (iter.eq.1) &
            .or. (iter.eq.num_iter)) lprint=.true.

       ldump=.false.
       if ( (num_dump_cycles.gt.0) .and. (mod(iter,num_dump_cycles).eq.0) ) ldump=.true.

       if(lprint) write(stdout,'(1x,a,i6)') 'Cycle: ',iter

       ! calculate gradient of omega
       call domega(csheet,sheet,rave,cdodq1,cdodq2,cdodq3,cdodq)

       ! gcnorm1 = Tr[gradient . gradient] -- NB gradient is anti-Hermitian
       gcnorm1=0.0_dp
       do nkp = 1, num_kpts  
          do n = 1, num_wann  
             do m = 1, num_wann  
                gcnorm1 = gcnorm1 + real(cdodq(m,n,nkp)*conjg(cdodq(m,n,nkp)),dp)
             enddo
          enddo
       enddo

       if ( lprint .and. iprint>2 ) &
            write(stdout,*) ' LINE --> Iteration                     :',iter

       ! calculate CG coefficient
       if ( (iter.eq.1) .or. (ncg.ge.num_cg_steps) ) then
          gcfac = 0.0_dp                 ! Steepest descents
          ncg   = 0
       else
          if (gcnorm0.gt.epsilon(1.0_dp)) then
             gcfac = gcnorm1/gcnorm0     ! Fletcher-Reeves CG coefficient
             ! prevent CG coefficient from getting too large
             if (gcfac.gt.3.0_dp) then
                if ( lprint .and. iprint>2 ) &
                     write(stdout,*) ' LINE --> CG coeff too large. Resetting CG:',gcfac
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
       doda0 = 0.0_dp
       do nkp = 1, num_kpts  
          do m = 1, num_wann  
             do n = 1, num_wann  
                doda0 = doda0 + real(cdq(m,n,nkp)*cdodq(n,m,nkp),dp)  
             enddo
          enddo
       enddo
       doda0 = doda0 / (4.0_dp*wbtot)

       ! check if search direction is uphill
       if (doda0.gt.0.0_dp) then
          ! if doing a CG step then reset CG
          if (ncg.gt.0) then
             if ( lprint .and. iprint>2 ) &
                  write(stdout,*) ' LINE --> Search direction uphill: resetting CG'
             cdq(:,:,:) = cdodq(:,:,:)
             ncg = 0
             gcfac = 0.0_dp
             ! re-calculate gradient along search direction
             doda0 = 0.0_dp
             do nkp = 1, num_kpts  
                do m = 1, num_wann  
                   do n = 1, num_wann  
                      doda0 = doda0 + real(cdq(m,n,nkp)*cdodq(n,m,nkp),dp)  
                   enddo
                enddo
             enddo
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

       ! save search direction 
       cdqkeep(:,:,:) = cdq(:,:,:)

       ! take trial step
       cdq(:,:,:)=cdqkeep(:,:,:)*( trial_step / (4.0_dp*wbtot) ) 

       ! update U and M
       call internal_new_u_and_m()

       ! calculate spread at trial step
       call omega(csheet,sheet,rave,r2ave,rave2,trial_spread)

       ! Calculate optimal parabolic step
       fac = abs(trial_spread%om_tot - wann_spread%om_tot)
       if ( fac.gt.tiny(1.0_dp) ) then
          fac = 1.0_dp/fac
          if ( trial_spread%om_tot .gt. wann_spread%om_tot ) then
             shift =  1.0_dp
          else
             shift = -1.0_dp
          endif
       else
          fac   = 1.0e5_dp
          shift = fac*trial_spread%om_tot - fac*wann_spread%om_tot
       endif
       eqb = fac*doda0  
       eqa = shift - eqb*trial_step
       if ( abs(eqa/(fac*wann_spread%om_tot)) .gt. epsilon(1.0_DP) ) then
          alphamin_quad = -eqb / (2.0_dp * eqa) * (trial_step**2)
          falphamin_quad = wann_spread%om_tot &
               - ( eqb*eqb / (4.0_dp * fac * eqa) ) * (trial_step**2)
       endif

       ! set line search coefficient
       alphamin=alphamin_quad
       falphamin=falphamin_quad

       ! print line search information
       if ( lprint .and. iprint>2 ) then
          write(stdout,*) ' LINE --> Spread at initial point       :',wann_spread%om_tot*lenconfac**2
          write(stdout,*) ' LINE --> Spread at trial step          :',trial_spread%om_tot*lenconfac**2
          write(stdout,*) ' LINE --> Slope along search direction  :',doda0*lenconfac**2
          write(stdout,*) ' LINE --> ||SD gradient||^2             :',gcnorm1*lenconfac**2
          write(stdout,*) ' LINE --> Trial step length             :',trial_step
          write(stdout,*) ' LINE --> Optimal parabolic step length :',alphamin
          write(stdout,*) ' LINE --> Spread at predicted minimum   :',falphamin*lenconfac**2
          write(stdout,*) ' LINE --> CG coefficient                :',gcfac
       endif

       ! bisection loop to try to avoid uphill moves
       bisection: do bis_loop=0,max_bis_iter

          ! take optimal step
          cdq(:,:,:) = cdqkeep(:,:,:) * ( alphamin / (4.0_dp*wbtot) ) 

          ! update U and M
          call internal_new_u_and_m()

          ! copy original spread to old_spread
          if (bis_loop.eq.0) call wann_spread_copy(wann_spread,old_spread)

          ! calculate the new centers and spread
          call omega(csheet,sheet,rave,r2ave,rave2,wann_spread)

          ! check for non-monotonic convergence
          if (max_bis_iter.gt.0 .and. bis_loop.lt.max_bis_iter) then
             if ( ((wann_spread%om_tot - old_spread%om_tot).gt.mono_thresh_omega) &
                  .and. (abs(doda0).lt.mono_thresh_grad) .and. (alphamin*doda0.lt.0.0_dp) )  then
                alphamin = alphamin * bis_factor
                if ( lprint .and. iprint>2 ) then 
                   write(stdout,*) ' LINE --> Bisection iteration             :',bis_loop+1
                   write(stdout,*) ' LINE --> Optimal step went uphill. Spread:',&
                        wann_spread%om_tot*lenconfac**2
                   write(stdout,*) ' LINE --> Reducing optimal step length to :',alphamin
                endif
             else
                exit bisection
             endif
          endif

       enddo bisection

       ! print the new centers and spreads
       if(lprint) then
          do loop_wann=1,num_wann
             write(stdout,1000) loop_wann,(rave(ind,loop_wann)*lenconfac,ind=1,3),&
                  (r2ave(loop_wann) - rave2(loop_wann))*lenconfac**2
          end do
          write(stdout,1001) (sum(rave(ind,:))*lenconfac,ind=1,3), (sum(r2ave)-sum(rave2))*lenconfac**2
          write(stdout,*)
          write(stdout,'(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
               iter,(wann_spread%om_tot-old_spread%om_tot)*lenconfac**2,sqrt(abs(gcnorm1))*lenconfac,&
               wann_spread%om_tot*lenconfac**2,io_time(),'<-- CONV'
          write(stdout,'(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
               'O_D=',wann_spread%om_d*lenconfac**2,' O_OD=',wann_spread%om_od*lenconfac**2,&
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



    deallocate(tmp_cdq,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating tmp_cdq in wann_main')

    write(stdout,'(1x,a)') 'Final State'
    do loop_wann=1,num_wann
       write(stdout,1000) loop_wann,(rave(ind,loop_wann)*lenconfac,ind=1,3),&
            (r2ave(loop_wann) - rave2(loop_wann))*lenconfac**2
    end do
    write(stdout,1001) (sum(rave(ind,:))*lenconfac,ind=1,3), (sum(r2ave)-sum(rave2))*lenconfac**2
    write(stdout,*)
!    write(stdout,'(3x,a,f15.9)') '     Spreads (Ang)         Omega I      = ',wann_spread%om_i  
    write(stdout,'(3x,a21,a,f15.9)') '     Spreads ('//trim(length_unit)//'^2)','       Omega I      = ',wann_spread%om_i*lenconfac**2
    write(stdout,'(3x,a,f15.9)') '     ================       Omega D      = ',wann_spread%om_d*lenconfac**2
    write(stdout,'(3x,a,f15.9)') '                            Omega OD     = ',wann_spread%om_od*lenconfac**2
!    write(stdout,'(3x,a,f15.9)') 'Final Spread (Ang)         Omega Total  = ',wann_spread%om_tot*lenconfac  
    write(stdout,'(3x,a21,a,f15.9)') 'Final Spread ('//trim(length_unit)//'^2)','       Omega Total  = ',wann_spread%om_tot*lenconfac**2  
    write(stdout,'(1x,a78)') repeat('-',78) 

    if (lrguide) call phases (  csheet, sheet, rguide, irguide)


    ! unitariety is checked
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


    if(iprint>2) then
       ! singular value decomposition
       omt1 = 0.d0  
       omt2 = 0.d0  
       omt3 = 0.d0
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
                omt1 = omt1 + wb (nkp, nn) * (1.d0 - singvd (nb) **2)  
                omt2 = omt2 - wb (nkp, nn) * (2.d0 * log (singvd (nb) ) )  
                omt3 = omt3 + wb (nkp, nn) * (acos (singvd (nb) ) **2)  
                !         write(*,'(i4,4x,f10.5,4x,3f10.5)')
                !    &     nb,singvd(nb),omt1,omt2,omt3
             enddo
          enddo
       enddo
       omt1 = omt1 / real(num_kpts,dp)  
       omt2 = omt2 / real(num_kpts,dp)  
       omt3 = omt3 / real(num_kpts,dp)  
       write ( stdout , * ) ' '  
       write(stdout,'(2x,a,f15.9,1x,a)') 'Omega Invariant:   1-s^2 = ',omt1*lenconfac**2,'('//trim(length_unit)//'^2)'
       write(stdout,'(2x,a,f15.9,1x,a)') '                 -2log s = ',omt2*lenconfac**2,'('//trim(length_unit)//'^2)'
       write(stdout,'(2x,a,f15.9,1x,a)') '                  acos^2 = ',omt3*lenconfac**2,'('//trim(length_unit)//'^2)'
       write ( stdout , * ) ' '  
    end if


    write ( stdout , * ) ' '  


    ! note that here I use formulas analogue to Eq. 23, and not to the
    ! shift-invariant Eq. 32 .
    open (20, file = 'wannier.r2_mn', form = 'formatted', status = 'unknown')
    do nw1 = 1, num_wann  
       do nw2 = 1, num_wann  
          r2ave_mn = 0.d0  
          delta = 0.d0  
          if (nw1.eq.nw2) delta = 1.d0  
          do nkp = 1, num_kpts  
             do nn = 1, nntot  
                r2ave_mn = r2ave_mn + wb (nkp, nn) * &
                     ( 2.d0 * delta - m_matrix(nw1,nw2,nn,nkp) - &
                     conjg(m_matrix(nw2,nw1,nn,nkp)) )
             enddo
          enddo
          r2ave_mn = r2ave_mn / real(num_kpts,dp)  
          write (20, '(2i4,f20.12)') nw1, nw2, r2ave_mn  
       enddo
    enddo
    close (20)  


    call param_write_um


    ! deallocate module data
    deallocate( rnkb,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating rnkb in wann_main')
    deallocate(  crt,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating crt in wann_main') 
    deallocate(  cr,stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating cr in wann_main') 


1000 format(2x,'WF centre and spread', &
         &       i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)

1001 format(2x,'Sum of centres and spreads', &
         &       1x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)


  contains


    !===============================================!
    subroutine internal_new_u_and_m()               
    !===============================================!
    !                                               !
    ! Update U and M matrices after a trial step    !
    !                                               !
    !===============================================!

      implicit none

      do nkp = 1, num_kpts  
         tmp_cdq(:,:) = cdq(:,:,nkp)
         call zgees ('V', 'N', select, num_wann, tmp_cdq, num_wann, nsdim, &
              cwschur1, cz, num_wann, cwschur2, 10 * num_wann, cwschur3, &
              cwschur4, info)
         if (info.ne.0) then  
            write ( stdout , * ) 'SCHUR: ', info  
            call io_error('wann_main: problem computing schur form 1') 
         endif
         cdq(:,:,nkp) = cmplx_0
         do j = 1, num_wann  
            do i = 1, num_wann  
               do m = 1, num_wann  
                  cdq (i, j, nkp) = cdq (i, j, nkp) + &
                       cz(i, m) * (exp(cwschur1(m))) * conjg(cz(j,m))
               enddo
            enddo
         enddo
      enddo

      ! the orbitals are rotated
      do nkp = 1, num_kpts  
         cmtmp = cmplx_0
         do j = 1, num_wann  
            do i = 1, num_wann  
               do m = 1, num_wann  
                  cmtmp (i, j) = cmtmp (i, j) + u_matrix (i, m, nkp) * cdq (m, j, nkp)  
               enddo
            enddo
         enddo
         u_matrix(:,:,nkp) = cmtmp(:,:)
      enddo

      ! and the M_ij are updated
      do nkp = 1, num_kpts  
         do nn = 1, nntot  
            nkp2 = nnlist (nkp, nn)  
            cmtmp = cmplx_0
            do j = 1, num_wann  
               do i = 1, num_wann  
                  do n = 1, num_wann  
                     do m = 1, num_wann  
                        cmtmp (i, j) = cmtmp (i, j) + conjg (cdq (m, i, nkp) ) * &
                             cdq (n, j, nkp2) * m_matrix (m, n, nn, nkp)
                     enddo
                  enddo
               enddo
            enddo
            m_matrix(:,:,nn,nkp) = cmtmp(:,:)
         enddo
      enddo

      return

    end subroutine internal_new_u_and_m


  end subroutine wann_main


  !==================================================================!
  subroutine phases (csheet, sheet, rguide, irguide)
    !==================================================================!
    !                                                                  !
    !                                                                  !
    !===================================================================  
    use parameters,     only : num_wann,m_matrix,nntot,neigh, &
         nnh,bk,bka,num_kpts
    use io,         only : stdout,io_error
    use utility,    only : utility_inv3

    implicit none

    complex(kind=dp), intent(out)   :: csheet (:,:,:)
    real(kind=dp)   , intent(out)   :: sheet (:,:,:)
    real(kind=dp)   , intent(inout) :: rguide (:,:)
    integer         , intent(in)    :: irguide

    !local
    complex(kind=dp)               :: csum (nntot/2)  
    real(kind=dp)                  ::  xx(nntot/2)
    real(kind=dp) ::  smat (3, 3), svec (3), sinv (3, 3)  
    real(kind=dp), parameter :: ci = (0.0_dp, 1.0_dp)  
    real(kind=dp) :: xx0,csumt,det,brn,pherr
    integer :: loop_wann,na,nkp,i,j,nn,ind,m,n

    ! AAM: variable pherr is given a value but never used

    csum=cmplx_0; xx=0.d0

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
       do j = 1, 3  
          do i = 1, 3  
             smat (j, i) = 0.d0  
          enddo
          svec (j) = 0.d0  
       enddo

       do nn = 1, nnh  
          if (nn.le.3) then  
             !         obtain xx with arbitrary branch cut choice
             xx (nn) = - aimag (log (csum (nn) ) )  
          else  
             !         obtain xx with branch cut choice guided by rguide
             xx0 = 0.d0  
             do j = 1, 3  
                xx0 = xx0 + bka (j, nn) * rguide (j, loop_wann)  
             enddo
             !         xx0 is expected value for xx
             csumt = exp (ci * xx0)  
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
             if (abs (det) .gt.1.e-06) then  
                !          to check that the first nn bka vectors are not
                !          linearly dependent - this is a change from original code
                if (irguide.ne.0) then  
                   do j = 1, 3  
                      rguide (j, loop_wann) = 0.d0  
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
             !           sheet (loop_wann, nn, nkp) = 0.d0
             do j = 1, 3  
                sheet (loop_wann, nn, nkp) = sheet (loop_wann, nn, nkp) + bk (j, nkp, nn) &
                     * rguide (j, loop_wann)
             enddo
             !           csheet (loop_wann, nn, nkp) = exp (ci * sheet (loop_wann, nn, nkp) )  
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
                brn = brn + bk (ind, nkp, nn) * rguide (ind, m)  
             enddo
             rnkb (m, nn, nkp) = rnkb (m, nn, nkp) + brn  
          enddo
       enddo
    enddo
    write ( stdout , * ) ' '  
    write ( stdout , * ) ' PHASES ARE SET USING THE GUIDING CENTERS'  
    write ( stdout , * ) ' '  
    do nkp = 1, num_kpts  
       do n = 1, num_wann  
          do nn = 1, nntot  
             pherr = aimag (log (csheet (n, nn, nkp) * m_matrix (n, n, nn, nkp) ) ) &
                  - sheet (n, nn, nkp) + rnkb (n, nn, nkp) - aimag (log (m_matrix (n, n, &
                  nn, nkp) ) )
          enddo
       enddo
    enddo

    return  

  end subroutine phases



  !==================================================================!
  subroutine omega ( csheet, sheet, rave, r2ave, rave2, wann_spread)
    !==================================================================!
    !                                                                  !
    !   Calculate the Wannier Function spread                          !
    !                                                                  !
    !===================================================================  
    use parameters,     only : num_wann,m_matrix,nntot,wb,bk,num_kpts,iprint
    use io,         only : stdout,io_error,io_stopwatch

    implicit none

    complex(kind=dp), intent(in)  :: csheet (:,:,:)
    real(kind=dp)   , intent(in)  :: sheet (:,:,:)
    real(kind=dp)   , intent(out) :: rave (:,:)
    real(kind=dp)   , intent(out) :: r2ave (:)
    real(kind=dp)   , intent(out) :: rave2 (:)
    type(localisation_vars)    , intent(out)  :: wann_spread

    !local variables
    real(kind=dp) :: bim2,sum
    real(kind=dp) :: brn,sqim, bim (3), rtot (3)  
    integer :: loop_wann,ind,nkp,nn,i,m,n




    rave  = 0.0_dp
    do loop_wann = 1, num_wann  
       do ind = 1, 3  
          !        rave (ind, loop_wann) = 0.0_dp  
          do nkp = 1, num_kpts  
             do nn = 1, nntot  
                rave (ind, loop_wann) = rave (ind, loop_wann) + wb (nkp, nn) * bk (ind, &
                     nkp, nn) * (aimag (log (csheet (loop_wann, nn, nkp) * m_matrix (loop_wann, &
                     loop_wann, nn, nkp) ) ) - sheet (loop_wann, nn, nkp) )
             enddo
          enddo
          !        rave (ind, loop_wann) = - rave (ind, loop_wann) / real(num_kpts,dp)  
       enddo
    enddo
    rave = -rave/real(num_kpts,dp)

    rave2 = 0.0_dp
    do loop_wann = 1, num_wann  
       !     rave2 (loop_wann) = 0.0_dp  
       do ind = 1, 3  
          rave2 (loop_wann) = rave2 (loop_wann) + rave (ind, loop_wann) **2  
       enddo
    enddo

    rtot=0.0_dp
    do ind = 1, 3  
       !     rtot (ind) = 0.0_dp  
       do loop_wann = 1, num_wann  
          rtot (ind) = rtot (ind) + rave (ind, loop_wann)  
       enddo
    enddo

    r2ave = 0.0_dp
    do loop_wann = 1, num_wann  
       !     r2ave (loop_wann) = 0.0_dp  
       do nkp = 1, num_kpts  
          do nn = 1, nntot  
             r2ave (loop_wann) = r2ave (loop_wann) + wb (nkp, nn) * (1.d0 - m_matrix (loop_wann, &
                  loop_wann, nn, nkp) * conjg (m_matrix (loop_wann, loop_wann, nn, nkp) ) + (aimag ( &
                  log (csheet (loop_wann, nn, nkp) * m_matrix (loop_wann, loop_wann, nn, nkp) ) ) &
                  - sheet (loop_wann, nn, nkp) ) **2)
          enddo
       enddo
       !     r2ave (loop_wann) = r2ave (loop_wann) / real(num_kpts,dp)  
    enddo
    r2ave = r2ave/real(num_kpts,dp)

    rave2=0.0_dp
    do loop_wann = 1, num_wann  
       !     rave2 (loop_wann) = 0.0_dp  
       do i = 1, 3  
          rave2 (loop_wann) = rave2 (loop_wann) + rave (i, loop_wann) * rave (i, loop_wann)  
       enddo
    enddo

    wann_spread%om_1 = 0.0_dp  
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          do loop_wann = 1, num_wann  
             wann_spread%om_1 = wann_spread%om_1 + wb (nkp, nn) * (1.0_dp - m_matrix (loop_wann, loop_wann, nn, &
                  nkp) * conjg (m_matrix (loop_wann, loop_wann, nn, nkp) ) )
          enddo
       enddo
    enddo
    wann_spread%om_1 = wann_spread%om_1 / real(num_kpts,dp)  

    wann_spread%om_2 = 0.0_dp  
    do loop_wann = 1, num_wann  
       sqim = 0.0_dp  
       do nkp = 1, num_kpts  
          do nn = 1, nntot  
             sqim = sqim + wb (nkp, nn) * ( (aimag (log (csheet (loop_wann, nn, &
                  nkp) * m_matrix (loop_wann, loop_wann, nn, nkp) ) ) - sheet (loop_wann, nn, nkp) ) ** &
                  2)
          enddo
       enddo
       sqim = sqim / real(num_kpts,dp)  
       wann_spread%om_2 = wann_spread%om_2 + sqim  
    enddo


    wann_spread%om_3 = 0.0_dp  
    do loop_wann = 1, num_wann  
       bim = 0.0_dp
       do ind = 1, 3  
          !        bim (ind) = 0.0_dp  
          do nkp = 1, num_kpts  
             do nn = 1, nntot  
                bim (ind) = bim (ind) + wb (nkp, nn) * bk (ind, nkp, nn) * &
                     (aimag (log (csheet (loop_wann, nn, nkp) * m_matrix (loop_wann, loop_wann, nn, nkp) &
                     ) ) - sheet (loop_wann, nn, nkp) )
             enddo
          enddo
          !        bim (ind) = bim (ind) / real(num_kpts,dp)  
       enddo
       bim = bim/real(num_kpts,dp)
       bim2 = 0.0_dp  
       do ind = 1, 3  
          bim2 = bim2 + bim (ind) * bim (ind)  
       enddo
       wann_spread%om_3 = wann_spread%om_3 - bim2  
    enddo



    !jry: Either the above (om1,2,3) or the following is redundant
    !     keep it in the code base for testing

    ! this only needs to be calculated on the first pass
    wann_spread%om_i = 0.0_dp  
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          sum = 0.0_dp  
          do m = 1, num_wann  
             do n = 1, num_wann  
                sum = sum + m_matrix (n, m, nn, nkp) * conjg (m_matrix (n, m, nn, nkp) )  
             enddo
          enddo
          wann_spread%om_i = wann_spread%om_i + wb (nkp, nn) * (real(num_wann,dp) - sum)  
       enddo
    enddo
    wann_spread%om_i = wann_spread%om_i / real(num_kpts,dp)  


    wann_spread%om_od = 0.0_dp  
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          sum = 0.0_dp  
          do m = 1, num_wann  
             do n = 1, num_wann  
                sum = sum + wb (nkp, nn) * m_matrix (n, m, nn, nkp) * conjg (m_matrix (n, m, &
                     nn, nkp) )
                if (m.eq.n) sum = sum - wb (nkp, nn) * m_matrix (n, m, nn, nkp) * conjg &
                     (m_matrix (n, m, nn, nkp) )
             enddo
          enddo
          wann_spread%om_od = wann_spread%om_od+sum  
       enddo
    enddo

    wann_spread%om_od = wann_spread%om_od / real(num_kpts,dp)  
    wann_spread%om_d = 0.0_dp  
    do nkp = 1, num_kpts  
       do nn = 1, nntot  
          sum = 0.0_dp  
          do n = 1, num_wann  
             brn = 0.0_dp  
             do ind = 1, 3  
                brn = brn + bk (ind, nkp, nn) * rave (ind, n)  
             enddo
             sum = sum + wb (nkp, nn) * (aimag (log (csheet (n, nn, nkp) &
                  * m_matrix (n, n, nn, nkp) ) ) - sheet (n, nn, nkp) + brn) **2
          enddo
          wann_spread%om_d = wann_spread%om_d+sum  
       enddo
    enddo

    wann_spread%om_d = wann_spread%om_d / real(num_kpts,dp)  
    wann_spread%om_tot = wann_spread%om_i + wann_spread%om_d+wann_spread%om_od  


    return  


  end subroutine omega





  !==================================================================!
  subroutine domega ( csheet, sheet, rave, cdodq1, cdodq2, cdodq3, &
       cdodq)
    !==================================================================!
    !                                                                  !
    !   Calculate the Gradient of the Wannier Function spread          !
    !                                                                  !
    !===================================================================  
    use parameters,     only : num_wann,wb,bk,nntot,m_matrix,num_kpts
    use io,         only : stdout,io_error,io_stopwatch
    implicit none


    complex(kind=dp), intent(in)  :: csheet (:,:,:)    
    complex(kind=dp), intent(out) :: cdodq (:,:,:)     
    complex(kind=dp), intent(out) :: cdodq1 (:,:,:)    
    complex(kind=dp), intent(out) :: cdodq2 (:,:,:)    
    complex(kind=dp), intent(out) :: cdodq3 (:,:,:)    
    real(kind=dp),    intent(in)  :: sheet (:,:,:)     
    real(kind=dp),    intent(out) :: rave (:,:)        


    ! local
    real(kind=dp) :: brn
    integer :: loop_wann,ind,nkp,nn,m,n
    complex(kind=dp) :: mnn



    ! recalculate rave
    rave = 0.0_dp
    do loop_wann = 1, num_wann  
       do ind = 1, 3  
          do nkp = 1, num_kpts  
             do nn = 1, nntot  
                rave(ind,loop_wann) = rave(ind,loop_wann) + &
                     wb(nkp,nn) * bk(ind,nkp,nn) * &
                     ( aimag(log(csheet(loop_wann,nn,nkp) * &
                     m_matrix(loop_wann,loop_wann,nn,nkp))) - sheet(loop_wann,nn,nkp) )
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
             do m=1,num_wann
                crt(m,n,nn,nkp) = m_matrix(m,n,nn,nkp) / mnn
                cr(m,n,nn,nkp)  = m_matrix(m,n,nn,nkp) * conjg(mnn)
             enddo
             brn = 0.0_dp  
             do ind = 1, 3  
                brn = brn + bk (ind, nkp, nn) * rave (ind, n)  
             enddo
             rnkb(n,nn,nkp) = brn
          enddo
       enddo
    enddo


    ! cd0dq(m,n,nkp) is calculated
    cdodq1 = cmplx_0
    cdodq2 = cmplx_0
    cdodq3 = cmplx_0
    do nkp = 1, num_kpts  
       do n = 1, num_wann  
          do m = 1, num_wann  
             do nn = 1, nntot  
                ! A[R^{k,b}]=(R-Rdag)/2
                cdodq1 (m, n, nkp) = cdodq1 (m, n, nkp) + wb (nkp, nn) * (cr (m, &
                     n, nn, nkp) / 2.0_dp - conjg (cr (n, m, nn, nkp) ) / 2.0_dp)
                ! -S[T^{k,b}]=-(T+Tdag)/2i ; T_mn = Rt_mn q_n
                cdodq2 (m, n, nkp) = cdodq2 (m, n, nkp) - wb (nkp, nn) * (crt (m, &
                     n, nn, nkp) * (aimag (log (csheet (n, nn, nkp) * m_matrix (n, n, nn, &
                     nkp) ) ) - sheet (n, nn, nkp) ) + conjg (crt (n, m, nn, nkp) &
                     * (aimag (log (csheet (m, nn, nkp) * m_matrix (m, m, nn, nkp) ) ) &
                     - sheet (m, nn, nkp) ) ) ) / (0.0_dp, 2.0_dp)
                cdodq3 (m, n, nkp) = cdodq3 (m, n, nkp) - wb (nkp, nn) * (crt (m, &
                     n, nn, nkp) * rnkb (n, nn, nkp) + conjg (crt (n, m, nn, nkp) &
                     * rnkb (m, nn, nkp) ) ) / (0.0_dp, 2.0_dp)
             enddo
          enddo
       enddo
    enddo
    cdodq1 = cdodq1 / real(num_kpts,dp) * 4.0_dp
    cdodq2 = cdodq2 / real(num_kpts,dp) * 4.0_dp
    cdodq3 = cdodq3 / real(num_kpts,dp) * 4.0_dp
    cdodq  = cdodq1 + cdodq2 + cdodq3



    return  


  end subroutine domega


  subroutine wann_spread_copy(orig,copy)

    implicit none

    type(localisation_vars), intent(in)  :: orig
    type(localisation_vars), intent(out) :: copy

    copy%om_i   =  orig%om_i  
    copy%om_d   =  orig%om_d  
    copy%om_od  =  orig%om_od 
    copy%om_tot =  orig%om_tot
    copy%om_1   =  orig%om_1  
    copy%om_2   =  orig%om_2  
    copy%om_3   =  orig%om_3  

    return

  end subroutine wann_spread_copy


!!$  subroutine wann_spread_write(wann_spread,old_spread)
!!$
!!$    use io, only : stdout
!!$    implicit none
!!$
!!$    type(localisation_vars),intent(in) :: wann_spread
!!$    type(localisation_vars),intent(in) :: old_spread
!!$
!!$    write ( stdout , * ) ' '  
!!$    write ( stdout, 1005) wann_spread%om_1,wann_spread%om_1-old_spread%om_1
!!$    write ( stdout, 1006) wann_spread%om_2,wann_spread%om_2-old_spread%om_2
!!$    write ( stdout, 1007) wann_spread%om_3,wann_spread%om_3-old_spread%om_3 
!!$    write ( stdout, 1008) wann_spread%om_2 + wann_spread%om_3,wann_spread%om_2-old_spread%om_2+&
!!$         wann_spread%om_3-old_spread%om_3  
!!$    write ( stdout , * ) ' '  
!!$
!!$    write ( stdout, 1010) wann_spread%om_i ,wann_spread%om_i-old_spread%om_i  
!!$    write ( stdout, 1011) wann_spread%om_d  ,wann_spread%om_d-old_spread%om_d
!!$    write ( stdout, 1012) wann_spread%om_od  ,wann_spread%om_od-old_spread%om_od 
!!$    write ( stdout, 1014) wann_spread%om_od+wann_spread%om_d  ,wann_spread%om_od-old_spread%om_od&
!!$         +wann_spread%om_d-old_spread%om_d 
!!$    write ( stdout , * ) ' '  
!!$    write ( stdout, 1013) wann_spread%om_tot  ,wann_spread%om_tot-old_spread%om_tot
!!$    write ( stdout , * ) ' '  
!!$
!!$1005 format(2x,'Omega 1    is   ',f18.12,3x,'Delta Omega 1    is  ',0pe16.8)  
!!$1006 format(2x,'Omega 2    is   ',f18.12,3x,'Delta Omega 2    is  ',0pe16.8)    
!!$1007 format(2x,'Omega 3    is   ',f18.12,3x,'Delta Omega 3    is  ',0pe16.8)  
!!$1008 format(2x,'Omega 2+3  is   ',f18.12,3x,'Delta Omega 2+3  is  ',0pe16.8)
!!$
!!$1010 format(2x,'Omega I    is   ',f18.12,3x,'Delta Omega I    is  ',0pe16.8)    
!!$1011 format(2x,'Omega D    is   ',f18.12,3x,'Delta Omega D    is  ',0pe16.8)  
!!$1012 format(2x,'Omega OD   is   ',f18.12,3x,'Delta Omega OD   is  ',0pe16.8)  
!!$1014 format(2x,'Omega D+OD is   ',f18.12,3x,'Delta Omega D+OD is  ',0pe16.8)  
!!$
!!$1013 format(2x,'Omega      is   ',f18.12,3x,'Delta Omega      is  ',0pe16.8)    
!!$
!!$  end subroutine wann_spread_write

end module wannierise
