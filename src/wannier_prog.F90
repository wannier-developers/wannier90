!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Wannier90 v2.0 authors:                                    !
!           Arash A. Mostofi   (Imperial College London)     !
!           Jonathan R. Yates  (University of Oxford)        !
!           Giovanni Pizzi     (EPFL, Switzerland)           !
!           Ivo Souza          (Universidad del Pais Vasco)  !
!                                                            !
! Contributorsr:                                             !
!          Young-Su Lee        (KIST, S. Korea)              !
!          Matthew Shelley     (Imperial College London)     !
!          Nicolas Poilvert    (Harvard)                     ! 
!                                                            !
!  Please cite                                               !
!                                                            !
!  [ref] A. A. Mostofi, J. R. Yates, Y.-S. Lee, I. Souza,    !
!        D. Vanderbilt and N. Marzari, "Wannier90: A Tool    !
!        for Obtaining Maximally Localised Wannier           !
!        Functions", Computer Physics Communications,        !
!        178, 685 (2008)                                     !
!                                                            !
!  in any publications arising from the use of this code.    !
!                                                            !
!  Wannier90 is based on Wannier77, written by N. Marzari,   !
!  I. Souza and D. Vanderbilt. For the method please cite    !
!                                                            !
!  [ref] N. Marzari and D. Vanderbilt,                       !
!        Phys. Rev. B 56 12847 (1997)                        !
!                                                            !
!  [ref] I. Souza, N. Marzari and D. Vanderbilt,             !
!        Phys. Rev. B 65 035109 (2001)                       !
!                                                            !
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

program wannier

  use w90_constants
  use w90_parameters
  use w90_io
  use w90_hamiltonian
  use w90_kmesh
  use w90_disentangle
  use w90_overlap
  use w90_wannierise
  use w90_plot
  use w90_transport
 
  implicit none

  real(kind=dp) time0,time1,time2
  character(len=9) :: stat,pos,cdate,ctime
  logical :: wout_found

  time0=io_time()

  library = .false.

  call io_get_seedname()

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.werr')
  call io_date(cdate,ctime)
  write(stdout,*)  'Wannier90: Execution started on ',cdate,' at ',ctime
  call param_read()
  close(stdout,status='delete')

  if (restart.eq.' ') then
     stat='replace'
     pos ='rewind'
  else
     inquire(file=trim(seedname)//'.wout',exist=wout_found)
     if (wout_found) then
        stat='old'
     else
        stat='replace'
     endif
     pos='append'
  endif

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.wout',status=trim(stat),position=trim(pos))
  call param_write_header()
  call param_write()

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to read parameters  ',time1-time0,' (sec)'

  if (transport .and. tran_read_ht) goto 3003

  call kmesh_get()
  call param_memory_estimate()

  ! Sort out restarts
  if (restart.eq.' ') then  ! start a fresh calculation
     write(stdout,'(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
     call param_read_chkpt()
!!$     call param_read_um
     select case (restart)
        case ('default')    ! continue from where last checkpoint was written
           write(stdout,'(/1x,a)',advance='no') 'Resuming a previous Wannier90 calculation '
           if (checkpoint.eq.'postdis') then 
              write(stdout,'(a/)') 'from wannierisation ...'
              goto 1001         ! go to wann_main
           elseif (checkpoint.eq.'postwann') then
              write(stdout,'(a/)') 'from plotting ...'
              goto 2002         ! go to plot_main
           else
              write(stdout,'(/a/)')
              call io_error('Value of checkpoint not recognised in wann_prog')
           endif
        case ('wannierise') ! continue from wann_main irrespective of value of last checkpoint
           write(stdout,'(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
           goto 1001
        case ('plot')       ! continue from plot_main irrespective of value of last checkpoint 
           write(stdout,'(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
           goto 2002       
        case ('transport')   ! continue from tran_main irrespective of value of last checkpoint 
           write(stdout,'(1x,a/)') 'Restarting Wannier90 from transport routines ...'
           goto 3003       
        case default        ! for completeness... (it is already trapped in param_read)
           call io_error('Value of restart not recognised in wann_prog')
     end select
  endif

  if (postproc_setup) then
     call kmesh_write()
     call kmesh_dealloc()
     call param_dealloc()
     write(stdout,'(1x,a25,f11.3,a)') 'Time to write kmesh      ',io_time(),' (sec)'
     write(stdout,'(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
     stop
  endif

  time2=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to get kmesh        ',time2-time1,' (sec)'

  call overlap_read()

  time1=io_time()
  write(stdout,'(/1x,a25,f11.3,a)') 'Time to read overlaps    ',time1-time2,' (sec)'

  have_disentangled = .false.

  if (disentanglement) then
     call dis_main()
     have_disentangled=.true.
     time2=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time to disentangle bands',time2-time1,' (sec)'     
  endif

  call param_write_chkpt('postdis')
!!$  call param_write_um

1001 time2=io_time()

  if (.not. gamma_only) then
       call wann_main()
  else
       call wann_main_gamma()
  end if

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time for wannierise      ',time1-time2,' (sec)'     

  call param_write_chkpt('postwann')

2002 time2=io_time()

  if (wannier_plot .or. bands_plot .or. fermi_surface_plot .or. hr_plot) then
     call plot_main()
     time1=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time for plotting        ',time1-time2,' (sec)'
  end if

3003 time2=io_time()

  if (transport) then
     call tran_main()
     time1=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time for transport       ',time1-time2,' (sec)'
     if (tran_read_ht) goto 4004
  end if

  call transport_dealloc()
  call hamiltonian_dealloc()
  call overlap_dealloc()
  call kmesh_dealloc()
  call param_dealloc()

4004 continue 

  write(stdout,'(1x,a25,f11.3,a)') 'Total Execution Time     ',io_time(),' (sec)'

  if (timing_level>0) call io_print_timings()

  write(stdout,*) 
  write(stdout,'(1x,a)') 'All done: wannier90 exiting'
 
  close(stdout)



end program wannier
  

