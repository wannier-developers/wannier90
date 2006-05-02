!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
!        Authors:                                            !
!              Arash A. Mostofi   (MIT)                      !
!              Jonathan R. Yates  (LBNL and UC Berkeley)     !
!                                                            !
!  Please cite                                               !
!                                                            !
!  [ref] A. A. Mostofi, J. R. Yates, N. Marzari, I. Souza    !
!        and D. Vanderbilt, http://www.wannier.org/          !
!                                                            !
!  in any publications arising from the use of this code.    !
!                                                            !
!  Wannier90 is based on routines written by N. Marzari,     !
!  I. Souza and D. Vanderbilt. For the method please cite    !
!                                                            !
!  [ref] N. Marzari and D. Vanderbilt,                       !
!        Phys. Rev. B 56 12847 (1997)                        !
!                                                            !
!  [ref] I. Souza, N. Marzari and D. Vanderbilt,             !
!        Phys. Rev. B 65 035109 (2001)                       !
!                                                            !
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

program wannier

  use constants
  use parameters
  use io
  use kmesh
  use disentangle
  use overlap
  use wannierise
  use plot
 
  implicit none

  real(kind=dp) time0,time1,time2
  character(len=9) :: stat,pos,cdate,ctime
  logical :: wout_found

  time0=io_time()

  call io_get_seedname

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.werr')
  call io_date(cdate,ctime)
  write(stdout,*)  'Wannier90: Execution started on ',cdate,' at ',ctime
  call param_read
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
  call param_write_header
  call param_write

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to read parameters  ',time1-time0,' (sec)'

  call io_stopwatch('kmesh_get',1)
  call kmesh_get
  call io_stopwatch('kmesh_get',2)

  ! Sort out restarts
  if (restart.eq.' ') then  ! start a fresh calculation
     write(stdout,'(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
     call param_read_chkpt
     call param_read_um
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
        case default        ! for completeness... (it is already trapped in param_read)
           call io_error('Value of restart not recognised in wann_prog')
     end select
  endif

  if (postproc_setup) then
     call kmesh_write
     call kmesh_dealloc
     call param_dealloc
     write(stdout,'(1x,a25,f11.3,a)') 'Time to write kmesh      ',io_time(),' (sec)'
     write(stdout,'(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
     stop
  endif

  time2=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to get kmesh        ',time2-time1,' (sec)'

  call io_stopwatch('overlap_read',1)

  call overlap_read

  call io_stopwatch('overlap_read',2)

  time1=io_time()
  write(stdout,'(/1x,a25,f11.3,a)') 'Time to read overlaps    ',time1-time2,' (sec)'

  have_disentangled = .false.

  if (disentanglement) then

     call io_stopwatch('dis_main',1)
     call dis_main
     call io_stopwatch('dis_main',2)

     have_disentangled=.true.

     time2=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time to disentangle bands',time2-time1,' (sec)'     

  endif

  call param_write_chkpt('postdis')
  call param_write_um

  time2=io_time()

  call io_stopwatch('wann_main',1)
1001 call wann_main
  call io_stopwatch('wann_main',2)

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time for wannierise      ',time1-time2,' (sec)'     

  call param_write_chkpt('postwann')

2002 continue

  if (wannier_plot .or. bands_plot .or. fermi_surface_plot) then
     call plot_main
     time2=io_time()
     write(stdout,'(1x,a25,f11.3,a)') 'Time for plotting        ',time2-time1,' (sec)'     
  end if

  call overlap_dealloc
  call kmesh_dealloc
  call param_dealloc

  write(stdout,'(1x,a25,f11.3,a)') 'Total Execution Time     ',io_time(),' (sec)'

  call io_print_timings()

  write(stdout,*) 
  write(stdout,'(1x,a)') 'All done: wannier90 exiting'
 
  close(stdout)


100 format(1x,'Execution started on ',a9,' at ',a9)
101 format(1x,'Reading parameters')

end program wannier
  

