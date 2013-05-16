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

module w90_io

  use w90_constants, only : dp

  implicit none

  private


#ifdef MPI
  include 'mpif.h'
#endif

  integer, public, save           :: stdout
  character(len=50), public, save :: seedname
  integer, parameter, public      :: maxlen = 120  ! Max column width of input file
  logical, public, save           :: post_proc_flag  ! Set post_processing from cmd line

  type timing_data
     integer :: ncalls           
     real(kind=DP) :: ctime      
     real(kind=DP) :: ptime     
     character(len=60) :: label      
  end type timing_data

  integer, parameter :: nmax = 100
  type(timing_data) :: clocks(nmax)
  integer, save     :: nnames=0


  public :: io_stopwatch
  public :: io_print_timings
  public :: io_get_seedname
  public :: io_time
  public :: io_date
  public :: io_error
  public :: io_file_unit

contains


  !==================================================================!
  subroutine io_stopwatch(tag,mode)
  !==================================================================!
  ! Stopwatch to time parts of the code                              !
  !==================================================================!

    implicit none

    character(len=*), intent(in) :: tag
    integer, intent(in)          :: mode

    integer :: i
    real(kind=dp) :: t
    
    call cpu_time(t)

    select case (mode)

       case (1)

          do i=1,nnames
             if (clocks(i)%label .eq. tag) then
                clocks(i)%ptime  = t
                clocks(i)%ncalls = clocks(i)%ncalls + 1
                return
             endif
          enddo

          nnames = nnames + 1
          if (nnames.gt.nmax) call io_error('Maximum number of calls to io_stopwatch exceeded')

          clocks(nnames)%label = tag
          clocks(nnames)%ctime = 0.0_dp
          clocks(nnames)%ptime = t
          clocks(nnames)%ncalls = 1

       case (2)

          do i=1,nnames
             if (clocks(i)%label .eq. tag) then
                clocks(i)%ctime = clocks(i)%ctime + t - clocks(i)%ptime
                return
             endif
          end do

          write(stdout,'(1x,3a)') 'WARNING: name = ',trim(tag),' not found in io_stopwatch' 

       case default

          write(stdout,*) ' Name = ',trim(tag),' mode = ',mode
          call io_error('Value of mode not recognised in io_stopwatch')

    end select

    return

  end subroutine io_stopwatch


  !==================================================================!
  subroutine io_print_timings()
  !==================================================================!
  ! Output timing information to stdout
  !==================================================================!

    implicit none

    integer :: i

    write(stdout,'(/1x,a)') '*===========================================================================*'
    write(stdout,'(1x,a)')  '|                             TIMING INFORMATION                            |'
    write(stdout,'(1x,a)')  '*===========================================================================*'
    write(stdout,'(1x,a)')  '|    Tag                                                Ncalls      Time (s)|'    
    write(stdout,'(1x,a)')  '|---------------------------------------------------------------------------|'    
    do i=1,nnames
       write(stdout,'(1x,"|",a50,":",i10,4x,f10.3,"|")') &
            clocks(i)%label,clocks(i)%ncalls,clocks(i)%ctime 
    enddo
    write(stdout,'(1x,a)')  '*---------------------------------------------------------------------------*'
    
    return

  end subroutine io_print_timings


    !==================================================================!
       subroutine io_get_seedname (  )
    !==================================================================!
    !                                                                  !
    ! Get the seedname from the commandline                            !
    ! Note iargc and getarg are not standard                           !
    ! Some platforms require them to be external or provide            !
    ! equivalent routines. Not a problem in f2003!                     !
    !===================================================================  


#ifdef NAG
    USE F90_UNIX_ENV, ONLY : IARGC,GETARG
#endif

         implicit none

         integer :: num_arg
#ifndef NAG
    integer :: iargc
#endif
         character(len=50) :: ctemp 

         post_proc_flag=.false.

         num_arg=iargc()
         if (num_arg==0) then
            seedname='wannier'
         elseif (num_arg==1) then
            call getarg(1,seedname)
            if( index(seedname,'-pp')>0 ) then
               post_proc_flag=.true.
               seedname='wannier'
            end if
         else
            call getarg(1,seedname)
            if( index(seedname,'-pp')>0 ) then
               post_proc_flag=.true.
               call getarg(2,seedname)
            else
               call getarg(2,ctemp)
               if( index(ctemp,'-pp')>0 ) post_proc_flag=.true.
            end if

         end if

         ! If on the command line the whole seedname.win was passed, I strip the last ".win"
         if (len(trim(seedname)).ge.5) then
            if (seedname(len(trim(seedname))-4+1:).eq.".win") then
               seedname = seedname(:len(trim(seedname))-4)
            end if
         end if

       end subroutine io_get_seedname



    !==================================================================!
       subroutine io_error ( error_msg )
    !==================================================================!
    !                                                                  !
    ! Aborts giving error message                                      !
    !                                                                  !
    !===================================================================  


         implicit none
         character(len=*), intent(in) :: error_msg

#ifdef MPI
         character(len=50) :: filename
         integer           :: stderr,ierr,whoami

         call mpi_comm_rank(mpi_comm_world, whoami, ierr)
         if(whoami>99999) then
            write(filename,'(a,a,I0,a)')trim(seedname),'.node_',whoami,'.werr'
         else
            write(filename,'(a,a,I5.5,a)')trim(seedname),'.node_',whoami,'.werr'
         endif
         stderr=io_file_unit()
         open(unit=stderr,file=trim(filename),form='formatted',err=105)
         write(stderr, '(1x,a)') trim(error_msg)
         close(stderr)

105      write(*,'(1x,a)') trim(error_msg)
106      write(*,'(1x,a,I0,a)') "Error on node ", &
              whoami, ": examine the output/error files for details"
         
         call MPI_abort(MPI_comm_world,1,ierr)

#else

         write(stdout,*)  'Exiting.......' 
         write(stdout, '(1x,a)') trim(error_msg)
         
         close(stdout)
         
         write(*, '(1x,a)') trim(error_msg)
         write(*,'(A)') "Error: examine the output/error file for details" 
#endif

#ifdef EXIT_FLAG
         call exit(1)
#else
         STOP
#endif
         
       end subroutine io_error


    !==================================================================!
      subroutine io_date(cdate, ctime)
    !==================================================================!
    !                                                                  !
    !     Returns two strings containing the date and the time         !
    !     in human-readable format. Uses a standard f90 call.          !
    !                                                                  !
    !===================================================================  
    implicit none
    character (len=9), intent(out) :: cdate
    character (len=9), intent(out) :: ctime

    character(len=3), dimension(12) :: months
    data months /'Jan','Feb','Mar','Apr','May','Jun',   &
         'Jul','Aug','Sep','Oct','Nov','Dec'/
    integer date_time(8)
    !
    call date_and_time(values=date_time)
    !
    write (cdate,'(i2,a3,i4)') date_time(3), months(date_time(2)), date_time(1)
    write (ctime,'(i2.2,":",i2.2,":",i2.2)') date_time(5), date_time(6), date_time(7)

  end subroutine io_date


    !==================================================================!
      function io_time()
    !==================================================================!
    !                                                                  !
    ! Returns elapsed CPU time in seconds since its first call         !
    ! uses standard f90 call                                           !
    !                                                                  !
    !===================================================================  
    use w90_constants, only : dp
    implicit none

    real(kind=dp) :: io_time

    ! t0 contains the time of the first call
    ! t1 contains the present time
    real(kind=dp) :: t0, t1
    logical :: first=.true.
    save first, t0
    !
    call cpu_time(t1)
    !
    if (first) then
       t0 = t1
       io_time = 0.0_dp
       first = .false.
    else
       io_time = t1 - t0
    endif
    return
  end function io_time

  !==================================================================!
  function io_file_unit()
  !==================================================================!
  !                                                                  !
  ! Returns an unused unit number                                    !
  ! (so we can open a file on that unit                              !
  !                                                                  !
  !=================================================================== 
  implicit none

  integer :: io_file_unit,unit
  logical :: file_open

  unit = 9
  file_open = .true.
  do while ( file_open )
     unit = unit + 1
     inquire( unit, OPENED = file_open )
  end do

  io_file_unit = unit


  return
end function io_file_unit


end module w90_io
