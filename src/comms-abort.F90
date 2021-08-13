module comms_abort_mod

  ! this module/routine corresponds to the MPI parts of io_error() (called everywhere)
  ! io_error previously called mpi_abort(mpi_comm_world) directly, but this means that the io
  ! module is dependent on the MPI compilation environment
  !
  ! it is anyway undesirable in a library to ever call mpi_abort (that is up to the calling program)
  ! and indeed the communicator may not be "world" as was previously assumed
  !
  ! routine here is temporary workaround until better error handling is introduced
  ! as soon as io_error no longer needs to stop/mpi_abort, this routine should be deleted, also
  !
  ! this needs to be it's own module because comms module already depends on io, therfore io cannot
  ! be made to depend on comms (otherwise module dependencies are circular)
  ! JJ Aug 2021

#ifdef MPI
#  if !(defined(MPI08) || defined(MPI90) || defined(MPIH))
#    error "You need to define which MPI interface you are using"
#  endif
#endif

#ifdef MPI08
  use mpi_f08 ! use f08 interface if possible
#endif
#ifdef MPI90
  use mpi ! next best, use fortran90 interface
#endif

  implicit none

#ifdef MPIH
  include 'mpif.h' ! worst case, use legacy interface
#endif

  type :: w90commtype
#ifdef MPI08
    type(mpi_comm) :: comm ! f08 mpi interface
#else
    integer :: comm ! f90 mpi or no mpi
#endif
  end type

  public :: comms_abort

contains
  subroutine comms_abort(seedname, error_msg, stdout)

    character(len=50), intent(in)  :: seedname
    character(len=*), intent(in) :: error_msg
    integer, intent(in) :: stdout

#ifdef MPI

    character(len=50) :: filename
    type(w90commtype) :: w90comm
    integer :: num_nodes, whoami, ierr, stderr

    w90comm%comm = MPI_COMM_WORLD !JJ not necessarily comm_world!!
    call mpi_comm_size(w90comm%comm, num_nodes, ierr)
    call mpi_comm_rank(w90comm%comm, whoami, ierr)

    if (num_nodes > 1) then
      if (whoami > 99999) then
        write (filename, '(a,a,I0,a)') trim(seedname), '.node_', whoami, '.werr'
      else
        write (filename, '(a,a,I5.5,a)') trim(seedname), '.node_', whoami, '.werr'
      endif
      open (newunit=stderr, file=trim(filename), form='formatted', err=105)
      write (stderr, '(1x,a)') trim(error_msg)
      close (stderr)
    end if
105 write (*, '(1x,a)') trim(error_msg)

    if (whoami == 0) then
      write (stdout, *) 'Exiting.......'
      write (stdout, '(1x,a)') trim(error_msg)
      close (stdout)
    end if

    call mpi_abort(w90comm%comm, 1, ierr)
#else
    write (stdout, *) 'Exiting.......'
    write (stdout, '(1x,a)') trim(error_msg)
#endif

  end subroutine comms_abort
end module comms_abort_mod

