!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!
!
! this routine contains the MPI parts of io_error() (called everywhere)
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

  subroutine comms_abort(seedname, error_msg, stdout)

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

    character(len=50), intent(in)  :: seedname
    character(len=*), intent(in) :: error_msg
    integer, intent(in) :: stdout

#ifdef MPI
    character(len=50) :: filename
    integer :: num_nodes, whoami, ierr, stderr

    ! this routine is not aware of any communicator other than WORLD
    call mpi_comm_size(MPI_COMM_WORLD, num_nodes, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, whoami, ierr)

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

    call mpi_abort(MPI_COMM_WORLD, 1, ierr)
#else
    write (stdout, *) 'Exiting.......'
    write (stdout, '(1x,a)') trim(error_msg)
#endif

  end subroutine comms_abort
