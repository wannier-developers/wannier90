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
!                                                            !
!  COMMS: set of MPI wrappers                                !
!  written 2006-2012 Jonathan R. Yates                       !
!    later additions Giovanni Pizzi                          !
!                                                            !
!------------------------------------------------------------!

module w90_comms
  !! This module handles all of the communications

  use w90_constants, only: dp
  use w90_error_base

#ifdef MPI
#  if !(defined(MPI08) || defined(MPI90) || defined(MPIH))
#    error "You need to define which MPI interface you are using"
#  endif
#else
#define MPI_COMM_NULL -1
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

  private

  integer, parameter :: code_mpi = 4 ! this is duplicated here to avoid circular dependency with error.F90

  integer, parameter :: mpi_send_tag = 77 ! arbitrary
  integer, parameter :: root_id = 0 ! not arbitrary

  public :: comms_allreduce  ! reduce data onto all nodes
  public :: comms_array_split
  public :: comms_barrier    ! puts a barrier so that the code goes on only when all nodes reach the barrier
  public :: comms_bcast      ! send data from the root node
  public :: comms_gatherv    ! gets chunks of an array from all nodes and gathers them on the root node
  !public :: comms_recv       ! accept data from one node to another
  public :: comms_reduce     ! reduce data onto root node (n.b. not allreduce); data is lost on all other nodes
  public :: comms_scatterv   ! sends chunks of an array to all nodes scattering them from the root node
  !public :: comms_send       ! send data from one node to another
  public :: mpirank
  public :: mpisize
  public :: comms_sync_error
  public :: valid_communicator
  !! test whether communicator is initialised; returns true always for serial build

  ! versions without error synchronisation, use at own risk
  public :: comms_no_sync_allreduce  ! reduce data onto all nodes
  public :: comms_no_sync_barrier    ! puts a barrier so that the code goes on only when all nodes reach the barrier
  public :: comms_no_sync_bcast      ! send data from the root node
  public :: comms_no_sync_gatherv    ! gets chunks of an array from all nodes and gathers them on the root node
  public :: comms_no_sync_recv       ! accept data from one node to another
  public :: comms_no_sync_reduce     ! reduce data onto root node (n.b. not allreduce); data is lost on all other nodes
  public :: comms_no_sync_scatterv   ! sends chunks of an array to all nodes scattering them from the root node
  public :: comms_no_sync_send       ! send data from one node to another

  type, public :: w90_comm_type
#ifdef MPI08
    type(mpi_comm) :: comm = MPI_COMM_NULL ! f08 mpi interface
#else
    integer :: comm = MPI_COMM_NULL ! f90 mpi or no mpi
#endif
  end type

  type, private :: w90stat_type
#ifdef MPI08
    type(mpi_status) :: stat ! f08 mpi interface
#elif MPI90
    integer :: stat(MPI_STATUS_SIZE)
#else
    integer :: stat ! not used
#endif
  end type

  interface comms_bcast
    module procedure comms_bcast_int
    module procedure comms_bcast_logical
    module procedure comms_bcast_real
    module procedure comms_bcast_cmplx
    module procedure comms_bcast_char
  end interface comms_bcast

  ! We don't currently have error synchronisation for point-to-point comms
  !interface comms_send
  !  module procedure comms_send_int
  !  module procedure comms_send_logical
  !  module procedure comms_send_real
  !  module procedure comms_send_cmplx
  !  module procedure comms_send_char
  !end interface comms_send

  !interface comms_recv
  !  module procedure comms_recv_int
  !  module procedure comms_recv_logical
  !  module procedure comms_recv_real
  !  module procedure comms_recv_cmplx
  !  module procedure comms_recv_char
  !end interface comms_recv

  interface comms_reduce
    module procedure comms_reduce_int
    module procedure comms_reduce_real
    module procedure comms_reduce_cmplx
  end interface comms_reduce

  interface comms_allreduce
!     module procedure comms_allreduce_int    ! to be done
    module procedure comms_allreduce_real
    module procedure comms_allreduce_cmplx
  end interface comms_allreduce

  interface comms_gatherv
!     module procedure comms_gatherv_int    ! to be done
    module procedure comms_gatherv_logical
    module procedure comms_gatherv_real_1
    module procedure comms_gatherv_real_2
    module procedure comms_gatherv_real_3
    module procedure comms_gatherv_real_2_3
    module procedure comms_gatherv_cmplx_1
    module procedure comms_gatherv_cmplx_2
    module procedure comms_gatherv_cmplx_3
    module procedure comms_gatherv_cmplx_3_4
    module procedure comms_gatherv_cmplx_4
  end interface comms_gatherv

  interface comms_scatterv
    module procedure comms_scatterv_int_1
    module procedure comms_scatterv_int_2
    module procedure comms_scatterv_int_3
    module procedure comms_scatterv_real_1
    module procedure comms_scatterv_real_2
    module procedure comms_scatterv_real_3
!     module procedure comms_scatterv_cmplx
    module procedure comms_scatterv_cmplx_3
    module procedure comms_scatterv_cmplx_4
  end interface comms_scatterv

  !! These routines do not synchronise errors. You should not normally use them.
  !! If you do use them, you are required to make sure that the error conditions
  !! between the MPI processes are correctly synced up on exit from your code.
  interface comms_no_sync_bcast
    module procedure comms_no_sync_bcast_int
    module procedure comms_no_sync_bcast_logical
    module procedure comms_no_sync_bcast_real
    module procedure comms_no_sync_bcast_cmplx
    module procedure comms_no_sync_bcast_char
  end interface comms_no_sync_bcast

  interface comms_no_sync_send
    module procedure comms_no_sync_send_int
    module procedure comms_no_sync_send_logical
    module procedure comms_no_sync_send_real
    module procedure comms_no_sync_send_cmplx
    module procedure comms_no_sync_send_char
  end interface comms_no_sync_send

  interface comms_no_sync_recv
    module procedure comms_no_sync_recv_int
    module procedure comms_no_sync_recv_logical
    module procedure comms_no_sync_recv_real
    module procedure comms_no_sync_recv_cmplx
    module procedure comms_no_sync_recv_char
  end interface comms_no_sync_recv

  interface comms_no_sync_reduce
    module procedure comms_no_sync_reduce_int
    module procedure comms_no_sync_reduce_real
    module procedure comms_no_sync_reduce_cmplx
  end interface comms_no_sync_reduce

  interface comms_no_sync_allreduce
!     module procedure comms_no_sync_allreduce_int    ! to be done
    module procedure comms_no_sync_allreduce_real
    module procedure comms_no_sync_allreduce_cmplx
  end interface comms_no_sync_allreduce

  interface comms_no_sync_gatherv
!     module procedure comms_no_sync_gatherv_int    ! to be done
    module procedure comms_no_sync_gatherv_logical
    module procedure comms_no_sync_gatherv_real_1
    module procedure comms_no_sync_gatherv_real_2
    module procedure comms_no_sync_gatherv_real_3
    module procedure comms_no_sync_gatherv_real_2_3
    module procedure comms_no_sync_gatherv_cmplx_1
    module procedure comms_no_sync_gatherv_cmplx_2
    module procedure comms_no_sync_gatherv_cmplx_3
    module procedure comms_no_sync_gatherv_cmplx_3_4
    module procedure comms_no_sync_gatherv_cmplx_4
  end interface comms_no_sync_gatherv

  interface comms_no_sync_scatterv
    module procedure comms_no_sync_scatterv_int_1
    module procedure comms_no_sync_scatterv_int_2
    module procedure comms_no_sync_scatterv_int_3
    module procedure comms_no_sync_scatterv_real_1
    module procedure comms_no_sync_scatterv_real_2
    module procedure comms_no_sync_scatterv_real_3
!     module procedure comms_no_sync_scatterv_cmplx
    module procedure comms_no_sync_scatterv_cmplx_3
    module procedure comms_no_sync_scatterv_cmplx_4
  end interface comms_no_sync_scatterv

contains

  logical function valid_communicator(comm)
    type(w90_comm_type), intent(in) :: comm
#ifdef MPI
    if (comm%comm == MPI_COMM_NULL) then
      valid_communicator = .false.
    else
      valid_communicator = .true.
    endif
#else
    valid_communicator = .true.
#endif
  end function

  ! mpi rank function for convenience
  integer function mpirank(comm)
    type(w90_comm_type), intent(in) :: comm
    integer :: ierr
#ifdef MPI
    call mpi_comm_rank(comm%comm, mpirank, ierr)
#else
    mpirank = 0
#endif
  end function

  ! mpi size function for convenience
  integer function mpisize(comm)
    type(w90_comm_type), intent(in) :: comm
    integer :: ierr
#ifdef MPI
    call mpi_comm_size(comm%comm, mpisize, ierr)
#else
    mpisize = 1
#endif
  end function

  ! synchronise error condition between MPI processess
  subroutine comms_sync_error(comm, error, ierr)
    implicit none
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(inout) :: error
    integer :: ierr, mpiierr, abserr

#if defined(MPI) && !defined(DISABLE_ERROR_SYNC)
    abserr = abs(ierr) ! possibility of -ve values, use abs for safety
    call mpi_allreduce(MPI_IN_PLACE, abserr, 1, MPI_INTEGER, MPI_SUM, comm%comm, mpiierr)
    ! you could check mpiierr here, but truly all bets are off in that case
    if (abserr > 0 .and. ierr == 0) call recv_error(error)
#endif
  end subroutine comms_sync_error

  subroutine comms_array_split(numpoints, counts, displs, comm)
    !! Given an array of size numpoints, we want to split on num_nodes nodes. This function returns
    !! two arrays: count and displs.
    !!
    !! The i-th element of the count array gives the number of elements
    !! that must be calculated by the process with id (i-1).
    !! The i-th element of the displs array gives the displacement of the array calculated locally on
    !! the process with id (i-1) with respect to the global array.
    !!
    !! These values are those to be passed to the functions MPI_Scatterv, MPI_Gatherv and MPI_Alltoallv.
    !!
    !! one can use the following do loop to run over the needed elements, if the full array is stored
    !! on all nodes:
    !! do i=displs(my_node_id)+1,displs(my_node_id)+counts(my_node_id)
    !!

    integer, intent(in) :: numpoints  !! Number of elements of the array to be scattered
    integer, intent(inout) :: counts(0:) !! Array (of size num_nodes) with the number of elements of the array on each node
    integer, intent(inout) :: displs(0:) !! Array (of size num_nodes) with the displacement relative to the global array
    type(w90_comm_type), intent(in) :: comm

    integer :: ratio, remainder, i
    integer :: num_nodes

    num_nodes = mpisize(comm)

    ratio = numpoints/num_nodes
    remainder = MOD(numpoints, num_nodes)

    do i = 0, num_nodes - 1
      if (i < remainder) then
        counts(i) = ratio + 1
        displs(i) = i*(ratio + 1)
      else
        counts(i) = ratio
        displs(i) = remainder*(ratio + 1) + (i - remainder)*ratio
      end if
    end do

  end subroutine comms_array_split

  subroutine comms_no_sync_barrier(comm)
    !! A barrier to synchronise all nodes
    implicit none
    type(w90_comm_type), intent(in) :: comm

#ifdef MPI
    integer :: ierr

    call mpi_barrier(comm%comm, ierr)
#endif

  end subroutine comms_no_sync_barrier

  subroutine comms_no_sync_bcast_int(array, size, error, comm)
    !! Send integar array from root node to all nodes
    implicit none

    integer, intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_bcast(array, size, MPI_INTEGER, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_bcast_int', code_mpi)
      return
    end if
#endif
  end subroutine comms_no_sync_bcast_int

  subroutine comms_no_sync_bcast_real(array, size, error, comm)
    !! Send real array from root node to all nodes
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_bcast(array, size, MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_bcast_real', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_bcast_real

  subroutine comms_no_sync_bcast_logical(array, size, error, comm)
    !! Send logical array from root node to all nodes
    implicit none

    logical, intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_bcast(array, size, MPI_LOGICAL, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_bcast_logical', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_bcast_logical

  subroutine comms_no_sync_bcast_char(array, size, error, comm)
    !! Send character array from root node to all nodes
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_bcast(array, size, MPI_CHARACTER, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_bcast_char', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_bcast_char

  subroutine comms_no_sync_bcast_cmplx(array, size, error, comm)
    !! Send character array from root node to all nodes

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_bcast(array, size, MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_bcast_cmplx', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_bcast_cmplx

  !--------- SEND ----------------

  subroutine comms_no_sync_send_logical(array, size, to, error, comm)
    !! Send logical array to specified node

    implicit none

    logical, intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: to
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_send(array, size, MPI_LOGICAL, to, mpi_send_tag, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_send_logical', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_send_logical

  subroutine comms_no_sync_send_int(array, size, to, error, comm)
    !! Send integer array to specified node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: to
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_send(array, size, MPI_INTEGER, to, mpi_send_tag, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_send_int', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_send_int

  subroutine comms_no_sync_send_char(array, size, to, error, comm)
    !! Send character array to specified node
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: to
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_send(array, size, MPI_CHARACTER, to, mpi_send_tag, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_send_char', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_send_char

  subroutine comms_no_sync_send_real(array, size, to, error, comm)
    !! Send real array to specified node
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: to
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_send(array, size, MPI_DOUBLE_PRECISION, to, mpi_send_tag, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_send_real', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_send_real

  subroutine comms_no_sync_send_cmplx(array, size, to, error, comm)
    !! Send complex array to specified node
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: to
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_send(array, size, MPI_DOUBLE_COMPLEX, to, mpi_send_tag, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_send_cmplx', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_send_cmplx

  !--------- RECV ----------------

  subroutine comms_no_sync_recv_logical(array, size, from, error, comm)
    !! Receive logical array from specified node
    implicit none

    logical, intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: from
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    type(w90stat_type) :: status
    integer :: ierr

    call mpi_recv(array, size, MPI_LOGICAL, from, mpi_send_tag, comm%comm, status%stat, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_recv_logical', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_recv_logical

  subroutine comms_no_sync_recv_int(array, size, from, error, comm)
    !! Receive integer array from specified node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: from
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    type(w90stat_type) :: status
    integer :: ierr

    call mpi_recv(array, size, MPI_INTEGER, from, mpi_send_tag, comm%comm, status%stat, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_recv_int', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_recv_int

  subroutine comms_no_sync_recv_char(array, size, from, error, comm)
    !! Receive character array from specified node
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: from
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    type(w90stat_type) :: status
    integer :: ierr

    call mpi_recv(array, size, MPI_CHARACTER, from, mpi_send_tag, comm%comm, status%stat, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_recv_char', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_recv_char

  subroutine comms_no_sync_recv_real(array, size, from, error, comm)
    !! Receive real array from specified node
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: from
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    type(w90stat_type) :: status
    integer :: ierr

    call mpi_recv(array, size, MPI_DOUBLE_PRECISION, from, mpi_send_tag, comm%comm, &
                  status%stat, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_recv_real', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_recv_real

  subroutine comms_no_sync_recv_cmplx(array, size, from, error, comm)
    !! Receive complex array from specified node
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    integer, intent(in) :: from
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    type(w90stat_type) :: status
    integer :: ierr

    call mpi_recv(array, size, MPI_DOUBLE_COMPLEX, from, mpi_send_tag, comm%comm, &
                  status%stat, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_recv_cmplx', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_recv_cmplx

  subroutine comms_no_sync_reduce_int(array, size, op, error, comm)
    !! Reduce integer data to root node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr
    integer :: rank
    rank = mpirank(comm)

    ! note, JJ 23/2/2021
    ! previously this routine alloc'd/used/dealloc'd a temp array
    ! to be used as receive buffer for MPI_reduce
    ! this temp array was then copied to argument "array"
    ! but: "array" needs to be of scalar type for the polymorphism to work
    ! so: need to copy array into a (fake) scalar
    ! previously: a subroutine my_icopy was used to help to do this.
    ! probably just reducing in place is better?

    select case (op)
    case ('SUM')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_INTEGER, MPI_SUM, root_id, comm%comm, &
                        ierr)
      else
        call mpi_reduce(array, array, size, MPI_INTEGER, MPI_SUM, root_id, comm%comm, ierr)
      endif
    case ('PRD')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_INTEGER, MPI_PROD, root_id, comm%comm, &
                        ierr)
      else
        call mpi_reduce(array, array, size, MPI_INTEGER, MPI_PROD, root_id, comm%comm, ierr)
      endif
    case default
      call set_base_error(error, 'Unknown operation in comms_reduce_int', code_mpi)
      return
    end select

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_reduce_int', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_reduce_int

  subroutine comms_no_sync_reduce_real(array, size, op, error, comm)
    !! Reduce real data to root node

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr
    integer :: rank
    rank = mpirank(comm)

    select case (op)

    case ('SUM')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_SUM, root_id, &
                        comm%comm, ierr)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_SUM, root_id, comm%comm, &
                        ierr)
      endif
    case ('PRD')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_PROD, root_id, &
                        comm%comm, ierr)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_PROD, root_id, comm%comm, &
                        ierr)
      endif
    case ('MIN')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MIN, root_id, &
                        comm%comm, ierr)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_MIN, root_id, comm%comm, &
                        ierr)
      endif
    case ('MAX')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MAX, root_id, &
                        comm%comm, ierr)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_MAX, root_id, comm%comm, &
                        ierr)
      endif
    case default
      call set_base_error(error, 'Unknown operation in comms_reduce_real', code_mpi)
      return

    end select

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_reduce_real', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_reduce_real

  subroutine comms_no_sync_reduce_cmplx(array, size, op, error, comm)
    !! Reduce complex data to root node

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr
    integer :: rank
    rank = mpirank(comm)

    select case (op)

    case ('SUM')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_SUM, root_id, &
                        comm%comm, ierr)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_COMPLEX, MPI_SUM, root_id, comm%comm, &
                        ierr)
      end if
    case ('PRD')
      if (rank == root_id) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_PROD, root_id, &
                        comm%comm, ierr)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_COMPLEX, MPI_PROD, root_id, comm%comm, &
                        ierr)
      end if
    case default
      call set_base_error(error, 'Unknown operation in comms_reduce_cmplx', code_mpi)
      return

    end select

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_reduce_cmplx', code_mpi)
      return
    end if

#endif

  end subroutine comms_no_sync_reduce_cmplx

  subroutine comms_no_sync_allreduce_real(array, size, op, error, comm)
    !! Reduce real data to all nodes

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    select case (op)

    case ('SUM')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_SUM, comm%comm, &
                         ierr)
    case ('PRD')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_PROD, comm%comm, &
                         ierr)
    case ('MIN')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MIN, comm%comm, &
                         ierr)
    case ('MAX')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MAX, comm%comm, &
                         ierr)
    case default
      call set_base_error(error, 'Unknown operation in comms_allreduce_real', code_mpi)
      return

    end select

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_allreduce_real', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_allreduce_real

  subroutine comms_no_sync_allreduce_cmplx(array, size, op, error, comm)
    !! Reduce complex data to all nodes
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    select case (op)

    case ('SUM')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_SUM, comm%comm, &
                         ierr)
    case ('PRD')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_PROD, comm%comm, &
                         ierr)
    case default
      call set_base_error(error, 'Unknown operation in comms_allreduce_cmplx', code_mpi)
      return

    end select

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_allreduce_cmplx', code_mpi)
      return
    end if
#endif

  end subroutine comms_no_sync_allreduce_cmplx

  subroutine comms_no_sync_gatherv_real_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 1)
    implicit none

    real(kind=dp), intent(inout) :: array(:) !! local array for sending data
    integer, intent(in) :: localcount !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, &
                     displs, MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_real_1', code_mpi)
      return
    end if
#else
    !call dcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif

  end subroutine comms_no_sync_gatherv_real_1

  subroutine comms_no_sync_gatherv_real_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 2)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :)           !! local array for sending data
    integer, intent(in) :: localcount                     !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:)                     !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, &
                     displs, MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_real_2', code_mpi)
      return
    end if
#else
    !call dcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif

  end subroutine comms_no_sync_gatherv_real_2

  subroutine comms_no_sync_gatherv_real_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 3)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :, :)           !! local array for sending data
    integer, intent(in) :: localcount                        !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:)                         !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, &
                     displs, MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_real_3', code_mpi)
      return
    end if

#else
    !call dcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif

  end subroutine comms_no_sync_gatherv_real_3

  subroutine comms_no_sync_gatherv_real_2_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 2 and 3, respectively)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :)              !! local array for sending data
    integer, intent(in) :: localcount                        !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:)                         !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_real_2_3', code_mpi)
      return
    end if

#else
    call dcopy(localcount, array, 1, rootglobalarray, 1)
    !rootglobalarray = array ! shapes don't match
#endif

  end subroutine comms_no_sync_gatherv_real_2_3

  ! Array: local array for sending data; localcount elements will be sent
  !        to the root node
  ! rootglobalarray: array on the root node to which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split

  subroutine comms_no_sync_gatherv_cmplx_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 1)
    implicit none

    complex(kind=dp), intent(inout) :: array(:)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_cmplx_1', code_mpi)
      return
    end if

#else
    !call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif

  end subroutine comms_no_sync_gatherv_cmplx_1

  subroutine comms_no_sync_gatherv_cmplx_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 2)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_cmplx_2', code_mpi)
      return
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif

  end subroutine comms_no_sync_gatherv_cmplx_2

  subroutine comms_no_sync_gatherv_cmplx_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 3)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_cmplx_3', code_mpi)
      return
    end if

#else
    !call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif

  end subroutine comms_no_sync_gatherv_cmplx_3

  subroutine comms_no_sync_gatherv_cmplx_3_4(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 3 and 4, respectively)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_cmplx_3_4', code_mpi)
      return
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
    !rootglobalarray = array ! shapes don't match
#endif

  end subroutine comms_no_sync_gatherv_cmplx_3_4

  subroutine comms_no_sync_gatherv_cmplx_4(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 4)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_cmplx_4', code_mpi)
      return
    end if

#else
    !call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif

  end subroutine comms_no_sync_gatherv_cmplx_4

  subroutine comms_no_sync_gatherv_logical(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node
    implicit none

    logical, intent(inout) :: array !! local array for sending data
    integer, intent(in) :: localcount !! localcount elements will be sent to the root node
    logical, intent(inout) :: rootglobalarray !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_gatherv(array, localcount, MPI_LOGICAL, rootglobalarray, counts, displs, &
                     MPI_LOGICAL, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_gatherv_logical', code_mpi)
      return
    end if
#else
!    rootglobalarray(1:localcount)=array(1:localcount)
#endif

  end subroutine comms_no_sync_gatherv_logical

  subroutine comms_no_sync_scatterv_real_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter real data from root node (array of rank 1)
    implicit none

    real(kind=dp), intent(inout) :: array(:) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_PRECISION, array, localcount, &
                      MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_real_1', code_mpi)
      return
    end if

#else
    !call dcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_real_1

  subroutine comms_no_sync_scatterv_real_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter real data from root node (array of rank 2)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_PRECISION, array, localcount, &
                      MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_real_2', code_mpi)
      return
    end if

#else
    !call dcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_real_2

  subroutine comms_no_sync_scatterv_real_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter real data from root node (array of rank 3)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_PRECISION, array, localcount, &
                      MPI_DOUBLE_PRECISION, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_real_3', code_mpi)
      return
    end if

#else
    !call dcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_real_3

  subroutine comms_no_sync_scatterv_cmplx_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter complex data from root node (array of rank 3)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_COMPLEX, array, localcount, &
                      MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_cmplx_3', code_mpi)
      return
    end if

#else
    !call zcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_cmplx_3

  subroutine comms_no_sync_scatterv_cmplx_4(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter complex data from root node (array of rank 4)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_COMPLEX, array, localcount, &
                      MPI_DOUBLE_COMPLEX, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_cmplx_4', code_mpi)
      return
    end if

#else
    !call zcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_cmplx_4

  subroutine comms_no_sync_scatterv_int_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter integer data from root node (array of rank 1)
    implicit none

    integer, intent(inout) :: array(:) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    integer, intent(inout) :: rootglobalarray(:) !!  array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_INTEGER, array, localcount, &
                      MPI_INTEGER, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_real', code_mpi)
      return
    end if

#else
    !call my_icopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_int_1

  subroutine comms_no_sync_scatterv_int_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter integer data from root node (array of rank 2)

    implicit none

    integer, intent(inout) :: array(:, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    integer, intent(inout) :: rootglobalarray(:, :) !!  array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_INTEGER, array, localcount, &
                      MPI_INTEGER, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_int_2', code_mpi)
      return
    end if

#else
    !call my_icopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_int_2

  subroutine comms_no_sync_scatterv_int_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter integer data from root node (array of rank 3)

    implicit none

    integer, intent(inout) :: array(:, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    integer, intent(inout) :: rootglobalarray(:, :, :) !!  array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

#ifdef MPI
    integer :: ierr

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_INTEGER, array, localcount, &
                      MPI_INTEGER, root_id, comm%comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call set_base_error(error, 'Error in comms_scatterv_int_3', code_mpi)
      return
    end if

#else
    !call my_icopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif

  end subroutine comms_no_sync_scatterv_int_3

  ! The routines that synchronise the errors
  subroutine comms_barrier(error, comm)
    !! A barrier to synchronise all nodes
    implicit none
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    call comms_sync_error(comm, error, 0)
    if (allocated(error)) return
    ! The barrier is almost redundant since the sync is global, unless DISABLE_ERROR_SYNC defined
    call comms_no_sync_barrier(comm)
  end subroutine comms_barrier

  subroutine comms_bcast_int(array, size, error, comm)
    !! Send integar array from root node to all nodes
    implicit none

    integer, intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0)
    if (allocated(error)) return

    call comms_no_sync_bcast_int(array, size, error, comm)
  end subroutine comms_bcast_int

  subroutine comms_bcast_real(array, size, error, comm)
    !! Send real array from root node to all nodes
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_bcast_real(array, size, error, comm)
  end subroutine comms_bcast_real

  subroutine comms_bcast_logical(array, size, error, comm)
    !! Send logical array from root node to all nodes
    implicit none

    logical, intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_bcast_logical(array, size, error, comm)
  end subroutine comms_bcast_logical

  subroutine comms_bcast_char(array, size, error, comm)
    !! Send character array from root node to all nodes
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_bcast_char(array, size, error, comm)
  end subroutine comms_bcast_char

  subroutine comms_bcast_cmplx(array, size, error, comm)
    !! Send character array from root node to all nodes

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_bcast_cmplx(array, size, error, comm)
  end subroutine comms_bcast_cmplx

  subroutine comms_reduce_int(array, size, op, error, comm)
    !! Reduce integer data to root node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_reduce_int(array, size, op, error, comm)
  end subroutine comms_reduce_int

  subroutine comms_reduce_real(array, size, op, error, comm)
    !! Reduce real data to root node

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_reduce_real(array, size, op, error, comm)
  end subroutine comms_reduce_real

  subroutine comms_reduce_cmplx(array, size, op, error, comm)
    !! Reduce complex data to root node

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_reduce_cmplx(array, size, op, error, comm)
  end subroutine comms_reduce_cmplx

  subroutine comms_allreduce_real(array, size, op, error, comm)
    !! Reduce real data to all nodes

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_allreduce_real(array, size, op, error, comm)
  end subroutine comms_allreduce_real

  subroutine comms_allreduce_cmplx(array, size, op, error, comm)
    !! Reduce complex data to all nodes
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in) :: size
    character(len=*), intent(in) :: op
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_allreduce_cmplx(array, size, op, error, comm)
  end subroutine comms_allreduce_cmplx

  subroutine comms_gatherv_real_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 1)
    implicit none

    real(kind=dp), intent(inout) :: array(:) !! local array for sending data
    integer, intent(in) :: localcount !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_real_1(array, localcount, rootglobalarray, counts, displs, &
                                      error, comm)
  end subroutine comms_gatherv_real_1

  subroutine comms_gatherv_real_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 2)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :)           !! local array for sending data
    integer, intent(in) :: localcount                     !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:)                     !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_real_2(array, localcount, rootglobalarray, counts, displs, &
                                      error, comm)
  end subroutine comms_gatherv_real_2

  subroutine comms_gatherv_real_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 3)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :, :)           !! local array for sending data
    integer, intent(in) :: localcount                        !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:)                         !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_real_3(array, localcount, rootglobalarray, counts, displs, &
                                      error, comm)
  end subroutine comms_gatherv_real_3

  subroutine comms_gatherv_real_2_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node (for arrays of rank 2 and 3, respectively)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :)              !! local array for sending data
    integer, intent(in) :: localcount                        !! localcount elements will be sent to the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:)                         !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_real_2_3(array, localcount, rootglobalarray, counts, displs, &
                                        error, comm)
  end subroutine comms_gatherv_real_2_3

  ! Array: local array for sending data; localcount elements will be sent
  !        to the root node
  ! rootglobalarray: array on the root node to which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split

  subroutine comms_gatherv_cmplx_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 1)
    implicit none

    complex(kind=dp), intent(inout) :: array(:)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_cmplx_1(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_gatherv_cmplx_1

  subroutine comms_gatherv_cmplx_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 2)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_cmplx_2(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_gatherv_cmplx_2

  subroutine comms_gatherv_cmplx_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 3)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_cmplx_3(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_gatherv_cmplx_3

  subroutine comms_gatherv_cmplx_3_4(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 3 and 4, respectively)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_cmplx_3_4(array, localcount, rootglobalarray, counts, displs, &
                                         error, comm)
  end subroutine comms_gatherv_cmplx_3_4

  subroutine comms_gatherv_cmplx_4(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather complex data to root node (for arrays of rank 4)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :, :)
    integer, intent(in) :: localcount
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :, :)
    integer, intent(in) :: counts(0:)
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_cmplx_4(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_gatherv_cmplx_4

  subroutine comms_gatherv_logical(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Gather real data to root node
    implicit none

    logical, intent(inout) :: array !! local array for sending data
    integer, intent(in) :: localcount !! localcount elements will be sent to the root node
    logical, intent(inout) :: rootglobalarray !! array on the root node to which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_gatherv_logical(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_gatherv_logical

  subroutine comms_scatterv_real_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter real data from root node (array of rank 1)
    implicit none

    real(kind=dp), intent(inout) :: array(:) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_real_1(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_scatterv_real_1

  subroutine comms_scatterv_real_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter real data from root node (array of rank 2)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_real_2(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_scatterv_real_2

  subroutine comms_scatterv_real_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter real data from root node (array of rank 3)
    implicit none

    real(kind=dp), intent(inout) :: array(:, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    real(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_real_3(array, localcount, rootglobalarray, counts, displs, &
                                       error, comm)
  end subroutine comms_scatterv_real_3

  subroutine comms_scatterv_cmplx_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter complex data from root node (array of rank 3)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_cmplx_3(array, localcount, rootglobalarray, counts, displs, &
                                        error, comm)
  end subroutine comms_scatterv_cmplx_3

  subroutine comms_scatterv_cmplx_4(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter complex data from root node (array of rank 4)
    implicit none

    complex(kind=dp), intent(inout) :: array(:, :, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    complex(kind=dp), intent(inout) :: rootglobalarray(:, :, :, :) !! array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_cmplx_4(array, localcount, rootglobalarray, counts, displs, &
                                        error, comm)
  end subroutine comms_scatterv_cmplx_4

  subroutine comms_scatterv_int_1(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter integer data from root node (array of rank 1)
    implicit none

    integer, intent(inout) :: array(:) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    integer, intent(inout) :: rootglobalarray(:) !!  array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_int_1(array, localcount, rootglobalarray, counts, displs, &
                                      error, comm)
  end subroutine comms_scatterv_int_1

  subroutine comms_scatterv_int_2(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter integer data from root node (array of rank 2)

    implicit none

    integer, intent(inout) :: array(:, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    integer, intent(inout) :: rootglobalarray(:, :) !!  array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_int_2(array, localcount, rootglobalarray, counts, displs, &
                                      error, comm)
  end subroutine comms_scatterv_int_2

  subroutine comms_scatterv_int_3(array, localcount, rootglobalarray, counts, displs, error, comm)
    !! Scatter integer data from root node (array of rank 3)

    implicit none

    integer, intent(inout) :: array(:, :, :) !! local array for getting data
    integer, intent(in) :: localcount !! localcount elements will be fetched from the root node
    integer, intent(inout) :: rootglobalarray(:, :, :) !!  array on the root node from which data will be sent
    integer, intent(in) :: counts(0:) !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, intent(in) :: displs(0:)
    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    call comms_sync_error(comm, error, 0) ! sync error state across comm
    if (allocated(error)) return

    call comms_no_sync_scatterv_int_3(array, localcount, rootglobalarray, counts, displs, &
                                      error, comm)
  end subroutine comms_scatterv_int_3
end module w90_comms
