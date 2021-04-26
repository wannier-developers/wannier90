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
!  COMMS: set of MPI wrapper                                 !
!  written 2006-2012 Jonathan R. Yates                       !
!    later additions Giovanni Pizzi                          !
!                                                            !
!------------------------------------------------------------!

! JJ 04.21 use MPI_IN_PLACE to avoid unnecessary array alloc/copy/deallocs

module w90_comms
  !! This module handles all of the communications

  use w90_constants, only: dp
  use w90_io, only: io_error

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

  logical, public, save :: on_root
  !! Are we the root node
  integer, public, save :: num_nodes
  !! Number of nodes
  integer, public, save :: my_node_id
  !! ID of this node
  integer, public, parameter :: root_id = 0
  !! ID of the root node

  integer, parameter :: mpi_send_tag = 77 !abitrary

  public :: comms_allreduce  ! reduce data onto all nodes
  public :: comms_array_split
  public :: comms_barrier    ! puts a barrier so that the code goes on only when all nodes reach the barrier
  public :: comms_bcast      ! send data from the root node
  public :: comms_end
  public :: comms_gatherv    ! gets chunks of an array from all nodes and gathers them on the root node
  public :: comms_recv       ! accept data from one node to another
  public :: comms_reduce     ! reduce data onto root node (n.b. not allreduce); data is lost on all other nodes
  public :: comms_scatterv   ! sends chunks of an array to all nodes scattering them from the root node
  public :: comms_send       ! send data from one node to another
  public :: comms_setup
  public :: comms_setup_vars

  type, public :: w90commtype
#ifdef MPI08
    type(mpi_comm) :: comm ! f08 mpi interface
#else
    integer :: comm ! f90 mpi or no mpi
#endif
  end type

  interface comms_bcast
    module procedure comms_bcast_int
    module procedure comms_bcast_logical
    module procedure comms_bcast_real
    module procedure comms_bcast_cmplx
    module procedure comms_bcast_char
  end interface comms_bcast

  interface comms_send
    module procedure comms_send_int
    module procedure comms_send_logical
    module procedure comms_send_real
    module procedure comms_send_cmplx
    module procedure comms_send_char
  end interface comms_send

  interface comms_recv
    module procedure comms_recv_int
    module procedure comms_recv_logical
    module procedure comms_recv_real
    module procedure comms_recv_cmplx
    module procedure comms_recv_char
  end interface comms_recv

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
    module procedure comms_scatterv_cmplx_4
  end interface comms_scatterv

contains

  subroutine comms_setup(w90comm)
    !! Set up communications
    !! comms_setup only invoked by main routine in executables wannier90.x and postw90.x
    implicit none
    type(w90commtype), intent(inout) :: w90comm

#ifdef MPI
    integer :: ierr

    w90comm%comm = MPI_COMM_WORLD
    call mpi_init(ierr)
    if (ierr .ne. 0) call io_error('MPI initialisation error')
#endif

    call comms_setup_vars(w90comm)

  end subroutine comms_setup

  subroutine comms_setup_vars(w90comm)
    !! Set up variables related to communicators
    !! This should be called also in library mode
    implicit none
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: ierr
    call mpi_comm_rank(w90comm%comm, my_node_id, ierr)
    call mpi_comm_size(w90comm%comm, num_nodes, ierr)
#else
    num_nodes = 1
    my_node_id = 0
#endif

    on_root = .false.
    if (my_node_id == root_id) on_root = .true.

  end subroutine comms_setup_vars

  subroutine comms_array_split(numpoints, counts, displs)
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
    use w90_io
    integer, intent(in) :: numpoints
    !! Number of elements of the array to be scattered
    integer, dimension(0:num_nodes - 1), intent(out) :: counts
    !! Array (of size num_nodes) with the number of elements of the array on each node
    integer, dimension(0:num_nodes - 1), intent(out) :: displs
    !! Array (of size num_nodes) with the displacement relative to the global array

    integer :: ratio, remainder, i

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

  subroutine comms_end
    !! Called to finalise the comms
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_finalize(ierr)
#endif

  end subroutine comms_end

  subroutine comms_barrier(w90comm)
    !! A barrier to synchronise all nodes
    implicit none
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: ierr

    call mpi_barrier(w90comm%comm, ierr)
#endif

  end subroutine comms_barrier

  subroutine comms_bcast_int(array, size, w90comm)
    !! Send integar array from root node to all nodes
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_bcast(array, size, MPI_INTEGER, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_bcast_int')
    end if
#endif
  end subroutine comms_bcast_int

  subroutine comms_bcast_real(array, size, w90comm)
    !! Send real array from root node to all nodes
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_bcast(array, size, MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_bcast_real')
    end if
#endif
  end subroutine comms_bcast_real

  subroutine comms_bcast_logical(array, size, w90comm)
    !! Send logical array from root node to all nodes
    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_bcast(array, size, MPI_LOGICAL, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_bcast_logical')
    end if
#endif
  end subroutine comms_bcast_logical

  subroutine comms_bcast_char(array, size, w90comm)
    !! Send character array from root node to all nodes
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_bcast(array, size, MPI_CHARACTER, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_bcast_char')
    end if
#endif
  end subroutine comms_bcast_char

  subroutine comms_bcast_cmplx(array, size, w90comm)
    !! Send character array from root node to all nodes

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_bcast(array, size, MPI_DOUBLE_COMPLEX, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_bcast_cmplx')
    end if
#endif
  end subroutine comms_bcast_cmplx

  !--------- SEND ----------------

  subroutine comms_send_logical(array, size, to, w90comm)
    !! Send logical array to specified node

    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_send(array, size, MPI_LOGICAL, to, mpi_send_tag, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_send_logical')
    end if
#endif
  end subroutine comms_send_logical

  subroutine comms_send_int(array, size, to, w90comm)
    !! Send integer array to specified node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_send(array, size, MPI_INTEGER, to, mpi_send_tag, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_send_int')
    end if
#endif
  end subroutine comms_send_int

  subroutine comms_send_char(array, size, to, w90comm)
    !! Send character array to specified node
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_send(array, size, MPI_CHARACTER, to, mpi_send_tag, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_send_char')
    end if
#endif
  end subroutine comms_send_char

  subroutine comms_send_real(array, size, to, w90comm)
    !! Send real array to specified node
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_send(array, size, MPI_DOUBLE_PRECISION, to, mpi_send_tag, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_send_real')
    end if
#endif
  end subroutine comms_send_real

  subroutine comms_send_cmplx(array, size, to, w90comm)
    !! Send complex array to specified node
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_send(array, size, MPI_DOUBLE_COMPLEX, to, mpi_send_tag, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_send_cmplx')
    end if
#endif
  end subroutine comms_send_cmplx

  !--------- RECV ----------------

  subroutine comms_recv_logical(array, size, from, w90comm)
    !! Receive logical array from specified node
    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
#ifdef MPIH
    integer :: status
#else
    type(mpi_status) :: status
#endif
    integer :: error

    call mpi_recv(array, size, MPI_LOGICAL, from, mpi_send_tag, w90comm%comm, status, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_recv_logical')
    end if
#endif
  end subroutine comms_recv_logical

  subroutine comms_recv_int(array, size, from, w90comm)
    !! Receive integer array from specified node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
#ifdef MPIH
    integer :: status
#else
    type(mpi_status) :: status
#endif
    integer :: error

    call mpi_recv(array, size, MPI_INTEGER, from, mpi_send_tag, w90comm%comm, status, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_recv_int')
    end if
#endif
  end subroutine comms_recv_int

  subroutine comms_recv_char(array, size, from, w90comm)
    !! Receive character array from specified node
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
#ifdef MPIH
    integer :: status
#else
    type(mpi_status) :: status
#endif
    integer :: error

    call mpi_recv(array, size, MPI_CHARACTER, from, mpi_send_tag, w90comm%comm, status, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_recv_char')
    end if
#endif
  end subroutine comms_recv_char

  subroutine comms_recv_real(array, size, from, w90comm)
    !! Receive real array from specified node
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
#ifdef MPIH
    integer :: status
#else
    type(mpi_status) :: status
#endif
    integer :: error

    call mpi_recv(array, size, MPI_DOUBLE_PRECISION, from, mpi_send_tag, w90comm%comm, status, &
                  error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_recv_real')
    end if
#endif
  end subroutine comms_recv_real

  subroutine comms_recv_cmplx(array, size, from, w90comm)
    !! Receive complex array from specified node
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
#ifdef MPIH
    integer :: status
#else
    type(mpi_status) :: status
#endif
    integer :: error

    call mpi_recv(array, size, MPI_DOUBLE_COMPLEX, from, mpi_send_tag, w90comm%comm, status, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_recv_cmplx')
    end if
#endif
  end subroutine comms_recv_cmplx

  subroutine comms_reduce_int(array, size, op, w90comm)
    !! Reduce integer data to root node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error, ierr

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
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_INTEGER, MPI_SUM, root_id, w90comm%comm, &
                        error)
      else
        call mpi_reduce(array, array, size, MPI_INTEGER, MPI_SUM, root_id, w90comm%comm, error)
      endif
    case ('PRD')
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_INTEGER, MPI_PROD, root_id, w90comm%comm, &
                        error)
      else
        call mpi_reduce(array, array, size, MPI_INTEGER, MPI_PROD, root_id, w90comm%comm, error)
      endif
    case default
      call io_error('Unknown operation in comms_reduce_int')
    end select

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_reduce_int')
    end if
#endif
  end subroutine comms_reduce_int

  subroutine comms_reduce_real(array, size, op, w90comm)
    !! Reduce real data to root node

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error, ierr

    select case (op)

    case ('SUM')
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_SUM, root_id, &
                        w90comm%comm, error)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_SUM, root_id, w90comm%comm, &
                        error)
      endif
    case ('PRD')
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_PROD, root_id, &
                        w90comm%comm, error)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_PROD, root_id, w90comm%comm, &
                        error)
      endif
    case ('MIN')
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MIN, root_id, &
                        w90comm%comm, error)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_MIN, root_id, w90comm%comm, &
                        error)
      endif
    case ('MAX')
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MAX, root_id, &
                        w90comm%comm, error)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_PRECISION, MPI_MAX, root_id, w90comm%comm, &
                        error)
      endif
    case default
      call io_error('Unknown operation in comms_reduce_real')

    end select

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_reduce_real')
    end if
#endif
  end subroutine comms_reduce_real

  subroutine comms_reduce_cmplx(array, size, op, w90comm)
    !! Reduce complex data to root node

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error, ierr

    select case (op)

    case ('SUM')
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_SUM, root_id, &
                        w90comm%comm, error)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_COMPLEX, MPI_SUM, root_id, w90comm%comm, &
                        error)
      end if
    case ('PRD')
      if (on_root) then
        call mpi_reduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_PROD, root_id, &
                        w90comm%comm, error)
      else
        call mpi_reduce(array, array, size, MPI_DOUBLE_COMPLEX, MPI_PROD, root_id, w90comm%comm, &
                        error)
      end if
    case default
      call io_error('Unknown operation in comms_reduce_cmplx')

    end select

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_reduce_cmplx')
    end if

#endif
  end subroutine comms_reduce_cmplx

  subroutine comms_allreduce_real(array, size, op, w90comm)
    !! Reduce real data to all nodes

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error, ierr

    select case (op)

    case ('SUM')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_SUM, w90comm%comm, &
                         error)
    case ('PRD')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_PROD, w90comm%comm, &
                         error)
    case ('MIN')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MIN, w90comm%comm, &
                         error)
    case ('MAX')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_PRECISION, MPI_MAX, w90comm%comm, &
                         error)
    case default
      call io_error('Unknown operation in comms_allreduce_real')

    end select

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_allreduce_real')
    end if
#endif
  end subroutine comms_allreduce_real

  subroutine comms_allreduce_cmplx(array, size, op, w90comm)
    !! Reduce complex data to all nodes
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error, ierr

    select case (op)

    case ('SUM')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_SUM, w90comm%comm, &
                         error)
    case ('PRD')
      call mpi_allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE_COMPLEX, MPI_PROD, w90comm%comm, &
                         error)
    case default
      call io_error('Unknown operation in comms_allreduce_cmplx')

    end select

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_allreduce_cmplx')
    end if
#endif
  end subroutine comms_allreduce_cmplx

  subroutine comms_gatherv_real_1(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather real data to root node (for arrays of rank 1)
    implicit none

    real(kind=dp), dimension(:), intent(inout)   :: array
    !! local array for sending data
    integer, intent(in)                          :: localcount
    !! localcount elements will be sent to the root node
    real(kind=dp), dimension(:), intent(inout)   :: rootglobalarray
    !! array on the root node to which data will be sent
    integer, dimension(num_nodes), intent(in)    :: counts
    !! how data should be partitioned, see MPI documentation or
    !! function comms_array_split
    integer, dimension(num_nodes), intent(in)    :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, &
                     displs, MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_real_1')
    end if
#else
    !call dcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif
  end subroutine comms_gatherv_real_1

  subroutine comms_gatherv_real_2(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather real data to root node (for arrays of rank 2)
    implicit none

    real(kind=dp), dimension(:, :), intent(inout) :: array
    !! local array for sending data
    integer, intent(in)                          :: localcount
    !! localcount elements will be sent to the root node
    real(kind=dp), dimension(:, :), intent(inout) :: rootglobalarray
    !! array on the root node to which data will be sent
    integer, dimension(num_nodes), intent(in)    :: counts
    !! how data should be partitioned, see MPI documentation or
    !! function comms_array_split
    integer, dimension(num_nodes), intent(in)    :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, &
                     displs, MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_real_2')
    end if
#else
    !call dcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif
  end subroutine comms_gatherv_real_2

  subroutine comms_gatherv_real_3(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather real data to root node (for arrays of rank 3)
    implicit none

    real(kind=dp), dimension(:, :, :), intent(inout) :: array
    !! local array for sending data
    integer, intent(in)                            :: localcount
    !! localcount elements will be sent to the root node
    real(kind=dp), dimension(:, :, :), intent(inout) :: rootglobalarray
    !! array on the root node to which data will be sent
    integer, dimension(num_nodes), intent(in)      :: counts
    !! how data should be partitioned, see MPI documentation or
    !! function comms_array_split
    integer, dimension(num_nodes), intent(in)      :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, &
                     displs, MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_real_3')
    end if

#else
    !call dcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif
  end subroutine comms_gatherv_real_3

  subroutine comms_gatherv_real_2_3(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather real data to root node (for arrays of rank 2 and 3, respectively)
    implicit none

    real(kind=dp), dimension(:, :), intent(inout) :: array
    !! local array for sending data
    integer, intent(in)                            :: localcount
    !! localcount elements will be sent to the root node
    real(kind=dp), dimension(:, :, :), intent(inout) :: rootglobalarray
    !! array on the root node to which data will be sent
    integer, dimension(num_nodes), intent(in)      :: counts
    !! how data should be partitioned, see MPI documentation or
    !! function comms_array_split
    integer, dimension(num_nodes), intent(in)      :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_PRECISION, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_real_2_3')
    end if

#else
    call dcopy(localcount, array, 1, rootglobalarray, 1)
    !rootglobalarray = array
    !JJ when/why should the shapes not match?
#endif
  end subroutine comms_gatherv_real_2_3

  ! Array: local array for sending data; localcount elements will be sent
  !        to the root node
  ! rootglobalarray: array on the root node to which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split

  subroutine comms_gatherv_cmplx_1(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather complex data to root node (for arrays of rank 1)
    implicit none

    complex(kind=dp), dimension(:), intent(inout)     :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:), intent(inout)     :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_cmplx_1')
    end if

#else
    !call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif
  end subroutine comms_gatherv_cmplx_1

  subroutine comms_gatherv_cmplx_2(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather complex data to root node (for arrays of rank 2)
    implicit none

    complex(kind=dp), dimension(:, :), intent(inout)   :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:, :), intent(inout)   :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_cmplx_2')
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif
  end subroutine comms_gatherv_cmplx_2

  subroutine comms_gatherv_cmplx_3(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather complex data to root node (for arrays of rank 3)
    implicit none

    complex(kind=dp), dimension(:, :, :), intent(inout) :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:, :, :), intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_cmplx_3')
    end if

#else
    !call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif
  end subroutine comms_gatherv_cmplx_3

  subroutine comms_gatherv_cmplx_3_4(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather complex data to root node (for arrays of rank 3 and 4, respectively)
    implicit none

    complex(kind=dp), dimension(:, :, :), intent(inout)   :: array
    integer, intent(in)                                 :: localcount
    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in)           :: counts
    integer, dimension(num_nodes), intent(in)           :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_cmplx_3_4')
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
    !rootglobalarray = array
    !JJ when/why should the shapes not match?
#endif
  end subroutine comms_gatherv_cmplx_3_4

  subroutine comms_gatherv_cmplx_4(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather complex data to root node (for arrays of rank 4)
    implicit none

    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_DOUBLE_COMPLEX, rootglobalarray, counts, displs, &
                     MPI_DOUBLE_COMPLEX, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_cmplx_4')
    end if

#else
    !call zcopy(localcount, array, 1, rootglobalarray, 1)
    rootglobalarray = array
#endif
  end subroutine comms_gatherv_cmplx_4

  subroutine comms_gatherv_logical(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Gather real data to root node
    implicit none

    logical, intent(inout)           :: array
    !! local array for sending data
    integer, intent(in)                       :: localcount
    !! localcount elements will be sent to the root node
    logical, intent(inout)           :: rootglobalarray
    !! array on the root node to which data will be sent
    integer, dimension(num_nodes), intent(in) :: counts
    !! how data should be partitioned, see MPI documentation or
    !! function comms_array_split
    integer, dimension(num_nodes), intent(in) :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_gatherv(array, localcount, MPI_LOGICAL, rootglobalarray, counts, displs, &
                     MPI_LOGICAL, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_gatherv_logical')
    end if
#else
!    rootglobalarray(1:localcount)=array(1:localcount)
#endif

  end subroutine comms_gatherv_logical

  subroutine comms_scatterv_real_1(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Scatter real data from root node (array of rank 1)
    implicit none

    real(kind=dp), dimension(:), intent(inout) :: array
    !! local array for getting data
    integer, intent(in)                        :: localcount
    !! localcount elements will be fetched from the root node
    real(kind=dp), dimension(:), intent(inout) :: rootglobalarray
    !! array on the root node from which data will be sent
    integer, dimension(num_nodes), intent(in)  :: counts
    !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, dimension(num_nodes), intent(in)  :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_PRECISION, array, localcount, &
                      MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_scatterv_real_1')
    end if

#else
    !call dcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif
  end subroutine comms_scatterv_real_1

  subroutine comms_scatterv_real_2(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Scatter real data from root node (array of rank 2)
    implicit none

    real(kind=dp), dimension(:, :), intent(inout) :: array
    !! local array for getting data
    integer, intent(in)                          :: localcount
    !! localcount elements will be fetched from the root node
    real(kind=dp), dimension(:, :), intent(inout) :: rootglobalarray
    !! array on the root node from which data will be sent
    integer, dimension(num_nodes), intent(in)    :: counts
    !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, dimension(num_nodes), intent(in)    :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_PRECISION, array, localcount, &
                      MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_scatterv_real_2')
    end if

#else
    !call dcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif
  end subroutine comms_scatterv_real_2

  subroutine comms_scatterv_real_3(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Scatter real data from root node (array of rank 3)
    implicit none

    real(kind=dp), dimension(:, :, :), intent(inout) :: array
    !! local array for getting data
    integer, intent(in)                            :: localcount
    !! localcount elements will be fetched from the root node
    real(kind=dp), dimension(:, :, :), intent(inout) :: rootglobalarray
    !! array on the root node from which data will be sent
    integer, dimension(num_nodes), intent(in)      :: counts
    !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, dimension(num_nodes), intent(in)      :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_PRECISION, array, localcount, &
                      MPI_DOUBLE_PRECISION, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_scatterv_real_3')
    end if

#else
    !call dcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif
  end subroutine comms_scatterv_real_3

  subroutine comms_scatterv_cmplx_4(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Scatter complex data from root node (array of rank 4)
    implicit none

    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: array
    !! local array for getting data
    integer, intent(in)                            :: localcount
    !! localcount elements will be fetched from the root node
    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: rootglobalarray
    !! array on the root node from which data will be sent
    integer, dimension(num_nodes), intent(in)      :: counts
    !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, dimension(num_nodes), intent(in)      :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_DOUBLE_COMPLEX, array, localcount, &
                      MPI_DOUBLE_COMPLEX, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_scatterv_cmplx_4')
    end if

#else
    !call zcopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif
  end subroutine comms_scatterv_cmplx_4

  subroutine comms_scatterv_int_1(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Scatter integer data from root node (array of rank 1)
    implicit none

    integer, dimension(:), intent(inout)      :: array
    !! local array for getting data
    integer, intent(in)                       :: localcount
    !! localcount elements will be fetched from the root node
    integer, dimension(:), intent(inout)      :: rootglobalarray
    !!  array on the root node from which data will be sent
    integer, dimension(num_nodes), intent(in) :: counts
    !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, dimension(num_nodes), intent(in) :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_INTEGER, array, localcount, &
                      MPI_INTEGER, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_scatterv_real')
    end if

#else
    !call my_icopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif
  end subroutine comms_scatterv_int_1

  subroutine comms_scatterv_int_2(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Scatter integer data from root node (array of rank 2)

    implicit none

    integer, dimension(:, :), intent(inout)    :: array
    !! local array for getting data
    integer, intent(in)                       :: localcount
    !! localcount elements will be fetched from the root node
    integer, dimension(:, :), intent(inout)    :: rootglobalarray
    !!  array on the root node from which data will be sent
    integer, dimension(num_nodes), intent(in) :: counts
    !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, dimension(num_nodes), intent(in) :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_INTEGER, array, localcount, &
                      MPI_INTEGER, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_scatterv_int_2')
    end if

#else
    !call my_icopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif
  end subroutine comms_scatterv_int_2

  subroutine comms_scatterv_int_3(array, localcount, rootglobalarray, counts, displs, w90comm)
    !! Scatter integer data from root node (array of rank 3)

    implicit none

    integer, dimension(:, :, :), intent(inout)  :: array
    !! local array for getting data
    integer, intent(in)                       :: localcount
    !! localcount elements will be fetched from the root node
    integer, dimension(:, :, :), intent(inout)  :: rootglobalarray
    !!  array on the root node from which data will be sent
    integer, dimension(num_nodes), intent(in) :: counts
    !! how data should be partitioned, see MPI documentation or function comms_array_split
    integer, dimension(num_nodes), intent(in) :: displs
    type(w90commtype), intent(in) :: w90comm

#ifdef MPI
    integer :: error

    call mpi_scatterv(rootglobalarray, counts, displs, MPI_INTEGER, array, localcount, &
                      MPI_INTEGER, root_id, w90comm%comm, error)

    if (error .ne. MPI_SUCCESS) then
      call io_error('Error in comms_scatterv_int_3')
    end if

#else
    !call my_icopy(localcount, rootglobalarray, 1, array, 1)
    array = rootglobalarray
#endif
  end subroutine comms_scatterv_int_3

end module w90_comms
