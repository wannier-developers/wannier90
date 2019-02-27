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

module w90_comms
  !! This module handles all of the communications

  use w90_constants, only: dp
  use w90_io, only: io_error

  implicit none

  private

#ifdef MPI
  include 'mpif.h'
#endif

  logical, public, save :: on_root
  !! Are we the root node
  integer, public, save :: num_nodes
  !! Number of nodes
  integer, public, save :: my_node_id
  !! ID of this node
  integer, public, parameter :: root_id = 0
  !! ID of the root node

  integer, parameter :: mpi_send_tag = 77 !abitrary

  public :: comms_setup
  public :: comms_setup_vars
  public :: comms_end
!  public :: comms_abort     ! [GP]: do not use, use io_error instead
  public :: comms_bcast      ! send data from the root node
  public :: comms_send       ! send data from one node to another
  public :: comms_recv       ! accept data from one node to another
  public :: comms_reduce     ! reduce data onto root node (n.b. not allreduce);
  ! note that on all other nodes, the data is lost
  public :: comms_allreduce  ! reduce data onto all nodes
  public :: comms_barrier    ! puts a barrier so that the code goes on only when all nodes reach the barrier
  public :: comms_gatherv    ! gets chunks of an array from all nodes and gathers them on the root node
  public :: comms_scatterv    ! sends chunks of an array to all nodes scattering them from the root node

  public :: comms_array_split

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

  subroutine comms_setup
    !! Set up communications
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_init(ierr)
    if (ierr .ne. 0) call io_error('MPI initialisation error')
#endif

    call comms_setup_vars

  end subroutine comms_setup

  subroutine comms_setup_vars
    !! Set up variables related to communicators
    !! This should be called also in library mode
    implicit none

#ifdef MPI
    integer :: ierr
    call mpi_comm_rank(mpi_comm_world, my_node_id, ierr)
    call mpi_comm_size(mpi_comm_world, num_nodes, ierr)
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

  subroutine comms_barrier
    !! A barrier to synchronise all nodes
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_barrier(mpi_comm_world, ierr)
#endif

  end subroutine comms_barrier

!  subroutine comms_abort
!
!    implicit none
!
!    integer :: ierr
!
!#ifdef MPI
!    call MPI_abort(MPI_comm_world,1,ierr)
!#else
!    STOP
!#endif
!
!  end subroutine comms_abort

  subroutine comms_bcast_int(array, size)
    !! Send integar array from root node to all nodes
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array, size, MPI_integer, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_bcast_int')
    end if
#endif

    return

  end subroutine comms_bcast_int

  subroutine comms_bcast_real(array, size)
    !! Send real array from root node to all nodes
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array, size, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_bcast_real')
    end if
#endif

    return

  end subroutine comms_bcast_real

  subroutine comms_bcast_logical(array, size)
    !! Send logical array from root node to all nodes
    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array, size, MPI_logical, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_bcast_logical')
    end if
#endif

    return

  end subroutine comms_bcast_logical

  subroutine comms_bcast_char(array, size)
    !! Send character array from root node to all nodes
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array, size, MPI_character, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_bcast_char')
    end if
#endif

    return

  end subroutine comms_bcast_char

  subroutine comms_bcast_cmplx(array, size)
    !! Send character array from root node to all nodes

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array, size, MPI_double_complex, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_bcast_cmplx')
    end if
#endif

    return

  end subroutine comms_bcast_cmplx

  !--------- SEND ----------------

  subroutine comms_send_logical(array, size, to)
    !! Send logical array to specified node

    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array, size, MPI_logical, to, &
                  mpi_send_tag, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_send_logical')
    end if
#endif

    return

  end subroutine comms_send_logical

  subroutine comms_send_int(array, size, to)
    !! Send integer array to specified node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array, size, MPI_integer, to, &
                  mpi_send_tag, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_send_int')
    end if
#endif

    return

  end subroutine comms_send_int

  subroutine comms_send_char(array, size, to)
    !! Send character array to specified node
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array, size, MPI_character, to, &
                  mpi_send_tag, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_send_char')
    end if
#endif

    return

  end subroutine comms_send_char

  subroutine comms_send_real(array, size, to)
    !! Send real array to specified node
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array, size, MPI_double_precision, to, &
                  mpi_send_tag, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_send_real')
    end if
#endif

    return

  end subroutine comms_send_real

  subroutine comms_send_cmplx(array, size, to)
    !! Send complex array to specified node
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array, size, MPI_double_complex, to, &
                  mpi_send_tag, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_send_cmplx')
    end if
#endif

    return

  end subroutine comms_send_cmplx

  !--------- RECV ----------------

  subroutine comms_recv_logical(array, size, from)
    !! Receive logical array from specified node
    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array, size, MPI_logical, from, &
                  mpi_send_tag, mpi_comm_world, status, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_recv_logical')
    end if
#endif

    return

  end subroutine comms_recv_logical

  subroutine comms_recv_int(array, size, from)
    !! Receive integer array from specified node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array, size, MPI_integer, from, &
                  mpi_send_tag, mpi_comm_world, status, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_recv_int')
    end if
#endif

    return

  end subroutine comms_recv_int

  subroutine comms_recv_char(array, size, from)
    !! Receive character array from specified node
    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array, size, MPI_character, from, &
                  mpi_send_tag, mpi_comm_world, status, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_recv_char')
    end if
#endif

    return

  end subroutine comms_recv_char

  subroutine comms_recv_real(array, size, from)
    !! Receive real array from specified node
    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array, size, MPI_double_precision, from, &
                  mpi_send_tag, mpi_comm_world, status, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_recv_real')
    end if
#endif

    return

  end subroutine comms_recv_real

  subroutine comms_recv_cmplx(array, size, from)
    !! Receive complex array from specified node
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error

    integer :: status(MPI_status_size)

    call MPI_recv(array, size, MPI_double_complex, from, &
                  mpi_send_tag, mpi_comm_world, status, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_recv_cmplx')
    end if

#endif

    return

  end subroutine comms_recv_cmplx

!  subroutine comms_error
!
!    implicit none
!
!#ifdef MPI
!    integer :: error
!
!    call MPI_abort(MPI_comm_world,1,error)
!
!#endif
!
!  end subroutine comms_error

  ! COMMS_REDUCE (collect data on the root node)

  subroutine comms_reduce_int(array, size, op)
    !! Reduce integer data to root node
    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error, ierr

    integer, allocatable :: array_red(:)

    allocate (array_red(size), stat=ierr)
    if (ierr /= 0) then
      call io_error('failure to allocate array_red in comms_reduce_int')
    end if

    select case (op)

    case ('SUM')
      call MPI_reduce(array, array_red, size, MPI_integer, MPI_sum, root_id, mpi_comm_world, error)
    case ('PRD')
      call MPI_reduce(array, array_red, size, MPI_integer, MPI_prod, root_id, mpi_comm_world, error)
    case default
      call io_error('Unknown operation in comms_reduce_int')

    end select

    call my_icopy(size, array_red, 1, array, 1)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_reduce_int')
    end if

    if (allocated(array_red)) deallocate (array_red)
#endif

    return

  end subroutine comms_reduce_int

  subroutine comms_reduce_real(array, size, op)
    !! Reduce real data to root node

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error, ierr

    real(kind=dp), allocatable :: array_red(:)

    allocate (array_red(size), stat=ierr)
    if (ierr /= 0) then
      call io_error('failure to allocate array_red in comms_reduce_real')
    end if

    select case (op)

    case ('SUM')
      call MPI_reduce(array, array_red, size, MPI_double_precision, MPI_sum, root_id, mpi_comm_world, error)
    case ('PRD')
      call MPI_reduce(array, array_red, size, MPI_double_precision, MPI_prod, root_id, mpi_comm_world, error)
    case ('MIN')
      call MPI_reduce(array, array_red, size, MPI_double_precision, MPI_MIN, root_id, mpi_comm_world, error)
    case ('MAX')
      call MPI_reduce(array, array_red, size, MPI_double_precision, MPI_max, root_id, mpi_comm_world, error)
    case default
      call io_error('Unknown operation in comms_reduce_real')

    end select

    call dcopy(size, array_red, 1, array, 1)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_reduce_real')
    end if

    if (allocated(array_red)) deallocate (array_red)
#endif

    return

  end subroutine comms_reduce_real

  subroutine comms_reduce_cmplx(array, size, op)
    !! Reduce complex data to root node

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error, ierr

    complex(kind=dp), allocatable :: array_red(:)

    allocate (array_red(size), stat=ierr)
    if (ierr /= 0) then
      call io_error('failure to allocate array_red in comms_reduce_cmplx')
    end if

    select case (op)

    case ('SUM')
      call MPI_reduce(array, array_red, size, MPI_double_complex, MPI_sum, root_id, mpi_comm_world, error)
    case ('PRD')
      call MPI_reduce(array, array_red, size, MPI_double_complex, MPI_prod, root_id, mpi_comm_world, error)
    case default
      call io_error('Unknown operation in comms_reduce_cmplx')

    end select

    call zcopy(size, array_red, 1, array, 1)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_reduce_cmplx')
    end if

    if (allocated(array_red)) deallocate (array_red)
#endif

    return

  end subroutine comms_reduce_cmplx

  subroutine comms_allreduce_real(array, size, op)
    !! Reduce real data to all nodes

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error, ierr

    real(kind=dp), allocatable :: array_red(:)

    allocate (array_red(size), stat=ierr)
    if (ierr /= 0) then
      call io_error('failure to allocate array_red in comms_allreduce_real')
    end if

    select case (op)

    case ('SUM')
      call MPI_allreduce(array, array_red, size, MPI_double_precision, MPI_sum, mpi_comm_world, error)
    case ('PRD')
      call MPI_allreduce(array, array_red, size, MPI_double_precision, MPI_prod, mpi_comm_world, error)
    case ('MIN')
      call MPI_allreduce(array, array_red, size, MPI_double_precision, MPI_MIN, mpi_comm_world, error)
    case ('MAX')
      call MPI_allreduce(array, array_red, size, MPI_double_precision, MPI_max, mpi_comm_world, error)
    case default
      call io_error('Unknown operation in comms_allreduce_real')

    end select

    call dcopy(size, array_red, 1, array, 1)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_allreduce_real')
    end if

    if (allocated(array_red)) deallocate (array_red)
#endif

    return

  end subroutine comms_allreduce_real

  subroutine comms_allreduce_cmplx(array, size, op)
    !! Reduce complex data to all nodes
    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error, ierr

    complex(kind=dp), allocatable :: array_red(:)

    allocate (array_red(size), stat=ierr)
    if (ierr /= 0) then
      call io_error('failure to allocate array_red in comms_allreduce_cmplx')
    end if

    select case (op)

    case ('SUM')
      call MPI_allreduce(array, array_red, size, MPI_double_complex, MPI_sum, mpi_comm_world, error)
    case ('PRD')
      call MPI_allreduce(array, array_red, size, MPI_double_complex, MPI_prod, mpi_comm_world, error)
    case default
      call io_error('Unknown operation in comms_allreduce_cmplx')

    end select

    call zcopy(size, array_red, 1, array, 1)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_allreduce_cmplx')
    end if

    if (allocated(array_red)) deallocate (array_red)
#endif

    return

  end subroutine comms_allreduce_cmplx

  subroutine comms_gatherv_real_1(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_precision, rootglobalarray, counts, &
                     displs, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_real_1')
    end if

#else
    call dcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_real_1

  subroutine comms_gatherv_real_2(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_precision, rootglobalarray, counts, &
                     displs, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_real_2')
    end if

#else
    call dcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_real_2

  subroutine comms_gatherv_real_3(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_precision, rootglobalarray, counts, &
                     displs, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_real_3')
    end if

#else
    call dcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_real_3

  subroutine comms_gatherv_real_2_3(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_precision, rootglobalarray, counts, &
                     displs, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_real_2_3')
    end if

#else
    call dcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_real_2_3

  ! Array: local array for sending data; localcount elements will be sent
  !        to the root node
  ! rootglobalarray: array on the root node to which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split

  subroutine comms_gatherv_cmplx_1(array, localcount, rootglobalarray, counts, displs)
    !! Gather complex data to root node (for arrays of rank 1)
    implicit none

    complex(kind=dp), dimension(:), intent(inout)     :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:), intent(inout)     :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_complex, rootglobalarray, counts, &
                     displs, MPI_double_complex, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_cmplx_1')
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_cmplx_1

  subroutine comms_gatherv_cmplx_2(array, localcount, rootglobalarray, counts, displs)
    !! Gather complex data to root node (for arrays of rank 2)
    implicit none

    complex(kind=dp), dimension(:, :), intent(inout)   :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:, :), intent(inout)   :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_complex, rootglobalarray, counts, &
                     displs, MPI_double_complex, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_cmplx_2')
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_cmplx_2

!!JRY  subroutine comms_gatherv_logical(array,localcount,rootglobalarray,counts,displs)
!!    !! Gather real data to root node
!!    implicit none

  subroutine comms_gatherv_cmplx_3(array, localcount, rootglobalarray, counts, displs)
    !! Gather complex data to root node (for arrays of rank 3)
    implicit none

    complex(kind=dp), dimension(:, :, :), intent(inout) :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:, :, :), intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_complex, rootglobalarray, counts, &
                     displs, MPI_double_complex, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_cmplx_3')
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_cmplx_3

  subroutine comms_gatherv_cmplx_3_4(array, localcount, rootglobalarray, counts, displs)
    !! Gather complex data to root node (for arrays of rank 3 and 4, respectively)
    implicit none

    complex(kind=dp), dimension(:, :, :), intent(inout)   :: array
    integer, intent(in)                                 :: localcount
    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in)           :: counts
    integer, dimension(num_nodes), intent(in)           :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_complex, rootglobalarray, counts, &
                     displs, MPI_double_complex, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_cmplx_3_4')
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_cmplx_3_4

  subroutine comms_gatherv_cmplx_4(array, localcount, rootglobalarray, counts, displs)
    !! Gather complex data to root node (for arrays of rank 4)
    implicit none

    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: array
    integer, intent(in)                               :: localcount
    complex(kind=dp), dimension(:, :, :, :), intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in)         :: counts
    integer, dimension(num_nodes), intent(in)         :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_double_complex, rootglobalarray, counts, &
                     displs, MPI_double_complex, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_cmplx_4')
    end if

#else
    call zcopy(localcount, array, 1, rootglobalarray, 1)
#endif

    return

  end subroutine comms_gatherv_cmplx_4

  subroutine comms_gatherv_logical(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_gatherv(array, localcount, MPI_logical, rootglobalarray, counts, &
                     displs, MPI_logical, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_gatherv_logical')
    end if
#else
!    rootglobalarray(1:localcount)=array(1:localcount)
#endif

  end subroutine comms_gatherv_logical

  subroutine comms_scatterv_real_1(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray, counts, displs, MPI_double_precision, &
                      array, localcount, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_scatterv_real_1')
    end if

#else
    call dcopy(localcount, rootglobalarray, 1, array, 1)
#endif

    return

  end subroutine comms_scatterv_real_1

  subroutine comms_scatterv_real_2(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray, counts, displs, MPI_double_precision, &
                      array, localcount, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_scatterv_real_2')
    end if

#else
    call dcopy(localcount, rootglobalarray, 1, array, 1)
#endif

    return

  end subroutine comms_scatterv_real_2

  subroutine comms_scatterv_real_3(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray, counts, displs, MPI_double_precision, &
                      array, localcount, MPI_double_precision, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_scatterv_real_3')
    end if

#else
    call dcopy(localcount, rootglobalarray, 1, array, 1)
#endif

    return

  end subroutine comms_scatterv_real_3

  subroutine comms_scatterv_cmplx_4(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray, counts, displs, MPI_double_complex, &
                      array, localcount, MPI_double_complex, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_scatterv_cmplx_4')
    end if

#else
    call zcopy(localcount, rootglobalarray, 1, array, 1)
#endif

    return

  end subroutine comms_scatterv_cmplx_4

  subroutine comms_scatterv_int_1(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray, counts, displs, MPI_Integer, &
                      Array, localcount, MPI_Integer, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_scatterv_real')
    end if

#else
    call my_icopy(localcount, rootglobalarray, 1, array, 1)
#endif

    return

  end subroutine comms_scatterv_int_1

  subroutine comms_scatterv_int_2(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray, counts, displs, MPI_Integer, &
                      Array, localcount, MPI_Integer, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_scatterv_int_2')
    end if

#else
    call my_icopy(localcount, rootglobalarray, 1, array, 1)
#endif

    return

  end subroutine comms_scatterv_int_2

  subroutine comms_scatterv_int_3(array, localcount, rootglobalarray, counts, displs)
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

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray, counts, displs, MPI_Integer, &
                      Array, localcount, MPI_Integer, root_id, mpi_comm_world, error)

    if (error .ne. MPI_success) then
      call io_error('Error in comms_scatterv_int_3')
    end if

#else
    call my_icopy(localcount, rootglobalarray, 1, array, 1)
#endif

    return

  end subroutine comms_scatterv_int_3

end module w90_comms

subroutine my_ICOPY(N, ZX, INCX, ZY, INCY)
  !     .. Scalar Arguments ..
  integer INCX, INCY, N
  !     ..
  !     .. Array Arguments ..
  integer ZX(*), ZY(*)
  !     ..
  !
  !  Purpose
  !  =======
  !
  !     copies a vector, x, to a vector, y.
  !     jack dongarra, linpack, 4/11/78.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  !
  !     .. Local Scalars ..
  integer I, IX, IY
  !     ..
  if (N .le. 0) return
  if (INCX .eq. 1 .and. INCY .eq. 1) GO TO 20
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
  IX = 1
  IY = 1
  if (INCX .lt. 0) IX = (-N + 1)*INCX + 1
  if (INCY .lt. 0) IY = (-N + 1)*INCY + 1
  do I = 1, N
    ZY(IY) = ZX(IX)
    IX = IX + INCX
    IY = IY + INCY
  end do
  return
  !
  !        code for both increments equal to 1
  !
20 do I = 1, N
    ZY(I) = ZX(I)
  end do
  return
end subroutine my_ICOPY
