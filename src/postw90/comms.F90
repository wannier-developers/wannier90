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
!------------------------------------------------------------!                                                            !
!  COMMS: set of MPI wrapper                                 !
!  written 2006-2012 Jonathan R. Yates                       !
!    later additions Giovanni Pizzi                          !
!------------------------------------------------------------!


module w90_comms

  use w90_constants, only : dp
  use w90_io, only: io_error

  implicit none


  private

#ifdef MPI
  include 'mpif.h'
#endif

  logical, public, save :: on_root
  integer, public, save :: num_nodes,my_node_id
  integer, public, parameter :: root_id=0

  integer, parameter :: mpi_send_tag=77 !abitrary

  public :: comms_setup
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
     module procedure comms_gatherv_real
!     module procedure comms_gatherv_cmplx
  end interface comms_gatherv

  interface comms_scatterv
     module procedure comms_scatterv_int  
     module procedure comms_scatterv_real
!     module procedure comms_scatterv_cmplx
  end interface comms_scatterv

contains

  subroutine comms_setup
 
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_init(ierr)
    if (ierr.ne.0) call io_error('MPI initialisation error')
    call mpi_comm_rank(mpi_comm_world, my_node_id, ierr)
    call mpi_comm_size(mpi_comm_world, num_nodes, ierr)
#else
    num_nodes=1
    my_node_id=0
#endif

    on_root=.false.
    if(my_node_id==root_id) on_root=.true.
    
  end subroutine comms_setup

  !> Given an array of size numpoints, we want to split on num_nodes nodes. This function returns
  !> two arrays: count and displs.
  !> The i-th element of the count array gives the number of elements
  !> that must be calculated by the process with id (i-1).
  !> The i-th element of the displs array gives the displacement of the array calculated locally on
  !> the process with id (i-1) with respect to the global array.
  !>
  !> \note These values are those to be passed to the functions MPI_Scatterv, MPI_Gatherv and MPI_Alltoallv.
  !>
  !> \note one can use the following do loop to run over the needed elements, if the full array is stored
  !> on all nodes:
  !> do i=displs(my_node_id)+1,displs(my_node_id)+counts(my_node_id)
  !> 
  !> \param numpoints Number of elements of the array to be scattered
  !> \param counts    Array (of size num_nodes) with the number of elements of the array on each node
  !> \param displs    Array (of size num_nodes) with the displacement relative to the global array
  subroutine comms_array_split(numpoints,counts,displs)
    integer, intent(in) :: numpoints
    integer, dimension(0:num_nodes-1), intent(out) :: counts
    integer, dimension(0:num_nodes-1), intent(out) :: displs

    integer :: ratio, remainder, i

    ratio = numpoints / num_nodes
    remainder = MOD(numpoints, num_nodes)

    do i=0,num_nodes-1
       if (i < remainder) then
          counts(i) = ratio+1
          displs(i) = i*(ratio+1)
       else
          counts(i) = ratio
          displs(i) = remainder*(ratio+1) + (i-remainder)*ratio
       end if
    end do

  end subroutine comms_array_split

  subroutine comms_end
 
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_finalize(ierr)
#endif
    
  end subroutine comms_end

  subroutine comms_barrier
 
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


  subroutine comms_bcast_int(array,size)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_integer,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_bcast_int')
    end if
#endif

    return

  end subroutine comms_bcast_int

  subroutine comms_bcast_real(array,size)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_double_precision,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_bcast_real')
    end if
#endif

    return

  end subroutine comms_bcast_real

  subroutine comms_bcast_logical(array,size)

    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_logical,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_bcast_logical')
    end if
#endif

    return

  end subroutine comms_bcast_logical

  subroutine comms_bcast_char(array,size)

    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size


#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_character,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_bcast_char')
    end if
#endif

    return

  end subroutine comms_bcast_char

  subroutine comms_bcast_cmplx(array,size)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size


#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_double_complex,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_bcast_cmplx')
    end if
#endif

    return

  end subroutine comms_bcast_cmplx


  !--------- SEND ----------------

  subroutine comms_send_logical(array,size,to)

    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array,size,MPI_logical,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_send_logical')
    end if
#endif

    return

  end subroutine comms_send_logical


  subroutine comms_send_int(array,size,to)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array,size,MPI_integer,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_send_int')
    end if
#endif

    return

  end subroutine comms_send_int


  subroutine comms_send_char(array,size,to)

    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array,size,MPI_character,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_send_char')
    end if
#endif

    return

  end subroutine comms_send_char


  subroutine comms_send_real(array,size,to)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array,size,MPI_double_precision,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_send_real')
    end if
#endif

    return

  end subroutine comms_send_real


  subroutine comms_send_cmplx(array,size,to)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to


#ifdef MPI
    integer :: error

    call MPI_send(array,size,MPI_double_complex,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_send_cmplx')
    end if
#endif

    return

  end subroutine comms_send_cmplx


  !--------- RECV ----------------

  subroutine comms_recv_logical(array,size,from)

    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_logical,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_recv_logical')
    end if
#endif

    return

  end subroutine comms_recv_logical


  subroutine comms_recv_int(array,size,from)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_integer,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_recv_int')
    end if
#endif

    return

  end subroutine comms_recv_int


  subroutine comms_recv_char(array,size,from)

    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_character,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_recv_char')
    end if
#endif

    return

  end subroutine comms_recv_char


  subroutine comms_recv_real(array,size,from)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_double_precision,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_recv_real')
    end if
#endif

    return

  end subroutine comms_recv_real


  subroutine comms_recv_cmplx(array,size,from)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

#ifdef MPI
    integer :: error

    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_double_complex,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
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

  subroutine comms_reduce_int(array,size,op)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error,ierr

    integer, allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       call io_error('failure to allocate array_red in comms_reduce_int')
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_integer,MPI_sum,root_id,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_integer,MPI_prod,root_id,mpi_comm_world,error)
    case default
       call io_error('Unknown operation in comms_reduce_int')

    end select

    call my_icopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_reduce_int')
    end if

    if (allocated(array_red)) deallocate(array_red)
#endif

    return

  end subroutine comms_reduce_int


  subroutine comms_reduce_real(array,size,op)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error,ierr

    real(kind=dp), allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       call io_error('failure to allocate array_red in comms_reduce_real')
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_sum,root_id,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_prod,root_id,mpi_comm_world,error)
    case ('MIN')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_MIN,root_id,mpi_comm_world,error)
    case ('MAX')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_max,root_id,mpi_comm_world,error)
    case default
       call io_error('Unknown operation in comms_reduce_real')

    end select

    call dcopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_reduce_real')
    end if

    if (allocated(array_red)) deallocate(array_red)
#endif

    return

  end subroutine comms_reduce_real


  subroutine comms_reduce_cmplx(array,size,op)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error,ierr

    complex(kind=dp), allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       call io_error('failure to allocate array_red in comms_reduce_cmplx')
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_double_complex,MPI_sum,root_id,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_double_complex,MPI_prod,root_id,mpi_comm_world,error)
    case default
       call io_error('Unknown operation in comms_reduce_cmplx')

    end select

    call zcopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_reduce_cmplx')
    end if

    if (allocated(array_red)) deallocate(array_red)
#endif

    return

  end subroutine comms_reduce_cmplx

  subroutine comms_allreduce_real(array,size,op)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error,ierr

    real(kind=dp), allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       call io_error('failure to allocate array_red in comms_allreduce_real')
    end if

    select case(op)

    case ('SUM')
       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_sum,mpi_comm_world,error)
    case ('PRD')
       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_prod,mpi_comm_world,error)
    case ('MIN')
       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_MIN,mpi_comm_world,error)
    case ('MAX')
       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_max,mpi_comm_world,error)
    case default
       call io_error('Unknown operation in comms_allreduce_real')

    end select

    call dcopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_allreduce_real')
    end if

    if (allocated(array_red)) deallocate(array_red)
#endif

    return

  end subroutine comms_allreduce_real

  subroutine comms_allreduce_cmplx(array,size,op)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error,ierr

    complex(kind=dp), allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       call io_error('failure to allocate array_red in comms_allreduce_cmplx')
    end if

    select case(op)

    case ('SUM')
       call MPI_allreduce(array,array_red,size,MPI_double_complex,MPI_sum,mpi_comm_world,error)
    case ('PRD')
       call MPI_allreduce(array,array_red,size,MPI_double_complex,MPI_prod,mpi_comm_world,error)
    case default
       call io_error('Unknown operation in comms_allreduce_cmplx')

    end select

    call zcopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_allreduce_cmplx')
    end if

    if (allocated(array_red)) deallocate(array_red)
#endif

    return

  end subroutine comms_allreduce_cmplx

  ! Array: local array for sending data; localcount elements will be sent
  !        to the root node
  ! rootglobalarray: array on the root node to which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split
  subroutine comms_gatherv_real(array,localcount,rootglobalarray,counts,displs)

    implicit none

    real(kind=dp), intent(inout)              :: array
    integer, intent(in)                       :: localcount
    real(kind=dp), intent(inout)              :: rootglobalarray
    integer, dimension(num_nodes), intent(in) :: counts
    integer, dimension(num_nodes), intent(in) :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array,localcount,MPI_double_precision,rootglobalarray,counts,&
         displs,MPI_double_precision,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_gatherv_real')
    end if

#else
    call dcopy(localcount,array,1,rootglobalarray,1)
#endif

    return

  end subroutine comms_gatherv_real


  ! Array: local array for getting data; localcount elements will be fetched
  !        from the root node
  ! rootglobalarray: array on the root node from which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split
  subroutine comms_scatterv_real(array,localcount,rootglobalarray,counts,displs)

    implicit none

    real(kind=dp), intent(inout)              :: array
    integer, intent(in)                       :: localcount
    real(kind=dp), intent(inout)              :: rootglobalarray
    integer, dimension(num_nodes), intent(in) :: counts
    integer, dimension(num_nodes), intent(in) :: displs

#ifdef MPI
    integer :: error

!    call MPI_scatterv(array,localcount,MPI_double_precision,rootglobalarray,counts,&
!         displs,MPI_double_precision,root_id,mpi_comm_world,error)
    call MPI_scatterv(rootglobalarray,counts,displs,MPI_double_precision,&
         array,localcount,MPI_double_precision,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_scatterv_real')
    end if

#else
    call dcopy(localcount,rootglobalarray,1,array,1)
#endif

    return

  end subroutine comms_scatterv_real

  ! Array: local array for getting data; localcount elements will be fetched
  !        from the root node
  ! rootglobalarray: array on the root node from which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split
  subroutine comms_scatterv_int(array,localcount,rootglobalarray,counts,displs)

    implicit none

    integer, intent(inout)                    :: array
    integer, intent(in)                       :: localcount
    integer, intent(inout)                    :: rootglobalarray
    integer, dimension(num_nodes), intent(in) :: counts
    integer, dimension(num_nodes), intent(in) :: displs

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray,counts,displs,MPI_Integer,&
         Array,localcount,MPI_Integer,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_scatterv_real')
    end if

#else
    call my_icopy(localcount,rootglobalarray,1,array,1)
#endif

    return

  end subroutine comms_scatterv_int


end module w90_comms


subroutine my_ICOPY(N,ZX,INCX,ZY,INCY)
  !     .. Scalar Arguments ..
  integer INCX,INCY,N
  !     ..
  !     .. Array Arguments ..
  integer ZX(*),ZY(*)
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
  integer I,IX,IY
  !     ..
  if (N.le.0) return
  if (INCX.eq.1 .and. INCY.eq.1) GO TO 20
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
  IX = 1
  IY = 1
  if (INCX.lt.0) IX = (-N+1)*INCX + 1
  if (INCY.lt.0) IY = (-N+1)*INCY + 1
  do I = 1,N
     ZY(IY) = ZX(IX)
     IX = IX + INCX
     IY = IY + INCY
  end do
  return
  !
  !        code for both increments equal to 1
  !
20 do I = 1,N
     ZY(I) = ZX(I)
  end do
  return
end subroutine my_ICOPY
