!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
!  COMMS: set of MPI wrapper                                 !
! (c) 2006-2007 Jonathan R. Yates                            !
!                                                            !
!------------------------------------------------------------!


module w90_comms

  use w90_constants, only : dp

  implicit none


  private

#ifdef MPI
  include 'mpif.h'
#endif

  logical, public :: on_root
  integer, public :: num_nodes,my_node_id
  integer, public, parameter :: root_id=0

  integer, parameter :: mpi_send_tag=77 !abitrary

  public :: comms_setup
  public :: comms_end
  public :: comms_bcast      ! send data from the root node
  public :: comms_send       ! send data from one node to another
  public :: comms_recv       ! accept data from one node to another
  public :: comms_reduce     ! reduce data onto root node (n.b. not allreduce) 

  interface comms_bcast
     module procedure comms_bcast_int
     module procedure comms_bcast_logical
     module procedure comms_bcast_real
     module procedure comms_bcast_cmplx
     module procedure comms_bcast_char
  end interface

  interface comms_send
     module procedure comms_send_int
     module procedure comms_send_logical
     module procedure comms_send_real
     module procedure comms_send_cmplx
     module procedure comms_send_char
  end interface

  interface comms_recv
     module procedure comms_recv_int
     module procedure comms_recv_logical
     module procedure comms_recv_real
     module procedure comms_recv_cmplx
     module procedure comms_recv_char
  end interface

  interface comms_reduce
!     module procedure comms_reduce_int !not working but trivial to fix
     module procedure comms_reduce_real
     module procedure comms_reduce_cmplx
  end interface


contains

  subroutine comms_setup
 
    implicit none

    integer :: ierr


#ifdef MPI
    call mpi_init(ierr)
    if (ierr.ne.0) stop 'MPI initialisation error'
    call mpi_comm_rank(mpi_comm_world, my_node_id, ierr)
    call mpi_comm_size(mpi_comm_world, num_nodes, ierr)
#else
    num_nodes=1
    my_node_id=0
#endif

    on_root=.false.
    if(my_node_id==root_id) on_root=.true.
    
  end subroutine comms_setup

  subroutine comms_end
 
    implicit none

    integer :: ierr


#ifdef MPI
    call mpi_finalize(ierr)
#else
    STOP
#endif
    
  end subroutine comms_end

  subroutine comms_bcast_int(array,size)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size

    integer :: error

#ifdef MPI

    call MPI_bcast(array,size,MPI_integer,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_bcast_int'
       call comms_error
    end if
#endif

    return

  end subroutine comms_bcast_int

  subroutine comms_bcast_real(array,size)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size

    integer :: error

#ifdef MPI

    call MPI_bcast(array,size,MPI_double_precision,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_bcast_real'
       call comms_error
    end if
#endif

    return

  end subroutine comms_bcast_real

  subroutine comms_bcast_logical(array,size)

    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size

    integer :: error

#ifdef MPI

    call MPI_bcast(array,size,MPI_logical,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_bcast_logical'
       call comms_error
    end if
#endif

    return

  end subroutine comms_bcast_logical

  subroutine comms_bcast_char(array,size)

    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size

    integer :: error

#ifdef MPI

    call MPI_bcast(array,size,MPI_character,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_bcast_char'
       call comms_error
    end if
#endif

    return

  end subroutine comms_bcast_char

  subroutine comms_bcast_cmplx(array,size)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size

    integer :: error

#ifdef MPI

    call MPI_bcast(array,size,MPI_double_complex,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_bcast_cmplx'
       call comms_error
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

    integer :: error

#ifdef MPI

    call MPI_send(array,size,MPI_logical,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_send_logical'
       call comms_error
    end if
#endif

    return

  end subroutine comms_send_logical


  subroutine comms_send_int(array,size,to)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

    integer :: error

#ifdef MPI

    call MPI_send(array,size,MPI_integer,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_send_int'
       call comms_error
    end if
#endif

    return

  end subroutine comms_send_int


  subroutine comms_send_char(array,size,to)

    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

    integer :: error

#ifdef MPI

    call MPI_send(array,size,MPI_character,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_send_char'
       call comms_error
    end if
#endif

    return

  end subroutine comms_send_char


  subroutine comms_send_real(array,size,to)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

    integer :: error

#ifdef MPI

    call MPI_send(array,size,MPI_double_precision,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_send_real'
       call comms_error
    end if
#endif

    return

  end subroutine comms_send_real


  subroutine comms_send_cmplx(array,size,to)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: to

    integer :: error

#ifdef MPI

    call MPI_send(array,size,MPI_double_complex,to, &
         mpi_send_tag,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_send_cmplx'
       call comms_error
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

    integer :: error

#ifdef MPI
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_logical,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_recv_logical'
       call comms_error
    end if
#endif

    return

  end subroutine comms_recv_logical


  subroutine comms_recv_int(array,size,from)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

    integer :: error

#ifdef MPI
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_integer,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_recv_int'
       call comms_error
    end if
#endif

    return

  end subroutine comms_recv_int


  subroutine comms_recv_char(array,size,from)

    implicit none

    character(len=*), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

    integer :: error

#ifdef MPI
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_character,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_recv_char'
       call comms_error
    end if
#endif

    return

  end subroutine comms_recv_char


  subroutine comms_recv_real(array,size,from)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

    integer :: error

#ifdef MPI
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_double_precision,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_recv_real'
       call comms_error
    end if
#endif

    return

  end subroutine comms_recv_real


  subroutine comms_recv_cmplx(array,size,from)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    integer, intent(in)    :: from

    integer :: error

#ifdef MPI

    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_double_complex,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if(error.ne.MPI_success) then
       print*,'Error in comms_recv_cmplx'
       call comms_error
    end if

#endif

    return

  end subroutine comms_recv_cmplx


  subroutine comms_error

    implicit none
    
    integer :: error
    
#ifdef MPI

    call MPI_abort(MPI_comm_world,1,error)

#endif

  end subroutine comms_error

  
  ! COMMS_REDUCE (collect data on the root node)

  subroutine comms_reduce_int(array,size,op)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    integer :: error,ierr

#ifdef MPI

    integer :: status(MPI_status_size)
    integer, allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       print*,'failure to allocate array_red in comms_reduce_int'
       call comms_error
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_integer,MPI_sum,0,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_integer,MPI_prod,0,mpi_comm_world,error)
    case default
       print*,'Unknown operation in comms_reduce_int'
       call comms_error

    end select

!    call icopy(size,array,1,array_red,1)

! For the moment I'm lazy and this routine doesn't work
! just need to write an integer BLAS copy routine (trivial)

    if(error.ne.MPI_success) then
       print*,'Error in comms_reduce_real'
       call comms_error
    end if
#endif

    return

  end subroutine comms_reduce_int


  subroutine comms_reduce_real(array,size,op)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    integer :: error,ierr

#ifdef MPI

    integer :: status(MPI_status_size)
    real(kind=dp), allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       print*,'failure to allocate array_red in comms_reduce_real'
       call comms_error
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_sum,0,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_prod,0,mpi_comm_world,error)
    case default
       print*,'Unknown operation in comms_reduce_real'
       call comms_error

    end select

    call dcopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       print*,'Error in comms_reduce_real'
       call comms_error
    end if
#endif

    return

  end subroutine comms_reduce_real


  subroutine comms_reduce_cmplx(array,size,op)

    implicit none

    complex(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op
    integer :: error,ierr

#ifdef MPI

    integer :: status(MPI_status_size)
    complex(kind=dp), allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       print*,'failure to allocate array_red in comms_reduce_cmplx'
       call comms_error
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_double_complex,MPI_sum,0,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_double_complex,MPI_prod,0,mpi_comm_world,error)
    case default
       print*,'Unknown operation in comms_reduce_cmplx'
       call comms_error

    end select

    call zcopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       print*,'Error in comms_reduce_cmplx'
       call comms_error
    end if
#endif

    return

  end subroutine comms_reduce_cmplx


end module w90_comms
