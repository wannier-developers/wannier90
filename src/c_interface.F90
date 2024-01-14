
module w90_library_c

!Fortran 2018: Assumed-length character dummy argument ‘keyword’ at (1) of procedure ‘cset_option_int’ with BIND(C) attribute
  use iso_c_binding
  use w90_library
  implicit none
  public
contains

!  type(c_ptr) function getglob(seedname) bind(c)
!    type(lib_common_type), pointer :: common_data
!    character(*, kind=c_char) :: seedname
!    allocate (common_data)
!    common_data%seedname = seedname
!    getglob = c_loc(common_data)
!  end function

  type(c_ptr) function getglob() bind(c)
    type(lib_common_type), pointer :: common_data
    allocate (common_data)
    getglob = c_loc(common_data)
  end function

  ! fixme needs "destructor"

  type(c_ptr) function getwann() bind(c)
    type(lib_wannier_type), pointer :: wannier_data
    allocate (wannier_data)
    getwann = c_loc(wannier_data)
  end function
  ! fixme needs "destructor"
  subroutine ccreate_kmesh(common_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call create_kmesh(common_fptr, istdout, istderr, ierr)
  end subroutine ccreate_kmesh

  subroutine cchkpt(common_cptr, wannier_cptr, label, seedname) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    character(*, kind=c_char) :: seedname, label
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call write_chkpt(common_fptr, wannier_fptr, label, seedname, istdout, istderr, ierr)
  end subroutine cchkpt

  subroutine cdisentangle(common_cptr, wannier_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call disentangle(common_fptr, wannier_fptr, istdout, istderr, ierr)
  end subroutine cdisentangle

  subroutine cinput_reader(common_cptr, wannier_cptr, seedname) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr
    character(*, kind=c_char) :: seedname
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call input_reader(common_fptr, wannier_fptr, seedname, istdout, istderr, ierr)
  end subroutine cinput_reader

  subroutine cinput_reader_special(common_cptr, wannier_cptr, seedname) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr
    character(*, kind=c_char) :: seedname
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call input_reader_special(common_fptr, wannier_fptr, seedname, istdout, istderr, ierr)
  end subroutine cinput_reader_special

  subroutine cinput_setopt(common_cptr, wannier_cptr, seedname, comm) bind(c)
#ifdef MPI08
    use mpi_f08
#endif
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr
    character(*, kind=c_char) :: seedname
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
#ifdef MPI08
    type(mpi_comm), value :: comm
#else
    integer(kind=c_int), value :: comm
#endif
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call input_setopt(common_fptr, wannier_fptr, seedname, comm, istdout, istderr, ierr)
  end subroutine cinput_setopt

  subroutine coverlaps(common_cptr, wannier_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call overlaps(common_fptr, wannier_fptr, istdout, istderr, ierr)
  end subroutine coverlaps

  subroutine cset_a_matrix(common_cptr, wannier_cptr, a_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr, a_cptr
    complex(kind=dp), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call c_f_pointer(a_cptr, a_fptr, [common_fptr%num_bands, common_fptr%num_wann, common_fptr%num_kpts])
    call set_a_matrix(wannier_fptr, a_fptr)
  end subroutine cset_a_matrix

  subroutine cset_exclude(common_cptr, exclude, nexc) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, exclude
    integer(kind=c_int), pointer :: exclude_fptr(:)
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), value :: nexc
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(exclude, exclude_fptr, [nexc])
    call set_exclude(common_fptr, exclude_fptr, nexc)
  end subroutine cset_exclude

  subroutine cset_eigval(common_cptr, eigval_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, eigval_cptr
    real(kind=dp), pointer :: eigval_fptr(:, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(eigval_cptr, eigval_fptr, [common_fptr%num_bands, common_fptr%num_kpts])
    call set_eigval(common_fptr, eigval_fptr)
  end subroutine cset_eigval

  subroutine cset_kpoint_distribution(common_cptr, dist_cptr) bind(c)
    ! this needs to be called after reading input in order to have num_kpts available
    ! fixme!!
    implicit none
    type(c_ptr), value :: common_cptr, dist_cptr
    integer(kind=c_int), pointer :: dist_fptr(:)
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(dist_cptr, dist_fptr, [common_fptr%num_kpts])
    call set_kpoint_distribution(common_fptr, dist_fptr, istdout, istderr, ierr)
  end subroutine cset_kpoint_distribution

  subroutine cset_m_matrix_local(common_cptr, wannier_cptr, m_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr, m_cptr
    complex(kind=dp), pointer :: m_fptr(:, :, :, :)
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr

    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call c_f_pointer(m_cptr, m_fptr, [common_fptr%num_wann, common_fptr%num_wann, &
                                      common_fptr%kmesh_info%nntot, common_fptr%num_kpts])
    !fixme jj call set_m_matrix_local(wannier_fptr, m_fptr)
  end subroutine cset_m_matrix_local

  subroutine cset_m_orig(common_cptr, wannier_cptr, m_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr, m_cptr
    complex(kind=dp), pointer :: m_fptr(:, :, :, :)
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call c_f_pointer(m_cptr, m_fptr, [common_fptr%num_bands, common_fptr%num_bands, &
                                      common_fptr%kmesh_info%nntot, common_fptr%num_kpts])
    call set_m_orig(wannier_fptr, m_fptr)
  end subroutine cset_m_orig

  subroutine cset_option_int(common_cptr, keyword, ival) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: ival
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call set_option(common_fptr, keyword, ival)
  end subroutine cset_option_int

  subroutine cset_option_int3(common_cptr, keyword, ival) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    type(c_ptr) ::  ival
    integer(kind=c_int), pointer :: m_fptr(:)
    character(*, kind=c_char) :: keyword
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(ival, m_fptr, [3])
    call set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_int3

  subroutine cset_option_intx(common_cptr, keyword, ival, ict) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    type(c_ptr), value ::  ival
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: ict
    integer(kind=c_int), pointer :: m_fptr(:)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(ival, m_fptr, [ict])
    call set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_intx

  subroutine cset_option_float(common_cptr, keyword, ival) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    real(kind=c_double), value  :: ival
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call set_option(common_fptr, keyword, ival)
  end subroutine cset_option_float

  subroutine cset_option_float33(common_cptr, keyword, ival) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    type(c_ptr)  :: ival
    type(lib_common_type), pointer :: common_fptr
    real(kind=dp), pointer :: m_fptr(:, :)
    call c_f_pointer(ival, m_fptr, [3, 3])
    call c_f_pointer(common_cptr, common_fptr)
    call set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_float33

  subroutine cset_option_floatxy(common_cptr, keyword, ival, x, y) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    type(c_ptr), value  :: ival
    integer(kind=c_int), value :: x, y
    type(lib_common_type), pointer :: common_fptr
    real(kind=dp), pointer :: m_fptr(:, :)
    call c_f_pointer(ival, m_fptr, [y, x]) ! haha check!! fixme jj
    call c_f_pointer(common_cptr, common_fptr)
    call set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_floatxy

  subroutine cset_parallel_comms(common_cptr, comm) bind(c)
    use w90_comms, only: w90_comm_type
#ifdef MPI08
    use mpi_f08
#endif
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
#ifdef MPI08
    type(mpi_comm), value :: comm
#else
    integer(kind=c_int), value :: comm
#endif
    call c_f_pointer(common_cptr, common_fptr)
    call set_parallel_comms(common_fptr, comm)
  end subroutine cset_parallel_comms

  subroutine cset_u_matrix(common_cptr, wannier_cptr, a_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr, a_cptr
    complex(kind=dp), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call c_f_pointer(a_cptr, a_fptr, [common_fptr%num_wann, common_fptr%num_wann, &
                                      common_fptr%num_kpts])
    call set_u_matrix(common_fptr, a_fptr)
  end subroutine cset_u_matrix

  subroutine cget_nn(common_cptr, n) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, n
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), pointer :: ndat
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(n, ndat)
    call get_nn(common_fptr, ndat)
  end subroutine cget_nn

  subroutine cget_nnkp(common_cptr, nnkp) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, nnkp
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), pointer :: nfptr(:, :)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(nnkp, nfptr, [common_fptr%num_kpts, common_fptr%kmesh_info%nntot])
    call get_nnkp(common_fptr, nfptr)
  end subroutine cget_nnkp

  subroutine cget_spreads(common_cptr, spreads) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, spreads
    type(lib_common_type), pointer :: common_fptr
    real(kind=dp), pointer :: fspreads(:)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(spreads, fspreads, [common_fptr%num_wann])
    call get_spreads(common_fptr, fspreads)
  end subroutine cget_spreads

  subroutine cget_centres(common_cptr, centres) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, centres
    type(lib_common_type), pointer :: common_fptr
    real(kind=dp), pointer :: fcentres(:, :)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(centres, fcentres, [3, common_fptr%num_wann])
    call get_centres(common_fptr, fcentres)
  end subroutine cget_centres

  subroutine cset_u_opt(common_cptr, wannier_cptr, a_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr, a_cptr
    complex(kind=dp), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call c_f_pointer(a_cptr, a_fptr, [common_fptr%num_bands, common_fptr%num_wann, &
                                      common_fptr%num_kpts])
    call set_u_opt(common_fptr, a_fptr)
  end subroutine cset_u_opt

  subroutine cread_eigvals(common_cptr, eig_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, eig_cptr
    type(lib_common_type), pointer :: common_fptr
    real(kind=dp), pointer :: eig_fptr(:, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(eig_cptr, eig_fptr, [common_fptr%num_bands, common_fptr%num_kpts])
    call read_eigvals(common_fptr, eig_fptr, istdout, istderr, ierr)
  end subroutine cread_eigvals

  subroutine cwannierise(common_cptr, wannier_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr, wannier_cptr
    type(lib_common_type), pointer :: common_fptr
    type(lib_wannier_type), pointer :: wannier_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(wannier_cptr, wannier_fptr)
    call wannierise(common_fptr, wannier_fptr, istdout, istderr, ierr)
  end subroutine cwannierise
end module w90_library_c
