
module w90_library_c

!Fortran 2018: Assumed-length character dummy argument ‘keyword’ at (1) of procedure ‘cset_option_int’ with BIND(C) attribute
  use iso_c_binding
  use w90_library
  implicit none
  public
contains
  type(c_ptr) function w90_create() bind(c)
! return a c-pointer to a instance of the wannier90 library data structure
    type(lib_common_type), pointer :: common_data
    allocate (common_data)
    w90_create = c_loc(common_data)
  end function

  subroutine w90_delete(common_cptr) bind(c)
! deallocates/clears a c-pointer to a instance of the wannier90 library data structure
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    deallocate (common_fptr)
    common_cptr = C_NULL_PTR
  end subroutine w90_delete

  subroutine cdisentangle(common_cptr, ierr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call w90_disentangle(common_fptr, istdout, istderr, ierr)
  end subroutine cdisentangle

  subroutine cwannierise(common_cptr, ierr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call w90_wannierise(common_fptr, istdout, istderr, ierr)
  end subroutine cwannierise

  subroutine cinput_setopt(common_cptr, seedname, ierr, comm) bind(c)
! specify parameters through the library interface
#ifdef MPI08
    use mpi_f08
#endif
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: seedname
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
#ifdef MPI08
    type(mpi_comm), value :: comm
#else
    integer(kind=c_int), value :: comm
#endif
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call w90_input_setopt(common_fptr, seedname, comm, istdout, istderr, ierr)
  end subroutine cinput_setopt

  subroutine cinput_reader(common_cptr, ierr) bind(c)
! read (optional) parameters from .win file
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call w90_input_reader(common_fptr, istdout, istderr, ierr)
  end subroutine cinput_reader

  subroutine cset_eigval(common_cptr, eigval_cptr) bind(c)
! copy a pointer to eigenvalue data
    implicit none
    type(c_ptr), value :: common_cptr, eigval_cptr
    real(8), pointer :: eigval_fptr(:, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(eigval_cptr, eigval_fptr, [common_fptr%num_bands, common_fptr%num_kpts])
    call w90_set_eigval(common_fptr, eigval_fptr)
  end subroutine cset_eigval

  subroutine cset_m_local(common_cptr, m_cptr) bind(c)
! copy a pointer to m-matrix data
    implicit none
    type(c_ptr), value :: common_cptr, m_cptr
    complex(8), pointer :: m_fptr(:, :, :, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(m_cptr, m_fptr, [common_fptr%num_bands, common_fptr%num_bands, &
                                      common_fptr%kmesh_info%nntot, common_fptr%num_kpts])
    call w90_set_m_local(common_fptr, m_fptr)
  end subroutine cset_m_local

  subroutine cset_u_matrix(common_cptr, a_cptr) bind(c)
! copy pointer to u-matrix
    implicit none
    type(c_ptr), value :: common_cptr, a_cptr
    complex(8), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(a_cptr, a_fptr, [common_fptr%num_wann, common_fptr%num_wann, &
                                      common_fptr%num_kpts]) ! these are reversed wrt c
    call w90_set_u_matrix(common_fptr, a_fptr)
  end subroutine cset_u_matrix

  subroutine cset_u_opt(common_cptr, a_cptr) bind(c)
! copy pointer to u-matrix (also used for initial projections)
    implicit none
    type(c_ptr), value :: common_cptr, a_cptr
    complex(kind=8), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(a_cptr, a_fptr, [common_fptr%num_bands, common_fptr%num_wann, &
                                      common_fptr%num_kpts]) ! these are reversed wrt c
    call w90_set_u_opt(common_fptr, a_fptr)
  end subroutine cset_u_opt

  subroutine cget_nn(common_cptr, n) bind(c)
! return the number of adjacent k-points in finite difference scheme
    implicit none
    type(c_ptr), value :: common_cptr, n
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), pointer :: ndat
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(n, ndat)
    call w90_get_nn(common_fptr, ndat, istdout, istderr, ierr)
  end subroutine cget_nn

  subroutine cget_nnkp(common_cptr, nnkp) bind(c)
! return the indexing of adjacent k-points in finite difference scheme
    implicit none
    type(c_ptr), value :: common_cptr, nnkp
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), pointer :: nfptr(:, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(nnkp, nfptr, [common_fptr%num_kpts, common_fptr%kmesh_info%nntot])
    call w90_get_nnkp(common_fptr, nfptr, istdout, istderr, ierr)
  end subroutine cget_nnkp

  subroutine cget_gkpb(common_cptr, gkpb) bind(c)
! return the g-offset of adjacent k-points in finite difference scheme
    implicit none
    type(c_ptr), value :: common_cptr, gkpb
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), pointer :: nfptr(:, :, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(gkpb, nfptr, [3, common_fptr%num_kpts, common_fptr%kmesh_info%nntot])
    call w90_get_gkpb(common_fptr, nfptr, istdout, istderr, ierr)
  end subroutine cget_gkpb

  subroutine cget_spreads(common_cptr, spreads) bind(c)
! returns the spreads of calulated mlwfs
    implicit none
    type(c_ptr), value :: common_cptr, spreads
    type(lib_common_type), pointer :: common_fptr
    real(kind=8), pointer :: fspreads(:)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(spreads, fspreads, [common_fptr%num_wann])
    call w90_get_spreads(common_fptr, fspreads)
  end subroutine cget_spreads

  subroutine cget_centres(common_cptr, centres) bind(c)
! returns the centres of calulated mlwfs
    implicit none
    type(c_ptr), value :: common_cptr, centres
    type(lib_common_type), pointer :: common_fptr
    real(kind=8), pointer :: fcentres(:, :)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(centres, fcentres, [3, common_fptr%num_wann])
    call w90_get_centres(common_fptr, fcentres)
  end subroutine cget_centres

  subroutine cset_option_int(common_cptr, keyword, ival) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: ival
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call w90_set_option(common_fptr, keyword, ival)
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
    call w90_set_option(common_fptr, keyword, m_fptr)
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
    call w90_set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_intx

  subroutine cset_option_float(common_cptr, keyword, ival) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    real(kind=c_double), value  :: ival
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call w90_set_option(common_fptr, keyword, ival)
  end subroutine cset_option_float

  subroutine cset_option_float33(common_cptr, keyword, ival) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    type(c_ptr)  :: ival
    type(lib_common_type), pointer :: common_fptr
    real(kind=8), pointer :: m_fptr(:, :)
    call c_f_pointer(ival, m_fptr, [3, 3])
    call c_f_pointer(common_cptr, common_fptr)
    call w90_set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_float33

  subroutine cset_option_floatxy(common_cptr, keyword, ival, x, y) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: keyword
    type(c_ptr), value  :: ival
    integer(kind=c_int), value :: x, y
    type(lib_common_type), pointer :: common_fptr
    real(kind=8), pointer :: m_fptr(:, :)
    call c_f_pointer(ival, m_fptr, [y, x]) ! these are reversed wrt c
    call c_f_pointer(common_cptr, common_fptr)
    call w90_set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_floatxy

end module w90_library_c
