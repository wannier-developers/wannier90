
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

  subroutine cdisentangle(common_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call disentangle(common_fptr, istdout, istderr, ierr)
  end subroutine cdisentangle

  subroutine cwannierise(common_cptr) bind(c)
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call wannierise(common_fptr, istdout, istderr, ierr)
  end subroutine cwannierise

  subroutine cinput_setopt(common_cptr, seedname, comm) bind(c)
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
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call input_setopt(common_fptr, seedname, comm, istdout, istderr, ierr)
  end subroutine cinput_setopt

  subroutine cinput_reader(common_cptr) bind(c)
! read (optional) parameters from .win file
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call input_reader(common_fptr, istdout, istderr, ierr)
  end subroutine cinput_reader

  subroutine cinput_reader_special(common_cptr, seedname) bind(c)
! read (non-optional/special) parameters from .win file
    implicit none
    type(c_ptr), value :: common_cptr
    character(*, kind=c_char) :: seedname
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call input_reader_special(common_fptr, seedname, istdout, istderr, ierr)
  end subroutine cinput_reader_special

  subroutine cset_eigval(common_cptr, eigval_cptr) bind(c)
! copy a pointer to eigenvalue data
    implicit none
    type(c_ptr), value :: common_cptr, eigval_cptr
    real(kind=dp), pointer :: eigval_fptr(:, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(eigval_cptr, eigval_fptr, [common_fptr%num_bands, common_fptr%num_kpts])
    call set_eigval(common_fptr, eigval_fptr)
  end subroutine cset_eigval

  subroutine cset_m_matrix_local(common_cptr, m_cptr) bind(c)
! copy a pointer to m-matrix data
    implicit none
    type(c_ptr), value :: common_cptr, m_cptr
    complex(kind=dp), pointer :: m_fptr(:, :, :, :)
    type(lib_common_type), pointer :: common_fptr

    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(m_cptr, m_fptr, [common_fptr%num_wann, common_fptr%num_wann, &
                                      common_fptr%kmesh_info%nntot, common_fptr%num_kpts])
  end subroutine cset_m_matrix_local

  subroutine cset_m_orig(common_cptr, m_cptr) bind(c)
! copy a pointer to m-matrix data
    implicit none
    type(c_ptr), value :: common_cptr, m_cptr
    complex(kind=dp), pointer :: m_fptr(:, :, :, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(m_cptr, m_fptr, [common_fptr%num_bands, common_fptr%num_bands, &
                                      common_fptr%kmesh_info%nntot, common_fptr%num_kpts])
    call set_m_orig(common_fptr, m_fptr)
  end subroutine cset_m_orig

  subroutine cset_u_matrix(common_cptr, a_cptr) bind(c)
! copy pointer to u-matrix
    implicit none
    type(c_ptr), value :: common_cptr, a_cptr
    complex(kind=dp), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(a_cptr, a_fptr, [common_fptr%num_wann, common_fptr%num_wann, &
                                      common_fptr%num_kpts]) ! these are reversed wrt c
    call set_u_matrix(common_fptr, a_fptr)
  end subroutine cset_u_matrix

  subroutine cset_u_opt(common_cptr, a_cptr) bind(c)
! copy pointer to u-matrix (also used for initial projections)
    implicit none
    type(c_ptr), value :: common_cptr, a_cptr
    complex(kind=dp), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: common_fptr
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(a_cptr, a_fptr, [common_fptr%num_bands, common_fptr%num_wann, &
                                      common_fptr%num_kpts]) ! these are reversed wrt c
    call set_u_opt(common_fptr, a_fptr)
  end subroutine cset_u_opt

  subroutine cget_nn(common_cptr, n) bind(c)
! return the number of adjacent k-points in finite difference scheme
    implicit none
    type(c_ptr), value :: common_cptr, n
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), pointer :: ndat
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(n, ndat)
    call get_nn(common_fptr, istdout, istderr, ndat, ierr)
  end subroutine cget_nn

  subroutine cget_nnkp(common_cptr, nnkp) bind(c)
! return the indexing of adjacent k-points in finite difference scheme
    implicit none
    type(c_ptr), value :: common_cptr, nnkp
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int), pointer :: nfptr(:, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(nnkp, nfptr, [common_fptr%num_kpts, common_fptr%kmesh_info%nntot])
    call get_nnkp(common_fptr, istdout, istderr, nfptr, ierr)
  end subroutine cget_nnkp

  subroutine cget_spreads(common_cptr, spreads) bind(c)
! returns the spreads of calulated mlwfs
    implicit none
    type(c_ptr), value :: common_cptr, spreads
    type(lib_common_type), pointer :: common_fptr
    real(kind=dp), pointer :: fspreads(:)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(spreads, fspreads, [common_fptr%num_wann])
    call get_spreads(common_fptr, fspreads)
  end subroutine cget_spreads

  subroutine cget_centres(common_cptr, centres) bind(c)
! returns the centres of calulated mlwfs
    implicit none
    type(c_ptr), value :: common_cptr, centres
    type(lib_common_type), pointer :: common_fptr
    real(kind=dp), pointer :: fcentres(:, :)
    call c_f_pointer(common_cptr, common_fptr)
    call c_f_pointer(centres, fcentres, [3, common_fptr%num_wann])
    call get_centres(common_fptr, fcentres)
  end subroutine cget_centres

  subroutine cread_eigvals(common_cptr, eig_cptr) bind(c)
! read eigenvalues from convetional .eig file (used by standalone wannier90.x)
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

  subroutine coverlaps(common_cptr) bind(c)
! read m matrix from conventional .mmn file (used by standalone wannier90.x)
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call overlaps(common_fptr, istdout, istderr, ierr)
  end subroutine coverlaps

  subroutine cchkpt(common_cptr, label) bind(c)
! writes a standard checkpoint file with specified label (used by standalone wannier90.x)
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    character(*, kind=c_char) :: label
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call write_chkpt(common_fptr, label, istdout, istderr, ierr)
  end subroutine cchkpt

  subroutine cset_parallel_comms(common_cptr, comm) bind(c)
! alternative method for specifying mpi communicator
! passing communicator to input_setopt() is anyway mandatory, so this is not needed
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

  subroutine cset_kpoint_distribution(common_cptr, dist_cptr) bind(c)
! alternative method for specifying k-point decompostion
! (using set_option is likely preferable)
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

  subroutine ccreate_kmesh(common_cptr) bind(c)
! cause the creation of the internal k-mesh data structures
! this is automatically done prior to use, so this call is generally unnecessary
    implicit none
    type(c_ptr), value :: common_cptr
    type(lib_common_type), pointer :: common_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call get_fortran_stderr(istderr)
    call get_fortran_stdout(istdout)
    call c_f_pointer(common_cptr, common_fptr)
    call create_kmesh(common_fptr, istdout, istderr, ierr)
  end subroutine ccreate_kmesh

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
    call c_f_pointer(ival, m_fptr, [y, x]) ! these are reversed wrt c
    call c_f_pointer(common_cptr, common_fptr)
    call set_option(common_fptr, keyword, m_fptr)
  end subroutine cset_option_floatxy

end module w90_library_c
