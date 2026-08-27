!!-*- mode: F90 -*-!
!------------------------------------------------------------!
! Copyright (C) 2026 Wannier Developer Group                 !
!                                                            !
! This library is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU Lesser General !
! Public License as published by the Free Software           !
! Foundation; either version 2.1 of the License, or (at your !
! option) any later version.                                 !
!                                                            !
! This library is distributed in the hope that it will be    !
! useful,but WITHOUT ANY WARRANTY; without even the implied  !
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR    !
! PURPOSE.  See the GNU Lesser General Public License for    !
! more details.                                              !
!                                                            !
! You should have received a copy of the GNU Lesser General  !
! Public License along with this library; if not, see        !
! <https://www.gnu.org/licenses/>.                           !
!                                                            !
! The webpage of the Wannier90 code is                       !
! <https://www.wannier.org>.                                 !
!                                                            !
! The Wannier90 code is hosted on GitHub                     !
! <https://github.com/wannier-developers/wannier90>          !
!------------------------------------------------------------!

module w90_library_c
!Fortran 2018: Assumed-length character dummy argument ‘keyword’ at (1) of procedure ‘cset_option_int’ with BIND(C) attribute
  use iso_c_binding
  use w90_library, only: w90_disentangle_f => w90_disentangle, &
                         w90_get_centres_f => w90_get_centres, w90_get_gkpb_f => w90_get_gkpb, &
                         w90_get_nn_f => w90_get_nn, w90_get_nnkp_f => w90_get_nnkp, &
                         w90_get_proj_f => w90_get_proj, &
                         w90_get_num_excl_bands_f => w90_get_num_excl_bands, &
                         w90_get_excl_bands_f => w90_get_excl_bands, &
                         w90_print_info_f => w90_print_info, &
                         w90_get_spreads_f => w90_get_spreads, w90_input_setopt_ff => w90_input_setopt, &
                         w90_input_reader_f => w90_input_reader, &
                         w90_project_overlap_f => w90_project_overlap, &
                         w90_set_eigval_f => w90_set_eigval, w90_set_m_local_f => w90_set_m_local, &
                         w90_set_u_matrix_f => w90_set_u_matrix, w90_set_u_opt_f => w90_set_u_opt, &
                         w90_wannierise_f => w90_wannierise, w90_set_option, &
                         w90_set_comm_ff => w90_set_comm, lib_common_type, &
                         w90_get_fortran_stderr, w90_get_fortran_stdout
  implicit none

  public

  type, bind(c) :: w90_data
    type(c_ptr) :: caddr
  end type

contains

  subroutine w90_create(w90_obj) bind(c)
    !! return a c-pointer to a instance of the wannier90 library data structure
    type(lib_common_type), pointer :: common_data
    type(w90_data) :: w90_obj
    !if (c_associated(w90_obj%caddr)) return ! we can't distinguish an uninitialised pointer (junk) from valid one here; test is reliable.
    ! -> policy: this function *always* creates a new object
    allocate (common_data)
    w90_obj%caddr = c_loc(common_data)
  end subroutine

  subroutine w90_delete(w90_obj) bind(c)
    !! deallocates/clears a c-pointer to a instance of the wannier90 library data structure
    implicit none
    type(w90_data) :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    if (.not. c_associated(w90_obj%caddr)) return
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    deallocate (w90_fptr)
    w90_obj%caddr = C_NULL_PTR
  end subroutine

  subroutine w90_get_nk(w90_obj, n) bind(c)
    ! return the number of k-points
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: n
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: ndat
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(n, ndat)
    ndat = w90_fptr%num_kpts
  end subroutine

  subroutine w90_get_nw(w90_obj, n) bind(c)
    ! return the number of wfs
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: n
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: ndat
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(n, ndat)
    ndat = w90_fptr%num_wann
  end subroutine

  subroutine w90_get_nn(w90_obj, n, ierr) bind(c)
    ! return the number of adjacent k-points in finite difference scheme
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: n
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: ndat
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(n, ndat)
    call w90_get_nn_f(w90_fptr, ndat, istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_gkpb(w90_obj, gkpb, ierr) bind(c)
    ! return the g-offset of adjacent k-points in finite difference scheme
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: gkpb
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: nfptr(:, :, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(gkpb, nfptr, [3, w90_fptr%num_kpts, w90_fptr%kmesh_info%nntot])
    call w90_get_gkpb_f(w90_fptr, nfptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_nnkp(w90_obj, nnkp, ierr) bind(c)
    ! return the indexing of adjacent k-points in finite difference scheme
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: nnkp
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: nfptr(:, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(nnkp, nfptr, [w90_fptr%num_kpts, w90_fptr%kmesh_info%nntot])
    call w90_get_nnkp_f(w90_fptr, nfptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_input_setopt_f(w90_obj, seedname, ierr) bind(c)
    ! specify parameters through the library interface
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: seedname
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_input_setopt_ff(w90_fptr, seedname, istdout, istderr, ierr)
  end subroutine

  subroutine w90_input_reader(w90_obj, ierr) bind(c)
    ! read (optional) parameters from .win file
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_input_reader_f(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_project_overlap(w90_obj, ierr) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_project_overlap_f(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_print_info(w90_obj, ierr) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_print_info_f(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_set_eigval(w90_obj, eigval_cptr) bind(c)
    ! copy a pointer to eigenvalue data
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: eigval_cptr
    real(kind=c_double), pointer :: eigval_fptr(:, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(eigval_cptr, eigval_fptr, [w90_fptr%num_bands, w90_fptr%num_kpts])
    call w90_set_eigval_f(w90_fptr, eigval_fptr)
  end subroutine

  subroutine w90_set_m_local(w90_obj, m_cptr) bind(c)
    ! copy a pointer to m-matrix data
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: m_cptr
    complex(kind=c_double_complex), pointer :: fptr(:, :, :, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(m_cptr, fptr, [w90_fptr%num_bands, w90_fptr%num_bands, w90_fptr%kmesh_info%nntot, w90_fptr%num_kpts])
    call w90_set_m_local_f(w90_fptr, fptr)
  end subroutine

  subroutine w90_set_u_matrix(w90_obj, a_cptr) bind(c)
    ! copy pointer to u-matrix
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: a_cptr
    complex(kind=c_double_complex), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(a_cptr, a_fptr, [w90_fptr%num_wann, w90_fptr%num_wann, w90_fptr%num_kpts]) ! these are reversed wrt c
    call w90_set_u_matrix_f(w90_fptr, a_fptr)
  end subroutine

  subroutine w90_set_u_opt(w90_obj, a_cptr) bind(c)
    ! copy pointer to u-matrix (also used for initial projections)
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: a_cptr
    complex(kind=c_double_complex), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(a_cptr, a_fptr, [w90_fptr%num_bands, w90_fptr%num_wann, w90_fptr%num_kpts]) ! these are reversed wrt c
    call w90_set_u_opt_f(w90_fptr, a_fptr)
  end subroutine

  subroutine w90_disentangle(w90_obj, ierr) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_disentangle_f(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_wannierise(w90_obj, ierr) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_wannierise_f(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_centres(w90_obj, centres) bind(c)
    ! returns the centres of calulated mlwfs
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: centres
    type(lib_common_type), pointer :: w90_fptr
    real(kind=c_double), pointer :: fcentres(:, :)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(centres, fcentres, [3, w90_fptr%num_wann])
    call w90_get_centres_f(w90_fptr, fcentres)
  end subroutine

  subroutine w90_get_spreads(w90_obj, spreads) bind(c)
    ! returns the spreads of calulated mlwfs
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: spreads
    type(lib_common_type), pointer :: w90_fptr
    real(kind=c_double), pointer :: fspreads(:)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(spreads, fspreads, [w90_fptr%num_wann])
    call w90_get_spreads_f(w90_fptr, fspreads)
  end subroutine

  subroutine w90_get_nproj(w90_obj, n) bind(c)
    ! return the number of projectors
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: n
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: n_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(n, n_fptr)
    n_fptr = size(w90_fptr%proj_input)
  end subroutine

  subroutine w90_get_proj(w90_obj, n, site, l, m, s, rad, x, z, sqa, zona, ierr) bind(c)
    ! probes projectors configured in the library object
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: n, site, l, m, s, rad, x, z, sqa, zona
    integer(kind=c_int) :: ierr
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: n_fptr, l_fptr(:), m_fptr(:), s_fptr(:), rad_fptr(:)
    real(kind=c_double), pointer :: site_fptr(:, :), x_fptr(:, :), z_fptr(:, :), sqa_fptr(:, :), zona_fptr(:)
    integer(kind=c_int) :: nproj
    integer(kind=c_int) :: istderr, istdout
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(n, n_fptr)
    nproj = n_fptr
    call c_f_pointer(site, site_fptr, [3, nproj])
    call c_f_pointer(l, l_fptr, [nproj])
    call c_f_pointer(m, m_fptr, [nproj])
    call c_f_pointer(s, s_fptr, [nproj])
    call c_f_pointer(rad, rad_fptr, [nproj])
    call c_f_pointer(x, x_fptr, [3, nproj])
    call c_f_pointer(z, z_fptr, [3, nproj])
    call c_f_pointer(sqa, sqa_fptr, [3, nproj])
    call c_f_pointer(zona, zona_fptr, [nproj])
    call w90_get_proj_f(w90_fptr, n_fptr, site_fptr, l_fptr, m_fptr, s_fptr, rad_fptr, x_fptr, z_fptr, sqa_fptr, zona_fptr, &
                        istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_num_excl_bands(w90_obj, num_excl_bands) bind(c)
    ! return the number of excluded bands
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: num_excl_bands
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: n_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(num_excl_bands, n_fptr)
    call w90_get_num_excl_bands_f(w90_fptr, n_fptr)
  end subroutine

  subroutine w90_get_excl_bands(w90_obj, excl_bands, ierr) bind(c)
    ! return the excluded bands
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: excl_bands
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: excl_bands_fptr(:)
    integer(kind=c_int), allocatable :: temp_excl(:)
    integer(kind=c_int) :: istderr, istdout, ierr
    integer :: nex
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_get_num_excl_bands_f(w90_fptr, nex)
    if (nex == 0) then
      ! no excluded bands, just return an empty array
      return
    end if
    allocate (temp_excl(nex))
    call w90_get_excl_bands_f(w90_fptr, temp_excl, istdout, istderr, ierr)
    ! copy to the output pointer
    call c_f_pointer(excl_bands, excl_bands_fptr, [nex])
    excl_bands_fptr = temp_excl
    deallocate (temp_excl)
  end subroutine

  subroutine w90_set_option_double_f(w90_obj, keyword, cdble) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    real(kind=c_double), value  :: cdble
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, cdble)
  end subroutine

  subroutine w90_set_option_double1d_f(w90_obj, keyword, arg_cptr, x) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    type(c_ptr), value  :: arg_cptr
    integer(kind=c_int), value :: x
    type(lib_common_type), pointer :: w90_fptr
    real(kind=c_double), pointer :: fptr(:)
    call c_f_pointer(arg_cptr, fptr, [x])
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_double2d_f(w90_obj, keyword, arg_cptr, x, y) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    type(c_ptr), value  :: arg_cptr
    integer(kind=c_int), value :: x, y
    type(lib_common_type), pointer :: w90_fptr
    real(kind=c_double), pointer :: fptr(:, :)
    call c_f_pointer(arg_cptr, fptr, [y, x]) ! these are reversed wrt c
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_int_f(w90_obj, keyword, cint) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: cint
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, cint)
  end subroutine

  subroutine w90_set_option_int1d_f(w90_obj, keyword, arg_cptr, x) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value ::  arg_cptr
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: x
    integer(kind=c_int), pointer :: fptr(:)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(arg_cptr, fptr, [x])
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_int2d_f(w90_obj, keyword, arg_cptr, x, y) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value ::  arg_cptr
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: x, y
    integer(kind=c_int), pointer :: fptr(:, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(arg_cptr, fptr, [y, x]) ! these are reversed wrt c
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_logical_f(w90_obj, keyword, boolarg) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    logical(kind=c_bool), value  :: boolarg
    logical :: fbool
    type(lib_common_type), pointer :: w90_fptr
    fbool = boolarg
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, fbool)
  end subroutine

  subroutine w90_set_option_text_f(w90_obj, keyword, text) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    character(*, kind=c_char) :: text
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, text)
  end subroutine

  subroutine w90_set_option_text2d_f(w90_obj, keyword, text) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char), intent(in) :: keyword
    ! len=* automatically picks up the 'max_len' we set in CFI_establish
    character(len=*, kind=c_char), intent(in) :: text(:)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, text)
  end subroutine

#ifdef MPI
  subroutine w90_set_comm_f(w90_obj, comm) bind(c)
#ifdef MPI08
    use mpi_f08
#endif
    implicit none
#ifdef MPIH
    include 'mpif.h'
#endif
    integer(kind=c_int), intent(in) :: comm
#ifdef MPI08
    type(mpi_comm) :: comm08
#endif
    type(w90_data), intent(in), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr

    call c_f_pointer(w90_obj%caddr, w90_fptr)

#ifdef MPI08
    ! Manually assign the integer to the type's internal handle
    comm08%MPI_VAL = comm
    call w90_set_comm_ff(w90_fptr, comm08)
#else
    call w90_set_comm_ff(w90_fptr, comm)
#endif
  end subroutine
#endif

end module w90_library_c
