module w90_lib
  implicit none

  private

  public :: read_overlaps_matrix

  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0, 1.0d0))

contains

  subroutine read_overlaps_matrix(mmn_in, overlaps_matrix, seedname, num_bands, num_kpts, nntot, &
                                  nnlist, nncell, mmn_tmp, error_handler)
    integer, intent(in) :: mmn_in
    complex(DPC), intent(out) :: overlaps_matrix(:, :, :, :)
    character(len=50), intent(in) :: seedname
    integer, intent(in) :: num_bands
    integer, intent(in) :: num_kpts
    integer, intent(in) :: nntot
    integer, intent(in) :: nnlist(:, :)   ! list of neighbours for each k-point
    integer, intent(in) :: nncell(:, :, :) ! gives BZ of each neighbour of each k-point
    complex(DPC), intent(out) :: mmn_tmp(:, :)

    interface
      subroutine error_handler(error_message)
        character(len=*), intent(in) :: error_message
      end subroutine
    end interface

    integer :: nkp, nkp2, inn, nn, n, m
    integer :: num_mmn
    integer :: nnl, nnm, nnn, ncount
    integer :: nb_tmp, nkp_tmp, nntot_tmp
    real(kind=dp) :: m_real, m_imag
    ! complex(kind=dp), allocatable :: mmn_tmp(:, :)
    ! character(len=50) :: dummy
    logical :: nn_found

    ! Read the number of bands, k-points and nearest neighbours
    read (mmn_in, *, err=103, end=103) nb_tmp, nkp_tmp, nntot_tmp

    ! Checks
    if (nb_tmp .ne. num_bands) &
      call error_handler(trim(seedname)//'.mmn has not the right number of bands')
    if (nkp_tmp .ne. num_kpts) &
      call error_handler(trim(seedname)//'.mmn has not the right number of k-points')
    if (nntot_tmp .ne. nntot) &
      call error_handler(trim(seedname)//'.mmn has not the right number of nearest neighbours')

    ! Read the overlaps
    num_mmn = num_kpts*nntot
    ! allocate (mmn_tmp(num_bands, num_bands), stat=ierr)
    ! if (ierr /= 0) call error_handler('Error in allocating mmn_tmp in overlap_read')
    do ncount = 1, num_mmn
      read (mmn_in, *, err=103, end=103) nkp, nkp2, nnl, nnm, nnn
      do n = 1, num_bands
        do m = 1, num_bands
          read (mmn_in, *, err=103, end=103) m_real, m_imag
          mmn_tmp(m, n) = cmplx(m_real, m_imag, kind=dp)
        enddo
      enddo
      nn = 0
      nn_found = .false.
      do inn = 1, nntot
        if ((nkp2 .eq. nnlist(nkp, inn)) .and. &
            (nnl .eq. nncell(1, nkp, inn)) .and. &
            (nnm .eq. nncell(2, nkp, inn)) .and. &
            (nnn .eq. nncell(3, nkp, inn))) then
          if (.not. nn_found) then
            nn_found = .true.
            nn = inn
          else
            call error_handler('Error reading '//trim(seedname)// &
                               '.mmn. More than one matching nearest neighbour found')
          endif
        endif
      end do
      if (nn .eq. 0) then
        ! if (on_root) write (stdout, '(/a,i8,2i5,i4,2x,3i3)') &
        !   ' Error reading '//trim(seedname)//'.mmn:', ncount, nkp, nkp2, nn, nnl, nnm, nnn
        call error_handler('Neighbour not found')
      end if

      overlaps_matrix(:, :, nn, nkp) = mmn_tmp(:, :)
    end do

    return

103 call error_handler('Error: Problem reading input file '//trim(seedname)//'.mmn')

  end subroutine read_overlaps_matrix

end module w90_lib
