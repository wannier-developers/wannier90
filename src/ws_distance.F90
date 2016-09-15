!-*- mode: F90 -*-!
!                                                            !
! Copyright (C) 2016 Lorenzo Paulatto                        !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

module w90_ws_distance
!     if( allocated( wdist_ndeg ) ) then
!        deallocate( wdist_ndeg, stat=ierr  )
!        if (ierr/=0) call io_error('Error in deallocating wdist_ndeg in param_dealloc')
!     end if
!     if( allocated( wdist_shiftj_wsi ) ) then
!        deallocate( wdist_shiftj_wsi, stat=ierr  )
!        if (ierr/=0) call io_error('Error in deallocating wdist_shiftj_wsi in param_dealloc')
!     end if

  use w90_constants, only : dp
  use w90_parameters, only : use_ws_distance

  implicit none

  private
  !
  public :: ws_translate_dist, clean_ws_translate
  !
  ! The number of unit cells to shift WF j to put its centre inside the Wigner-Seitz
  ! of wannier function i. If several shifts are equivalent (i.e. they take the function
  ! on the edge of the WS) they are all listed
  integer, public, save, allocatable :: wdist_shiftj_wsi(:,:,:,:,:)!(3,ndegenx,num_wann,num_wann,nrpts)
  ! The number of equivalent shifts (see above)
  integer, public, save, allocatable :: wdist_ndeg(:,:,:)!(num_wann,num_wann,nrpts)
  !
  ! next parameter moved to parameters, used here
  !logical, save, public :: use_ws_distance = .false.
  logical, save :: done_ws_distance = .false.
  integer, parameter :: ndegenx = 8 ! max number of unit cells that can touch
                                    ! in a single point (i.e.  vertex of cube)

   ! How far shold we look? It depends on how far the WFs have been wandering
   ! from their unit-cell, or if the cell is very slanted. 2 is almost always enough,
   ! 3 is enough in every case I have ever met
   integer,parameter :: far = 3

  contains

! The next three subroutines find the supercell translation (i.e. the translation
! by a integer number of supercell) That minimizes the distance between two given funtions,
! i and j, the first in unit cell 0, the other in unit cell R.
! I.e. we put  WF j in the Wigner-Seitz of WF i.
! We also look for the number of equivalent translation, meaning that j is on the edge of the WS
! The results are stored in global arrays wdist_ndeg and wdist_shiftj_wsi
!====================================================!
subroutine ws_translate_dist(nrpts, irvec, force_recompute)
    !====================================================!

    use w90_parameters, only : num_wann,wannier_centres, real_lattice, &
                               recip_lattice, iprint
                                !translation_centre_frac, automatic_translation,lenconfac
    use w90_io,         only : stdout, io_error
    use w90_utility,    only : utility_cart_to_frac, utility_frac_to_cart

    implicit none

    integer, intent(in) :: nrpts
    integer, intent(in) :: irvec(3,nrpts)
    logical, optional, intent(in):: force_recompute ! set to true to force recomputing everything

    ! <<<local variables>>>
    integer  :: iw, jw, ideg, ir, ierr
    real(DP) :: wdist_wssc_frac(3,num_wann,num_wann,nrpts)
    real(DP) :: wdist_deg_cart(3,ndegenx,num_wann,num_wann,nrpts)
    real(DP) :: wdist_deg_frac(3,ndegenx,num_wann,num_wann,nrpts)
    real(DP) :: irvec_cart(3)

    ! The subroutine does nothing if called more than once, which may
    ! not be the best thing if you invoke it while the WFs are moving
    if(present(force_recompute))then
      if(force_recompute) then
        call clean_ws_translate()
      endif
    endif
    if(done_ws_distance) return
    done_ws_distance = .true.

   if (ndegenx*num_wann*nrpts<=0) call io_error("unexpected dimensions in ws_translate_dist")

   allocate(wdist_shiftj_wsi(3,ndegenx,num_wann,num_wann,nrpts),stat=ierr)
   if (ierr/=0) call io_error('Error in allocating wdist_shiftj_wsi in ws_translate_dist')
   allocate(wdist_ndeg(num_wann,num_wann,nrpts),stat=ierr)
   if (ierr/=0) call io_error('Error in allocating wcenter_ndeg in ws_translate_dist')

    !translation_centre_frac = 0._dp
    wdist_ndeg   = 0
    wdist_shiftj_wsi   = 0
    wdist_deg_cart = 0
    wdist_deg_frac = 0
    ! take to WS cell
!     write(stdout,'(1x,a)') 'Translated centres'
!     write(stdout,'(4x,a,3f10.6)') 'translation centre in fractional coordinate:',translation_centre_frac(:)
    do ir =1,nrpts
    do jw=1,num_wann
        do iw=1,num_wann
            call utility_frac_to_cart(DBLE(irvec(:,ir)),irvec_cart,real_lattice)
            ! function IW translated in the Wigner-Size around function JW
            wdist_wssc_frac(:,iw,jw,ir) = R_wz_sc( -wannier_centres(:,iw)&
                        +(irvec_cart+wannier_centres(:,jw)), (/0._dp,0._dp,0._dp/) )
            !find its degeneracy
            CALL R_wz_sc_equiv(wdist_wssc_frac(:,iw,jw,ir), (/0._dp,0._dp,0._dp/), &
                            wdist_ndeg(iw,jw,ir), wdist_deg_cart(:,:,iw,jw,ir))
            IF(wdist_ndeg(iw,jw,ir)>ndegenx) call io_error('surprising ndeg')
            do ideg = 1,wdist_ndeg(iw,jw,ir)
            wdist_deg_cart(:,ideg,iw,jw,ir) = wdist_deg_cart(:,ideg,iw,jw,ir)&
                            +wannier_centres(:,iw)-wannier_centres(:,jw)
            call utility_cart_to_frac(wdist_deg_cart(:,ideg,iw,jw,ir),&
                            wdist_deg_frac(:,ideg,iw,jw,ir),recip_lattice)
            enddo
        enddo
    enddo
    enddo
    ! apply translation
    wdist_shiftj_wsi = NINT(wdist_deg_frac)
    !
    IF(iprint>3)then
      write(stdout,'(1x,a78)') repeat('-',78)
      do ir=1,nrpts
      write(stdout,'(i5)') ir
      do iw=1,num_wann
          write(stdout,'("deg:",100i2)') wdist_ndeg(:,iw,ir)
      enddo
      enddo
      write(stdout,'(1x,a78)') repeat('-',78)
    endif
    !
    IF(ANY(ABS(DBLE(wdist_shiftj_wsi)-wdist_deg_frac)>1.d-6)) &
    call io_error('wrong wdist_shiftj_wsi')

    return
end subroutine ws_translate_dist

! puts R_in in the Wigner-Seitz cell centered around R0
!====================================================!
function R_wz_sc(R_in, R0) result (R_bz)
    use w90_parameters, only : real_lattice,recip_lattice, mp_grid
    use w90_utility,    only : utility_cart_to_frac,utility_frac_to_cart
    use w90_io,         only : stdout
    implicit none
    real(DP),intent(in) :: R_in(3), R0(3)
    real(DP) :: R(3), R_bz(3), R_f(3), R_in_f(3), mod2_R_bz
    integer :: i,j,k

    R_bz = R_in
    mod2_R_bz = SUM((R_bz-R0)**2)
    !
    ! take R_bz to cryst(frac) coord for translating
    call utility_cart_to_frac(R_in,R_in_f,recip_lattice)

    do i = -far, far
    do j = -far, far
        do k = -far, far

        R_f = R_in_f + REAL( (/i*mp_grid(1),j*mp_grid(2),k*mp_grid(3)/), kind=DP)
        call utility_frac_to_cart(R_f,R,real_lattice)

        if(SUM((R-R0)**2)<mod2_R_bz) then
            R_bz = R
            mod2_R_bz = SUM((R_bz-R0)**2)
        endif

        enddo
    enddo
    enddo
end function R_wz_sc
!
! Find the list list of R_out that differ from R_in by a lattice vector
! and are equally distant from R0 (i.e. that are on the edges of the WS
! cell centered on R0)
!====================================================!
subroutine R_wz_sc_equiv(R_in, R0, ndeg, R_out)
    use w90_parameters, only : real_lattice,recip_lattice, mp_grid
    use w90_utility,    only : utility_cart_to_frac,utility_frac_to_cart
    use w90_io,         only : stdout
    implicit none
    real(DP),intent(in)  :: R_in(3), R0(3)
    real(DP),intent(out) :: R_out(3,ndegenx)
    integer,intent(out)  :: ndeg

    real(DP) :: R(3), R_f(3), R_in_f(3), mod2_R_bz
    integer :: i,j,k
    integer,parameter :: far = 3
    real(DP),parameter :: eps = 1.e-6_dp !d-1

    ! init
    ndeg=0
    R_out = 0._dp

    mod2_R_bz = SUM((R_in-R0)**2)
    if(mod2_R_bz<eps)then
    ndeg=1
    return
    endif
    !
    ! take R_bz to cryst(frac) coord for translating
    call utility_cart_to_frac(R_in,R_in_f,recip_lattice)

    do i = -far, far
    do j = -far, far
        do k = -far, far

        R_f = R_in_f + REAL( (/i*mp_grid(1),j*mp_grid(2),k*mp_grid(3)/), kind=DP)
        call utility_frac_to_cart(R_f,R,real_lattice)

        if(  ABS(SUM((R-R0)**2)-mod2_R_bz)/mod2_R_bz<eps) then
            ndeg=ndeg+1
            R_out(:,ndeg) = R
        endif

        enddo
    enddo
    enddo
    !====================================================!
end subroutine R_wz_sc_equiv
!====================================================!
!====================================================!
subroutine clean_ws_translate()
    !====================================================!
    implicit none
    done_ws_distance = .false.
    if(allocated(wdist_shiftj_wsi)) deallocate(wdist_shiftj_wsi)
    if(allocated(wdist_ndeg)) deallocate(wdist_ndeg)
  !====================================================!
end subroutine

end module w90_ws_distance
