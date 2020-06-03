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

module w90_spin
  !! Module to compute spin

  use w90_constants, only: dp

  implicit none

  private

  public :: spin_get_moment, spin_get_nk, spin_get_S

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       !
  !===========================================================!

  subroutine spin_get_moment
    !============================================================!
    !                                                            !
    !! Computes the spin magnetic moment by Wannier interpolation
    !                                                            !
    !============================================================!

    use w90_constants, only: dp, pi, cmplx_i
    use w90_comms, only: on_root, my_node_id, num_nodes, comms_reduce
    use w90_io, only: io_error, stdout
    use w90_postw90_common, only: num_int_kpts_on_node, int_kpts, weight
    use w90_parameters, only: spin_kmesh, wanint_kpoint_file, &
      nfermi, fermi_energy_list
    use w90_get_oper, only: get_HH_R, get_SS_R

    integer       :: loop_x, loop_y, loop_z, loop_tot
    real(kind=dp) :: kweight, kpt(3), spn_k(3), spn_all(3), &
                     spn_mom(3), magnitude, theta, phi, conv

    if (nfermi > 1) call io_error('Routine spin_get_moment requires nfermi=1')

    call get_HH_R
    call get_SS_R

    if (on_root) then
      write (stdout, '(/,/,1x,a)') '------------'
      write (stdout, '(1x,a)') 'Calculating:'
      write (stdout, '(1x,a)') '------------'
      write (stdout, '(/,3x,a)') '* Spin magnetic moment'
    end if

    spn_all = 0.0_dp
    if (wanint_kpoint_file) then

      if (on_root) then
        write (stdout, '(/,1x,a)') 'Sampling the irreducible BZ only'
        write (stdout, '(5x,a)') &
          'WARNING: - IBZ implementation is currently limited to simple cases:'
        write (stdout, '(5x,a)') &
          '               Check results against a full BZ calculation!'
      end if
      !
      ! Loop over k-points on the irreducible wedge of the Brillouin zone,
      ! read from file 'kpoint.dat'
      !
      do loop_tot = 1, num_int_kpts_on_node(my_node_id)
        kpt(:) = int_kpts(:, loop_tot)
        kweight = weight(loop_tot)
        call spin_get_moment_k(kpt, fermi_energy_list(1), spn_k)
        spn_all = spn_all + spn_k*kweight
      end do

    else

      if (on_root) &
        write (stdout, '(/,1x,a)') 'Sampling the full BZ (not using symmetry)'
      kweight = 1.0_dp/real(PRODUCT(spin_kmesh), kind=dp)
      do loop_tot = my_node_id, PRODUCT(spin_kmesh) - 1, num_nodes
        loop_x = loop_tot/(spin_kmesh(2)*spin_kmesh(3))
        loop_y = (loop_tot - loop_x*(spin_kmesh(2)*spin_kmesh(3)))/spin_kmesh(3)
        loop_z = loop_tot - loop_x*(spin_kmesh(2)*spin_kmesh(3)) &
                 - loop_y*spin_kmesh(3)
        kpt(1) = (real(loop_x, dp)/real(spin_kmesh(1), dp))
        kpt(2) = (real(loop_y, dp)/real(spin_kmesh(2), dp))
        kpt(3) = (real(loop_z, dp)/real(spin_kmesh(3), dp))
        call spin_get_moment_k(kpt, fermi_energy_list(1), spn_k)
        spn_all = spn_all + spn_k*kweight
      end do

    end if

    ! Collect contributions from all nodes
    !
    call comms_reduce(spn_all(1), 3, 'SUM')

    ! No factor of g=2 because the spin variable spans [-1,1], not
    ! [-1/2,1/2] (i.e., it is really the Pauli matrix sigma, not S)
    !
    spn_mom(1:3) = -spn_all(1:3)

    if (on_root) then
      write (stdout, '(/,1x,a)') 'Spin magnetic moment (Bohr magn./cell)'
      write (stdout, '(1x,a,/)') '===================='
      write (stdout, '(1x,a18,f11.6)') 'x component:', spn_mom(1)
      write (stdout, '(1x,a18,f11.6)') 'y component:', spn_mom(2)
      write (stdout, '(1x,a18,f11.6)') 'z component:', spn_mom(3)

      ! Polar and azimuthal angles of the magnetization (defined as in pwscf)
      !
      conv = 180.0_dp/pi
      magnitude = sqrt(spn_mom(1)**2 + spn_mom(2)**2 + spn_mom(3)**2)
      theta = acos(spn_mom(3)/magnitude)*conv
      phi = atan(spn_mom(2)/spn_mom(1))*conv
      write (stdout, '(/,1x,a18,f11.6)') 'Polar theta (deg):', theta
      write (stdout, '(1x,a18,f11.6)') 'Azim. phi (deg):', phi
    end if

  end subroutine spin_get_moment

! =========================================================================

  subroutine spin_get_nk(kpt, spn_nk)
    !=============================================================!
    !                                                             !
    !! Computes <psi_{mk}^(H)|S.n|psi_{mk}^(H)> (m=1,...,num_wann)
    !! where S.n = n_x.S_x + n_y.S_y + n_z.Z_z
    !!
    !! S_i are the Pauli matrices and n=(n_x,n_y,n_z) is the unit
    !! vector along the chosen spin quantization axis
    !                                                             !
    !============================================================ !

    use w90_constants, only: dp, pi, cmplx_0, cmplx_i
    use w90_io, only: io_error
    use w90_utility, only: utility_diagonalize, utility_rotate_diag
    use w90_parameters, only: num_wann, spin_axis_polar, &
      spin_axis_azimuth
    use w90_postw90_common, only: pw90common_fourier_R_to_k
    use w90_get_oper, only: HH_R, SS_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: spn_nk(num_wann)

    ! Physics
    !
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: SS(:, :, :), SS_n(:, :)

    ! Misc/Dummy
    !
    integer          :: is
    real(kind=dp)    :: eig(num_wann), alpha(3), conv

    allocate (HH(num_wann, num_wann))
    allocate (UU(num_wann, num_wann))
    allocate (SS(num_wann, num_wann, 3))
    allocate (SS_n(num_wann, num_wann))

    call pw90common_fourier_R_to_k(kpt, HH_R, HH, 0)
    call utility_diagonalize(HH, num_wann, eig, UU)

    do is = 1, 3
      call pw90common_fourier_R_to_k(kpt, SS_R(:, :, :, is), SS(:, :, is), 0)
    enddo

    ! Unit vector along the magnetization direction
    !
    conv = 180.0_dp/pi
    alpha(1) = sin(spin_axis_polar/conv)*cos(spin_axis_azimuth/conv)
    alpha(2) = sin(spin_axis_polar/conv)*sin(spin_axis_azimuth/conv)
    alpha(3) = cos(spin_axis_polar/conv)

    ! Vector of spin matrices projected along the quantization axis
    !
    SS_n(:, :) = alpha(1)*SS(:, :, 1) + alpha(2)*SS(:, :, 2) + alpha(3)*SS(:, :, 3)

    spn_nk(:) = real(utility_rotate_diag(SS_n, UU, num_wann), dp)

  end subroutine spin_get_nk

  !===========================================================!
  !                   PRIVATE PROCEDURES                      !
  !===========================================================!

  subroutine spin_get_moment_k(kpt, ef, spn_k)
    !! Computes the spin magnetic moment by Wannier interpolation
    !! at the specified k-point
    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_io, only: io_error
    use w90_utility, only: utility_diagonalize, utility_rotate_diag
    use w90_parameters, only: num_wann
    use w90_postw90_common, only: pw90common_fourier_R_to_k, pw90common_get_occ
    use w90_get_oper, only: HH_R, SS_R
    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(in)  :: ef
    real(kind=dp), intent(out) :: spn_k(3)

    ! Physics
    !
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: SS(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    real(kind=dp)                 :: spn_nk(num_wann, 3)

    ! Misc/Dummy
    !
    integer          :: i, is
    real(kind=dp)    :: eig(num_wann), occ(num_wann)

    allocate (HH(num_wann, num_wann))
    allocate (UU(num_wann, num_wann))
    allocate (SS(num_wann, num_wann, 3))

    call pw90common_fourier_R_to_k(kpt, HH_R, HH, 0)
    call utility_diagonalize(HH, num_wann, eig, UU)
    call pw90common_get_occ(eig, occ, ef)

    spn_k(1:3) = 0.0_dp
    do is = 1, 3
      call pw90common_fourier_R_to_k(kpt, SS_R(:, :, :, is), SS(:, :, is), 0)
      spn_nk(:, is) = aimag(cmplx_i*utility_rotate_diag(SS(:, :, is), UU, num_wann))
      do i = 1, num_wann
        spn_k(is) = spn_k(is) + occ(i)*spn_nk(i, is)
      end do
    enddo

  end subroutine spin_get_moment_k

  subroutine spin_get_S(kpt, S)
    !===========================================================!
    !                                                           !
    ! Computes <psi_{nk}^(H)|S|psi_{nk}^(H)> (n=1,...,num_wann) !
    ! where S = (S_x,S_y,S_z) is the vector of Pauli matrices   !
    !                                                           !
    !========================================================== !

    use w90_constants, only: dp, pi, cmplx_0, cmplx_i
    use w90_io, only: io_error
    use w90_utility, only: utility_diagonalize, utility_rotate_diag
    use w90_parameters, only: num_wann
    use w90_postw90_common, only: pw90common_fourier_R_to_k
    use w90_get_oper, only: HH_R, SS_R

    ! Arguments
    !
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), intent(out) :: S(num_wann, 3)

    ! Physics
    !
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: SS(:, :, :)
    real(kind=dp)                 :: eig(num_wann)

    ! Misc/Dummy
    !
    integer :: i

    allocate (HH(num_wann, num_wann))
    allocate (UU(num_wann, num_wann))
    allocate (SS(num_wann, num_wann, 3))

    call pw90common_fourier_R_to_k(kpt, HH_R, HH, 0)
    call utility_diagonalize(HH, num_wann, eig, UU)

    do i = 1, 3
      call pw90common_fourier_R_to_k(kpt, SS_R(:, :, :, i), SS(:, :, i), 0)
      S(:, i) = real(utility_rotate_diag(SS(:, :, i), UU, num_wann), dp)
    enddo

  end subroutine spin_get_S

end module w90_spin
