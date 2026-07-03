!-*- mode: F90 -*-!
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
!                                                            !
!  w90_decompose: decomposition of Wannier function          !
!  densities onto an orthonormalised Gaussian-radial x        !
!  real-spherical-harmonic basis (Himanen et al. 2020)        !
!                                                            !
!------------------------------------------------------------!

module w90_decompose

  !! Decompose each Wannier function's density onto a basis of
  !! Löwdin-orthonormalised Gaussian radial functions times real
  !! spherical harmonics, following Himanen et al., Adv. Sci. 7,
  !! 1902333 (2020) (SOAP-style). This module currently provides only
  !! the basis mathematics; density projection, power spectra and file
  !! I/O are added by later work.

  use w90_constants, only: dp, pi

  implicit none

  private

  public :: decompose_radial_params
  public :: decompose_real_ylm

contains

  !================================================================!
  subroutine decompose_radial_params(n_max, l_max, r_min, r_max, alphas, betas)
    !! Precompute the parameters that define the radial basis
    !! \(\phi_{nl}(r) = r^{l}\,\exp(-\alpha_{nl} r^{2})\) and its
    !! Löwdin-orthonormalised counterpart
    !! \(g_{nl}(r) = \sum_{n'} \beta_{n'nl}\,\phi_{n'l}(r)\).
    !!
    !! The decay coefficients follow the analytic decay condition of
    !! Himanen et al. (2020): demanding \(\phi_{nl}(r_{\mathrm{thr},n})
    !! = \mathrm{thr}\) with \(\mathrm{thr} = 10^{-3}\) and
    !! \(r_{\mathrm{thr},n}\) equally spaced on \([r_{\min},r_{\max}]\)
    !! gives \(\alpha_{nl} = -\ln(\mathrm{thr}/r_{\mathrm{thr},n}^{l})
    !! / r_{\mathrm{thr},n}^{2}\).
    !!
    !! The overlap of two radial functions is evaluated analytically,
    !! \(S_{nn'l} = \tfrac12\,\Gamma(l+\tfrac32)\,
    !! (\alpha_{nl}+\alpha_{n'l})^{-(l+3/2)}\), and orthonormalised by
    !! Löwdin symmetric orthogonalisation
    !! \(\beta = V\,\mathrm{diag}(1/\sqrt{e})\,V^{T}\), where \(V\), \(e\)
    !! are the eigenvectors/eigenvalues of \(S\) (LAPACK dsyev).

    implicit none

    integer, intent(in) :: n_max
    !! Number of radial functions
    integer, intent(in) :: l_max
    !! Maximum angular momentum
    real(kind=dp), intent(in) :: r_min
    !! Smallest threshold radius (Angstrom)
    real(kind=dp), intent(in) :: r_max
    !! Largest threshold radius (Angstrom)
    real(kind=dp), intent(out) :: alphas(n_max, 0:l_max)
    !! Gaussian decay coefficients \(\alpha_{nl}\)
    real(kind=dp), intent(out) :: betas(n_max, n_max, 0:l_max)
    !! Löwdin orthonormalisation coefficients \(\beta_{nn'l}\)

    real(kind=dp), parameter :: thr = 1.0e-3_dp

    real(kind=dp) :: r_thr(n_max)
    real(kind=dp) :: smat(n_max, n_max)
    real(kind=dp) :: evec(n_max, n_max)
    real(kind=dp) :: eval(n_max)
    real(kind=dp) :: work(3*n_max)
    real(kind=dp) :: gamma_half, exponent, cc
    integer :: n, np, l, info

    ! Threshold radii: linspace(r_min, r_max, n_max)
    if (n_max == 1) then
      r_thr(1) = r_min
    else
      do n = 1, n_max
        r_thr(n) = r_min + (r_max - r_min)*real(n - 1, dp)/real(n_max - 1, dp)
      end do
    end if

    ! Decay coefficients from the analytic decay condition
    do l = 0, l_max
      do n = 1, n_max
        alphas(n, l) = -(log(thr) - real(l, dp)*log(r_thr(n)))/r_thr(n)**2
      end do
    end do

    ! For each l build the analytic overlap matrix and Löwdin-orthonormalise
    do l = 0, l_max
      ! Gamma(l + 3/2) = sqrt(pi) * prod_{k=0}^{l} (k + 1/2)
      gamma_half = sqrt(pi)
      do n = 0, l
        gamma_half = gamma_half*(real(n, dp) + 0.5_dp)
      end do
      exponent = real(l, dp) + 1.5_dp

      do np = 1, n_max
        do n = 1, n_max
          cc = alphas(n, l) + alphas(np, l)
          smat(n, np) = 0.5_dp*gamma_half*cc**(-exponent)
        end do
      end do

      ! Symmetric eigendecomposition S = V diag(eval) V^T
      evec = smat
      call dsyev('V', 'U', n_max, evec, n_max, eval, work, 3*n_max, info)
      if (info /= 0) then
        ! Leave betas unset on failure; caller-side error handling is
        ! added when this routine is wired into the workflow.
        betas(:, :, l) = 0.0_dp
        cycle
      end if

      ! beta = V diag(1/sqrt(eval)) V^T  (symmetric, sign-independent)
      do np = 1, n_max
        do n = 1, n_max
          betas(n, np, l) = sum(evec(n, :)*evec(np, :)/sqrt(eval(:)))
        end do
      end do
    end do

  end subroutine decompose_radial_params

  !================================================================!
  subroutine decompose_real_ylm(l_max, theta, phi, ylm)
    !! Real spherical harmonics \(Y_{lm}(\theta,\phi)\) for all
    !! \(0 \le l \le l_{\max}\), \(-l \le m \le l\), evaluated by an
    !! associated-Legendre recurrence.
    !!
    !! \(\theta\) is the polar (colatitude) angle in \([0,\pi]\) and
    !! \(\phi\) the azimuthal angle in \([0,2\pi]\), matching the
    !! spherical convention of the python reference
    !! (`koopmans.ml._misc.cart2sph_array`). The real harmonics use the
    !! scipy real convention with the \(\sqrt{2}\,(-1)^{m}\)
    !! Condon-Shortley factor:
    !! \(Y_{l0} = N_{l}^{0} P_{l}^{0}(\cos\theta)\),
    !! \(Y_{lm} = \sqrt{2}\,(-1)^{m} N_{l}^{m} P_{l}^{m}(\cos\theta)
    !! \cos(m\phi)\) for \(m>0\), and
    !! \(Y_{lm} = \sqrt{2}\,(-1)^{|m|} N_{l}^{|m|} P_{l}^{|m|}(\cos\theta)
    !! \sin(|m|\phi)\) for \(m<0\), with
    !! \(N_{l}^{m} = \sqrt{(2l+1)/(4\pi)\,(l-m)!/(l+m)!}\) and
    !! \(P_{l}^{m}\) the Condon-Shortley associated Legendre functions.
    !!
    !! The result is indexed `ylm(l, i)` with `i = 0..2*l` mapping to
    !! \(m = i - l\), i.e. `i=0` is \(m=-l\) and `i=2l` is \(m=+l\).
    !! Entries with `i > 2*l` are left as zero.

    implicit none

    integer, intent(in) :: l_max
    !! Maximum angular momentum
    real(kind=dp), intent(in) :: theta
    !! Polar (colatitude) angle in [0, pi]
    real(kind=dp), intent(in) :: phi
    !! Azimuthal angle in [0, 2*pi]
    real(kind=dp), intent(out) :: ylm(0:l_max, 0:2*l_max)
    !! Real spherical harmonics; ylm(l, m+l)

    real(kind=dp) :: x, plm, norm, ang, sqrt2
    integer :: l, m

    x = cos(theta)
    sqrt2 = sqrt(2.0_dp)
    ylm = 0.0_dp

    do l = 0, l_max
      do m = 0, l
        plm = norm_assoc_legendre(l, m, x)   ! N_l^m * P_l^m(x), CS phase included
        if (m == 0) then
          ylm(l, l) = plm
        else
          ! (-1)^m
          norm = merge(-1.0_dp, 1.0_dp, mod(m, 2) == 1)
          ang = real(m, dp)*phi
          ! m > 0  ->  index l + m
          ylm(l, l + m) = sqrt2*norm*plm*cos(ang)
          ! m < 0  ->  index l - m
          ylm(l, l - m) = sqrt2*norm*plm*sin(ang)
        end if
      end do
    end do

  end subroutine decompose_real_ylm

  !================================================================!
  pure function norm_assoc_legendre(l, m, x) result(res)
    !! Normalised associated Legendre function
    !! \(N_{l}^{m} P_{l}^{m}(x)\) for \(0 \le m \le l\), with
    !! \(N_{l}^{m} = \sqrt{(2l+1)/(4\pi)\,(l-m)!/(l+m)!}\) and the
    !! Condon-Shortley phase included in \(P_{l}^{m}\), so that the
    !! result equals the modulus part of scipy's complex spherical
    !! harmonic \(Y_{l}^{m}\). Uses the standard stable upward
    !! recurrence in \(l\).

    implicit none

    integer, intent(in) :: l
    !! Degree
    integer, intent(in) :: m
    !! Order (0 <= m <= l)
    real(kind=dp), intent(in) :: x
    !! Argument, x = cos(theta), |x| <= 1
    real(kind=dp) :: res

    real(kind=dp) :: pmm, pmmp1, pll, somx2, fact, norm
    integer :: i, ll

    ! P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2)   (Condon-Shortley phase)
    pmm = 1.0_dp
    if (m > 0) then
      somx2 = sqrt(max(0.0_dp, (1.0_dp - x)*(1.0_dp + x)))
      fact = 1.0_dp
      do i = 1, m
        pmm = -pmm*fact*somx2
        fact = fact + 2.0_dp
      end do
    end if

    if (l == m) then
      pll = pmm
    else
      ! P_{m+1}^m(x) = x (2m+1) P_m^m(x)
      pmmp1 = x*real(2*m + 1, dp)*pmm
      if (l == m + 1) then
        pll = pmmp1
      else
        ! Upward recurrence in l
        do ll = m + 2, l
          pll = (x*real(2*ll - 1, dp)*pmmp1 - real(ll + m - 1, dp)*pmm)/real(ll - m, dp)
          pmm = pmmp1
          pmmp1 = pll
        end do
      end if
    end if

    ! Normalisation N_l^m = sqrt((2l+1)/(4 pi) * (l-m)!/(l+m)!)
    norm = real(2*l + 1, dp)/(4.0_dp*pi)
    do i = l - m + 1, l + m
      norm = norm/real(i, dp)
    end do
    res = sqrt(norm)*pll

  end function norm_assoc_legendre

end module w90_decompose
