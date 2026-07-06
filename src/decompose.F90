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
  !! 1902333 (2020) (SOAP-style). The module provides the basis
  !! mathematics (\ref decompose_radial_params, \ref decompose_real_ylm),
  !! the density projection (\ref decompose_project), the orbital-orbital
  !! power spectrum (\ref decompose_power_orb) and a root-only driver
  !! (\ref decompose_main) that decomposes a set of Wannier densities and
  !! writes the coefficient/power files.

  use w90_constants, only: dp, pi

  implicit none

  private

  public :: decompose_radial_params
  public :: decompose_real_ylm
  public :: decompose_project
  public :: decompose_power_orb
  public :: decompose_main

contains

  !================================================================!
  subroutine decompose_radial_params(n_max, l_max, r_min, r_max, alphas, betas, error, comm)
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

    use w90_comms, only: w90_comm_type
    use w90_error, only: w90_error_type, set_error_fatal

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
    type(w90_error_type), allocatable, intent(out) :: error
    !! Error handler
    type(w90_comm_type), intent(in) :: comm
    !! Communicator

    real(kind=dp), parameter :: thr = 1.0e-3_dp

    real(kind=dp) :: r_thr(n_max)
    real(kind=dp) :: smat(n_max, n_max)
    real(kind=dp) :: evec(n_max, n_max)
    real(kind=dp) :: eval(n_max)
    real(kind=dp) :: work(3*n_max)
    real(kind=dp) :: gamma_half, exponent, cc
    integer :: n, np, l, info
    character(len=128) :: mesg

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
        write (mesg, '(a,i0,a,i0)') &
          'decompose_radial_params: dsyev failed for l=', l, ', info=', info
        call set_error_fatal(error, trim(mesg), comm)
        return
      end if
      if (any(eval <= 0.0_dp)) then
        write (mesg, '(a,i0,a)') &
          'decompose_radial_params: overlap matrix for l=', l, &
          ' is not positive definite (basis too ill-conditioned)'
        call set_error_fatal(error, trim(mesg), comm)
        return
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

  !================================================================!
  pure function decompose_num_coeff(n_max, l_max) result(nc)
    !! Number of expansion coefficients,
    !! \(n_{\max}\,(l_{\max}+1)^{2} = n_{\max}\sum_{l=0}^{l_{\max}}(2l+1)\).
    implicit none
    integer, intent(in) :: n_max, l_max
    integer :: nc
    nc = n_max*(l_max + 1)**2
  end function decompose_num_coeff

  !================================================================!
  pure function decompose_num_power(n_max, l_max) result(np_)
    !! Number of orbital-orbital power-spectrum entries: one per
    !! \((n_1 \le n_2, l)\) triple, i.e.
    !! \((l_{\max}+1)\,n_{\max}(n_{\max}+1)/2\).
    implicit none
    integer, intent(in) :: n_max, l_max
    integer :: np_
    np_ = (l_max + 1)*n_max*(n_max + 1)/2
  end function decompose_num_power

  !================================================================!
  subroutine decompose_project(rho, ngx, ngy, ngz, nxx_lo, nxx_hi, nyy_lo, nyy_hi, &
                               nzz_lo, nzz_hi, real_lattice, centre, n_max, l_max, &
                               r_cut, alphas, betas, coeff)
    !! Project a real-space density onto the orthonormalised
    !! Gaussian-radial x real-spherical-harmonic basis about a single
    !! centre, returning the flat coefficient vector
    !! \(c_{nlm} = \int \rho(\mathbf{r})\,g_{nl}(|\mathbf{r}-\mathbf{c}|)\,
    !! Y_{lm}(\theta,\phi)\,\mathrm{d}V\).
    !!
    !! The density is sampled on the regular grid that spans the
    !! Born-von-Karman supercell defined by the index bounds
    !! (`nxx_lo:nxx_hi`, ...) with `ngx/ngy/ngz` points per unit cell.
    !! Grid point \((i_x,i_y,i_z)\) sits at fractional coordinate
    !! \((i_x/n_{gx}, i_y/n_{gy}, i_z/n_{gz})\) of the (primitive) cell.
    !! The Cartesian displacement from the centre is taken with the
    !! minimum-image convention over the supercell, and only points with
    !! \(|\mathbf{r}-\mathbf{c}| \le r_{\mathrm{cut}}\) contribute
    !! (spherical mask). The volume element is
    !! \(\mathrm{d}V = |\det(A)| / (n_{gx} n_{gy} n_{gz})\) with \(A\) the
    !! primitive real lattice, which equals |det(supercell)|/(total grid
    !! points). No normalisation of `rho` is applied here.
    !!
    !! Angles follow the python reference convention
    !! (`koopmans.ml._misc.cart2sph_array`):
    !! \(\theta = \mathrm{atan2}(z, \sqrt{x^2+y^2}) + \pi/2\),
    !! \(\phi = \mathrm{atan2}(y, x) + \pi\). The flat ordering is outer
    !! \(n\) (0..n_max-1), then \(l\) (0..l_max), then inner \(m\) (0..2l):
    !! `coeff((n-1)*(l_max+1)**2 + l**2 + m + 1)`.

    use w90_utility, only: utility_inverse_mat, utility_cart_to_frac, &
                           utility_recip_lattice_base

    implicit none

    integer, intent(in) :: nxx_lo, nxx_hi, nyy_lo, nyy_hi, nzz_lo, nzz_hi
    !! Supercell grid index bounds
    real(kind=dp), intent(in) :: rho(nxx_lo:nxx_hi, nyy_lo:nyy_hi, nzz_lo:nzz_hi)
    !! Real-space density on the supercell grid (x contiguous)
    integer, intent(in) :: ngx, ngy, ngz
    !! Grid points per (primitive) unit cell along each lattice vector
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    !! Primitive real lattice; `real_lattice(i,:)` is the i-th vector (Angstrom)
    real(kind=dp), intent(in) :: centre(3)
    !! Cartesian coordinates of the projection centre (Angstrom)
    integer, intent(in) :: n_max
    !! Number of radial functions
    integer, intent(in) :: l_max
    !! Maximum angular momentum
    real(kind=dp), intent(in) :: r_cut
    !! Spherical cutoff radius (Angstrom)
    real(kind=dp), intent(in) :: alphas(n_max, 0:l_max)
    !! Gaussian decay coefficients (from decompose_radial_params)
    real(kind=dp), intent(in) :: betas(n_max, n_max, 0:l_max)
    !! Löwdin coefficients (from decompose_radial_params)
    real(kind=dp), intent(out) :: coeff(n_max*(l_max + 1)**2)
    !! Flat expansion coefficients

    real(kind=dp) :: inv_lat(3, 3), cfrac(3)
    real(kind=dp) :: ylm(0:l_max, 0:2*l_max)
    real(kind=dp) :: phinl(n_max), rpow(0:l_max)
    real(kind=dp) :: dvol, det_cell, r_cut2, r2, rr, xyproj, theta, phi
    real(kind=dp) :: fx, fy, fz, dcart(3), gnl, weight, rhoval, recip(3, 3)
    integer :: nsc(3), ix, iy, iz, n, np, l, m, ibase
    real(kind=dp) :: volume

    coeff = 0.0_dp

    ! Volume element: |det(primitive lattice)| / (points per cell)
    call utility_inverse_mat(real_lattice, inv_lat)
    call utility_recip_lattice_base(real_lattice, recip, volume)
    det_cell = volume
    dvol = det_cell/real(ngx*ngy*ngz, dp)

    ! Supercell replication factors implied by the grid bounds
    nsc(1) = (nxx_hi - nxx_lo + 1)/ngx
    nsc(2) = (nyy_hi - nyy_lo + 1)/ngy
    nsc(3) = (nzz_hi - nzz_lo + 1)/ngz

    ! Centre in (primitive-cell) fractional coordinates
    call utility_cart_to_frac(centre, cfrac, inv_lat)

    r_cut2 = r_cut*r_cut

    do iz = nzz_lo, nzz_hi
      fz = real(iz, dp)/real(ngz, dp) - cfrac(3)
      fz = fz - real(nsc(3), dp)*anint(fz/real(nsc(3), dp))
      do iy = nyy_lo, nyy_hi
        fy = real(iy, dp)/real(ngy, dp) - cfrac(2)
        fy = fy - real(nsc(2), dp)*anint(fy/real(nsc(2), dp))
        do ix = nxx_lo, nxx_hi
          rhoval = rho(ix, iy, iz)
          fx = real(ix, dp)/real(ngx, dp) - cfrac(1)
          fx = fx - real(nsc(1), dp)*anint(fx/real(nsc(1), dp))

          ! Minimum-image Cartesian displacement from the centre
          dcart(1) = fx*real_lattice(1, 1) + fy*real_lattice(2, 1) + fz*real_lattice(3, 1)
          dcart(2) = fx*real_lattice(1, 2) + fy*real_lattice(2, 2) + fz*real_lattice(3, 2)
          dcart(3) = fx*real_lattice(1, 3) + fy*real_lattice(2, 3) + fz*real_lattice(3, 3)

          r2 = dcart(1)**2 + dcart(2)**2 + dcart(3)**2
          if (r2 > r_cut2) cycle
          rr = sqrt(r2)

          ! Spherical angles (python cart2sph convention)
          xyproj = sqrt(dcart(1)**2 + dcart(2)**2)
          theta = atan2(dcart(3), xyproj) + 0.5_dp*pi
          phi = atan2(dcart(2), dcart(1)) + pi
          call decompose_real_ylm(l_max, theta, phi, ylm)

          ! Powers r^l
          rpow(0) = 1.0_dp
          do l = 1, l_max
            rpow(l) = rpow(l - 1)*rr
          end do

          do l = 0, l_max
            ! Primitive radial functions phi_{n'l}(r) = r^l exp(-alpha_{n'l} r^2)
            do np = 1, n_max
              phinl(np) = rpow(l)*exp(-alphas(np, l)*r2)
            end do
            do n = 1, n_max
              ! g_{nl}(r) = sum_{n'} beta_{n'nl} phi_{n'l}(r)
              gnl = 0.0_dp
              do np = 1, n_max
                gnl = gnl + betas(np, n, l)*phinl(np)
              end do
              weight = rhoval*dvol*gnl
              ibase = (n - 1)*(l_max + 1)**2 + l**2
              do m = 0, 2*l
                coeff(ibase + m + 1) = coeff(ibase + m + 1) + weight*ylm(l, m)
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine decompose_project

  !================================================================!
  subroutine decompose_power_orb(coeff, n_max, l_max, power)
    !! Orbital-orbital rotationally-invariant power spectrum
    !! \(p_{n_1 n_2 l} = \sum_{m} c_{n_1 l m}\,c_{n_2 l m}\) for
    !! \(n_1 \le n_2\) and all \(l\), following
    !! `koopmans.ml._compute_power.compute_power_mat` (orb block only).
    !! Ordering: outer \(n_1\), then \(n_2 \ge n_1\), then \(l\).

    implicit none

    integer, intent(in) :: n_max
    !! Number of radial functions
    integer, intent(in) :: l_max
    !! Maximum angular momentum
    real(kind=dp), intent(in) :: coeff(n_max*(l_max + 1)**2)
    !! Flat expansion coefficients (from decompose_project)
    real(kind=dp), intent(out) :: power((l_max + 1)*n_max*(n_max + 1)/2)
    !! Power-spectrum entries

    integer :: n1, n2, l, m, ib1, ib2, ip
    real(kind=dp) :: s

    ip = 0
    do n1 = 1, n_max
      do n2 = n1, n_max
        do l = 0, l_max
          ib1 = (n1 - 1)*(l_max + 1)**2 + l**2
          ib2 = (n2 - 1)*(l_max + 1)**2 + l**2
          s = 0.0_dp
          do m = 0, 2*l
            s = s + coeff(ib1 + m + 1)*coeff(ib2 + m + 1)
          end do
          ip = ip + 1
          power(ip) = s
        end do
      end do
    end do

  end subroutine decompose_power_orb

  !================================================================!
  subroutine decompose_main(wann_func, ngx, ngy, ngz, nxx_lo, nxx_hi, nyy_lo, nyy_hi, &
                            nzz_lo, nzz_hi, n_wf, wf_index, centres, n_out, out_slot, &
                            real_lattice, n_max, l_max, r_min, r_max, r_cut, n_ext, &
                            ext_centres, seedname, stdout, error, comm)
    !! Root-only driver for the Wannier-density decomposition.
    !!
    !! Each of the `n_wf` Wannier functions supplied in `wann_func`
    !! (already reduced onto root by `plot_build_wannier_grid`) has its
    !! density \(\rho_n = |w_n|^2\) normalised to unit integral over the
    !! supercell. The `n_out` functions selected by `out_slot` (indices
    !! into the last dimension of `wann_func`) are decomposed about their
    !! own centres `centres(:,n)`; the coefficients are written to
    !! `<seed>_NNNNN.coeff` and the orbital-orbital power spectrum to
    !! `<seed>_NNNNN.power`, where NNNNN is the global index
    !! `wf_index(n)`. If `n_ext > 0`, the group density \(\sum_n \rho_n\)
    !! summed over ALL `n_wf` functions (not just the `out_slot` subset)
    !! is additionally decomposed about each external centre
    !! `ext_centres(:,e)` and written to `<seed>_gc_NNNNN.coeff`
    !! (NNNNN = e). The group-density coefficients are linear in
    !! \(\rho\), so a caller can sum them across runs to build
    !! total-density coefficients (see the design notes).
    !!
    !! The radial-basis precomputation and the `r_cut` validation are
    !! collective (all ranks) so their error synchronisation is safe; the
    !! projection and file I/O run on root only.

    use w90_comms, only: w90_comm_type, mpirank
    use w90_error, only: w90_error_type, set_error_input
    use w90_utility, only: utility_recip_lattice_base

    implicit none

    integer, intent(in) :: nxx_lo, nxx_hi, nyy_lo, nyy_hi, nzz_lo, nzz_hi
    !! Supercell grid index bounds
    integer, intent(in) :: n_wf
    !! Number of Wannier functions on the grid
    complex(kind=dp), intent(in) :: wann_func(nxx_lo:nxx_hi, nyy_lo:nyy_hi, nzz_lo:nzz_hi, n_wf)
    !! Complex Wannier functions on the supercell grid (root-valid)
    integer, intent(in) :: ngx, ngy, ngz
    !! Grid points per unit cell along each lattice vector
    integer, intent(in) :: wf_index(n_wf)
    !! Global index of each Wannier function (used in file names)
    real(kind=dp), intent(in) :: centres(3, n_wf)
    !! Cartesian centres of the Wannier functions (Angstrom)
    integer, intent(in) :: n_out
    !! Number of Wannier functions to decompose about their own centres
    integer, intent(in) :: out_slot(n_out)
    !! Indices (into the last dimension of `wann_func`) of the functions
    !! to decompose about their own centres
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    !! Primitive real lattice; `real_lattice(i,:)` is the i-th vector (Angstrom)
    integer, intent(in) :: n_max
    !! Number of radial functions
    integer, intent(in) :: l_max
    !! Maximum angular momentum
    real(kind=dp), intent(in) :: r_min
    !! Smallest threshold radius (Angstrom)
    real(kind=dp), intent(in) :: r_max
    !! Largest threshold radius (Angstrom)
    real(kind=dp), intent(in) :: r_cut
    !! Spherical cutoff radius (Angstrom)
    integer, intent(in) :: n_ext
    !! Number of external centres for the group-density channel (0 if none)
    real(kind=dp), intent(in) :: ext_centres(3, n_ext)
    !! Cartesian external centres (Angstrom)
    character(len=*), intent(in) :: seedname
    !! Seed name for output files
    integer, intent(in) :: stdout
    !! Output unit
    type(w90_error_type), allocatable, intent(out) :: error
    !! Error handler
    type(w90_comm_type), intent(in) :: comm
    !! Communicator

    real(kind=dp), allocatable :: alphas(:, :), betas(:, :, :)
    real(kind=dp), allocatable :: rho(:, :, :), group_rho(:, :, :)
    real(kind=dp), allocatable :: coeff(:), power(:)
    real(kind=dp) :: asc(3, 3), cross(3), det_sc, area, dist, bound
    real(kind=dp) :: rsum, dvol, det_cell, recip(3, 3), volume
    integer :: nsc(3), n_coeff, n_power, n, e, ierr
    logical :: do_own(n_wf)
    character(len=5) :: numstr
    character(len=256) :: fname

    ! --- Collective section (grid-independent): safe error synchronisation ---

    nsc(1) = (nxx_hi - nxx_lo + 1)/ngx
    nsc(2) = (nyy_hi - nyy_lo + 1)/ngy
    nsc(3) = (nzz_hi - nzz_lo + 1)/ngz

    ! Supercell lattice built from the grid bounds
    asc(1, :) = real(nsc(1), dp)*real_lattice(1, :)
    asc(2, :) = real(nsc(2), dp)*real_lattice(2, :)
    asc(3, :) = real(nsc(3), dp)*real_lattice(3, :)

    ! Signed cell volume via triple product a1 . (a2 x a3)
    cross(1) = asc(2, 2)*asc(3, 3) - asc(2, 3)*asc(3, 2)
    cross(2) = asc(2, 3)*asc(3, 1) - asc(2, 1)*asc(3, 3)
    cross(3) = asc(2, 1)*asc(3, 2) - asc(2, 2)*asc(3, 1)
    det_sc = abs(asc(1, 1)*cross(1) + asc(1, 2)*cross(2) + asc(1, 3)*cross(3))

    ! Inscribed-sphere radius = 0.5 * min over faces of |det|/area(face)
    bound = huge(1.0_dp)
    call cross_product_local(asc(2, :), asc(3, :), cross)
    area = sqrt(dot_product(cross, cross))
    dist = det_sc/area
    bound = min(bound, 0.5_dp*dist)
    call cross_product_local(asc(3, :), asc(1, :), cross)
    area = sqrt(dot_product(cross, cross))
    dist = det_sc/area
    bound = min(bound, 0.5_dp*dist)
    call cross_product_local(asc(1, :), asc(2, :), cross)
    area = sqrt(dot_product(cross, cross))
    dist = det_sc/area
    bound = min(bound, 0.5_dp*dist)

    if (r_cut <= 0.0_dp) then
      call set_error_input(error, 'decompose_main: r_cut must be positive', comm)
      return
    end if
    if (r_cut > bound) then
      write (fname, '(a,f0.6,a)') &
        'decompose_main: r_cut exceeds the inscribed-sphere radius of the &
        &Born-von-Karman supercell (', bound, ' Angstrom); reduce r_cut or &
        &increase the k-point mesh'
      call set_error_input(error, trim(fname), comm)
      return
    end if

    ! Radial basis (collective error handling)
    allocate (alphas(n_max, 0:l_max), betas(n_max, n_max, 0:l_max), stat=ierr)
    if (ierr /= 0) then
      call set_error_input(error, 'decompose_main: allocation of basis arrays failed', comm)
      return
    end if
    call decompose_radial_params(n_max, l_max, r_min, r_max, alphas, betas, error, comm)
    if (allocated(error)) return

    ! --- Root-only section: projection and I/O ---
    if (mpirank(comm) /= 0) return

    n_coeff = decompose_num_coeff(n_max, l_max)
    n_power = decompose_num_power(n_max, l_max)

    ! Volume element on the primitive-cell grid
    call utility_recip_lattice_base(real_lattice, recip, volume)
    det_cell = volume
    dvol = det_cell/real(ngx*ngy*ngz, dp)

    allocate (rho(nxx_lo:nxx_hi, nyy_lo:nyy_hi, nzz_lo:nzz_hi), stat=ierr)
    if (ierr /= 0) then
      call set_error_local(error, 'decompose_main: allocation of rho failed')
      return
    end if
    allocate (coeff(n_coeff), power(n_power), stat=ierr)
    if (ierr /= 0) then
      call set_error_local(error, 'decompose_main: allocation of coefficient arrays failed')
      return
    end if
    if (n_ext > 0) then
      allocate (group_rho(nxx_lo:nxx_hi, nyy_lo:nyy_hi, nzz_lo:nzz_hi), stat=ierr)
      if (ierr /= 0) then
        call set_error_local(error, 'decompose_main: allocation of group density failed')
        return
      end if
      group_rho = 0.0_dp
    end if

    write (stdout, '(/,1x,a)') 'Decomposing Wannier function densities'
    write (stdout, '(3x,a,i0,a,i0)') 'n_max = ', n_max, ',  l_max = ', l_max
    write (stdout, '(3x,a,f0.4,a,f0.4,a,f0.4,a)') &
      'r_min = ', r_min, ',  r_max = ', r_max, ',  r_cut = ', r_cut, ' (Angstrom)'
    write (stdout, '(3x,a,i0,a,i0)') 'Number of Wannier functions: ', n_wf, &
      ', decomposed about their own centres: ', n_out
    if (n_ext > 0) write (stdout, '(3x,a,i0)') &
      'Group-density external centres: ', n_ext

    if (any(out_slot < 1) .or. any(out_slot > n_wf)) then
      call set_error_local(error, 'decompose_main: out_slot index out of range')
      return
    end if
    do_own = .false.
    do n = 1, n_out
      do_own(out_slot(n)) = .true.
    end do

    do n = 1, n_wf
      ! Density and normalisation to unit integral over the supercell
      rho = real(wann_func(:, :, :, n)*conjg(wann_func(:, :, :, n)), dp)
      rsum = sum(rho)
      if (rsum <= 0.0_dp) then
        call set_error_local(error, 'decompose_main: non-positive density norm')
        return
      end if
      rho = rho/(rsum*dvol)
      if (n_ext > 0) group_rho = group_rho + rho

      ! Own-centre orbital channel (only for the selected subset)
      if (.not. do_own(n)) cycle

      call decompose_project(rho, ngx, ngy, ngz, nxx_lo, nxx_hi, nyy_lo, nyy_hi, &
                             nzz_lo, nzz_hi, real_lattice, centres(:, n), n_max, l_max, &
                             r_cut, alphas, betas, coeff)
      call decompose_power_orb(coeff, n_max, l_max, power)

      write (numstr, '(i5.5)') wf_index(n)
      fname = trim(seedname)//'_'//numstr//'.coeff'
      call decompose_write_coeff(fname, coeff, n_coeff, n_max, l_max, r_min, r_max, r_cut, error)
      if (allocated(error)) return
      fname = trim(seedname)//'_'//numstr//'.power'
      call decompose_write_power(fname, power, n_power, n_max, l_max, r_min, r_max, r_cut, error)
      if (allocated(error)) return
    end do

    ! Group-density channel about each external centre
    do e = 1, n_ext
      call decompose_project(group_rho, ngx, ngy, ngz, nxx_lo, nxx_hi, nyy_lo, nyy_hi, &
                             nzz_lo, nzz_hi, real_lattice, ext_centres(:, e), n_max, l_max, &
                             r_cut, alphas, betas, coeff)
      write (numstr, '(i5.5)') e
      fname = trim(seedname)//'_gc_'//numstr//'.coeff'
      call decompose_write_coeff(fname, coeff, n_coeff, n_max, l_max, r_min, r_max, r_cut, error)
      if (allocated(error)) return
    end do

    write (stdout, '(3x,a,/)') 'Decomposition complete'

  end subroutine decompose_main

  !================================================================!
  subroutine decompose_write_coeff(fname, coeff, n_coeff, n_max, l_max, &
                                   r_min, r_max, r_cut, error)
    !! Write a flat coefficient vector to an ASCII, self-describing file.

    use w90_error, only: w90_error_type

    implicit none

    character(len=*), intent(in) :: fname
    !! Output file name
    integer, intent(in) :: n_coeff
    !! Number of coefficients
    real(kind=dp), intent(in) :: coeff(n_coeff)
    !! Coefficients
    integer, intent(in) :: n_max, l_max
    !! Basis sizes
    real(kind=dp), intent(in) :: r_min, r_max, r_cut
    !! Basis radii (Angstrom)
    type(w90_error_type), allocatable, intent(out) :: error
    !! Error handler

    integer :: unit, ierr, i

    open (newunit=unit, file=trim(fname), status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
      call set_error_local(error, 'decompose_write_coeff: cannot open '//trim(fname))
      return
    end if

    write (unit, '(a)') '# Wannier density decomposition coefficients'
    write (unit, '(a,i0)') '# n_max = ', n_max
    write (unit, '(a,i0)') '# l_max = ', l_max
    write (unit, '(a,es24.15e3)') '# r_min = ', r_min
    write (unit, '(a,es24.15e3)') '# r_max = ', r_max
    write (unit, '(a,es24.15e3)') '# r_cut = ', r_cut
    write (unit, '(a)') '# ordering: outer n (0..n_max-1), then l (0..l_max), &
                        &then inner m (0..2l)'
    write (unit, '(a,i0)') '# number of coefficients = ', n_coeff
    do i = 1, n_coeff
      write (unit, '(es24.15e3)') coeff(i)
    end do
    close (unit)

  end subroutine decompose_write_coeff

  !================================================================!
  subroutine decompose_write_power(fname, power, n_power, n_max, l_max, &
                                   r_min, r_max, r_cut, error)
    !! Write an orbital-orbital power spectrum to an ASCII file.

    use w90_error, only: w90_error_type

    implicit none

    character(len=*), intent(in) :: fname
    !! Output file name
    integer, intent(in) :: n_power
    !! Number of power entries
    real(kind=dp), intent(in) :: power(n_power)
    !! Power spectrum
    integer, intent(in) :: n_max, l_max
    !! Basis sizes
    real(kind=dp), intent(in) :: r_min, r_max, r_cut
    !! Basis radii (Angstrom)
    type(w90_error_type), allocatable, intent(out) :: error
    !! Error handler

    integer :: unit, ierr, i

    open (newunit=unit, file=trim(fname), status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
      call set_error_local(error, 'decompose_write_power: cannot open '//trim(fname))
      return
    end if

    write (unit, '(a)') '# Wannier density decomposition power spectrum (orbital-orbital)'
    write (unit, '(a,i0)') '# n_max = ', n_max
    write (unit, '(a,i0)') '# l_max = ', l_max
    write (unit, '(a,es24.15e3)') '# r_min = ', r_min
    write (unit, '(a,es24.15e3)') '# r_max = ', r_max
    write (unit, '(a,es24.15e3)') '# r_cut = ', r_cut
    write (unit, '(a)') '# ordering: outer n1 (0..n_max-1), then n2 (n1..n_max-1), &
                        &then l (0..l_max); entry = sum_m c(n1,l,m) c(n2,l,m)'
    write (unit, '(a,i0)') '# number of entries = ', n_power
    do i = 1, n_power
      write (unit, '(es24.15e3)') power(i)
    end do
    close (unit)

  end subroutine decompose_write_power

  !================================================================!
  subroutine set_error_local(error, mesg)
    !! Set a fatal error object without collective synchronisation, for
    !! use inside the root-only section where the other ranks have already
    !! returned (a collective `set_error_*` would deadlock).

    use w90_error, only: w90_error_type, set_base_error, code_fatal

    implicit none

    type(w90_error_type), allocatable, intent(out) :: error
    !! Error handler
    character(len=*), intent(in) :: mesg
    !! Error message

    call set_base_error(error, mesg, code_fatal)

  end subroutine set_error_local

  !================================================================!
  subroutine cross_product_local(a, b, c)
    !! Cartesian cross product c = a x b.
    implicit none
    real(kind=dp), intent(in) :: a(3), b(3)
    real(kind=dp), intent(out) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end subroutine cross_product_local

end module w90_decompose
