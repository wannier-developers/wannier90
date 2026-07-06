#!/usr/bin/env python3
"""Sanity checker for the synthetic wannier_decompose physics-pinning tests
(testw90_decompose_sphere / testw90_decompose_l2). See
generate_decompose_synthetic.py for the systems.

Given a directory in which wannier90.x has been run, this script

  1. parses <seed>.wout and compares the WF centre wannier90 reports with
     the true (analytically chosen) orbital centre;
  2. parses <seed>_00001.coeff (flat ordering: outer n, then l, then inner
     m = 0..2l) and reports max |c| per l and the l-purity ratios:
       - sphere: max |c(l>=1)| / max |c(l=0)|   (must be ~0; the density is
         pure l=0)
       - dz2   : max |c(odd l)| / max |c(l=0)| and max |c(even l, m/=0)| /
         max |c(l=0)|  (must be ~0), while c(l=2,m=0) must be clearly nonzero;
  3. compares every symmetry-allowed coefficient against its analytic
     continuum value, computed by 1D radial quadrature with the same
     Gaussian-radial basis (alphas/betas replicated from
     decompose_radial_params, analytic overlap + Loewdin).

Exits nonzero if any check fails.

Usage:
    python3 check_decompose_synthetic.py <rundir> sphere|dz2
"""

import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as gamma_fn

RATIO_TOL = 1.0e-4      # required l-purity
CENTRE_TOL = 1.0e-3     # Angstrom, wout centre vs true centre
ANALYTIC_RTOL = 1.0e-3  # continuum vs grid quadrature (discretisation)

# Must match generate_decompose_synthetic.py
SYSTEMS = {
    'sphere': {'sigma': 0.8, 'centre': (3.2117, 2.8331, 3.4569)},
    'dz2': {'sigma': 0.8, 'centre': (3.3579, 3.1259, 2.9873)},
}


def radial_params(n_max, l_max, r_min, r_max):
    """Replicates decompose_radial_params (analytic overlap + Loewdin)."""
    thr = 1.0e-3
    if n_max == 1:
        r_thr = np.array([r_min])
    else:
        r_thr = np.linspace(r_min, r_max, n_max)
    alphas = np.empty((n_max, l_max + 1))
    betas = np.empty((n_max, n_max, l_max + 1))
    for l in range(l_max + 1):
        alphas[:, l] = -(np.log(thr) - l * np.log(r_thr)) / r_thr ** 2
        a = alphas[:, l]
        smat = 0.5 * gamma_fn(l + 1.5) * (a[:, None] + a[None, :]) ** -(l + 1.5)
        eval_, evec = np.linalg.eigh(smat)
        betas[:, :, l] = (evec / np.sqrt(eval_)) @ evec.T
    return alphas, betas


def parse_coeff(fname):
    header = {}
    vals = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                if '=' in line:
                    k = line[1:].split('=')[0].strip()
                    v = line.split('=')[1].strip()
                    header[k] = v
            elif line.strip():
                vals.append(float(line))
    n_max = int(header['n_max'])
    l_max = int(header['l_max'])
    r_min = float(header['r_min'])
    r_max = float(header['r_max'])
    r_cut = float(header['r_cut'])
    c = np.array(vals).reshape(n_max, (l_max + 1) ** 2)
    return c, n_max, l_max, r_min, r_max, r_cut


def parse_wout_centre(fname):
    with open(fname) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if 'Final State' in line:
            w = lines[i + 1]
            assert 'WF centre and spread' in w
            xyz = w.split('(')[1].split(')')[0].split(',')
            return np.array([float(x) for x in xyz])
    raise RuntimeError('Final State block not found in ' + fname)


def per_l_max(c, n_max, l_max):
    """max |c| over n and m for each l, plus max |c| for m != 0 per l."""
    out = np.zeros(l_max + 1)
    out_moff = np.zeros(l_max + 1)
    for l in range(l_max + 1):
        block = c[:, l ** 2:(l + 1) ** 2]           # m index 0..2l, m = idx - l
        out[l] = np.abs(block).max()
        moff = np.delete(block, l, axis=1)
        out_moff[l] = np.abs(moff).max() if moff.size else 0.0
    return out, out_moff


def analytic_coeffs(kind, sigma, n_max, l_max, r_min, r_max, r_cut):
    """Continuum c_{n,l,0} for the symmetry-allowed channels."""
    alphas, betas = radial_params(n_max, l_max, r_min, r_max)

    def gnl(r, n, l):
        return sum(betas[np_, n, l] * r ** l * np.exp(-alphas[np_, l] * r * r)
                   for np_ in range(n_max))

    out = {}
    if kind == 'sphere':
        # rho = C exp(-r^2/sigma^2), C = 1/(pi^{3/2} sigma^3), only (l,m)=(0,0):
        # c_n00 = (1/sqrt(4pi)) * 4pi * int rho g_n0 r^2 dr
        C = np.pi ** -1.5 * sigma ** -3
        for n in range(n_max):
            val = quad(lambda r: C * np.exp(-(r / sigma) ** 2) * gnl(r, n, 0) * r * r,
                       0.0, r_cut)[0]
            # c_n00 = Y00 * int rho g_n0 dV = (1/sqrt(4pi)) * 4pi * val
            out[(n, 0)] = 4.0 * np.pi * val / np.sqrt(4.0 * np.pi)
    elif kind == 'dz2':
        # rho = A^2 r^4 (3cos^2(th)-1)^2 exp(-r^2/sigma^2),
        # (3t^2-1)^2 = a0 P0 + a2 P2 + a4 P4, a = (4/5, 8/7, 72/35)
        # c_nl0 = A^2 sqrt(4pi/(2l+1)) a_l int r^6 exp(-r^2/s^2) g_nl dr
        A2 = 1.0 / (3.0 * np.pi ** 1.5 * sigma ** 7)
        a_l = {0: 4.0 / 5.0, 2: 8.0 / 7.0, 4: 72.0 / 35.0}
        for l in (0, 2, 4):
            if l > l_max:
                continue
            for n in range(n_max):
                val = quad(lambda r: r ** 6 * np.exp(-(r / sigma) ** 2) * gnl(r, n, l),
                           0.0, r_cut)[0]
                out[(n, l)] = A2 * np.sqrt(4.0 * np.pi / (2 * l + 1)) * a_l[l] * val
    return out


def main():
    rundir, kind = sys.argv[1], sys.argv[2]
    sysdef = SYSTEMS[kind]
    seed = kind if kind == 'sphere' else 'dz2'
    ok = True

    # --- centre ---
    centre = parse_wout_centre(os.path.join(rundir, seed + '.wout'))
    true_centre = np.array(sysdef['centre'])
    dc = np.abs(centre - true_centre).max()
    print('W90 centre   : {0}'.format(centre))
    print('true centre  : {0}'.format(true_centre))
    print('max |diff|   : {0:.3e} Angstrom (tol {1:.0e})'.format(dc, CENTRE_TOL))
    if dc > CENTRE_TOL:
        ok = False
        print('  ** CENTRE CHECK FAILED **')

    # --- coefficients ---
    c, n_max, l_max, r_min, r_max, r_cut = parse_coeff(
        os.path.join(rundir, seed + '_00001.coeff'))
    lmax_abs, lmax_moff = per_l_max(c, n_max, l_max)
    print('max |c| per l        :', ' '.join('{0:.3e}'.format(v) for v in lmax_abs))
    print('max |c| per l (m!=0) :', ' '.join('{0:.3e}'.format(v) for v in lmax_moff))

    if kind == 'sphere':
        ratio = lmax_abs[1:].max() / lmax_abs[0]
        print('l-purity ratio max|c(l>=1)|/max|c(l=0)| = {0:.3e} (tol {1:.0e})'
              .format(ratio, RATIO_TOL))
        if ratio > RATIO_TOL:
            ok = False
            print('  ** L-PURITY CHECK FAILED **')
    else:
        odd = max(lmax_abs[l] for l in range(1, l_max + 1, 2))
        moff = max(lmax_moff[l] for l in range(0, l_max + 1, 2))
        r_odd = odd / lmax_abs[0]
        r_moff = moff / lmax_abs[0]
        r_l2 = lmax_abs[2] / lmax_abs[0]
        print('odd-l  ratio max|c(odd l)|/max|c(l=0)|      = {0:.3e} (tol {1:.0e})'
              .format(r_odd, RATIO_TOL))
        print('m!=0   ratio max|c(even l,m!=0)|/max|c(l=0)| = {0:.3e} (tol {1:.0e})'
              .format(r_moff, RATIO_TOL))
        print('signal ratio max|c(l=2)|/max|c(l=0)|         = {0:.3e} (must be O(0.1-1))'
              .format(r_l2))
        if r_odd > RATIO_TOL or r_moff > RATIO_TOL:
            ok = False
            print('  ** L-PURITY CHECK FAILED **')
        if r_l2 < 1.0e-2:
            ok = False
            print('  ** L=2 SIGNAL CHECK FAILED (no d-character?) **')

    # --- analytic continuum values ---
    ana = analytic_coeffs(kind, sysdef['sigma'], n_max, l_max, r_min, r_max, r_cut)
    print('symmetry-allowed coefficients vs analytic continuum values:')
    for (n, l), val in sorted(ana.items()):
        got = c[n, l ** 2 + l]   # m = 0
        rel = abs(got - val) / max(abs(val), 1e-30)
        flag = ''
        if rel > ANALYTIC_RTOL:
            ok = False
            flag = '   ** ANALYTIC CHECK FAILED **'
        print('  c(n={0},l={1},m=0): w90 = {2: .9e}  analytic = {3: .9e}  '
              'rel.diff = {4:.2e}{5}'.format(n, l, got, val, rel, flag))

    print('RESULT:', 'OK' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
