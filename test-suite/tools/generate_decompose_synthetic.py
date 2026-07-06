#!/usr/bin/env python3
"""Generate the synthetic physics-pinning inputs for the wannier_decompose
test-suite tests

  testw90_decompose_sphere : a single "Wannier function" that is an isotropic
      Gaussian u(r) = N exp(-|r - c|^2 / (2 sigma^2)) about a deliberately
      non-grid-aligned centre c, in a deliberately IRREGULAR TRICLINIC cell
      (unequal edges, all angles well off 90 degrees). Its density is pure
      l=0, so every l >= 1 coefficient of the decomposition is ~0 (up to a
      tiny, deterministic grid-discretisation residue). The triclinic cell is
      the point of this test: it exercises the non-orthorhombic geometry path
      of decompose_project (fractional<->Cartesian transforms through the full
      lattice matrix, minimum-image over a skewed Born-von-Karman cell, general
      volume element). A cubic cell cannot expose lattice-convention errors
      (e.g. a transposed real_lattice) because its lattice matrix is diagonal;
      the triclinic cell can.

  testw90_decompose_l2 : a single d_z2-like "Wannier function"
      u(r) = N (3 z'^2 - r'^2) exp(-r'^2 / (2 sigma^2)), r' = r - c, in a cubic
      cell. Its density contains ONLY (l, m) = (0,0), (2,0) and (4,0)
      components: all odd-l and all m /= 0 coefficients are ~0, and the l=2 /
      l=4 values are analytically known (see check_decompose_synthetic.py).
      This test stays cubic on purpose (it pins the l-selectivity of the
      angular basis, which the sphere test cannot); the triclinic sphere pins
      the geometry.

Both systems are gamma-only (mp_grid = 1 1 1), num_bands = num_wann = 1,
with an identity .amn and num_iter = 0 so that the U matrix stays 1 and the
"Wannier function" is exactly the synthetic u(r) written to UNK00001.1.
The .mmn overlaps M(b) = <u| e^{-i b.r} |u> are computed numerically on the
grid for exactly the b-vectors wannier90 selects (obtained by running
`wannier90.x -pp`; for a triclinic cell this may be MORE than the cubic six,
in several shells). Because |u|^2 is even about c, M(b) = e^{-i b.c} |M(b)|
exactly, hence Im ln M(b) = -b.c and the WF centre wannier90 reports equals
the true Gaussian centre (the least-squares recovery is printed as a check).

Everything is deterministic (no random numbers, no timestamps): regenerating
the inputs reproduces them bit-for-bit. The generated files are nevertheless
committed, so the tests are self-contained and this script only needs to be
re-run if the synthetic systems are changed.

Usage:
    python3 generate_decompose_synthetic.py --w90 /path/to/wannier90.x \
            [--testdir /path/to/test-suite/tests]
"""

import argparse
import os
import subprocess

import numpy as np


# ----------------------------------------------------------------------
# Lattice helpers
# ----------------------------------------------------------------------
def lattice_from_params(a, b, c, alpha, beta, gamma):
    """Row-stacked lattice vectors (Angstrom) from cell lengths (Angstrom)
    and angles (degrees), in the wannier90 real_lattice convention
    real_lattice(i, :) = a_i. Standard crystallographic setting: a_1 along x,
    a_2 in the xy-plane."""
    al, be, ga = np.deg2rad([alpha, beta, gamma])
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([b * np.cos(ga), b * np.sin(ga), 0.0])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = c * np.sqrt(1.0 - (cx / c) ** 2 - (cy / c) ** 2)
    return np.array([a1, a2, [cx, cy, cz]])


def reciprocal(A):
    """Reciprocal lattice (rows) b_i with a_i . b_j = 2 pi delta_ij, via cross
    products (exact zeros for a diagonal cell, so the cubic l2 system is
    reproduced bit-for-bit)."""
    vol = np.dot(A[0], np.cross(A[1], A[2]))
    b1 = 2.0 * np.pi * np.cross(A[1], A[2]) / vol
    b2 = 2.0 * np.pi * np.cross(A[2], A[0]) / vol
    b3 = 2.0 * np.pi * np.cross(A[0], A[1]) / vol
    return np.array([b1, b2, b3])


def inscribed_bound(A):
    """r_cut upper bound = 0.5 * min over faces of |det A| / area(face)
    (the decompose_main validation, for mp_grid = 1)."""
    vol = abs(np.dot(A[0], np.cross(A[1], A[2])))
    bounds = []
    for (j, k) in ((1, 2), (2, 0), (0, 1)):
        area = np.linalg.norm(np.cross(A[j], A[k]))
        bounds.append(0.5 * vol / area)
    return min(bounds)


# ----------------------------------------------------------------------
# System definitions (shared)
# ----------------------------------------------------------------------
# Sphere: irregular triclinic cell a=7,b=8,c=9 Ang; alpha=70,beta=80,gamma=75 deg.
# Inscribed bound (0.5*min interplanar spacing) = 3.3658 Ang, so r_cut = 3.0 is
# comfortably inside; the Gaussian density (sigma=0.8) is negligible
# (~8e-7) at r = r_cut, hence fully contained. Grid 28^3 -> worst-direction
# spacing (min interplanar / 28) = 0.24 Ang, at least as fine as the old
# cubic test (0.33 Ang). Centre placed off-grid at fractional
# (0.4123, 0.5237, 0.4871) (none a multiple of 1/28).
_SPHERE_CELL = lattice_from_params(7.0, 8.0, 9.0, 70.0, 80.0, 75.0)
_SPHERE_CENTRE = tuple((np.array([0.4123, 0.5237, 0.4871]) @ _SPHERE_CELL).tolist())

# Cubic l2 cell (unchanged; kept bit-identical to the committed inputs)
_L2_CELL = 8.0 * np.eye(3)
_L2_CELL_LINES = [
    ' 8.0  0.0  0.0',
    '  0.0 8.0  0.0',
    '  0.0  0.0 8.0',
]

SYSTEMS = {
    'testw90_decompose_sphere': {
        'seed': 'sphere',
        'kind': 'sphere',
        'cell': _SPHERE_CELL,
        'ng': (28, 28, 28),
        'sigma': 0.8,                          # Angstrom
        'centre': _SPHERE_CENTRE,              # Angstrom, off-grid
        'n_max': 3,
        'l_max': 2,
        'r_min': 2.5,
        'r_max': 4.0,
        'r_cut': 3.0,
        'comment': [
            '! Synthetic physics-pinning test for wannier_decompose.',
            '! Single "WF" = isotropic Gaussian about a non-grid-aligned centre',
            '! in an irregular TRICLINIC cell (a=7,b=8,c=9 Ang; angles',
            '! 70,80,75 deg): the density is pure l=0, so all l>=1 coefficients',
            '! must vanish (up to a tiny grid-discretisation residue). The',
            '! skewed cell exercises the non-orthorhombic geometry path',
            '! (frac<->cart, minimum-image, general volume element). Identity',
            '! amn + num_iter=0 keep U=1, so the WF is exactly the UNK orbital.',
            '! Inputs generated by tools/generate_decompose_synthetic.py.',
        ],
    },
    'testw90_decompose_l2': {
        'seed': 'dz2',
        'kind': 'dz2',
        'cell': _L2_CELL,
        'cell_lines': _L2_CELL_LINES,
        'ng': (24, 24, 24),
        'sigma': 0.8,
        'centre': (3.3579, 3.1259, 2.9873),
        'n_max': 2,
        'l_max': 4,
        'r_min': 2.5,
        'r_max': 4.0,
        'r_cut': 3.8,
        'comment': [
            '! Synthetic physics-pinning test for wannier_decompose.',
            '! Single "WF" = d_z2-like orbital u ~ (3z\'^2 - r\'^2) exp(-r\'^2/2s^2)',
            '! about a non-grid-aligned centre: |u|^2 contains only',
            '! (l,m) = (0,0), (2,0), (4,0), so all odd-l and m/=0 coefficients',
            '! must vanish and the (2,0)/(4,0) values are analytically known.',
            '! Inputs generated by tools/generate_decompose_synthetic.py.',
        ],
    },
}

UNK_FMT = '{0: .10E}  0.0'   # real part, imaginary part is exactly zero


def grid_coords(cell, ng):
    """Cartesian coordinates of the 1-based grid point (i) at fractional
    coordinate (i-1)/ng along each axis (the wannier90 UNK /
    plot_build_wannier_grid convention). Returns X, Y, Z arrays indexed
    [ix, iy, iz]."""
    A = np.asarray(cell)
    fx = np.arange(ng[0]) / ng[0]
    fy = np.arange(ng[1]) / ng[1]
    fz = np.arange(ng[2]) / ng[2]
    Fx, Fy, Fz = np.meshgrid(fx, fy, fz, indexing='ij')
    X = Fx * A[0, 0] + Fy * A[1, 0] + Fz * A[2, 0]
    Y = Fx * A[0, 1] + Fy * A[1, 1] + Fz * A[2, 1]
    Z = Fx * A[0, 2] + Fy * A[1, 2] + Fz * A[2, 2]
    return X, Y, Z


def build_u(kind, sigma, centre, cell, ng):
    """Synthetic orbital on the grid, periodised over 3^3 lattice images."""
    A = np.asarray(cell)
    X, Y, Z = grid_coords(cell, ng)
    u = np.zeros_like(X)
    for Rx in (-1, 0, 1):
        for Ry in (-1, 0, 1):
            for Rz in (-1, 0, 1):
                sx = Rx * A[0, 0] + Ry * A[1, 0] + Rz * A[2, 0]
                sy = Rx * A[0, 1] + Ry * A[1, 1] + Rz * A[2, 1]
                sz = Rx * A[0, 2] + Ry * A[1, 2] + Rz * A[2, 2]
                dx = X - centre[0] - sx
                dy = Y - centre[1] - sy
                dz = Z - centre[2] - sz
                r2 = dx * dx + dy * dy + dz * dz
                g = np.exp(-r2 / (2.0 * sigma ** 2))
                if kind == 'sphere':
                    u += g
                elif kind == 'dz2':
                    u += (3.0 * dz * dz - r2) * g
                else:
                    raise ValueError(kind)
    # Normalize sum |u|^2 dV = 1 (cosmetic: wannier90 renormalizes anyway).
    # dV = |det(cell)| / (ngx*ngy*ngz), the general (triclinic) volume element.
    dv = abs(np.dot(A[0], np.cross(A[1], A[2]))) / (ng[0] * ng[1] * ng[2])
    u /= np.sqrt(np.sum(u * u) * dv)
    return u


def roundtrip_unk(u, ng):
    """Round u through the exact text written to UNK00001.1, so that all
    subsequent numbers (mmn overlaps, analytic checks) are computed from
    the values wannier90 actually reads."""
    flat = u.transpose(2, 1, 0).ravel()          # x contiguous
    lines = [UNK_FMT.format(v) for v in flat]
    rounded = np.array([float(l.split()[0]) for l in lines])
    u_r = rounded.reshape(ng[2], ng[1], ng[0]).transpose(2, 1, 0)
    return u_r, lines


def write_unk(dirname, lines, ng):
    with open(os.path.join(dirname, 'UNK00001.1'), 'w') as f:
        f.write('  {0}  {1}  {2}    1    1\n'.format(ng[0], ng[1], ng[2]))
        f.write('\n'.join(lines))
        f.write('\n')


def cell_block_lines(sysdef):
    """The three unit_cell_cart rows for the .win. Uses the verbatim
    'cell_lines' if the system provides them (so the cubic l2 win stays
    bit-identical); otherwise formats the numeric cell."""
    if 'cell_lines' in sysdef:
        return list(sysdef['cell_lines'])
    return ['  {0: .10f}  {1: .10f}  {2: .10f}'.format(*row)
            for row in np.asarray(sysdef['cell'])]


def write_win(dirname, sysdef):
    seed = sysdef['seed']
    c = sysdef['centre']
    lines = list(sysdef['comment'])
    lines += [
        '',
        'num_wann    = 1',
        'num_iter    = 0',
        '',
        'begin unit_cell_cart',
        'ang',
    ]
    lines += cell_block_lines(sysdef)
    lines += [
        'end unit_cell_cart',
        '',
        '! Dummy atom marking the true orbital centre',
        'begin atoms_cart',
        'ang',
        'H  {0:.4f}  {1:.4f}  {2:.4f}'.format(*c),
        'end atoms_cart',
        '',
        'begin projections',
        'random',
        'end projections',
        '',
        'mp_grid : 1 1 1',
        '',
        'begin kpoints',
        '0.0 0.0 0.0',
        'end kpoints',
        '',
        'wvfn_formatted = .true.',
        '',
        'wannier_decompose = true',
        'decompose_n_max   = {0}'.format(sysdef['n_max']),
        'decompose_l_max   = {0}'.format(sysdef['l_max']),
        'decompose_r_min   = {0}'.format(sysdef['r_min']),
        'decompose_r_max   = {0}'.format(sysdef['r_max']),
        'decompose_r_cut   = {0}'.format(sysdef['r_cut']),
    ]
    with open(os.path.join(dirname, seed + '.win'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


def read_nnkp_bvectors(dirname, seed):
    """Parse the nnkpts block of seed.nnkp: records (ik, ikp, G)."""
    fname = os.path.join(dirname, seed + '.nnkp')
    with open(fname) as f:
        lines = [l.strip() for l in f]
    i0 = lines.index('begin nnkpts')
    nntot = int(lines[i0 + 1])
    recs = []
    for l in lines[i0 + 2:i0 + 2 + nntot]:
        w = l.split()
        recs.append((int(w[0]), int(w[1]), int(w[2]), int(w[3]), int(w[4])))
    return recs


def bvector(rec, recip):
    """Cartesian b-vector for an nnkp record: b = G . B (B rows = recip)."""
    g = np.array(rec[2:5], dtype=float)
    return g @ recip


def overlap(rho, norm, X, Y, Z, b):
    """M(b) = sum_r rho(r) e^{-i b.r} dV / sum_r rho(r) dV."""
    phase = np.exp(-1j * (b[0] * X + b[1] * Y + b[2] * Z))
    return (rho * phase).sum() / norm


def write_mmn(dirname, seed, u, recs, X, Y, Z, recip):
    """M(b) for each nnkp record on the grid (dV and normalization cancel)."""
    rho = u * u
    norm = rho.sum()
    header = 'Synthetic gamma-only overlaps (generate_decompose_synthetic.py)'
    out = [header, '    1    1    {0}'.format(len(recs))]
    for rec in recs:
        m = overlap(rho, norm, X, Y, Z, bvector(rec, recip))
        out.append('    {0}    {1}    {2}    {3}    {4}'.format(*rec))
        out.append('  {0: .12f}  {1: .12f}'.format(m.real, m.imag))
    with open(os.path.join(dirname, seed + '.mmn'), 'w') as f:
        f.write('\n'.join(out) + '\n')


def write_amn(dirname, seed):
    with open(os.path.join(dirname, seed + '.amn'), 'w') as f:
        f.write('Identity projection (generate_decompose_synthetic.py)\n')
        f.write('           1           1           1\n')
        f.write('    1    1    1    1.000000000000    0.000000000000\n')


def report_bvectors_and_centre(sysdef, u, recs, X, Y, Z, recip):
    """Print the b-vector shells and recover the WF centre directly from the
    phases Im ln M(b) = -b.c (least-squares over all listed b), as a sanity
    check that is valid for a general (triclinic) cell."""
    rho = u * u
    norm = rho.sum()
    bmat = np.array([bvector(rec, recip) for rec in recs])
    bmag = np.linalg.norm(bmat, axis=1)
    shells = sorted(set(np.round(bmag, 6)))
    print('  b-vectors: {0} in {1} shell(s) |b| = {2}'.format(
        len(recs), len(shells), ['{0:.4f}'.format(s) for s in shells]))
    phases = np.array([np.angle(overlap(rho, norm, X, Y, Z, b)) for b in bmat])
    # b.c = -Im ln M(b); solve the overdetermined system for c. For a skewed
    # cell the phases b.c can exceed pi and wrap, so the recovered centre (like
    # wannier90's own) may sit a lattice vector away from c; compare modulo the
    # lattice (a WF centre is only defined up to a lattice translation).
    centre, *_ = np.linalg.lstsq(bmat, -phases, rcond=None)
    true_centre = np.array(sysdef['centre'])
    frac = (centre - true_centre) @ np.linalg.inv(np.asarray(sysdef['cell']))
    frac -= np.round(frac)
    dc = np.abs(frac @ np.asarray(sysdef['cell'])).max()
    print('  true centre        :', true_centre)
    print('  centre from M(b)   :', centre)
    print('  max|diff| (min-img): {0:.3e} Angstrom'.format(dc))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--w90', required=True, help='path to a wannier90.x binary '
                   '(used only for -pp, to obtain the b-vector list)')
    p.add_argument('--testdir', default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'tests'))
    args = p.parse_args()
    w90 = os.path.abspath(args.w90)

    for dirname, sysdef in SYSTEMS.items():
        seed = sysdef['seed']
        cell = np.asarray(sysdef['cell'])
        ng = sysdef['ng']
        recip = reciprocal(cell)
        d = os.path.abspath(os.path.join(args.testdir, dirname))
        os.makedirs(d, exist_ok=True)
        print('=== {0} (seed {1}) ==='.format(dirname, seed))
        print('  inscribed r_cut bound = {0:.6f} Ang (r_cut = {1})'.format(
            inscribed_bound(cell), sysdef['r_cut']))

        write_win(d, sysdef)

        # b-vectors exactly as wannier90 will choose them
        subprocess.run([w90, '-pp', seed], cwd=d, check=True,
                       stdout=subprocess.DEVNULL)
        recs = read_nnkp_bvectors(d, seed)
        os.remove(os.path.join(d, seed + '.nnkp'))
        for junk in (seed + '.wout',):
            path = os.path.join(d, junk)
            if os.path.exists(path):
                os.remove(path)

        u = build_u(sysdef['kind'], sysdef['sigma'], sysdef['centre'], cell, ng)
        u, lines = roundtrip_unk(u, ng)
        write_unk(d, lines, ng)
        X, Y, Z = grid_coords(cell, ng)
        write_mmn(d, seed, u, recs, X, Y, Z, recip)
        write_amn(d, seed)

        report_bvectors_and_centre(sysdef, u, recs, X, Y, Z, recip)


if __name__ == '__main__':
    main()
