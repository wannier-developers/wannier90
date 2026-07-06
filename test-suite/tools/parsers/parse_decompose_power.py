"""
Parser function parse() to parse the <seedname>_NNNNN.power output file of
wannier90.x written by the `wannier_decompose` feature: the (orbital-orbital)
power spectrum built from the density-projection coefficients, together with
the self-describing header.
"""
from __future__ import print_function, unicode_literals

import inspect
from collections import defaultdict

from . import show_output


def parse(fname):
    """
    Open the file, parse the header (n_max, l_max, r_min, r_max, r_cut,
    number of entries) and the flat list of power-spectrum entries, and
    return the values.
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('# n_max'):
            retdict['n_max'] = [int(line.split('=')[1])]
        elif line.startswith('# l_max'):
            retdict['l_max'] = [int(line.split('=')[1])]
        elif line.startswith('# r_min'):
            retdict['r_min'] = [float(line.split('=')[1])]
        elif line.startswith('# r_max'):
            retdict['r_max'] = [float(line.split('=')[1])]
        elif line.startswith('# r_cut'):
            retdict['r_cut'] = [float(line.split('=')[1])]
        elif line.startswith('# number of entries'):
            retdict['nentries'] = [int(line.split('=')[1])]

    retdict['power'] = [float(line) for line in lines
                         if line.strip() and not line.startswith('#')]

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
