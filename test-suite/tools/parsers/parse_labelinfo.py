"""
Parser function parse() to parse the _band.labelinfo.dat output file of Wannier90 (info on high-symmetry labels and k-points).
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

def parse(fname):
    """
    Open the file, parses it and return the values

    For now, I just check that the size of the file is correct, but
    I don't check the actual content
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    labels_dict = {}
    for line_idx, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        label, idx, xval, kx, ky, kz = line.split()
        idx = int(idx)
        xval = float(xval)
        kx = float(kx)
        ky = float(ky)
        kz = float(kz)

    #read the values on the second line, that are the size of the matrix
    retdict['labels'].append(label)
    retdict['indices'].append(idx)
    retdict['xvals'].append(xval)
    retdict['kx'].append(kx)
    retdict['ky'].append(ky)
    retdict['kz'].append(kz)

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
