"""
Parser function parse() to parse the <seedname>_kubo*.dat output file of postw90.x.

"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

def parse(fname):
    """
    Open the file, parses it and return the values.
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for lno, l in enumerate(lines):
        
        if l.strip().startswith('#'):
            # Skip headers
            continue

        pieces = l.split()

        if len(pieces) == 0 :
            # skip blank line
            continue

        if len(pieces) == 11 : 
            retdict['efermi'].append(float(pieces[0]))
            retdict['omega'].append(float(pieces[1]))
            retdict['gyro_xx'].append(float(pieces[2]))
            retdict['gyro_yy'].append(float(pieces[3]))
            retdict['gyro_zz'].append(float(pieces[4]))
            retdict['gyro_xy'].append(float(pieces[5]))
            retdict['gyro_xz'].append(float(pieces[6]))
            retdict['gyro_yz'].append(float(pieces[7]))
            retdict['gyro_x'].append(float(pieces[8]))
            retdict['gyro_y'].append(float(pieces[9]))
            retdict['gyro_z'].append(float(pieces[10]))
        else:
            raise ValueError("Wrong line length ({}, instead of 3); line content: {}".format(
                len(pieces)), l)


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
