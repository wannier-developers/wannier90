"""
Parser function parse() to parse the <seedname>_tdf2.dat output file of postw90.x.
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
        print( "[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for lno, l in enumerate(lines):
        
        if l.strip().startswith('#'):
            # Skip headers
            continue

        pieces = l.split()
        if len(pieces) == 10 : # No derivatives
            retdict['energy'].append(float(pieces[0]))
            retdict['tdf2_xxz'].append(float(pieces[1]))
            retdict['tdf2_xyz'].append(float(pieces[2]))
            retdict['tdf2_yyz'].append(float(pieces[3]))
            retdict['tdf2_xzz'].append(float(pieces[4]))
            retdict['tdf2_yzz'].append(float(pieces[5]))
            retdict['tdf2_zzz'].append(float(pieces[6]))
            retdict['tdf2_yxz'].append(float(pieces[7]))
            retdict['tdf2_zxz'].append(float(pieces[8]))
            retdict['tdf2_zyz'].append(float(pieces[9]))
        else:
            raise ValueError("Wrong line length ({}, instead of 10 ); line content: {}".format(
                len(pieces)), l)


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
