"""
Parser function parse() to parse the <seedname>_geninerp.dat output file of postw90.x.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

def parse_elcond(fname):
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
        if len(pieces) != 8:
            raise ValueError("Wrong line length ({}, instead of 5 or 8); line content: {}".format(
                len(pieces)), l)
        retdict['mu'].append(float(pieces[0]))
        retdict['temp'].append(float(pieces[1]))
        retdict['elcond_xx'].append(float(pieces[2]))
        retdict['elcond_xy'].append(float(pieces[3]))
        retdict['elcond_yy'].append(float(pieces[4]))
        retdict['elcond_xz'].append(float(pieces[5]))
        retdict['elcond_yz'].append(float(pieces[6]))
        retdict['elcond_zz'].append(float(pieces[7]))

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
