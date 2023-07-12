"""
Parser function parse() to parse the <seedname>-shg*.dat output file of postw90.x.

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

        if len(pieces) == 5 :
            retdict['energy'].append(float(pieces[0]))
            retdict['re_chi'].append(float(pieces[1]))
            retdict['im_chi'].append(float(pieces[2]))
            retdict['re_sigma'].append(float(pieces[3]))
            retdict['im_sigma'].append(float(pieces[4]))
        else:
            raise ValueError("Wrong line length ({}, instead of 3); line content: {}".format(
                len(pieces)), l)


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
