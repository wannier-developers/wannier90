"""
Parser function parse() to parse the <seedname>_geninerp.dat output file of postw90.x.
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
        if len(pieces) == 5 or len(pieces) == 8: # No derivatives
            retdict['bandidx'].append(int(pieces[0]))
            retdict['bandkptx'].append(float(pieces[1]))
            retdict['bandkpty'].append(float(pieces[2]))
            retdict['bandkptz'].append(float(pieces[3]))
            retdict['bandenergy'].append(float(pieces[4]))
            if len(pieces) == 8:
                retdict['bandderivx'].append(float(pieces[5]))
                retdict['bandderivy'].append(float(pieces[6]))
                retdict['bandderivz'].append(float(pieces[7]))
        else:
            raise ValueError("Wrong line length ({}, instead of 5 or 8); line content: {}".format(
                len(pieces)), l)


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
