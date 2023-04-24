"""
Parser function parse() to parse the <seedname>-kslice-shc.dat output file of postw90.x.

"""
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

        pieces = l.split()

        if len(pieces) == 0:
            # skip blank line
            continue
        elif len(pieces) == 1:
            retdict['shc'].append(float(pieces[0]))
        else:
            raise ValueError("Wrong line length ({}, instead of 1); line content: {}".format(
                len(pieces), l))

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
