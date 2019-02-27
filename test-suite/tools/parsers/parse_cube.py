"""
Parser function parse() to parse the <seedname>_00001.cube output file of Wannier90 (Gaussian cube format).
"""
from __future__ import print_function, unicode_literals

import inspect
import re
from collections import defaultdict

from . import show_output

def parse(fname):
    """
    Open the file, parses it and return the values

    For now, I just check that the origins of the small grid are the same
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    #read the values on the second line, that is the origin
    retdict['origin'] = [float(_) for _ in lines[2].split()]

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
