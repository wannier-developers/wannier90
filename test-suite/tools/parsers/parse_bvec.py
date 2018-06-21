"""
Parser function parse() to parse the .bvec output file of Wannier90 (bvec matrices).
"""
from __future__ import print_function, unicode_literals
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

    #read the values on the second line, that are the size of the matrix
    retdict['size'] = [int(_) for _ in lines[1].split()]

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
