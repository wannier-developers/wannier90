"""
Parser function parse() to parse the <seedname>_band.dat output file of wannier90.x.

"""
import inspect
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

            if len(pieces) == 0: # blank line
                continue

            if len(pieces) == 2: # kpath eigval
                retdict["kpath"].append(float(pieces[0]))
                retdict["eigval"].append(float(pieces[1]))
            else:
                raise ValueError("Wrong line length ({}, instead of 2); line content: {}".format(
                    len(pieces), l))

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict



