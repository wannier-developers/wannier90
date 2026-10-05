"""
Parser function parse() to parse the <seedname>_band.kpt output file of wannier90.x.

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
            if lno == 0:
                retdict["nkp"] = [int(l)]
                continue
            pieces = l.split()

            if len(pieces) == 4: # three kpoint coordinates and weight
                retdict["kpath"].append([float(pieces[0]), float(pieces[1]), float(pieces[2])])
                retdict["weight"].append(float(pieces[3]))
            else:
                raise ValueError("Wrong line length ({}, instead of 4); line content: {}".format(
                    len(pieces), l))

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
