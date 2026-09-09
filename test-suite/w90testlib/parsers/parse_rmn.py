"""
Parser function parse() to parse the seedname_r.dat output file of Wannier90 (<i|r|j> matrix elements).
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

    # skip the first line, which contains only the date
    retdict['num_wann'] = [int(_) for _ in lines[1].split()]
    retdict['nrpts'] = [int(_) for _ in lines[2].split()]
    retdict['irvec_a'] = [float(line.split()[0]) for line in lines[3:]]
    retdict['irvec_b'] = [float(line.split()[1]) for line in lines[3:]]
    retdict['irvec_c'] = [float(line.split()[2]) for line in lines[3:]]
    retdict['index_i'] = [float(line.split()[3]) for line in lines[3:]]
    retdict['index_j'] = [float(line.split()[4]) for line in lines[3:]]
    retdict['real_x'] = [float(line.split()[5]) for line in lines[3:]]
    retdict['imag_x'] = [float(line.split()[6]) for line in lines[3:]]
    retdict['real_y'] = [float(line.split()[7]) for line in lines[3:]]
    retdict['imag_y'] = [float(line.split()[8]) for line in lines[3:]]
    retdict['real_z'] = [float(line.split()[9]) for line in lines[3:]]
    retdict['imag_z'] = [float(line.split()[10]) for line in lines[3:]]

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
