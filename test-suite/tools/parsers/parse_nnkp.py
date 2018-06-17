"""
Parser function parse() to parse the .nnkp output file of Wannier90 -pp.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

# Match the lines in the nnkpts block
# Groups: 
# 0: idx of first kpt
# 1: idx of second kpt
# 2,3,4: g1,g2,g3 (reciprocal lattice vector needed to bring them together)
nnkpts_re = re.compile("^\s*\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$")
# Match the lines in the kpoints block
# Groups
# 0: k1
# 0: k2
# 0: k3
kpts_re = re.compile("^\s*([\d\.-]+)\s+([\d\.-]+)\s+([\d\.-]+)\s*$")

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

    # The file is case insensitive, change everything to lowercase
    lines = [l.lower() for l in lines]

    for lno, l in enumerate(lines):
        
        ###############################################################
        # Look for 'begin nnkpts'
        # Start from the second line after (the immediate next one is just 
        # a count) and stop at 'end nnkpts'
        if "begin nnkpts" in l:
            for l2 in lines[lno+2:]: # Skip 1 line
                match = nnkpts_re.search(l2)
                if not match or 'end nnkpts' in l2:
                    break
                k1, k2, g1, g2, g3 = match.groups() 
                retdict["nnkpts_k1"].append(int(k1))
                retdict["nnkpts_k2"].append(int(k2))
                retdict["nnkpts_g1"].append(int(g1))
                retdict["nnkpts_g2"].append(int(g2))
                retdict["nnkpts_g3"].append(int(g3))
            continue

        ###############################################################
        # Look for 'begin kpoints'
        # Start from the second line after (the immediate next one is just 
        # a count) and stop at 'end kpoints'
        if "begin kpoints" in l:
            for l2 in lines[lno+2:]: # Skip 1 line
                match = kpts_re.search(l2)
                if not match or 'end kpoints' in l2:
                    break
                k1, k2, k3 = match.groups() 
                retdict["kpoints_k1"].append(float(k1))
                retdict["kpoints_k2"].append(float(k2))
                retdict["kpoints_k3"].append(float(k3))
            continue

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
