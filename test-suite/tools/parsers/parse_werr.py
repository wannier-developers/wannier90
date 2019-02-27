"""
Parser function parse() to parse the .werr output file of Wannier90.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

def parse(fname):
    """
    Open the file, parses it and return the values
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for lno, l in enumerate(lines):

        ###############################################################
        # Look for 'Exiting...', then take the text of the next line
        if "Exiting..." in l:
            for l2 in lines[lno+1:]: 
                # I store the stripped, lower-case string
                retdict['warning_msg'].append(l2.strip().lower())
                # For now, I just take one line and stop.
                # I still do it in a for loop in case 'Exiting...' is on the
                # last row, to avoid exceptions. Remove the next 'break' if 
                # you want to get all lines
                break

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
