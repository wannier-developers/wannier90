"""
Parser function parse() to parse the <seedname>-shc-{fermiscan,freqscan}.dat output file of postw90.x.

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

        if lno == 0:
            if len(pieces) < 2:
                raise ValueError("Wrong SHC dat file: header is wrong")
            else:
                if pieces[1] == 'Fermi':
                    freqscan = False
                    continue
                elif pieces[1] == 'Frequency(eV)':
                    freqscan = True
                    continue
                else:
                    raise ValueError("Wrong SHC dat file: header is wrong")

        if len(pieces) == 0:
            # skip blank line
            continue

        if freqscan:
            if len(pieces) == 4 : 
                retdict['frequency'].append(float(pieces[1]))
                retdict['shc_re'].append(float(pieces[2]))
                retdict['shc_im'].append(float(pieces[3]))
            else:
                raise ValueError("Wrong line length ({}, instead of 4); line content: {}".format(
                    len(pieces), l))
        else:
            if len(pieces) == 3 : 
                retdict['energy'].append(float(pieces[1]))
                retdict['shc'].append(float(pieces[2]))
            else:
                raise ValueError("Wrong line length ({}, instead of 3); line content: {}".format(
                    len(pieces), l))


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
