"""
Parser function parse() to parse the seedname_hr.dat output file of Wannier90 (<0i|H|Rj> matrix elements).
"""
import inspect
from collections import defaultdict

from . import show_output

def parse(fname):
    """
    Open the file, parse it and return the values.

    The layout is: comment line, num_wann, nrpts, the ndegen block wrapped at
    15 entries per line, then num_wann^2 * nrpts matrix-element lines.
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    # skip the first line, which contains only the date
    num_wann = int(lines[1].split()[0])
    nrpts = int(lines[2].split()[0])
    retdict['num_wann'] = [num_wann]
    retdict['nrpts'] = [nrpts]

    ndegen = []
    iline = 3
    while len(ndegen) < nrpts:
        ndegen += [int(_) for _ in lines[iline].split()]
        iline += 1
    retdict['ndegen'] = [float(_) for _ in ndegen]

    for line in lines[iline:iline + num_wann * num_wann * nrpts]:
        fields = line.split()
        retdict['irvec_a'].append(float(fields[0]))
        retdict['irvec_b'].append(float(fields[1]))
        retdict['irvec_c'].append(float(fields[2]))
        retdict['index_i'].append(float(fields[3]))
        retdict['index_j'].append(float(fields[4]))
        retdict['real_h'].append(float(fields[5]))
        retdict['imag_h'].append(float(fields[6]))

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
