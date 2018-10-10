"""
Parser function parse() to parse the .wout output file of Wannier90.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

# Match the lines describing the nearest-neighbour Shells
# Groups:
# 0: Shell Index
# 1: distance (ang^-1)
# 2: multiplicity
near_neigh_re = re.compile("^\s*\|\s+(\d+)\s+([\d\.]+)\s*(\d+)\s*")


# Match the lines describing the b_k vectors for the completeness relation
# Groups:
# 0: Index
# 1: b_k(x)
# 2: b_k(y)
# 3: b_k(z)
# 4: w_b
completeness_re = re.compile("^\s*\|\s+(\d+)\s+([\d\.-]+)\s+([\d\.-]+)\s+([\d\.-]+)\s*([\d\.]+)\s*")

# Match the 'WF centre and spread' line. 
# Groups:
# 0: idx
# 1: centre_x
# 2: centre_y
# 3: centre_z
# 4: spread
spread_re = re.compile("^\s*WF centre and spread\s+(\d+)\s+\(\s*([0-9\.-]+)\s*,\s*([0-9\.-]+)\s*,\s*([0-9\.-]+)\s*\)\s*([0-9\.-]+)\s*$")

# Match the lines with the Omegas
# Groups:
# 0: Omega_*
omegaI_re = re.compile("Omega\ I\s+=\s*([0-9\.-]+)\s*$")
omegaD_re = re.compile("Omega\ D\s+=\s*([0-9\.-]+)\s*$")
omegaOD_re = re.compile("Omega\ OD\s+=\s*([0-9\.-]+)\s*$")
omegaTotal_re = re.compile("Omega\ Total\s+=\s*([0-9\.-]+)\s*$")
omegaIOD_C_re = re.compile("Omega\ IOD_C\s+=\s*([0-9\.-]+)\s*$")
omegaRest_re = re.compile("Omega\ Rest\s+=\s*([0-9\.-]+)\s*$")
penaltyfunc_re = re.compile("Penalty\ func\s+=\s*([0-9\.-]+)\s*$")
omegaTotal_C_re = re.compile("Omega\ Total_C\s+=\s*([0-9\.-]+)\s*$")

## A comment on regexps: re.match only checks the beginning of the line, while
## re.search anywhere in the string (like perl)

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
        # Nearest-neighbour Shells
        # Start from the fourth line after 
        # 'Distance to Nearest-Neighbour Shells',
        # then stop at the line with ------------------
        if "Distance to Nearest-Neighbour Shells" in l:
            for l2 in lines[lno+4:]: # Skip 4 lines
                match = near_neigh_re.search(l2)
                if not match or '--------------------------------------' in l2:
                    break
                _, dist, mult = match.groups() 
                retdict["near_neigh_dist"].append(float(dist))
                retdict["near_neigh_mult"].append(int(mult))
            continue
        
        ###############################################################
        # Completeness relation
        # Start from the sixth line after 
        # 'Completeness relation is fully satisfied',
        # then stop at the line with ------------------
        if "Completeness relation is fully satisfied" in l:
            for l2 in lines[lno+6:]: # Skip 6 lines
                match = completeness_re.search(l2)
                if not match or '--------------------------------------' in l2:
                    break
                _, bkx, bky, bkz, bkw = match.groups() 
                retdict["completeness_x"].append(float(bkx))
                retdict["completeness_y"].append(float(bky))
                retdict["completeness_z"].append(float(bkz))
                retdict["completeness_weight"].append(float(bkw))
            continue

        ###############################################################
        # Final state spreads and centres: get all lines after
        # 'Final state' that contain the spreads
        if "Final State" in l:
            for l2 in lines[lno+1:]:
                match = spread_re.search(l2)
                if not match:
                    break
                _, x, y, z, s = match.groups()
                retdict["final_centres_x"].append(float(x))
                retdict["final_centres_y"].append(float(y))
                retdict["final_centres_z"].append(float(z))
                retdict["final_spreads"].append(float(s))
            continue
        ###############################################################
        # various Omegas (four numbers, typically at the end)
        match = omegaI_re.search(l)
        if match:
            retdict["omegaI"].append(float(match.groups()[0]))
            continue
        match = omegaD_re.search(l)
        if match:
            retdict["omegaD"].append(float(match.groups()[0]))
            continue
        match = omegaOD_re.search(l)
        if match:
            retdict["omegaOD"].append(float(match.groups()[0]))
            continue
        match = omegaTotal_re.search(l)
        if match:
            retdict["omegaTotal"].append(float(match.groups()[0]))
            continue
        match = omegaIOD_C_re.search(l)
        if match:
            retdict["omegaIOD_C"].append(float(match.groups()[0]))
            continue
        match = omegaRest_re.search(l)
        if match:
            retdict["omegaRest"].append(float(match.groups()[0]))
            continue
        match = penaltyfunc_re.search(l)
        if match:
            retdict["penaltyfunc"].append(float(match.groups()[0]))
            continue
        match = omegaTotal_C_re.search(l)
        if match:
            retdict["omegaTotal_C"].append(float(match.groups()[0]))
            continue
        ###############################################################
        

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
