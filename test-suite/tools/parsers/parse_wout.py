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
near_neigh_re = re.compile(r"^\s*\|\s+(\d+)\s+([\d\.]+)\s*(\d+)\s*")

# Match the 'WF centre and spread' line. 
# Groups:
# 0: idx
# 1: centre_x
# 2: centre_y
# 3: centre_z
# 4: spread
spread_re = re.compile(r"^\s*WF centre and spread\s+(\d+)\s+\(\s*([0-9\.-]+)\s*,\s*([0-9\.-]+)\s*,\s*([0-9\.-]+)\s*\)\s*([0-9\.-]+)\s*$")

# Match the lines with the Omegas
# Groups:
# 0: Omega_*
omegaI_re = re.compile(r"Omega\ I\s+=\s*([0-9\.-]+)\s*$")
omegaD_re = re.compile(r"Omega\ D\s+=\s*([0-9\.-]+)\s*$")
omegaOD_re = re.compile(r"Omega\ OD\s+=\s*([0-9\.-]+)\s*$")
omegaTotal_re = re.compile(r"Omega\ Total\s+=\s*([0-9\.-]+)\s*$")
omegaIOD_C_re = re.compile(r"Omega\ IOD_C\s+=\s*([0-9\.-]+)\s*$")
omegaRest_re = re.compile(r"Omega\ Rest\s+=\s*([0-9\.-]+)\s*$")
penaltyfunc_re = re.compile(r"Penalty\ func\s+=\s*([0-9\.-]+)\s*$")
omegaTotal_C_re = re.compile(r"Omega\ Total_C\s+=\s*([0-9\.-]+)\s*$")

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
        # b_k vectors
        # Start from the fourth line after
        # 'b_k Vectors (Ang^-1) and Weights (Ang^2)',
        # then stop at the line with ------------------
        if "b_k Vectors (Ang^-1) and Weights (Ang^2)" in l:
            for l2 in lines[lno+4:]: # Skip 6 lines
                data = l2.split()
                if len(data) != 7 or '-' * 76 in l2:
                    break
                _, bkx, bky, bkz, bkw = data[1:-1]
                retdict["b_k_x"].append(float(bkx))
                retdict["b_k_y"].append(float(bky))
                retdict["b_k_z"].append(float(bkz))
                retdict["w_b"].append(float(bkw))
            continue

        ###############################################################
        # b_k directions
        # Start from the fourth line after
        # 'b_k Directions (Ang^-1)',
        # then stop at the line with ------------------
        if "b_k Directions (Ang^-1)" in l:
            for l2 in lines[lno+4:]: # Skip 6 lines
                data = l2.split()
                if len(data) != 6 or '-' * 76 in l2:
                    break
                _, bkx, bky, bkz = data[1:-1]
                retdict["b_k_direction_x"].append(float(bkx))
                retdict["b_k_direction_y"].append(float(bky))
                retdict["b_k_direction_z"].append(float(bkz))
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
