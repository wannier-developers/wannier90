#!/bin/bash
# This script is used to extract the projectability from
# the output of projwfc.x (proj.out) and it works only
# for this specific system with 21 bands.

# Remove temporary files if exist
rm -f e.dat
rm -f p.dat
rm -f tmp.dat
rm -f p_vs_e.dat

# Check proj.out exists
[[ -f "proj.out" ]] || { echo "proj.out not found!"; echo "Aborting!"; exit 1; }

# Get energies and projectability in the correct order
cat proj.out |grep '=='|awk '{print $5}' > e.dat
cat proj.out |grep '|psi|^2'|awk '{print $3}' > p.dat
paste e.dat p.dat > tmp.dat

sort -k1n tmp.dat > p_vs_e.dat 

# Clean workspace
rm e.dat p.dat tmp.dat
