#!/bin/bash
# This script is used to extract the projectability from
# the output of projwfc.x (proj.out) and it works only
# for this specific system with 21 bands.

# Remove temporary files if exist
if [ -f "e.dat" ]; then
	rm e.dat;
fi

if [ -f "p.dat" ]; then
        rm p.dat;
fi

if [ -f "p_vs_e.dat" ]; then
        rm p_vs_e.dat;
fi

# Check proj.out exists
[[ -f "proj.out" ]] || { echo "proj.out not found!"; echo "Aborting!"; exit; }

# Start loop over bands
for i in {1..21}; do
        echo $i 
	if [ $i -lt 10 ]; then 
           line=`echo "e(   $i)"`;
        else
           line=`echo "e(  $i)"`;
	fi;
	grep "$line" proj.out | awk '{print $5}' >> e.dat; 
        grep -A4 "$line" proj.out | grep "|psi|^2" | awk '{print $3}' >> p.dat; 
	paste e.dat p.dat >> tmp.dat;
done
sort -k1n tmp.dat > p_vs_e.dat 

# Clean workspace
rm e.dat p.dat tmp.dat
