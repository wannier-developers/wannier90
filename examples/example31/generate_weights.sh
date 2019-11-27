#!/bin/bash
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
rm e.dat p.dat tmp.dat
