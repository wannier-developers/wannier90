#!/bin/sh
#
# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainers: Samuel Ponce
#
# SP: This can lead to issue if you reach the OS pipe buffer 
#     You can increase the buffer in /proc/sys/fs/pipe-max-size  

fname=$1
args=$(echo $fname | awk -F= '{print $NF}')

if [ "$fname" == "" ]
then
    echo "No file name passed, stopping..." >&2
    exit 1
fi

# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`

# NSCF
ef1=`grep "the Fermi energy is" $fname | awk '{print $5}'`
eh1=`grep "highest occupied" $fname | awk '{print $5}'`
el1=`grep "highest occupied" $fname | awk '{print $8}'`
tf1=`grep " P = " $fname | head -1 | awk '{printf "%7.5f", $3}'`
eigenval=`sed -n '/   k =/{n;n;p;}' $fname | awk '{print $1}; {print $2}; {print $3}; {print $4}'`

# WANNIER
wfcenter=`sed -n '/Final State/{n;p;}' $fname | awk '{print $7; {print $8}; {print $9}}'`
spread=`sed -n '/Final State/{n;p;}' $fname | awk '{print $11}'`

omegaI=`grep "  Omega Total " $fname | awk '{print $6}'`
omegaD=`grep "  Omega Total " $fname | awk '{print $5}'`
omegaOD=`grep "  Omega Total " $fname | awk '{print $4}'`
omegaT=`grep "  Omega Total " $fname | awk '{print $7}'`

# Wannier -pp
nearn=`sed -n '/ Distance (Ang^-1)/{n;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p}'\
 $fname | awk '{print $2}; {print $3}; {print $4}'`
compl=`sed -n '/ b_k(x) /{n;n;p;n;p;n;p;n;p;n;p;n;p;}' \
 $fname | awk '{print $2}; {print $3}; {print $4}; {print $5}; {print $6}'`

# PW2WANNIER
proji=`sed -n '/ Projections/{n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $1}; {print $2}; {print $3}; {print $4}; {print $5}; {print $6}; {print $7};'`



if test "$e1" != ""; then
        echo e1
        echo $e1
fi

if test "$n1" != ""; then
        echo n1
        echo $n1
fi

if test "$f1" != ""; then
        echo f1
        echo $f1
fi

if test "$p1" != ""; then
        echo p1
        echo $p1
fi

# NSCF

if test "$ef1" != ""; then
        echo ef1
        for x in $ef1; do echo $x; done
fi

if test "$eh1" != ""; then
        echo eh1
        for x in $eh1; do echo $x; done
fi

if test "$el1" != ""; then
        echo el1
        for x in $el1; do echo $x; done
fi

if test "$tf1" != ""; then
        echo tf1
        for x in $tf1; do echo $x; done
fi

if test "$eigenval" != ""; then
        echo eigenval
        for x in $eigenval; do echo $x; done
fi

# Wannier

if test "$wfcenter" != ""; then
        echo wfcenter
        for x in $wfcenter; do echo $x; done
fi

if test "$omegaI" != ""; then
        echo omegaI
        echo $omegaI
fi

if test "$omegaD" != ""; then
        echo omegaD
        echo $omegaD
fi


if test "$omegaOD" != ""; then
        echo omegaOD
        echo $omegaOD
fi


if test "$omegaT" != ""; then
        echo omegaT
        echo $omegaT
fi

# Wannier -pp

if test "$nearn" != ""; then
        echo nearn
        for x in $nearn; do echo $x; done
fi

if test "$compl" != ""; then
        echo compl
        for x in $compl; do echo $x; done
fi

# PW2WANNIER

if test "$proji" != ""; then
        echo proji
        for x in $proji; do echo $x; done
fi


