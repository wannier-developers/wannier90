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

# WANNIER
wfcenter=`sed -n '/Final State/{n;p;}' $fname | awk '{print $7; {print $8}; {print $9}}'`
spread=`sed -n '/Final State/{n;p;}' $fname | awk '{print $11}'`

omegaI=`grep "  Omega Total " $fname | awk '{print $6}'`
omegaD=`grep "  Omega Total " $fname | awk '{print $5}'`
omegaOD=`grep "  Omega Total " $fname | awk '{print $4}'`
omegaT=`grep "  Omega Total " $fname | awk '{print $7}'`

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




