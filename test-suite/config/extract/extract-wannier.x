#!/bin/bash
#
# Copyright (C) 2001-2016 Wannier90 Developers Group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainers: Samuel Ponce
#

fname=$1
args=$(echo $fname | awk -F= '{print $NF}')

if [ "$fname" == "" ]
then
    echo "No file name passed, stopping..." >&2
    exit 1
fi

# WANNIER
# center and spread, but only for the first WF (this is suboptimal, should
# be improved...)
wfcenter=`sed -n '/Final State/{n;p;}' $fname | awk '{print $7; print $8; print $9}'| tr , " " `
spread=`sed -n '/Final State/{n;p;}' $fname | awk '{print $11}'`

omegaI=`grep "  Omega I " $fname | awk '{print $6}'`
omegaD=`grep "  Omega D " $fname | awk '{print $5}'`
omegaOD=`grep "  Omega OD " $fname | awk '{print $4}'`
omegaT=`grep "  Omega Total " $fname | awk '{print $7}'`

# Wannier -pp
nearn=`sed -n '/ Distance (Ang^-1)/{n;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $2}; {print $3}; {print $4}'`
compl=`sed -n '/ b_k(x) /{n;n;p;n;p;n;p;n;p;n;p;n;p;}' \
 $fname | awk '{print $2; print $3; print $4; print $5; print $6}'`



proji=`sed -n '/ Projections/{n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $1; print $2; print $3; print $4; print $5; print $6; print $7};'`

nnkpt_kpt=`sed -n '/begin kpoints/{n;n;p;n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $1; print $2; print $3};'`

nnkpt=`sed -n '/begin nnkpts/{n;n;p;n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $1; print $2; print $3; print $4; print $5};'`

wrmsg=`grep -A1 'Exiting.......'  "${fname}" | wc -l`



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


if test "$proji" != ""; then
        echo proji
        for x in $proji; do echo $x; done
fi

if test "$nnkpt_kpt" != ""; then
        echo nnkpt_kpt
        for x in $nnkpt_kpt; do echo $x; done
fi

if test "$nnkpt" != ""; then
        echo nnkpt
        for x in $nnkpt; do echo $x; done
fi

if test "$wrmsg" != ""; then
        echo wrmsg
        echo $wrmsg
fi


