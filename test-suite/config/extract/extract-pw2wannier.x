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

proji=`sed -n '/ Projections/{n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $1; print $2; print $3; print $4; print $5; print $6; print $7};'`

nnkpt_kpt=`sed -n '/begin kpoints/{n;n;p;n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $1; print $2; print $3};'`

nnkpt=`sed -n '/begin nnkpts/{n;n;p;n;p;n;p;n;p;n;p;}'\
 $fname | awk '{print $1; print $2; print $3; print $4; print $5};'`

wrmsg=`grep -A1 'Exiting.......'  "${fname}" | wc -l`


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
