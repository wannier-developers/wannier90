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
# Maintainers: Giovanni Pizzi

fname=$1

if [ "$fname" == "" ]
then
    echo "No file name passed, stopping..." >&2
    exit 1
fi

# POSTW90 - BANDS
bandidx=`grep -v '#' $fname | awk '{print $1}'`
bandkpt=`grep -v '#' $fname | awk '{print $2, $3, $4}'`
bandenergy=`grep -v '#' $fname | awk '{print $5}'`
bandderiv=`grep -v '#' $fname | awk '{print $6, $7, $8}'`

if test "$bandidx" != ""; then
        echo bandidx
        echo $bandidx | tr ' ' '\n'
fi

if test "$bandkpt" != ""; then
        echo bandkpt
        echo $bandkpt | tr ' ' '\n'
fi

if test "$bandenergy" != ""; then
        echo bandenergy
        echo $bandenergy | tr ' ' '\n'
fi

if test "$bandderiv" != ""; then
        echo bandderiv
        echo $bandderiv | tr ' ' '\n'
fi

