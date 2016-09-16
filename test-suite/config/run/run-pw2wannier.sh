#!/bin/bash
#
# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainer: Samuel Ponce

#include ${ESPRESSO_ROOT}/test-suite/ENVIRONMENT
#bash ../ENVIRONMENT

if [[ $QE_USE_MPI == 1 ]]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
  export PARA_SUFFIX="-npool ${TESTCODE_NPROCS}"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

echo $0" "$@
if [[ "$1" == "1" ]]
then
  echo "Running pw2wannier90 ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw2wannier90.x < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw2wannier90.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

#rm -f input_tmp.in
