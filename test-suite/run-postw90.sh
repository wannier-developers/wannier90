#!/bin/bash
#
# Copyright (C) 2016 Wannier90 developers
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#

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
  echo "Running w90+postw90 for geninterp and checking for energies and derivatives ..."
  export TMP=$2
  export OUTPUT="${TMP%????}_geninterp.dat"
  echo "${WANNIER_ROOT}/wannier90.x $2 2> $4 # $3"
  ${WANNIER_ROOT}/wannier90.x $2 2> $4
  echo "${PARA_PREFIX} ${WANNIER_ROOT}/postw90.x $2 2> $4 # $3"
  ${PARA_PREFIX} ${WANNIER_ROOT}/postw90.x $2 2> $4
  cp $OUTPUT $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
fi


