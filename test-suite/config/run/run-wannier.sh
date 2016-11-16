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
  echo "Running Wannier90 ..."
  export TMP1=$2
  export OUTPUT="${TMP1%??}out"
  echo "${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x $2 $3 2> $4"
  ${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x $2 $3 2> $4
  cp $OUTPUT $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
elif [[ "$1" == "2" ]]
then
  echo "Running PW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "3" ]]
then
  echo "Running pw2wannier90 ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw2wannier90.x < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw2wannier90.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "4" ]]
then
  echo "Running PP wannier ..."
  export TMP1=$2
  export OUTPUT="${TMP1%??}out"
  echo "${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x -pp $2 $3 2> $4"
  ${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x -pp $2 $3 2> $4
  cp $OUTPUT $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
elif [[ "$1" == "5" ]]
then
  echo "Running PP wannier and checking for nnkp ..."
  export TMP1=$2
  export OUTPUT="${TMP1%???}nnkp"
  echo "${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x -pp $2 $3 2> $4"
  ${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x -pp $2 $3 2> $4
  cp $OUTPUT $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
elif [[ "$1" == "6" ]]
then
  echo "Running PP wannier and checking for crash ..."
  export TMP1=$2
  export OUTPUT="${TMP1%???}werr"
  echo "${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x -pp $2 $3 2> $4"
  ${PARA_PREFIX} ${WANNIER_ROOT}/wannier90.x $2 $3 2> $4
  cp $OUTPUT $3
#  if [[ -e CRASH ]]
#  then
#    cat $3
#  fi  
elif [[ "$1" == "7" ]]
then
  echo "Running w90+postw90 for geninterp and checking for energies and derivatives ..."
  # SP: DO NOT set this to TMP or it will clash with mpirun
  export TMP1=$2
  export OUTPUT="${TMP1%????}_geninterp.dat"
  echo " ${WANNIER_ROOT}/wannier90.x $2 2> $4 # $3"
  ${WANNIER_ROOT}/wannier90.x $2 2> ${4}_wan
  echo "${PARA_PREFIX} ${WANNIER_ROOT}postw90.x $2 2> $4 # $3"
  ${PARA_PREFIX} ${WANNIER_ROOT}postw90.x $2 2> $4
  cp $OUTPUT $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
fi

#rm -f input_tmp.in
