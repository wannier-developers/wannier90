#!/bin/bash

## Set here, if needed, the location of the executables
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export ESPRESSO_ROOT="${THISDIR}/external-codes/espresso/"

cd "$THISDIR"

## Get the external codes, but only if the environment variable
## W90TESTSWITHINTERFACE is set to true
if [ "$W90TESTSWITHINTERFACE" == "true" ]
then
    # Only tests involving interfaces
    make run-tests-interface
elif [ "$W90TESTSWITHINTERFACE" == "false" ]
then
    # Only wannier tests
    make run-tests
else
    # By default: run both
    make run-tests-all
fi
