#!/bin/bash

## Set here, if needed, the location of the executables
TESTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
export ESPRESSO_ROOT="${TESTDIR}/tools/external-codes/espresso/"

cd "$TESTDIR"

pwd

## Get the external codes, but only if the environment variable
## W90TESTSWITHINTERFACE is set to true
if [ "$W90TESTSWITHINTERFACE" == "true" ]
then
    # Only tests involving interfaces
    ./run_tests --category=run-tests-interface
elif [ "$W90TESTSWITHINTERFACE" == "false" ]
then
    # Only wannier tests
    ./run_tests --category=run-tests
else
    # By default: run default tests
    ./run_tests -d
fi
