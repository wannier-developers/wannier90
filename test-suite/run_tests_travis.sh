#!/bin/bash

THEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$THEDIR"

## Get the external codes, but only if the environment variable
## W90TESTSWITHINTERFACE is set to true
if [ "$W90TESTSWITHINTERFACE" == "true" ]
then
    # Only wannier tests
    make run-tests
elif [ "$W90TESTSWITHINTERFACE" == "false" ]
then
    # Only tests involving interfaces
    make run-tests-interface
else
    # By default: run both
    make run-tests-all
fi
