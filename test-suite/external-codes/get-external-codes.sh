#!/bin/bash
#Stop the full script if one line crashes
set -e

THEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$THEDIR"

## Get the external codes, but only if the environment variable
## W90TESTSWITHINTERFACE is set to true
if [ "$W90TESTSWITHINTERFACE" == "true" ]
then
    ./get-quantum-espresso.sh
fi
