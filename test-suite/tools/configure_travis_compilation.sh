#!/bin/bash

#Stop the full script if one line crashes
set -e

## Set here, if needed, the location of the executables
TESTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
cd "$TESTDIR"

# Move the correct make.inc to prepare compilation
if [ "$W90BINARYPARALLEL" == "true" ]
then
    cp config/TravisCI/make.inc.gfort+openmpi ../make.inc
else
    cp config/TravisCI/make.inc.gfort ../make.inc 
fi
