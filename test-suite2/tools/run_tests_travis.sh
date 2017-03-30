#!/bin/bash

set -e

## Set here, if needed, the location of the executables
TESTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
cd "$TESTDIR"

# Default: serial, no mpirun. Run these in any case
./run_tests --category=default 
if [ "$W90BINARYPARALLEL" == "true" ]
then
  # If running in parallel: run also the tests in parallel
  # I hardcode this to four, in case change it or set it as an ENV
  # var in the .travis.yml
  make ./run_tests --category=default --num-procs=4
fi
