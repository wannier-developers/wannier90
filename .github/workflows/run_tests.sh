#!/bin/bash

set -e

TESTDIR=$(pwd)/test-suite/
echo $TESTDIR
cd "$TESTDIR"

# Default: serial, no mpirun. Run these in any case
echo "************************"
echo "* RUNNING SERIAL TESTS *"
echo "************************"
./run_tests --category=default

if [ "$W90BINARYPARALLEL" == "true" ]
then
    # If running in parallel: run also the tests in parallel
    echo ""
    echo "**************************"
    echo "* RUNNING PARALLEL TESTS *"
    echo "**************************"
    # I hardcode the numprocs to four, in case change it or set it as an ENV
    # var in the .travis.yml
    ./run_tests --category=par --numprocs=2
fi
