#!/bin/bash

set -e

## Set here, if needed, the location of the executables
TESTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
cd "$TESTDIR"

echo "****************************************"
echo "* RUNNING SERIAL TEST FOR LIBRARY MODE *"
echo "****************************************"
cd library-mode-test
# Create the binary for testing
make
# remove the output file from old runs, if there
rm -f results.dat
# Run it - if it crashes it should give a non-zero error
./test_library_serial.x
# Check the output values
./compare_results.py

## Going back to the test dir
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
    ./run_tests --category=default --numprocs=4
fi
