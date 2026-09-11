#!/bin/bash

set -e

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO_ROOT"

# Artefacts land outside the source tree so a failed run can be uploaded wholesale and
# `git status` stays clean.
ARTIFACTS="$REPO_ROOT/test-artifacts"

echo "************************"
echo "* RUNNING SERIAL TESTS *"
echo "************************"
python -m pytest test-suite/tests --workdir="$ARTIFACTS/serial" -ra

if [ "$W90BINARYPARALLEL" == "true" ]
then
    echo ""
    echo "**************************"
    echo "* RUNNING PARALLEL TESTS *"
    echo "**************************"
    python -m pytest test-suite/tests --nprocs=2 \
        -m "wannier90 or postw90 or checkpoint or parallel" \
        --workdir="$ARTIFACTS/parallel" -ra
fi
