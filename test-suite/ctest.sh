#!/bin/bash
# vim: dict+=/usr/share/beakerlib/dictionary.vim cpt=.,w,b,u,t,i,k
# shellcheck disable=all
source /usr/share/beakerlib/beakerlib.sh || exit 1

rlJournalStart
  rlPhaseStartSetup
    [[ -z "${MPI_VARIANT}" ]] && rlDie "MPI_VARIANT was not provided"
    rlRun 'tmp=$(mktemp -d)' 0 "Create tmp directory" || rlDie "Could not create tmp directory"
  rlPhaseEnd

  rlPhaseStartTest
    if [[ "${MPI_VARIANT}" != "serial" ]]; then
      rlRun 'module load "mpi/${MPI_VARIANT}-$(uname -m)"' 0 "Load MPI" || rlDie "Could not load MPI"
    fi
    rlRun 'cmake -S . -B "$tmp/build" -DWannier90_ROOT="${MPI_HOME:-/usr}"' 0 "CMake configure" &&
      rlRun 'cmake --build "$tmp/build"' 0 "CMake build" &&
      rlRun 'TMPDIR="$tmp" PYTHONDONTWRITEBYTECODE=1 ctest --test-dir "$tmp/build" --output-on-failure' 0 "Run CTest"
  rlPhaseEnd

  rlPhaseStartCleanup
    rlRun 'rm -r -- "$tmp"' 0 "Remove tmp directory"
  rlPhaseEnd
rlJournalEnd
