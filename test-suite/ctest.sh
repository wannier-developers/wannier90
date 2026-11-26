#!/bin/bash
# vim: dict+=/usr/share/beakerlib/dictionary.vim cpt=.,w,b,u,t,i,k
# shellcheck disable=all
source /usr/share/beakerlib/beakerlib.sh || exit 1

rlJournalStart
  rlPhaseStartSetup
    [[ -n "${MPI_VARIANT}" ]] || rlDie "MPI_VARIANT was not provided"
    rlRun 'tmp=$(mktemp -d)' 0 "Create temporary build directory" || rlDie "mktemp failed"
    prefix=/usr
    if [[ "${MPI_VARIANT}" != "serial" ]]; then
      rlRun 'module load "mpi/${MPI_VARIANT}-$(uname -m)"' 0 "Load MPI" || rlDie "MPI module unavailable"
      prefix=${MPI_HOME:?MPI module did not set MPI_HOME}
    fi
  rlPhaseEnd

  rlPhaseStartTest
    rlRun 'cmake -S . -B "$tmp/build" -DWannier90_ROOT="$prefix" -DWANNIER90_TEST_KEEP_WORKDIR=OFF' 0 "CMake configure" &&
      rlRun 'cmake --build "$tmp/build"' 0 "CMake build" &&
      rlRun 'TMPDIR="$tmp" PYTHONDONTWRITEBYTECODE=1 ctest --test-dir "$tmp/build" --output-on-failure' 0 "Run CTest"
  rlPhaseEnd

  rlPhaseStartCleanup
    rlRun 'rm -r -- "$tmp"' 0 "Remove temporary build directory"
  rlPhaseEnd
rlJournalEnd
