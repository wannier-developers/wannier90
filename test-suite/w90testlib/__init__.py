"""Helpers for the Wannier90 regression suite.

The suite is driven by pytest; this package holds the pieces that do the actual work, so
they can be unit-tested independently of any test run:

* :mod:`compare`     -- the tolerance engine
* :mod:`profiles`    -- profiles.yaml: program, parser and tolerances
* :mod:`case`        -- test.yaml: one regression test
* :mod:`prepare`     -- declarative input preparation (replaces the per-test Makefiles)
* :mod:`runner`      -- subprocess / mpirun invocation
* :mod:`executables` -- locating wannier90.x, postw90.x and w90chk2chk.x
* :mod:`parsers`     -- output parsers, one module per file format
"""

from .case import TestCase, discover_cases, load_case
from .compare import ComparisonResult, Tolerance, compare, format_failure
from .executables import Executables, MissingExecutableError
from .prepare import run_prepare
from .profiles import Profile, load_profiles
from .runner import RunResult, build_argv, run

__all__ = [
    "ComparisonResult",
    "Executables",
    "MissingExecutableError",
    "Profile",
    "RunResult",
    "TestCase",
    "Tolerance",
    "build_argv",
    "compare",
    "discover_cases",
    "format_failure",
    "load_case",
    "load_profiles",
    "run",
    "run_prepare",
]
