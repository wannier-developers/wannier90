"""Running one test case out of tree.

Every run happens in a fresh work directory, never in the source tree: no dirty git status,
no contamination between runs, and the last few runs stay on disk for post-mortem debugging.
"""

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .case import TestCase
from .executables import Executables
from .prepare import run_prepare
from .profiles import Profile
from .runner import DEFAULT_MPI_LAUNCHER, RunResult, build_argv, run

#: Never copied into the work directory: the reference is read from the source tree, and
#: copying it in would let a test accidentally read its own answer key.
EXCLUDED_FROM_WORKDIR = ("benchmark", "test.yaml", "Makefile", ".gitignore")


@dataclass(frozen=True)
class RunOptions:
    """How to run, independent of which test is running."""

    nprocs: int | None = None
    launcher: str = DEFAULT_MPI_LAUNCHER
    timeout: int | None = None  # None -> use each case's own timeout


def populate_work_dir(case: TestCase, work_dir: Path) -> None:
    """Copy a test's inputs into ``work_dir``.

    ``symlinks=False`` materialises symlinks into real files.  The 156-odd symlinks under
    tests/ point at ../../checkpoints/... and at sibling test directories, so copying them
    as links would leave them dangling.  This mirrors what ctest.sh does with `rsync -rL`.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        case.directory, work_dir,
        symlinks=False,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns(*EXCLUDED_FROM_WORKDIR),
    )


def execute_case(
    case: TestCase,
    profile: Profile,
    executables: Executables,
    work_dir: Path,
    options: RunOptions = RunOptions(),
    *,
    resolve_dependency: Callable[[str], Path] | None = None,
) -> list[RunResult]:
    """Populate a work directory, run the prepare steps, then run the program."""
    populate_work_dir(case, work_dir)

    if case.prepare:
        run_prepare(
            case.prepare, work_dir,
            w90chk2chk=executables.require("w90chk2chk") if _needs_chk(case) else None,
            resolve_dependency=resolve_dependency,
        )

    exe = executables.require(profile.program)
    results = []
    for spec in case.runs:
        argv = build_argv(
            exe, spec.input, spec.args,
            nprocs=case.runs_at(options.nprocs),
            launcher=options.launcher,
        )
        results.append(run(argv, work_dir, timeout=options.timeout or case.timeout))
    return results


def _needs_chk(case: TestCase) -> bool:
    return any("chk_from_bz2" in step for step in case.prepare)
