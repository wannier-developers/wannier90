"""Declarative input preparation -- the replacement for the 55 per-test Makefiles.

Large inputs are stored bzip2-compressed, and checkpoints are stored *formatted* and
compressed because the binary ``.chk`` format is compiler- and machine-dependent.  Keep it
that way: never commit a binary ``.chk``.

All work happens inside the per-test work directory, never in the source tree.  Everything
is done with Python's ``bz2`` module rather than a shell pipe, so it also works where
``bunzip2`` is not on PATH.
"""

import bz2
import shutil
import subprocess
from pathlib import Path
from typing import Any, Callable

#: Extensions that are simply decompressed in place (``<name>.<ext>.bz2 -> <name>.<ext>``).
#: Derived from the wildcard rules in the old Makefiles.
BUNZIP2_SUFFIX = ".bz2"

#: Suffix identifying a compressed formatted checkpoint.
CHK_FMT_SUFFIX = ".chk.fmt.bz2"


class PrepareError(Exception):
    """Raised when a prepare step is malformed or fails."""


def _decompress(src: Path, dest: Path) -> None:
    """Stream-decompress ``src`` to ``dest`` (files here reach hundreds of MB)."""
    with bz2.open(src, "rb") as fh_in, open(dest, "wb") as fh_out:
        shutil.copyfileobj(fh_in, fh_out)


def step_bunzip2(work_dir: Path, pattern: str) -> list[Path]:
    """Decompress every file matching ``pattern``, dropping the ``.bz2`` suffix."""
    produced = []
    for src in sorted(work_dir.glob(pattern)):
        if src.suffix != BUNZIP2_SUFFIX:
            raise PrepareError(f"bunzip2 step matched {src.name}, which is not a .bz2 file")
        dest = src.with_suffix("")
        _decompress(src, dest)
        produced.append(dest)
    return produced


def step_chk_from_bz2(work_dir: Path, pattern: str, w90chk2chk: Path) -> list[Path]:
    """Turn ``<seed>.chk.fmt.bz2`` into a binary ``<seed>.chk``.

    Decompress to ``<seed>.chk.fmt``, convert with ``w90chk2chk.x -f2u <seed>``, then remove
    the formatted intermediate.  This is what the Makefiles' ``%.chk: %.chk.fmt.bz2`` rule did.
    """
    produced = []
    for src in sorted(work_dir.glob(pattern)):
        if not src.name.endswith(CHK_FMT_SUFFIX):
            raise PrepareError(
                f"chk_from_bz2 step matched {src.name}, which is not a {CHK_FMT_SUFFIX} file"
            )
        seed = src.name[: -len(CHK_FMT_SUFFIX)]
        formatted = work_dir / f"{seed}.chk.fmt"
        _decompress(src, formatted)
        result = subprocess.run(
            [str(w90chk2chk), "-f2u", seed],
            cwd=work_dir, capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise PrepareError(
                f"w90chk2chk.x -f2u {seed} failed with exit code {result.returncode} "
                f"in {work_dir}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        formatted.unlink()
        produced.append(work_dir / f"{seed}.chk")
    return produced


def step_copy_from_dependency(
    work_dir: Path,
    spec: dict[str, Any],
    resolve_dependency: Callable[[str], Path] | None,
) -> list[Path]:
    """Copy artefacts produced by another test's run.

    Used only by ``checkpoint1_read``, which consumes the ``copper.chk`` that
    ``checkpoint0_write`` writes.  The old Makefile symlinked into the sibling *source*
    directory, which only worked because tests ran in-tree and in order.
    """
    if resolve_dependency is None:
        raise PrepareError(
            "a copy_from_dependency step was requested but no dependency resolver was "
            "supplied; this test cannot run standalone"
        )
    dependency = spec.get("test")
    files = spec.get("files")
    if not dependency or not files:
        raise PrepareError("copy_from_dependency needs both 'test' and 'files'")

    source_dir = resolve_dependency(str(dependency))
    produced = []
    for name in files:
        src = source_dir / str(name)
        if not src.is_file():
            raise PrepareError(
                f"dependency {dependency!r} did not produce {name!r} (looked in {source_dir})"
            )
        dest = work_dir / str(name)
        shutil.copy2(src, dest)
        produced.append(dest)
    return produced


def run_prepare(
    steps,
    work_dir: Path,
    *,
    w90chk2chk: Path | None = None,
    resolve_dependency: Callable[[str], Path] | None = None,
) -> list[Path]:
    """Execute a test's ``prepare:`` list in order."""
    produced: list[Path] = []
    for step in steps:
        if not isinstance(step, dict) or len(step) != 1:
            raise PrepareError(
                f"each prepare step must be a single-key mapping, got {step!r}"
            )
        (kind, argument), = step.items()
        match kind:
            case "bunzip2":
                produced += step_bunzip2(work_dir, str(argument))
            case "chk_from_bz2":
                if w90chk2chk is None:
                    raise PrepareError(
                        "a chk_from_bz2 step needs w90chk2chk.x; pass --w90chk2chk-exe or "
                        "build it first"
                    )
                produced += step_chk_from_bz2(work_dir, str(argument), w90chk2chk)
            case "copy_from_dependency":
                produced += step_copy_from_dependency(work_dir, argument, resolve_dependency)
            case _:
                raise PrepareError(f"unknown prepare step {kind!r}")
    return produced
