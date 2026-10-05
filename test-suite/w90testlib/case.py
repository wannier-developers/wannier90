"""The TestCase model: one ``tests/<name>/test.yaml``, loaded and validated.

Migrated from the ``[section]`` blocks of the old ``tests/jobconfig``.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

#: Default wall-clock limit for a single run, in seconds.  Overridable per test.
DEFAULT_TIMEOUT = 600

#: Valid values of the ``status`` key.
STATUSES = ("enabled", "skip", "xfail")

#: Directory-name prefix -> pytest marker.  Order matters: ``testpostw90_`` must be tried
#: before ``testw90_`` would ever be ambiguous (it is not today, but keep it explicit).
MARKER_PREFIXES = (
    ("partest", "parallel"),
    ("testpostw90_", "postw90"),
    ("testw90_", "wannier90"),
    ("checkpoint", "checkpoint"),
)


class CaseError(Exception):
    """Raised when a test.yaml is malformed."""


@dataclass(frozen=True)
class Run:
    """A single invocation of the program under test.

    Every test currently has exactly one of these; the list form exists so a future test
    can declare several without a schema change.
    """

    input: str
    output: str
    benchmark: str
    args: tuple[str, ...] = ()


@dataclass(frozen=True)
class TestCase:
    """One regression test."""

    # Not a pytest test class, despite the name -- stop pytest trying to collect it.
    __test__ = False

    name: str
    directory: Path
    profile: str
    runs: tuple[Run, ...]
    description: str = ""
    serial_only: bool = False
    min_nprocs: int | None = None
    expect_failure: bool = False
    timeout: int = DEFAULT_TIMEOUT
    tags: tuple[str, ...] = ()
    depends_on: tuple[str, ...] = ()
    prepare: tuple[dict[str, Any], ...] = ()
    status: str = "enabled"
    status_reason: str = ""

    @property
    def markers(self) -> tuple[str, ...]:
        """Markers for this test: derived from the directory name, plus any ``tags``."""
        marks = [marker for prefix, marker in MARKER_PREFIXES if self.name.startswith(prefix)]
        if self.serial_only:
            marks.append("serial")
        marks.extend(self.tags)
        # dict.fromkeys keeps first-seen order while removing duplicates.
        return tuple(dict.fromkeys(marks))

    def benchmark_path(self, run: Run) -> Path:
        """Absolute path of a run's reference file, in the source tree."""
        return self.directory / run.benchmark

    def runs_at(self, nprocs: int | None) -> int | None:
        """How many processes this test should actually use.

        ``serial_only`` tests always run serially, even when --nprocs is given: they are the
        Gamma-only tests, which the old harness pinned with ``max_nprocs = 0``.
        """
        if self.serial_only:
            return None
        return nprocs

    def skip_reason_for(self, nprocs: int | None) -> str | None:
        """Why this test cannot run at this process count, if it cannot."""
        if self.status == "skip":
            return self.status_reason or "disabled in test.yaml"
        if self.min_nprocs is not None and (nprocs or 0) < self.min_nprocs:
            return (
                f"needs at least {self.min_nprocs} MPI processes "
                f"(run with --nprocs={self.min_nprocs} or more)"
            )
        return None


def _as_str_tuple(value: Any, what: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        raise CaseError(f"{what}: expected a list, got the string {value!r}")
    return tuple(str(v) for v in value)


def _run_from_dict(name: str, raw: dict[str, Any]) -> Run:
    unknown = set(raw) - {"input", "args", "output", "benchmark"}
    if unknown:
        raise CaseError(f"{name}: unknown field(s) in run: {', '.join(sorted(unknown))}")
    for required in ("input", "output"):
        if not raw.get(required):
            raise CaseError(f"{name}: run is missing required field {required!r}")
    output = str(raw["output"])
    return Run(
        input=str(raw["input"]),
        output=output,
        # The reference is a saved copy of the output file, so it is named after it.
        benchmark=str(raw.get("benchmark") or f"benchmark/{output}"),
        args=_as_str_tuple(raw.get("args"), f"{name}: run args"),
    )


def load_case(directory: Path) -> TestCase:
    """Load ``<directory>/test.yaml``."""
    path = directory / "test.yaml"
    if not path.is_file():
        raise CaseError(f"no test.yaml in {directory}")
    raw = yaml.safe_load(path.read_text()) or {}
    if not isinstance(raw, dict):
        raise CaseError(f"{path}: expected a mapping")

    known = {
        "description", "profile", "runs", "serial_only", "min_nprocs", "expect_failure",
        "timeout", "tags", "depends_on", "prepare", "status", "status_reason",
    }
    unknown = set(raw) - known
    if unknown:
        raise CaseError(f"{path}: unknown field(s): {', '.join(sorted(unknown))}")

    name = directory.name
    if not raw.get("profile"):
        raise CaseError(f"{path}: 'profile' is required")
    runs_raw = raw.get("runs")
    if not runs_raw:
        raise CaseError(f"{path}: 'runs' is required and must list at least one run")

    status = str(raw.get("status", "enabled"))
    if status not in STATUSES:
        raise CaseError(f"{path}: status must be one of {', '.join(STATUSES)}, got {status!r}")
    status_reason = str(raw.get("status_reason") or "")
    if status != "enabled" and not status_reason:
        raise CaseError(f"{path}: status is {status!r}, so status_reason must be non-empty")

    min_nprocs = raw.get("min_nprocs")
    return TestCase(
        name=name,
        directory=directory,
        profile=str(raw["profile"]),
        runs=tuple(_run_from_dict(name, r) for r in runs_raw),
        description=str(raw.get("description") or "").strip(),
        serial_only=bool(raw.get("serial_only", False)),
        min_nprocs=int(min_nprocs) if min_nprocs is not None else None,
        expect_failure=bool(raw.get("expect_failure", False)),
        timeout=int(raw.get("timeout") or DEFAULT_TIMEOUT),
        tags=_as_str_tuple(raw.get("tags"), f"{path}: tags"),
        depends_on=_as_str_tuple(raw.get("depends_on"), f"{path}: depends_on"),
        prepare=tuple(raw.get("prepare") or ()),
        status=status,
        status_reason=status_reason,
    )


def discover_cases(tests_root: Path) -> list[TestCase]:
    """Every test under ``tests_root``, ordered by name."""
    return [
        load_case(path.parent)
        for path in sorted(tests_root.glob("*/test.yaml"))
    ]
