"""Tolerance engine.

This is a deliberate, literal reimplementation of the comparison semantics of
``testcode2.validation`` (the harness this suite used until 2026).  Every rule below was
checked against that source before being written here; see ``README.md`` and the notes in
each function.  Do not "improve" the arithmetic: the reference outputs in ``tests/*/benchmark/``
were accepted under exactly these rules, and changing them would silently reclassify tests.
"""

import math
from dataclasses import dataclass, field
from typing import Any

#: Tolerance applied to any key a profile does not mention.  From testcode's built-in
#: default user option ``tolerance = '(1.e-10,None)'`` (testcode2/config.py).  Note this
#: means unlisted keys are still checked, and checked strictly.
DEFAULT_TOLERANCE_ABS = 1.0e-10
DEFAULT_TOLERANCE_REL = None

#: Cap on how many individual deviations get rendered in a failure message.
MAX_REPORTED_DEVIATIONS = 50


@dataclass(frozen=True)
class Tolerance:
    """Thresholds for comparing one key.

    ``None`` means "do not check this at all" -- distinct from ``0.0``, which is a real
    threshold that nothing can pass (the test is ``err < tol``).  ``rel_tol=None`` is used
    deliberately by some profiles for quantities that pass through zero.

    ``strict`` (testcode's default, and the only value the old ``userconfig`` ever used)
    requires *both* thresholds to be met when both are set; with ``strict=False`` either
    one passing is enough.
    """

    abs_tol: float | None = DEFAULT_TOLERANCE_ABS
    rel_tol: float | None = DEFAULT_TOLERANCE_REL
    strict: bool = True

    def describe(self) -> str:
        parts = []
        parts.append(f"abs={self.abs_tol:.3e}" if self.abs_tol is not None else "abs=off")
        parts.append(f"rel={self.rel_tol:.3e}" if self.rel_tol is not None else "rel=off")
        if not self.strict:
            parts.append("non-strict")
        return ", ".join(parts)


#: The tolerance used for keys a profile does not name.
DEFAULT_TOLERANCE = Tolerance()


@dataclass
class Deviation:
    """One element that failed its tolerance check."""

    key: str
    index: int
    test: Any
    benchmark: Any
    reason: str
    abs_err: float | None = None
    rel_err: float | None = None
    tolerance: Tolerance | None = None

    def render(self) -> str:
        head = f"  {self.key}[{self.index}]: test={self.test!r} benchmark={self.benchmark!r}"
        if self.tolerance is None:
            return f"{head}\n      {self.reason}"
        bits = []
        if self.abs_err is not None:
            limit = "off" if self.tolerance.abs_tol is None else f"{self.tolerance.abs_tol:.3e}"
            bits.append(f"abs err {self.abs_err:.6e} (limit {limit})")
        if self.rel_err is not None:
            limit = "off" if self.tolerance.rel_tol is None else f"{self.tolerance.rel_tol:.3e}"
            bits.append(f"rel err {self.rel_err:.6e} (limit {limit})")
        detail = "; ".join(bits) if bits else self.reason
        return f"{head}\n      {detail}"


@dataclass
class ComparisonResult:
    """Outcome of comparing a parsed output against a parsed benchmark."""

    deviations: list[Deviation] = field(default_factory=list)
    structural_errors: list[str] = field(default_factory=list)
    compared: int = 0

    @property
    def ok(self) -> bool:
        return not self.deviations and not self.structural_errors


def _is_number(value: Any) -> bool:
    """True for values testcode would have done arithmetic on.

    testcode distinguished numbers from everything else by catching ``TypeError`` from
    ``test - benchmark``; this is the same partition, spelled explicitly.  ``bool`` is a
    subclass of ``int`` and would have taken the numeric path there, so it does here too.
    """
    return isinstance(value, (int, float))


def absolute_error(test: float, benchmark: float) -> float:
    """``|test - benchmark|``."""
    return abs(test - benchmark)


def relative_error(test: float, benchmark: float) -> float:
    """Error normalised by the *benchmark* (not the mean, not the max).

    A zero benchmark with a zero difference is an exact match; a zero benchmark with any
    difference is infinitely wrong.  This is testcode's rule and is *not* the same as
    ``math.isclose``, which normalises by the larger magnitude -- do not substitute it.
    """
    diff = test - benchmark
    if benchmark == 0 and diff == 0:
        return 0.0
    if benchmark == 0:
        return math.inf
    return abs(diff / benchmark)


def check_value(test: Any, benchmark: Any, tol: Tolerance) -> Deviation | None:
    """Compare one element.  Returns ``None`` when it passes."""
    if not (_is_number(test) and _is_number(benchmark)):
        # Non-numeric (parsers return some strings and ints): plain equality, as testcode
        # did via its TypeError fallback.  A number vs a string also lands here and fails.
        if test == benchmark:
            return None
        return Deviation(key="", index=-1, test=test, benchmark=benchmark,
                         reason="values differ (compared for exact equality)")

    if math.isnan(test) or math.isnan(benchmark):
        # NaN never compares equal, so it can never be "within tolerance".
        return Deviation(key="", index=-1, test=test, benchmark=benchmark,
                         reason="cannot compare NaN")

    abs_err = absolute_error(test, benchmark)
    rel_err = relative_error(test, benchmark)

    # Strict inequality throughout: a value passes when err < tol, never when err == tol.
    abs_ok = True if tol.abs_tol is None else abs_err < tol.abs_tol
    rel_ok = True if tol.rel_tol is None else rel_err < tol.rel_tol

    if tol.abs_tol is not None and tol.rel_tol is not None and not tol.strict:
        passed = abs_ok or rel_ok
    else:
        # One threshold set (only the active one matters), or both set under strict mode.
        passed = abs_ok and rel_ok

    if passed:
        return None
    return Deviation(
        key="", index=-1, test=test, benchmark=benchmark,
        reason="outside tolerance",
        abs_err=abs_err if tol.abs_tol is not None else None,
        rel_err=rel_err if tol.rel_tol is not None else None,
        tolerance=tol,
    )


def compare(
    test: dict[str, list[Any]],
    benchmark: dict[str, list[Any]],
    tolerances: dict[str, Tolerance] | None = None,
    default_tolerance: Tolerance = DEFAULT_TOLERANCE,
) -> ComparisonResult:
    """Compare two parsed dictionaries key-by-key, element-by-element."""
    tolerances = tolerances or {}
    result = ComparisonResult()

    # An empty parse is a failure, not a pass.  testcode reported "no data" as success,
    # which let a test keep "passing" forever after a regex stopped matching.
    if not benchmark or all(len(v) == 0 for v in benchmark.values()):
        result.structural_errors.append(
            "the parser extracted no data from the benchmark file: either the reference is "
            "empty or the parser no longer matches its format"
        )
        return result

    # Key sets and per-key lengths must match exactly.  A mismatch is a failure, not a
    # silently truncated comparison.
    test_keys, bench_keys = set(test), set(benchmark)
    if test_keys != bench_keys:
        only_bench = sorted(bench_keys - test_keys)
        only_test = sorted(test_keys - bench_keys)
        if only_bench:
            result.structural_errors.append(
                f"keys present in benchmark but missing from test output: {', '.join(only_bench)}")
        if only_test:
            result.structural_errors.append(
                f"keys present in test output but missing from benchmark: {', '.join(only_test)}")

    for key in sorted(test_keys & bench_keys):
        if len(test[key]) != len(benchmark[key]):
            result.structural_errors.append(
                f"key {key!r}: test output has {len(test[key])} value(s), "
                f"benchmark has {len(benchmark[key])}"
            )

    if result.structural_errors:
        return result

    for key in sorted(bench_keys):
        tol = tolerances.get(key, default_tolerance)
        for index, (test_value, bench_value) in enumerate(zip(test[key], benchmark[key])):
            result.compared += 1
            deviation = check_value(test_value, bench_value, tol)
            if deviation is not None:
                deviation.key = key
                deviation.index = index
                result.deviations.append(deviation)

    return result


def format_failure(result: ComparisonResult, *, work_dir=None, output=None, benchmark=None) -> str:
    """Render a comparison failure for a pytest assertion message."""
    lines: list[str] = []
    if result.structural_errors:
        lines.append("Parsed data does not line up with the benchmark:")
        lines.extend(f"  {e}" for e in result.structural_errors)

    if result.deviations:
        total = len(result.deviations)
        lines.append(f"{total} value(s) outside tolerance (of {result.compared} compared):")
        for deviation in result.deviations[:MAX_REPORTED_DEVIATIONS]:
            lines.append(deviation.render())
        if total > MAX_REPORTED_DEVIATIONS:
            lines.append(f"  ... and {total - MAX_REPORTED_DEVIATIONS} more")

    if output is not None:
        lines.append(f"output:    {output}")
    if benchmark is not None:
        lines.append(f"benchmark: {benchmark}")
    if work_dir is not None:
        lines.append(f"work dir:  {work_dir}")
    return "\n".join(lines)
