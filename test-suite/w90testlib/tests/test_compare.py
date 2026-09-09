"""Unit tests for the tolerance engine.

These pin down the semantics inherited from testcode2.  Each test names the rule it
protects; if one of these fails after a refactor, the refactor is wrong, not the test.
"""

import math

import pytest

from w90testlib.compare import (
    DEFAULT_TOLERANCE,
    MAX_REPORTED_DEVIATIONS,
    Tolerance,
    absolute_error,
    check_value,
    compare,
    format_failure,
    relative_error,
)


# --- the default tolerance -------------------------------------------------------------

def test_default_tolerance_is_abs_1e10_with_no_relative_check():
    assert DEFAULT_TOLERANCE.abs_tol == 1.0e-10
    assert DEFAULT_TOLERANCE.rel_tol is None
    assert DEFAULT_TOLERANCE.strict is True


def test_keys_absent_from_the_profile_are_still_checked_strictly():
    """Unlisted keys get the default tolerance -- they are NOT skipped."""
    result = compare({"x": [1.0 + 1e-8]}, {"x": [1.0]}, tolerances={})
    assert not result.ok
    assert result.deviations[0].key == "x"


# --- absolute and relative error -------------------------------------------------------

def test_absolute_error_is_the_magnitude_of_the_difference():
    assert absolute_error(1.5, 1.0) == pytest.approx(0.5)
    assert absolute_error(1.0, 1.5) == pytest.approx(0.5)


def test_relative_error_is_normalised_by_the_benchmark_not_the_mean():
    # 1.0 vs 2.0 -> 0.5 when normalised by the benchmark, 0.5 either way here, so use a
    # case where the choice matters: normalising by the test value would give 1.0.
    assert relative_error(2.0, 1.0) == pytest.approx(1.0)
    assert relative_error(1.0, 2.0) == pytest.approx(0.5)


def test_zero_benchmark_with_zero_difference_is_an_exact_match():
    assert relative_error(0.0, 0.0) == 0.0


def test_zero_benchmark_with_any_difference_is_infinitely_wrong():
    assert relative_error(1e-30, 0.0) == math.inf


def test_relative_check_is_not_math_isclose():
    """isclose normalises by the larger magnitude; we must normalise by the benchmark."""
    test, benchmark = 1e-9, 0.0
    assert math.isclose(test, benchmark, abs_tol=1e-8)  # isclose would pass this
    assert relative_error(test, benchmark) == math.inf  # our rule fails it


# --- strict inequality -----------------------------------------------------------------

def test_error_exactly_equal_to_the_tolerance_fails():
    """The test is err < tol, never err <= tol."""
    tol = Tolerance(abs_tol=0.5, rel_tol=None)
    assert check_value(1.5, 1.0, tol) is not None


def test_error_just_below_the_tolerance_passes():
    tol = Tolerance(abs_tol=0.5, rel_tol=None)
    assert check_value(1.0 + 0.4999, 1.0, tol) is None


def test_a_zero_tolerance_can_never_be_passed():
    """0.0 is a real threshold (nothing satisfies err < 0), distinct from None."""
    assert check_value(1.0, 1.0, Tolerance(abs_tol=0.0, rel_tol=None)) is not None


# --- None means "do not check" ---------------------------------------------------------

def test_none_relative_tolerance_disables_the_relative_check():
    """Used deliberately for quantities that pass through zero."""
    tol = Tolerance(abs_tol=1e-3, rel_tol=None)
    assert check_value(1e-9, 0.0, tol) is None  # rel err is inf, but rel is switched off


def test_none_absolute_tolerance_disables_the_absolute_check():
    tol = Tolerance(abs_tol=None, rel_tol=1e-3)
    assert check_value(1000.0, 1000.01, tol) is None  # abs err 0.01, but abs is off


def test_none_is_distinct_from_zero():
    assert check_value(5.0, 5.0, Tolerance(abs_tol=None, rel_tol=None)) is None
    assert check_value(5.0, 5.0, Tolerance(abs_tol=0.0, rel_tol=None)) is not None


# --- combining the two thresholds ------------------------------------------------------

def test_strict_mode_requires_both_thresholds():
    tol = Tolerance(abs_tol=1.0, rel_tol=1e-12, strict=True)
    # abs err 0.5 passes, rel err 0.5 fails -> overall fail
    assert check_value(1.5, 1.0, tol) is not None


def test_non_strict_mode_accepts_either_threshold():
    tol = Tolerance(abs_tol=1.0, rel_tol=1e-12, strict=False)
    assert check_value(1.5, 1.0, tol) is None


def test_only_the_set_threshold_is_checked():
    assert check_value(1.5, 1.0, Tolerance(abs_tol=1.0, rel_tol=None)) is None
    assert check_value(1.5, 1.0, Tolerance(abs_tol=None, rel_tol=1e-12)) is not None


# --- NaN -------------------------------------------------------------------------------

@pytest.mark.parametrize("test,benchmark", [
    (math.nan, 1.0), (1.0, math.nan), (math.nan, math.nan),
])
def test_nan_always_fails(test, benchmark):
    deviation = check_value(test, benchmark, Tolerance(abs_tol=1e9, rel_tol=None))
    assert deviation is not None
    assert "NaN" in deviation.reason


# --- non-numeric values ----------------------------------------------------------------

def test_equal_strings_pass():
    assert check_value("Gamma", "Gamma", DEFAULT_TOLERANCE) is None


def test_different_strings_fail():
    assert check_value("Gamma", "X", DEFAULT_TOLERANCE) is not None


def test_a_number_against_a_string_fails():
    assert check_value(1.0, "1.0", DEFAULT_TOLERANCE) is not None


# --- structural checks -----------------------------------------------------------------

def test_a_key_only_in_the_benchmark_is_a_failure():
    result = compare({"a": [1.0]}, {"a": [1.0], "b": [2.0]})
    assert not result.ok
    assert any("missing from test output" in e for e in result.structural_errors)


def test_a_key_only_in_the_test_output_is_a_failure():
    result = compare({"a": [1.0], "b": [2.0]}, {"a": [1.0]})
    assert not result.ok
    assert any("missing from benchmark" in e for e in result.structural_errors)


def test_differing_value_counts_are_a_failure_not_a_truncated_comparison():
    result = compare({"a": [1.0, 2.0, 3.0]}, {"a": [1.0, 2.0]})
    assert not result.ok
    assert any("3 value(s)" in e and "2" in e for e in result.structural_errors)
    assert result.compared == 0  # nothing was silently compared


def test_matching_data_passes():
    result = compare({"a": [1.0, 2.0], "b": ["x"]}, {"a": [1.0, 2.0], "b": ["x"]})
    assert result.ok
    assert result.compared == 3


# --- the empty-parse guard (an improvement over testcode) ------------------------------

def test_an_empty_benchmark_dict_fails():
    """testcode reported this as a pass, masking parsers that stopped matching."""
    result = compare({}, {})
    assert not result.ok
    assert any("extracted no data" in e for e in result.structural_errors)


def test_a_benchmark_with_only_empty_lists_fails():
    result = compare({"a": []}, {"a": []})
    assert not result.ok
    assert any("extracted no data" in e for e in result.structural_errors)


# --- failure reporting -----------------------------------------------------------------

def test_every_out_of_tolerance_element_is_reported_not_just_the_first():
    test = {"a": [1.0, 2.0, 3.0]}
    benchmark = {"a": [9.0, 9.0, 9.0]}
    result = compare(test, benchmark)
    assert len(result.deviations) == 3
    assert [d.index for d in result.deviations] == [0, 1, 2]


def test_failure_message_names_key_index_both_values_and_both_errors():
    result = compare({"omegaD": [1.5]}, {"omegaD": [1.0]},
                     tolerances={"omegaD": Tolerance(abs_tol=1e-6, rel_tol=5e-6)})
    message = format_failure(result, work_dir="/tmp/w", output="a.wout", benchmark="b/a.wout")
    assert "omegaD[0]" in message
    assert "1.5" in message and "1.0" in message
    assert "abs err" in message and "rel err" in message
    assert "1.000e-06" in message and "5.000e-06" in message


def test_failure_message_reports_the_work_directory():
    result = compare({"a": [1.0]}, {"a": [2.0]})
    assert "/tmp/somewhere" in format_failure(result, work_dir="/tmp/somewhere")


def test_long_failures_are_capped_with_a_tail_count():
    n = MAX_REPORTED_DEVIATIONS + 17
    result = compare({"a": [float(i) for i in range(n)]},
                     {"a": [float(i) + 1.0 for i in range(n)]})
    message = format_failure(result)
    assert message.count("a[") == MAX_REPORTED_DEVIATIONS
    assert "... and 17 more" in message
