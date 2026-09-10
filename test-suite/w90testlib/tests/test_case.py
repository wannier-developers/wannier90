"""Unit tests for loading and validating test.yaml."""

import pytest

from w90testlib.case import DEFAULT_TIMEOUT, CaseError, load_case

MINIMAL = """
profile: wannier90_wout
runs:
  - input: gaas.win
    output: gaas.wout
"""


def write_case(tmp_path, name, body):
    directory = tmp_path / name
    directory.mkdir()
    (directory / "test.yaml").write_text(body)
    return directory


def test_minimal_case_loads_with_sensible_defaults(tmp_path):
    case = load_case(write_case(tmp_path, "testw90_example01", MINIMAL))
    assert case.name == "testw90_example01"
    assert case.profile == "wannier90_wout"
    assert len(case.runs) == 1
    assert case.runs[0].args == ()
    assert case.serial_only is False
    assert case.expect_failure is False
    assert case.timeout == DEFAULT_TIMEOUT
    assert case.status == "enabled"


def test_benchmark_defaults_to_benchmark_slash_output(tmp_path):
    case = load_case(write_case(tmp_path, "testw90_example01", MINIMAL))
    assert case.runs[0].benchmark == "benchmark/gaas.wout"


def test_explicit_benchmark_overrides_the_default(tmp_path):
    body = MINIMAL + "    benchmark: benchmark/other.wout\n"
    case = load_case(write_case(tmp_path, "t", body))
    assert case.runs[0].benchmark == "benchmark/other.wout"


def test_args_are_kept_as_a_list(tmp_path):
    body = """
profile: wannier90_nnkp
runs:
  - input: wannier.win
    args: ["-pp"]
    output: wannier.nnkp
"""
    case = load_case(write_case(tmp_path, "testw90_nnkpt2", body))
    assert case.runs[0].args == ("-pp",)


def test_missing_profile_is_rejected(tmp_path):
    with pytest.raises(CaseError, match="'profile' is required"):
        load_case(write_case(tmp_path, "t", "runs:\n  - {input: a.win, output: a.wout}\n"))


def test_missing_runs_is_rejected(tmp_path):
    with pytest.raises(CaseError, match="'runs' is required"):
        load_case(write_case(tmp_path, "t", "profile: p\n"))


def test_unknown_top_level_field_is_rejected(tmp_path):
    with pytest.raises(CaseError, match="unknown field"):
        load_case(write_case(tmp_path, "t", MINIMAL + "\nnonsense: 1\n"))


def test_a_non_enabled_status_requires_a_reason(tmp_path):
    with pytest.raises(CaseError, match="status_reason must be non-empty"):
        load_case(write_case(tmp_path, "t", MINIMAL + "\nstatus: skip\n"))


def test_skip_with_a_reason_is_accepted(tmp_path):
    body = MINIMAL + "\nstatus: skip\nstatus_reason: flaky\n"
    case = load_case(write_case(tmp_path, "t", body))
    assert case.skip_reason_for(nprocs=None) == "flaky"


# --- markers ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,expected", [
    ("testw90_example01", "wannier90"),
    ("testpostw90_fe_ahc", "postw90"),
    ("checkpoint0_write", "checkpoint"),
    ("partestw90_mpierr", "parallel"),
])
def test_markers_are_derived_from_the_directory_name(tmp_path, name, expected):
    case = load_case(write_case(tmp_path, name, MINIMAL))
    assert expected in case.markers


def test_testpostw90_is_not_also_marked_wannier90(tmp_path):
    case = load_case(write_case(tmp_path, "testpostw90_fe_ahc", MINIMAL))
    assert "wannier90" not in case.markers


def test_serial_only_adds_the_serial_marker(tmp_path):
    case = load_case(write_case(tmp_path, "testw90_example07", MINIMAL + "\nserial_only: true\n"))
    assert "serial" in case.markers


def test_tags_are_appended_to_the_derived_markers(tmp_path):
    case = load_case(write_case(tmp_path, "testw90_x", MINIMAL + "\ntags: [slow]\n"))
    assert "wannier90" in case.markers and "slow" in case.markers


# --- process counts --------------------------------------------------------------------

def test_serial_only_tests_stay_serial_even_when_nprocs_is_given(tmp_path):
    case = load_case(write_case(tmp_path, "t", MINIMAL + "\nserial_only: true\n"))
    assert case.runs_at(4) is None


def test_ordinary_tests_honour_nprocs(tmp_path):
    case = load_case(write_case(tmp_path, "t", MINIMAL))
    assert case.runs_at(4) == 4
    assert case.runs_at(None) is None


def test_min_nprocs_skips_when_run_with_too_few_processes(tmp_path):
    case = load_case(write_case(tmp_path, "t", MINIMAL + "\nmin_nprocs: 2\n"))
    assert "at least 2" in case.skip_reason_for(nprocs=None)
    assert "at least 2" in case.skip_reason_for(nprocs=1)
    assert case.skip_reason_for(nprocs=2) is None
