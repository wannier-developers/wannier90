"""Unit tests for profiles.yaml loading and tolerance resolution."""

import pytest

from w90testlib.compare import DEFAULT_TOLERANCE
from w90testlib.profiles import ProfileError, load_profiles, profile_from_dict


def write_profiles(tmp_path, body):
    path = tmp_path / "profiles.yaml"
    path.write_text(body)
    return path


def test_a_profile_loads_with_its_tolerances(tmp_path):
    path = write_profiles(tmp_path, """
wannier90_wout:
  program: wannier90
  parser: parse_wout
  tolerances:
    final_spreads: {abs: 3.0e-6, rel: 3.0e-6}
""")
    profile = load_profiles(path)["wannier90_wout"]
    assert profile.program == "wannier90"
    assert profile.parser == "parse_wout"
    assert profile.function == "parse"
    tol = profile.tolerance_for("final_spreads")
    assert tol.abs_tol == pytest.approx(3.0e-6)
    assert tol.rel_tol == pytest.approx(3.0e-6)


def test_unlisted_keys_get_the_default_tolerance(tmp_path):
    path = write_profiles(tmp_path, """
p:
  program: wannier90
  parser: parse_wout
  tolerances:
    a: {abs: 1.0e-5, rel: 1.0e-5}
""")
    assert load_profiles(path)["p"].tolerance_for("not_mentioned") == DEFAULT_TOLERANCE


def test_explicit_null_relative_tolerance_disables_the_relative_check(tmp_path):
    """Distinct from omitting 'rel', which leaves it at the default (also None here)."""
    path = write_profiles(tmp_path, """
p:
  program: postw90
  parser: parse_geninterp_dat
  tolerances:
    near_zero: {abs: 1.0e-8, rel: null}
""")
    tol = load_profiles(path)["p"].tolerance_for("near_zero")
    assert tol.abs_tol == pytest.approx(1.0e-8)
    assert tol.rel_tol is None


def test_omitting_abs_leaves_it_at_the_default(tmp_path):
    path = write_profiles(tmp_path, "p:\n  program: wannier90\n  parser: x\n"
                                    "  tolerances:\n    k: {rel: 1.0e-4}\n")
    tol = load_profiles(path)["p"].tolerance_for("k")
    assert tol.abs_tol == DEFAULT_TOLERANCE.abs_tol
    assert tol.rel_tol == pytest.approx(1.0e-4)


def test_a_non_default_parser_function_is_honoured(tmp_path):
    path = write_profiles(tmp_path, """
postw90_boltzwann_elcond:
  program: postw90
  parser: parse_boltzwann
  function: parse_elcond
""")
    assert load_profiles(path)["postw90_boltzwann_elcond"].function == "parse_elcond"


def test_an_unknown_program_is_rejected():
    with pytest.raises(ProfileError, match="must be one of"):
        profile_from_dict("p", {"program": "quantum_espresso"})


def test_an_unknown_profile_field_is_rejected():
    with pytest.raises(ProfileError, match="unknown field"):
        profile_from_dict("p", {"program": "wannier90", "exe": "../../wannier90.x"})


def test_an_unknown_tolerance_field_is_rejected():
    with pytest.raises(ProfileError, match="unknown field"):
        profile_from_dict("p", {"program": "wannier90",
                                "tolerances": {"k": {"abs": 1e-6, "absolute": 1e-6}}})


def test_a_profile_without_a_parser_explains_itself_when_used():
    profile = profile_from_dict("postw90_higher", {"program": "postw90"})
    with pytest.raises(ProfileError, match="declares no parser"):
        profile.load_parser()


def test_a_missing_parser_module_is_reported_clearly():
    profile = profile_from_dict("p", {"program": "wannier90", "parser": "parse_nonexistent"})
    with pytest.raises(ProfileError, match="cannot import parser module"):
        profile.load_parser()


def test_a_tolerance_yaml_parsed_as_a_string_is_refused(tmp_path):
    """`1e-6` without a decimal point loads as a str in YAML 1.1; that must not pass.

    Silently accepting it would disable the check entirely, which is exactly the kind of
    failure this suite exists to catch.
    """
    path = write_profiles(tmp_path, "p:\n  program: wannier90\n  parser: x\n"
                                    "  tolerances:\n    k: {abs: 1e-6, rel: null}\n")
    with pytest.raises(ProfileError, match="not a number"):
        load_profiles(path)


def test_the_error_says_how_to_write_it_correctly(tmp_path):
    path = write_profiles(tmp_path, "p:\n  program: wannier90\n  parser: x\n"
                                    "  tolerances:\n    k: {abs: 1e-6}\n")
    with pytest.raises(ProfileError, match=r"1\.0e-6"):
        load_profiles(path)
