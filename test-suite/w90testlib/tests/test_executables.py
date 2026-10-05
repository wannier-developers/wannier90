"""Unit tests for locating the binaries under test."""

import pytest

from w90testlib.executables import MissingExecutableError, resolve


def test_defaults_point_at_the_makefile_build_locations(tmp_path):
    exes = resolve(tmp_path)
    assert exes.wannier90.name == "wannier90.x"
    assert exes.postw90.name == "postw90.x"
    assert exes.w90chk2chk.name == "w90chk2chk.x"
    # ../wannier90.x relative to test-suite/ is the repository root.
    assert exes.wannier90.parent == tmp_path.parent


def test_an_explicit_path_wins(tmp_path):
    explicit = tmp_path / "build" / "wannier90.x"
    assert resolve(tmp_path, wannier90=str(explicit)).wannier90 == explicit.resolve()


def test_resolution_does_not_require_the_file_to_exist(tmp_path):
    """A wannier90-only run must not fail because postw90.x was never built."""
    resolve(tmp_path)  # no exception


def test_require_names_the_missing_binary_and_how_to_build_it(tmp_path):
    exes = resolve(tmp_path)
    with pytest.raises(MissingExecutableError) as excinfo:
        exes.require("postw90")
    message = str(excinfo.value)
    assert "postw90.x not found" in message
    assert "make post" in message
    assert "--postw90-exe" in message


def test_require_returns_the_path_when_the_binary_is_present(tmp_path):
    binary = tmp_path / "wannier90.x"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)
    exes = resolve(tmp_path, wannier90=str(binary))
    assert exes.require("wannier90") == binary.resolve()


def test_a_non_executable_file_is_rejected(tmp_path):
    binary = tmp_path / "wannier90.x"
    binary.write_text("not executable")
    binary.chmod(0o644)
    exes = resolve(tmp_path, wannier90=str(binary))
    with pytest.raises(MissingExecutableError, match="not executable"):
        exes.require("wannier90")
