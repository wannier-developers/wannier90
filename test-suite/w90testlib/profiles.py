"""Profiles: which program runs a test, which parser reads its output, and to what tolerance.

Migrated from the ``[SECTION]`` blocks of the old ``tests/userconfig``.  One profile is
shared by many tests; see ``profiles.yaml``.
"""

import importlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import yaml

from .compare import DEFAULT_TOLERANCE, Tolerance

#: Programs a profile may name, mapped to the pytest option that supplies the binary.
PROGRAMS = ("wannier90", "postw90")

#: Parser modules live here and export a function (``parse`` unless stated otherwise)
#: taking a filename and returning ``dict[str, list]``.
PARSER_PACKAGE = "w90testlib.parsers"

DEFAULT_PARSER_FUNCTION = "parse"


class ProfileError(Exception):
    """Raised when profiles.yaml is malformed or names something that does not exist."""


@dataclass(frozen=True)
class Profile:
    """How to run and check one class of test."""

    name: str
    program: str
    parser: str | None = None
    function: str = DEFAULT_PARSER_FUNCTION
    tolerances: dict[str, Tolerance] = field(default_factory=dict)
    default_tolerance: Tolerance = DEFAULT_TOLERANCE

    def tolerance_for(self, key: str) -> Tolerance:
        """Tolerance for one parsed key.

        Keys the profile does not name fall back to the default (abs 1e-10, no relative
        check).  The old harness additionally tried to *regex*-match unlisted keys against
        the tolerance names and would borrow the first match; that path never fired for any
        key produced by any parser against any current benchmark, so it is deliberately not
        reproduced here.  See README.md.
        """
        return self.tolerances.get(key, self.default_tolerance)

    def load_parser(self) -> Callable[[str], dict[str, list[Any]]]:
        """Import the parser function this profile names."""
        if self.parser is None:
            raise ProfileError(
                f"profile {self.name!r} declares no parser, so its output cannot be checked"
            )
        try:
            module = importlib.import_module(f"{PARSER_PACKAGE}.{self.parser}")
        except ImportError as exc:
            raise ProfileError(
                f"profile {self.name!r}: cannot import parser module "
                f"{PARSER_PACKAGE}.{self.parser}: {exc}"
            ) from exc
        try:
            return getattr(module, self.function)
        except AttributeError as exc:
            raise ProfileError(
                f"profile {self.name!r}: parser module {self.parser} has no function "
                f"{self.function!r}"
            ) from exc


def _tolerance_from_yaml(profile_name: str, key: str, raw: Any) -> Tolerance:
    if not isinstance(raw, dict):
        raise ProfileError(
            f"profile {profile_name!r}, tolerance {key!r}: expected a mapping with "
            f"'abs' and/or 'rel', got {raw!r}"
        )
    unknown = set(raw) - {"abs", "rel", "strict"}
    if unknown:
        raise ProfileError(
            f"profile {profile_name!r}, tolerance {key!r}: unknown field(s) "
            f"{', '.join(sorted(unknown))}"
        )
    # A key present with the value null means "this check is switched off"; a key that is
    # absent entirely means "leave it at the default".  These are different things.
    abs_tol = raw["abs"] if "abs" in raw else DEFAULT_TOLERANCE.abs_tol
    rel_tol = raw["rel"] if "rel" in raw else DEFAULT_TOLERANCE.rel_tol
    for field, value in (("abs", abs_tol), ("rel", rel_tol)):
        # A tolerance written as `1e-6` rather than `1.0e-6` loads as a *string*, because
        # YAML 1.1 wants a decimal point and a signed exponent.  Silently accepting that
        # would disable the check, so refuse it loudly and say how to fix it.
        if value is not None and not isinstance(value, (int, float)):
            raise ProfileError(
                f"profile {profile_name!r}, tolerance {key!r}: {field} is {value!r}, which "
                f"YAML parsed as {type(value).__name__}, not a number. Write it with a "
                f"decimal point and a signed exponent, e.g. 1.0e-6 rather than 1e-6."
            )
    return Tolerance(abs_tol=abs_tol, rel_tol=rel_tol, strict=bool(raw.get("strict", True)))


def profile_from_dict(name: str, raw: dict[str, Any]) -> Profile:
    unknown = set(raw) - {"program", "parser", "function", "tolerances", "comment"}
    if unknown:
        raise ProfileError(f"profile {name!r}: unknown field(s) {', '.join(sorted(unknown))}")

    program = raw.get("program")
    if program not in PROGRAMS:
        raise ProfileError(
            f"profile {name!r}: 'program' must be one of {', '.join(PROGRAMS)}, got {program!r}"
        )

    tolerances = {
        key: _tolerance_from_yaml(name, key, value)
        for key, value in (raw.get("tolerances") or {}).items()
    }
    return Profile(
        name=name,
        program=program,
        parser=raw.get("parser"),
        function=raw.get("function", DEFAULT_PARSER_FUNCTION),
        tolerances=tolerances,
    )


def load_profiles(path: Path) -> dict[str, Profile]:
    """Load every profile from ``profiles.yaml``."""
    if not path.is_file():
        raise ProfileError(f"profiles file not found: {path}")
    raw = yaml.safe_load(path.read_text()) or {}
    if not isinstance(raw, dict):
        raise ProfileError(f"{path}: expected a mapping of profile name -> definition")
    return {name: profile_from_dict(name, body or {}) for name, body in raw.items()}
