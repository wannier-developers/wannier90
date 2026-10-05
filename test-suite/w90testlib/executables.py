"""Locating the binaries under test.

Never fall back to ``PATH``: that risks silently testing a Wannier90 installed system-wide
instead of the one just built, which is the most confusing possible outcome for a
regression suite.
"""

import os
from dataclasses import dataclass
from pathlib import Path

#: Default locations, relative to ``test-suite/`` -- where a plain ``make`` build leaves them.
DEFAULT_RELATIVE_PATHS = {
    "wannier90": "../wannier90.x",
    "postw90": "../postw90.x",
    "w90chk2chk": "../w90chk2chk.x",
}

#: How to build each one, quoted back at the user when it is missing.
BUILD_HINTS = {
    "wannier90": "make wannier",
    "postw90": "make post",
    "w90chk2chk": "make w90chk2chk",
}


class MissingExecutableError(Exception):
    """Raised when a binary a test needs has not been built."""


@dataclass(frozen=True)
class Executables:
    """Resolved paths to the three binaries the suite drives."""

    wannier90: Path | None = None
    postw90: Path | None = None
    w90chk2chk: Path | None = None

    def require(self, which: str) -> Path:
        """Return the path for ``which``, or explain how to build it."""
        path = getattr(self, which, None)
        if path is None:
            raise MissingExecutableError(
                f"no path configured for {which}.x; pass --{which.replace('_', '-')}-exe"
            )
        if not path.is_file():
            raise MissingExecutableError(
                f"{which}.x not found at {path}\n"
                f"Build it first, e.g. `{BUILD_HINTS.get(which, 'make')}` from the "
                f"repository root, or point at it with "
                f"--{which.replace('_', '-')}-exe=/path/to/{which}.x"
            )
        if not os.access(path, os.X_OK):
            raise MissingExecutableError(f"{path} exists but is not executable")
        return path



def resolve(
    test_suite_root: Path,
    *,
    wannier90: str | None = None,
    postw90: str | None = None,
    w90chk2chk: str | None = None,
) -> Executables:
    """Resolve each binary from an explicit path, else the default build location.

    Paths are resolved but *not* checked for existence here -- a test that never runs
    postw90.x should not fail because postw90.x is missing.  The check happens in
    ``Executables.require`` at the point of use.
    """
    def pick(explicit: str | None, key: str) -> Path:
        if explicit:
            return Path(explicit).expanduser().resolve()
        return (test_suite_root / DEFAULT_RELATIVE_PATHS[key]).resolve()

    return Executables(
        wannier90=pick(wannier90, "wannier90"),
        postw90=pick(postw90, "postw90"),
        w90chk2chk=pick(w90chk2chk, "w90chk2chk"),
    )
