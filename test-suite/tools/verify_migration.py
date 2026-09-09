#!/usr/bin/env python3
"""Prove the pytest configuration means exactly what the testcode configuration meant.

For all 97 tests this derives, independently and from scratch:

    (program, parser module, parser function, input, args, output, benchmark bytes,
     fully-resolved tolerance for every key the parser emits, serial_only, min_nprocs,
     expect_failure)

once from the old tests/jobconfig + tests/userconfig, and once from the new
tests/*/test.yaml + profiles.yaml, then diffs them.  It must report zero differences.

The old side reimplements testcode's resolution rules directly, *including* the regex
fallback in testcode2.validation.compare_data (an unlisted key borrows the tolerance of the
first profile entry whose name regex-matches it) so that any place where that fires shows up
as a difference rather than being assumed away.

Run from test-suite/:  python3 tools/verify_migration.py
Deleted once the migration lands.
"""

import ast
import configparser
import importlib
import re
import subprocess
import sys
from pathlib import Path

TEST_SUITE = Path(__file__).resolve().parents[1]
TESTS = TEST_SUITE / "tests"
sys.path.insert(0, str(TEST_SUITE))

import yaml  # noqa: E402

from w90testlib.case import load_case  # noqa: E402
from w90testlib.profiles import load_profiles  # noqa: E402

DEFAULT_TOL = (1.0e-10, None, True)
EXE_TO_PROGRAM = {"../../wannier90.x": "wannier90", "../../postw90.x": "postw90"}


def eval_nested_tuple(raw):
    value = ast.literal_eval(raw)
    if isinstance(value[0], (list, tuple)):
        return value
    return ast.literal_eval(f"{raw},")


def old_side():
    userconfig = configparser.RawConfigParser(); userconfig.optionxform = str
    userconfig.read(TESTS / "userconfig")
    jobconfig = configparser.RawConfigParser(); jobconfig.optionxform = str
    jobconfig.read(TESTS / "jobconfig")

    profiles = {}
    for section in userconfig.sections():
        if section == "user":
            continue
        entry = {"program": EXE_TO_PROGRAM[userconfig.get(section, "exe")],
                 "parser": None, "function": None, "tolerances": {},
                 "can_fail": userconfig.getboolean(section, "can_fail", fallback=False)}
        if userconfig.has_option(section, "extract_fn"):
            module_path, function = userconfig.get(section, "extract_fn").split()[-1].rsplit(".", 1)
            entry["parser"] = module_path.split(".")[-1]
            entry["function"] = function
        if userconfig.has_option(section, "tolerance"):
            for item in eval_nested_tuple(userconfig.get(section, "tolerance")):
                strict = item[3] if len(item) >= 4 else True
                entry["tolerances"][item[2]] = (item[0], item[1], strict)
        profiles[section] = entry

    cases = {}
    for section in jobconfig.sections():
        if section == "categories":
            continue
        name = section.rstrip("/")
        program = jobconfig.get(section, "program")
        input_file, raw_args = eval_nested_tuple(jobconfig.get(section, "inputs_args"))[0]
        args = raw_args.split() if raw_args else []
        benchmark = "benchmark.out.default.inp=" + input_file
        if args:
            benchmark += ".args=" + " ".join(args)
        cases[name] = {
            "program": profiles[program]["program"],
            "parser": profiles[program]["parser"],
            "function": profiles[program]["function"],
            "input": input_file,
            "args": args,
            "output": jobconfig.get(section, "output"),
            "benchmark_old_path": f"tests/{name}/{benchmark}",
            "serial_only": jobconfig.getint(section, "max_nprocs", fallback=None) == 0,
            "min_nprocs": jobconfig.getint(section, "min_nprocs", fallback=None),
            "expect_failure": profiles[program]["can_fail"],
            "_profile": profiles[program],
        }
    return cases


def resolve_old_tolerance(profile, key):
    """testcode's rule, regex fallback included."""
    if key in profile["tolerances"]:
        return profile["tolerances"][key]
    matches = [tol for name, tol in profile["tolerances"].items()
               if name and re.match(name, key)]
    if matches:
        return matches[0]
    return DEFAULT_TOL


def new_side():
    profiles = load_profiles(TEST_SUITE / "profiles.yaml")
    cases = {}
    for path in sorted(TESTS.glob("*/test.yaml")):
        case = load_case(path.parent)
        profile = profiles[case.profile]
        run = case.runs[0]
        cases[case.name] = {
            "program": profile.program,
            "parser": profile.parser,
            "function": profile.function,
            "input": run.input,
            "args": list(run.args),
            "output": run.output,
            "benchmark_new_path": f"tests/{case.name}/{run.benchmark}",
            "serial_only": case.serial_only,
            "min_nprocs": case.min_nprocs,
            "expect_failure": case.expect_failure,
            "_profile": profile,
        }
    return cases


def main():
    old, new = old_side(), new_side()
    problems = []

    if set(old) != set(new):
        problems.append(f"test sets differ: only old={sorted(set(old)-set(new))} "
                        f"only new={sorted(set(new)-set(old))}")

    checked_keys = 0
    for name in sorted(set(old) & set(new)):
        o, n = old[name], new[name]
        for field in ("program", "parser", "function", "input", "args", "output",
                      "serial_only", "min_nprocs", "expect_failure"):
            if o[field] != n[field]:
                problems.append(f"{name}: {field}: old={o[field]!r} new={n[field]!r}")

        # The reference file must be byte-identical to the one it was renamed from.
        old_bytes = subprocess.run(
            ["git", "show", f"HEAD:test-suite/{o['benchmark_old_path']}"],
            cwd=TEST_SUITE.parent, capture_output=True, check=False).stdout
        new_path = TEST_SUITE / n["benchmark_new_path"]
        if not old_bytes:
            problems.append(f"{name}: could not read the pre-rename benchmark from git")
        elif not new_path.is_file():
            problems.append(f"{name}: renamed benchmark missing at {n['benchmark_new_path']}")
        elif new_path.read_bytes() != old_bytes:
            problems.append(f"{name}: benchmark content changed during the rename")

        # Resolved tolerance for every key this test's parser actually emits.
        if n["parser"] is None:
            continue
        module = importlib.import_module(f"w90testlib.parsers.{n['parser']}")
        parsed = getattr(module, n["function"])(str(new_path))
        for key in parsed:
            checked_keys += 1
            old_tol = resolve_old_tolerance(o["_profile"], key)
            new_tol_obj = n["_profile"].tolerance_for(key)
            new_tol = (new_tol_obj.abs_tol, new_tol_obj.rel_tol, new_tol_obj.strict)
            if old_tol != new_tol:
                problems.append(
                    f"{name}: tolerance for {key!r}: old={old_tol} new={new_tol}")

    print(f"tests compared:            {len(set(old) & set(new))}")
    print(f"tolerance lookups checked: {checked_keys}")
    print(f"differences:               {len(problems)}")
    for problem in problems:
        print(f"  {problem}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
