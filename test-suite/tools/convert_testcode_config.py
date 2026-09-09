#!/usr/bin/env python3
"""One-shot converter: tests/jobconfig + tests/userconfig -> profiles.yaml + test.yaml.

Run once, from ``test-suite/``:

    python3 tools/convert_testcode_config.py [--apply]

Without ``--apply`` it only reports what it would do.  With it, it writes ``profiles.yaml``
and one ``tests/<name>/test.yaml`` per test, and ``git mv``s each benchmark file into
``tests/<name>/benchmark/<output filename>``.

The script asserts loudly on anything it does not understand rather than guessing: an
unknown option key, a test with more than one (input, args) pair, a profile carrying an
option other than exe/extract_fn/tolerance/can_fail, a missing benchmark, or a benchmark
filename that does not match the expected pattern.

This file is deleted once the migration lands; it is kept in history for reference.
"""

import argparse
import ast
import configparser
import re
import subprocess
import sys
from pathlib import Path

TEST_SUITE = Path(__file__).resolve().parents[1]
TESTS = TEST_SUITE / "tests"
JOBCONFIG = TESTS / "jobconfig"
USERCONFIG = TESTS / "userconfig"
PROFILES_YAML = TEST_SUITE / "profiles.yaml"

# --- what we accept ---------------------------------------------------------------------

JOBCONFIG_KEYS = {"program", "inputs_args", "output", "max_nprocs", "min_nprocs"}
USERCONFIG_KEYS = {"exe", "extract_fn", "tolerance", "can_fail"}

EXE_TO_PROGRAM = {"../../wannier90.x": "wannier90", "../../postw90.x": "postw90"}

BENCHMARK_RE = re.compile(r"^benchmark\.out\.default\.inp=(?P<input>[^.]+\.win)"
                          r"(?:\.args=(?P<args>.+))?$")

# Profile names are lowercased with the redundant _OK / _FAIL suffix dropped (expected
# failure is a per-test property now).  These few also get word breaks inserted, because
# "postw90_shckpathbandsdat" is unreadable.
PROFILE_NAME_OVERRIDES = {
    "POSTW90_GENINTERPDAT_OK": "postw90_geninterp_dat",
    "POSTW90_MORBDAT_OK": "postw90_morb_dat",
    "POSTW90_CURVDAT_OK": "postw90_curv_dat",
    "POSTW90_SHCFERMIDAT_OK": "postw90_shc_fermi_dat",
    "POSTW90_SHCFREQDAT_OK": "postw90_shc_freq_dat",
    "POSTW90_SHCKPATHBANDSDAT_OK": "postw90_shc_kpath_bands_dat",
    "POSTW90_SHCKPATHDAT_OK": "postw90_shc_kpath_dat",
    "POSTW90_SHCKSLICEDAT_OK": "postw90_shc_kslice_dat",
}

# Compressed input families, from the wildcard rules in the per-test Makefiles.
PLAIN_BUNZIP2_EXTENSIONS = ("mmn", "uHu", "amn", "spn", "sHu", "sIu")

# Carried over verbatim from the comments in the old tests/CMakeLists.txt, where these were
# commented out of the test list.
DISABLED = {
    "partestw90_mpierr": "test is intended to fail; disabled in CI",
    "testw90_example11_2": ("CI's old/unfortunate intel/flexiblas versions give too large "
                            "numerical fuzz"),
    "testw90_na_chain_gamma": "flaky",
}

# The one genuine inter-test dependency: checkpoint1_read consumes a .chk that
# checkpoint0_write produces at run time (its Makefile symlinked the sibling directory).
DEPENDENCIES = {
    "checkpoint1_read": {"test": "checkpoint0_write", "files": ["copper.chk"]},
}


class ConversionError(Exception):
    """Raised on anything the converter does not understand."""


def profile_name(section: str) -> str:
    if section in PROFILE_NAME_OVERRIDES:
        return PROFILE_NAME_OVERRIDES[section]
    name = section.lower()
    for suffix in ("_ok", "_fail"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def eval_nested_tuple(raw: str):
    """testcode's tolerance syntax: a tuple of tuples, or a single bare tuple."""
    value = ast.literal_eval(raw)
    if isinstance(value[0], (list, tuple)):
        return value
    return ast.literal_eval(f"{raw},")


def harvest_descriptions(path: Path) -> dict[str, str]:
    """Collect the comment block immediately preceding each [section].

    These comments are the only documentation of what each test does, so they become the
    `description` field rather than being lost.
    """
    descriptions: dict[str, str] = {}
    pending: list[str] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            pending.append(stripped.lstrip("#").strip())
        elif stripped.startswith("["):
            descriptions[stripped[1:-1].rstrip("/")] = " ".join(pending).strip()
            pending = []
        elif not stripped:
            pending = []  # a blank line ends the block
    return descriptions


# --- profiles ---------------------------------------------------------------------------

def convert_profiles(userconfig: configparser.RawConfigParser) -> tuple[dict, dict[str, str]]:
    profiles: dict = {}
    mapping: dict[str, str] = {}

    for section in userconfig.sections():
        if section == "user":
            continue
        unknown = set(userconfig.options(section)) - USERCONFIG_KEYS
        if unknown:
            raise ConversionError(
                f"profile [{section}] has unexpected option(s): {', '.join(sorted(unknown))}")

        exe = userconfig.get(section, "exe")
        if exe not in EXE_TO_PROGRAM:
            raise ConversionError(f"profile [{section}]: unrecognised exe {exe!r}")

        name = profile_name(section)
        if name in profiles:
            raise ConversionError(f"profile name collision: {section} -> {name}")
        mapping[section] = name

        body: dict = {"program": EXE_TO_PROGRAM[exe]}

        if userconfig.has_option(section, "extract_fn"):
            extract = userconfig.get(section, "extract_fn").split()
            # "tools parsers.parse_wout.parse" -> module parse_wout, function parse
            module_path, function = extract[-1].rsplit(".", 1)
            body["parser"] = module_path.split(".")[-1]
            if function != "parse":
                body["function"] = function
        else:
            # POSTW90_HIGHER_OK is a stub with no parser and no test using it.  Keep it,
            # flagged, rather than silently dropping a profile the maintainers wrote.
            body["parser"] = None
            body["comment"] = ("TODO(review): the old userconfig said "
                               "'To do: add extract_fn, tolerance'. No test uses this "
                               "profile; it cannot check anything until a parser is set.")

        if userconfig.has_option(section, "tolerance"):
            tolerances = {}
            for item in eval_nested_tuple(userconfig.get(section, "tolerance")):
                if len(item) < 3:
                    raise ConversionError(
                        f"profile [{section}]: tolerance {item!r} has no key name")
                if len(item) > 4:
                    raise ConversionError(
                        f"profile [{section}]: tolerance {item!r} has too many fields")
                abs_tol, rel_tol, key = item[0], item[1], item[2]
                entry = {"abs": abs_tol, "rel": rel_tol}
                if len(item) == 4 and not item[3]:
                    entry["strict"] = False
                tolerances[key] = entry
            body["tolerances"] = tolerances

        profiles[name] = body

    return profiles, mapping


# --- tests ------------------------------------------------------------------------------

def find_benchmark(directory: Path, input_file: str, args: list[str], output: str) -> Path:
    already_moved = directory / "benchmark" / output
    if already_moved.is_file() and not any(directory.glob("benchmark.out.*")):
        # A previous --apply run already moved it; the script is re-runnable.
        return already_moved
    candidates = sorted(directory.glob("benchmark.out.*"))
    if len(candidates) != 1:
        raise ConversionError(
            f"{directory.name}: expected exactly one benchmark file, found {len(candidates)}")
    benchmark = candidates[0]
    match = BENCHMARK_RE.match(benchmark.name)
    if not match:
        raise ConversionError(f"{directory.name}: benchmark {benchmark.name!r} does not match "
                              "the expected pattern")
    if match.group("input") != input_file:
        raise ConversionError(
            f"{directory.name}: benchmark names input {match.group('input')!r} but jobconfig "
            f"says {input_file!r}")
    benchmark_args = match.group("args")
    expected_args = args[0] if args else None
    if benchmark_args != expected_args:
        raise ConversionError(
            f"{directory.name}: benchmark names args {benchmark_args!r} but jobconfig says "
            f"{expected_args!r}")
    return benchmark


def derive_prepare(directory: Path) -> list[dict]:
    """Work out the prepare steps from the compressed inputs actually present.

    The Makefiles used `$(wildcard *.chk.fmt.bz2)` etc., so "what is in the directory" is
    exactly what they acted on.  The result is cross-checked against the Makefile below.
    """
    steps: list[dict] = []
    if any(directory.glob("*.chk.fmt.bz2")):
        steps.append({"chk_from_bz2": "*.chk.fmt.bz2"})
    for ext in PLAIN_BUNZIP2_EXTENSIONS:
        if any(directory.glob(f"*.{ext}.bz2")):
            steps.append({"bunzip2": f"*.{ext}.bz2"})
    if directory.name in DEPENDENCIES:
        steps.append({"copy_from_dependency": DEPENDENCIES[directory.name]})
    return steps


def check_prepare_against_makefile(directory: Path, steps: list[dict]) -> None:
    """Confirm the derived steps cover exactly what the Makefile declared."""
    makefile = directory / "Makefile"
    if not makefile.is_file():
        if steps:
            raise ConversionError(
                f"{directory.name}: derived prepare steps but there is no Makefile")
        return
    text = makefile.read_text()
    if "ln -sf" in text:
        # checkpoint1_read's symlink rule; handled by DEPENDENCIES.
        if directory.name not in DEPENDENCIES:
            raise ConversionError(
                f"{directory.name}: Makefile creates a symlink but no dependency is declared")
        return
    declared = set(re.findall(r"%\.?([A-Za-z.]*?)\.bz2", text))
    declared = {d for d in declared if d}
    derived = set()
    for step in steps:
        pattern = next(iter(step.values()))
        if isinstance(pattern, str):
            derived.add(pattern.removeprefix("*.").removesuffix(".bz2"))
    # The Makefiles are generic templates driven by $(wildcard ...), so one may declare a
    # family the directory does not actually contain -- that rule simply never fires.  What
    # must not happen is the reverse: a compressed input present but never prepared.
    if not derived <= declared:
        raise ConversionError(
            f"{directory.name}: directory contains {sorted(derived - declared)} which the "
            f"Makefile never prepares (Makefile declares {sorted(declared)})")


def convert_tests(jobconfig, descriptions, profile_map) -> dict[str, dict]:
    cases: dict[str, dict] = {}
    for section in jobconfig.sections():
        if section == "categories":
            continue
        name = section.rstrip("/")
        directory = TESTS / name
        if not directory.is_dir():
            raise ConversionError(f"jobconfig names [{section}] but {directory} does not exist")

        unknown = set(jobconfig.options(section)) - JOBCONFIG_KEYS
        if unknown:
            raise ConversionError(
                f"test [{section}] has unexpected option(s): {', '.join(sorted(unknown))}")

        program = jobconfig.get(section, "program")
        if program not in profile_map:
            raise ConversionError(f"test [{section}] names unknown program {program!r}")

        inputs_args = eval_nested_tuple(jobconfig.get(section, "inputs_args"))
        if len(inputs_args) != 1:
            raise ConversionError(
                f"test [{section}] has {len(inputs_args)} (input, args) pairs; the benchmark "
                "rename assumes exactly one")
        input_file, raw_args = inputs_args[0]
        args = raw_args.split() if raw_args else []
        output = jobconfig.get(section, "output")

        benchmark = find_benchmark(directory, input_file, args, output)
        prepare = derive_prepare(directory)
        check_prepare_against_makefile(directory, prepare)

        case: dict = {}
        if descriptions.get(name):
            case["description"] = descriptions[name]
        case["profile"] = profile_map[program]
        run: dict = {"input": input_file}
        if args:
            run["args"] = args
        run["output"] = output
        case["runs"] = [run]

        if jobconfig.has_option(section, "max_nprocs"):
            max_nprocs = jobconfig.getint(section, "max_nprocs")
            if max_nprocs != 0:
                raise ConversionError(
                    f"test [{section}]: max_nprocs={max_nprocs}; only 0 (serial only) is handled")
            case["serial_only"] = True
        if jobconfig.has_option(section, "min_nprocs"):
            case["min_nprocs"] = jobconfig.getint(section, "min_nprocs")
        if program == "WANNIER90_WERR_FAIL":
            case["expect_failure"] = True
        if prepare:
            case["prepare"] = prepare
        if name in DEPENDENCIES:
            case["depends_on"] = [DEPENDENCIES[name]["test"]]
        if name in DISABLED:
            case["status"] = "skip"
            case["status_reason"] = DISABLED[name]

        case["_benchmark_source"] = benchmark
        case["_benchmark_target"] = directory / "benchmark" / output
        cases[name] = case
    return cases


# --- emitting ----------------------------------------------------------------------------

def yaml_float(value: float) -> str:
    """Format a float so YAML parses it back as a float.

    YAML 1.1 (what PyYAML implements) only recognises a float when it has a decimal point
    *and*, if an exponent is present, an explicitly signed one.  Plain `1e-06` loads as the
    string "1e-06", which would silently disable every comparison using it.
    """
    text = f"{value:.12g}"
    if "e" in text.lower():
        mantissa, _, exponent = text.lower().partition("e")
        if "." not in mantissa:
            mantissa += ".0"
        sign = "-" if exponent.startswith("-") else "+"
        digits = exponent.lstrip("+-").lstrip("0") or "0"
        return f"{mantissa}e{sign}{digits}"
    return text if "." in text else f"{text}.0"


def yaml_scalar(value) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return yaml_float(value)
    if isinstance(value, int):
        return str(value)
    text = str(value)
    if text != text.strip() or not text or any(c in text for c in ":#{}[],&*?|>'\"%@`"):
        escaped = text.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped}"'
    return text


def dump_profiles(profiles: dict) -> str:
    lines = [
        "# Test profiles: which program runs a test, which parser reads its output, and to",
        "# what tolerance its values are compared.",
        "#",
        "# Migrated from the old tests/userconfig. Tolerances are reproduced exactly; do not",
        "# adjust one to make a test pass. A key not listed here is still checked, against the",
        "# default tolerance of abs=1.0e-10 with no relative check.",
        "#",
        "# 'rel: null' means the relative check is switched off for that key (used where values",
        "# pass through zero); it is not the same as omitting the key.",
        "",
    ]
    for name in sorted(profiles):
        body = profiles[name]
        lines.append(f"{name}:")
        if "comment" in body:
            lines.append(f"  # {body['comment']}")
        lines.append(f"  program: {body['program']}")
        lines.append(f"  parser: {yaml_scalar(body.get('parser'))}")
        if "function" in body:
            lines.append(f"  function: {body['function']}")
        if body.get("tolerances"):
            lines.append("  tolerances:")
            for key, entry in body["tolerances"].items():
                fields = ", ".join(f"{k}: {yaml_scalar(v)}" for k, v in entry.items())
                lines.append(f"    {key}: {{{fields}}}")
        lines.append("")
    return "\n".join(lines)


def dump_case(name: str, case: dict) -> str:
    lines = [f"# {name}", ""]
    if "description" in case:
        lines.append(f"description: {yaml_scalar(case['description'])}")
    lines.append(f"profile: {case['profile']}")
    lines.append("runs:")
    for run in case["runs"]:
        lines.append(f"  - input: {run['input']}")
        if "args" in run:
            args = ", ".join(yaml_scalar(a) for a in run["args"])
            lines.append(f"    args: [{args}]")
        lines.append(f"    output: {run['output']}")
    for key in ("serial_only", "min_nprocs", "expect_failure"):
        if key in case:
            lines.append(f"{key}: {yaml_scalar(case[key])}")
    if "depends_on" in case:
        lines.append(f"depends_on: [{', '.join(case['depends_on'])}]")
    if "prepare" in case:
        lines.append("prepare:")
        for step in case["prepare"]:
            (kind, argument), = step.items()
            if isinstance(argument, dict):
                files = ", ".join(yaml_scalar(f) for f in argument["files"])
                lines.append(f"  - {kind}: {{test: {argument['test']}, files: [{files}]}}")
            else:
                lines.append(f"  - {kind}: {yaml_scalar(argument)}")
    if "status" in case:
        lines.append(f"status: {case['status']}")
        lines.append(f"status_reason: {yaml_scalar(case['status_reason'])}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true",
                        help="write the files and perform the git mv (default: dry run)")
    args = parser.parse_args()

    for path in (JOBCONFIG, USERCONFIG):
        if not path.is_file():
            raise ConversionError(f"{path} not found; has the migration already run?")

    userconfig = configparser.RawConfigParser()
    userconfig.optionxform = str
    userconfig.read(USERCONFIG)
    jobconfig = configparser.RawConfigParser()
    jobconfig.optionxform = str
    jobconfig.read(JOBCONFIG)

    profiles, profile_map = convert_profiles(userconfig)
    descriptions = harvest_descriptions(JOBCONFIG)
    cases = convert_tests(jobconfig, descriptions, profile_map)

    print(f"profiles: {len(profiles)}")
    print(f"tests:    {len(cases)}")
    missing = [n for n, c in cases.items() if "description" not in c]
    print(f"tests without a description comment: {len(missing)} {missing or ''}")

    if not args.apply:
        print("\ndry run; pass --apply to write profiles.yaml, the test.yaml files and "
              "git mv the benchmarks")
        return 0

    PROFILES_YAML.write_text(dump_profiles(profiles))
    print(f"wrote {PROFILES_YAML.relative_to(TEST_SUITE)}")

    for name, case in sorted(cases.items()):
        (TESTS / name / "test.yaml").write_text(dump_case(name, case))
        source = case["_benchmark_source"]
        target = case["_benchmark_target"]
        target.parent.mkdir(exist_ok=True)
        if source == target or (target.is_file() and not source.exists()):
            continue  # already moved on a previous run
        subprocess.run(["git", "mv", str(source), str(target)], cwd=TEST_SUITE, check=True)
    print(f"wrote {len(cases)} test.yaml files and moved {len(cases)} benchmarks")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ConversionError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
