"""The regression suite: one generic test, parametrised over every tests/*/test.yaml.

Discovery, options and fixtures live in ../conftest.py.
"""

import shutil

import pytest

from conftest import UPDATED_BENCHMARKS
from w90testlib.compare import compare, format_failure
from w90testlib.execute import execute_case


def _describe_run_failure(case, result, expected_nonzero: bool = False) -> str:
    what = ("was expected to fail but exited 0" if expected_nonzero
            else f"exited with code {result.returncode}")
    lines = [f"{case.name}: {result.command}", f"  {what}"]
    if result.timed_out:
        lines.append(f"  timed out after {result.timeout}s and was killed")
    tail = result.stderr_tail()
    if tail:
        lines.append("  last lines of stderr.log:")
        lines.extend(f"    {line}" for line in tail.splitlines())
    lines.append(f"  work dir: {result.work_dir}")
    return "\n".join(lines)


def test_regression(case, profiles, executables, run_options, work_dir,
                    dependency_artifacts, pytestconfig):
    skip_reason = case.skip_reason_for(run_options.nprocs)
    if skip_reason:
        pytest.skip(skip_reason)

    profile = profiles[case.profile]
    updating = pytestconfig.getoption("--update-benchmarks")

    results = execute_case(
        case, profile, executables, work_dir, run_options,
        resolve_dependency=dependency_artifacts,
    )

    for spec, result in zip(case.runs, results):
        # An `expect_failure` test is supposed to abort. TODO(review): the migration brief
        # asked for a strictly stronger rule than testcode's `can_fail` -- treat a zero exit
        # code as a failure, on the grounds that the code was meant to abort and did not.
        # That rule is not implementable as written: serial wannier90.x ends a fatal error
        # with a bare Fortran `stop` (wannier_prog.F90), which exits 0. Only the MPI path
        # aborts nonzero, via MPI_Abort. Asserting a nonzero exit would make nnkpt4 and
        # nnkpt5 fail permanently on every serial build.
        #
        # The intent is preserved by a more reliable signal: a run that aborted writes its
        # <seedname>.werr, and a run that wrongly succeeded writes a .wout instead. The
        # output-existence check below therefore catches "was supposed to abort and did
        # not", and the comparison then checks the error message itself.
        if not case.expect_failure:
            assert result.returncode == 0, _describe_run_failure(case, result, False)

        output = work_dir / spec.output
        missing_hint = (" -- an expect_failure test that did not abort writes its .wout "
                        "instead" if case.expect_failure else "")
        assert output.is_file(), (
            f"{case.name}: the run produced no {spec.output!r}{missing_hint}\n"
            f"  command:  {result.command}\n"
            f"  exit code: {result.returncode}\n"
            f"  stderr:   {result.stderr_log}\n"
            f"  work dir: {work_dir}"
        )

        parse = profile.load_parser()
        parsed = parse(str(output))

        # Refuse to promote an unparseable output to a reference, and refuse to let an
        # empty parse count as a pass.
        assert parsed and any(len(v) for v in parsed.values()), (
            f"{case.name}: parser {profile.parser}.{profile.function} extracted nothing "
            f"from {spec.output}\n"
            f"  either the run produced no useful output or the parser no longer matches "
            f"its format\n"
            f"  work dir: {work_dir}"
        )

        benchmark_path = case.benchmark_path(spec)

        if updating:
            benchmark_path.parent.mkdir(parents=True, exist_ok=True)
            changed = (not benchmark_path.is_file()
                       or benchmark_path.read_bytes() != output.read_bytes())
            shutil.copyfile(output, benchmark_path)
            pytestconfig.stash.setdefault(UPDATED_BENCHMARKS, {})[case.name] = changed
            continue

        assert benchmark_path.is_file(), (
            f"{case.name}: no reference file at {benchmark_path}\n"
            f"  generate one with: pytest --update-benchmarks -k {case.name}"
        )
        expected = parse(str(benchmark_path))
        result_of_comparison = compare(
            parsed, expected,
            tolerances=profile.tolerances,
            default_tolerance=profile.default_tolerance,
        )
        assert result_of_comparison.ok, format_failure(
            result_of_comparison,
            work_dir=work_dir, output=output, benchmark=benchmark_path,
        )
