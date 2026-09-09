# Wannier90 test suite

Functional/regression tests for `wannier90.x` and `postw90.x`. Each test runs the code on a
small input and compares selected values from one output file against a stored reference,
within per-quantity tolerances.

The suite is driven by [pytest](https://docs.pytest.org/).

## Requirements

- Python **3.10** or newer
- `pytest` and `PyYAML`; `pytest-xdist` and `pytest-timeout` are optional but recommended:

  ```bash
  pip install -r test-suite/requirements.txt
  ```

- The binaries under test, built from the repository root:

  ```bash
  make wannier post w90chk2chk
  ```

  `w90chk2chk.x` is needed by the `postw90.x` tests, which start from a checkpoint file.

By default the harness looks for `wannier90.x`, `postw90.x` and `w90chk2chk.x` in the
repository root. Point it elsewhere with `--wannier90-exe`, `--postw90-exe` and
`--w90chk2chk-exe`, or the environment variables `W90_EXE`, `POSTW90_EXE` and
`W90CHK2CHK_EXE`. It never falls back to `PATH`, so it cannot accidentally test a
system-installed Wannier90 instead of your build.

## Running the tests

Run everything from the `test-suite` directory:

| Command | What it does |
|---|---|
| `pytest` | run the whole suite, serially |
| `pytest -k testw90_example01` | run tests whose name matches |
| `pytest tests/testw90_example01` | run one test by directory |
| `pytest -m postw90` | run only the `postw90.x` tests |
| `pytest -m "wannier90 or checkpoint"` | combine markers |
| `pytest --nprocs=2` | run each test under `mpirun -np 2` |
| `pytest -n auto` | run tests in parallel across CPUs (needs `pytest-xdist`) |
| `pytest --workdir=/tmp/w90 ` | run in a fixed location and keep the results |
| `pytest -x --ff` | stop at the first failure, failed tests first |
| `ctest --test-dir build -L wannier90` | the same tests through CTest |

Markers are derived from the directory name: `testw90_*` are `wannier90`, `testpostw90_*`
are `postw90`, `checkpoint*` are `checkpoint`, `partest*` are `parallel`. Tests that must
always run serially are additionally marked `serial`.

`--nprocs=N` selects a parallel run; without it, the binary is invoked directly rather than
through `mpirun -np 1`, which is not the same thing. Change the launcher for your MPI with
`--mpi-launcher`, e.g. `--mpi-launcher="srun -n {nprocs}"`.

### Where tests run

Every test runs in its own **work directory outside the source tree**, into which the test's
inputs are copied. Nothing is written back into `tests/`, so `git status` stays clean and
runs cannot contaminate one another.

By default that is a pytest temporary directory; pytest keeps the last three runs under
`/tmp/pytest-of-<user>/`, which is usually what you want for a post-mortem. Use
`--workdir=PATH` to put them somewhere predictable (CI does this so it can upload them).

## Adding a test for `wannier90.x`

1. Create `tests/<name>/`. The name **must** start with `testw90_` so it is grouped
   correctly.
2. Put the inputs in it (`.win`, `.amn`, `.mmn`, …) and `git add` them. Large inputs are
   stored bzip2-compressed and decompressed at run time; see `prepare:` below.
3. Write `tests/<name>/test.yaml`:

   ```yaml
   description: |
     Gallium Arsenide, valence bands
   profile: wannier90_wout
   runs:
     - input: gaas.win
       output: gaas.wout
   ```

   `profile` picks the program, parser and tolerances from `profiles.yaml` (see below).
   `output` is the file the code writes and the harness parses.
4. Generate the reference and inspect it:

   ```bash
   pytest --update-benchmarks -k <name>
   git diff --stat
   ```

5. Run `pytest -k <name>` and confirm it passes. Commit the inputs, the `test.yaml` and
   `benchmark/<output>`.

**There is no second place to register the test.** CMake discovers tests by globbing
`tests/*/test.yaml`, so creating the directory is enough.

**Strong suggestion.** Once the test passes, deliberately change one of the checked values
in the reference file and confirm the test *fails*. A test that cannot fail is worse than no
test. Put the correct value back afterwards.

### The full `test.yaml` schema

Everything except `profile` and `runs` is optional.

```yaml
description: |                     # free text; shown in the file, not by pytest
  What this test covers
profile: wannier90_wout            # required; key into profiles.yaml
runs:                              # required; a list, currently always of length one
  - input: gaas.win
    args: ["-pp"]                  # command-line arguments, default none
    output: gaas.wout              # file the code writes and we parse
    benchmark: benchmark/gaas.wout # default: benchmark/<output>
serial_only: false                 # never run under mpirun, even with --nprocs
min_nprocs: 2                      # skip unless --nprocs is at least this
expect_failure: false              # the run is expected to abort
timeout: 600                       # seconds for one run
tags: []                           # extra pytest markers
depends_on: []                     # tests whose artefacts this one consumes
prepare: []                        # input preparation, see below
status: enabled                    # enabled | skip | xfail
status_reason: ""                  # required when status is not "enabled"
```

`prepare:` replaces the per-test `Makefile`s that used to do this work:

```yaml
prepare:
  - chk_from_bz2: "*.chk.fmt.bz2"  # decompress, w90chk2chk.x -f2u, drop the intermediate
  - bunzip2: "*.mmn.bz2"           # plain decompression
```

A `status: skip` test is skipped by pytest with its reason shown in the summary, and marked
`DISABLED` in CTest, so the two agree.

## Adding a test for `postw90.x`

Follow the steps above, with `testpostw90_` as the prefix and a `postw90_*` profile. A
`postw90.x` run starts from a checkpoint, which needs preparing first.

**Never commit a binary `.chk` file.** The format is compiler- and machine-dependent, so it
would not be usable on another computer. Commit the *formatted and compressed* `.chk.fmt.bz2`
instead, and let `prepare:` convert it at run time.

1. Create `checkpoints/<name>/` and put the Wannier90 inputs in it, with exactly one `.win`.
2. Copy the `Makefile` from a sibling (e.g. `checkpoints/si_geninterp`) and run `make` there.
   That runs `wannier90.x`, converts the checkpoint to formatted form and bzips it.
   (`checkpoints/` still uses Makefiles; only the per-test ones under `tests/` were removed.)
3. In your test directory, symlink the result:

   ```bash
   ln -s ../../checkpoints/<name>/<seedname>.chk.fmt.bz2
   ```

   You will usually need the `.amn`, `.mmn` and similar files too.
4. Declare the preparation in `test.yaml`:

   ```yaml
   prepare:
     - chk_from_bz2: "*.chk.fmt.bz2"
     - bunzip2: "*.mmn.bz2"
   ```

Symlinks are materialised into real files when the work directory is populated, so they may
point outside the test directory.

## Regenerating the reference files

```bash
pytest --update-benchmarks              # regenerate everything
git diff --stat                         # see what moved
git diff -- tests/testw90_example01     # inspect one
```

This runs each test exactly as normal, then copies the produced output over the reference
instead of comparing. It refuses to write a reference when the parser extracts nothing, and
it still checks the return code, so a crashed or unparseable run cannot be promoted. At the
end it reports how many references were updated, unchanged, or failed to run.

It works with any subset: `pytest --update-benchmarks -k testw90_example01`,
`pytest --update-benchmarks -m postw90`.

Two caveats:

- **Regenerate serially.** MPI runs can differ in ordering and in the low-order digits, so
  `--update-benchmarks` together with `--nprocs` is refused unless you also pass
  `--allow-parallel-update`.
- **Expect cosmetic diffs.** A `.wout` contains timings, dates and paths, so a regenerated
  reference will differ in those regions even when nothing physical changed. The parsers
  extract specific keys and ignore the rest, so those diffs are harmless — but do read the
  diff rather than committing it blind.

## Profiles

`profiles.yaml` says, for each class of test, which program to run, which parser reads its
output, and to what tolerance values are compared:

```yaml
wannier90_wout:
  program: wannier90            # wannier90 | postw90
  parser: parse_wout            # module in w90testlib/parsers/
  tolerances:
    final_spreads: {abs: 3.0e-6, rel: 3.0e-6}
    omegaD:        {abs: 1.0e-6, rel: 5.0e-6}
```

A value passes when its error is **strictly less** than the tolerance. The relative error is
normalised by the benchmark value. With both tolerances set, both must be satisfied.

A key not listed is **still checked**, against a default of `abs: 1.0e-10` with no relative
check — so listing a key can only loosen the comparison, never tighten it. Listing every key
your parser returns is good practice anyway, because it documents what the test is for.

`rel: null` switches the relative check off for that key, which is what you want for
quantities that pass through zero. Note that tolerances must be written so YAML reads them
as numbers: `1.0e-6`, not `1e-6` — the latter parses as a *string* and would be rejected.

To parse a file no existing parser handles, add `w90testlib/parsers/parse_<something>.py`
exporting

```python
def parse(filename) -> dict[str, list]:
    ...
```

returning `{'key': [value, ...]}` for the values to check, then reference it from a new
profile. Take inspiration from the existing parsers, including how they handle verbose
output (`W90VERBOSETESTS=true`).

Choosing tolerances: start from a similar existing profile. Values that are large in
magnitude usually need a loose absolute tolerance and a tight relative one; values near zero
need the opposite, and often `rel: null`.

## Debugging a failure

The assertion message names the test, every value that fell outside tolerance — with its
key, index, both values, and both errors against both thresholds — and the **work
directory**.

```bash
pytest -k testw90_example01 --workdir=/tmp/w90   # keep the run around
cd /tmp/w90/testw90_example01
cat stderr.log            # why the run died, if it did
cat stdout.log
diff gaas.wout <path-to>/tests/testw90_example01/benchmark/gaas.wout
```

Useful flags: `-x` stops at the first failure, `-ra` (on by default) lists skip reasons,
`-v` shows each test name, `--tb=short` shortens tracebacks, and `-n auto` speeds up a full
run once you are not debugging.

If a test fails with "the parser extracted nothing", the run probably produced no useful
output, or an output-format change stopped the parser matching. That case is deliberately a
failure rather than a pass.

## Running the library-mode tests

`library-mode-test*/` are self-contained CMake projects, not part of the pytest suite. They
are built and run by CTest as part of a normal `cmake`/`ctest` cycle.

## Acknowledgements

We acknowledge S. Poncé for the first implementation of the test-suite in Wannier90.
