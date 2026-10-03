# LayNii test suite

The tests run the compiled LayNii programs on the example data in
`test_data/` and check their outputs. They run on every pull request via
GitHub Actions (`.github/workflows/tests.yml`).

## Running the tests

```bash
make -j4 all wip                        # build released and work in progress programs
pip install -r tests/requirements.txt
python -m pytest tests -n auto          # or: make tests
```

Useful options:

```bash
python -m pytest tests -k LN2_LAYERS    # only tests matching a name
python -m pytest tests -m "not slow"    # skip the slowest cases
python -m pytest tests -m "not wip"     # skip work in progress programs
LAYNII_BIN=/path/to/build python -m pytest tests   # binaries from elsewhere
```

## What is tested

| File | Content |
| --- | --- |
| `test_cli.py` | Every program in the Makefile is built, prints its help, fails cleanly on a missing input file and on unknown options. |
| `test_programs.py` | Every program runs on the example data. All expected outputs must be written and must be readable niftis with the right shape, without NaN/Inf and not empty. |
| `test_correctness.py` | Numerical checks against numpy/scipy and against synthetic inputs with a known answer (e.g. layers of a flat cortical slab). |
| `compare_outputs.py` | Not a test: compares the outputs of two test runs, see below. |

The program lists are read from the Makefile, so a new program added to the
Makefile fails `test_every_program_has_a_case` until a case is added to
`CASES` in `test_programs.py`.

## Known bugs

Known bugs are marked as `xfail(strict=True)` with the reason. Such a test
is expected to fail; once the bug is fixed the test passes, pytest reports it
as `XPASS` and the run fails, as a reminder to remove the marker.

Out-of-bounds memory accesses are marked with `memory_bug=` in
`test_programs.py`. They are only reliably detected by the sanitizer build,
where those cases are expected to fail (`LAYNII_SANITIZED=1`).

## Continuous integration

`.github/workflows/tests.yml` has three jobs:

- **test**: build and run the suite on Ubuntu (g++ and clang++) and macOS.
- **sanitizers**: build with `-fsanitize=address,undefined` and run the suite.
- **regression** (pull requests only): build both the pull request and its
  base branch, run the tests with each set of binaries
  (`LAYNII_SAVE_OUTPUTS=<dir>` stores all outputs) and compare them with
  `compare_outputs.py`. The job fails when outputs change or disappear. If the
  change is intended, add the `output-change` label to the pull request.

To run the comparison locally:

```bash
LAYNII_BIN=/path/to/old/build LAYNII_SAVE_OUTPUTS=/tmp/old python -m pytest tests -n auto
LAYNII_SAVE_OUTPUTS=/tmp/new python -m pytest tests -n auto
python tests/compare_outputs.py /tmp/old /tmp/new
```
