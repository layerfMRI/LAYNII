"""Shared helpers for the LayNii test suite.

The test suite treats every LayNii program as a black box: it runs the
compiled binary on the example data in `test_data/` inside a temporary
directory and inspects the files that were written.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import nibabel as nib
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DATA = REPO_ROOT / "test_data"

# Set by the CI job that builds with -fsanitize=address,undefined
SANITIZED = os.environ.get("LAYNII_SANITIZED") == "1"

# Default per-call timeout in seconds. Sanitizer builds are much slower, so
# CI can raise this via the environment.
DEFAULT_TIMEOUT = float(os.environ.get("LAYNII_TEST_TIMEOUT", "300"))


def bin_dir():
    """Directory with the compiled LayNii programs (default: repo root)."""
    return Path(os.environ.get("LAYNII_BIN", REPO_ROOT)).resolve()


def exe_path(program):
    path = bin_dir() / program
    if os.name == "nt":
        path = path.with_suffix(".exe")
    return path


# =============================================================================
# Program lists, parsed from the Makefile so that the tests never go stale
# =============================================================================
def _makefile_variable(text, name):
    match = re.search(r"^%s\s*=(.*?)(?:\n\s*\n|\n(?=\S))" % re.escape(name),
                      text, flags=re.M | re.S)
    if match is None:
        return []
    body = match.group(1).replace("\\", " ")
    return [tok for tok in body.split() if not tok.startswith("$(")]


def makefile_programs():
    """Return (released_programs, work_in_progress_programs)."""
    text = (REPO_ROOT / "Makefile").read_text()
    released = []
    for group in ("LAYNII2", "HIGH_PRIORITY", "LOW_PRIORITY", "DERIVATIVES"):
        released += _makefile_variable(text, group)
    wip = _makefile_variable(text, "WORK_IN_PROGRESS")
    # Preserve order, drop duplicates (e.g. LN2_PHASE_JOLT is listed twice)
    released = list(dict.fromkeys(released))
    wip = [p for p in dict.fromkeys(wip) if p not in released]
    return released, wip


RELEASED_PROGRAMS, WIP_PROGRAMS = makefile_programs()
ALL_PROGRAMS = RELEASED_PROGRAMS + WIP_PROGRAMS


def memory_bug(reason, flaky=False):
    """Marks for a known out-of-bounds access.

    Only the sanitizer build reliably detects these, so the test is expected
    to fail there. Normal builds usually survive the bug, unless it is flaky.
    """
    if SANITIZED:
        return [pytest.mark.xfail(reason=reason, strict=True)]
    if flaky:
        return [pytest.mark.xfail(reason=reason, strict=False)]
    return []


def known_memory_bug(reason, flaky=False):
    """Decorator version of memory_bug() for test functions."""
    def decorate(func):
        for mark in memory_bug(reason, flaky):
            func = mark(func)
        return func
    return decorate


# =============================================================================
# Running programs
# =============================================================================
class RunResult:
    def __init__(self, cmd, returncode, stdout, stderr, workdir, new_files):
        self.cmd = cmd
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.workdir = workdir
        self.new_files = new_files

    @property
    def output(self):
        return self.stdout + self.stderr

    def describe(self, tail=40):
        lines = self.output.splitlines()[-tail:]
        return ("command: %s\nexit code: %s\nworkdir: %s\n--- last %d lines of "
                "output ---\n%s" % (" ".join(self.cmd), self.returncode,
                                     self.workdir, tail, "\n".join(lines)))


def link_inputs(workdir, names, source_dirs):
    """Make input files available in `workdir` (symlink, copy as fallback)."""
    for name in names:
        for src_dir in source_dirs:
            src = Path(src_dir) / name
            if src.exists():
                break
        else:
            raise FileNotFoundError("Test input '%s' not found in %s"
                                    % (name, [str(d) for d in source_dirs]))
        dst = Path(workdir) / name
        if dst.exists():
            continue
        try:
            os.symlink(src, dst)
        except OSError:
            shutil.copy(src, dst)


def run_program(program, args, workdir, timeout=None, check=True,
                stdin=subprocess.DEVNULL):
    """Run a LayNii program in `workdir` and report which files it created."""
    exe = exe_path(program)
    if not exe.exists():
        raise FileNotFoundError(
            "%s not found in %s. Build it first (make all wip) or set "
            "LAYNII_BIN." % (program, bin_dir()))
    workdir = Path(workdir)
    before = {p.name for p in workdir.iterdir()}
    cmd = [str(exe)] + [str(a) for a in args]
    try:
        proc = subprocess.run(cmd, cwd=str(workdir), stdin=stdin,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              timeout=timeout or DEFAULT_TIMEOUT)
    except subprocess.TimeoutExpired as err:
        raise AssertionError("%s timed out after %ss: %s"
                             % (program, err.timeout, " ".join(cmd)))
    after = {p.name for p in workdir.iterdir()}
    result = RunResult(cmd, proc.returncode,
                       proc.stdout.decode(errors="replace"),
                       proc.stderr.decode(errors="replace"),
                       workdir, sorted(after - before))
    if check:
        assert proc.returncode == 0, "%s failed\n%s" % (program,
                                                       result.describe())
    return result


# =============================================================================
# Inspecting outputs
# =============================================================================
def load(path):
    """Load a nifti and return (image, data array with scaling applied)."""
    img = nib.load(str(path))
    return img, np.asanyarray(img.dataobj)


def check_nifti(path, shape=None, disk_dtype=None, allow_nonfinite=False,
                allow_all_zero=False, value_range=None):
    """Assert that `path` is a readable, plausible nifti; return its data."""
    path = Path(path)
    assert path.exists(), "expected output %s was not written" % path.name
    img, data = load(path)
    # Reading the full array catches truncated or inconsistent files
    data = np.asarray(data)
    if shape is not None:
        assert data.shape[:len(shape)] == tuple(shape), (
            "%s: shape %s, expected %s" % (path.name, data.shape, shape))
        extra = data.shape[len(shape):]
        assert all(s == 1 for s in extra), (
            "%s: unexpected extra dimensions %s" % (path.name, data.shape))
    if disk_dtype is not None:
        assert img.get_data_dtype() == np.dtype(disk_dtype), (
            "%s: stored as %s, expected %s"
            % (path.name, img.get_data_dtype(), disk_dtype))
    if data.dtype.kind == "f" and not allow_nonfinite:
        n_bad = int(np.count_nonzero(~np.isfinite(data)))
        assert n_bad == 0, "%s: %d NaN/Inf voxels" % (path.name, n_bad)
    if not allow_all_zero:
        assert np.any(data), "%s: output is all zeros" % path.name
    if value_range is not None:
        lo, hi = value_range
        finite = data[np.isfinite(data)] if data.dtype.kind == "f" else data
        assert finite.min() >= lo and finite.max() <= hi, (
            "%s: values in [%s, %s], expected within [%s, %s]"
            % (path.name, finite.min(), finite.max(), lo, hi))
    return data


def spatial_shape(name, source_dirs=(TEST_DATA,)):
    for d in source_dirs:
        p = Path(d) / name
        if p.exists():
            return nib.load(str(p)).shape[:3]
    raise FileNotFoundError(name)
