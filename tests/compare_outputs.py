#!/usr/bin/env python3
"""Compare LayNii outputs written by two test runs.

Run the test suite twice with LAYNII_SAVE_OUTPUTS pointing to different
directories (e.g. once with binaries built from the base branch and once from
a pull request), then:

    python tests/compare_outputs.py base_outputs/ head_outputs/

Writes a markdown report to stdout (and to $GITHUB_STEP_SUMMARY when set).
Exits with 1 if any output changed or disappeared, or if there is nothing to
compare. New outputs (e.g. from new test cases) are only reported.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import nibabel as nib


# Tests whose outputs differ between identical runs (uninitialized memory)
NONDETERMINISTIC = [
    "test_program_LN3_NOLAD_",
    "test_program_LN2_DIRECTIONALITY_BIN_",
]


def is_nifti(path):
    return path.name.endswith((".nii", ".nii.gz"))


def compare_nifti(a, b, rtol, atol):
    img_a, img_b = nib.load(str(a)), nib.load(str(b))
    da = np.asanyarray(img_a.dataobj)
    db = np.asanyarray(img_b.dataobj)
    if da.shape != db.shape:
        return "shape %s -> %s" % (da.shape, db.shape)
    if img_a.get_data_dtype() != img_b.get_data_dtype():
        return "datatype %s -> %s" % (img_a.get_data_dtype(),
                                      img_b.get_data_dtype())
    if not np.allclose(img_a.affine, img_b.affine, atol=1e-6):
        return "affine changed"
    da = da.astype(np.float64)
    db = db.astype(np.float64)
    nan_a, nan_b = np.isnan(da), np.isnan(db)
    if not np.array_equal(nan_a, nan_b):
        return "NaN voxels %d -> %d" % (nan_a.sum(), nan_b.sum())
    ok = ~nan_a
    close = np.isclose(da[ok], db[ok], rtol=rtol, atol=atol)
    if close.all():
        return None
    diff = np.abs(da[ok] - db[ok])
    return ("%d of %d voxels differ (%.3g%%), max abs diff %.4g"
            % ((~close).sum(), close.size, 100.0 * (~close).mean(),
               diff.max()))


def compare_file(a, b, rtol, atol):
    if a.read_bytes() == b.read_bytes():
        return None
    if is_nifti(a):
        try:
            return compare_nifti(a, b, rtol, atol)
        except Exception as err:  # unreadable output in one of the runs
            return "could not compare: %s" % err
    return "content differs"


def files(root):
    return {str(p.relative_to(root)) for p in Path(root).rglob("*")
            if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("base", type=Path)
    parser.add_argument("head", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--ignore", action="append",
                        default=list(NONDETERMINISTIC),
                        help="test directory prefix to skip (repeatable)")
    parser.add_argument("--no-fail", action="store_true",
                        help="always exit with 0")
    args = parser.parse_args()

    base, head = files(args.base), files(args.head)

    def ignored(rel):
        return any(rel.startswith(prefix) for prefix in args.ignore)

    changed, added, removed = [], [], []
    for rel in sorted(base & head):
        if ignored(rel):
            continue
        msg = compare_file(args.base / rel, args.head / rel, args.rtol,
                           args.atol)
        if msg:
            changed.append((rel, msg))
    added = [r for r in sorted(head - base) if not ignored(r)]
    removed = [r for r in sorted(base - head) if not ignored(r)]

    lines = ["## LayNii output regression check", ""]
    lines.append("Compared %d output files (rtol=%g, atol=%g)."
                 % (len(base & head), args.rtol, args.atol))
    lines.append("")
    common = base & head
    if not common:
        lines.append(":x: No outputs to compare. Did the base run fail?")
        lines.append("")
    elif not (changed or removed):
        lines.append("All outputs are identical to the base branch. "
                     ":white_check_mark:")
        lines.append("")
    else:
        lines.append(":warning: Outputs differ from the base branch. If this "
                     "is intended, add the `output-change` label to the pull "
                     "request.")
        lines.append("")
    if changed:
        lines += ["### Changed", "", "| Output | Difference |",
                  "| --- | --- |"]
        lines += ["| `%s` | %s |" % (rel, msg) for rel, msg in changed]
        lines.append("")
    if removed:
        lines += ["### Missing outputs", ""]
        lines += ["- `%s`" % rel for rel in removed]
        lines.append("")
    if added:
        lines += ["<details><summary>%d new outputs (not in base)</summary>"
                  % len(added), ""]
        lines += ["- `%s`" % rel for rel in added]
        lines += ["", "</details>", ""]
    report = "\n".join(lines) + "\n"
    sys.stdout.write(report)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as fh:
            fh.write(report)

    if (changed or removed or not common) and not args.no_fail:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
