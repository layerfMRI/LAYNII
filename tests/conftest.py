"""pytest configuration and fixtures for the LayNii test suite."""
import os
import re
import shutil
from pathlib import Path

import numpy as np
import nibabel as nib
import pytest

import laynii_testing as lt


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "wip: program is listed under WORK_IN_PROGRESS in the "
        "Makefile")
    config.addinivalue_line(
        "markers", "slow: takes more than ~10 seconds with an optimized build")


# =============================================================================
# Saving outputs for regression comparisons (see tests/compare_outputs.py)
# =============================================================================
SAVE_DIR = os.environ.get("LAYNII_SAVE_OUTPUTS")


def _save_outputs(nodeid, result):
    if not SAVE_DIR:
        return
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", nodeid.split("::", 1)[-1])
    dst = Path(SAVE_DIR) / name
    dst.mkdir(parents=True, exist_ok=True)
    for fname in result.new_files:
        src = Path(result.workdir) / fname
        if src.is_file():
            shutil.copy(src, dst / fname)


# =============================================================================
# Fixtures
# =============================================================================
class Runner:
    """Runs LayNii programs inside a per-test scratch directory."""

    def __init__(self, workdir, nodeid, sources):
        self.workdir = Path(workdir)
        self.nodeid = nodeid
        self.sources = sources

    def path(self, name):
        return self.workdir / name

    def link(self, *names):
        lt.link_inputs(self.workdir, names, self.sources)

    def save(self, name, data, affine=None, dtype=None):
        """Write a synthetic nifti into the scratch directory."""
        data = np.asarray(data if dtype is None else data.astype(dtype))
        img = nib.Nifti1Image(data, np.eye(4) if affine is None else affine)
        img.header.set_data_dtype(data.dtype)
        nib.save(img, str(self.path(name)))
        return self.path(name)

    def run(self, program, *args, inputs=(), check=True, timeout=None):
        """Run `program`; arguments that name test data are linked in."""
        needed = set(inputs)
        for arg in args:
            arg = str(arg)
            if (not arg.startswith("-") and not self.path(arg).exists()
                    and any((Path(d) / arg).is_file() for d in self.sources)):
                needed.add(arg)
        self.link(*sorted(needed))
        result = lt.run_program(program, args, self.workdir, timeout=timeout,
                                check=check)
        _save_outputs(self.nodeid, result)
        return result

    def load(self, name):
        return lt.load(self.path(name))[1]


@pytest.fixture
def laynii(tmp_path, request, derived):
    return Runner(tmp_path, request.node.nodeid, [lt.TEST_DATA, derived])


@pytest.fixture(scope="session")
def derived(tmp_path_factory):
    """Inputs that are derived from test_data once per session.

    Some programs need outputs of other programs (e.g. UV coordinates from
    LN2_MULTILATERATE) or smaller inputs to run in a reasonable time.
    """
    out = tmp_path_factory.mktemp("derived")
    lt.link_inputs(out, [p.name for p in lt.TEST_DATA.glob("*.nii*")],
                   [lt.TEST_DATA])

    def crop(src, dst, slices):
        img = nib.load(str(lt.TEST_DATA / src))
        data = np.asanyarray(img.dataobj)[slices]
        new = nib.Nifti1Image(data, img.affine, img.header)
        new.set_data_dtype(img.get_data_dtype())
        nib.save(new, str(out / dst))

    # Small cut-outs of the columnar test data (LN_COLUMNAR_DIST is slow)
    box = (slice(265, 365), slice(265, 365), slice(6, 9))
    crop("sc_layers_3dcolumns.nii.gz", "crop_layers.nii.gz", box)
    crop("sc_landmarks.nii.gz", "crop_landmarks.nii.gz", box)
    box = (slice(270, 365), slice(265, 365), slice(6, 9))
    crop("sc_layers_3dcolumns.nii.gz", "crop_layers_rect.nii.gz", box)
    crop("sc_landmarks.nii.gz", "crop_landmarks_rect.nii.gz", box)

    # Single 2D slices for the 2D programs
    plane = (slice(None), slice(None), slice(7, 8))
    for name in ("sc_rim", "sc_midGM", "sc_VASO_act"):
        crop(name + ".nii.gz", "slice_%s.nii.gz" % name, plane)

    # Cortical depth, UV coordinates and a flat patch of the occipital data
    lt.run_program("LN2_LAYERS", ["-rim", "Ding2016_occip_rim.nii.gz",
                                  "-nr_layers", "3", "-output", "ding"], out)
    lt.run_program("LN2_MULTILATERATE", [
        "-rim", "Ding2016_occip_rim.nii.gz",
        "-control_points",
        "Ding2016_occipital_rim_midGM_equidist_control_point0.nii.gz",
        "-radius", "10", "-output", "ding"], out)
    lt.run_program("LN2_PATCH_FLATTEN", [
        "-values", "Ding2016_occip_T2starweighted_filtered_for_tests.nii.gz",
        "-coord_uv", "ding_UV_coordinates.nii",
        "-coord_d", "ding_metric_equidist.nii",
        "-domain", "ding_perimeter_chunk.nii",
        "-bins_u", "30", "-bins_v", "30", "-bins_d", "5", "-voronoi",
        "-output", "ding"], out)
    lt.run_program("LN2_LAYERS", ["-rim", "slice_sc_rim.nii.gz",
                                  "-nr_layers", "3", "-output", "slice"], out)
    lt.run_program("LN2_GEODISTANCE", ["-domain", "slice_sc_rim.nii.gz",
                                       "-init", "slice_sc_midGM.nii.gz",
                                       "-output", "slice_geodist"], out)
    lt.run_program("LN2_IFPOINTS", ["-domain", "sc_midGM.nii.gz",
                                    "-nr_points", "10"], out)

    # Synthetic SIEMENS pulse log: 5 header values, interleaved channel
    # values with trigger markers (5000) and the end-of-data marker (5003)
    values = [1, 2, 40, 280, 5]
    rng = np.random.default_rng(0)
    for i in range(400):
        if i % 50 == 10:
            values.append(5000)
        values += [10240 + int(rng.integers(0, 200)),
                   2048 + int(rng.integers(0, 200))]
    values.append(5003)
    (out / "synthetic.puls").write_text(" ".join(map(str, values)) + "\n")
    return out
