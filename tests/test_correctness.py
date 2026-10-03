"""Numerical checks of LayNii outputs.

Outputs are compared against independent numpy/scipy implementations or
against synthetic inputs with a known answer.
"""
import re

import numpy as np
import pytest
from scipy import ndimage

import laynii_testing as lt


def test_data(name):
    return lt.load(lt.TEST_DATA / name)[1]


test_data.__test__ = False  # not a test, despite the name


# =============================================================================
# Data type conversion
# =============================================================================
def test_float_me_preserves_values(laynii):
    laynii.run("LN_FLOAT_ME", "-input", "lo_BOLD_intemp.nii.gz")
    out = laynii.load("lo_BOLD_intemp_float.nii.gz")
    expected = test_data("lo_BOLD_intemp.nii.gz").astype(np.float32)
    np.testing.assert_array_equal(out, expected)


def test_int_me_truncates(laynii):
    laynii.run("LN_INT_ME", "-input", "lo_BOLD_act.nii.gz")
    out = laynii.load("lo_BOLD_act_int16.nii.gz")
    np.testing.assert_array_equal(out, np.trunc(test_data("lo_BOLD_act.nii.gz")))


def test_short_me_keeps_three_decimals(laynii):
    laynii.run("LN_SHORT_ME", "-input", "lo_VASO_act.nii.gz",
               "-output", "short.nii.gz")
    out = laynii.load("short.nii.gz")
    np.testing.assert_allclose(out, test_data("lo_VASO_act.nii.gz"),
                               rtol=0, atol=1.001e-3)


# =============================================================================
# Time series statistics
# =============================================================================
def test_zscore_mean_std(laynii):
    laynii.run("LN2_ZSCORE", "-input", "lo_BOLD_intemp.nii.gz", "-mean",
               "-std")
    data = test_data("lo_BOLD_intemp.nii.gz").astype(np.float64)
    mean, std = data.mean(axis=-1), data.std(axis=-1)

    out_mean = laynii.load("lo_BOLD_intemp_mean.nii.gz")[..., 0]
    out_std = laynii.load("lo_BOLD_intemp_std.nii.gz")[..., 0]
    np.testing.assert_allclose(out_mean, mean, rtol=1e-5, atol=1e-3)
    np.testing.assert_allclose(out_std, std, rtol=1e-5, atol=1e-3)

    valid = std > 0
    z = laynii.load("lo_BOLD_intemp_zscore.nii.gz")
    expected = (data[valid] - mean[valid, None]) / std[valid, None]
    np.testing.assert_allclose(z[valid], expected, rtol=1e-4, atol=1e-4)
    assert np.all(np.isnan(z[~valid]))


@pytest.mark.parametrize("option", ["-box", "-gaus"])
def test_tempsmooth_preserves_mean_reduces_variance(laynii, option):
    laynii.run("LN_TEMPSMOOTH", "-input", "lo_BOLD_intemp.nii.gz", option,
               "1", "-output", "smooth.nii.gz")
    data = test_data("lo_BOLD_intemp.nii.gz")
    out = laynii.load("smooth.nii.gz")
    assert out.shape == data.shape
    np.testing.assert_allclose(out.mean(-1), data.mean(-1), rtol=0,
                               atol=1e-3 * np.abs(data).max())
    assert out.var(-1).sum() < 0.6 * data.var(-1).sum()


def test_noiseme_adds_requested_noise(laynii):
    laynii.run("LN_NOISEME", "-input", "lo_VASO_act.nii.gz", "-std", "1")
    diff = laynii.load("lo_VASO_act_noised.nii.gz") - test_data(
        "lo_VASO_act.nii.gz")
    assert abs(diff.mean()) < 0.05
    assert 0.95 < diff.std() < 1.05


# =============================================================================
# Simple voxel-wise and geometric operations
# =============================================================================
def test_reciprocal(laynii):
    laynii.run("LN2_RECIPROCAL", "-input", "lo_T1EPI.nii.gz")
    data = test_data("lo_T1EPI.nii.gz")
    out = laynii.load("lo_T1EPI_recip.nii.gz")
    nonzero = data != 0
    np.testing.assert_allclose(out[nonzero],
                               1e6 / np.clip(data[nonzero], 1.0, None),
                               rtol=1e-6)
    assert np.all(out[~nonzero] == 0)


def test_circshift_matches_numpy_roll(laynii):
    laynii.run("LN2_CIRCSHIFT", "-input", "sc_UNI.nii.gz", "-shift_neg",
               "10", "-axis", "2")
    out = laynii.load("sc_UNI_circshift-y.nii.gz")
    np.testing.assert_array_equal(out, np.roll(test_data("sc_UNI.nii.gz"),
                                               -10, axis=1))


@pytest.mark.xfail(reason="window is [i-range, i+range) and skips index 0, "
                          "instead of the documented range on either side",
                   strict=True)
def test_intpro_max_uses_symmetric_window(laynii):
    laynii.run("LN2_INTPRO", "-input", "sc_UNI.nii.gz", "-range", "3",
               "-type", "max")
    out = laynii.load("sc_UNI_maxip-z_range-3.nii.gz")
    data = test_data("sc_UNI.nii.gz").astype(np.float32)
    expected = ndimage.maximum_filter1d(data, size=7, axis=2, mode="nearest")
    np.testing.assert_array_equal(out, expected)


def test_zoom_crops_to_mask_bounding_box(laynii):
    laynii.run("LN_ZOOM", "-mask", "sc_layers_3dcolumns.nii.gz", "-input",
               "sc_UNI.nii.gz")
    mask = test_data("sc_layers_3dcolumns.nii.gz")
    lo = np.argwhere(mask > 0).min(0)
    hi = np.argwhere(mask > 0).max(0) + 1
    expected = test_data("sc_UNI.nii.gz")[lo[0]:hi[0], lo[1]:hi[1],
                                         lo[2]:hi[2]]
    np.testing.assert_array_equal(laynii.load("sc_UNI_zoomed.nii.gz"),
                                  expected)


def test_direct_smooth_keeps_constant_image(laynii):
    laynii.save("const.nii.gz", np.full((20, 20, 20), 3.0, np.float32))
    laynii.run("LN_DIRECT_SMOOTH", "-input", "const.nii.gz", "-FWHM", "2",
               "-direction", "3")
    np.testing.assert_allclose(laynii.load("const_smooth.nii.gz"), 3.0,
                               atol=1e-5)


def test_layer_smooth_keeps_constant_image(laynii):
    layers = test_data("sc_layers.nii.gz")
    laynii.save("const.nii.gz", np.where(layers > 0, 5.0, 0.0), dtype=np.float32,
                affine=lt.load(lt.TEST_DATA / "sc_layers.nii.gz")[0].affine)
    laynii.run("LN2_LAYER_SMOOTH", "-input", "const.nii.gz", "-layer_file",
               "sc_layers.nii.gz", "-FWHM", "1")
    out = laynii.load("const_layer_smoothed.nii.gz")
    np.testing.assert_allclose(out[layers > 0], 5.0, atol=1e-3)
    assert np.all(out[layers == 0] == 0)


def test_gradient_magnitude_matches_gradients(laynii):
    laynii.run("LN2_GRADIENTS", "-input", "Smagn_t-001.nii.gz")
    laynii.run("LN2_GRAMAG", "-input", "Smagn_t-001.nii.gz")
    gx, gy, gz = (laynii.load("Smagn_t-001_gradient_%s.nii.gz" % a)
                  for a in "xyz")
    np.testing.assert_allclose(laynii.load("Smagn_t-001_gramag.nii.gz"),
                               np.sqrt(gx**2 + gy**2 + gz**2), rtol=1e-5,
                               atol=1e-3)


def test_info_reports_dimensions(laynii):
    result = laynii.run("LN_INFO", "-input", "lo_T1EPI.nii.gz", "-NoPlot")
    assert re.search(r"162 X \| 162 Y \| 3 Z \| 1 T", result.output), \
        result.describe()


# =============================================================================
# Segmentation utilities
# =============================================================================
def test_copy_geometry(laynii):
    data = np.arange(4 * 5 * 6, dtype=np.float32).reshape(4, 5, 6)
    laynii.save("input.nii.gz", data)
    affine = np.diag([0.8, 0.8, 1.2, 1.0])
    affine[:3, 3] = [-10, 5, 3]
    laynii.save("reference.nii.gz", np.zeros((4, 5, 6, 2), np.int16),
                affine=affine)
    laynii.run("LN2_COPY_GEOMETRY", "-input", "input.nii.gz", "-geometry",
               "reference.nii.gz")
    img, out = lt.load(laynii.path("input_geomcopied.nii.gz"))
    np.testing.assert_array_equal(out, data)
    assert img.get_data_dtype() == np.float32
    np.testing.assert_allclose(img.affine, affine, atol=1e-6)
    np.testing.assert_allclose(img.header.get_zooms()[:3], (0.8, 0.8, 1.2),
                               atol=1e-6)


def test_copy_geometry_rejects_different_grid(laynii):
    laynii.save("input.nii.gz", np.zeros((4, 5, 6), np.float32))
    laynii.save("reference.nii.gz", np.zeros((4, 5, 7), np.float32))
    result = laynii.run("LN2_COPY_GEOMETRY", "-input", "input.nii.gz",
                        "-geometry", "reference.nii.gz", check=False)
    assert result.returncode > 0, result.describe()
    assert result.new_files == []


def test_rimify_relabels(laynii):
    # Swap the meaning of 1 and 2: input 1 becomes inner GM (2) and vice versa
    laynii.run("LN2_RIMIFY", "-input", "sc_rim.nii.gz", "-innergm", "1",
               "-outergm", "2", "-gm", "3", "-output", "out.nii.gz")
    rim = test_data("sc_rim.nii.gz")
    expected = np.choose(rim, [0, 2, 1, 3])
    np.testing.assert_array_equal(laynii.load("out.nii.gz"), expected)


def test_borderize_matches_6_connected_erosion(laynii):
    laynii.run("LN2_BORDERIZE", "-input", "sc_rim.nii.gz")
    rim = test_data("sc_rim.nii.gz")
    expected = np.zeros_like(rim)
    for label in (1, 2, 3):
        region = rim == label
        inner = ndimage.binary_erosion(region, border_value=1)
        expected[region & ~inner] = label
    np.testing.assert_array_equal(laynii.load("sc_rim_borders.nii.gz"),
                                  expected)


def test_connected_clusters_synthetic(laynii):
    blobs = np.zeros((30, 30, 30), np.int16)
    blobs[2:6, 2:6, 2:6] = 1        # 64 voxels
    blobs[10:15, 10:12, 10:20] = 1  # 100 voxels
    blobs[20:25, 20:25, 20:25] = 1  # 125 voxels
    blobs[27, 27, 27] = 1           # 1 voxel
    laynii.save("blobs.nii.gz", blobs)
    result = laynii.run("LN2_CONNECTED_CLUSTERS", "-input", "blobs.nii.gz")
    assert result.new_files == ["blobs_connected_clusters4.nii.gz"]
    out = laynii.load("blobs_connected_clusters4.nii.gz")
    assert np.array_equal(out > 0, blobs > 0)
    sizes = sorted(int(np.sum(out == k)) for k in range(1, 5))
    assert sizes == [1, 64, 100, 125]


def test_connected_clusters_uses_26_neighbourhood(laynii):
    diagonal = np.zeros((10, 10, 10), np.int16)
    diagonal[2, 2, 2] = diagonal[3, 3, 3] = 1
    laynii.save("diag.nii.gz", diagonal)
    result = laynii.run("LN2_CONNECTED_CLUSTERS", "-input", "diag.nii.gz")
    assert result.new_files == ["diag_connected_clusters1.nii.gz"]


def test_connected_clusters_real_data(laynii):
    result = laynii.run("LN2_CONNECTED_CLUSTERS", "-input", "sc_midGM.nii.gz")
    mid = test_data("sc_midGM.nii.gz") > 0
    _, n = ndimage.label(mid, ndimage.generate_binary_structure(3, 3))
    assert result.new_files == ["sc_midGM_connected_clusters%d.nii.gz" % n]


def test_ifpoints_and_voronoi(laynii):
    laynii.run("LN2_IFPOINTS", "-domain", "sc_midGM.nii.gz", "-nr_points",
               "10")
    domain = test_data("sc_midGM.nii.gz") > 0
    points = laynii.load("sc_midGM_points10.nii.gz")
    assert np.count_nonzero(points) == 10
    assert sorted(np.unique(points[points > 0])) == list(range(1, 11))
    assert np.all(domain[points > 0])

    laynii.run("LN2_VORONOI", "-domain", "sc_midGM.nii.gz", "-init",
               "sc_midGM_points10.nii.gz")
    cells = laynii.load("sc_midGM_voronoi.nii.gz")
    assert np.all(domain[cells > 0])
    np.testing.assert_array_equal(cells[points > 0], points[points > 0])
    # Every cell is one connected region containing its seed point
    for label in range(1, 11):
        _, n = ndimage.label(cells == label,
                             ndimage.generate_binary_structure(3, 3))
        assert n == 1, "cell %d is split into %d parts" % (label, n)


def test_geodistance_on_flat_domain(laynii):
    domain = np.zeros((41, 41, 3), np.int16)
    domain[:, :, 1] = 1
    init = np.zeros_like(domain)
    init[20, 20, 1] = 1
    laynii.save("domain.nii.gz", domain)
    laynii.save("init.nii.gz", init)
    laynii.run("LN2_GEODISTANCE", "-domain", "domain.nii.gz", "-init",
               "init.nii.gz", "-no_smooth", "-output", "geo")
    dist = laynii.load("geo.nii")
    assert np.all(dist[:, :, (0, 2)] == 0)
    x, y = np.meshgrid(np.arange(41), np.arange(41), indexing="ij")
    euclid = np.hypot(x - 20, y - 20)
    # The distance of the initial voxel is defined as half a voxel. Distances
    # are propagated through the voxel grid, so they are exact along axes and
    # diagonals, never shorter than euclidean, and at most ~8% longer.
    geo = dist[:, :, 1] - 0.5
    exact = (x == 20) | (y == 20) | (np.abs(x - 20) == np.abs(y - 20))
    np.testing.assert_allclose(geo[exact], euclid[exact], atol=1e-3)
    assert np.all(geo >= euclid - 1e-3)
    assert np.all(geo <= 1.09 * euclid + 1e-3)


@lt.known_memory_bug("thr_exceed[column - 1] is read before checking that "
                     "column > 0")
def test_mask_selects_columns_by_mean_score(laynii):
    laynii.run("LN2_MASK", "-scores", "lo_BOLD_act.nii.gz", "-columns",
               "lo_columns.nii.gz", "-mean_thr", "1", "-output", "mask.nii.gz",
               "-abs")
    scores = test_data("lo_BOLD_act.nii.gz")
    columns = test_data("lo_columns.nii.gz").astype(int)
    selected = [c for c in np.unique(columns[columns > 0])
                if np.abs(scores[columns == c]).mean() >= 1]
    expected = np.where(np.isin(columns, selected), columns, 0)
    np.testing.assert_array_equal(laynii.load("mask.nii.gz"), expected)


# =============================================================================
# Layers, profiles and columns
# =============================================================================
@pytest.fixture
def slab(laynii):
    """Flat synthetic cortex: WM border at z=5, GM z=6..25, CSF at z=26."""
    rim = np.zeros((30, 30, 40), np.int16)
    rim[:, :, 5] = 2
    rim[:, :, 6:26] = 3
    rim[:, :, 26] = 1
    laynii.save("slab.nii.gz", rim)
    return rim


def test_layers_on_flat_slab(laynii, slab):
    laynii.run("LN2_LAYERS", "-rim", "slab.nii.gz", "-nr_layers", "5",
               "-equivol")
    gm_z = np.arange(6, 26)

    metric = laynii.load("slab_metric_equidist.nii.gz")
    plane_mean = metric[:, :, gm_z].mean(axis=(0, 1))
    assert metric[:, :, gm_z].std(axis=(0, 1)).max() < 1e-3
    assert np.all(np.diff(plane_mean) > 0), "depth is not monotonic"
    assert np.corrcoef(plane_mean, gm_z)[0, 1] > 0.999
    assert 0 < plane_mean[0] < 0.1 and 0.9 < plane_mean[-1] < 1

    # 20 GM planes split into 5 layers of 4 planes each, layer 1 at the WM
    expected = np.repeat(np.arange(1, 6), 4)
    for name in ("slab_layers_equidist.nii.gz", "slab_layers_equivol.nii.gz"):
        layers = laynii.load(name)
        assert np.all(layers[slab != 3] == 0), name
        for z, layer in zip(gm_z, expected):
            assert np.all(layers[:, :, z] == layer), (name, z)

    mid = laynii.load("slab_midGM_equidist.nii.gz")
    assert sorted(np.unique(np.argwhere(mid > 0)[:, 2])) == [15, 16]


def test_layers_on_example_data(laynii):
    laynii.run("LN2_LAYERS", "-rim", "sc_rim.nii.gz", "-nr_layers", "10")
    rim = test_data("sc_rim.nii.gz")
    layers = laynii.load("sc_rim_layers_equidist.nii.gz")
    metric = laynii.load("sc_rim_metric_equidist.nii.gz")
    gm = rim == 3
    assert np.all(layers[~gm] == 0)
    assert np.all(layers[gm] > 0)
    assert sorted(np.unique(layers[gm])) == list(range(1, 11))
    assert metric.min() >= 0 and metric.max() <= 1

    near_wm = ndimage.binary_dilation(rim == 2) & gm
    near_csf = ndimage.binary_dilation(rim == 1) & gm
    assert layers[near_wm].mean() < 3 < 8 < layers[near_csf].mean()


def test_profile_matches_numpy(laynii):
    laynii.run("LN2_PROFILE", "-input", "sc_VASO_act.nii.gz", "-layers",
               "sc_layers.nii.gz", "-plot")
    values = test_data("sc_VASO_act.nii.gz")
    layers = test_data("sc_layers.nii.gz")
    rows = [line.split() for line in
            laynii.path("sc_VASO_act_profile.txt").read_text().splitlines()
            if line.strip()]
    assert len(rows) == layers.max()
    for row in rows:
        layer, mean, std, count = int(row[0]), *map(float, row[1:4])
        voxels = values[layers == layer]
        assert count == voxels.size
        np.testing.assert_allclose(mean, voxels.mean(), rtol=1e-4, atol=1e-5)
        np.testing.assert_allclose(std, voxels.std(), rtol=1e-4, atol=1e-5)


def test_columns(laynii):
    laynii.run("LN2_COLUMNS", "-rim", "sc_rim.nii.gz", "-midgm",
               "sc_midGM.nii.gz", "-nr_columns", "300")
    rim = test_data("sc_rim.nii.gz")
    columns = laynii.load("sc_rim_columns300.nii.gz")
    centroids = laynii.load("sc_rim_centroids300.nii.gz")
    assert np.all(rim[columns > 0] == 3)
    assert sorted(np.unique(columns[columns > 0])) == list(range(1, 301))
    assert np.count_nonzero(centroids) == 300
    assert sorted(np.unique(centroids[centroids > 0])) == list(range(1, 301))
