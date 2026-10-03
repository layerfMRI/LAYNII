"""Run every LayNii program on the example data and validate its outputs.

Each case lists the command line and the files it is expected to write. Every
nifti output is checked to be readable (catches truncated files), to have the
expected shape, to contain no NaN/Inf (unless allowed) and not to be empty.
Programs with dedicated numerical checks are additionally covered in
test_correctness.py.
"""
import fnmatch

import pytest

import laynii_testing as lt


class Case:
    """One program invocation and the outputs it must produce.

    xfail:       known bug that always makes the case fail.
    memory_bug:  known out-of-bounds access, see laynii_testing.memory_bug.
    """

    def __init__(self, id, program, args, outputs, like=None, xfail=None,
                 memory_bug=None, flaky=False, slow=False):
        self.id = id
        self.program = program
        self.args = args.split()
        self.outputs = outputs  # {filename or glob: check_nifti kwargs}
        self.like = like        # input whose spatial shape outputs share
        self.xfail = xfail
        self.memory_bug = memory_bug
        self.flaky = flaky
        self.slow = slow

    def param(self):
        marks = []
        if self.program in lt.WIP_PROGRAMS:
            marks.append(pytest.mark.wip)
        if self.slow:
            marks.append(pytest.mark.slow)
        if self.xfail:
            marks.append(pytest.mark.xfail(reason=self.xfail, strict=True))
        elif self.memory_bug:
            marks += lt.memory_bug(self.memory_bug, self.flaky)
        return pytest.param(self, id=self.id, marks=marks)


NONFINITE = {"allow_nonfinite": True}
ZEROS_OK = {"allow_all_zero": True}
TEXT = {"text": True}

SC = "sc_rim.nii.gz"
LO = "lo_BOLD_act.nii.gz"
DING = "Ding2016_occip_rim.nii.gz"
PHASE = "Sphase_t-001.nii.gz"
MAGN = "Smagn_t-001.nii.gz"
T2S = "Ding2016_occip_T2starweighted_filtered_for_tests.nii.gz"
CP0 = "Ding2016_occipital_rim_midGM_equidist_control_point0.nii.gz"

CASES = [
    # ------------------------------------------------------------------------
    # Layers and columns
    # ------------------------------------------------------------------------
    Case("LN2_LAYERS-equivol", "LN2_LAYERS",
         "-rim sc_rim.nii.gz -nr_layers 10 -equivol", {
             "sc_rim_layers_equidist.nii.gz": {"value_range": (0, 10)},
             "sc_rim_layers_equivol.nii.gz": {"value_range": (0, 10)},
             "sc_rim_metric_equidist.nii.gz": {"value_range": (0, 1)},
             "sc_rim_metric_equivol.nii.gz": {"value_range": (0, 1)},
             "sc_rim_midGM_equidist.nii.gz": {"value_range": (0, 1)},
             "sc_rim_midGM_equivol.nii.gz": {"value_range": (0, 1)},
         }, like=SC),
    Case("LN2_LAYERS-curvature-thickness-streamlines", "LN2_LAYERS",
         "-rim rim_M.nii.gz -nr_layers 3 -curvature -thickness -streamlines "
         "-output rimM", {
             "rimM_layers_equidist.nii": {"value_range": (0, 3)},
             "rimM_metric_equidist.nii": {"value_range": (0, 1)},
             "rimM_midGM_equidist.nii": {"value_range": (0, 1)},
             "rimM_curvature.nii": {"value_range": (-1, 1)},
             "rimM_curvature_binned.nii": {},
             "rimM_thickness.nii": {"value_range": (0, 100)},
             "rimM_streamline_vectors.nii": {"shape": (100, 100, 100, 3),
                                             "value_range": (-1, 1)},
         }, like="rim_M.nii.gz"),
    Case("LN2_LAYERS-incl_borders-equal_counts", "LN2_LAYERS",
         "-rim sc_rim.nii.gz -nr_layers 3 -incl_borders -equal_counts "
         "-output out", {
             "out_layers_equidist.nii": {"value_range": (0, 3)},
             "out_layers_equicount.nii": {"value_range": (0, 3)},
             "out_metric_equidist.nii": {"value_range": (0, 1)},
             "out_midGM_equidist.nii": {},
         }, like=SC),
    Case("LN2_COLUMNS", "LN2_COLUMNS",
         "-rim sc_rim.nii.gz -midgm sc_midGM.nii.gz -nr_columns 300", {
             "sc_rim_columns300.nii.gz": {"value_range": (0, 300)},
             "sc_rim_centroids300.nii.gz": {"value_range": (0, 300)},
         }, like=SC),
    Case("LN2_MULTILATERATE", "LN2_MULTILATERATE",
         "-rim %s -control_points %s -radius 10 -output out" % (DING, CP0), {
             "out_UV_coordinates.nii": {"shape": (130, 114, 107, 2)},
             "out_UV_axes.nii": {},
             "out_perimeter.nii": {},
             "out_perimeter_chunk.nii": {},
         }, like=DING),
    Case("LN2_HEXBIN", "LN2_HEXBIN",
         "-coord_uv ding_UV_coordinates.nii -radius 2", {
             "ding_UV_coordinates_hexbins2.nii": {},
         }, like=DING),
    Case("LN2_PATCH_FLATTEN", "LN2_PATCH_FLATTEN",
         "-values %s -coord_uv ding_UV_coordinates.nii "
         "-coord_d ding_metric_equidist.nii -domain ding_perimeter_chunk.nii "
         "-bins_u 50 -bins_v 50 -bins_d 21 -voronoi -density -output out"
         % T2S, {
             "out_flat_50x50x21_voronoi.nii": {"shape": (50, 50, 21)},
             "out_flat_50x50x21_density_voronoi.nii": {"shape": (50, 50, 21)},
             "out_flat_50x50x21_foldedcoords_voronoi.nii": {
                 "shape": (50, 50, 21, 3)},
         }),
    Case("LN2_PATCH_UNFLATTEN", "LN2_PATCH_UNFLATTEN",
         "-values ding_flat_30x30x5_voronoi.nii "
         "-coord_xyz ding_flat_30x30x5_foldedcoords_voronoi.nii "
         "-ref Ding2016_occip_ROI.nii.gz -output out", {
             "out_unflattened.nii": {},
             "out_unflattened_density.nii": {},
         }, like=DING),
    Case("LN2_PATCH_FLATTEN_2D", "LN2_PATCH_FLATTEN_2D",
         "-values slice_sc_VASO_act.nii.gz -coord_tan slice_geodist.nii "
         "-coord_rad slice_metric_equidist.nii -domain slice_sc_rim.nii.gz "
         "-bins_rad 21 -bins_tan 50 -output out", {
             "out_flat_21x50.nii": {"shape": (50, 21)},
         }),
    Case("LN2_UVD_FILTER", "LN2_UVD_FILTER",
         "-values %s -coord_uv ding_UV_coordinates.nii "
         "-coord_d ding_metric_equidist.nii -domain ding_perimeter_chunk.nii "
         "-radius 3 -height 0.25 -output out" % T2S, {
             "out_UVD_median_filter.nii": {},
         }, like=DING),
    Case("LN2_GEODISTANCE", "LN2_GEODISTANCE",
         "-domain Ding2016_occip_ROI.nii.gz -init %s -init_val 2 -output out"
         % CP0, {"out.nii": {"value_range": (0, 1000)}}, like=DING),
    Case("LN2_IFPOINTS", "LN2_IFPOINTS",
         "-domain sc_midGM.nii.gz -nr_points 10", {
             "sc_midGM_points10.nii.gz": {"value_range": (0, 10)},
             "sc_midGM_cells10.nii.gz": {"value_range": (0, 10)},
         }, like=SC),
    Case("LN2_VORONOI", "LN2_VORONOI",
         "-domain sc_midGM.nii.gz -init sc_midGM_points10.nii.gz", {
             "sc_midGM_voronoi.nii.gz": {"value_range": (0, 10)},
         }, like=SC),
    Case("LN2_ZERO_CROSSING", "LN2_ZERO_CROSSING",
         "-values lo_BOLD_act.nii.gz -domain lo_layers.nii.gz -output out", {
             "out_zero_crossing.nii": {"value_range": (0, 1)},
         }, like=LO),
    Case("LN2_CHOLMO", "LN2_CHOLMO",
         "-layers sc_layers.nii.gz -outer -nr_layers 3 -layer_thickness 0.4 "
         "-output padded.nii.gz", {
             "padded.nii.gz": {"value_range": (0, 23)},
         }, like=SC),
    Case("LN2_PROFILE", "LN2_PROFILE",
         "-input sc_VASO_act.nii.gz -layers sc_layers.nii.gz -plot", {
             "sc_VASO_act_profile.txt": TEXT,
         }),
    Case("LN2_LAYERDIMENSION", "LN2_LAYERDIMENSION",
         "-values lo_BOLD_act.nii.gz -layers lo_layers.nii.gz "
         "-columns lo_columns.nii.gz", {
             "lo_BOLD_act_layerdim.nii.gz": {"shape": (162, 162, 3, 10)},
         }),
    Case("LN2_MASK", "LN2_MASK",
         "-scores lo_BOLD_act.nii.gz -columns lo_columns.nii.gz -mean_thr 1 "
         "-output mask.nii.gz -abs", {"mask.nii.gz": {}}, like=LO,
         memory_bug="thr_exceed[column - 1] is read before checking that "
                    "column > 0"),
    Case("LN2_DIRECTIONALITY_BIN", "LN2_DIRECTIONALITY_BIN",
         "-input sc_midGM.nii.gz -layers sc_layers.nii.gz "
         "-columns sc_columns_3dcolumns.nii.gz -output out.nii.gz", {
             "out.nii.gz": {"value_range": (0, 1)},
         }, like=SC, flaky=True,
         memory_bug="background voxels (layer 0, column 0) are used as parcel "
                    "index -1, writing before the start of vec_nrVox_pacels"),
    Case("LN2_LAYER_SMOOTH", "LN2_LAYER_SMOOTH",
         "-input sc_VASO_act.nii.gz -layer_file sc_layers.nii.gz -FWHM 1", {
             "sc_VASO_act_layer_smoothed.nii.gz": {},
         }, like=SC),
    Case("LN_LAYER_SMOOTH-NoKissing", "LN_LAYER_SMOOTH",
         "-input sc_VASO_act.nii.gz -layer_file sc_layers.nii.gz -FWHM 0.3 "
         "-NoKissing", {
             "smoothed_sc_VASO_act.nii.gz": {},
             "hairy_brain.nii": {"value_range": (0, 1)},
         }, like=SC),
    Case("LN_3DCOLUMNS", "LN_3DCOLUMNS",
         "-layers sc_layers_3dcolumns.nii.gz "
         "-landmarks sc_landmarks_3dcolumns.nii.gz", {
             "sc_layers_3dcolumns_column_coordinates.nii.gz": {},
             "sc_layers_3dcolumns_finding_leaks.nii.gz": {},
             "sc_layers_3dcolumns_grow_from_left.nii.gz": {},
             "sc_layers_3dcolumns_grow_from_right.nii.gz": {},
             "sc_layers_3dcolumns_lateral_cord.nii.gz": {},
         }, like=SC),
    Case("LN_COLUMNAR_DIST", "LN_COLUMNAR_DIST",
         "-layers crop_layers.nii.gz -landmarks crop_landmarks.nii.gz", {
             "crop_layers_coordinates_final.nii.gz": {},
         }, like="crop_layers.nii.gz", slow=True),
    Case("LN_COLUMNAR_DIST-non_square", "LN_COLUMNAR_DIST",
         "-layers crop_layers_rect.nii.gz "
         "-landmarks crop_landmarks_rect.nii.gz", {
             "crop_layers_rect_coordinates_final.nii.gz": {},
         }, like="crop_layers_rect.nii.gz", slow=True,
         memory_bug="the x loop runs to size_y instead of size_x"),
    Case("LN_IMAGIRO", "LN_IMAGIRO",
         "-layers sc_layers_3dcolumns.nii.gz "
         "-columns sc_columns_3dcolumns.nii.gz -data sc_BOLD_act.nii.gz", {
             "sc_BOLD_act_unfolded.nii.gz": {"shape": (45, 15, 20)},
             "sc_BOLD_act_nr_voxels.nii.gz": {"shape": (45, 15, 20)},
         }, memory_bug="jx_stop is clamped with size_y_imagiro instead of "
                       "size_x_imagiro"),
    Case("LN_GROW_LAYERS", "LN_GROW_LAYERS", "-rim sc_rim.nii.gz", {
        "sc_rim_layers.nii.gz": {"value_range": (0, 20)},
    }, like=SC),
    Case("LN_LEAKY_LAYERS", "LN_LEAKY_LAYERS", "-rim lo_rim_LL.nii.gz", {
        "lo_rim_LL_leaky_layers.nii.gz": {"value_range": (0, 20)},
    }, like=LO),
    Case("LN_LOITUMA", "LN_LOITUMA",
         "-equidist sc_distlay_1000.nii.gz -leaky sc_leakylay_1000.nii.gz "
         "-FWHM 1 -nr_layers 10", {
             "equi_distance_layers.nii": {"value_range": (0, 10)},
             "equi_volume_layers.nii": {"value_range": (0, 10)},
             "leaky_layers.nii": {"value_range": (0, 10)},
         }, like=SC),
    Case("LN_CONLAY", "LN_CONLAY",
         "-layers lo_sc_layers.nii.gz -ref lo_T1EPI.nii.gz -subsample "
         "-output out.nii.gz", {"out.nii.gz": {"value_range": (0, 10)}},
         like=LO),
    Case("LN2_DEVEIN-ALF", "LN2_DEVEIN",
         "-layer_file lo_layers.nii.gz -column_file lo_columns.nii.gz "
         "-input lo_BOLD_act.nii.gz -ALF lo_ALF.nii.gz", {
             "lo_BOLD_act_deveinDeconv.nii.gz": {},
         }, like=LO),
    Case("LN2_DEVEIN-linear", "LN2_DEVEIN",
         "-layer_file lo_layers.nii.gz -input lo_BOLD_act.nii.gz -linear "
         "-output out.nii.gz", {"out_deveinLinear.nii.gz": {}}, like=LO),
    # ------------------------------------------------------------------------
    # Segmentation utilities
    # ------------------------------------------------------------------------
    Case("LN2_RIMIFY", "LN2_RIMIFY",
         "-input sc_rim.nii.gz -innergm 2 -outergm 1 -gm 3 -output out.nii.gz",
         {"out.nii.gz": {"value_range": (0, 3)}}, like=SC),
    Case("LN2_BORDERIZE", "LN2_BORDERIZE", "-input sc_rim.nii.gz", {
        "sc_rim_borders.nii.gz": {"value_range": (0, 3)},
    }, like=SC),
    Case("LN2_RIM_BORDERIZE", "LN2_RIM_BORDERIZE", "-rim sc_rim.nii.gz", {
        "sc_rim_borderized.nii.gz": {"value_range": (0, 3)},
    }, like=SC),
    Case("LN2_RIM_POLISH", "LN2_RIM_POLISH", "-rim sc_rim.nii.gz", {
        "sc_rim_polished.nii.gz": {"value_range": (0, 3)},
    }, like=SC),
    Case("LN2_CONNECTED_CLUSTERS", "LN2_CONNECTED_CLUSTERS",
         "-input sc_midGM.nii.gz", {
             "sc_midGM_connected_clusters*.nii.gz": {},
         }, like=SC),
    Case("LN2_NEIGHBORS", "LN2_NEIGHBORS", "-input lo_columns.nii.gz", {
        "lo_columns_neighbors.csv": TEXT,
    }),
    Case("LN_RAGRUG", "LN_RAGRUG", "-input sc_rim.nii.gz", {
        "sc_rim_ragrug.nii.gz": {"value_range": (1, 8)},
    }, like=SC),
    Case("LN2_COPY_GEOMETRY", "LN2_COPY_GEOMETRY",
         "-input lo_BOLD_intemp.nii.gz -geometry lo_Nulled_intemp.nii.gz", {
             "lo_BOLD_intemp_geomcopied.nii.gz": {"shape": (162, 162, 3, 80)},
         }),
    Case("LN_ZOOM", "LN_ZOOM",
         "-mask sc_layers_3dcolumns.nii.gz -input sc_UNI.nii.gz", {
             "sc_UNI_zoomed.nii.gz": {"shape": (79, 92, 15)},
         }),
    # ------------------------------------------------------------------------
    # Data type conversion and information
    # ------------------------------------------------------------------------
    Case("LN_FLOAT_ME", "LN_FLOAT_ME", "-input lo_BOLD_intemp.nii.gz", {
        "lo_BOLD_intemp_float.nii.gz": {"disk_dtype": "float32",
                                        "shape": (162, 162, 3, 80)},
    }),
    Case("LN_SHORT_ME", "LN_SHORT_ME",
         "-input lo_VASO_act.nii.gz -output short.nii.gz", {
             "short.nii.gz": {"disk_dtype": "int16"},
         }, like=LO),
    Case("LN_INT_ME", "LN_INT_ME", "-input lo_BOLD_act.nii.gz", {
        "lo_BOLD_act_int16.nii.gz": {"disk_dtype": "int16"},
    }, like=LO),
    # ------------------------------------------------------------------------
    # Time series
    # ------------------------------------------------------------------------
    Case("LN_BOCO-trialBOCO-shift", "LN_BOCO",
         "-Nulled lo_Nulled_intemp.nii.gz -BOLD lo_BOLD_intemp.nii.gz "
         "-trialBOCO 40 -shift", {
             "VASO_LN.nii": {"shape": (162, 162, 3, 80)},
             "BOLD_trialAV_LN.nii": {"shape": (162, 162, 3, 40)},
             "VASO_trialAV_LN.nii": dict(NONFINITE, shape=(162, 162, 3, 40)),
             "_shift_correlated.nii": {"shape": (162, 162, 3, 7)},
         }, memory_bug="NaN replacement loops over nr_voxels (all time points)"
                       " but the correlation image only has 7 volumes"),
    Case("LN_BOCO", "LN_BOCO",
         "-Nulled lo_Nulled_intemp.nii.gz -BOLD lo_BOLD_intemp.nii.gz "
         "-output out.nii.gz", {
             "out_VASO_LN.nii.gz": {"shape": (162, 162, 3, 80)},
         }),
    Case("LN_CORREL2FILES", "LN_CORREL2FILES",
         "-file1 lo_Nulled_intemp.nii.gz -file2 lo_BOLD_intemp.nii.gz", {
             "lo_Nulled_intemp_correlated.nii.gz": {"value_range": (-1, 1)},
         }, like=LO),
    Case("LN_EXTREMETR", "LN_EXTREMETR", "-input lo_BOLD_intemp.nii.gz", {
        "lo_BOLD_intemp_MaxTR.nii.gz": {},
        "lo_BOLD_intemp_MinTR.nii.gz": {},
    }, like=LO),
    Case("LN_SKEW", "LN_SKEW", "-input lo_BOLD_intemp.nii.gz", {
        "lo_BOLD_intemp_%s.nii.gz" % s: {} for s in (
            "autocorr", "imageSNR", "kurt", "local_gradient", "mean", "noise",
            "overall_correl", "skew", "stdev", "tSNR")
    }, like=LO),
    Case("LN_TEMPSMOOTH-box", "LN_TEMPSMOOTH",
         "-input lo_BOLD_intemp.nii.gz -box 1", {
             "lo_BOLD_intemp_tempsmooth.nii.gz": {"shape": (162, 162, 3, 80)},
         }),
    Case("LN_TEMPSMOOTH-gaus", "LN_TEMPSMOOTH",
         "-input lo_BOLD_intemp.nii.gz -gaus 1 -output out.nii.gz", {
             "out.nii.gz": {"shape": (162, 162, 3, 80)},
         }),
    Case("LN_TRIAL", "LN_TRIAL", "-input lo_BOLD_intemp.nii.gz -trialdur 20", {
        "lo_BOLD_intemp_TrialAverage.nii.gz": {"shape": (162, 162, 3, 20)},
    }, xfail="output nvox is computed as nr_voxels / nr_trials instead of "
             "nr_voxels * trial_dur, so the file is truncated (heap overflow)"),
    Case("LN_NOISE_KERNEL", "LN_NOISE_KERNEL",
         "-input lo_Nulled_intemp.nii.gz -kernel_size 7", {
             "lo_Nulled_intemp_fPSF.nii.gz": {"shape": (7, 7, 7)},
         }),
    Case("LN2_FRISGO-simple", "LN2_FRISGO",
         "-input lo_BOLD_intemp.nii.gz -simple -output out.nii.gz", {
             "out_frisgo-simple.nii.gz": {"shape": (162, 162, 3, 80)},
         }),
    Case("LN2_FRISGO-spline", "LN2_FRISGO",
         "-input lo_BOLD_intemp.nii.gz -spline -output out.nii.gz", {
             "out_frisgo-spline.nii.gz": {"shape": (162, 162, 3, 80)},
         }),
    Case("LN2_FRISGO-lpass", "LN2_FRISGO",
         "-input lo_BOLD_intemp.nii.gz -lpass 1.0 -output out.nii.gz", {
             "out_frisgo-lpass.nii.gz": {"shape": (162, 162, 3, 80)},
         }),
    Case("LN2_FRISGO-box", "LN2_FRISGO",
         "-input lo_BOLD_intemp.nii.gz -box 1 -output out.nii.gz", {
             "out_frisgo-box.nii.gz": {"shape": (162, 162, 3, 80)},
         }, xfail="-box is documented in the help but rejected as invalid"),
    Case("LN2_DESPIKE", "LN2_DESPIKE", "-input lo_BOLD_intemp.nii.gz", {
        "lo_BOLD_intemp_despike-simple.nii.gz": {"shape": (162, 162, 3, 80)},
        "lo_BOLD_intemp_despike-outlier_counts.nii.gz": {},
    }),
    Case("LN2_ZSCORE", "LN2_ZSCORE", "-input lo_BOLD_intemp.nii.gz -mean -std",
         {
             # Voxels with zero variance have an undefined z-score
             "lo_BOLD_intemp_zscore.nii.gz": dict(NONFINITE,
                                                  shape=(162, 162, 3, 80)),
             "lo_BOLD_intemp_mean.nii.gz": {},
             "lo_BOLD_intemp_std.nii.gz": {},
         }, like=LO),
    Case("LN2_SENSITIVITY", "LN2_SENSITIVITY", "-input lo_BOLD_intemp.nii.gz",
         {"lo_BOLD_intemp_sensitivity.nii.gz": {}}, like=LO),
    Case("LN2_SPECIFICITY", "LN2_SPECIFICITY", "-input lo_BOLD_intemp.nii.gz",
         {"lo_BOLD_intemp_specificity.nii.gz": {"value_range": (0, 1)}},
         like=LO),
    Case("LN_PHYSIO_PARS", "LN_PHYSIO_PARS", "synthetic.puls out.txt", {
        "out.txt": TEXT,
    }),
    # ------------------------------------------------------------------------
    # Image processing
    # ------------------------------------------------------------------------
    Case("LN_MP2RAGE_DNOISE", "LN_MP2RAGE_DNOISE",
         "-INV1 sc_INV1.nii.gz -INV2 sc_INV2.nii.gz -UNI sc_UNI.nii.gz", {
             # Background voxels with INV1 = INV2 = 0 divide by zero
             "sc_UNI_denoised.nii.gz": NONFINITE,
             "sc_UNI_border_enhance.nii.gz": NONFINITE,
         }, like=SC),
    Case("LN_DIRECT_SMOOTH", "LN_DIRECT_SMOOTH",
         "-input sc_UNI.nii.gz -FWHM 2 -direction 3", {
             "sc_UNI_smooth.nii.gz": {},
         }, like=SC),
    Case("LN_GRADSMOOTH", "LN_GRADSMOOTH",
         "-gradfile lo_gradT1.nii.gz -input lo_VASO_act.nii.gz -FWHM 1 "
         "-within -selectivity 0.1", {"lo_VASO_act_smoothed.nii.gz": {}},
         like=LO),
    Case("LN_GFACTOR", "LN_GFACTOR",
         "-input sc_INV2.nii.gz -variance 1 -direction 1 -grappa 2 "
         "-cutoff 200", {
             "sc_INV2_Amplified_GRAPPA.nii.gz": {},
             "sc_INV2_Gfactormap.nii.gz": {},
             "sc_INV2_Gfactormap_binary.nii.gz": {"value_range": (0, 1)},
         }, like=SC, slow=True),
    Case("LN_NOISEME", "LN_NOISEME", "-input lo_VASO_act.nii.gz -std 1", {
        "lo_VASO_act_noised.nii.gz": {},
    }, like=LO),
    Case("LN_INTPRO", "LN_INTPRO",
         "-image sc_UNI.nii.gz -min -direction 2 -range 3", {
             "sc_UNI_collapsed.nii.gz": {},
         }, like=SC),
    Case("LN2_INTPRO", "LN2_INTPRO",
         "-input sc_UNI.nii.gz -range 3 -type max", {
             "sc_UNI_maxip-z_range-3.nii.gz": {},
         }, like=SC),
    Case("LN2_RECIPROCAL", "LN2_RECIPROCAL", "-input lo_T1EPI.nii.gz", {
        "lo_T1EPI_recip.nii.gz": {"value_range": (0, 1e6)},
    }, like=LO),
    Case("LN2_SNAPCAST-mean", "LN2_SNAPCAST", "-input %s" % T2S, {
        "Ding2016_occip_T2starweighted_filtered_for_tests_snapcast-mean_"
        "steps-5.nii.gz": {"shape": (130, 130, 6)},
    }),
    Case("LN2_SNAPCAST-max", "LN2_SNAPCAST",
         "-input %s -type max -steps 3 -output out.nii.gz" % T2S, {
             "out_snapcast-max_steps-3.nii.gz": {"shape": (130, 130, 6)},
         }),
    Case("LN_INFO", "LN_INFO", "-input lo_T1EPI.nii.gz -NoPlot", {}),
    # ------------------------------------------------------------------------
    # Spatial derivatives
    # ------------------------------------------------------------------------
    Case("LN2_GRADIENTS", "LN2_GRADIENTS", "-input %s" % MAGN, {
        "Smagn_t-001_gradient_%s.nii.gz" % ax: {} for ax in "xyz"
    }, like=MAGN),
    Case("LN2_GRAMAG", "LN2_GRAMAG", "-input %s" % MAGN, {
        "Smagn_t-001_gramag.nii.gz": {},
    }, like=MAGN),
    Case("LN2_LAPLACIAN", "LN2_LAPLACIAN", "-input %s" % MAGN, {
        "Smagn_t-001_laplacian.nii.gz": {},
    }, like=MAGN),
    Case("LN2_PHASE_GRADIENTS", "LN2_PHASE_GRADIENTS", "-input %s -int13"
         % PHASE, {
             "Sphase_t-001_phase_gradient_%s.nii.gz" % ax: {
                 "value_range": (-3.1416, 3.1416)} for ax in "xyz"
         }, like=PHASE),
    Case("LN2_PHASE_JOLT", "LN2_PHASE_JOLT", "-input %s -int13" % PHASE, {
        "Sphase_t-001_phase_jolt.nii.gz": {},
    }, like=PHASE),
    Case("LN2_PHASE_LAPLACIAN", "LN2_PHASE_LAPLACIAN", "-input %s -int13"
         % PHASE, {"Sphase_t-001_phase_laplacian.nii.gz": {}}, like=PHASE),
    # ------------------------------------------------------------------------
    # Work in progress programs
    # ------------------------------------------------------------------------
    Case("LN2_CIRCSHIFT", "LN2_CIRCSHIFT",
         "-input sc_UNI.nii.gz -shift_neg 10 -axis 2", {
             "sc_UNI_circshift-y.nii.gz": {},
         }, like=SC),
    Case("LN2_PEAK_DETECT", "LN2_PEAK_DETECT",
         "-values lo_BOLD_act.nii.gz -max", {
             "lo_BOLD_act_peaks.nii.gz": {"value_range": (0, 1)},
         }, like=LO),
    Case("LN2_SKELETONIZE", "LN2_SKELETONIZE", "-input sc_midGM.nii.gz", {
        "sc_midGM_*.nii.gz": {},
    }, like=SC),
    Case("LN2_WINDOWED_COUNTER_2D", "LN2_WINDOWED_COUNTER_2D",
         "-input sc_midGM_points10.nii.gz -radius 5", {
             "sc_midGM_points10_counts_rad-5.nii.gz": {},
         }, like=SC),
    Case("LN2_REGRESS_OUT", "LN2_REGRESS_OUT",
         "-input1 lo_BOLD_intemp.nii.gz -input2 lo_Nulled_intemp.nii.gz", {
             "lo_BOLD_intemp_fitted.nii.gz": dict(NONFINITE,
                                                  shape=(162, 162, 3, 80)),
             "lo_BOLD_intemp_residuals.nii.gz": dict(NONFINITE,
                                                     shape=(162, 162, 3, 80)),
             "lo_BOLD_intemp_slope.nii.gz": NONFINITE,
             "lo_BOLD_intemp_intercept.nii.gz": {},
         }, like=LO),
    Case("LN3_LAYERS", "LN3_LAYERS", "-rim rim_M.nii.gz", {
        "rim_M_*layers.nii.gz": {"value_range": (0, 3)},
    }, like="rim_M.nii.gz",
        memory_bug="neighbour lookup through voi_id_inv reads past voi_rim"),
    Case("LN3_NOLAD", "LN3_NOLAD", "-input %s" % MAGN, {
        "Smagn_t-001_*.nii.gz": NONFINITE,
    }, like=MAGN),
    Case("LN2_UVD_LSTSQR", "LN2_UVD_LSTSQR",
         "-values %s -coord_uv ding_UV_coordinates.nii "
         "-coord_d ding_metric_equidist.nii -radius 3 -height 0.25 "
         "-output out" % T2S, {"out_*.nii*": {}}, like=DING,
         xfail="aborts (SIGABRT) on the example data"),
]


def test_every_program_has_a_case():
    covered = {case.program for case in CASES}
    missing = sorted(set(lt.ALL_PROGRAMS) - covered)
    assert not missing, "add test cases for: %s" % ", ".join(missing)


@pytest.mark.parametrize("case", [c.param() for c in CASES])
def test_program(case, laynii, derived):
    result = laynii.run(case.program, *case.args)

    shape = None
    if case.like is not None:
        shape = lt.spatial_shape(case.like, [lt.TEST_DATA, derived])

    for pattern, checks in case.outputs.items():
        checks = dict(checks)
        matches = fnmatch.filter(result.new_files, pattern)
        assert matches, "%s did not write %s; new files: %s\n%s" % (
            case.program, pattern, result.new_files, result.describe())
        for name in matches:
            if checks.pop("text", False):
                text = laynii.path(name).read_text()
                assert text.strip(), "%s is empty" % name
                continue
            if "shape" not in checks and shape is not None:
                checks["shape"] = shape
            lt.check_nifti(laynii.path(name), **checks)
