
// TODO(Faruk): Curvature shows some artifacts in rim_circles test case. Needs
// further investiation.
// TODO(Faruk): First implementation of equi-volume works but it needs
// extensive testing. I might need to smooth previous derivatives
// (hotspots, curvature etc.) a bit to make layering smoother on empirical
// data overall.
// TODO(Faruk): Memory usage is a bit sloppy for now. Low priority but might
// need to have a look at it in the future if we start hitting ram limits.
// NOTE(Faruk): Might put neighbour visits into a function.
// NOTE(Faruk): Might replace the smoothing code with a faster method later.
// NOTE(Faruk): Might use bspline weights to smooth curvature maps a bit.
// NOTE(Faruk): Might be better to use step 1 id's to define columns.


#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_LAYERS: Generates equi-distant cortical gray matter layers with\n"
    "            an option to also generate equi-volume layers.\n"
    "\n"
    "Usage:\n"
    "    LN2_LAYERS -rim rim.nii\n"
    "    LN2_LAYERS -rim rim.nii -nr_layers 3\n"
    "    LN2_LAYERS -rim rim.nii -nr_layers 3 -equivol\n"
    "    LN2_LAYERS -rim rim.nii -nr_layers 3 -equivol -iter_smooth 1000\n"
    "    ../LN2_LAYERS -rim sc_rim.nii -nr_layers 10 -equivol \n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -rim          : Specify input dataset. Use 1 to code outer gray\n"
    "                    matter surface (facing mostly CSF), 2 to code inner\n"
    "                    gray matter surdafe (facing mostly white matter),\n"
    "                    and 3 to code pure gray matter voxels.\n"
    "                    note that values 1 and 2 will not be included in the\n"
    "                    layerification, this is in contrast to the programs\n"
    "                    LN_GROW_LAYERS and LN_LEAKY LAYERS \n"
    "    -nr_layers    : Number of layers. Default is 3.\n"
    "    -equivol      : (Optional) Create equi-volume layers. We do not\n"
    "                    recommend this option if your rim file is above 0.3mm\n"
    "                    resolution. You can always upsample your rim file to\n"
    "                    a higher resolution first (<0.3mm) and then use this\n"
    "                    option.\n"
    "    -iter_smooth  : (Optional) Number of smoothing iterations. Default\n"
    "                    is 100. Only used together with '-equivol' flag. Use\n"
    "                    larger values when equi-volume layers are jagged.\n"
    "    -debug        : (Optional) Save extra intermediate outputs.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n"
    "Notes:\n"
    "    - You can find further explanation of this algorithm at:\n"
    "      <https://thingsonthings.org/ln2_layers>\n"
    "    - For a general discussion on equi-volume layering:\n"
    "      <https://layerfmri.com/equivol>\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

bool use_outpath = false;
    nifti_image *nii1 = NULL;
    char *fin = NULL, *fout = NULL;
    uint16_t ac, nr_layers = 3;
    uint16_t iter_smooth = 100;
    bool mode_equivol = false, mode_debug = false;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -rim\n");
                return 1;
            }
            fin = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-nr_layers")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -nr_layers\n");
            } else {
                nr_layers = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-iter_smooth")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -iter_smooth\n");
            } else {
                iter_smooth = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-equivol")) {
            mode_equivol = true;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-rim'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN2_LAYERS");
    log_nifti_descriptives(nii1);

    cout << "\n  Nr. layers: " << nr_layers << endl;

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    const float dX = nii1->pixdim[1];
    const float dY = nii1->pixdim[2];
    const float dZ = nii1->pixdim[3];

    // Short diagonals
    const float dia_xy = sqrt(dX * dX + dY * dY);
    const float dia_xz = sqrt(dX * dX + dZ * dZ);
    const float dia_yz = sqrt(dY * dY + dZ * dZ);
    // Long diagonals
    const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_rim = copy_nifti_as_int16(nii1);
    int16_t* nii_rim_data = static_cast<int16_t*>(nii_rim->data);

    // Prepare required nifti images
    nifti_image* nii_layers  = copy_nifti_as_int16(nii_rim);
    int16_t* nii_layers_data = static_cast<int16_t*>(nii_layers->data);
    // Setting zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_layers_data + i) = 0;
    }

    nifti_image* innerGM_step = copy_nifti_as_int16(nii_layers);
    int16_t* innerGM_step_data = static_cast<int16_t*>(innerGM_step->data);
    nifti_image* innerGM_dist = copy_nifti_as_float32(nii_layers);
    float* innerGM_dist_data = static_cast<float*>(innerGM_dist->data);

    nifti_image* outerGM_step = copy_nifti_as_int16(nii_layers);
    int16_t* outerGM_step_data = static_cast<int16_t*>(outerGM_step->data);
    nifti_image* outerGM_dist = copy_nifti_as_float32(nii_layers);
    float* outerGM_dist_data = static_cast<float*>(outerGM_dist->data);

    // nifti_image* err_dist = copy_nifti_as_float16(nii_layers);
    // short* err_dist_data = static_cast<short*>(err_dist->data);

    nifti_image* innerGM_id = copy_nifti_as_int32(nii_layers);
    int32_t* innerGM_id_data = static_cast<int32_t*>(innerGM_id->data);
    nifti_image* outerGM_id = copy_nifti_as_int32(nii_layers);
    int32_t* outerGM_id_data = static_cast<int32_t*>(outerGM_id->data);

    nifti_image* innerGM_prevstep_id = copy_nifti_as_int32(nii_layers);
    int32_t* innerGM_prevstep_id_data =
        static_cast<int32_t*>(innerGM_prevstep_id->data);
    nifti_image* outerGM_prevstep_id = copy_nifti_as_int32(nii_layers);
    int32_t* outerGM_prevstep_id_data =
        static_cast<int32_t*>(outerGM_prevstep_id->data);
    nifti_image* normdist = copy_nifti_as_float32(nii_layers);
    float* normdist_data = static_cast<float*>(normdist->data);
    nifti_image* normdistdiff = copy_nifti_as_float32(nii_layers);
    float* normdistdiff_data = static_cast<float*>(normdistdiff->data);

    nifti_image* nii_columns = copy_nifti_as_int32(nii_layers);
    int32_t* nii_columns_data = static_cast<int32_t*>(nii_columns->data);

    nifti_image* midGM = copy_nifti_as_int16(nii_layers);
    int16_t* midGM_data = static_cast<int16_t*>(midGM->data);
    nifti_image* midGM_id = copy_nifti_as_int32(nii_layers);
    int32_t* midGM_id_data = static_cast<int32_t*>(midGM_id->data);

    nifti_image* hotspots = copy_nifti_as_int32(nii_layers);
    int32_t* hotspots_data = static_cast<int32_t*>(hotspots->data);
    nifti_image* curvature = copy_nifti_as_float32(nii_layers);
    float* curvature_data = static_cast<float*>(curvature->data);

    // ========================================================================
    // Grow from WM
    // ========================================================================
    cout << "\n  Start growing from inner GM (WM-facing border)..." << endl;

    // Initialize grow volume
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 2) {  // WM boundary voxels within GM
            *(innerGM_step_data + i) = 1;
            *(innerGM_dist_data + i) = 0.;
            *(innerGM_id_data + i) = i;
        } else {
            *(innerGM_step_data + i) = 0.;
            *(innerGM_dist_data + i) = 0.;
        }
    }

    uint16_t grow_step = 1;
    uint32_t voxel_counter = nr_voxels;
    uint32_t ix, iy, iz, j, k;
    float d;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(innerGM_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dX;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dX;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dY;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dY;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dZ;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dZ;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xy;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_yz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(innerGM_dist_data + i) + dia_xyz;
                        if (d < *(innerGM_dist_data + j)
                            || *(innerGM_dist_data + j) == 0) {
                            *(innerGM_dist_data + j) = d;
                            *(innerGM_step_data + j) = grow_step + 1;
                            *(innerGM_id_data + j) = *(innerGM_id_data + i);
                            *(innerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }
    if (mode_debug) {
        save_output_nifti(fout, "innerGM_step", innerGM_step, false);
        save_output_nifti(fout, "innerGM_dist", innerGM_dist, false);
        save_output_nifti(fout, "innerGM_id", innerGM_id, false);
    }

    // ========================================================================
    // Grow from CSF
    // ========================================================================
    cout << "\n  Start growing from outer GM..." << endl;

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 1) {
            *(outerGM_step_data + i) = 1.;
            *(outerGM_dist_data + i) = 0.;
            *(outerGM_id_data + i) = i;
        } else {
            *(outerGM_step_data + i) = 0.;
            *(outerGM_dist_data + i) = 0.;
        }
    }

    grow_step = 1, voxel_counter = nr_voxels;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(outerGM_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dX;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dX;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dY;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dY;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dZ;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dZ;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xy;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_yz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(outerGM_dist_data + i) + dia_xyz;
                        if (d < *(outerGM_dist_data + j)
                            || *(outerGM_dist_data + j) == 0) {
                            *(outerGM_dist_data + j) = d;
                            *(outerGM_step_data + j) = grow_step + 1;
                            *(outerGM_id_data + j) = *(outerGM_id_data + i);
                            *(outerGM_prevstep_id_data + j) = i;
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }
    if (mode_debug) {
        save_output_nifti(fout, "outerGM_step", outerGM_step, false);
        save_output_nifti(fout, "outerGM_dist", outerGM_dist, false);
        save_output_nifti(fout, "outerGM_id", outerGM_id, false);
    }

    // ========================================================================
    // Layers
    // ========================================================================
    cout << "\n  Start layering (equi-distant)..." << endl;
    float x, y, z, wm_x, wm_y, wm_z, gm_x, gm_y, gm_z;

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            tie(wm_x, wm_y, wm_z) = ind2sub_3D(*(innerGM_id_data + i),
                                               size_x, size_y);
            tie(gm_x, gm_y, gm_z) = ind2sub_3D(*(outerGM_id_data + i),
                                               size_x, size_y);

            // // Normalize distance
            // float dist1 = dist(x, y, z, wm_x, wm_y, wm_z, dX, dY, dZ);
            // float dist2 = dist(x, y, z, gm_x, gm_y, gm_z, dX, dY, dZ);
            // float dist_normalized = dist1 / (dist1 + dist2);

            // Normalize distance (completely discrete)
            float dist1 = *(innerGM_dist_data + i);
            float dist2 = *(outerGM_dist_data + i);
            float total_dist = dist1 + dist2;;
            float dist_normalized = dist1 / total_dist;

            // To export equi-distant metric in a simple 0-1 range
            *(normdist_data + i) = dist_normalized;
            // Difference of normalized distances
            *(normdistdiff_data + i) = (dist1 - dist2) / total_dist;

            // Cast distances to integers as number of desired layers
            if (dist_normalized != 0) {
                *(nii_layers_data + i) = ceil(nr_layers * dist_normalized);
            } else {
                *(nii_layers_data + i) = 1;
            }

            // Count inner and outer GM anchor voxels
            j = *(innerGM_id_data + i);
            *(hotspots_data + j) += 1;
            j = *(outerGM_id_data + i);
            *(hotspots_data + j) -= 1;
        }
    }
    save_output_nifti(fout, "metric_equidist", normdist);
    save_output_nifti(fout, "layers_equidist", nii_layers, true);

    if (mode_debug) {
        save_output_nifti(fout, "hotspots", hotspots, false);
        save_output_nifti(fout, "normdistdiff_equidist", normdistdiff, false);
    }

    // ========================================================================
    // Middle gray matter
    // ========================================================================
    cout << "\n  Start finding middle gray matter (equi-distant)..." << endl;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            // Check sign changes in normalized distance differences between
            // neighbouring voxels on a column path (a.k.a. streamline)
            if (*(normdistdiff_data + i) == 0) {
                *(midGM_data + i) = 1;
                *(midGM_id_data + i) = i;
            } else {
                float m = *(normdistdiff_data + i);
                float n;

                // Inner neighbour
                j = *(innerGM_prevstep_id_data + i);
                if (*(nii_rim_data + j) == 3) {
                    n = *(normdistdiff_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if ( m*m < n*n) {
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = j;
                        } else {  // Equal +/- normalized distance
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = i;  // On purpose
                        }
                    }
                }

                // Outer neighbour
                j = *(outerGM_prevstep_id_data + i);
                if (*(nii_rim_data + j) == 3) {
                    n = *(normdistdiff_data + j);
                    if (signbit(m) - signbit(n) != 0) {
                        if (m*m < n*n) {
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                        } else if (m*m > n*n) {  // Closer to prev. step
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = j;
                        } else {  // Equal +/- normalized distance
                            *(midGM_data + i) = 1;
                            *(midGM_id_data + i) = i;
                            *(midGM_data + j) = 1;
                            *(midGM_id_data + j) = i;  // On purpose
                        }
                    }
                }
            }
        }
    }
    save_output_nifti(fout, "midGM_equidist", midGM, true);

    // ========================================================================
    // Columns
    // ========================================================================
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            // Approximate curvature measurement per column/streamline
            j = *(innerGM_id_data + i);
            k = *(outerGM_id_data + i);  // These values are negative
            *(curvature_data + i) = *(hotspots_data + j) + *(hotspots_data + k);
            *(curvature_data + i) /=
                max(*(hotspots_data + j), -*(hotspots_data + k));  // normalize

            // Re-assign mid-GM id based on curvature
            if (*(curvature_data + i) >= 0) {  // Gyrus
                *(nii_columns_data + i) = j;
            } else {
                *(nii_columns_data + i) = k;
            }

            // MiddleGM ids are used to find centroids in the next step
            if (*(midGM_data + i) == 1) {
                *(midGM_id_data + i) = *(nii_columns_data + i);
            }
        }
    }
    if (mode_debug) {
        save_output_nifti(fout, "curvature_init", curvature, false);
    }

    // ========================================================================
    // Find Middle Gray Matter centroids
    // ========================================================================
    // NOTE(Faruk): I am a bit sluggish with memory usage here. Might optimize
    // later by switching nifti images to vectors.
    nifti_image* coords_x = copy_nifti_as_float32(nii_rim);
    float* coords_x_data = static_cast<float*>(coords_x->data);
    nifti_image* coords_y = copy_nifti_as_float32(nii_rim);
    float* coords_y_data = static_cast<float*>(coords_y->data);
    nifti_image* coords_z = copy_nifti_as_float32(nii_rim);
    float* coords_z_data = static_cast<float*>(coords_z->data);
    nifti_image* coords_count = copy_nifti_as_int32(nii_rim);
    int32_t* coords_count_data = static_cast<int32_t*>(coords_count->data);
    nifti_image* centroid = copy_nifti_as_int32(nii_rim);
    int32_t* centroid_data = static_cast<int32_t*>(centroid->data);
    nifti_image* midGM_centroid_id = copy_nifti_as_int32(nii_columns);
    int32_t* midGM_centroid_id_data = static_cast<int32_t*>(midGM_centroid_id->data);

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(coords_x_data + i) = 0;
        *(coords_y_data + i) = 0;
        *(coords_z_data + i) = 0;
        *(coords_count_data + i) = 0;
        *(centroid_data + i) = 0;
    }

    // Sum x, y, z coordinates of same-column middle GM voxels
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(midGM_data + i) == 1) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            j = *(midGM_id_data + i);  // used to determine storage voxel
            *(coords_x_data + j) += x;
            *(coords_y_data + j) += y;
            *(coords_z_data + j) += z;
            *(coords_count_data + j) += 1;
        }
    }
    // Divide summed coordinates to find centroid
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(coords_count_data + i) != 0) {
            // Assign centroid id in place of inner/outer border voxel id
            x = floor(*(coords_x_data + i) / *(coords_count_data + i));
            y = floor(*(coords_y_data + i) / *(coords_count_data + i));
            z = floor(*(coords_z_data + i) / *(coords_count_data + i));
            j = sub2ind_3D(x, y, z, size_x, size_y);
            *(centroid_data + i) = j;
        }
    }
    // Map new centroid IDs to columns/streamlines
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            j = *(nii_columns_data + i);
            *(midGM_centroid_id_data + i) = *(centroid_data + j);

            if (*(midGM_data + i) == 1) {  // Update Mid GM id
                *(midGM_id_data + i) = *(centroid_data + j);
            }
        }
    }
    if (mode_debug) {
        save_output_nifti(fout, "midGM_equidist_id", midGM_id, false);
        save_output_nifti(fout, "columns", midGM_centroid_id, false);
    }

    // ========================================================================
    // Update curvature along column/streamline based on midGM curvature
    // ========================================================================
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            j = *(midGM_centroid_id_data + i);
            *(curvature_data + i) = *(curvature_data + j);
        }
    }
    if (mode_debug) {
        save_output_nifti(fout, "curvature", curvature, true);
    }

    // ========================================================================
    // Equi-volume layers
    // ========================================================================
    if (mode_equivol) {
        cout << "\n  Start equi-volume stage..." << endl;

        nifti_image* hotspots_i = copy_nifti_as_float32(nii_rim);
        float* hotspots_i_data = static_cast<float*>(hotspots_i->data);
        nifti_image* hotspots_o = copy_nifti_as_float32(nii_rim);
        float* hotspots_o_data = static_cast<float*>(hotspots_o->data);

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(hotspots_i_data + i) = 0;
            *(hotspots_o_data + i) = 0;
        }

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_rim_data + i) == 3) {
                // Find inner/outer anchors
                j = *(innerGM_id_data + i);
                k = *(outerGM_id_data + i);

                // Count how many voxels fall inner and outer shells from MidGM
                if (*(curvature_data + i) < 0) {
                    if (*(normdistdiff_data + i) <= 0) {
                        *(hotspots_i_data + k) += 1;
                    }
                    if (*(normdistdiff_data + i) >= 0) {
                        *(hotspots_o_data + k) += 1;
                    }
                }
                if (*(curvature_data + i) > 0) {
                    if (*(normdistdiff_data + i) <= 0) {
                        *(hotspots_i_data + j) += 1;
                    }
                    if (*(normdistdiff_data + i) >= 0) {
                        *(hotspots_o_data + j) += 1;
                    }
                }
            }
        }
        if (mode_debug) {
            save_output_nifti(fout, "hotspots_in", hotspots_i, false);
            save_output_nifti(fout, "hotspots_out", hotspots_o, false);
        }

        // --------------------------------------------------------------------
        // Compute equi-volume factors
        // --------------------------------------------------------------------
        cout << "\n  Start computing equi-volume factors..." << endl;
        nifti_image* equivol_factors = copy_nifti_as_float32(nii_rim);
        float* equivol_factors_data = static_cast<float*>(equivol_factors->data);
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(equivol_factors_data + i) = 0;
        }

        float w;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_rim_data + i) == 3) {
                // Find mass at each end of the given column
                j = *(innerGM_id_data + i);
                k = *(outerGM_id_data + i);
                if (*(curvature_data + i) == 0) {
                    w = 0.5;
                } else if (*(curvature_data + i) < 0) {
                    w = *(hotspots_i_data + k)
                        / (*(hotspots_i_data + k) + *(hotspots_o_data + k));
                } else if (*(curvature_data + i) > 0) {
                    w = *(hotspots_i_data + j)
                        / (*(hotspots_i_data + j) + *(hotspots_o_data + j));
                }
                *(equivol_factors_data + i) = w;
            }
        }

        if (mode_debug) {
            save_output_nifti(fout, "equivol_factors", equivol_factors, false);
        }

        // --------------------------------------------------------------------
        // Smooth equi-volume factors for seamless transitions
        // --------------------------------------------------------------------
        cout << "\n  Start smoothing transitions..." << endl;
        nifti_image* smooth = copy_nifti_as_float32(curvature);
        float* smooth_data = static_cast<float*>(smooth->data);
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(smooth_data + i) = 0;
        }

        // Pre-compute weights
        float FWHM_val = 1;  // TODO(Faruk): Might tweak this one
        float w_0 = gaus(0, FWHM_val);
        float w_dX = gaus(dX, FWHM_val);
        float w_dY = gaus(dY, FWHM_val);
        float w_dZ = gaus(dZ, FWHM_val);

        for (uint16_t n = 0; n != iter_smooth; ++n) {
            cout << "\r    Iteration: " << n+1 << "/" << iter_smooth << flush;
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                if (*(nii_rim_data + i) == 3) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    float new_val = 0, total_weight = 0;

                    // Start with the voxel itself
                    new_val += *(equivol_factors_data + i) * w_0;
                    total_weight += w_0;

                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            new_val += *(equivol_factors_data + j) * w_dX;
                            total_weight += w_dX;
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            new_val += *(equivol_factors_data + j) * w_dX;
                            total_weight += w_dX;
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            new_val += *(equivol_factors_data + j) * w_dY;
                            total_weight += w_dY;
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            new_val += *(equivol_factors_data + j) * w_dY;
                            total_weight += w_dY;
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            new_val += *(equivol_factors_data + j) * w_dZ;
                            total_weight += w_dZ;
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            new_val += *(equivol_factors_data + j) * w_dZ;
                            total_weight += w_dZ;
                        }
                    }
                    *(smooth_data + i) = new_val / total_weight;
                }
            }
            // Swap image data
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                *(equivol_factors_data + i) = *(smooth_data + i);
            }
        }
        cout << endl;
        if (mode_debug) {
            save_output_nifti(fout, "equivol_factors_smooth", smooth, false);
        }

        // --------------------------------------------------------------------
        // Apply equi-volume factors
        // --------------------------------------------------------------------
        cout << "\n  Start final layering..." << endl;
        float d1_new, d2_new, a, b;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_rim_data + i) == 3) {
                // Find normalized distances from a given point on a column
                float dist1 = *(innerGM_dist_data + i);
                float dist2 = *(outerGM_dist_data + i);
                float total_dist = dist1 + dist2;;
                dist1 /= total_dist;
                dist2 /= total_dist;

                a = *(equivol_factors_data + i);
                b = 1 - a;

                // Perturb using masses to modify distances in simplex space
                tie(d1_new, d2_new) = simplex_perturb_2D(dist1, dist2, a, b);

                // Difference of normalized distances (used in finding midGM)
                *(normdistdiff_data + i) = d1_new - d2_new;

                // Cast distances to integers as number of desired layers
                if (d1_new != 0 && isfinite(d1_new)) {
                    *(nii_layers_data + i) =  ceil(nr_layers * d1_new);
                } else {
                    *(nii_layers_data + i) = 1;
                }
            }
        }

        save_output_nifti(fout, "layers_equivol", nii_layers, true);

        // Save equi-volume metric in a simple 0-1 range form.
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_rim_data + i) == 3) {
                *(normdistdiff_data + i) /= 2;
                *(normdistdiff_data + i) += 0.5;
            }
        }
        save_output_nifti(fout, "metric_equivol", normdistdiff);

        if (mode_debug) {
            save_output_nifti(fout, "normdistdiff_equivol", normdistdiff, false);
        }

        // ====================================================================
        // Middle gray matter for equi-volume
        // ====================================================================
        cout << "\n  Start finding middle gray matter (equi-volume)..." << endl;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(midGM_data + i) = 0;
            *(midGM_id_data + i) = 0;
            // Change back normalized dist. differences after saves above
            if (*(nii_rim_data + i) == 3) {
                *(normdistdiff_data + i) -= 0.5;
                *(normdistdiff_data + i) *= 2;
            }
        }

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_rim_data + i) == 3) {
                // Check sign changes in normalized distance differences between
                // neighbouring voxels on a column path (a.k.a. streamline)
                if (*(normdistdiff_data + i) == 0) {
                    *(midGM_data + i) = 1;
                    *(midGM_id_data + i) = i;
                } else {
                    float m = *(normdistdiff_data + i);
                    float n;

                    // Inner neighbour
                    j = *(innerGM_prevstep_id_data + i);
                    if (*(nii_rim_data + j) == 3) {
                        n = *(normdistdiff_data + j);
                        if (signbit(m) - signbit(n) != 0) {
                            if (m*m < n*n) {
                                *(midGM_data + i) = 1;
                                *(midGM_id_data + i) = i;
                            } else if (m*m > n*n) {  // Closer to prev. step
                                *(midGM_data + j) = 1;
                                *(midGM_id_data + j) = j;
                            } else {  // Equal +/- normalized distance
                                *(midGM_data + i) = 1;
                                *(midGM_id_data + i) = i;
                                *(midGM_data + j) = 1;
                                *(midGM_id_data + j) = i;  // On purpose
                            }
                        }
                    }

                    // Outer neighbour
                    j = *(outerGM_prevstep_id_data + i);
                    if (*(nii_rim_data + j) == 3) {
                        n = *(normdistdiff_data + j);
                        if (signbit(m) - signbit(n) != 0) {
                            if (m*m < n*n) {
                                *(midGM_data + i) = 1;
                                *(midGM_id_data + i) = i;
                            } else if (m*m > n*n) {  // Closer to prev. step
                                *(midGM_data + j) = 1;
                                *(midGM_id_data + j) = j;
                            } else {  // Equal +/- normalized distance
                                *(midGM_data + i) = 1;
                                *(midGM_id_data + i) = i;
                                *(midGM_data + j) = 1;
                                *(midGM_id_data + j) = i;  // On purpose
                            }
                        }
                    }
                }
            }
        }
        save_output_nifti(fout, "midGM_equivol", midGM, true);
    }

    // ========================================================================
    // Cortical thickness
    // ========================================================================
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(innerGM_dist_data + i) += *(outerGM_dist_data + i);
    }
    save_output_nifti(fout, "thickness", innerGM_dist, true);
    // ========================================================================
    // Streamline angles
    // ========================================================================
    // NOTE(Faruk): I have implemented this quickly for B0 related before
    // ISMRM2021 abstract submission.

    // Prepare a 4D nifti for streamline vectors
    nifti_image* svec = nifti_copy_nim_info(normdist);
    svec->dim[0] = 4;  // For proper 4D nifti
    svec->dim[1] = size_x;
    svec->dim[2] = size_y;
    svec->dim[3] = size_z;
    svec->dim[4] = 3;
    nifti_update_dims_from_array(svec);
    svec->nvox = normdist->nvox * 3;
    svec->nbyper = sizeof(float);
    svec->data = calloc(svec->nvox, svec->nbyper);
    svec->scl_slope = 1;
    float* svec_data = static_cast<float*>(svec->data);

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            tie(wm_x, wm_y, wm_z) = ind2sub_3D(*(innerGM_id_data + i),
                                               size_x, size_y);
            tie(gm_x, gm_y, gm_z) = ind2sub_3D(*(outerGM_id_data + i),
                                               size_x, size_y);

            // Compute approx. surface angles from streamline start-end points
            float vec_x = wm_x - gm_x;
            float vec_y = wm_y - gm_y;
            float vec_z = wm_z - gm_z;
            float vec_norm = std::sqrt(vec_x * vec_x + vec_y * vec_y + vec_z * vec_z);

            *(svec_data + nr_voxels*0 + i) = vec_x / vec_norm;
            *(svec_data + nr_voxels*1 + i) = vec_y / vec_norm;
            *(svec_data + nr_voxels*2 + i) = vec_z / vec_norm;

            // Angular difference
            // float ref_x = 0, ref_y = 0, ref_z = 1;
            // float temp_dot = ref_x * vec_x + ref_y * vec_y + ref_z * vec_z;
            // float term1 = ref_x * ref_x + ref_y * ref_y + ref_z * ref_z;
            // float term2 = vec_x * vec_x + vec_y * vec_y + vec_z * vec_z;
            // float temp_angle = std::acos(temp_dot / std::sqrt(term1 * term2));
            // *(curvature_data + i) = temp_angle * 180 / PI;
            // ----------------------------------------------------------------
        }
    }
    save_output_nifti(fout, "streamline_vectors", svec, true);

    cout << "\n  Finished." << endl;
    return 0;
}
