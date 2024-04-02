
#include "../dep/laynii_lib.h"
#include <limits>
#include <set>


int show_help(void) {
    printf(
    "LN3_LAYERS: !!!EXPERIMENTAL - OPTIMUM RAM USAGE TESTS!!!"
    "\n"
    "Usage:\n"
    "    LN3_LAYERS -rim rim.nii\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -rim          : A segmented image. This image must use 1 to code outer\n"
    "                    gray matter border voxels (facing mostly CSF), 2 to code inner\n"
    "                    gray matter border voxels (facing mostly white matter), and\n"
    "                    3 to code pure gray matter voxels.\n"
    "    -nr_layers    : Number of layers. Default is 3.\n"
    "    -equivol      : (Optional) Create equi-volume layers. We do not\n"
    "                    recommend this option if your rim file is above 0.3mm\n"
    "                    resolution. You can always upsample your rim file to\n"
    "                    a higher resolution first (<0.3mm) and then use this\n"
    "                    option.\n"
    "    -iter_smooth  : (Optional) Number of smoothing iterations. Default\n"
    "                    is 100. Only used together with '-equivol' flag. Use\n"
    "                    larger values when equi-volume layers are jagged.\n"
    "    -curvature    : (Optional) Compute curvature. Uses -iter_smooth value\n"
    "                    for smoothing the curvature estimates. Off by default.\n"
    "    -streamlines  : (Optional) Export streamline vectors. Useful for e.g.\n"
    "                    computing B0 angular differences. Off by default.\n"
    "    -thickness    : (Optional) Export cortical thickness. Uses -iter_smooth\n"
    "                    value for smoothing the thickness. Off by default.\n"
    "    -incl_borders : (Optional) Include inner and outer gray matter borders\n"
    "                    into the layering. This treats the borders as \n"
    "                    a part of gray matter. Off by default.\n"
    "    -equal_counts : (Optional) Equalize number of voxels for each layer.\n"
    "                    This option inherently includes the borders.\n"
    "                    output is given with file name addition `*layers_equicount*.\n"
    "                    Useful for ~0.8 mm inputs where no upsampling is done.\n"
    "    -no_smooth    : (Optional) Disable smoothing on cortical depth metric.\n"
    "    -debug        : (Optional) Save extra intermediate outputs.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin = NULL, *fout = NULL;
    uint16_t ac, nr_layers = 3;
    uint16_t iter_smooth = 100;
    bool mode_equivol = false, mode_debug = false, mode_incl_borders = false;
    bool mode_curvature =false, mode_streamlines = false, mode_smooth = true;
    bool mode_thickness = false, mode_equal_counts = false;

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

    log_welcome("LN3_LAYERS");
    log_nifti_descriptives(nii1);

    cout << "  Nr. layers: " << nr_layers << endl;

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
    // ========================================================================
    nifti_image* nii_rim = copy_nifti_as_int8(nii1);
    int8_t* nii_rim_data = static_cast<int8_t*>(nii_rim->data);
    free(nii1);

    // ------------------------------------------------------------------------
    // Prepare for RAM optimization RAM by allocating minimal data
    // ------------------------------------------------------------------------
    // First find all gray matter voxels
    int nr_voi_gm = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {  // Only gray matter voxels
            nr_voi_gm += 1;
        }
    }

    // Allocate memory to only the voxel of interest
    int* voi_id_gm = (int*) malloc(nr_voi_gm*sizeof(int));

    // Fill in indices to be able to remap from subset to full set of voxels
    int ii = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            *(voi_id_gm + ii) = i;
            ii += 1;
        }
    }

    // Second, find all gray matter boundaries
    std::set<int> voi_id_borders;  // Need sets as j's can overlap over loop
    int ix, iy, iz, j;
    for (int ii = 0; ii != nr_voi_gm; ++ii) {
        int i = *(voi_id_gm + ii);
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

        // --------------------------------------------------------------------
        // 1-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            if (*(nii_rim_data + j) == 1 || *(nii_rim_data + j) == 2) {
                voi_id_borders.insert(j);
            }
        }
        if (ix < end_x) {
            j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            if (*(nii_rim_data + j) == 1 || *(nii_rim_data + j) == 2) {
                voi_id_borders.insert(j);
            }
        }
        if (iy > 0) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            if (*(nii_rim_data + j) == 1 || *(nii_rim_data + j) == 2) {
                voi_id_borders.insert(j);
            }
        }
        if (iy < end_y) {
            j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            if (*(nii_rim_data + j) == 1 || *(nii_rim_data + j) == 2) {
                voi_id_borders.insert(j);
            }
        }
        if (iz > 0) {
            j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
            if (*(nii_rim_data + j) == 1 || *(nii_rim_data + j) == 2) {
                voi_id_borders.insert(j);
            }
        }
        if (iz < end_z) {
            j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
            if (*(nii_rim_data + j) == 1 || *(nii_rim_data + j) == 2) {
                voi_id_borders.insert(j);
            }
        }
    }
    int nr_voi_borders = voi_id_borders.size();

    // ------------------------------------------------------------------------
    // Allocate memory to gray matter + border voxels of interest
    int nr_voi = (nr_voi_gm + nr_voi_borders);
    int* voi_id = (int*)malloc(nr_voi * sizeof(int));

    // Merge gray matter and border voxel indices
    for (int i = 0; i != nr_voi_gm; ++i) {
        *(voi_id + i) = *(voi_id_gm + i);
    }
    int temp_counter = 0;
    for (const int& i : voi_id_borders) {
        *(voi_id + nr_voi_gm + temp_counter) = i;
        temp_counter += 1;
    }

    // ------------------------------------------------------------------------
    printf("  Sparsity: %.1f %% (%.2fM out of %.2fM voxels are gray matter + borders)", 
           static_cast<float>(nr_voi) / nr_voxels * 100, 
           static_cast<float>(nr_voi) / 1000000, 
           static_cast<float>(nr_voxels) / 1000000
           );

    // ------------------------------------------------------------------------
    // Reduce input to voxels of interest
    // ------------------------------------------------------------------------
    int8_t* voi_rim = (int8_t*)malloc(nr_voi * sizeof(int8_t));
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        *(voi_rim + ii) = *(nii_rim_data + i);
    }

    // DEBUG
    nifti_image* nii_out = copy_nifti_as_int8(nii1);
    int8_t* nii_out_data = static_cast<int8_t*>(nii_out->data);
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_out_data + i) = 0;
    }
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        *(nii_out_data + i) = *(voi_rim + ii);
    }
    save_output_nifti(fout, "DEBUG_voi", nii_out, false);

    // // ------------------------------------------------------------------------
    // // Allocate memory for neighbor indices
    // // ------------------------------------------------------------------------
    // int* voi_id_11 = (int*)malloc(nr_voi * sizeof(int));
    // int* voi_id_12 = (int*)malloc(nr_voi * sizeof(int));
    // int* voi_id_13 = (int*)malloc(nr_voi * sizeof(int));
    // int* voi_id_14 = (int*)malloc(nr_voi * sizeof(int));
    // int* voi_id_15 = (int*)malloc(nr_voi * sizeof(int));
    // int* voi_id_16 = (int*)malloc(nr_voi * sizeof(int));

    // // Populate neighbors
    // int ix, iy, iz, j, k;
    // for (int ii = 0; ii != nr_voi; ++ii) {
    //     int i = *(voi_id + ii);
    //     tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

    //     // --------------------------------------------------------------------
    //     // 1-jump neighbours
    //     // --------------------------------------------------------------------
    //     if (ix > 0) {
    //         j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
    //         *(voi_id_11 + ii) = j;
    //     }
    //     if (ix < end_x) {
    //         j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
    //         *(voi_id_12 + ii) = j;
    //     }
    //     if (iy > 0) {
    //         j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
    //         *(voi_id_13 + ii) = j;
    //     }
    //     if (iy < end_y) {
    //         j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
    //         *(voi_id_14 + ii) = j;
    //     }
    //     if (iz > 0) {
    //         j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
    //         *(voi_id_15 + ii) = j;
    //     }
    //     if (iz < end_z) {
    //         j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
    //         *(voi_id_16 + ii) = j;
    //     }
    // }

    // // ------------------------------------------------------------------------
    // // Allocate memory for only voxels of interest
    // int16_t* innerGM_step_data;
    // innerGM_step_data = (int16_t*) malloc(nr_voi * sizeof(int16_t));

    // float* innerGM_dist_data;
    // innerGM_dist_data = (float*) malloc(nr_voi * sizeof(float));

    // int32_t* innerGM_id_data;
    // innerGM_id_data = (int32_t*) malloc(nr_voi * sizeof(int32_t));

    // // ========================================================================
    // // Grow from WM
    // // ========================================================================
    // cout << "\n  Start growing from inner GM (WM-facing border)..." << endl;

    // // Initialize grow volume
    // for (int ii = 0; ii != nr_voi; ++ii) {
    //     int i = *(voi_id + ii);

    //     if (*(nii_rim_data + i) == 2) {  // WM boundary voxels within GM
    //         *(innerGM_step_data + ii) = 1;
    //         *(innerGM_dist_data + ii) = 0.;
    //         *(innerGM_id_data + ii) = i;
    //     } else {
    //         *(innerGM_step_data + ii) = 0.;
    //         *(innerGM_dist_data + ii) = 0.;
    //     }
    // }

    // uint16_t grow_step = 1;
    // uint32_t voxel_counter = nr_voxels;
    // float d;
    // while (voxel_counter != 0) {
    //     voxel_counter = 0;
    //     for (int ii = 0; ii != nr_voi; ++ii) {
    //         int i = *(voi_id + ii);

    //         if (*(innerGM_step_data + ii) == grow_step) {
    //             voxel_counter += 1;

    //             // ------------------------------------------------------------
    //             // 1-jump neighbours
    //             // ------------------------------------------------------------
    //             if (ix > 0) {
    //                 j = *(voi_id_11 + ii);
    //                 if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //                     d = *(innerGM_dist_data + ii) + dX;
    //                     if (d < *(innerGM_dist_data + j)
    //                         || *(innerGM_dist_data + j) == 0) {
    //                         *(innerGM_dist_data + j) = d;
    //                         *(innerGM_step_data + j) = grow_step + 1;
    //                         *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //                         *(innerGM_prevstep_id_data + j) = i;
    //                     }
    //                 }
    //             }
    //             if (ix < end_x) {
    //                 j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //                     d = *(innerGM_dist_data + i) + dX;
    //                     if (d < *(innerGM_dist_data + j)
    //                         || *(innerGM_dist_data + j) == 0) {
    //                         *(innerGM_dist_data + j) = d;
    //                         *(innerGM_step_data + j) = grow_step + 1;
    //                         *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //                         *(innerGM_prevstep_id_data + j) = i;
    //                     }
    //                 }
    //             }
    //             if (iy > 0) {
    //                 j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //                     d = *(innerGM_dist_data + i) + dY;
    //                     if (d < *(innerGM_dist_data + j)
    //                         || *(innerGM_dist_data + j) == 0) {
    //                         *(innerGM_dist_data + j) = d;
    //                         *(innerGM_step_data + j) = grow_step + 1;
    //                         *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //                         *(innerGM_prevstep_id_data + j) = i;
    //                     }
    //                 }
    //             }
    //             if (iy < end_y) {
    //                 j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //                     d = *(innerGM_dist_data + i) + dY;
    //                     if (d < *(innerGM_dist_data + j)
    //                         || *(innerGM_dist_data + j) == 0) {
    //                         *(innerGM_dist_data + j) = d;
    //                         *(innerGM_step_data + j) = grow_step + 1;
    //                         *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //                         *(innerGM_prevstep_id_data + j) = i;
    //                     }
    //                 }
    //             }
    //             if (iz > 0) {
    //                 j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //                     d = *(innerGM_dist_data + i) + dZ;
    //                     if (d < *(innerGM_dist_data + j)
    //                         || *(innerGM_dist_data + j) == 0) {
    //                         *(innerGM_dist_data + j) = d;
    //                         *(innerGM_step_data + j) = grow_step + 1;
    //                         *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //                         *(innerGM_prevstep_id_data + j) = i;
    //                     }
    //                 }
    //             }
    //             if (iz < end_z) {
    //                 j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

    //                 if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //                     d = *(innerGM_dist_data + i) + dZ;
    //                     if (d < *(innerGM_dist_data + j)
    //                         || *(innerGM_dist_data + j) == 0) {
    //                         *(innerGM_dist_data + j) = d;
    //                         *(innerGM_step_data + j) = grow_step + 1;
    //                         *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //                         *(innerGM_prevstep_id_data + j) = i;
    //                     }
    //                 }
    //             }

    //             // // ------------------------------------------------------------
    //             // // 2-jump neighbours
    //             // // ------------------------------------------------------------
    //             // if (ix > 0 && iy > 0) {
    //             //     j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xy;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix > 0 && iy < end_y) {
    //             //     j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xy;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iy > 0) {
    //             //     j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xy;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iy < end_y) {
    //             //     j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xy;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (iy > 0 && iz > 0) {
    //             //     j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_yz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (iy > 0 && iz < end_z) {
    //             //     j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_yz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (iy < end_y && iz > 0) {
    //             //     j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_yz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (iy < end_y && iz < end_z) {
    //             //     j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_yz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix > 0 && iz > 0) {
    //             //     j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iz > 0) {
    //             //     j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix > 0 && iz < end_z) {
    //             //     j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iz < end_z) {
    //             //     j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }

    //             // // ------------------------------------------------------------
    //             // // 3-jump neighbours
    //             // // ------------------------------------------------------------
    //             // if (ix > 0 && iy > 0 && iz > 0) {
    //             //     j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix > 0 && iy > 0 && iz < end_z) {
    //             //     j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix > 0 && iy < end_y && iz > 0) {
    //             //     j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iy > 0 && iz > 0) {
    //             //     j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix > 0 && iy < end_y && iz < end_z) {
    //             //     j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iy > 0 && iz < end_z) {
    //             //     j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iy < end_y && iz > 0) {
    //             //     j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //             // if (ix < end_x && iy < end_y && iz < end_z) {
    //             //     j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

    //             //     if (*(nii_rim_data + j) == 3 || *(nii_rim_data + j) == 1 ) {
    //             //         d = *(innerGM_dist_data + i) + dia_xyz;
    //             //         if (d < *(innerGM_dist_data + j)
    //             //             || *(innerGM_dist_data + j) == 0) {
    //             //             *(innerGM_dist_data + j) = d;
    //             //             *(innerGM_step_data + j) = grow_step + 1;
    //             //             *(innerGM_id_data + j) = *(innerGM_id_data + i);
    //             //             *(innerGM_prevstep_id_data + j) = i;
    //             //         }
    //             //     }
    //             // }
    //         }
    //     }
    //     grow_step += 1;
    // }

    // if (mode_debug) {
    //     save_output_nifti(fout, "innerGM_step", innerGM_step, false);
    //     save_output_nifti(fout, "innerGM_dist", innerGM_dist, false);
    //     save_output_nifti(fout, "innerGM_id", innerGM_id, false);
    // }

    cout << "\n  Finished." << endl;
    return 0;
}
