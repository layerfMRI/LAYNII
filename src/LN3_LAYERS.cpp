
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
    const float dia_XY = sqrt(dX * dX + dY * dY);
    const float dia_XZ = sqrt(dX * dX + dZ * dZ);
    const float dia_YZ = sqrt(dY * dY + dZ * dZ);
    // Long diagonals
    const float dia_XYZ = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ------------------------------------------------------------------------
    // Fix input datatype issues
    nifti_image* nii_rim = copy_nifti_as_int32(nii1);
    int* nii_rim_data = static_cast<int*>(nii_rim->data);

    // ================================================================================================================
    // Prepare for RAM optimization RAM by allocating minimal data
    // ================================================================================================================
    cout << "\n  Start memory optimization..." << endl;

    // First find all gray matter voxels
    int nr_voi_gm = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {  // Only gray matter voxels
            nr_voi_gm += 1;
        }
    }

    // Allocate memory to only the voxel of interest
    int* voi_id_gm = (int*)malloc(nr_voi_gm*sizeof(int));

    // Fill in indices to be able to remap from subset to full set of voxels
    int ii = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            *(voi_id_gm + ii) = i;
            ++ii;
        }
    }

    // Copy rim file
    int8_t* temp_rim = (int8_t*)malloc(nr_voxels*sizeof(int8_t));
    for (int i = 0; i != nr_voxels; ++i) {
        *(temp_rim + i) = *(nii_rim_data + i);
    }

    // Allocate memory to non-unique border voxels (assumes that border voxels 
    // can not be more than twice of gray matter voxels)
    int* voi_id_borders = (int*)malloc(nr_voi_gm*2*sizeof(int));

    // Loop through gm voxels, find neighbors, and erase each encounter
    int nr_voi_borders = 0;
    int ix, iy, iz, j;
    int k = 0;
    for (int ii = 0; ii != nr_voi_gm; ++ii) {
        int i = *(voi_id_gm + ii);
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

        // --------------------------------------------------------------------
        // 1-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x) {
            j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iy > 0) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iy < end_y) {
            j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iz > 0) {
            j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iz < end_z) {
            j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        // --------------------------------------------------------------------
        // 2-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && iy > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix > 0 && iy < end_y) {
            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iy > 0) {
            j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iy < end_y) {
            j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iy > 0 && iz > 0) {
            j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iy < end_y && iz > 0) {
            j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iz > 0) {
            j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix > 0 && iz < end_z) {
            j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iz < end_z) {
            j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        // --------------------------------------------------------------------
        // 3-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && iy > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix > 0 && iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix > 0 && iy < end_y && iz > 0) {
            j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iy > 0 && iz > 0) {
            j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix > 0 && iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iy < end_y && iz > 0) {
            j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
        if (ix < end_x && iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
            if (*(temp_rim + j) == 1 || *(temp_rim + j) == 2) {
                nr_voi_borders += 1; *(voi_id_borders + k) = j; ++k; *(temp_rim + j) = 0;
            }
        }
    }
    free(temp_rim);

    int nr_voi = nr_voi_gm + nr_voi_borders;
    printf("    Nr. voxels of gray matter : %d\n", nr_voi_gm);
    printf("    Nr. voxels of borders     : %d\n", nr_voi_borders);
    printf("    Nr. voxels of all         : %d\n", nr_voxels);
    printf("    Sparsity : %.1f %% (%.2fM out of %.2fM voxels are gray matter + borders)\n", 
           static_cast<float>(nr_voi) / nr_voxels * 100, 
           static_cast<float>(nr_voi) / 1000000, 
           static_cast<float>(nr_voxels) / 1000000
           );

    // ------------------------------------------------------------------------
    // Allocate memory to gray matter + border voxels of interest
    int* voi_id = (int*)malloc(nr_voi * sizeof(int));

    // Merge gray matter and border voxel indices
    for (int i = 0; i != nr_voi_gm; ++i) {
        *(voi_id + i) = *(voi_id_gm + i);
    }
    for (int i = 0; i != nr_voi_borders; ++i) {
        *(voi_id + nr_voi_gm + i) = *(voi_id_borders + i);
    }
    free(voi_id_borders);
    free(voi_id_gm);

    // ------------------------------------------------------------------------
    // Reduce input to voxels of interest
    // ------------------------------------------------------------------------
    int8_t* voi_rim = (int8_t*)malloc((nr_voi+1) * sizeof(int8_t));  // +1 is very important to hold references to 0
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        *(voi_rim + ii) = *(nii_rim_data + i);
    }
    *(voi_rim + nr_voi) = 0;  // Reserve last data point to 0 label in rim (very important)
    free(nii_rim_data);
    free(nii_rim);

    // // DEBUG
    // nifti_image* nii_out1 = copy_nifti_as_int32(nii1);
    // int* nii_out1_data = static_cast<int*>(nii_out1->data);
    // for (int i = 0; i != nr_voxels; ++i) {
    //     *(nii_out1_data + i) = 0;
    // }
    // for (int ii = 0; ii != nr_voi; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out1_data + i) = *(voi_rim + ii);
    // }
    // save_output_nifti(fout, "DEBUG_voi", nii_out1, false);
    // free(nii_out1);

    // ------------------------------------------------------------------------
    // Prepare inverse mapping
    // ------------------------------------------------------------------------
    int* voi_id_inv = (int*)malloc(nr_voxels * sizeof(int));  // 3D
    // All voxels point to 0 in voi_id by default if the following assignment is not done. Therefore, I decided to
    // allocate one extra voxels in voi_id to point to 0.
    for (int i = 0; i != nr_voxels; ++i) {
        *(voi_id_inv + i) = nr_voi + 1;  
    }
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        *(voi_id_inv + i) = ii;
    }

    // ------------------------------------------------------------------------
    // Prepare output nifti
    // ------------------------------------------------------------------------
    nifti_image* nii_out_float32 = copy_nifti_as_float32(nii1);
    float* nii_out_float32_data = static_cast<float*>(nii_out_float32->data);
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_out_float32_data + i) = 0;
    }

    // ================================================================================================================
    // Grow from inner GM border
    // ================================================================================================================
    printf("\n  Start growing from inner GM border (facing WM)...\n");

    // Allocate memory
    int16_t* voi_step = (int16_t*)malloc(nr_voi * sizeof(int16_t));
    float* innerGM_dist = (float*)malloc(nr_voi * sizeof(float));

    // Initialize grow volume
    for (int ii = 0; ii != nr_voi; ++ii) {
        if (*(voi_rim + ii) == 2) {  // WM boundary voxels within GM
            *(voi_step + ii) = 1;
            *(innerGM_dist + ii) = 0.;
        } else {
            *(voi_step + ii) = 0;
            *(innerGM_dist + ii) = 0.;
        }
    }

    // Initialize unprocessed and processed sets
    int* vox_id_undetermined = (int*)malloc(nr_voi * sizeof(int));
    int size_undetermined = nr_voi, size_determined = 0;
    for (int i = 0; i != nr_voi; ++i) {
        *(vox_id_undetermined + i) = *(voi_id + i);
    } 

    // Loop over gray matter voxels, compute distance if neighbor is border
    printf("    Processing voxels [undetermined -> determined]\n");
    int step = 1;
    int jj;
    int temp_counter = -1;
    while (temp_counter != 0) {
        std::cout << "\r      Step: " << step << " | " << size_undetermined << " -> " << size_determined << "        ";
        std::cout.flush();

        // Loop over undetermined voxels only
        temp_counter = 0;
        for (int ii = 0; ii != size_undetermined; ++ii) {
            int i = *(vox_id_undetermined + ii);
            int kk = *(voi_id_inv + i);

            if (*(voi_step + kk) == step) {

                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                bool border_detected = false;
                float d = std::numeric_limits<float>::max();

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dX;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dX;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dY;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dY;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XY;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XY;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XY;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XY;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_YZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_YZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_YZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_YZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 1 ) {
                        d = *(innerGM_dist + kk) + dia_XYZ;
                        if (d < *(innerGM_dist + jj)
                            || *(innerGM_dist + jj) == 0) {
                            *(innerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }

                ++size_determined;
                ++temp_counter;
            } else {
                *(vox_id_undetermined + ii - temp_counter) = i;
            }
        }
        size_undetermined = nr_voi - size_determined;
        ++step;
    }
    std::cout << std::endl;
    printf("    Processed all voxels [%d -> %d]\n", size_undetermined, size_determined);

    // // DEBUG
    // for (int ii = 0; ii != nr_voi; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out_float32_data + i) = *(innerGM_dist + ii);
    // }
    // save_output_nifti(fout, "DEBUG_innerGM_dist", nii_out_float32, false);

    // ================================================================================================================
    // Grow from outer GM border
    // ================================================================================================================
    printf("\n  Start growing from outer GM border...\n");

    // Allocate memory
    float* outerGM_dist = (float*)malloc(nr_voi * sizeof(float));

    // Initialize grow volume
    for (int ii = 0; ii != nr_voi; ++ii) {
        if (*(voi_rim + ii) == 1) {  // WM boundary voxels within GM
            *(voi_step + ii) = 1;
            *(outerGM_dist + ii) = 0.;
        } else {
            *(voi_step + ii) = 0;
            *(outerGM_dist + ii) = 0.;
        }
    }

    // Initialize unprocessed and processed sets
    size_undetermined = nr_voi, size_determined = 0;
    for (int i = 0; i != nr_voi; ++i) {
        *(vox_id_undetermined + i) = *(voi_id + i);
    } 

    // Loop over gray matter voxels, compute distance if neighbor is border
    printf("    Processing voxels [undetermined -> determined]\n");
    step = 1;
    temp_counter = -1;
    while (temp_counter != 0) {
        std::cout << "\r      Step: " << step << " | " << size_undetermined << " -> " << size_determined << "        ";
        std::cout.flush();

        // Loop over undetermined voxels only
        temp_counter = 0;
        for (int ii = 0; ii != size_undetermined; ++ii) {
            int i = *(vox_id_undetermined + ii);
            int kk = *(voi_id_inv + i);

            if (*(voi_step + kk) == step) {

                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                bool border_detected = false;
                float d = std::numeric_limits<float>::max();

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dX;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dX;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dY;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dY;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XY;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XY;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XY;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XY;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_YZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_YZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_YZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_YZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    jj = *(voi_id_inv + j);
                    if (*(voi_rim + jj) == 3 || *(voi_rim + jj) == 2) {
                        d = *(outerGM_dist + kk) + dia_XYZ;
                        if (d < *(outerGM_dist + jj)
                            || *(outerGM_dist + jj) == 0) {
                            *(outerGM_dist + jj) = d;
                            *(voi_step + jj) = step + 1;
                        }
                    }
                }

                ++size_determined;
                ++temp_counter;
            } else {
                *(vox_id_undetermined + ii - temp_counter) = i;
            }
        }
        size_undetermined = nr_voi - size_determined;
        ++step;
    }
    std::cout << std::endl;
    printf("    Processed all voxels [%d -> %d]\n", size_undetermined, size_determined);

    // // DEBUG
    // for (int ii = 0; ii != nr_voi; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out_float32_data + i) = *(outerGM_dist + ii);
    // }
    // save_output_nifti(fout, "DEBUG_outerGM_dist", nii_out_float32, false);

    // ================================================================================================================
    // Compute Thickness
    // ================================================================================================================
    float* voi_thickness = (float*)malloc(nr_voi * sizeof(float));
    for (int ii = 0; ii != nr_voi; ++ii) {
        *(voi_thickness + ii) = *(innerGM_dist + ii) + *(outerGM_dist + ii);
    }

    // DEBUG
    for (int ii = 0; ii != nr_voi_gm; ++ii) {
        int i = *(voi_id + ii);
        *(nii_out_float32_data + i) = *(voi_thickness + ii);
    }
    save_output_nifti(fout, "DEBUG_thickness", nii_out_float32, false);

    // ================================================================================================================
    // Compute Equidistant Metric
    // ================================================================================================================
    printf("  Saving metric...\n");

    float* voi_metric = (float*)malloc(nr_voi * sizeof(float));
    for (int ii = 0; ii != nr_voi; ++ii) {
        *(voi_metric + ii) = *(innerGM_dist + ii) / *(voi_thickness + ii);
    }

    // DEBUG
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        *(nii_out_float32_data + i) = *(voi_metric + ii);
    }
    save_output_nifti(fout, "DEBUG_equidist_metric", nii_out_float32, false);

    // ================================================================================================================
    // Compute Equidistant Layers
    // ================================================================================================================

    // TODO

    cout << "\n  Finished.\n" << endl;
    return 0;
}
