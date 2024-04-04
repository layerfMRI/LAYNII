
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
    nifti_image* nii_rim = copy_nifti_as_int8(nii1);
    int8_t* nii_rim_data = static_cast<int8_t*>(nii_rim->data);

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
            ii += 1;
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
    // Prepare inverse mapping
    // ------------------------------------------------------------------------
    int* voi_id_inv = (int*)malloc(nr_voxels * sizeof(int));  // 3D
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        *(voi_id_inv + i) = ii;
    }

    // ------------------------------------------------------------------------
    // Reduce input to voxels of interest
    // ------------------------------------------------------------------------
    int8_t* voi_rim = (int8_t*)malloc(nr_voi * sizeof(int8_t));
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);
        *(voi_rim + ii) = *(nii_rim_data + i);
    }
    free(nii_rim_data);
    free(nii_rim);

    // // DEBUG
    // nifti_image* nii_out1 = copy_nifti_as_int8(nii1);
    // int8_t* nii_out1_data = static_cast<int8_t*>(nii_out1->data);
    // for (int i = 0; i != nr_voxels; ++i) {
    //     *(nii_out1_data + i) = 0;
    // }
    // for (int ii = 0; ii != nr_voi; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out1_data + i) = *(voi_rim + ii);
    // }
    // save_output_nifti(fout, "DEBUG_voi", nii_out1, false);
    // free(nii_out1);

    // ================================================================================================================
    // Grow from inner GM border
    // ================================================================================================================
    cout << "\n  Start growing from inner GM (WM-facing border)..." << endl;

    // Allocate memory
    int16_t* innerGM_step = (int16_t*)malloc(nr_voi * sizeof(int16_t));
    float* innerGM_dist = (float*)malloc(nr_voi * sizeof(float));

    // Initialize grow volume
    for (int ii = 0; ii != nr_voi; ++ii) {
        if (*(voi_rim + ii) == 2) {  // WM boundary voxels within GM
            *(innerGM_step + ii) = 1;
            *(innerGM_dist + ii) = 0.;
        } else {
            *(innerGM_step + ii) = 0;
            *(innerGM_dist + ii) = 0.;
        }
    }

    // Initialize unprocessed and processed sets
    int* vox_id_undetermined = (int*)malloc(nr_voi_gm * sizeof(int));
    int* temp = (int*)malloc(nr_voi_gm * sizeof(int));

    int size_undetermined = nr_voi_gm;
    int size_determined = 0;

    for (int i = 0; i != nr_voi_gm; ++i) {
        *(vox_id_undetermined + i) = *(voi_id + i);
    } 

    // Loop over gray matter voxels, compute distance if neighbor is border
    int step = 1;
    int jj;
    int counter = 666;
    while (counter != 0) {
    // for (int iter = 0; iter != 12; ++iter) {  // DEBUG
        cout << "    Processing voxels: " << size_undetermined << " -> " << size_determined << endl;

        // Loop over undetermined voxels only
        counter = 0;
        for (int ii = 0; ii != size_undetermined; ++ii) {
            int i = *(vox_id_undetermined + ii);
            int kk = *(voi_id_inv + i);
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
            bool border_detected = false;
            float d = std::numeric_limits<float>::max();

            // ----------------------------------------------------------------
            // 1-jump neighbours
            // ----------------------------------------------------------------
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dX;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                    *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dX;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                    *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dY;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;                }
            }

            j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dY;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                }
            }

            j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;                }
            }

            j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                }
            }

            // ----------------------------------------------------------------
            // 2-jump neighbours
            // ----------------------------------------------------------------
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XY;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XY;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XY;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XY;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_YZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_YZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_YZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_YZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            // ----------------------------------------------------------------
            // 3-jump neighbours
            // ----------------------------------------------------------------
            j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
            jj = *(voi_id_inv + j);
            if (*(innerGM_step + jj) == step) {
                border_detected = true;
                d = *(innerGM_dist + jj) + dia_XYZ;
                if (*(innerGM_dist + kk) == 0. || d < *(innerGM_dist + kk)) {
                *(innerGM_dist + kk) = d;
                }
            }

            // ----------------------------------------------------------------
            if (border_detected) {
                ++size_determined;
                ++counter;
                *(innerGM_step + kk) = step + 1;
            } else {
                *(vox_id_undetermined + ii-counter) = i;

            }
            // ----------------------------------------------------------------
        }
        ++step;
        size_undetermined = nr_voi_gm - size_determined;  // Update this outside
    }

    // // DEBUG
    // nifti_image* nii_out3 = copy_nifti_as_int16(nii1);
    // int16_t* nii_out3_data = static_cast<int16_t*>(nii_out3->data);
    // for (int i = 0; i != nr_voxels; ++i) {
    //     *(nii_out3_data + i) = 0;
    // }
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out3_data + i) = *(innerGM_step + ii);
    // }
    // save_output_nifti(fout, "DEBUG_innerGM_step", nii_out3, false);


    // DEBUG
    nifti_image* nii_out = copy_nifti_as_float32(nii1);
    float* nii_out_data = static_cast<float*>(nii_out->data);
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_out_data + i) = 0;
    }
    for (int ii = 0; ii != nr_voi_gm; ++ii) {
        int i = *(voi_id + ii);
        *(nii_out_data + i) = *(innerGM_dist + ii);
    }
    save_output_nifti(fout, "DEBUG_innerGM_dist", nii_out, false);

    // // ================================================================================================================
    // // Grow from outer GM border
    // // ================================================================================================================
    // cout << "\n  Start growing from outer GM..." << endl;

    // // Allocate memory
    // int16_t*    outerGM_step           = (int16_t*)malloc(nr_voi * sizeof(int16_t));
    // float*      outerGM_dist           = (float*)malloc(nr_voi * sizeof(float));

    // // Initialize grow volume
    // for (int ii = 0; ii != nr_voi; ++ii) {
    //     if (*(voi_rim + ii) == 1) {  // WM boundary voxels within GM
    //         *(outerGM_step + ii) = 1;
    //         *(outerGM_dist + ii) = 0.;
    //     } else {
    //         *(outerGM_step + ii) = 0.;
    //         *(outerGM_dist + ii) = 0.;
    //     }
    // }

    // // ------------------------------------------------------------------------
    // // Fill in first iteration of neighbors
    // // ------------------------------------------------------------------------
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     // --------------------------------------------------------------------
    //     // 1-jump neighbours
    //     // --------------------------------------------------------------------
    //     j = *(voi_id3D_n11 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dX; }

    //     j = *(voi_id3D_n12 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dX; }

    //     j = *(voi_id3D_n13 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dX; }

    //     j = *(voi_id3D_n14 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dX; }

    //     j = *(voi_id3D_n15 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dX; }

    //     j = *(voi_id3D_n16 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dX; }

    //     // --------------------------------------------------------------------
    //     // 2-jump neighbours
    //     // --------------------------------------------------------------------
    //     j = *(voi_id3D_n201 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xy; }

    //     j = *(voi_id3D_n202 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xy; }

    //     j = *(voi_id3D_n203 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xy; }

    //     j = *(voi_id3D_n204 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xy; }

    //     j = *(voi_id3D_n205 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_yz; }

    //     j = *(voi_id3D_n206 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_yz; }

    //     j = *(voi_id3D_n207 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_yz; }

    //     j = *(voi_id3D_n208 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_yz; }

    //     j = *(voi_id3D_n209 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xz; }

    //     j = *(voi_id3D_n210 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xz; }

    //     j = *(voi_id3D_n211 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xz; }

    //     j = *(voi_id3D_n212 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xz; }

    //     // --------------------------------------------------------------------
    //     // 3-jump neighbours
    //     // --------------------------------------------------------------------
    //     j = *(voi_id3D_n31 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    //     j = *(voi_id3D_n32 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    //     j = *(voi_id3D_n33 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    //     j = *(voi_id3D_n34 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    //     j = *(voi_id3D_n35 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    //     j = *(voi_id3D_n36 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    //     j = *(voi_id3D_n37 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    //     j = *(voi_id3D_n38 + ii); jj = *(voi_id_inv + j);
    //     if (*(voi_rim + jj) == 1) { *(outerGM_step + ii) = 2; *(outerGM_dist + ii) = dia_xyz; }

    // }

    // // ------------------------------------------------------------------------
    // // Fill rest of gray matter
    // // ------------------------------------------------------------------------
    // grow_step = 2;
    // voxel_counter = nr_voi_gm;
    // while (voxel_counter != 0) {
    //     voxel_counter = 0;
    //     for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //         if (*(outerGM_step + ii) == grow_step) {
    //             voxel_counter += 1;

    //             // ------------------------------------------------------------
    //             // 1-jump neighbours
    //             // ------------------------------------------------------------
    //             j = *(voi_id3D_n11 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dX;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n12 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dX;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n13 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dX;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n14 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dX;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n15 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dX;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n16 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dX;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             // ------------------------------------------------------------
    //             // 2-jump neighbours
    //             // ------------------------------------------------------------
    //             j = *(voi_id3D_n201 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xy;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n202 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xy;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n203 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xy;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n204 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xy;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n205 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_yz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n206 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_yz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n207 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_yz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n208 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_yz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n209 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n210 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n211 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n212 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             // ------------------------------------------------------------
    //             // 3-jump neighbours
    //             // ------------------------------------------------------------
    //             j = *(voi_id3D_n31 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n32 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n33 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n34 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n35 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n36 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n37 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }

    //             j = *(voi_id3D_n38 + ii);
    //             jj = *(voi_id_inv + j);
    //             d = *(outerGM_dist + ii) + dia_xyz;
    //             if (d < *(outerGM_dist + jj) || *(outerGM_step + jj) == 0) {
    //                 *(outerGM_dist + jj) = d;
    //                 *(outerGM_step + jj) = grow_step + 1;
    //             }
    //         }
    //     }
    //     grow_step += 1;
    // }


    // // DEBUG
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out_data + i) = *(outerGM_dist + ii);
    // }
    // save_output_nifti(fout, "DEBUG_outerGM_dist", nii_out, false);

    // // ================================================================================================================
    // // Compute Thickness
    // // ================================================================================================================
    // float* voi_thickness = (float*)malloc(nr_voi_gm * sizeof(float));
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     *(voi_thickness + ii) = *(innerGM_dist + ii) + *(outerGM_dist + ii);
    // }

    // // DEBUG
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out_data + i) = *(voi_thickness + ii);
    // }
    // save_output_nifti(fout, "DEBUG_thickness", nii_out, false);

    // // ================================================================================================================
    // // Compute Equidistant Metric
    // // ================================================================================================================
    // float* voi_metric = (float*)malloc(nr_voi_gm * sizeof(float));
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     *(voi_metric + ii) = *(innerGM_dist + ii) / *(voi_thickness + ii);
    // }

    // // DEBUG
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out_data + i) = *(voi_metric + ii);
    // }
    // save_output_nifti(fout, "DEBUG_equidist_metric", nii_out, false);

    // // ================================================================================================================
    // // Compute Equidistant Layers
    // // ================================================================================================================
    // int* voi_layers = (int*)malloc(nr_voi_gm * sizeof(int));
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     *(voi_layers + ii) = static_cast<int>(std::ceil(*(voi_metric + ii) * nr_layers));
    // }

    // // DEBUG
    // nifti_image* nii_out_int = copy_nifti_as_float32(nii_out);
    // float* nii_out_int_data = static_cast<float*>(nii_out_int->data);
    // for (int ii = 0; ii != nr_voi_gm; ++ii) {
    //     int i = *(voi_id + ii);
    //     *(nii_out_int_data + i) = *(voi_layers + ii);
    // }
    // save_output_nifti(fout, "DEBUG_equidist_layers", nii_out_int, false);


    cout << "\n  Finished.\n" << endl;
    return 0;
}
