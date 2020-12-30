
#include "../dep/laynii_lib.h"
#include <limits>

int show_help(void) {
    printf(
    "LN2_FLATTEN:\n"
    "\n"
    "!!! WORK IN PROGRESS !!!\n"
    "\n"
    "Usage:\n"
    "    LN2_FLATTEN -rim rim.nii -midgm rim_midgm_equidist.nii\n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -rim       : Segmentation input. Use 3 to code pure gray matter \n"
    "                 voxels. This program only generates columns in the \n"
    "                 voxels coded with 3.\n"
    "    -midgm     : Middle gray matter file (from LN2_LAYERS output).\n"
    "    -intensity : Intensity image that will be flattened.\n"
    "    -debug     : (Optional) Save extra intermediate outputs.\n"
    "    -output    : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2 = NULL, *fin3 = NULL;
    uint16_t ac;
    bool mode_debug = false;
    uint32_t q_factor = 10;

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
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-midgm")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -midgm\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-intensity")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -intensity\n");
                return 1;
            }
            fin3 = argv[ac];
        } else if (!strcmp(argv[ac], "-q_factor")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -q_factor\n");
            } else {
                q_factor = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-rim'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-midgm'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-intensity'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }
    nii2 = nifti_image_read(fin2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }
    nii3 = nifti_image_read(fin3, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin3);
        return 2;
    }

    log_welcome("LN2_FLATTEN");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);

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
    nifti_image* nii_midgm = copy_nifti_as_int16(nii2);
    int16_t* nii_midgm_data = static_cast<int16_t*>(nii_midgm->data);
    nifti_image* nii_intensity = copy_nifti_as_float32(nii3);
    float* nii_intensity_data = static_cast<float*>(nii_intensity->data);

    // Prepare required nifti images
    nifti_image* nii_columns  = copy_nifti_as_int32(nii_rim);
    int32_t* nii_columns_data = static_cast<int32_t*>(nii_columns->data);
    // Setting zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_columns_data + i) = 0;
    }

    nifti_image* flood_step = copy_nifti_as_int32(nii_columns);
    int32_t* flood_step_data = static_cast<int32_t*>(flood_step->data);
    nifti_image* flood_dist = copy_nifti_as_float32(nii_columns);
    float* flood_dist_data = static_cast<float*>(flood_dist->data);

    // Prepare 4D nifti output
    nifti_image* nii_coords = nifti_copy_nim_info(flood_dist);
    nii_coords->dim[0] = 4;  // For proper 4D nifti
    nii_coords->dim[1] = size_x;
    nii_coords->dim[2] = size_y;
    nii_coords->dim[3] = size_z;
    nii_coords->dim[4] = 2;  // X and Y coordinates
    nifti_update_dims_from_array(nii_coords);
    nii_coords->nvox = nr_voxels * 2;
    nii_coords->data = calloc(nii_coords->nvox, nii_coords->nbyper);
    float* nii_coords_data = static_cast<float*>(nii_coords->data);

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 1){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 1){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Find column centers through farthest flood distance
    // ========================================================================
    cout << "  Start finding origins..." << flush;
    // Find the initial voxel
    uint32_t start_voxel;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 1) {
            start_voxel = i;
        }
    }
    *(nii_midgm_data + start_voxel) = 2;

    // Initialize new voxel
    uint32_t new_voxel_id;

    // First and third points are enough to parametrize
    for (uint32_t n = 0; n != 3; ++n) {

        uint16_t grow_step = 1;
        uint32_t voxel_counter = nr_voxels;
        uint32_t ix, iy, iz, i, j;
        float d;

        // Initialize grow volume
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_midgm_data + i) == 2) {
                *(flood_step_data + i) = 1.;
                *(flood_dist_data + i) = 0.;
            } else {
                *(flood_step_data + i) = 0.;
                *(flood_dist_data + i) = 0.;
            }
        }

        while (voxel_counter != 0) {
            voxel_counter = 0;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                // Map subset to full set
                i = *(voi_id + ii);
                if (*(flood_step_data + i) == grow_step) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    voxel_counter += 1;

                    // ------------------------------------------------------------
                    // 1-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 2-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 3-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                }
            }
            grow_step += 1;
        }
        *(nii_columns_data + new_voxel_id) = n+1;
        *(nii_midgm_data + new_voxel_id) = 2;

        // Remove the initial voxel (reduces arbitrariness of the 1st point)
        // NOTE(Faruk): This step guarantees to start from extrememums. The
        // initial point is only used to determine an extremum distance.
        if (n == 0) {
            *(nii_midgm_data + start_voxel) = 1;
            // Also reset distances
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                i = *(voi_id + ii);
                *(flood_step_data + i) = 0.;
                *(flood_dist_data + i) = 0.;
            }
        }
    }
    cout << endl;

    if (mode_debug) {
        save_output_nifti(fout, "flood_step", flood_step, false);
        save_output_nifti(fout, "flood_dist", flood_dist, false);
    }
    // save_output_nifti(fout, "origins", nii_columns, true);

    // ========================================================================
    // Flood distance for parametrization
    // ========================================================================
    cout << "\n  Start parametrization..." << endl;

    for (uint32_t n = 0; n != 2; ++n) {

        uint16_t grow_step = 1;
        uint32_t voxel_counter = nr_voxels;
        uint32_t ix, iy, iz, i, j;
        int32_t k;
        float d, d_max;

        // Determine points that will be used as origins
        if (n == 0) {
            k = 1;
        } else {
            k = 3;
        }
        // Initialize grow volume
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_columns_data + i) == k) {
                *(flood_step_data + i) = 1.;
                *(flood_dist_data + i) = 0.;
            } else {
                *(flood_step_data + i) = 0.;
                *(flood_dist_data + i) = 0.;
            }
        }

        while (voxel_counter != 0) {
            voxel_counter = 0;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                // Map subset to full set
                i = *(voi_id + ii);
                if (*(flood_step_data + i) == grow_step) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    voxel_counter += 1;

                    // ------------------------------------------------------------
                    // 1-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 2-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 3-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                            }
                        }
                    }
                }
                if (d_max < d) {
                    d_max = d;
                }
            }
            grow_step += 1;
        }
        // Store distances as coordinates output
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_coords_data + n * nr_voxels + i) = *(flood_dist_data + i) / d_max;
        }
    }
    save_output_nifti(fout, "coords", nii_coords, true);

    // ========================================================================
    // Quantize
    // ========================================================================
    cout << "\n  Start quantization..." << endl;
    float coord_x, coord_y;
    int32_t q_coord_x, q_coord_y;
    uint32_t i;
    int32_t q_factor_int = static_cast<int32_t>(q_factor);
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);
        // Get coordinates
        coord_x = *(nii_coords_data + 0 * nr_voxels + i);
        coord_y = *(nii_coords_data + 1 * nr_voxels + i);
        // Quantize
        q_coord_x = static_cast<int32_t>(coord_x * q_factor);
        q_coord_y = static_cast<int32_t>(coord_y * q_factor);
        // Assign unique zone id
        *(nii_columns_data + i) = (q_coord_y * q_factor_int + q_coord_x);
    }
    save_output_nifti(fout, "coords_quantized", nii_columns, true);

    // ========================================================================
    // Generate flatmap nifti
    // ========================================================================
    cout << "  Start flattening..." << endl;

    // Prepare new output
    nifti_image* nii_flat = nifti_copy_nim_info(flood_dist);
    nii_flat->dim[0] = 3;
    nii_flat->dim[1] = q_factor_int;
    nii_flat->dim[2] = q_factor_int;
    nii_flat->dim[3] = 1;
    nifti_update_dims_from_array(nii_flat);
    nii_flat->nvox = q_factor_int * q_factor_int;
    nii_flat->data = calloc(nii_flat->nvox, nii_flat->nbyper);
    float* nii_flat_data = static_cast<float*>(nii_flat->data);

    for (uint32_t i = 0; i != q_factor_int * q_factor_int; ++i) {
        *(nii_flat_data + i) = 0;
    }

    // Prep required niftis
    nifti_image* nii_density = copy_nifti_as_float32(nii_flat);
    float* nii_density_data = static_cast<float*>(nii_density->data);

    // ------------------------------------------------------------------------
    // Transform 3D data onto flatmap
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);
        // Get coordinates
        coord_x = *(nii_coords_data + 0 * nr_voxels + i);
        coord_y = *(nii_coords_data + 1 * nr_voxels + i);
        // Quantize
        q_coord_x = static_cast<int32_t>(coord_x * q_factor);
        q_coord_y = static_cast<int32_t>(coord_y * q_factor);

        // *(nii_flat_data + q_factor_int * q_coord_y + q_coord_x) = *(nii_columns_data + i);
        *(nii_flat_data + q_factor_int * q_coord_y + q_coord_x) += *(nii_intensity_data + i);
        *(nii_density_data + q_factor_int * q_coord_y + q_coord_x) += 1;
    }

    // Take mean of pooled intensities
    for (uint32_t i = 0; i != q_factor_int * q_factor; ++i) {
        if (*(nii_density_data + i) > 0) {
            *(nii_flat_data + i) /= *(nii_density_data + i);
        }
    }

    save_output_nifti(fout, "density", nii_density, false);
    save_output_nifti(fout, "flatmap", nii_flat, true);

    // // ========================================================================
    // // Voronoi cell from MidGM cells to rest of the GM (gray matter)
    // // ========================================================================
    // cout << "\n  Start Voronoi..." << endl;
    //
    // // ------------------------------------------------------------------------
    // // Reduce number of looped-through voxels
    // // TODO[Faruk]: Put this into a function to reduce code repetition.
    // nr_voi = 0;  // Voxels of interest
    // for (uint32_t i = 0; i != nr_voxels; ++i) {
    //     if (*(nii_rim_data + i) == 3){
    //         nr_voi += 1;
    //     }
    // }
    // // Allocate memory to only the voxel of interest
    // free(voi_id);
    // voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));
    //
    // // Fill in indices to be able to remap from subset to full set of voxels
    // ii = 0;
    // for (uint32_t i = 0; i != nr_voxels; ++i) {
    //     if (*(nii_rim_data + i) == 3){
    //         *(voi_id + ii) = i;
    //         ii += 1;
    //     }
    // }
    // // ------------------------------------------------------------------------
    //
    // // Initialize grow volume
    // for (uint32_t i = 0; i != nr_voxels; ++i) {
    //     if (*(nii_columns_data + i) != 0) {
    //         *(flood_step_data + i) = 1.;
    //         *(flood_dist_data + i) = 0.;
    //     } else {
    //         *(flood_step_data + i) = 0.;
    //         *(flood_dist_data + i) = 0.;
    //     }
    // }
    //
    // uint16_t grow_step = 1;
    // uint32_t voxel_counter = nr_voxels;
    // uint32_t ix, iy, iz, j;
    // float d;
    // voxel_counter = nr_voxels;
    // while (voxel_counter != 0) {
    //     voxel_counter = 0;
    //     for (uint32_t ii = 0; ii != nr_voi; ++ii) {
    //         i = *(voi_id + ii);
    //         if (*(flood_step_data + i) == grow_step) {
    //             tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
    //             voxel_counter += 1;
    //
    //             // ------------------------------------------------------------
    //             // 1-jump neighbours
    //             // ------------------------------------------------------------
    //             if (ix > 0) {
    //                 j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dX;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);                        }
    //                 }
    //             }
    //             if (ix < end_x) {
    //                 j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dX;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //
    //                     }
    //                 }
    //             }
    //             if (iy > 0) {
    //                 j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dY;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (iy < end_y) {
    //                 j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dY;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (iz > 0) {
    //                 j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dZ;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (iz < end_z) {
    //                 j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dZ;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //
    //             // ------------------------------------------------------------
    //             // 2-jump neighbours
    //             // ------------------------------------------------------------
    //             if (ix > 0 && iy > 0) {
    //                 j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xy;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix > 0 && iy < end_y) {
    //                 j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xy;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iy > 0) {
    //                 j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xy;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iy < end_y) {
    //                 j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xy;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (iy > 0 && iz > 0) {
    //                 j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_yz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (iy > 0 && iz < end_z) {
    //                 j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_yz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (iy < end_y && iz > 0) {
    //                 j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_yz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (iy < end_y && iz < end_z) {
    //                 j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_yz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix > 0 && iz > 0) {
    //                 j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iz > 0) {
    //                 j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix > 0 && iz < end_z) {
    //                 j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iz < end_z) {
    //                 j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //
    //             // ------------------------------------------------------------
    //             // 3-jump neighbours
    //             // ------------------------------------------------------------
    //             if (ix > 0 && iy > 0 && iz > 0) {
    //                 j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix > 0 && iy > 0 && iz < end_z) {
    //                 j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix > 0 && iy < end_y && iz > 0) {
    //                 j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iy > 0 && iz > 0) {
    //                 j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix > 0 && iy < end_y && iz < end_z) {
    //                 j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iy > 0 && iz < end_z) {
    //                 j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iy < end_y && iz > 0) {
    //                 j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //             if (ix < end_x && iy < end_y && iz < end_z) {
    //                 j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
    //
    //                 if (*(nii_rim_data + j) == 3) {
    //                     d = *(flood_dist_data + i) + dia_xyz;
    //                     if (d < *(flood_dist_data + j)
    //                         || *(flood_dist_data + j) == 0) {
    //                         *(flood_dist_data + j) = d;
    //                         *(flood_step_data + j) = grow_step + 1;
    //                         *(nii_columns_data + j) = *(nii_columns_data + i);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     grow_step += 1;
    // }
    //
    // // ========================================================================
    // save_output_nifti(fout, "cells", nii_columns, true);
    // if (mode_debug) {
    //     save_output_nifti(fout, "voronoi_flood_step", flood_step, false);
    //     save_output_nifti(fout, "voronoi_flood_dist", flood_dist, false);
    // }

    cout << "\n  Finished." << endl;
    return 0;
}
