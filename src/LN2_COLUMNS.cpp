
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_COLUMNS: Generate columns using the outputs of LN2_LAYERS.\n"
    "\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "!!! WORK IN PROGRESS !!!\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!"
    "\n"
    "Usage:\n"
    "    LN2_COLUMNS\n"
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
    "    -midgm        : Middle gray matter file (from LN2_LAYERS output).\n"
    "    -nr_columns   : Number of columns.\n"
    "    -debug        : (Optional) Save extra intermediate outputs.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2;
    uint16_t ac, nr_columns = 5;
    bool mode_debug = true;

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
        } else if (!strcmp(argv[ac], "-nr_columns")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -nr_columns\n");
            } else {
                nr_columns = atof(argv[ac]);
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

    log_welcome("LN2_COLUMNS");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    cout << "\n  Nr. columns: " << nr_columns << endl;

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

    // Prepare required nifti images
    nifti_image* nii_columns  = copy_nifti_as_int16(nii_rim);
    int16_t* nii_columns_data = static_cast<int16_t*>(nii_columns->data);
    // Setting zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_columns_data + i) = 0;
    }

    nifti_image* flood_step = copy_nifti_as_int16(nii_columns);
    int16_t* flood_step_data = static_cast<int16_t*>(flood_step->data);
    nifti_image* flood_dist = copy_nifti_as_float32(nii_columns);
    float* flood_dist_data = static_cast<float*>(flood_dist->data);

    nifti_image* anchor_id = copy_nifti_as_int32(nii_columns);
    int32_t* anchor_id_data = static_cast<int32_t*>(anchor_id->data);

    // ========================================================================
    // Grow from an initial point
    // ========================================================================
    // Find the first voxel
    uint32_t start_voxel;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 1) {
            start_voxel = i;
        }
    }
    *(anchor_id_data + start_voxel) = 1;
    *(nii_midgm_data + start_voxel) = 2;

    // Loop until desired number of columns reached
    for (uint32_t n = 1; n != nr_columns; ++n) {
        cout << "\n  Start flooding..." << n << endl;

        // Initialize new point
        int32_t new_voxel_id;

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

        uint16_t grow_step = 1;
        uint32_t voxel_counter = nr_voxels;
        uint32_t ix, iy, iz, j, k;
        float d;
        while (voxel_counter != 0) {
            voxel_counter = 0;
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                if (*(flood_step_data + i) == grow_step) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    voxel_counter += 1;

                    // ------------------------------------------------------------
                    // 1-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dX;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dY;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 2-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 3-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                        if (*(nii_midgm_data + j) == 1) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                //*(nii_columns_data + j) = *(anchor_id_data + i);
                                new_voxel_id = j;
                            }
                        }
                    }
                }
            }
            grow_step += 1;
        }
        *(nii_midgm_data + new_voxel_id) = 2;
    }

    if (mode_debug) {
        save_output_nifti(fout, "flood_step", flood_step, false);
        save_output_nifti(fout, "flood_dist", flood_dist, false);
        save_output_nifti(fout, "anchors", nii_midgm, false);
    }

    // ========================================================================
    // Voronoi cell from the points
    // ========================================================================
    cout << "\n  Start Voronoi..." << endl;

    // Prepare centroid ids
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_midgm_data + i) == 2) {
            *(nii_columns_data + i) = i;
            *(nii_rim_data + i) = 0;
        }
    }

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

    uint16_t grow_step = 1;
    uint32_t voxel_counter = nr_voxels;
    uint32_t ix, iy, iz, j, k;
    float d;
    voxel_counter = nr_voxels;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(flood_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);
                        }
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);


                        }
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xy;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_yz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);

                    if (*(nii_rim_data + j) == 3) {
                        d = *(flood_dist_data + i) + dia_xyz;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_columns_data + j) = *(nii_columns_data + i);

                        }
                    }
                }
            }
        }
        grow_step += 1;
    }

    // ========================================================================
    save_output_nifti(fout, "columns", nii_columns, true);
    if (mode_debug) {
        save_output_nifti(fout, "voronoi_flood_step", flood_step, false);
        save_output_nifti(fout, "voronoi_flood_dist", flood_dist, false);
    }

    cout << "\n  Finished." << endl;
    return 0;
}
