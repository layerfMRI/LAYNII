
#include <vector>
#include <algorithm>
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
        "LN2_LLOYD: Lloyd's algorithm.\n"
        "\n"
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        "!! BEWARE! WORK IN PROGRESS... USE WITH CAUTION !!\n"
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        "\n"
        "Usage:\n"
        "    LN2_LAYERS -rim rim.nii -init rim_initial_spots.nii\n"
        "\n"
        "Options:\n"
        "    -help         : Show this help. \n"
        "    -rim          : Specify input dataset.\n"
        "    -init         : Labels for initializing centroids.\n"
        "    -nr_iter      : Number of iterations for centroid updates.\n"
        "\n");
    return 0;
}

// ============================================================================
// Local functions

void run_updates(int i, int j, float d, float grow_step,
                 float * dist_data, float * step_data, int32_t * cell_data) {
    // Procedure for updating nifti data in neighbour grow section
    // d = pow(d, 2);
    if (d < *(dist_data + j) || *(step_data + j) == 0) {
        *(dist_data + j) = d;
        *(step_data + j) = grow_step + 1;
        *(cell_data + j) = *(cell_data + i);
    }
}

void run_updates(int i, int j, float d, float grow_step,
                 float * dist_data, float * step_data, int32_t * cell_data);

// ============================================================================

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fin2 = NULL;
    uint16_t ac;
    int nr_iter = 10;

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
        } else if (!strcmp(argv[ac], "-init")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -init\n");
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-nr_iter")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -nr_iter\n");
            } else {
                nr_iter = atof(argv[ac]);
            }
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
        fprintf(stderr, "** missing option '-init'\n");
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

    log_welcome("LN2_LLOYD");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;

    const int end_x = size_x - 1;
    const int end_y = size_y - 1;
    const int end_z = size_z - 1;

    const int nr_voxels = size_z * size_y * size_x;

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
    nifti_image* nii_rim = copy_nifti_as_int32(nii1);
    int32_t* nii_rim_data = static_cast<int32_t*>(nii_rim->data);
    nifti_image* nii_init = copy_nifti_as_int32(nii2);
    int32_t* nii_init_data = static_cast<int32_t*>(nii_init->data);

    // Prepare required nifti images
    nifti_image* nii_centroid = copy_nifti_as_int32(nii_rim);
    int32_t* nii_centroid_data = static_cast<int32_t*>(nii_centroid->data);
    for (int i = 0; i < nr_voxels; ++i) {
        *(nii_centroid_data + i) = 0;
    }

    nifti_image* nii_step = copy_nifti_as_float32(nii_centroid);
    float* nii_step_data = static_cast<float*>(nii_step->data);
    nifti_image* nii_dist = copy_nifti_as_float32(nii_centroid);
    float* nii_dist_data = static_cast<float*>(nii_dist->data);
    nifti_image* nii_cell = copy_nifti_as_int32(nii_centroid);
    int32_t* nii_cell_data = static_cast<int32_t*>(nii_cell->data);
    nifti_image* nii_border = copy_nifti_as_int32(nii_centroid);
    int32_t* nii_border_data = static_cast<int32_t*>(nii_border->data);

    // 4D output
    nifti_image* nii_new = nifti_copy_nim_info(nii_centroid);
    nii_new->datatype = NIFTI_TYPE_FLOAT32;
    nii_new->dim[0] = 4;
    nii_new->dim[1] = size_x;
    nii_new->dim[2] = size_y;
    nii_new->dim[3] = size_z;
    nii_new->dim[4] = nr_iter;
    nifti_update_dims_from_array(nii_new);
    nii_new->nvox = nr_iter * nr_voxels;
    nii_new->nbyper = sizeof(float);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    float* nii_new_data = static_cast<float*>(nii_new->data);

    nifti_image* nii_new_cells = copy_nifti_as_int32(nii_new);
    int32_t* nii_new_cell_data = static_cast<int32_t*>(nii_new_cells->data);

    // ------------------------------------------------------------------------
    int nr_valid_voxels = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) == 3) {
            nr_valid_voxels += 1;
        }
    }

    // ========================================================================
    // Find centroid ids
    // ========================================================================
    std::vector<int> v(nr_voxels);
    for (int i = 0; i != nr_voxels; ++i) {
        v[i] = *(nii_init_data + i);
    }
    // Remove duplicate elements
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
    int nr_cells = v.size();

    // ========================================================================
    // Find centroids of initial cells
    // ========================================================================
    std::vector<float> c_x(nr_cells, 0);
    std::vector<float> c_y(nr_cells, 0);
    std::vector<float> c_z(nr_cells, 0);
    std::vector<float> counts(nr_cells, 0);

    // Sum x, y, z coordinates of same-column middle GM voxels
    int j, x, y, z;
    for (int i = 0; i != nr_voxels; ++i) {
        j = *(nii_init_data + i);  // used to determine storage voxel
        if (j > 0) {
            tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
            c_x[j] += x;
            c_y[j] += y;
            c_z[j] += z;
            counts[j] += 1;
        }
    }
    // Divide summed coordinates to find centroids
    for (int i = 0; i != nr_cells; ++i) {
        if (counts[i] != 0) {
            c_x[i] = floor(c_x[i] / counts[i]);
            c_y[i] = floor(c_y[i] / counts[i]);
            c_z[i] = floor(c_z[i] / counts[i]);
            j = sub2ind_3D(c_x[i], c_y[i], c_z[i], size_x, size_y);
            *(nii_centroid_data + j) = i;
        }
    }
    save_output_nifti(fin1, "centroids", nii_centroid, true);

    // ====================================================================
    // Grow from centroids
    // ====================================================================
    cout << "\n  Start growing from centroids..." << endl;

    // Initialize grow volume
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_centroid_data + i) != 0) {
            *(nii_step_data + i) = 1.;
            *(nii_dist_data + i) = 0.;
            *(nii_cell_data + i) = *(nii_centroid_data + i);
        } else {
            *(nii_step_data + i) = 0.;
            *(nii_dist_data + i) = 0.;
        }
    }

    for (int n = 0; n < nr_iter; ++n) {
        cout << "\r    Iteration: " << n+1 << "/" << nr_iter << flush;

        int grow_step = 1, voxel_counter = nr_voxels;
        int ix, iy, iz;
        float d;

        while (voxel_counter != 0) {
            voxel_counter = 0;
            for (int i = 0; i != nr_voxels; ++i) {
                if (*(nii_step_data + i) == grow_step
                    && *(nii_rim_data + i) == 3) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    voxel_counter += 1;

                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dX;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dX;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dY;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dY;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dZ;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dZ;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }

                    // --------------------------------------------------------
                    // 2-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xy;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xy;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xy;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xy;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_yz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_yz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_yz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_yz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }

                    // --------------------------------------------------------
                    // 3-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_rim_data + j) == 3) {
                            d = *(nii_dist_data + i) + dia_xyz;
                            run_updates(i, j, d, grow_step,
                                        nii_dist_data,
                                        nii_step_data,
                                        nii_cell_data);
                        }
                    }
                }
            }
            grow_step += 1;
        }

        // Update 4D image
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + n * nr_voxels + i) = *(nii_dist_data + i);
            *(nii_new_cell_data + n * nr_voxels + i) = *(nii_cell_data + i);
        }

        // // ----------------------------------------------------------------
        // // Find borders
        // // ----------------------------------------------------------------
        // for (int i = 0; i < nr_voxels; ++i) {
        //     if (*(nii_rim_data + i) == 3) {
        //         tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
        //         // --------------------------------------------------------
        //         // 1-jump neighbours
        //         // --------------------------------------------------------
        //         int count_n = 0;
        //         if (ix > 0) {
        //             j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
        //             if (*(nii_cell_data + j) == *(nii_cell_data + i)) {
        //                 count_n += 1;
        //             }
        //         }
        //         if (ix < end_x) {
        //             j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
        //             if (*(nii_cell_data + j) == *(nii_cell_data + i)) {
        //                 count_n += 1;
        //             }
        //         }
        //         if (iy > 0) {
        //             j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
        //             if (*(nii_cell_data + j) == *(nii_cell_data + i)) {
        //                 count_n += 1;
        //             }
        //         }
        //         if (iy < end_y) {
        //             j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
        //             if (*(nii_cell_data + j) == *(nii_cell_data + i)) {
        //                 count_n += 1;
        //             }
        //         }
        //         if (iz > 0) {
        //             j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
        //             if (*(nii_cell_data + j) == *(nii_cell_data + i)) {
        //                 count_n += 1;
        //             }
        //         }
        //         if (iz < end_z) {
        //             j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
        //             if (*(nii_cell_data + j) == *(nii_cell_data + i)) {
        //                 count_n += 1;
        //             }
        //         }
        //
        //         if (count_n < 6) {
        //             *(nii_border_data + i) = *(nii_cell_data + i);
        //         } else {
        //             *(nii_border_data + i) = 0;
        //         }
        //     }
        // }

        if (n == 0) {  // For debugging
            save_output_nifti(fin1, "lloyd_step", nii_step, false);
            save_output_nifti(fin1, "lloyd_dist", nii_dist, false);
            save_output_nifti(fin1, "lloyd_cell", nii_cell, false);
            save_output_nifti(fin1, "lloyd_border", nii_border, false);
        }

        // ====================================================================
        // Find centroids of new cells
        // ====================================================================
        // NOTE(Faruk): This part can be move at the top of the loop
        if (n != nr_iter - 1) {
            for (int i = 0; i != nr_cells; ++i) {
                c_x[i] = 0;
                c_y[i] = 0;
                c_z[i] = 0;
                counts[i] = 0;
            }

            // Sum x, y, z coordinates
            for (int i = 0; i != nr_voxels; ++i) {
                // j = *(nii_border_data + i);
                j = *(nii_cell_data + i);
                if (j > 0) {
                    tie(x, y, z) = ind2sub_3D(i, size_x, size_y);
                    c_x[j] += x;
                    c_y[j] += y;
                    c_z[j] += z;
                    counts[j] += 1;
                }
            }
            // Reset centroids
            for (int i = 0; i != nr_voxels; ++i) {
                *(nii_centroid_data + i) = 0;
            }
            // Divide summed coordinates to find centroids
            for (int i = 0; i != nr_cells; ++i) {
                if (counts[i] != 0) {
                    c_x[i] = floor(c_x[i] / counts[i]);
                    c_y[i] = floor(c_y[i] / counts[i]);
                    c_z[i] = floor(c_z[i] / counts[i]);
                    j = sub2ind_3D(c_x[i], c_y[i], c_z[i], size_x, size_y);
                    *(nii_centroid_data + j) = i;
                }
            }

            // // TODO(Faruk): To penalize unequal volume cells
            // for (int i = 0; i != nr_cells; ++i) {
            //     counts[i] /= nr_valid_voxels / (nr_cells-1);
            //     cout << "  " << counts[i] << " | " << flush;
            // }
            // cout << endl;

            // Reset other images
            for (int i = 0; i != nr_voxels; ++i) {
                if (*(nii_centroid_data + i) != 0) {
                    *(nii_step_data + i) = 1.;
                    *(nii_dist_data + i) = 0.;
                    *(nii_cell_data + i) = *(nii_centroid_data + i);
                } else {
                    *(nii_step_data + i) = 0.;
                    *(nii_dist_data + i) = 0.;
                }
            }
        }
    }
    cout << "\n" << endl;

    save_output_nifti(fin1, "lloyd_dist", nii_new, true);
    save_output_nifti(fin1, "lloyd", nii_new_cells, true);

    // For debugging
    save_output_nifti(fin1, "lloyd_step_final", nii_step, false);
    save_output_nifti(fin1, "lloyd_dist_final", nii_dist, false);
    save_output_nifti(fin1, "lloyd_cell_final", nii_cell, false);
    save_output_nifti(fin1, "lloyd_border_final", nii_border, false);


    cout << "\n  Finished." << endl;
    return 0;
}
