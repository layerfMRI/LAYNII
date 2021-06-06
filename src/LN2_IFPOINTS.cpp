
#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_IFPOINTS: Iterative farthest points finder.\n"
    "\n"
    "Usage:\n"
    "    LN2_IFPOINTS -domain domain.nii -nr_points 10\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -domain       : Set of voxels in which points will be generated.\n"
    "                    All non-zero voxels will be considered.\n"
    "    -nr_points    : Number of points that will be generated\n"
    "    -init         : (Optional) New points will be added based on these\n"
    "                    initial points.\n"
    "    -debug        : (Optional) Save extra intermediate outputs.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n"
    "Notes:\n"
    "    - This program is a stripped-down version of LN2_COLUMNS that work\n"
    "      with fewer inputs and constraints.\n"
    "    - You can find further explanation of LN2_COLUMNS at:\n"
    "      <https://thingsonthings.org/ln2_columns>\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    bool use_outpath = false;
    nifti_image *nii1 = NULL, *nii3 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin3=NULL;
    int ac;
    int32_t nr_points = 3;
    bool mode_debug = false, mode_initialize_with_centroids = false;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-domain")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -domain\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-init")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -init\n");
                return 1;
            }
            fin3 = argv[ac];
            mode_initialize_with_centroids = true;
        } else if (!strcmp(argv[ac], "-nr_points")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -nr_points\n");
            } else {
                nr_points = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
            use_outpath = true;
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-domain'\n");
        return 1;
    }
    if (mode_initialize_with_centroids) {
        if (!fin3) {
            fprintf(stderr, "** missing option '-init'\n");
            return 1;
        }
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }
    if (mode_initialize_with_centroids) {
        nii3 = nifti_image_read(fin3, 1);
        if (!nii3) {
            fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin3);
            return 2;
        }
    }

    log_welcome("LN2_IFPOINTS");
    log_nifti_descriptives(nii1);

    if (mode_initialize_with_centroids) {
        log_nifti_descriptives(nii3);
    }

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
    nifti_image* nii_domain = copy_nifti_as_int32(nii1);
    int32_t* nii_domain_data = static_cast<int32_t*>(nii_domain->data);
    // Binarize
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_domain_data + i) != 0) {
            *(nii_domain_data + i) = 1;
        } else {
            *(nii_domain_data + i) = 0;
        }
    }

    // Prepare required nifti images
    nifti_image* nii_points  = copy_nifti_as_int32(nii_domain);
    int32_t* nii_points_data = static_cast<int32_t*>(nii_points->data);
    // Setting zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_points_data + i) = 0;
    }

    nifti_image* flood_step = copy_nifti_as_int32(nii_points);
    int32_t* flood_step_data = static_cast<int32_t*>(flood_step->data);
    nifti_image* flood_dist = copy_nifti_as_float32(nii_points);
    float* flood_dist_data = static_cast<float*>(flood_dist->data);

    // ------------------------------------------------------------------------
    // Find initial number of points if the optional init input is given
    int32_t max_point_id = 0;
    if (mode_initialize_with_centroids) {
        nifti_image* nii_init_points = copy_nifti_as_int32(nii3);
        int32_t* nii_init_points_data = static_cast<int32_t*>(nii_init_points->data);

        // Find maximum column id
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_init_points_data + i) > max_point_id) {
                max_point_id = *(nii_init_points_data + i);
            }
        }

        // Remove centroids if the desired number of columns is less than
        // initially given centroids.
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_init_points_data + i) > nr_points) {
                *(nii_points_data + i) = 0;
            } else if (*(nii_init_points_data + i) < 0) {  // for signed ids error
                *(nii_points_data + i) = 0;
            } else {
                *(nii_points_data + i) = *(nii_init_points_data + i);
            }
        }
        if (mode_debug) {
            save_output_nifti(fout, "initial_points", nii_points, false);
        }
    }
    cout << "  Initial number of points: " << max_point_id << endl;
    cout << "  Desired number of points: " << nr_points << endl;

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_domain_data + i) != 0) {
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_domain_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Find connected clusters to initialize one voxel in each
    // ========================================================================
    cout << "  Start finding connected clusters..." << endl;

    // Loop until all clusters have one initial voxel
    uint32_t voxel_counter = 0, prev_voxel_counter = 0;
    int32_t init_voxel_id = 1;
    bool terminate_switch1 = true;
    while (terminate_switch1) {
        uint32_t ix, iy, iz, i, j;

        if (voxel_counter == nr_voi) {
            // Indicates all clusters are reached. Terminate condition.
            terminate_switch1 = false;
            cout << "    Nr. of connected clusters within domain: "
                << init_voxel_id - 1 << endl;
        } else if (voxel_counter == prev_voxel_counter) {
            // Find the initial voxel for each disconnected cluster
            uint32_t start_voxel;
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                if (*(nii_domain_data + i) != 0) {
                    start_voxel = i;
                }
            }
            voxel_counter += 1;
            init_voxel_id += 1;
            *(nii_domain_data + start_voxel) = init_voxel_id;
        }

        while (prev_voxel_counter != voxel_counter) {
            prev_voxel_counter = voxel_counter;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                i = *(voi_id + ii);  // Map subset to full set
                if (*(nii_domain_data + i) == init_voxel_id) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    // --------------------------------------------------------
                    // 2-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }

                    // --------------------------------------------------------
                    // 3-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) == 1) {
                            *(nii_domain_data + j) = init_voxel_id;
                        }
                    }
                }
            }

            // Count cluster assigned voxels
            voxel_counter = 0;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                i = *(voi_id + ii);  // Map subset to full set
                if (*(nii_domain_data + i) > 1) {
                    voxel_counter += 1;
                }
            }
        }
    }
    if (mode_debug) {
        save_output_nifti(fout, "connected_clusters", nii_domain, false);
    }

    // ========================================================================
    // Find points through farthest flood distance
    // ========================================================================
    cout << "  Start generating points..." << endl;
    // Find the initial voxel
    uint32_t start_voxel;
    for (int32_t n = 2; n <= init_voxel_id; ++n) {
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            uint32_t i = *(voi_id + ii);  // Map subset to full set
            if (*(nii_domain_data + i) == n) {
                start_voxel = i;
                *(nii_domain_data + i) = 1;  // Reset midgm
            }
        }
        *(nii_domain_data + start_voxel) = 2;  // Reduce to single initial voxel
    }

    // Initialize new voxel
    uint32_t new_voxel_id;
    float flood_dist_thr = std::numeric_limits<float>::infinity();

    // Loop until desired number of points reached
    for (int32_t n = max_point_id; n < nr_points; ++n) {
        cout << "\r    Point [" << n+1 << "/" << nr_points << "]" << flush;

        int32_t grow_step = 1;
        voxel_counter = 1;
        uint32_t ix, iy, iz, i, j;
        float d;

        // Initialize grow volume
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            if (*(nii_domain_data + i) == 2) {
                *(flood_step_data + i) = 1.;
                *(flood_dist_data + i) = 0.;
            } else if (*(flood_dist_data + i) >= flood_dist_thr
                       && *(flood_dist_data + i) > 0) {
                *(flood_step_data + i) = 0.;
                *(flood_dist_data + i) = 0.;
                *(nii_domain_data + i) = 1;
            } else if (*(flood_dist_data + i) < flood_dist_thr
                       && *(flood_dist_data + i) > 0) {
                *(nii_domain_data + i) = 0;  // no need to recompute
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

                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
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
                        if (*(nii_domain_data + j) != 0) {
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
                        if (*(nii_domain_data + j) != 0) {
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
                        if (*(nii_domain_data + j) != 0) {
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
                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dZ;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }
                    // --------------------------------------------------------
                    // 2-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                new_voxel_id = j;
                            }
                        }
                    }

                    // --------------------------------------------------------
                    // 3-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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

                        if (*(nii_domain_data + j) != 0) {
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
        flood_dist_thr = *(flood_dist_data + new_voxel_id) / 2.;
        *(nii_domain_data + new_voxel_id) = 2;
        *(nii_points_data + new_voxel_id) = n + 1;
    }
    cout << endl;

    // Add number of points into the output tag
    std::ostringstream tag;
    tag << nr_points;
    save_output_nifti(fout, "points"+tag.str(), nii_points, true, use_outpath);

    // ========================================================================
    // Grow Voronoi cells from points towards the rest of the domain
    // ========================================================================
    cout << "\n  Start growing Voronoi cells..." << endl;

    // Reset domain
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        uint32_t i = *(voi_id + ii);
        *(nii_domain_data + i) = 1;
    }

    // Initialize grow volume
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_points_data + i) != 0) {
            *(flood_step_data + i) = 1;
            *(flood_dist_data + i) = 0;
        } else {
            *(flood_step_data + i) = 0;
            *(flood_dist_data + i) = 0;
        }
    }

    int32_t grow_step = 1;
    uint32_t ix, iy, iz, i, j;
    float d;
    voxel_counter = nr_voxels;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            i = *(voi_id + ii);
            if (*(flood_step_data + i) == grow_step) {
                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
                voxel_counter += 1;
                bool jump_lock = false;
                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_domain_data + j) != 0) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_points_data + j) = *(nii_points_data + i);
                        }
                    } else {
                        jump_lock = true;
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_domain_data + j) != 0) {
                        d = *(flood_dist_data + i) + dX;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_points_data + j) = *(nii_points_data + i);
                        }
                    } else {
                        jump_lock = true;
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_domain_data + j) != 0) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_points_data + j) = *(nii_points_data + i);
                        }
                    } else {
                        jump_lock = true;
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_domain_data + j) != 0) {
                        d = *(flood_dist_data + i) + dY;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_points_data + j) = *(nii_points_data + i);
                        }
                    } else {
                        jump_lock = true;
                    }
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    if (*(nii_domain_data + j) != 0) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_points_data + j) = *(nii_points_data + i);
                        }
                    } else {
                        jump_lock = true;
                    }
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    if (*(nii_domain_data + j) != 0) {
                        d = *(flood_dist_data + i) + dZ;
                        if (d < *(flood_dist_data + j)
                            || *(flood_dist_data + j) == 0) {
                            *(flood_dist_data + j) = d;
                            *(flood_step_data + j) = grow_step + 1;
                            *(nii_points_data + j) = *(nii_points_data + i);
                        }
                    } else {
                        jump_lock = true;
                    }
                }

                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
                if (jump_lock == false) {

                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xy;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_yz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }

                    // ------------------------------------------------------------
                    // 3-jump neighbours
                    // ------------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_domain_data + j) != 0) {
                            d = *(flood_dist_data + i) + dia_xyz;
                            if (d < *(flood_dist_data + j)
                                || *(flood_dist_data + j) == 0) {
                                *(flood_dist_data + j) = d;
                                *(flood_step_data + j) = grow_step + 1;
                                *(nii_points_data + j) = *(nii_points_data + i);
                            }
                        }
                    }
                }
            }
        }
        grow_step += 1;
    }

    if (mode_debug) {
        save_output_nifti(fout, "flood_step", flood_step, false);
        save_output_nifti(fout, "flood_dist", flood_dist, false);
    }
    // Add number of points into the output tag
    save_output_nifti(fout, "cells"+tag.str(), nii_points, true, use_outpath);

    cout << "\n  Finished." << endl;
    return 0;
}
