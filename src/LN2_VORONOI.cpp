
#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_VORONOI: Voronoi propagation from a set of voxels.\n"
    "\n"
    "Usage:\n"
    "    LN2_VORONOI -domain domain.nii -init initial_voxels.nii\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -domain       : Set of voxels in which points will be used in Voronoi.\n"
    "                    propagation.\n"
    "    -init         : Initial voxels.\n"
    "    -max_dist     : (Optional) Maximum distance from the initial voxels\n"
    "                    where Voronoi cells will be propagated.\n"
    "    -iter_smooth  : (Optional) Number of smoothing iterations. Default\n"
    "                    is 0 (no smoothing).\n"
    "    -debug        : (Optional) Save extra intermediate outputs.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL;
    int ac;
    bool mode_debug = false, mode_initialize_with_centroids = false;
    float max_dist = std::numeric_limits<float>::max();
    int iter_smooth = 0;

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
            fin2 = argv[ac];
            mode_initialize_with_centroids = true;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-max_dist")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -max_dist\n");
            } else {
                max_dist = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-iter_smooth")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -iter_smooth\n");
            } else {
                iter_smooth = atof(argv[ac]);
            }
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
        if (!fin2) {
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
    nii2 = nifti_image_read(fin2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }

    log_welcome("LN2_VORONOI");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    if (max_dist < std::numeric_limits<float>::max()) {
        cout << "    Maximum distance is: " << max_dist << endl;
    } else
        cout << "    No maximum distance is selected." << endl;

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
    nifti_image* nii_init  = copy_nifti_as_int32(nii2);
    int32_t* nii_init_data = static_cast<int32_t*>(nii_init->data);

    nifti_image* flood_step = copy_nifti_as_int32(nii_init);
    int32_t* flood_step_data = static_cast<int32_t*>(flood_step->data);
    nifti_image* flood_dist = copy_nifti_as_float32(nii_init);
    float* flood_dist_data = static_cast<float*>(flood_dist->data);

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
    // Grow Voronoi cells from points towards the rest of the domain
    // ========================================================================
    cout << "\n  Start growing Voronoi cells..." << endl;

    // Reset domain
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        uint32_t i = *(voi_id + ii);
        *(nii_domain_data + i) = 1;
    }

    // Initialize grow volume
    uint32_t i;
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);
        if (*(nii_init_data + i) != 0) {
            *(flood_step_data + i) = 1;
            *(flood_dist_data + i) = 0;
        } else {
            *(flood_step_data + i) = 0;
            *(flood_dist_data + i) = 0;
        }
    }

    int32_t grow_step = 1;
    uint32_t ix, iy, iz, j;
    float d;
    int voxel_counter = nr_voxels;
    while (voxel_counter != 0) {
        voxel_counter = 0;
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            i = *(voi_id + ii);
            if (*(flood_step_data + i) == grow_step && *(flood_dist_data + i) < max_dist) {
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
                            *(nii_init_data + j) = *(nii_init_data + i);
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
                            *(nii_init_data + j) = *(nii_init_data + i);
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
                            *(nii_init_data + j) = *(nii_init_data + i);
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
                            *(nii_init_data + j) = *(nii_init_data + i);
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
                            *(nii_init_data + j) = *(nii_init_data + i);
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
                            *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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
                                *(nii_init_data + j) = *(nii_init_data + i);
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

    // Smooth
    if (iter_smooth > 0) {
        nii_init = iterative_smoothing(nii_init, iter_smooth, nii_domain, 1);
    }

    // Threshold
    if (iter_smooth > 0) {
        cout << "\n  Start mildly smoothing distances before thresholding..." << endl;
        flood_dist = iterative_smoothing(flood_dist, 3, nii_domain, 1);
        float* flood_dist_data = static_cast<float*>(flood_dist->data);

        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            i = *(voi_id + ii);
            if (*(flood_dist_data + i) > max_dist) {
                *(nii_init_data + i) = 0;
            }
        }
    }

    // Add number of points into the output tag
    save_output_nifti(fout, "voronoi", nii_init, true);

    cout << "\n  Finished." << endl;
    return 0;
}
