
#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_CONNECTED_CLUSTERS: Find connected clusters in a binary image.\n"
    "\n"
    "Usage:\n"
    "    LN2_CONNECTED_CLUSTERS -input input.nii\n"
    "    ../LN2_CONNECTED_CLUSTERS -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -input        : Binary nifti image (only consists of 0s and 1s).\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    log_welcome("LN2_CONNECTED_CLUSTERS");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_int32(nii1);
    int32_t* nii_input_data = static_cast<int32_t*>(nii_input->data);

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_input_data + i) != 0){
            *(nii_input_data + i) = 1;
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_input_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Find connected clusters
    // ========================================================================
    cout << "  Start finding connected clusters (3-jump neighbourhood)..." << endl;

    // Loop until all clusters have one initial voxel
    uint32_t voxel_counter = 0, prev_voxel_counter = 0;
    int32_t init_voxel_id = 1;
    bool terminate_switch1 = true;
    while (terminate_switch1) {
        uint32_t ix, iy, iz, i, j;
        cout << "  " << voxel_counter << "/" << nr_voi << flush;

        if (voxel_counter == nr_voi) {
            // Indicates all clusters are reached. Terminate condition.
            terminate_switch1 = false;
            cout << "    Nr. of connected clusters within midgm input: "
                << init_voxel_id - 1 << endl;
        } else if (voxel_counter == prev_voxel_counter) {
            // Find the initial voxel for each disconnected cluster
            uint32_t start_voxel;
            for (uint32_t i = 0; i != nr_voxels; ++i) {
                if (*(nii_input_data + i) == 1) {
                    start_voxel = i;
                }
            }
            voxel_counter += 1;
            init_voxel_id += 1;
            *(nii_input_data + start_voxel) = init_voxel_id;
        }

        while (prev_voxel_counter != voxel_counter) {
            prev_voxel_counter = voxel_counter;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                // Map subset to full set
                i = *(voi_id + ii);
                if (*(nii_input_data + i) == init_voxel_id) {
                    tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                    // --------------------------------------------------------
                    // 1-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0) {
                        j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x) {
                        j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iy > 0) {
                        j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iy < end_y) {
                        j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iz > 0) {
                        j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iz < end_z) {
                        j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    // --------------------------------------------------------
                    // 2-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy < end_y) {
                        j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy < end_y) {
                        j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iz > 0) {
                        j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }

                    // --------------------------------------------------------
                    // 3-jump neighbours
                    // --------------------------------------------------------
                    if (ix > 0 && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy > 0 && iz > 0) {
                        j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix > 0 && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy > 0 && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy < end_y && iz > 0) {
                        j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                    if (ix < end_x && iy < end_y && iz < end_z) {
                        j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                        if (*(nii_input_data + j) == 1) {
                            *(nii_input_data + j) = init_voxel_id;
                        }
                    }
                }
            }

            // Count cluster assigned voxels
            voxel_counter = 0;
            for (uint32_t ii = 0; ii != nr_voi; ++ii) {
                i = *(voi_id + ii);  // Map subset to full set
                if (*(nii_input_data + i) > 1) {
                    voxel_counter += 1;
                }
            }
        }
    }
    cout << endl;
    cout << "  Nr. connected clusters = " << init_voxel_id - 1 << endl;

    // Make uniqute cluster ids start from 1 instead of 2
    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        uint32_t i = *(voi_id + ii);  // Map subset to full set
        if (*(nii_input_data + i) != 0) {
            *(nii_input_data + i) -= 1;
        }
    }

    // Add number of clusters into the output tag
    std::ostringstream tag;
    tag << init_voxel_id - 1;
    save_output_nifti(fout, "connected_clusters" + tag.str(), nii_input, true);

    cout << "\n  Finished." << endl;
    return 0;
}
