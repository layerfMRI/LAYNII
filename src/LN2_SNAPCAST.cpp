#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_SNAPCAST: Create snapshots of isometric projections over the cardinal axes.\n"
    "              It is intended to be used as a fast way to visualize subdural and\n"
    "              leptomeningeal brain surfaces.\n"
    "\n"
    "Usage:\n"
    "    LN2_SNAPCAST -input T2star.nii.gz -mask brain_mask.nii.gz\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : A 3D or 4D nifti file that contains values to project.\n"
    "              Only non zero voxels will be used.\n"
    "    -type   : Projection type. Options are: 'mean', 'min', 'max'. Default is 'mean'.\n"
    "    -steps  : A positive integer. Determines the number of voxels for\n"
    "              ray penetration. Higher values give smoother results. Default is '5'.\n"
    // "    -mask     : (Optional) 3D binary nifti file that will be used to mask the input.\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    int64_t steps = 5, step_count_all = 0, step_count = 0;
    std::string projection_type = "mean";

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
        } else if (!strcmp(argv[ac], "-steps")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -steps\n");
            } else {
                steps = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-type")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -type\n");
                return 1;
            }
            projection_type = argv[ac];
            for (auto &c : projection_type) {  // Convert to lowercase
                c = std::tolower(c);
            }
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

    log_welcome("LN2_SNAPCAST");
    // Check projection type input
    if (projection_type == "mean" || projection_type == "min" || projection_type == "max" ) {
        std::cout << "  Type: '" << projection_type << "' intensity projection.\n" << std::endl;
    } else {
        fprintf(stderr, "** invalid projection type: '%s'\n", projection_type.c_str());
        return 1;          
    }
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const int64_t size_x = nii1->nx;
    const int64_t size_y = nii1->ny;
    const int64_t size_z = nii1->nz;
    const int64_t size_time = nii1->nt;

    const int64_t end_x = size_x - 1;
    const int64_t end_y = size_y - 1;
    const int64_t end_z = size_z - 1;

    const int64_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // ========================================================================
    // Compute the new dimensions of output
    // ========================================================================
    const int64_t size_max = std::max(size_z, std::max(size_x, size_y));
    const int64_t size_x_out = size_max;
    const int64_t size_y_out = size_max;
    const int64_t size_z_out = 6;  // For 6 cube faces (cardinal data axes)
    const int64_t size_time_out = size_time;

    const int64_t nr_voxels_out = size_z_out * size_y_out * size_x_out;

    const int64_t offset_x = (size_max - size_x) / 2;
    const int64_t offset_y = (size_max - size_y) / 2;
    const int64_t offset_z = (size_max - size_z) / 2;

    // ========================================================================
    // Allocating new nifti for output 4D nifti
    // ========================================================================
    // Prepare the output nifti
    nifti_image* nii_output = nifti_copy_nim_info(nii_input);
    nii_output->datatype = NIFTI_TYPE_FLOAT32;
    nii_output->dim[0] = 4;  // For proper 4D nifti
    nii_output->dim[1] = size_x_out;
    nii_output->dim[2] = size_y_out;
    nii_output->dim[3] = size_z_out;
    nii_output->dim[4] = size_time_out;
    nifti_update_dims_from_array(nii_output);
    nii_output->nvox = nr_voxels_out * size_time_out;
    nii_output->nbyper = sizeof(float);
    nii_output->data = calloc(nii_output->nvox, nii_output->nbyper);
    nii_output->scl_slope = 1;
    float* nii_output_data = static_cast<float*>(nii_output->data);

    log_nifti_descriptives(nii_output);

    // Prepare for mean intensity projection
    for (int64_t i = 0; i != nr_voxels_out * size_time_out; ++i) {
        *(nii_output_data + i) = 0;
    }

    // ================================================================================================================
    // Compute rays with mean intensity projection
    // ================================================================================================================
    if ( projection_type == "mean" ) {
        for (int64_t t = 0; t != size_time; ++t) {  // Loop across time points
            cout << "\r    Volume: " << t+1 << "/" << size_time << flush;

            for (int64_t z = 0; z != size_z; ++z) { // Loop across rays on x plane
                for (int64_t y = 0; y != size_y; ++y) {

                    // Loop across a ray
                    for (int64_t x = 0; x != size_x; ++x) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            int64_t k = sub2ind_4D_64(y + offset_y, z + offset_z, 0, t, size_max, size_max, 6);
                            *(nii_output_data + k) += *(nii_input_data + i) * 1. / steps;
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t x = size_x-1; x >= 0; --x) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            int64_t k = sub2ind_4D_64(y + offset_y, z + offset_z, 1, t, size_max, size_max, 6);
                            *(nii_output_data + k) += *(nii_input_data + i) * 1. / steps;
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }

            // --------------------------------------------------------------------
            // Loop across rays on y plane
            for (int64_t z = 0; z != size_z; ++z) {
                for (int64_t x = 0; x != size_x; ++x) {

                    // Loop across a ray
                    for (int64_t y = 0; y != size_y; ++y) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            if ( step_count == steps || step_count_all == size_y ) {
                                break;
                            }
                            int64_t k = sub2ind_4D_64(x + offset_x, z + offset_z, 2, t, size_max, size_max, 6);
                            *(nii_output_data + k) += *(nii_input_data + i) * 1. / steps;
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t y = size_y-1; y >= 0; --y) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            if ( step_count == steps || step_count_all == size_y ) {
                                break;
                            }
                            int64_t k = sub2ind_4D_64(x + offset_x, z + offset_z, 3, t, size_max, size_max, 6);
                            *(nii_output_data + k) += *(nii_input_data + i) * 1. / steps;
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }

            // --------------------------------------------------------------------
            // Loop across rays on z plane
            for (int64_t y = 0; y != size_y; ++y) {
                for (int64_t x = 0; x != size_x; ++x) {

                    // Loop across a ray
                    for (int64_t z = 0; z != size_z; ++z) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            if ( step_count == steps || step_count_all == size_z ) {
                                break;
                            }
                            int64_t k = sub2ind_4D_64(x + offset_x, y + offset_y, 4, t, size_max, size_max, 6);
                            *(nii_output_data + k) += *(nii_input_data + i) * 1. / steps;
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t z = size_z-1; z >= 0; --z) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            if ( step_count == steps || step_count_all == size_z ) {
                                break;
                            }
                            int64_t k = sub2ind_4D_64(x + offset_x, y + offset_y, 5, t, size_max, size_max, 6);
                            *(nii_output_data + k) += *(nii_input_data + i) * 1. / steps;
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }
        }
        cout << endl;
    }

    // ================================================================================================================
    // Compute rays with minimum intensity projection
    // ================================================================================================================
    else if ( projection_type == "min" ) {

        for (int64_t t = 0; t != size_time; ++t) {  // Loop across time points
            cout << "\r    Volume: " << t+1 << "/" << size_time << flush;

            for (int64_t z = 0; z != size_z; ++z) { // Loop across rays on x plane
                for (int64_t y = 0; y != size_y; ++y) {

                    // Loop across a ray
                    for (int64_t x = 0; x != size_x; ++x) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(y + offset_y, z + offset_z, 0, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) < *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t x = size_x-1; x >= 0; --x) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(y + offset_y, z + offset_z, 1, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) < *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }

            // --------------------------------------------------------------------
            // Loop across rays on y plane
            for (int64_t z = 0; z != size_z; ++z) {
                for (int64_t x = 0; x != size_x; ++x) {

                    // Loop across a ray
                    for (int64_t y = 0; y != size_y; ++y) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, z + offset_z, 2, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) < *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t y = size_y-1; y >= 0; --y) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, z + offset_z, 3, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) < *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }

            // --------------------------------------------------------------------
            // Loop across rays on z plane
            for (int64_t y = 0; y != size_y; ++y) {
                for (int64_t x = 0; x != size_x; ++x) {

                    // Loop across a ray
                    for (int64_t z = 0; z != size_z; ++z) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, y + offset_y, 4, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) < *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t z = size_z-1; z >= 0; --z) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, y + offset_y, 5, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) < *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }
        }
        cout << endl;
    }

    // ================================================================================================================
    // Compute rays with maximum intensity projection
    // ================================================================================================================
    else if ( projection_type == "max" ) {

        for (int64_t t = 0; t != size_time; ++t) {  // Loop across time points
            cout << "\r    Volume: " << t+1 << "/" << size_time << flush;

            for (int64_t z = 0; z != size_z; ++z) { // Loop across rays on x plane
                for (int64_t y = 0; y != size_y; ++y) {

                    // Loop across a ray
                    for (int64_t x = 0; x != size_x; ++x) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(y + offset_y, z + offset_z, 0, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) > *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t x = size_x-1; x >= 0; --x) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(y + offset_y, z + offset_z, 1, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) > *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }

            // --------------------------------------------------------------------
            // Loop across rays on y plane
            for (int64_t z = 0; z != size_z; ++z) {
                for (int64_t x = 0; x != size_x; ++x) {

                    // Loop across a ray
                    for (int64_t y = 0; y != size_y; ++y) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, z + offset_z, 2, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) > *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t y = size_y-1; y >= 0; --y) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, z + offset_z, 3, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) > *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }

            // --------------------------------------------------------------------
            // Loop across rays on z plane
            for (int64_t y = 0; y != size_y; ++y) {
                for (int64_t x = 0; x != size_x; ++x) {

                    // Loop across a ray
                    for (int64_t z = 0; z != size_z; ++z) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, y + offset_y, 4, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) > *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;

                    // Loop across a ray in reverse
                    for (int64_t z = size_z-1; z >= 0; --z) {
                        int64_t i = sub2ind_4D_64(x, y, z, t, size_x, size_y, size_z);
                        if ( *(nii_input_data + i) != 0 ) {
                            step_count += 1;
                            int64_t k = sub2ind_4D_64(x + offset_x, y + offset_y, 5, t, size_max, size_max, 6);
                            if ( step_count == 1 ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                            else if ( step_count == steps || step_count_all == size_x ) {
                                break;
                            }
                            if ( *(nii_input_data + i) > *(nii_output_data + k) ) {
                                *(nii_output_data + k) = *(nii_input_data + i);
                            }
                        }
                    }
                    step_count_all = 0;
                    step_count = 0;
                }
            }
        }
        cout << endl;
    }

    // ========================================================================
    // Save
    // ========================================================================
    // Prepare tags for the output
    std::ostringstream tag_1;
    tag_1 << steps;
    cout << "  Saving output..." << endl;
    save_output_nifti(fout, "snapcast-"+projection_type+"_steps-"+tag_1.str(), nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
