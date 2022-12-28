
#include <limits>
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_INTPRO: Do maximum and minimum intensity projections. For example,\n"
    "      this is useful for visualizing vessels.\n"
    "\n"
    "Usage:\n"
    "     LN_INTPRO -image file.nii -min -direction 3 \n"
    "     ../LN_INTPRO -image sc_UNI.nii -min -direction 2 -range 3 \n"
    "\n"
    "Options:\n"
    "    -help      : Show this help\n"
    "    -image     : Nifti (.nii) for intensity projections. \n"
    "    -max       : Maximum intensity projection. Do not combine with min.\n"
    "    -min       : Minimum intensity projection. Do not combine with max.\n"
    "    -direction : Direction in which the dimension is collapsted. \n"
    "                 1 for x, 2 for y, and 3 for z.\n"
    "    -range     : (Optional) Range of neighbouring voxels included.\n"
    "                 Default is all aslices. \n"
    "    -output    : (Optional) Output filename, including .nii or\n"
    "                 .nii.gz, and path if needed. Overwrites existing files.\n"
    "\n"
    "Notes:\n"
    "    - If the input is a time series, the entire time domain is also\n"
    "      collapsed. For example, the program chooses the minimum/maximmum \n"
    "      value across time and slices.\n"
    "    - An example application in a blog post is here: \n"
    "      <https://layerfmri.comintensity-projections-in-laynii>\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false ;
    char *fout = NULL ;
    char *fin_1 = NULL;
    int ac, is_max = 0, is_min = 0;
    int is_direction = 0, is_range = 0;
    if (argc < 2) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-direction")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -direction\n");
                return 1;
            }
            is_direction = atoi(argv[ac]);
        } else if (!strcmp(argv[ac], "-range")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -range\n");
                return 1;
            }
            is_range = atoi(argv[ac]);
        } else if (!strcmp(argv[ac], "-min")) {
            is_min = 1;
            cout << "Doing Min. intensity projections." << endl;
        } else if (!strcmp(argv[ac], "-max")) {
            is_max = 1;
            cout << "Doing Max. intensity projections." << endl;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-image")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -timeseries\n");
                return 1;
            }
            fin_1 = argv[ac];
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-timeseries'\n");
        return 1;
    }
    // Read input dataset
    nifti_image * nii_input = nifti_image_read(fin_1, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_1);
        return 2;
    }

    if (is_min+is_max !=1) {
        fprintf(stderr, "** Do either maximum or minimum. Not both.\n");
        return 1;
    }
    if (is_direction-2 >1) {
        cout << "  Invalid direction. ";
        return 1;
    }

    cout << "The signal is collapsed along: ";
    if (is_direction == 1) {
        cout << " x-direction " << endl;
    }
    if (is_direction == 2) {
        cout << " y-direction " << endl;
    }
    if (is_direction == 3) {
        cout << " z-direction " << endl;
    }

    log_welcome("LN_INTPRO");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_time = nii_input->nt;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    if (is_range > 0) {
        cout << "  Only looking for the Max/Min value in a proximity of "
             << is_range << " voxels." << endl;
    }
    if (is_range == 0) {
        is_range = max(max(size_z, size_x), size_y);
    }

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocate new nifti images
    nifti_image * nii_collapse = nifti_copy_nim_info(nii);
    nii_collapse->nt = 1;
    nii_collapse->nvox = nii->nvox / size_time;
    nii_collapse->datatype = NIFTI_TYPE_FLOAT32;
    nii_collapse->nbyper = sizeof(float);
    nii_collapse->data = calloc(nii_collapse->nvox, nii_collapse->nbyper);
    float* nii_collapse_data = static_cast<float*>(nii_collapse->data);

    double vec_file1[size_time];

    // ========================================================================
    cout << "  Starting with dimensionality collapse = " << endl;

    float extreme_val = 0.0;
    if (is_min == 1) {
        extreme_val = numeric_limits<float>::max();
    }

    for (int ix = 0; ix < size_y; ++ix) {
        for (int iy = 0; iy < size_x; ++iy) {
            for (int iz = 0; iz < size_z; ++iz) {
                int voxel_i, voxel_j;
                if (is_min == 1) extreme_val = numeric_limits<float>::max();
                if (is_max == 1) extreme_val = 0.0;
                // ------------------------------------------------------------
                if (is_direction == 1) {
                    for (int ix_i = max(ix - is_range, 0);
                         ix_i < min(ix + is_range, size_y); ++ix_i) {
                        for (int it = 0; it < size_time; ++it) {
                            voxel_i = nxyz * it + nxy * iz + nx * ix + iy;
                            voxel_j = nxyz * it + nxy * iz + nx * ix_i + iy;
                            if (is_min == 1
                                && *(nii_data + voxel_j) < extreme_val
                                && *(nii_data + voxel_j) > 0.0) {
                                extreme_val = *(nii_data + voxel_j);
                            }
                            if (is_max == 1
                                && *(nii_data + voxel_j) > extreme_val) {
                                extreme_val = *(nii_data + voxel_j);
                            }
                        }
                    }
                    for (int it = 0; it < size_time; ++it) {
                        if (is_min == 1
                            && *(nii_data + voxel_i) >= extreme_val) {
                            *(nii_collapse_data + voxel_i) = extreme_val;
                        }
                        if (is_max == 1
                            && *(nii_data + voxel_i) <= extreme_val) {
                            *(nii_collapse_data + voxel_i) = extreme_val;
                        }
                    }
                    // ------------------------------------------------------------
                } else if (is_direction == 2) {
                    for (int iy_i = max(iy - is_range, 0);
                         iy_i < min(iy + is_range, size_x); ++iy_i) {
                        for (int it = 0; it < size_time; ++it) {
                            voxel_i = nxyz * it + nxy * iz + nx * ix + iy;
                            voxel_j = nxyz * it + nxy * iz + nx * ix + iy_i;

                            if (is_min == 1
                                && *(nii_data + voxel_j) < extreme_val
                                && *(nii_data + voxel_j) > 0.0) {
                                extreme_val = *(nii_data + voxel_j);
                            }
                            if (is_max == 1
                                && *(nii_data + voxel_j) > extreme_val) {
                                extreme_val = *(nii_data + voxel_j);
                            }
                        }
                    }
                    for (int it = 0; it < size_time; ++it) {
                        if (is_min == 1
                            && *(nii_data + voxel_i) >= extreme_val) {
                            *(nii_collapse_data + voxel_i) = extreme_val;
                        }
                        if (is_max == 1
                            && *(nii_data + voxel_i) <= extreme_val) {
                            *(nii_collapse_data + voxel_i) = extreme_val;
                        }
                    }
                    // ------------------------------------------------------------
                } else if (is_direction == 3) {
                    for (int iz_i = max(iz - is_range, 0);
                         iz_i < min(iz + is_range, size_z); ++iz_i) {
                        for (int it = 0; it < size_time; ++it) {
                            voxel_i = nxyz * it + nxy * iz + nx * ix + iy;
                            voxel_j = nxyz * it + nxy * iz_i + nx * ix + iy;
                            if (is_min == 1
                                && *(nii_data + voxel_j) < extreme_val
                                && *(nii_data + voxel_j) > 0.0) {
                                extreme_val = *(nii_data + voxel_j);
                            }
                            if (is_max == 1
                                && *(nii_data + voxel_j) > extreme_val) {
                                extreme_val = *(nii_data + voxel_j);
                            }
                        }
                    }
                    for (int it = 0; it < size_time; ++it) {
                        if (is_min == 1
                            && *(nii_data + voxel_i) >= extreme_val) {
                            *(nii_collapse_data + voxel_i) = extreme_val;
                        }
                        if (is_max == 1
                            && *(nii_data + voxel_i) <= extreme_val) {
                            *(nii_collapse_data + voxel_i) = extreme_val;
                        }
                    }
                }
            }
        }
    }


    if (!use_outpath) fout = fin_1;
    save_output_nifti(fout, "collapsed", nii_collapse, true, use_outpath);

    return 0;
}
