
// TODO(Faruk): I think I have resolved a logical error while tidying up.
// Need to ask Renzo about this to make sure it was a bug before.

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_DIRECT_SMOOTH : Smoothing in specific directions only.\n"
    "\n"
    "Usage:\n"
    "    LN_DIRECT_SMOOTH -input activty_map.nii -FWHM 1 -direction 1 \n"
    "    LN_DIRECT_SMOOTH -input bico_VASO.Mean.nii -FWHM 0.5 -direction 3 -laurenzian \n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -input         : Nifti (.nii) file that will be smooth. It \n"
    "                     should have same dimensions as layer file.\n"
    "    -FWHM          : Amount of smoothing in unts of voxels.\n"
    "    -direction     : Axis of smoothing. 1 for x, 2 for y or 3 for z. \n"
    "    -laurenzian    : Use Laurenzian smoothing. Default is Gaussian \n"
    "                   : only for division images.\n"
    "    -Anonymous_sri : You know what you did (no FWHM).\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char* fin = NULL;
    int ac, direction = 0, option = 0;
    float FWHM_val = 10, strength = 1;
    float laur(float distance, float sigma);
    float ASLFt(float distance, float strength);
    if (argc < 3) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-direction")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -direction\n");
                return 1;
            }
            direction = atoi(argv[ac]);
        } else if (!strcmp(argv[ac], "-laurenzian")) {
            option = 1;
        } else if (!strcmp(argv[ac], "-Anonymous_sri")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -Anonymous_sri\n");
                // return 1;
            }
            option = 2;
            strength = atof(argv[ac]);
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(fin, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }
    if (direction == 0) {
        fprintf(stderr, "** invalid direction '%i'\n", direction);
        return 2;
    }

    log_welcome("LN_DIRECT_SMOOTH");
    log_nifti_descriptives(nii1);

    if (option == 0) {
        cout << "  Gaussian smoothing is selected." << endl;;
    } else if (option == 1) {
        cout << "  Use Laurenzian smoothing instead of Gaussian." << endl;;
    } else if (option == 2) {
        cout << "  Sri option is selected." << endl;;
        cout << "    No FWHM_val is used. Strength is " << strength << endl;
    }

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;
    const int size_t = nii1->nt;
    const int nx = nii1->nx;
    const int nxy = nii1->nx * nii1->ny;
    const int nxyz = nii1->nx * nii1->ny * nii1->nz;
    // TODO(Faruk): Need to ask to Renzo about kernel symmetry
    const float dX = 1;  // nii1->pxdim[1];
    const float dY = 1;  // nii1->pxdim[2];
    const float dZ = 1;  // nii1->pxdim[3];

    // ========================================================================
    // Fx datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Allocate new nifti images
    nifti_image* smooth = copy_nifti_as_float32(nii1);
    float* smooth_data = static_cast<float*>(smooth->data);

    // ========================================================================
    int vic = max(1., 2. * FWHM_val / dX);  // Vicinity, ignore far voxels
    cout << "    vic = " << vic << endl;
    cout << "    FWHM = " << FWHM_val << endl;

    // ========================================================================
    // Smoothing loop
    // ========================================================================
    cout << "  Smoothing dimension = " << direction << endl;

    for (int t = 0; t < size_t; ++t) {
        for (int z = 0; z < size_z; ++z) {
            for (int y = 0; y < size_y; ++y) {
                for (int x = 0; x < size_x; ++x) {
                    int voxel_i = nxyz * t + nxy * z + nx * y + x;
                    *(smooth_data + voxel_i) = 0;

                    float total_weight = 0;
                    int start_j, stop_j;

                    if (direction == 1) {
                        start_j = max(0, x - vic);
                        stop_j = min(x + vic, size_x);
                    } else if (direction == 2) {
                        start_j = max(0, y - vic);
                        stop_j = min(y + vic, size_y);
                    } else if (direction == 3) {
                        start_j = max(0, z - vic);
                        stop_j = min(z + vic, size_z);
                    }

                    for (int j = start_j; j <= stop_j; ++j) {
                        int voxel_j;
                        float d, w;
                        if (direction == 1) {
                            voxel_j = nxyz * t + nxy * z + nx * y + j;
                            d = dist((float)x, (float)y, (float)z,
                                     (float)j, (float)y, (float)z,
                                     dX, dY, dZ);
                        } else if (direction == 2) {
                            voxel_j = nxyz * t + nxy * z + nx * j + x;
                            d = dist((float)x, (float)y, (float)z,
                                     (float)x, (float)j, (float)z,
                                     dX, dY, dZ);
                        } else if (direction == 3) {
                            voxel_j = nxyz * t + nxy * j + nx * y + x;
                            d = dist((float)x, (float)y, (float)z,
                                     (float)x, (float)y, (float)j,
                                     dX, dY, dZ);
                        }
                        // TODO(Faruk): Need to ask to Renzo about j masking
                        // Potentially problematic with sulci + big FWHM
                        if (*(nii_input_data + voxel_j) != 0) {
                            if (option == 0) {
                                w = gaus(d, FWHM_val);
                            } else if (option == 1) {
                                w = laur(d, FWHM_val);
                            } else if (option == 2) {
                                w = ASLFt(d, strength);
                            }
                            *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * w;
                            total_weight += w;
                        }
                    }
                    if (total_weight != 0) {
                        *(smooth_data + voxel_i) /= total_weight;
                    }
                }
            }
        }
    }
    save_output_nifti(fin, "smooth", smooth, true);

    cout << "  Finished." << endl;
    return 0;
}

float laur(float distance, float sigma) {
    // Note: For consistency wth Gaus's sigma, I am using a scaled version of
    // the FWHM.
    // sigma = sigma / sqrt(2 * log (2));
    float result = 0;
    if ((int)distance % 2 == 0) {
        // sigma = 2 * sigma;
        result = 1 * 1. / (3.141592 * sigma) * 1
                 / (1 + distance * distance / (sigma * sigma));
    }
    if ((int)distance % 2 == 1) {
        // sigma = 2 * sigma;
        result = 1.6 * 1. / (3.141592 * sigma) * 1
                 / (1 + distance * distance / (sigma * sigma));
    }
    return result;
// return 1 / 3.141592 * 1 / ((distance) * (distance) + (0.5 * sigma) * (0.5 * sigma));
}

float ASLFt(float distance, float strength) {
    float value;
    if (distance >= 4.0) value = 0.0;
    if (distance <= 4.0) value = 0.0245;
    if (distance <= 3.0) value = 0.024;
    if (distance <= 2.0) value = 0.06;
    if (distance <= 1.0) value = 0.18;
    if (distance == 0.0) value = 1.0;
    return value * strength;
}
