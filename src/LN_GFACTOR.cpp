

#include "../dep/laynii_lib.h"

int N_rand = 1000;
double_t lower = -10;
double_t upper = 10;

double verteilung(double x);
typedef double ( * Functions)(double);
Functions pFunc = verteilung;

double_t arb_pdf_num(int N_rand, double ( * pFunc)(double), double_t lower,
                     double_t upper);
double adjusted_rand_numbers(double mean, double stdev, double value);

int show_help(void) {
    printf(
        "LN_GFACTOR: Simulating where the g-factor penalty would be largest.\n"
        "\n"
        "Usage:\n"
        "    LN_GFACTOR -input MEAN.nii -variance 1 -direction 1 -grappa 2 -cutoff 150 \n"
        "\n"
        "Options:\n"
        "    -help       : Show this help.\n"
        "    -input      : Specify input dataset.\n"
        "    -variance   : How much noise there will be.\n"
        "    -direction  : Phase encoding direction [0=x, 1=y, 2=z].\n"
        "    -grappa     : GRAPPA factor."
        "    -cutoff     : Value to seperate noise from signal.\n"
        "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    nifti_image* nii = NULL;
    char * fin = NULL, * fout = NULL;
    int grappa_int, direction_int, ac;
    float cutoff, variance_val;
    if (argc < 2) return show_help();

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
        } else if (!strcmp(argv[ac], "-variance")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            variance_val = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-direction")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -direction\n");
                return 1;
            }
            direction_int = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-grappa")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -grappa\n");
                return 1;
            }
            grappa_int = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-cutoff")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -cutoff\n");
                return 1;
            }
            cutoff = atof(argv[ac]);
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset
    nii = nifti_image_read(fin, 1);
    if (!nii) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_GFACTOR");
    log_nifti_descriptives(nii);

    cout << "  Variance  = " << variance_val << endl;
    cout << "  Direction = " << direction_int << endl;
    cout << "  GRAPPA    = " << grappa_int << endl;
    cout << "  Cut-off   = " << cutoff << endl;

    // Get dimensions of input
    int size_x = nii->nx;
    int size_y = nii->ny;
    int size_z = nii->nz;
    int size_t = nii->nt;
    int nx = nii->nx;
    int nxy = nii->nx * nii->ny;
    int nxyz = nii->nx * nii->ny * nii->nz;

    if (direction_int == 0) {
        size_x = size_x - (size_x % grappa_int);
    }
    if (direction_int == 1) {
        size_y = size_y - (size_y % grappa_int);
    }
    if (direction_int == 2) {
        size_z = size_z - (size_z % grappa_int);
    }

    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32(nii);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Allocating additional images
    nifti_image* nii_gfactormap = copy_nifti_as_float32(nii);
    float* nii_gfactormap_data = static_cast<float*>(nii_gfactormap->data);
    nifti_image* nii_binary = copy_nifti_as_float32(nii);
    float* nii_binary_data = static_cast<float*>(nii_binary->data);
    nifti_image* nii_noise = copy_nifti_as_float32(nii);
    float* nii_noise_data = static_cast<float*>(nii_noise->data);

    // ========================================================================

    // for (int it = 0; it < size_t; ++it) {
    //     for (int iz = 0; iz < size_z; ++iz) {
    //         for (int iy = 0; iy < size_x; ++iy) {
    //             for (int ix = 0; ix < size_y; ++ix) {
    //                 *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) = *(nii_input_data + nxyz * it + nxy * iz + nx * ix + iy) + adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper));
    //                 //cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
    //             }
    //         }
    //     }
    // }

    for (int it = 0; it < size_t; ++it) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    if (cutoff <= *(nii_input_data + nxyz * it + nxy * iz + nx * ix + iy)) {
                        *(nii_binary_data + nxyz * it + nxy * iz + nx * ix + iy) = 1;
                    } else {
                        *(nii_binary_data + nxyz * it + nxy * iz + nx * ix + iy) = 0;
                    }
                    // cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
                }
            }
        }
    }
    save_output_nifti(fin, "Gfactormap_binary", nii_binary, false);

    // Setting initial condition
    for (int it = 0; it < size_t; ++it) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) = 0.;
                    // cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
                }
            }
        }
    }

    ////////////////////////////////////
    // Adding additional GRAPPA noise //
    ////////////////////////////////////

    if (direction_int == 0) {
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int it = 0; it < size_t; ++it) {
                for (int iz = 0; iz < size_z; ++iz) {
                    for (int iy = 0; iy < size_x / grappa_int; ++iy) {
                        for (int ix = 0; ix < size_y; ++ix) {
                            *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) = *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nii_binary_data + nxyz * it + nxy * iz + nx * ix + iy + grappa_seg * size_x / grappa_int) / grappa_int;
                        }
                    }
                }
            }
        }
        // Completing dataset to full slice
        for (int grappa_seg=1; grappa_seg < grappa_int; ++grappa_seg) {
            for (int it = 0; it < size_t; ++it) {
                for (int iz = 0; iz < size_z; ++iz) {
                    for (int iy = 0; iy < size_x / grappa_int; ++iy) {
                        for (int ix = 0; ix < size_y; ++ix) {
                            *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy+ grappa_seg * size_x / grappa_int) = *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy);
                        }
                    }
                }
            }
        }
    }

    if (direction_int == 1) {
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int it = 0; it < size_t; ++it) {
                for (int iz = 0; iz < size_z; ++iz) {
                    for (int iy = 0; iy < size_x; ++iy) {
                        for (int ix = 0; ix < size_y / grappa_int; ++ix) {
                            *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) = *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nii_binary_data + nxyz * it + nxy * iz + nx * (ix+ grappa_seg * size_y / grappa_int) + iy) / grappa_int;
                        }
                    }
                }
            }
        }
        // Completing dataset to full slice
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int it = 0; it < size_t; ++it) {
                for (int iz = 0; iz < size_z; ++iz) {
                    for (int iy = 0; iy < size_x; ++iy) {
                        for (int ix = 0; ix < size_y / grappa_int; ++ix) {
                            *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * (ix+ grappa_seg * size_y / grappa_int) + iy) = *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy);
                        }
                    }
                }
            }
        }
    }

    if (direction_int == 2) {
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int it = 0; it < size_t; ++it) {
                for (int iz = 0; iz < size_z / grappa_int; ++iz) {
                    for (int iy = 0; iy < size_x; ++iy) {
                        for (int ix = 0; ix < size_y; ++ix) {
                            *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) = *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nii_binary_data + nxyz * it + nxy * (iz+ grappa_seg * size_z / grappa_int) + nx * ix + iy) / grappa_int;
                        }
                    }
                }
            }
        }
        // Completing dataset to full slice
        for (int grappa_seg = 0; grappa_seg < grappa_int; ++grappa_seg) {
            for (int it = 0; it < size_t; ++it) {
                for (int iz = 0; iz < size_z / grappa_int; ++iz) {
                    for (int iy = 0; iy < size_x; ++iy) {
                        for (int ix = 0; ix < size_y; ++ix) {
                            *(nii_gfactormap_data + nxyz * it +nxy * (iz+ grappa_seg * size_z / grappa_int) + nx * ix + iy) = *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy);
                        }
                    }
                }
            }
        }
    }

    for (int it = 0; it < size_t; ++it) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) = *(nii_binary_data + nxyz * it + nxy * iz + nx * ix + iy) **(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy);
                }
            }
        }
    }
    for (int it = 0; it < size_t; ++it) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_x; ++iy) {
                for (int ix = 0; ix < size_y; ++ix) {
                    *(nii_noise_data + nxyz * it + nxy * iz + nx * ix + iy) = *(nii_input_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nii_gfactormap_data + nxyz * it + nxy * iz + nx * ix + iy) * cutoff * adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper));
                    // cout << adjusted_rand_numbers(0, variance_val, arb_pdf_num(N_rand, pFunc, lower, upper)) << " noise    " << endl;
                }
            }
        }
    }

    save_output_nifti(fin, "Gfactormap", nii_gfactormap, true);
    save_output_nifti(fin, "Amplified_GRAPPA", nii_noise, true);

    cout << "  Finished." << endl;
    return 0;
}

// Gauss lower = -5, upper = 5
double verteilung(double z) {
    return exp(-z * z / (2.)) * 1. / sqrt(2. * M_PI);
}

double_t arb_pdf_num(int N_rand, double ( * pFunc)(double), double_t lower,
                     double_t upper) {
    double_t binwidth = (upper - lower)/(double_t)N_rand;
    double_t integral = 0.0;
    double_t rand_num = rand() / (double_t)RAND_MAX;
    int i;

    for (i = 0; integral < rand_num; i++) {
        integral += pFunc(lower + (double_t) i * binwidth) * binwidth;
        if ((lower + (double_t) i * binwidth) > upper) {
            cout << "  Upper limit, vielleicht sollte da limit angepasst werden." << i << endl;
            return lower + (double_t) i * binwidth;
        }
    }
    return lower + (double_t) i * binwidth;
}

double adjusted_rand_numbers(double mean, double stdev, double value) {
    return value * stdev + mean;
}
