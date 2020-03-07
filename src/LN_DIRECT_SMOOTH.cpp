

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_DIRECT_SMOOTH : Smoothing in specific directions only.\n"
    "\n"
    "    This program smooths data within layer or columns. In order to\n"
    "    avoid smoothing across masks a crawler smooths only across connected "
    "    voxels.\n"
    "\n"
    "Usage:\n"
    "    LN_DIRECT_SMOOTH -input activity_map.nii -FWHM 1 -direction 1 \n"
    "    LN_DIRECT_SMOOTH -input bico_VASO.Mean.nii -FWHM 0.5 -direction 3 -laurenzian \n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -input         : Nifti (.nii) file that will be smooth. It \n"
    "                     should have same dimensions as layer file.\n"
    "    -FWHM          : Amount of smoothing in units of voxels.\n"
    "    -direction     : Axis of smoothing. 1 for x, 2 for y or 3 for z. \n"
    "    -laurenzian    : Use Laurenzian smoothing. Default is Gaussian \n"
    "                   : only for division images.\n"
    "    -Anonymous_sri : You know what you did (no FWHM).\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char* fin_1 = NULL;
    int ac, direction_i = 0, do_laurenz = 0, do_gauss = 0, do_sri = 0;
    float FWHM_val = 10, strength = 1;
    if (argc < 3) {
        return show_help();  // Typing '-help' is sooo much work
    }
    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);  // Assign pointer, no string copy.
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin_1 = argv[ac];  // Assign pointer, no string copy.
        } else if (!strcmp(argv[ac], "-direction")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -direction\n");
                return 1;
            }
            direction_i = atoi(argv[ac]);  // Assign pointer, no string copy.
        } else if (!strcmp(argv[ac], "-laurenzian")) {
            do_laurenz = 1;
            fprintf(stderr, "Use Laurenzian smoothing instead of Gaussian.");
        } else if (!strcmp(argv[ac], "-Anonymous_sri")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -Anonymous_sri\n");
                // return 1;
            }
            do_sri = 1;
            strength = atof(argv[ac]);  // Assign pointer, no string copy.
            fprintf(stderr, "Yes Sri I am doing you.");
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", fin_1);
        return 2;
    }
    if (direction_i == 0) {
        fprintf(stderr, "** failed to read direction '%i'\n", direction_i);
        return 2;
    }

    log_welcome("LN_DIRECT_SMOOTH");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;
    const int size_t = nii1->nt;
    const int nx = nii1->nx;
    const int nxy = nii1->nx * nii1->ny;
    const int nxyz = nii1->nx * nii1->ny * nii1->nz;
    const float dX = 1;  // nii1->pixdim[1];
    const float dY = 1;  // nii1->pixdim[2];
    const float dZ = 1;  // nii1->pixdim[3];
    const int nr_voxels = size_t * size_z * size_y * size_x;


    // ========================================================================
    // Fix datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Allocate new nifti images
    nifti_image* smooth = copy_nifti_as_float32(nii_input);
    float* smooth_data = static_cast<float*>(smooth->data);
    nifti_image* gaussw = copy_nifti_as_float32(nii_input);
    float* gaussw_data = static_cast<float*>(gaussw->data);

    // ========================================================================

    // float* layer_data = (float*) layer->data;
    // float* leak_layer_data = (float*) leak_layer->data;
    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);
    float gaus(float distance, float sigma);
    float laur(float distance, float sigma);
    float ASLFt(float distance, float strength);

    // Float kernel_size = 10;  // corresponds to one voxel sice.
    int vinc = max(1., 2. * FWHM_val / dX);  // Ignore far voxels
    float dist_i = 0.;
    cout << "  vinc " << vinc<< endl;
    cout << "  FWHM_val " << FWHM_val<< endl;

    if (do_sri == 0 && do_laurenz == 0) {
        do_gauss = 1;
    }
    if (do_sri == 1) {
        cout << "  No FWHM_val is used. Strength is " << strength << endl;
    }

    ////////////////////
    // Smoothing loop //
    ////////////////////
    // cout << " DEBUG " << dist(1.,1.,1.,1.,2.,1.,dX,dY,dZ) << endl;
    cout << "  Smoothing in dimension " << direction_i << endl;
    cout << "  Time step... " << flush;

    for (int i = 0; i < 5; i++) {
        cout << "  " << ASLFt(i, strength) << endl;
    }

    for (int it = 0; it < size_t; ++it) {
        for (int iz = 0; iz < size_z; ++iz) {
            for (int iy = 0; iy < size_y; ++iy) {
                for (int ix = 0; ix < size_x; ++ix) {
                    int voxel_i = nxyz * it + nxy * iz + nx * iy + ix;
                    *(gaussw_data + voxel_i) = 0;
                    if (direction_i == 1) {
                        for (int ix_i = max(0, ix - vinc); ix_i < min(ix + vinc + 1, size_x); ++ix_i) {
                            int voxel_j = nxyz * it + nxy * iz + nx * iy + ix_i;
                            if (*(nii_input_data + voxel_j) != 0) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy, (float)iz, dX, dY, dZ);
                                if (do_gauss == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * gaus(dist_i, FWHM_val);
                                    *(gaussw_data + voxel_i) += gaus(dist_i, FWHM_val);
                                }
                                if (do_laurenz == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * laur(dist_i, FWHM_val);
                                    *(gaussw_data + voxel_i) += laur(dist_i, FWHM_val);
                                }
                                if (do_sri == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * ASLFt(dist_i, strength);
                                    *(gaussw_data + voxel_i) += ASLFt(dist_i, strength);
                                }
                            }
                        }
                    } else if (direction_i == 2) {
                        for (int iy_i = max(0, iy - vinc); iy_i < min(iy + vinc + 1, size_y); ++iy_i) {
                            int voxel_j = nxyz * it + nxy * iz + nx * iy_i + ix;
                            if (*(nii_input_data + voxel_j) != 0) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix, (float)iy_i, (float)iz, dX, dY, dZ);

                                if (do_gauss == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * gaus(dist_i, FWHM_val);
                                    *(gaussw_data + voxel_i) += gaus(dist_i, FWHM_val);
                                }
                                if (do_laurenz == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * laur(dist_i, FWHM_val);
                                    *(gaussw_data + voxel_i) += laur(dist_i, FWHM_val);
                                }
                                if (do_sri == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * ASLFt(dist_i, strength);
                                    *(gaussw_data + voxel_i) += ASLFt(dist_i, strength);
                                }
                            }
                        }
                    } else if (direction_i == 3) {
                        for (int iz_i = max(0, iz - vinc); iz_i < min(iz + vinc + 1, size_z); ++iz_i) {
                            int voxel_j = nxyz * it + nxy * iz_i + nx * iy + ix;
                            if (*(nii_input_data + voxel_j) != 0) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix, (float)iy, (float)iz_i, dX, dY, dZ);

                                if (do_gauss == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * gaus(dist_i, FWHM_val);
                                    *(gaussw_data + voxel_i) += gaus(dist_i, FWHM_val);
                                }
                                if (do_laurenz == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * laur(dist_i, FWHM_val);
                                    *(gaussw_data + voxel_i) += laur(dist_i, FWHM_val);
                                }
                                if (do_sri == 1) {
                                    *(smooth_data + voxel_i) += *(nii_input_data + voxel_j) * ASLFt(dist_i, strength);
                                    *(gaussw_data + voxel_i) += ASLFt(dist_i, strength);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ///////////////////////////////
    // Correcting for edge error //
    ///////////////////////////////
    for (int i = 0; i < nr_voxels; ++i) {
        *(smooth_data + i) /= *(gaussw_data + i);
    }

    smooth->scl_slope = nii1->scl_slope;

    if (nii1->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    save_output_nifti(fin_1, "smooth", smooth, true);

    cout << "  Finished." << endl;
    return 0;
}

float laur(float distance, float sigma) {
    // Note: For consistency with Gaus's sigma, I am using a scaled version of
    // the FWHM.
    // sigma = sigma / sqrt(2 * log (2));
    float result = 0;
    if ((int)distance%2 == 0) {
        // sigma = 2 * sigma;
        result = 1 * 1. / (3.141592 * sigma) * 1
                 / (1 + distance * distance / (sigma * sigma));
    }
    if ((int)distance%2 == 1) {
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
