

#include "./laynii_lib.h"

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
    "    -input         : Nifti (.nii) file that will be smoothed. It \n"
    "                     should have same dimensions as layer file.\n"
    "    -FWHM          : Amount of smoothing in units of voxels.\n"
    "    -laurenzian    : Use Laurenzian smoothing. Default is Gaussian \n"
    "                   : only for division images.\n"
    "    -direction     : Axis of smoothing. 1 for x, 2 for y or 3 for z. \n"
    "    -Anonymous_sri : You know what you did (no FWHM).\n"
    "\n"
    "Note: This program supports INT16, INT32 and FLOAT32.\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char* finfi = NULL;
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
            finfi = argv[ac];  // Assign pointer, no string copy.
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

    if (!finfi) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(finfi, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", finfi);
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

    // ========================================================================
    // Fix datatype issues

    nifti_image* nii1_temp = copy_nifti_header_as_float(nii1);
    float* nii1_temp_data = static_cast<float*>(nii1_temp->data);

    // ========================================================================

    /////////////////////////////////////
    // Make allocating necessary files //
    /////////////////////////////////////
    nifti_image* smoothed = nifti_copy_nim_info(nii1_temp);
    nifti_image* gausweight = nifti_copy_nim_info(nii1_temp);
    // nifti_image* layer = nifti_copy_nim_info(nim_input);
    // nifti_image* leak_layer = nifti_copy_nim_info(nim_input);

    smoothed->datatype = NIFTI_TYPE_FLOAT32;
    gausweight->datatype = NIFTI_TYPE_FLOAT32;
    // layer->datatype = NIFTI_TYPE_FLOAT32;
    // leak_layer->datatype = NIFTI_TYPE_FLOAT32;

    smoothed->nbyper = sizeof(float);
    gausweight->nbyper = sizeof(float);
    // layer->nbyper = sizeof(float);
    // leak_layer->nbyper = sizeof(float);

    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);
    // layer->data = calloc(layer->nvox, layer->nbyper);
    // leak_layer->data = calloc(leak_layer->nvox, leak_layer->nbyper);

    float* smoothed_data = (float*) smoothed->data;
    float* gausweight_data = (float*) gausweight->data;
    // float* layer_data = (float*) layer->data;
    // float* leak_layer_data = (float*) leak_layer->data;
    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);
    float gaus(float distance, float sigma);
    float laur(float distance, float sigma);
    float ASLFt(float distance, float strength);

    // Float kernal_size = 10;  // corresponds to one voxel sice.
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

    if (direction_i == 1) {
        FOR_EACH_VOXEL_TZYX
            *(gausweight_data + VOXEL_ID) = 0;

            for (int ix_i = max(0, ix - vinc); ix_i < min(ix + vinc + 1, size_x); ++ix_i) {
                if (*(nii1_temp_data + nxyz * it + nxy * iz + nx * iy + ix_i) != 0) {
                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy, (float)iz, dX, dY, dZ);

                    if (do_gauss == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz + nx * iy + ix_i) * gaus(dist_i, FWHM_val);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + gaus(dist_i, FWHM_val);
                    }
                    if (do_laurenz == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz + nx * iy + ix_i) * laur(dist_i, FWHM_val);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + laur(dist_i, FWHM_val);
                    }
                    if (do_sri == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz + nx * iy + ix_i) * ASLFt(dist_i, strength);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + ASLFt(dist_i, strength);
                    }
                }
            }
        END_FOR_EACH_VOXEL_TZYX
    }
    if (direction_i == 2) {
        FOR_EACH_VOXEL_TZYX
            *(gausweight_data + VOXEL_ID) = 0;

            for (int iy_i = max(0, iy - vinc); iy_i < min(iy + vinc + 1, size_y); ++iy_i) {
                if (*(nii1_temp_data + nxyz * it + nxy * iz + nx * iy_i + ix) != 0) {
                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix, (float)iy_i, (float)iz, dX, dY, dZ);

                    if (do_gauss == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz + nx * iy_i + ix) * gaus(dist_i, FWHM_val);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + gaus(dist_i, FWHM_val);
                    }
                    if (do_laurenz == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz + nx * iy_i + ix) * laur(dist_i, FWHM_val);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + laur(dist_i, FWHM_val);
                    }
                    if (do_sri == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz + nx * iy_i + ix) * ASLFt(dist_i, strength);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + ASLFt(dist_i, strength);
                    }
                }
            }
        END_FOR_EACH_VOXEL_TZYX
    }
    if (direction_i == 3) {
        FOR_EACH_VOXEL_TZYX
            *(gausweight_data + VOXEL_ID) = 0;

            for (int iz_i = max(0, iz - vinc); iz_i < min(iz + vinc + 1, size_z); ++iz_i) {
                if (*(nii1_temp_data + nxyz * it + nxy * iz_i + nx * iy + ix) != 0) {
                    dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix, (float)iy, (float)iz_i, dX, dY, dZ);

                    if (do_gauss == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz_i + nx * iy + ix) * gaus(dist_i, FWHM_val);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + gaus(dist_i, FWHM_val);
                    }
                    if (do_laurenz == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz_i + nx * iy + ix) * laur(dist_i, FWHM_val);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + laur(dist_i, FWHM_val);
                    }
                    if (do_sri == 1) {
                        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) + *(nii1_temp_data + nxyz * it + nxy * iz_i + nx * iy + ix) * ASLFt(dist_i, strength);
                        *(gausweight_data + VOXEL_ID) = *(gausweight_data + VOXEL_ID) + ASLFt(dist_i, strength);
                    }
                }
            }
        END_FOR_EACH_VOXEL_TZYX
    }

    ///////////////////////////////
    // Correcting for edge error //
    ///////////////////////////////
    FOR_EACH_VOXEL_TZYX
        *(smoothed_data + VOXEL_ID) = *(smoothed_data + VOXEL_ID) / *(gausweight_data + VOXEL_ID);
    END_FOR_EACH_VOXEL_TZYX

    smoothed->scl_slope = nii1->scl_slope;

    if (nii1->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    // output file name
    // const char*fout_4 ="leaky_layers.nii";
    // if (nifti_set_filenames(leak_layer, fout_4, 1, 1)) return 1;
    // nifti_image_write(leak_layer);
    //
    // const char*fout_5 = "input_file.nii";
    // if (nifti_set_filenames(nii1_temp, fout_5, 1, 1)) return 1;
    // nifti_image_write(nii1_temp);
    //
    // const char*fout_2 = "mask.nii";
    // if (nifti_set_filenames(nim_mask, fout_2, 1, 1)) return 1;
    // nifti_image_write(nim_mask);

    string prefix = "smoothed_";
    string filename = (string) (finfi);
    string outfilename = prefix + filename;
    log_output(outfilename.c_str());

    const char*fout_1 = outfilename.c_str();
    if (nifti_set_filenames(smoothed, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(smoothed);

    // const char*fout_1 = "layer.nii";
    // if (nifti_set_filenames(layer, fout_1, 1, 1)) return 1;
    // nifti_image_write(layer);

    cout << "  Finished." << endl;
    return 0;
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ) {
    return sqrt((x1 - x2) * (x1-x2) * dX * dX
                + (y1 - y2) * (y1 - y2) * dY * dY
                + (z1 - z2) * (z1 - z2) * dZ * dZ);
}

float gaus(float distance, float sigma) {
    return 1. / (sigma * sqrt(2. * 3.141592))
           * exp(-0.5 * distance * distance / (sigma * sigma));
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
