
#include "./common.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_LEAKY_LAYERS: Layering algorithm based on iterative smothing.\n"
    "\n"
    "    This program derives 20 layers from GM/WM and GM/CSF border lines.\n"
    "\n"
    "Usage:\n"
    "    LN_LEAKY_LAYERS -rim rim.nii -dim 2 \n"
    "\n"
    "Options:\n"
    "    -help               : Show this help.\n"
    "    -disp_float_example : Show some voxel's data.\n"
    "    -rim border         : Specify input dataset INT16.\n"
    "    -dim value          : Specify value (2 or 3) layer algorithm. \n"
    "                          Default is 3 (3D).\n"
    "\n"
    "Notes:\n"
    "    - This program now supports INT16, INT32 and FLOAT32. For other \n"
    "      data types use: \n"
    "          3dcalc -a rim.nii -datum short -expr 'a' -prefix rim_short.nii\n"
    "    - This can be 3D. Hence the rim file should be dmsmooth in all three\n"
    "      dimensions.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL, *fout = NULL;
    int ac, disp_float_eg = 0, dim;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options.
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-disp_float_example")) {
            disp_float_eg = 1;
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-dim")) {
            if (++ac >= argc) {
                fprintf(stderr, " * * missing argument for -input\n");
                return 1;
            }
            dim = atof(argv[ac]);  // Assign pointer, no string copy
        } else {
            fprintf(stderr, " * * invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, " * * missing option '-rim'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_inputr = nifti_image_read(fin, 1);
    if (!nim_inputr) {
        fprintf(stderr, " * * failed to read NIfTI image from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_LEAKY_LAYERS");
    log_nifti_descriptives(nim_inputr);

    // Get dimensions of input
    int sizeSlice = nim_inputr->nz;
    int sizePhase = nim_inputr->nx;
    int sizeRead = nim_inputr->ny;
    int nrep = nim_inputr->nt;
    int nx = nim_inputr->nx;
    int nxy = nim_inputr->nx * nim_inputr->ny;
    int nxyz = nim_inputr->nx * nim_inputr->ny * nim_inputr->nz;
    float dX = nim_inputr->pixdim[1];
    float dY = nim_inputr->pixdim[2];
    float dZ = nim_inputr->pixdim[3];

    nifti_image * nim_input = nifti_copy_nim_info(nim_inputr);
    nim_input->datatype = NIFTI_TYPE_FLOAT32;
    nim_input->nbyper = sizeof(float);
    nim_input->data = calloc(nim_input->nvox, nim_input->nbyper);
    float * nim_input_data = (float *) nim_input->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////
    if (nim_inputr->datatype == NIFTI_TYPE_FLOAT32) {
        float * nim_inputr_data = (float *) nim_inputr->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputr_data + nxyz * it + nxy * islice + nx * ix + iy));
                        if ((*(nim_inputr_data + nxyz * it + nxy * islice + nx * ix + iy)) != 0  && (*(nim_inputr_data + nxyz * it + nxy * islice + nx * ix + iy)) <= 4) {
                            *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputr_data + nxyz * it + nxy * islice + nx * ix + iy));
                        }
                        else *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = 0;
                    }
                }
            }
        }
    }

    if (nim_inputr->datatype == NIFTI_TYPE_INT16) {
        short *nim_inputr_data = (short *) nim_inputr->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputr_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    if (nim_inputr->datatype == NIFTI_TYPE_INT32) {
        int *nim_inputr_data = (int *) nim_inputr->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputr_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    // cout << "  nim_input->intent_code " << nim_input->intent_code << endl;

    ////////////////////////////
    // Make it 2D if you want //
    ////////////////////////////
    if (dim == 2) {
        dZ = 1000.;
        cout << "  Calculating layers only in 2D." << endl;
    }

    nifti_image *smoothed = nifti_copy_nim_info(nim_input);
    nifti_image *gausweight = nifti_copy_nim_info(nim_input);
    nifti_image *layer = nifti_copy_nim_info(nim_input);
    nifti_image *leak_layer = nifti_copy_nim_info(nim_input);

    smoothed->datatype = NIFTI_TYPE_FLOAT32;
    gausweight->datatype = NIFTI_TYPE_FLOAT32;
    layer->datatype = NIFTI_TYPE_FLOAT32;
    leak_layer->datatype = NIFTI_TYPE_FLOAT32;

    smoothed->nbyper = sizeof(float);
    gausweight->nbyper = sizeof(float);
    layer->nbyper = sizeof(float);
    leak_layer->nbyper = sizeof(float);

    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);
    layer->data = calloc(layer->nvox, layer->nbyper);
    leak_layer->data = calloc(leak_layer->nvox, leak_layer->nbyper);

    float *smoothed_data = (float *) smoothed->data;
    float *gausweight_data = (float *) gausweight->data;
    float *layer_data = (float *) layer->data;
    float *leak_layer_data = (float *) leak_layer->data;

    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);
    float gaus(float distance, float sigma);

    cout << "  Debug 2 " << endl;
    float kernel_size = 2;  // Corresponds to one voxel sice.
    int vinc = max(1., 2. * kernel_size);  // Ignore If voxel is too far
    float dist_i = 0.;
    cout << "  vinc " << vinc << endl;
    cout << "  kernel_size " << kernel_size << endl;

    int number_of_layers = 20;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 1) {
                     *(layer_data + nxy * iz + nx * ix + iy) = -200.;
                }
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 2) {
                    *(layer_data + nxy * iz + nx * ix + iy) = 200.;
                }
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) {
                    *(layer_data + nxy * iz + nx * ix + iy) = 0.;
                }
            }
        }
    }

    ///////////////////////////////
    // Start iterative loop here //
    ///////////////////////////////
    int N_iteratiosn = 400;

    float dummy = 0;
    for (int iteration = 0; iteration < N_iteratiosn; ++iteration) {
        cout << "\r  Iteration: " << iteration << " of " << N_iteratiosn << flush;
        // cout << "  Debug 3 = " << gaus(1 , 3) << endl;
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    *(gausweight_data + nxy * iz + nx * ix + iy) = 0;
                    *(smoothed_data + nxy * iz + nx * ix + iy) = 0;
                    if (*(nim_input_data + nxy * iz + nx * ix + iy)  > 0) {
                        for (int iz_i = max(0, iz - vinc); iz_i < min(iz + vinc + 1, sizeRead); ++iz_i) {
                            for (int iy_i = max(0, iy - vinc); iy_i < min(iy + vinc + 1, sizePhase); ++iy_i) {
                                for (int ix_i = max(0, ix - vinc); ix_i < min(ix + vinc + 1, sizeRead); ++ix_i) {
                                    if (*(nim_input_data + nxy * iz_i + nx * ix_i + iy_i)  > 0) {
                                        dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                        // cout << "debug  4 " << gaus(dist_i , kernel_size)  endl;
                                        // if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) cout << "  Debug 4b " << endl;
                                        if (dist_i < vinc && *(nim_input_data + nxy * iz + nx * ix + iy) == 3) {
                                            //dummy = * (layer_data + nxy * iz_i + nx * ix_i + iy_i);
                                            *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy) + *(layer_data + nxy * iz_i + nx * ix_i + iy_i) * gaus(dist_i, kernel_size);
                                            *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy) + gaus(dist_i, kernel_size);
                                        }
                                    }
                                }
                            }
                        }
                        if (*(gausweight_data + nxy * iz + nx * ix + iy) >= 0) {
                            *(smoothed_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy)/ *(gausweight_data + nxy * iz + nx * ix + iy);
                        }
                    }
                }
            }
        }
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 1) {
                        *(layer_data + nxy * iz + nx * ix + iy) = -200;
                    }
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 2) {
                        *(layer_data + nxy * iz + nx * ix + iy) = 200;
                    }
                    if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) {
                        *(layer_data + nxy * iz + nx * ix + iy) = *(smoothed_data + nxy * iz + nx * ix + iy);
                    }
                }
            }
        }
    }  // Iteration loop closed

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                *(leak_layer_data + nxy * iz + nx * ix + iy) = 2+19 * (*(smoothed_data + nxy * iz + nx * ix + iy) -(-200))/(200-(-200));
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 1) {
                    *(leak_layer_data + nxy * iz + nx * ix + iy) = 1;
                }
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 2) {
                    *(leak_layer_data + nxy * iz + nx * ix + iy) = 21;
                }
                if (*(nim_input_data + nxy * iz + nx * ix + iy) == 0) {
                    *(leak_layer_data + nxy * iz + nx * ix + iy) = 0;
                }
                if (*(leak_layer_data + nxy * iz + nx * ix + iy) < 0) {
                    *(leak_layer_data + nxy * iz + nx * ix + iy) = 0;
                }
            }
        }
    }

    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (*(nim_input_data + nxy * iz + nx * ix + iy) != 0) {
                    *(leak_layer_data + nxy * iz + nx * ix + iy) = (int) (22- *(leak_layer_data + nxy * iz + nx * ix + iy));
                }
                if (*(leak_layer_data + nxy * iz + nx * ix + iy) == 20) {
                    *(leak_layer_data + nxy * iz + nx * ix + iy) = 19;
                }
                if (*(leak_layer_data + nxy * iz + nx * ix + iy) == 21) {
                    *(leak_layer_data + nxy * iz + nx * ix + iy) = 20;
                }
            }
        }
    }
    cout << "  Running also until here 5..." << endl;

    // Output file name
    const char *fout_4 = "leaky_layers.nii";
    if (nifti_set_filenames(leak_layer, fout_4, 1, 1)) {
        return 1;
    }
    nifti_image_write(leak_layer);

    // const char *fout_5 = "input.nii";
    // if (nifti_set_filenames(nim_input, fout_5 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(nim_input);
    //
    // const char *fout_2 = "gausweights.nii";
    // if (nifti_set_filenames(gausweight, fout_2 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(gausweight);
    //
    // const char *fout_1 = "layer.nii";
    // if (nifti_set_filenames(layer, fout_1 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(layer);

    cout << "  Finished." << endl;
    return 0;
}

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
            float dX, float dY, float dZ) {
    return sqrt((x1 - x2) * (x1 - x2) * dX * dX
                + (y1 - y2) * (y1 - y2) * dY * dY
                + (z1 - z2) * (z1 - z2) * dZ * dZ);
}

float gaus(float distance, float sigma) {
    return (1. / (sigma * sqrt(2. * 3.141592))
            * exp(-0.5 * distance * distance / (sigma * sigma)));
}
