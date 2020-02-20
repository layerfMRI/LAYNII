

#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_TEMPSMOOTH : Temporal filtering.\n"
    "\n"
    "    Smooths data within the time domain. It removes high frequencies \n"
    "    spikes like a low-pass filter.\n"
    "\n"
    "Usage:\n"
    "    LN_TEMPSMOOTH -timeseries file.nii -gaus 1.0 \n"
    "    LN_TEMPSMOOTH -timeseries file.nii -box 1 \n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -timeseries : Nifti (.nii) file that with the series that will be \n"
    "                  smoothed. Only the first time point is used. \n"
    "    -gaus value : Doing the smoothing with a Gaussian weight function. \n"
    "                  A travelling window of averaging. Specify the value \n"
    "                  of the Gaussian size (float values) in units of TR. \n"
    "    -box value  : Doing the smoothing with a box-var. Specify the value \n"
    "                  of the box sice (integer value). This is like a \n"
    "                  running average sliding window.\n"
    "\n"
    "Note: This program now supports INT16, INT32 and FLOAT32.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fmaski = NULL, *fout = NULL, *finfi = NULL;
    int ac, do_gaus = 0, do_box = 0, bFWHM_val = 0;
    float gFWHM_val = 0.0;
    if (argc  <  3) {  // Typing '-help' is sooo much work
        return show_help();
    }
    // Process user options
    for (ac = 1; ac  <  argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-gaus")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -gaus\n");
                return 1;
            }
            gFWHM_val = atof(argv[ac]);  // Assign pointer, no string copy
            do_gaus = 1;
            cout << "Selected temporal smoothing: Gaussian" << endl;
        } else if (!strcmp(argv[ac], "-box")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -box\n");
                return 1;
            }
            bFWHM_val = atoi(argv[ac]);  // Assign pointer, no string copy
            do_box = 1;
            cout << "Selected temporal smoothing: Box-car" << endl;
        } else if (!strcmp(argv[ac], "-timeseries")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            finfi = argv[ac];  // Assign pointer, no string copy
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!finfi) {
        fprintf(stderr, "** missing option '-timeseries'\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image * nim_inputfi = nifti_image_read(finfi, 1);
    if (!nim_inputfi) {
        fprintf(stderr, "** failed to read layer NIfTI image from '%s'\n", finfi);
        return 2;
    }
    if (do_box + do_gaus !=1) {
        cout << "  Invalid smoothing option. Please select Gaussian or box-car." << endl;
        return 2;
    }

    log_welcome("LN_TEMPSMOOTH");
    log_nifti_descriptives(nim_inputfi);

    // Get dimensions of input
    int sizeSlice = nim_inputfi->nz;
    int sizePhase = nim_inputfi->nx;
    int sizeRead = nim_inputfi->ny;
    int nrep = nim_inputfi->nt;
    int nx = nim_inputfi->nx;
    int nxy = nim_inputfi->nx * nim_inputfi->ny;
    int nxyz = nim_inputfi->nx * nim_inputfi->ny * nim_inputfi->nz;
    float dX = nim_inputfi->pixdim[1];
    // float dY = nim_inputfi->pixdim[2];
    // float dZ = nim_inputfi->pixdim[3];

    // Note: If you are running the smoothing in 2D, it will still go through
    // the entire pipeline. The only difference is that the weights in a
    // certain direction are suppressed. Doing it in 2D, will be faster.

    nifti_image * nim_inputf = nifti_copy_nim_info(nim_inputfi);
    nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
    nim_inputf->nbyper = sizeof(float);
    nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
    float* nim_inputf_data = (float* ) nim_inputf->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    // here, I am loading them in their native datatype //////////
    // and translate them to the datatime I like best  ///////////
    //////////////////////////////////////////////////////////////
    if (nim_inputfi->datatype == NIFTI_TYPE_FLOAT32) {
        float* nim_inputfi_data = (float* ) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_inputfi->datatype == NIFTI_TYPE_INT16) {
        short *nim_inputfi_data = (short *) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz *it + nxy*islice + nx*ix + iy) = (float) (*(nim_inputfi_data + nxyz *it + nxy*islice + nx*ix + iy));
                    }
                }
            }
        }
    }
    if (nim_inputfi->datatype == NIFTI_TYPE_INT32) {
        int *nim_inputfi_data = (int *) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }

    /////////////////////////////////////
    // Make allocating necessary files //
    /////////////////////////////////////
    nifti_image * smoothed = nifti_copy_nim_info(nim_inputf);
    nifti_image * gausweight = nifti_copy_nim_info(nim_inputf);
    smoothed->datatype = NIFTI_TYPE_FLOAT32;
    gausweight->datatype = NIFTI_TYPE_FLOAT32;
    smoothed->nbyper = sizeof(float);
    gausweight->nbyper = sizeof(float);
    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);

    // Note: Gauss weight is a geometry factor. Only needs to be estimated once.
    // So with the next lines I am saving RAM.
    gausweight->nt = 1;
    gausweight->nvox = gausweight->nvox / nrep;
    float* smoothed_data = (float* ) smoothed->data;
    float* gausweight_data = (float* ) gausweight->data;

    // nifti_image *debug = nifti_copy_nim_info(nim_inputf);
    // debug->datatype = NIFTI_TYPE_FLOAT32;
    // debug->nbyper = sizeof(float);
    // debug->data = calloc(debug->nvox, debug->nbyper);
    // float* debug_data = (float* ) debug->data;

    // float dist(float x1, float y1, float z1, float x2, float y2, float z2,
    //            float dX, float dY, float dZ);
    float gaus(float distance, float sigma);

    cout << "  Debug 2 " << endl;

    ////////////////////////////////////////////
    // Smoothing loop (for Gaussian smoothing //
    ////////////////////////////////////////////
    if (do_gaus) {
        // float_kernel_size = 10;  // Corresponds to one voxel size.
        int vinc = max(1., 2. * gFWHM_val / dX);  // Ignore if voxel is too far
        float dist_i = 0.;
        cout << "  vinc " << vinc << endl;
        cout << "  FWHM_val " << gFWHM_val << endl;

        // To estimate how much longer the program will take.
        int nvoxels_to_go_across = sizeSlice * sizePhase * sizeRead;
        int running_index = 0;
        int pref_ratio = 0;

        /////////////////////
        /// Smoothing loop //
        /////////////////////
        // cout << "  DEBUG " << dist(1., 1., 1., 1., 2., 1., dX, dY, dZ) << endl;
        cout << "  Smoothing with Gaussian..." << endl;
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    running_index++;
                    if ((running_index * 100) / nvoxels_to_go_across != pref_ratio) {
                        cout << "\r  Progress: %" << (running_index * 100) / nvoxels_to_go_across << flush;
                        pref_ratio = (running_index * 100) / nvoxels_to_go_across;
                    }
                    for (int it = 0; it < nrep; ++it) {
                        *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = 0;
                        // *(smoothed_data + nxy * iz + nx * ix + iy) = 0;
                        for (int it_i = max(0, it - vinc); it_i < min(it + vinc + 1, nrep); ++it_i) {
                            if (*(nim_inputf_data + nxyz * it_i + nxy * iz + nx * ix + iy) != 0) {
                                dist_i = abs (it-it_i);
                                // cout << "  Debug 4 " << gaus(dist_i , FWHM_val) << endl;
                                // cout << "  Debug 5 " << dist_i << endl;
                                // if (*(nim_input_data + nxy * iz + nx * ix + iy) == 3) cout << "  Debug 4b " << endl;

                                // dummy = *(layer_data + nxy * iz_i + nx * ix_i + iy_i);
                                *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it_i + nxy * iz + nx * ix + iy) * gaus(dist_i, gFWHM_val);
                                *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + gaus(dist_i, gFWHM_val);
                            }
                        }
                    }
                }
            }
        }
        cout << endl;
    }  // Close smoothing loop (for Gaussian smoothing)

    ////////////////////////////////////////////
    // Smoothing loop (for box-car smoothing) //
    ////////////////////////////////////////////
    if (do_box) {
        // float_kernel_size = 10;  // Corresponds to one voxel size.
        int vinc = bFWHM_val;  // Ignore if voxel is too far.
        cout << "  vinc " << vinc << endl;
        float dist_i = 0.;
        // Estimate how long this program will take.
        int nvoxels_to_go_across = sizeSlice * sizePhase * sizeRead;
        int running_index = 0;
        int pref_ratio = 0;

        ////////////////////
        // Smoothing loop //
        ////////////////////
        // cout << " DEBUG " << dist(1., 1., 1., 1., 2., 1., dX, dY, dZ) << endl;
        cout << "  Smoothing with box-car function." << endl;
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    running_index++;
                    if ((running_index * 100) / nvoxels_to_go_across != pref_ratio) {
                        cout << "\r  Progress: %" << (running_index * 100) / nvoxels_to_go_across << flush;
                        pref_ratio = (running_index * 100) / nvoxels_to_go_across;
                    }
                    for (int it = 0; it < nrep; ++it) {
                        *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = 0;
                        // *(smoothed_data + nxy * iz + nx * ix + iy) = 0;
                        for (int it_i=max(0, it - vinc); it_i < min(it + vinc + 1, nrep); ++it_i) {
                            if (*(nim_inputf_data + nxyz * it_i + nxy * iz + nx * ix + iy) != 0) {
                                // No distance here. Just averaging.
                                // dist_i = abs (it-it_i);
                                *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) + *(nim_inputf_data + nxyz * it_i + nxy * iz + nx * ix + iy);
                                *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) = *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) + 1;
                            }
                        }
                    }
                }
            }
        }
        cout << endl;
    }  // Closed smoothing loop (for box-car smoothing)

    ///////////////////////////////
    // Correcting for edge error //
    ///////////////////////////////
    for (int it = 0; it < nrep; ++it) {
        for (int iz = 0; iz < sizeSlice; ++iz) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead-0; ++ix) {
                    if (*(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy) != 0) {
                        *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) / *(gausweight_data + nxyz * it + nxy * iz + nx * ix + iy);
                    }
                }
            }
        }
    }
    // cout << "  Running also until here 5..." << endl;
    // cout << "  Slope " << smoothed->scl_slope << " " << nim_inputfi->scl_slope << endl;
    smoothed->scl_slope = nim_inputfi->scl_slope;

    if (nim_inputfi->scl_inter != 0) {
        cout << "  ########################################## " << endl;
        cout << "  #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << "  ## The NIFTI scale factor is asymmetric ## " << endl;
        cout << "  ## Why would you do such a thing????    ## " << endl;
        cout << "  #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << "  ########################################## " << endl;
    }

    string prefix = "smoothed_";
    string filename = (string) (finfi);
    string outfilename = prefix+filename;
    cout << "  Writing as = " << outfilename.c_str() << endl;

    const char *fout_1 = outfilename.c_str();
    if (nifti_set_filenames(smoothed, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(smoothed);

    // const char *fout_2 = "debug.nii";
    // if (nifti_set_filenames(debug, fout_2 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(debug);

    cout << "  Finished." << endl;
    return 0;
}

// float dist(float x1, float y1, float z1, float x2, float y2, float z2,
//            float dX, float dY, float dZ) {
// return sqrt((x1 - x2) * (x1 - x2) * dX * dX + (y1 - y2) * (y1 - y2) * dY * dY + (z1 - z2) * (z1 - z2) * dZ * dZ);
// }

float gaus(float distance, float sigma) {
    return 1. / (sigma * sqrt(2. * 3.141592)) * exp(-0.5 * distance * distance / (sigma * sigma));
}
