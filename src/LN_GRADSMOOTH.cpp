
#include "./common.h"
#include "./renzo_stat.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_GRADSMOOTH : Layering algorithm based on iterative smoothing.\n"
    "\n"
    "    This program smooths data within layer or columns. In order to \n"
    "    avoid smoothing across masks, a crawler smooths only across \n"
    "    connected voxels.\n"
    "\n"
    "Usage:\n"
    "    LN_GRADSMOOTH -gradfile gradfile.nii -input activity_map.nii -FWHM 1 -within -selectivity 0.1 \n"
    "\n"
    "Options:\n"
    "    -help        : Show this help.\n"
    "    -gradfile    : Nifti (.nii) that is used to estimate local gradients\n"
    "                   only the first time point of this file is used. It \n"
    "                   should have the same spatial dimensions as the input.\n"
    "    -input       : Nifti (.nii) that should be smoothed. It should have\n"
    "                   same dimensions as layer file.\n"
    "    -FWHM        : Amount of smoothing in mm.\n"
    "    -twodim      : (Optional) Smooth in 2 dimensions only. \n"
    "    -mask        : (Optional) Nifti (.nii) that is used mask activity \n"
    "                   outside. This option can speed up processing.\n"
    "    -within      : (Optional) Determines that smoothing should happen \n"
    "                   within similar values, not across different values.\n"
    "    -acros       : (Optional) Determines that smoothing should happen \n"
    "                   across different values, not within similar values.\n"
    "                   NOTE: This option is not working yet.\n"
    "    -selectivity : (Optional) Makes the smoothing more or less specific \n"
    "                   to a certain gradient range. 0.05 is only within \n"
    "                   very similar values. 0.9 is almost independent of \n"
    "                   the gradient file 0.1 is default.\n"
    "\n"
    "Notes:\n"
    "     - If you run this on EPI-T1 data consider making them pretty, E.g:  \n"
    "       start_bias_field.sh T1.nii  \n"
    "       denoise_me.sh bico_T1.nii \n"
    "       short_me.sh denoised_bico_T1.nii \n"
    "       smooth_me.sh denoised_bico_T1.nii 0.5 \n"
    "       mv smoothed_denoised_bico_T1.nii new_T1.nii \n"
    "     - This program now supports INT16, INT32 and FLOAT32.\n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * fmaski = NULL,  * fout = NULL,  * finfi = NULL,  * froii = NULL;
    int ac, twodim = 0, do_masking = 0, within = 0, acros = 0;
    float FWHM_val = 0, selectivity = 0.1;
    if (argc < 3) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-gradfile")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -gradfile\n");
                return 1;
            }
            fmaski = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            finfi = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-twodim")) {
            twodim = 1;
            cout << "Smooth only in 2D." << endl;
        } else if (!strcmp(argv[ac], "-within")) {
            within = 1;
            cout << "Within similar values." << endl;
        } else if (!strcmp(argv[ac], "-acros")) {
            acros = 1;
            cout << "Across different values." << endl;
        } else if (!strcmp(argv[ac], "-mask")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -mask\n");
                return 1;
            }
            froii = argv[ac];  // Assign pointer, no string copy
            do_masking = 1;
            cout << "Set every thing to zero outside the layers (masking option)." << endl;
        } else if (!strcmp(argv[ac], "-selectivity")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -selectivity\n");
                return 1;
            }
            selectivity = atof(argv[ac]);  // Assign pointer, no string copy
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
    nifti_image * nim_inputfi = nifti_image_read(finfi, 1);
    if (!nim_inputfi) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", finfi);
        return 2;
    }
    if (!fmaski) {
        fprintf(stderr, "** missing option '-gradfile'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image * nim_maski = nifti_image_read(fmaski, 1);
    if (!nim_maski) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", fmaski);
        return 2;
    }
    if (acros + within !=1) {
        cout << "Please choose within or across." << endl;
        return 2;
    }
    if (acros ==1) {
        cout << "Smoothing across gradients is not implemented yet. Use -within instead  " << endl;
        return 2;
    }

    // Get dimensions of input
    int sizeSlice = nim_maski->nz;
    int sizePhase = nim_maski->nx;
    int sizeRead = nim_maski->ny;
    int nrep = nim_inputfi->nt;
    int nx = nim_maski->nx;
    int nxy = nim_maski->nx  * nim_maski->ny;
    int nxyz = nim_maski->nx  * nim_maski->ny  * nim_maski->nz;
    float dX = nim_maski->pixdim[1];
    float dY = nim_maski->pixdim[2];
    float dZ = nim_maski->pixdim[3];


    // If you are running the smoothing in 2D, it will still go thought the
    // entire pipeline the only difference is that the weights in a certain
    // direction are suppressed doing it in 2D, will not speed up the program
    if  (twodim == 1) {
        dZ = 1000  * dZ;
    }

    nifti_image * nim_inputf = nifti_copy_nim_info(nim_inputfi);
    nim_inputf->datatype = NIFTI_TYPE_FLOAT32;
    nim_inputf->nbyper = sizeof(float);
    nim_inputf->data = calloc(nim_inputf->nvox, nim_inputf->nbyper);
    float * nim_inputf_data = (float * ) nim_inputf->data;

    nifti_image * nim_mask = nifti_copy_nim_info(nim_maski);
    nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    nim_mask->nbyper = sizeof(float);
    nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);
    float * nim_mask_data = (float * ) nim_mask->data;

    nifti_image * nim_roi = nifti_copy_nim_info(nim_maski);
    nim_roi->datatype = NIFTI_TYPE_FLOAT32;
    nim_roi->nbyper = sizeof(float);
    nim_roi->data = calloc(nim_roi->nvox, nim_roi->nbyper);
    float * nim_roi_data = (float * ) nim_roi->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    // here, I am loading them in their native datatype //////////
    // and translate them to the datatime I like best ////////////
    //////////////////////////////////////////////////////////////
    if (nim_inputfi->datatype == NIFTI_TYPE_FLOAT32 ||  nim_inputfi->datatype == NIFTI_TYPE_INT32) {
        float * nim_inputfi_data = (float * ) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz  * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz  * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_inputfi->datatype == NIFTI_TYPE_INT16 || nim_inputfi->datatype == DT_UINT16) {
        short  * nim_inputfi_data = (short  * ) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz  * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz  * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_inputfi->datatype == DT_FLOAT64 || nim_inputfi->datatype == NIFTI_TYPE_FLOAT64) {
        double  * nim_inputfi_data = (double  * ) nim_inputfi->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_inputf_data + nxyz  * it + nxy * islice + nx * ix + iy) = (float) (*(nim_inputfi_data + nxyz  * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_maski->datatype == NIFTI_TYPE_FLOAT32 || nim_maski->datatype == NIFTI_TYPE_INT32) {
        float * nim_maski_data = (float * ) nim_maski->data;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_mask_data + nxy * islice + nx * ix + iy) = (float) (*(nim_maski_data + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (nim_maski->datatype == NIFTI_TYPE_INT16 || nim_maski->datatype == DT_UINT16) {
        short  * nim_maski_data = (short  * ) nim_maski->data;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_mask_data + nxy * islice + nx * ix + iy) = (float) (*(nim_maski_data + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (nim_maski->datatype == DT_FLOAT64 || nim_maski->datatype == NIFTI_TYPE_FLOAT64) {
        double  * nim_maski_data = (double  * ) nim_maski->data;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    *(nim_mask_data + nxy * islice + nx * ix + iy) = (float) (*(nim_maski_data + nxy * islice + nx * ix + iy));
                }
            }
        }
    }
    if (do_masking == 1) {
        // Read input dataset, including data
        nifti_image * nim_roii = nifti_image_read(froii, 1);
        if (!nim_roii) {
            fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", froii);
            return 2;
        }
        if (nim_roii->datatype == NIFTI_TYPE_FLOAT32 ||  nim_roii->datatype == NIFTI_TYPE_INT32) {
            float * nim_roii_data = (float * ) nim_roii->data;
            for (int it = 0; it < nrep; ++it) {
                for (int islice = 0; islice < sizeSlice; ++islice) {
                    for (int iy = 0; iy < sizePhase; ++iy) {
                        for (int ix = 0; ix < sizeRead; ++ix) {
                            *(nim_roi_data + nxyz  * it + nxy * islice + nx * ix + iy) = (float) (*(nim_roii_data + nxyz  * it + nxy * islice + nx * ix + iy));
                        }
                    }
                }
            }
        }
        if (nim_roii->datatype == NIFTI_TYPE_INT16 || nim_roii->datatype == DT_UINT16) {
            short  * nim_roii_data = (short  * ) nim_roii->data;
            for (int it = 0; it < nrep; ++it) {
                for (int islice = 0; islice < sizeSlice; ++islice) {
                    for (int iy = 0; iy < sizePhase; ++iy) {
                        for (int ix = 0; ix < sizeRead; ++ix) {
                            *(nim_roi_data + nxyz  * it + nxy * islice + nx * ix + iy) = (float) (*(nim_roii_data + nxyz  * it + nxy * islice + nx * ix + iy));
                        }
                    }
                }
            }
        }
    }
    // Write out some stuff that might be good to know, if you want to debug
    cout << sizeSlice << " Slices | " << sizePhase << " Phase_steps | " << sizeRead << " Read_steps | " << nrep << " Time_steps " << endl;
    cout << "  Voxel size = " << dX << " x " << dY << " x " << dZ << endl;
    cout << "  Datatype 1 = " << nim_inputf->datatype << endl;
    cout << "  Datatype 2 = " << nim_mask->datatype << endl;

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

    // Gauss weight is a gemoetry factor and only needs to be estimated once
    // So with the next lines I am saving RAM
    gausweight->nt = 1;
    gausweight->nvox = gausweight->nvox / nrep;
    float * smoothed_data = (float * ) smoothed->data;
    float * gausweight_data = (float * ) gausweight->data;

    // nifti_image * debug = nifti_copy_nim_info(nim_inputf);
    // debug->datatype = NIFTI_TYPE_FLOAT32;
    // debug->nbyper = sizeof(float);
    // debug->data = calloc(debug->nvox, debug->nbyper);
    // float * debug_data = (float * ) debug->data;

    // if (do_masking == 1) {
    //     for (int islice = 0; islice < sizeSlice; ++islice) {
    //         for (int iy = 0; iy < sizePhase; ++iy) {
    //             for (int ix = 0; ix < sizeRead; ++ix) {
    //                 *(debug_data + nxy * islice + nx * ix + iy) = 0;
    //             }
    //         }
    //     }
    // }

    float dist(float x1, float y1, float z1, float x2, float y2, float z2,
               float dX, float dY, float dZ);
    float gaus(float distance, float sigma);

    cout << "  Time dimension of smoothed. Output file: " << smoothed->nt << endl;

    //  float kernel_size = 10;  // corresponds to one voxel sice.
    int vinc = max(1., 2.  * FWHM_val/dX);  // Ignore if voxel is too far
    float dist_i = 0.;
    cout << "  vinc " << vinc << endl;
    cout << "  FWHM_val " << FWHM_val << endl;

    // To store temp values, so I don't make the same computations repeatedly.
    float temp_wight_factor = 0;

    //////////////////////////////////////////
    // Finding the range of gradient values //
    //////////////////////////////////////////

    // valued that I need to characterize the local signals in the vicinity.
    float local_val = 0;
    int NvoxInVinc = (2 * vinc+1) * (2 * vinc+1) * (2 * vinc+1);
    double vec1[NvoxInVinc];
    for (int it = 0; it < NvoxInVinc; it++) vec1[it] = 0;
    float grad_stdev = 0;
    float value_dist = 0;


    // Estimate and output of program process and how much longer it will take.
    int nvoxels_to_go_across = sizeSlice  * sizePhase  * sizeRead;
    int running_index = 0;
    int pref_ratio = 0;

    if (sizeSlice  * sizePhase  * sizeRead > 32767) {
        cout << "  The number of voxels is bigger than the range of int the time estimation will be wrong." << endl;
    }

    if (do_masking == 1) {
        nvoxels_to_go_across = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    if (*(nim_roi_data + nxy * islice + nx * ix + iy) > 0) {
                        nvoxels_to_go_across = nvoxels_to_go_across +1;
                    }
                }
            }
        }
    }
    cout << "  The number of voxels to go across = " << nvoxels_to_go_across << endl;

    ////////////////////
    // Smoothing loop //
    ////////////////////

    cout << "  Big smoothing loop is being done..." << endl;
    for (int iz = 0; iz < sizeSlice; ++iz) {
        for (int iy = 0; iy < sizePhase; ++iy) {
            for (int ix = 0; ix < sizeRead; ++ix) {
                if (!(!(*(nim_roi_data + nxy * iz + nx * ix + iy) > 0) && (do_masking == 1))) {
                    // if (iz==sizeSlice/2 && iy == sizePhase/2-4 && ix == sizeRead/2-4) { // debug loop open

                    // This is to write out how many more voxels I have to go through.
                    running_index++;
                    if ((running_index * 100)/nvoxels_to_go_across != pref_ratio) {
                        cout << "\r  Progress: " << (running_index * 100)/nvoxels_to_go_across << "%" << flush;
                        pref_ratio = (running_index * 100)/nvoxels_to_go_across;
                    }
                    // I am cooking in a clean kitchen.
                    // This might not be neccessary, just to be on the safe side
                    *(gausweight_data + nxy * iz + nx * ix + iy) = 0;
                    *(smoothed_data + nxy * iz + nx * ix + iy) = 0;
                    NvoxInVinc = 0;
                    local_val = *(nim_mask_data + nxy * iz + nx * ix + iy);

                    // Examining the environment and determining what the
                    // signal intensities are and what its distribution are
                    for (int iz_i = max(0, iz - vinc); iz_i <= min(iz + vinc, sizeSlice - 1); ++iz_i) {
                        for (int iy_i = max(0, iy - vinc); iy_i <= min(iy + vinc, sizePhase - 1); ++iy_i) {
                            for (int ix_i = max(0, ix - vinc); ix_i <= min(ix + vinc, sizeRead - 1); ++ix_i) {
                                vec1[NvoxInVinc] = (double)  *(nim_mask_data + nxy * iz_i + nx * ix_i + iy_i);
                                NvoxInVinc++;
                            }
                        }
                    }

                    // Standard deviation of the signal valued in the vicinity.
                    // This is necessary to normalize how many voxels are
                    // contributing to the local smoothing.
                    // grad_stdev = (float)  gsl_stats_sd (vec1, 1, NvoxInVinc);
                    grad_stdev = (float) ren_stdev (vec1, NvoxInVinc);
                    for (int iz_i = max(0, iz-vinc); iz_i <= min(iz+vinc, sizeSlice-1); ++iz_i) {
                        for (int iy_i = max(0, iy-vinc); iy_i <= min(iy+vinc, sizePhase-1); ++iy_i) {
                            for (int ix_i = max(0, ix-vinc); ix_i <= min(ix+vinc, sizeRead-1); ++ix_i) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i, dX, dY, dZ);
                                value_dist = fabs(local_val - *(nim_mask_data + nxy * iz_i + nx * ix_i + iy_i) );

                                // * (debug_data + nxy * iz_i + nx * ix_i + iy_i) = gaus(dist_i , FWHM_val)/gaus(0, FWHM_val)
                                //                                                   * gaus(value_dist, grad_stdev * 0.1) /gaus(0, grad_stdev * 0.1) ;

                                temp_wight_factor = gaus(dist_i, FWHM_val)  * gaus(value_dist, grad_stdev * selectivity)/gaus(0, grad_stdev * selectivity);

                                // Gauss data are important to avoid local
                                // scaling differences. When the kernel size
                                // changes (e.g. at edge of images).
                                // This is a geometric parameter and only need
                                // to be calculated for one time point. This
                                // might be avoidable, if the Gauss fucnction
                                // is better normalized.
                                *(gausweight_data + nxy * iz + nx * ix + iy) = *(gausweight_data + nxy * iz + nx * ix + iy)
                                                                               + temp_wight_factor;

                                for (int it = 0; it < nrep; ++it) {  // loop across all time steps
                                    *(smoothed_data + nxyz  * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz  * it + nxy * iz + nx * ix + iy  )
                                                                                              + *(nim_inputf_data + nxyz  * it + nxy * iz_i + nx * ix_i + iy_i)  * temp_wight_factor;
                                }
                            }
                        }
                    }
                    // Scaling signal intensity with the overall Gauss leakage
                    if (*(gausweight_data + nxy * iz + nx * ix + iy) > 0) {
                        for (int it = 0; it < nrep; ++it) {
                            *(smoothed_data + nxyz  * it + nxy * iz + nx * ix + iy) = *(smoothed_data + nxyz  * it + nxy * iz + nx * ix + iy)/  *(gausweight_data + nxy * iz + nx * ix + iy);
                        }
                    }
                    // if (* (nim_mask_data + nxy * iz + nx * ix + iy) < = 0) {
                    //     * (smoothed_data + nxy * iz + nx * ix + iy) = *(nim_inputf_data + nxy * iz + nx * ix + iy) ;
                    // }
                    // }  //debug loop closed
                }
            }
        }
    }
    cout << endl;
    // Note(Renzo): I am not sure if there is a case there masking makes sense?
    // I just leave it in.
    // if (do_masking == 1) {
    //     for (int it = 0; it < nrep; ++it) {
    //         for (int islice = 0; islice < sizeSlice; ++islice) {
    //             for (int iy = 0; iy < sizePhase; ++iy) {
    //                 for (int ix = 0; ix < sizeRead; ++ix) {
    //                     if (*(nim_mask_data + nxy * islice + nx * ix + iy) == 0) {
    //                         *(smoothed_data + nxyz  * it + nxy * islice + nx * ix + iy) = 0;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // cout << "  Runing also until here 5... " << endl;
    // cout << "  Slope " << smoothed->scl_slope << " " << nim_inputfi->scl_slope << endl;

    smoothed->scl_slope = nim_inputfi->scl_slope;

    if (nim_inputfi->scl_inter != 0) {
        cout << " ########################################### " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ###### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric  ## " << endl;
        cout << " ## Why would you do such a thing????     ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ###### " << endl;
        cout << " ########################################### " << endl;
    }

    string prefix = "smoothed_";
    string filename = (string) (finfi);
    string outfilename = prefix+filename;

    cout << "  Writing as = " << outfilename.c_str() << endl;

    const char * fout_1 = outfilename.c_str();
    if (nifti_set_filenames(smoothed, fout_1, 1, 1)) {
        return 1;
    }
    nifti_image_write(smoothed);

    // const char * fout_2="debug.nii" ;
    // if (nifti_set_filenames(debug, fout_2 , 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(debug);

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
    return (1./(sigma * sqrt(2. * 3.141592))
            * exp(-0.5 * distance * distance/(sigma * sigma)));
}
