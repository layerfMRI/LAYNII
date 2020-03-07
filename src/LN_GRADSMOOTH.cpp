

#include "../dep/laynii_lib.h"

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
    "    -input       : Nifti (.nii) that should be smooth. It should have\n"
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
    "    - If you run this on EPI-T1 data consider making them pretty, E.g: \n"
    "        start_bias_field.sh T1.nii \n"
    "        denoise_me.sh bico_T1.nii \n"
    "        short_me.sh denoised_bico_T1.nii \n"
    "        smooth_me.sh denoised_bico_T1.nii 0.5 \n"
    "        mv smooth_denoised_bico_T1.nii new_T1.nii \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char * fin_2 = NULL, * fin_1 = NULL, * fin_3 = NULL;
    int ac, twodim = 0, do_masking = 0, within = 0, acros = 0;
    float FWHM_val = 0, selectivity = 0.1;
    if (argc < 3) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-gradfile")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -gradfile\n");
                return 1;
            }
            fin_2 = argv[ac];
        } else if (!strcmp(argv[ac], "-FWHM")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin_1 = argv[ac];
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
            fin_3 = argv[ac];
            do_masking = 1;
            cout << "Set voxels to 0 outside layers (mask option)." << endl;
        } else if (!strcmp(argv[ac], "-selectivity")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -selectivity\n");
                return 1;
            }
            selectivity = atof(argv[ac]);
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-gradfile'\n");
        return 1;
    }
    if (acros + within !=1) {
        cout << "Please choose within or across." << endl;
        return 2;
    }
    if (acros ==1) {
        cout << "Smoothing across gradients is not implemented yet. Use -within instead  " << endl;
        return 2;
    }

    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", fin_1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", fin_2);
        return 2;
    }

    log_welcome("LN_GRADSMOOTH");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    int sizeSlice = nii2->nz;
    int sizePhase = nii2->nx;
    int sizeRead = nii2->ny;
    int nrep = nii1->nt;
    int nx = nii2->nx;
    int nxy = nii2->nx * nii2->ny;
    int nxyz = nii2->nx * nii2->ny * nii2->nz;
    float dX = nii2->pixdim[1];
    float dY = nii2->pixdim[2];
    float dZ = nii2->pixdim[3];
    // If you are running the smoothing in 2D, it will still go thought the
    // entire pipeline the only difference is that the weights in a certain
    // direction are suppressed doing it in 2D, will not speed up the program
    if  (twodim == 1) {
        dZ = 1000 * dZ;
    }

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* nii_mask = copy_nifti_as_float32(nii2);
    float* nii_mask_data = static_cast<float*>(nii_mask->data);

    // Allocate new nifti images
    nifti_image* smooth = copy_nifti_as_float32(nii_input);
    float* smooth_data = static_cast<float*>(smooth->data);
    nifti_image* gaussw = copy_nifti_as_float32(nii_input);
    // NOTE(Renzo): Gauss weight is a gemoetry factor and only needs to be
    // estimated once. So with the next lines I am saving RAM.
    gaussw->nt = 1;
    gaussw->nvox = gaussw->nvox / nrep;
    float* gaussw_data = static_cast<float*>(gaussw->data);

    // ------------------------------------------------------------------------
    // Mask nifti related part
    nifti_image* nii_roi;
    float* nii_roi_data;
    if (do_masking == 1) {
        nifti_image* nii3 = nifti_image_read(fin_3, 1);
        if (!nii3) {
            fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", fin_3);
            return 2;
        }
        log_nifti_descriptives(nii3);
        nii_roi = copy_nifti_as_float32(nii3);
        nii_roi_data = static_cast<float*>(nii_roi->data);
    }
    // ========================================================================

    cout << "  Time dimension of smooth. Output file: " << smooth->nt << endl;

    //  float kernel_size = 10;  // corresponds to one voxel sice.
    int vinc = max(1., 2. * FWHM_val/dX);  // Ignore if voxel is too far
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
    int nvoxels_to_go_across = sizeSlice * sizePhase * sizeRead;
    int running_index = 0;
    int pref_ratio = 0;

    if (sizeSlice * sizePhase * sizeRead > 32767) {
        cout << "  The number of voxels is bigger than the range of int the time estimation will be wrong." << endl;
    }

    if (do_masking == 1) {
        nvoxels_to_go_across = 0;
        for (int islice = 0; islice < sizeSlice; ++islice) {
            for (int iy = 0; iy < sizePhase; ++iy) {
                for (int ix = 0; ix < sizeRead; ++ix) {
                    if (*(nii_roi_data + nxy * islice + nx * ix + iy) > 0) {
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
                if (!(!(*(nii_roi_data + nxy * iz + nx * ix + iy) > 0) && (do_masking == 1))) {
                    // if (iz==sizeSlice/2 && iy == sizePhase/2-4 && ix == sizeRead/2-4) { // debug loop open

                    // This is to write out how many more voxels I have to go through.
                    running_index++;
                    if ((running_index * 100)/nvoxels_to_go_across != pref_ratio) {
                        cout << "\r  Progress: " << (running_index * 100)/nvoxels_to_go_across << "%" << flush;
                        pref_ratio = (running_index * 100)/nvoxels_to_go_across;
                    }
                    // I am cooking in a clean kitchen.
                    // This might not be neccessary, just to be on the safe side
                    *(gaussw_data + nxy * iz + nx * ix + iy) = 0;
                    *(smooth_data + nxy * iz + nx * ix + iy) = 0;
                    NvoxInVinc = 0;
                    local_val = *(nii_mask_data + nxy * iz + nx * ix + iy);

                    // Examining the environment and determining what the
                    // signal intensities are and what its distribution are
                    for (int iz_i = max(0, iz - vinc); iz_i <= min(iz + vinc, sizeSlice - 1); ++iz_i) {
                        for (int iy_i = max(0, iy - vinc); iy_i <= min(iy + vinc, sizePhase - 1); ++iy_i) {
                            for (int ix_i = max(0, ix - vinc); ix_i <= min(ix + vinc, sizeRead - 1); ++ix_i) {
                                vec1[NvoxInVinc] = (double) *(nii_mask_data + nxy * iz_i + nx * ix_i + iy_i);
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
                                value_dist = fabs(local_val - *(nii_mask_data + nxy * iz_i + nx * ix_i + iy_i));

                                // * (debug_data + nxy * iz_i + nx * ix_i + iy_i) = gaus(dist_i , FWHM_val)/gaus(0, FWHM_val)
                                //                                               * gaus(value_dist, grad_stdev * 0.1) /gaus(0, grad_stdev * 0.1) ;

                                temp_wight_factor = gaus(dist_i, FWHM_val) * gaus(value_dist, grad_stdev * selectivity)/gaus(0, grad_stdev * selectivity);

                                // Gauss data are important to avoid local
                                // scaling differences. When the kernel size
                                // changes (e.g. at edge of images).
                                // This is a geometric parameter and only need
                                // to be calculated for one time point. This
                                // might be avoidable, if the Gauss fucnction
                                // is better normalized.
                                *(gaussw_data + nxy * iz + nx * ix + iy) = *(gaussw_data + nxy * iz + nx * ix + iy)
                                                                               + temp_wight_factor;

                                for (int it = 0; it < nrep; ++it) {  // loop across all time steps
                                    *(smooth_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smooth_data + nxyz * it + nxy * iz + nx * ix + iy )
                                                                                             + *(nii_input_data + nxyz * it + nxy * iz_i + nx * ix_i + iy_i) * temp_wight_factor;
                                }
                            }
                        }
                    }
                    // Scaling signal intensity with the overall Gauss leakage
                    if (*(gaussw_data + nxy * iz + nx * ix + iy) > 0) {
                        for (int it = 0; it < nrep; ++it) {
                            *(smooth_data + nxyz * it + nxy * iz + nx * ix + iy) = *(smooth_data + nxyz * it + nxy * iz + nx * ix + iy) / *(gaussw_data + nxy * iz + nx * ix + iy);
                        }
                    }
                    // if (* (nii_mask_data + nxy * iz + nx * ix + iy) < = 0) {
                    // * (smooth_data + nxy * iz + nx * ix + iy) = *(nii_input_data + nxy * iz + nx * ix + iy) ;
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
    //                     if (*(nii_mask_data + nxy * islice + nx * ix + iy) == 0) {
    //                     *(smooth_data + nxyz * it + nxy * islice + nx * ix + iy) = 0;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // cout << "  Runing also until here 5... " << endl;
    // cout << "  Slope " << smooth->scl_slope << " " << nii1->scl_slope << endl;

    smooth->scl_slope = nii1->scl_slope;

    if (nii1->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " ## Why would you do such a thing????    ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    save_output_nifti(fin_1, "smooth", smooth, true);

    cout << "  Finished." << endl;
    return 0;
}
