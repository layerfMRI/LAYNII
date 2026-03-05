
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_FRISGO : This program corrects/mitigates the Fuzzy Ripple artifacts\n"
    "             in dual polarity 3D-EPI fMRI data.\n"
    "              \n"
    "\n"
    "Usage:\n"
    "    LN2_FRISGO -input timeseries.nii -simple\n"
    "    LN2_FRISGO -input timeseries.nii -spline\n"
    "    LN2_FRISGO -input timeseries.nii -lpass 1.0 \n"
    "    LN2_FRISGO -input timeseries.nii -box 1 \n"
    "    ../LN2_FRISGO -input lo_BOLD_intemp.nii -spline \n" 
    "    ../LN2_FRISGO -input lo_BOLD_intemp.nii -lpass 1 \n" 
    "    ../LN2_FRISGO -input lo_BOLD_intemp.nii -box 1 \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii) file with time series data. \n"
    "    -simple : Weighed average with 3 time points [0.25*t-1 + 0.5*t0 + 0.25*t+1].\n"
    "              Particularly useful for despiked data.\n"
    "    -spline : Estimates the fMRI signal right between TRs. Similar to slice scan \n"
    "              timing correction, but for 3D-EPI. This way the trigger timing is\n"
    "              directly aligned with the k-space center. Here we use a third order\n"
    "              interpolation across time to minimize temporal blurring.\n"
    "    -lpass  : Estimates slow changes of the fuzzy ripples across time. A Gaussian\n"
    "              sliding window is used to estimate the fuzzy ripples. Specify the value\n"
    "              of the Gaussian size in units of TR (float precision). Make sure to use\n"
    "              values bigger than your trial timing. E.g. '40'.\n"
    "    -runwise: Assumes that the fuzzy ripples are identical for the entire run.\n"
    // "    -calib  : (WIP) This option assumes that the first 4 volumes have alternating \n"
    // "              read polarity, while the remainder of the run has one polarity.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "    -output : (Optional) Output filename, including .nii or .nii.gz, and path if needed.\n"
    "              Overwrites existing files.\n"    
    "\n"
    "Notes:\n"
    "    This program was featured in the following ISMRM 2025 presentation:\n"
    "      <https://youtu.be/a-X1JW7-Pk4?si=QRG_wgqPiO74Zphe>\n"
    "\n"
    "    More on fuzzy ripples: \n"
    "      - Huber, L., Stirnberg, R., Morgan, A.T., Feinberg, D.A., Ehses, P., Knudsen, L., \n"
    "        Gulban, O.F., Koiso, K., Gephart, I., Swegle, S., Wardle, S.G., Persichetti, A.S.,\n"
    "        Beckett, A.J.S., Stöcker, T., Boulant, N., Poser, B.A., Bandettini, P.A., 2025.\n"
    "        Short‐term gradient imperfections in high‐resolution EPI lead to Fuzzy Ripple artifacts.\n"
    "        Magnetic Resonance in Med. <https://doi.org/10.1002/mrm.30489>\n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false, mode_debug = false;
    bool mode_simple=false, mode_spline=false, mode_lpass=false, mode_runwise=false;
    bool mode_oldname=false;
    char *fin = NULL, *fout = NULL;
    int ac;

    float gFWHM_val;

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
            fin = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-simple")) {
            mode_simple = true;

        // --------------------------------------------------------------------
        } else if (!strcmp(argv[ac], "-spline")) {
            mode_spline = true;
        } else if (!strcmp(argv[ac], "-tshift")) {  // Old convention 
            mode_spline = true;
            mode_oldname = true;

        } else if (!strcmp(argv[ac], "-runwise")) {
            mode_runwise = true;
        } else if (!strcmp(argv[ac], "-run")) {  // Old convention
            mode_runwise = true;
            mode_oldname = true;

        // --------------------------------------------------------------------
        } else if (!strcmp(argv[ac], "-lpass")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -lpass\n");
                return 1;
            } 
            gFWHM_val = atof(argv[ac]);
            mode_lpass = true;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    if (mode_simple==false && mode_spline==false && mode_lpass==false && mode_runwise==false) {
        fprintf(stderr, "** select at least one correction method.\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN2_FRISGO");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    const uint64_t size_x = nii_input->nx;
    const uint64_t size_y = nii_input->ny;
    const uint64_t size_z = nii_input->nz;
    const uint64_t size_time = nii_input->nt;

    const uint64_t end_time = size_time - 1;

    const uint64_t nxyz = size_x * size_y * size_z;
    const uint64_t nxyzt = size_x * size_y * size_z * size_time;

    const float dT = 1.0;

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocating necessary files
    nifti_image* nii_smooth = copy_nifti_as_float32(nii);
    float* nii_smooth_data = static_cast<float*>(nii_smooth->data);

    nifti_image* nii_weight = nifti_copy_nim_info(nii);
    nii_weight->nt = 1;
    nii_weight->datatype = NIFTI_TYPE_FLOAT32;
    nii_weight->nbyper = sizeof(float);
    nii_weight->nvox = size_x * size_y * size_z;
    nii_weight->data = calloc(nii_weight->nvox, nii_weight->nbyper);
    float* nii_weight_data = static_cast<float*>(nii_weight->data);

    // ========================================================================
    // Fuzzy ripple correction using weighted averaging of a symmetric window
    // ========================================================================
    // Set to zero
    for (int t = 0; t < size_time; ++t) {
        for (int i = 0; i < nxyz; ++i) {
            *(nii_smooth_data + t*nxyz + i) = 0.0;
        }
    }

    if (mode_simple) {
        // NOTE[Faruk]: I think this is a better combination for despiked data
        // where the spike volume is imputed with time neighbors. For example
        // 3rd order spline with 4 neighbors gives asymmetric smoothing, as one
        // of the 4 volumes is made of other 2.
        cout << "  Doing simple correction..." << endl;

        uint64_t ix, iy, iz, it;
        for (uint64_t i = 0; i < nxyzt; ++i) {
            std::tie(ix, iy, iz, it) = ind2sub_4D_64(i, size_x, size_y, size_z);

            // Weighted average of neighboring time points
            if (it > 0 && it < end_time) {
                uint64_t j = sub2ind_4D_64(ix, iy, iz, it-1, size_x, size_y, size_z);
                uint64_t k = sub2ind_4D_64(ix, iy, iz, it+1, size_x, size_y, size_z);
                *(nii_smooth_data + i) = *(nii_data + i) * 0.5 + ( (*(nii_data + j) + *(nii_data + k) ) * 0.25);
            } else if (it == end_time) {
                uint64_t j = sub2ind_4D_64(ix, iy, iz, it-1, size_x, size_y, size_z);
                *(nii_smooth_data + i) = ( *(nii_data + i) + *(nii_data + j) ) * 0.5;
            } else if (it == 0) {
                uint64_t k = sub2ind_4D_64(ix, iy, iz, it+1, size_x, size_y, size_z);
                *(nii_smooth_data + i) = ( *(nii_data + i) + *(nii_data + k) ) * 0.5;
            }
        }
        std::cout << "    Saving..." << std::endl;
        save_output_nifti(fout, "frisgo-simple", nii_smooth, false);
    }

    // ========================================================================
    // Fuzzy ripple using 4 TRs with 3rd order spline
    // ========================================================================    
    if (mode_spline) {
        cout << "  Doing spline correction..." << endl;

        // NOTE[Renzo]: Convention is to fist forward. See assumed time table below
        //   t is trigger
        //   k time point of largest power of fmri signal: k-space center
        // input data are like this:  
        // t       t       t       t       t       t       t       ...
        //     k1      k2      k3      k4      k5      k6      k7  ...
        //
        // output data are like this: 
        // t       t       t       t       t       t       t       ...
        // k1      k1-2    k2-3    k3-4    k4-5    k5-6    k6-7    ...

        if (size_time < 4) {
            std::cout << "  ERROR: Less than 4 TRs, 3rd order spline fitting is not possible." << std::endl;
            return 1; 
        }
        
        // Deal with first two time points, interpolating back
        for (uint64_t i = 0; i < nxyz; ++i) {
            uint64_t t0 = nxyz*0, t1 = nxyz*1, t2 = nxyz*2;
            *(nii_smooth_data + t0 + i) = 0.5 * ( *(nii_data + t0 + i) + *(nii_data + t1 + i) ) ;
            *(nii_smooth_data + t1 + i) = 0.5 * ( *(nii_data + t1 + i) + *(nii_data + t2 + i) ) ;
        }

        // Deal with last two time points copying averages
        for (uint64_t i = 0; i < nxyz; ++i) {
            uint64_t t0 = nxyz*(size_time - 1), t1 = nxyz*(size_time - 2), t2 = nxyz*(size_time - 3);
           *(nii_smooth_data + t0 + i) = 0.5 * ( *(nii_data + t0 + i) + *(nii_data + t1 + i) );
           *(nii_smooth_data + t1 + i) = 0.5 * ( *(nii_data + t1 + i) + *(nii_data + t2 + i) );
        }

        // Deal with all the time points in between
        // 
        // NOTE[Renzo]: y = a + b x + c x^2 + d x^3, for 4 data points, 
        //     forward matrix is y = M A, 
        //     with y = (y1, y2, y3, y4) the signal strengths across 4 time point
        //     with A = a, b, c, d, the prefactors of the cubic polinomial 
        //     y (X=0) is desired, so the other x values are -1.5, -0.5, 0.5, 1.5
        //     so the forward matrix is M = (1, -1.5, 2.25, -3.375; 
        //                                   1, -0.5, 0.25, -0.125; 
        //                                   1, 0.5,  0.25,  0.125; 
        //                                   1, 1.5,  2.25,  3.375)
        //     so the inverse matrix is M^-1 = (- 0.0625, + 0.5625, + 0.5625, - 0.0625 )... (we only need th first row).
        //     and the output is y = M^-1 y, so we can just multiply the input data with the inverse matrix.
        for (uint64_t t = 2; t != size_time-2; ++t) {
            for (uint64_t i = 0; i < nxyz; ++i) {
                *(nii_smooth_data + nxyz*t + i) = 
                    - 0.0625 * ( *(nii_data + nxyz*(t-2) + i))
                    + 0.5625 * ( *(nii_data + nxyz*(t-1) + i))
                    + 0.5625 * ( *(nii_data + nxyz*(t)   + i))
                    - 0.0625 * ( *(nii_data + nxyz*(t+1) + i));
            }
        }
        std::cout << "    Saving..." << std::endl;
        if (mode_oldname) {
            save_output_nifti(fout, "Tshift", nii_smooth, false);  // Old convention 
        } else {
            save_output_nifti(fout, "frisgo-spline", nii_smooth, false);
        }
    }

    // ========================================================================
    // Prepare for low pass or run-wise Fuzzy Ripple correction
    // ========================================================================
    if (mode_lpass || mode_runwise) {
        cout << "  Preparing for low pass or run-wise correction..." << endl;

        // Allocate memory for Gaussian weights
        nifti_image* nii_weight = nifti_copy_nim_info(nii);
        nii_weight->nt = 1;
        nii_weight->datatype = NIFTI_TYPE_FLOAT32;
        nii_weight->nbyper = sizeof(float);
        nii_weight->nvox = size_x * size_y * size_z;
        nii_weight->data = calloc(nii_weight->nvox, nii_weight->nbyper);
        float* nii_weight_data = static_cast<float*>(nii_weight->data);

        // Error time series
        nifti_image* nii_error = copy_nifti_as_float32(nii);
        float* nii_error_data = static_cast<float*>(nii_error->data);

        // Set all error data to zero
        for (uint64_t i = 0; i < nxyzt ; ++i) {
            *(nii_error_data + i) = 0;
        }

        // Calculate error of first time point
        for (uint64_t i = 0; i < nxyz; ++i) {
            *(nii_error_data + i) = *(nii_data + nxyz * 0 + i) 
                                    - *(nii_data + nxyz * 1 + i);
        }

        // Calculate error of last time point
        for (uint64_t i = 0; i < nxyz; ++i) {
            *(nii_error_data + nxyz * (size_time-1) + i) = *(nii_data + nxyz * (size_time-1) + i)
                                                           - *(nii_data + nxyz * (size_time-2) + i);
        }
        
        // Calculate error for the rest of the time points
        for (uint64_t t = 1; t < size_time-1; ++t) {
            for (uint64_t i = 0; i < nxyz; ++i) {
                *(nii_error_data + t * nxyz + i) =  *(nii_data + t * nxyz + i) 
                                                  - 0.5 * ( *(nii_data + (t-1) * nxyz + i))
                                                  - 0.5 * ( *(nii_data + (t+1) * nxyz + i));
            }
        }

        // Inverting error every other TR
            for (uint64_t t = 1; t < size_time; t += 2) {
                for (uint64_t i = 0; i < nxyz; ++i) {
                        *(nii_error_data + t * nxyz + i) *= -1;
                    }
            }

        if (mode_debug) {
            save_output_nifti(fout, "error_maps", nii_error, false);
        }

        // ========================================================================
        // Gaussian smoothing loop
        // ========================================================================
        if (mode_lpass) {
            int vic = static_cast<int>(max(1., 2. * gFWHM_val / dT));  // Ignore if voxel is too far

            std::cout << "  Doing low pass correction (This might take a while)..." << std::endl;
            std::cout << "    FWHM_val = " << gFWHM_val << std::endl;
            std::cout << "    Vicinity = " << vic << std::endl;

            // Allocate memory for smoothed error time series.
            nifti_image* nii_smooth_error = copy_nifti_as_float32(nii);
            float* nii_smooth_error_data = static_cast<float*>(nii_smooth_error->data);

            for (uint64_t i = 0; i < nxyz; ++i) {
                // Progress
                uint64_t progress = (i * 100) / nxyz;
                if (i % (nxyz / 100) == 0) {
                    std::cout << "    Progress: " << progress << "%\r" << std::flush;
                }

                *(nii_weight_data + i) = 0;
                if (*(nii_error_data + i) != 0) {
                    for (uint64_t t = 0; t < size_time; ++t) {
                        // Handle edge cases
                        uint64_t t_start, t_stop;
                        if ( (static_cast<int>(t) - vic ) < 0 ) {
                            t_start = 0;
                        } else {
                            t_start = t - static_cast<uint64_t>(vic);
                        }
                        if ( (static_cast<int>(t) + vic + 1 ) > static_cast<int>(size_time) ) {
                            t_stop = size_time;
                        } else {
                            t_stop = t + static_cast<uint64_t>(vic + 1);
                        }

                        uint64_t j = nxyz * t + i;
                        *(nii_smooth_error_data + j) = 0;
                        float weight = 0;
                        for (uint64_t jt = t_start; jt < t_stop; ++jt) {
                            uint64_t k = nxyz * jt + i;
                            float dist = static_cast<float>( std::abs(static_cast<int>(t) - static_cast<int>(jt)) );
                            float g = gaus(dist, gFWHM_val);
                            *(nii_smooth_error_data + j) += *(nii_error_data + k) * g ;
                            weight += g;
                        }
                        *(nii_smooth_error_data + j) /= weight;
                    }
                }
            }
            std::cout << std::endl;

            // Write out smoothed error time series   
            if (mode_debug) {
                save_output_nifti(fout, "smooth_error_maps", nii_smooth_error, false);
            }

            // Correct original data with smoothed error time series
            for (uint64_t t = 0; t < size_time; ++t) {
                for (uint64_t i = 0; i < nxyz; ++i) {
                    uint64_t j = nxyz * t + i;
                    if (t % 2 == 0) {
                        *(nii_smooth_data + j) = *(nii_data + j) - *(nii_smooth_error_data + j) * 0.5;
                    } else {
                        *(nii_smooth_data + j) = *(nii_data + j) + *(nii_smooth_error_data + j) * 0.5;
                    }
                }
            }

            std::cout << "    Saving..." << endl;
            save_output_nifti(fout, "frisgo-lpass", nii_smooth, false);
        }

        // ====================================================================
        // Run wise estimation of Fuzzy Ripple
        // ====================================================================
        if (mode_runwise) {
           std::cout << "  Doing run-wise correction..." << std::endl;
    
            // Allocating memory for average error 
            nifti_image* nii_mean_error = copy_nifti_as_float32(nii_weight);
            float* nii_mean_error_data = static_cast<float*>(nii_mean_error->data);
            for (uint64_t i = 0; i < nxyz ; ++i) {
                *(nii_mean_error_data + i) = 0;
            }

            // Average errors across time
            for (uint64_t i = 0; i < nxyz; ++i) {
                float sum = 0;
                for (uint64_t it = 0; it < size_time; ++it) {
                    sum += *(nii_error_data + nxyz * it + i);
                }
                *(nii_mean_error_data + i) = sum / static_cast<float>(size_time);
            }
        
            if (mode_debug) {
                save_output_nifti(fout, "mean_error_maps", nii_mean_error, false);
            }

            // Subtracting error from original data with inverse sign of every other TR
            for (uint64_t it = 0; it < size_time; ++it) {
                for (uint64_t i = 0; i < nxyz; ++i) {
                    uint64_t j = nxyz * it + i;
                    if (it % 2 == 0) {
                        *(nii_smooth_data + j) = *(nii_data + j) - ( *(nii_mean_error_data + i) ) * 0.5;
                    } else {
                        *(nii_smooth_data + j) = *(nii_data + j) + ( *(nii_mean_error_data + i) ) * 0.5;
                    }
                }
            }
            std::cout << "    Saving..." << std::endl;
            if (mode_oldname) {
                save_output_nifti(fout, "run", nii_smooth, false);  // Old convention 
            } else {
                save_output_nifti(fout, "frisgo-runwise", nii_smooth, false);
            }
        }
    }
    std::cout << "  Finished." << std::endl;
    return 0;
}
