

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_FRISGO :  This program amins to correct for Fuzzy Ripple \n"
    "              artifacts in dual polarity EPI data.\n"
    "              \n"
    "\n"
    "Usage:\n"
    "    LN2_FRISGO -input timeseries.nii -lpass 1.0 \n"
    "    LN2_FRISGO -input timeseries.nii -box 1 \n"
    "    ../LN2_FRISGO -input lo_BOLD_intemp.nii -box 1 \n" 
    "    ../LN2_FRISGO -input lo_BOLD_intemp.nii -lpass 1 \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii) file with time series data that will be \n"
    "              nii_smooth. Only the first time point is used. \n"
    "    -tshift : This option of dual polarity correction estimates \n"
    "              the fMRI signal rigth between TRs. Just like slice  \n"
    "              timing correction, but for 3D-EPI This way the triger \n"
    "              timing is directly aligned with the k-space center.\n"
    "              Here we use a third order interpolation across time \n"
    "              to minimize temporal blurring.\n"
    "              \n"
    "    -lpass  : This option of dual polarity correction estimates slow \n"
    "              changes of the Fruzzy ripple across time. \n"
    "              A Gaussian travelling window of averaging is used to \n"
    "              to estimate the Fuzzy Ripple modulations. Specify the value \n"
    "              of the Gaussian size (float values) in units of TR. \n"
    "              Make sure to use values that are bigger than you trial timing.\n"
    "              The default value is 40 TRs.\n"
    "    -run    : Doing the dual-polarity correction by assuming that the \n"
    "              Fuzzy Ripples are indenital for the entire run. \n"
    "    -calib  : This option assumes that the first 4 volumes have alternating \n"
    "              read polarity, while the remainer of the run has one polarity (WIP). \n"
    "    -output : (Optional) Output filename, including .nii or\n"
    "              .nii.gz, and path if needed. Overwrites existing files.\n"    
     "    -verb  : wrtiting out all the inbetween steps, for debugging. \n"
    "\n"
    "Notes:\n"
    "    This program was featured in the following ISMRM presentation:\n"
    "    <https://youtu.be/a-X1JW7-Pk4?si=QRG_wgqPiO74Zphe> \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false ;
    char  *fout = NULL ;
    char* fin = NULL;
    int ac, do_gaus = 0, do_run = 0, do_tshift = 0, bFWHM_val = 0, do_verb = 0;
    float gFWHM_val = 0.0;
    if (argc  <  3) return show_help();

    // Process user options
    for (ac = 1; ac  <  argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-lpass")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -lpass\n");
                return 1;
            }
            gFWHM_val = atof(argv[ac]);
            do_gaus = 1;
        } else if (!strcmp(argv[ac], "-run")) {
            do_run = 1;
        } else if (!strcmp(argv[ac], "-verb")) {
            do_verb = 1;
        } else if (!strcmp(argv[ac], "-tshift")) {
            do_tshift = 1;
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    if (!use_outpath) fout = fin;

    // Read input dataset
    nifti_image * nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    if ( gFWHM_val ==0 & do_gaus == 1) {
        cout << "  *******************************************************" << endl;
        cout << "  You did not specify a value for -lpass, so I will use the default value of 40 TRs." << endl;
        cout << "  If you want to change this, please specify a value for -lpass." << endl;
        cout << "  *******************************************************" << endl;
        gFWHM_val = 40.0; // Default value
    } else {
        cout << "  Using Gaussian low pass filter with FWHM of " << gFWHM_val << " TRs." << endl; 
    }

    if (do_run + do_gaus + do_tshift == 0) {
        cout << "  *******************************************************" << endl;
        cout << "  Invalid smoothing option. Select at least one: lpass OR run OR tshift," << endl;
        cout << "  of you selected neither, I dont know what to do." << endl;
        cout << "  *******************************************************" << endl;
        return 2;
    }


    log_welcome("LN2_FRISGO");
    log_nifti_descriptives(nii_input);
    if (do_gaus) {
        cout << "Selected temporal low pass filter: Gaussian" << endl;
    } else if (do_run) {
        cout << "Selected run wise estimation and correction of Fuzzy Ripples" << endl;
    }   else if (do_tshift) {
        cout << "Selected slice timing correction for Fuzzy Ripples" << endl;
    }

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_t = nii_input->nt;
    int nr_voxels = size_x * size_y * size_z;
    // int nx = nii_input->nx;
    // int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;
    int nxyzt = nii_input->nx * nii_input->ny * nii_input->nz * nii_input->nt;
    float dT = 1;

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

    if (do_verb) cout << "done allocating the main stuff, not going though the spefics algorythm as desired " << endl;

    // ========================================================================
    // Fuzzy ripple correction via slice time correction
    // ========================================================================
cout << "run until here -2 " << endl;

    for (int t = 0; t < size_t; ++t) {
            for (int i = 0; i < nr_voxels; ++i) {
                *(nii_smooth_data + t*nxyz + i) = 0.0 ;
            }
        }
cout << "run until here -1 " << endl;

    if (do_tshift) {
        
        // convention is to fist forward
        // see assumed time table below
        // t is trigger
        // k time point of largest power of fmri signal: k-scpace center for acquisition   // 
        // input data are like this:  
        // t       t       t       t       t       t       t       ....
        //     k1      k2      k3      k4      k5      k6      k7      ....
        //
        // output data are like this: 
        // t       t       t       t       t       t       t       ....
        // k1      k1-2    k2-3    k3-4    k4-5    k5-6    k6-7    ....

        // if there are less than 4 TRs, there is nothing to do for me.
        if (size_t<4) {
            cout << "###############################" << endl;
            cout << "you have less than 4 TRs, so I canot reallt do 3rd order spline fitting with this" << endl;
            cout << "go back to the scanner and collect longer time series" << endl;
            return 1; 
        }
cout << "run until here0 " << endl;
        // dealing with first two time point, interpolating back
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_smooth_data + nxyz * 0 + i) = 0.5 * ( *(nii_data + nxyz * 0 + i) + *(nii_data + nxyz * 1 + i) ) ;
            *(nii_smooth_data + nxyz * 1 + i) = 0.5 * ( *(nii_data + nxyz * 1 + i) + *(nii_data + nxyz * 2 + i) ) ;
        }
        cout << "run until here1 " << endl;
        // dealing with last two time points copying averages
        // removed for now because I had segmegtnation errors here
        for (int i = 0; i < nr_voxels; ++i) {
         //   *(nii_smooth_data + nxyz * (size_t-3) + i) = 0.5 * ( *(nii_data + nxyz * (size_t-3) + i) + *(nii_data + nxyz * (size_t-4) + i) );
         //   *(nii_smooth_data + nxyz * (size_t-4) + i) = 0.5 * ( *(nii_data + nxyz * (size_t-3) + i) + *(nii_data + nxyz * (size_t-2) + i) ) ;
        }
        cout << "run until here2 " << endl;
        // dealing with all the time points in between
        // y = a + b x + c x^2 + d x^3, for 4 data points, 
        //forward matrix is y = M A, 
        //with y = (y1, y2, y3, y4) the signal strenghts across 4 time point
        //with A = a, b, c, d, the prefactors of the cubic polinomial 
        //y (X=0) is desired, so the other x values are -1.5, -0.5, 0.5, 1.5
        //so the forward matrix is M = (1, -1.5, 2.25, -3.375; 
        //                              1, -0.5, 0.25, -0.125; 
        //                              1, 0.5,  0.25,  0.125; 
        //                              1, 1.5,  2.25,  3.375)
        //so the inverse matrix is M^-1 = (- 0.0625, + 0.5625, + 0.5625, - 0.0625 ) .... (we only need th first row).
        // and the output is y = M^-1 y, so we can just multiply the input data with the inverse matrix.

        for (int t = 2; t < size_t-2; ++t) {
            for (int i = 0; i < nr_voxels; ++i) {
                *(nii_smooth_data + t*nxyz + i) = // 0.5 * ( *(nii_data + nxyz * (t-2) + i) + *(nii_data + nxyz * (t-1) + i) ) ;
                                              - 0.0625 *  ( *(nii_data + (t-2)*nxyz + i))
                                              + 0.5625 *  ( *(nii_data + (t-1)*nxyz + i))
                                              + 0.5625 *  ( *(nii_data + (t)*nxyz + i))
                                              - 0.0625 *  ( *(nii_data + (t+1)*nxyz + i));
            }
        }
        save_output_nifti(fout, "tshift_corrected", nii_smooth, true, use_outpath);
    }



    // ========================================================================
    // Preparing things I need for runwise and low pass correcfted ways of Fuzzy Ripple correction
    // ========================================================================
    if(do_gaus || do_run) {
        if (do_verb) cout << "preparing for runwise and low pass correction" << endl;

        // allocating memory for weight time series
        nifti_image* nii_weight = nifti_copy_nim_info(nii);
        nii_weight->nt = 1;
        nii_weight->datatype = NIFTI_TYPE_FLOAT32;
        nii_weight->nbyper = sizeof(float);
        nii_weight->nvox = size_x * size_y * size_z;
        nii_weight->data = calloc(nii_weight->nvox, nii_weight->nbyper);
        float* nii_weight_data = static_cast<float*>(nii_weight->data);

    
        // error time series
        nifti_image* nii_error = copy_nifti_as_float32(nii);
        float* nii_error_data = static_cast<float*>(nii_error->data);
        if (do_verb) cout << "done allocating error time series" << endl;
        
        // checking if the allocation was successful
        if (!nii_weight_data) {
            std::cerr << "Allocation failed for nii_weight_data\n";
            return 1;
        }
        if (!nii_smooth_data) {
            std::cerr << "Allocation failed for nii_smooth_data\n";
            return 1;
        }

        // setting all error data to zero
        for (int i = 0; i < nxyzt ; ++i)   *(nii_error_data + i) = 0;

        if (do_verb) cout << "calculating error across time: first " << endl;
        // getting error of first time point
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_error_data + i) = *(nii_data + nxyz * 0 + i) - *(nii_data + nxyz * 1 + i);
        }

        if (do_verb) cout << "calculating error across time: last" << endl;
        // getting error of last time point
        // removed - was causing segmentation faults
        for (int i = 0; i < nr_voxels; ++i) {
        // *(nii_error_data + (nxyz*(size_t-1)) + i) = *(nii_data + (nxyz*(size_t-1)) + i) - *(nii_data + (nxyz*(size_t-2)) + i);
        }
        
        if (do_verb) cout << "calculating error across time: all the rest" << endl;
        // getting error for all the other time points in between
        for (int t = 1; t < size_t-1; ++t) {
            for (int i = 0; i < nr_voxels; ++i) {
                *(nii_error_data + t*nxyz + i) =  *(nii_data + t*nxyz + i) 
                                                - 0.5 *  ( *(nii_data + (t-1)*nxyz + i))
                                                - 0.5 *  ( *(nii_data + (t+1)*nxyz + i));
            }
        }

        // invertign error every other TR
            for (int t = 1; t < size_t; t += 2) {
            for (int i = 0; i < nr_voxels; ++i) {
                    *(nii_error_data + t*nxyz + i) *= -1;
                }
            }

        if (do_verb) save_output_nifti(fout, "error_maps", nii_error, true, use_outpath);


        // ========================================================================
        // Run wise estimation of Fuzzy Ripple
        // ========================================================================
        
        if (do_run) {
            if (do_verb) cout << "starting run wise Fuzzy Ripple correction" << endl;
            // allocating memory for average error 
            nifti_image* nii_mean_error = copy_nifti_as_float32(nii_weight);
            float* nii_mean_error_data = static_cast<float*>(nii_mean_error->data);
            for (int i = 0; i < nxyz ; ++i)   *(nii_mean_error_data + i) = 0;

            // getting average of error across time
            for (int i = 0; i < nr_voxels; ++i) {
                float sum = 0;
                for (int it = 0; it < size_t; ++it) {
                    sum += *(nii_error_data + nxyz * it + i);
                }
                *(nii_mean_error_data + i) = sum / size_t;
            }
        
            if (do_verb) save_output_nifti(fout, "mean_error_maps", nii_mean_error, true, use_outpath);

            // subtracting error from original data with inverse sign of every other TR
            for (int it = 0; it < size_t; ++it) {
                for (int i = 0; i < nr_voxels; ++i) {
                    int j = nxyz * it + i;
                    if (it % 2 == 0) {
                        *(nii_smooth_data + j) = *(nii_data + j) - ( *(nii_mean_error_data + i)) * 0.5 ;
                    } else {
                        *(nii_smooth_data + j) = *(nii_data + j) + ( *(nii_mean_error_data + i)) * 0.5 ;
                    }
                }
            }
            save_output_nifti(fout, "runwise_corrected", nii_smooth, true, use_outpath);


        }


        // ========================================================================
        // Smoothing loop
        // ========================================================================
        int vic;
        if (do_gaus) {
            vic = max(1., 2. * gFWHM_val / dT);  // Ignore if voxel is too far

            cout << "    vic " << vic << endl;
            cout << "    FWHM_val " << gFWHM_val << endl;

            // allocating memory for smoothed error time series.
            nifti_image* nii_smooth_error = copy_nifti_as_float32(nii);
            float* nii_smooth_error_data = static_cast<float*>(nii_smooth_error->data);

            cout << "starting smoothing loop....This might take a while" << endl;
            cout << "progress is shown in percent" << endl;
            cout << "  " << flush;
            for (int i = 0; i < nr_voxels; ++i) {

                if (i % (nr_voxels / 100) == 0) {
                    cout << "\b\b\b\b" << i * 100 / nr_voxels << "%" << flush;
                }

                *(nii_weight_data + i) = 0;
                if (*(nii_error_data + i) != 0) {
                    for (int it = 0; it < size_t; ++it) {
                        int j = nxyz * it + i;
                        *(nii_smooth_error_data + j) = 0;

                            float weight = 0;
                            int jt_start = max(0, it - vic);
                            int jt_stop = min(it + vic + 1, size_t);
                            for (int jt = jt_start; jt < jt_stop; ++jt) {
                                int k = nxyz * jt + i;
                                float dist = abs(it - jt);
                                float g = gaus(dist, gFWHM_val);
                                *(nii_smooth_error_data + j) += (*(nii_error_data + k) * g);
                                weight += g;
                            }
                            *(nii_smooth_error_data + j) /= weight;
                    }
                }
            }

            // wrtitng out smoothed error time series   
            if (do_verb) save_output_nifti(fout, "smooth_error_maps", nii_smooth_error, true, use_outpath);

             // correcting original data with smoothed error time series
            for (int it = 0; it < size_t; ++it) {
                for (int i = 0; i < nr_voxels; ++i) {
                    int j = nxyz * it + i;
                    if (it % 2 == 0) {
                        *(nii_smooth_data + j) = *(nii_data + j) - (*(nii_smooth_error_data + j)) * 0.5 ;
                    } else {
                        *(nii_smooth_data + j) = *(nii_data + j) + (*(nii_smooth_error_data + j)) * 0.5 ;
                    }
                }
            }
            
            save_output_nifti(fout, "lpass_corrected", nii_smooth, true, use_outpath);

        } // smoothing loop closed

    } 
    cout << "  Finished." << endl;
    return 0;
}
