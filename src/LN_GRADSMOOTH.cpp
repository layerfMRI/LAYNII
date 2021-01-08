#include "../dep/laynii_lib.h"


int show_help(void){
    printf(
    "LN_GRADSMOOTH : Local gradient based smoothing.\n"
    "\n"
    "Usage:\n"
    "    LN_GRADSMOOTH -input lo_VASO_act.nii -gradfile lo_gradT1.nii -FWHM 1 -within -selectivity 0.1 \n"
    "    ../LN_GRADSMOOTH -input lo_VASO_act.nii -gradfile lo_gradT1.nii -FWHM 1 -within -selectivity 0.1 \n"
    "\n"
    "Options:\n"
    "    -help        : Show this help\n"
    "    -input       : Nifti (.nii) that will be smoothed.\n"
    "    -gradfile    : Nifti (.nii) used to estimate local gradients.\n"
    "                   Only the first time point of this file is used. It \n"
    "                   should have the same spatial dimensions as the input.\n"
    "    -FWHM        : Amount of smoothing in mm.\n"
    "    -twodim      : (Optional) Smoothing in 2 Dim only \n"
    "    -mask        : (Optional) Nifti (.nii) that is used mask activity \n"
    "                   outside. This option can speed up processing.\n"
    "    -keep_masked : (Optional) Keep masked-out voxels in the output nifti.\n"
    "                   Useful for having a composite image of smoothed and\n"
    "                   un-smoothed voxels. Only used with -mask option.\n"
    "    -selectivity : (Optional) Makes the smoothing more or less specific \n"
    "                   to a certain gradient range. 0.05 is only within \n"
    "                   very similar values. 0.9 is almost independent of \n"
    "                   the gradient file 0.1 is default.\n"
    "    -within      : (Optional) Determines that smoothing should happen \n"
    "                   within similar values, not across different values.\n"
    "    -acros       : (Optional) Determines that smoothing should happen \n"
    "                   across different values, not within similar values.\n"
    "                   NOTE: This option is not working yet.\n"
    "    -output      : (Optional) Output filename, including .nii or\n"
    "                   .nii.gz, and path if needed. Overwrites existing files.\n"
    "\n"
    "Notes:\n"
    "    - If you run this on EPI-T1 data consider making them pretty, E.g:  \n"
    "      start_bias_field.sh T1.nii  \n"
    "      denoise_me.sh bico_T1.nii \n"
    "      short_me.sh denoised_bico_T1.nii \n"
    "      smooth_me.sh denoised_bico_T1.nii 0.5 \n"
    "      mv smoothed_denoised_bico_T1.nii new_T1.nii \n"
    "    - An example application is mentioned on the blog post here:\n"
    "      <https://layerfmri.com/anatomically-informed-spatial-smoothing> \n"
    "\n");
    return 0;
}

int main(int argc, char * argv[])
{
    bool use_outpath = false, keep_masked_voxels = false;
    char *fout = NULL;
    char *fgradi=NULL, *finfi=NULL, *fmaski=NULL;
    int ac, twodim=0, do_masking=0, within = 0, across = 0;
    float FWHM_val=0, selectivity=0.1;
    if( argc < 3 ) return show_help();

    // Process user options
    for( ac = 1; ac < argc; ac++ ) {
        if( !strncmp(argv[ac], "-h", 2) ) {
            return show_help();
        }
        else if( !strcmp(argv[ac], "-gradfile") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -gradfile\n");
                return 1;
            }
            fgradi = argv[ac];
        }
        else if( !strcmp(argv[ac], "-FWHM") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -FWHM\n");
                return 1;
            }
            FWHM_val = atof(argv[ac]);
        }
        else if( !strcmp(argv[ac], "-input") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            finfi = argv[ac];
        }
        else if( !strcmp(argv[ac], "-twodim") ) {
            twodim = 1;
            cout << "I will do smoothing only in 2D"  << endl;
        }
        else if( !strcmp(argv[ac], "-within") ) {
            within = 1;
            cout << "I will within similar values"  << endl;
        }
        else if( !strcmp(argv[ac], "-across") ) {
            across = 1;
            cout << "I will across different values"  << endl;
        }
        else if( !strcmp(argv[ac], "-mask") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -mask\n");
                return 1;
            }
            fmaski = argv[ac];
            do_masking = 1;
        }
        else if( !strcmp(argv[ac], "-keep_masked") ) {
            keep_masked_voxels = true;
        }
        else if( !strcmp(argv[ac], "-selectivity") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -selectivity\n");
                return 1;
            }
            selectivity = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else {
            fprintf(stderr,"** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!finfi) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fgradi) {
        fprintf(stderr, "** missing option '-gradfile'\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image *nim_inputfi = nifti_image_read(finfi, 1);
    if (!nim_inputfi) {
        fprintf(stderr,"** failed to read NIfTI from '%s'\n", finfi);
        return 2;
    }

    nifti_image *nim_gradi = nifti_image_read(fgradi, 1);
    if (!nim_gradi) {
        fprintf(stderr,"** failed to read NIfTI from '%s'\n", fgradi);
        return 2;
    }

    log_welcome("LN_GRADSMOOTH");
    log_nifti_descriptives(nim_inputfi);
    log_nifti_descriptives(nim_gradi);

    if (across + within !=1) {
        cout << " Please select either -within or -across" << endl;
        return 2;
    }
    if (across ==1) {
        cout << " -across is not implemented. Use -within instead." << endl;
        return 2;
    }

    // Get dimensions of input
    int size_x = nim_gradi->nx;
    int size_y = nim_gradi->ny;
    int size_z = nim_gradi->nz;
    int size_t = nim_inputfi->nt;
    int nx = nim_gradi->nx;
    int nxy = nim_gradi->nx * nim_gradi->ny;
    int nxyz = nim_gradi->nx * nim_gradi->ny * nim_gradi->nz;
    float dX = nim_gradi->pixdim[1];
    float dY = nim_gradi->pixdim[2];
    float dZ = nim_gradi->pixdim[3];

    // NOTE(Renzo): If you are running the smoothing in 2D, it will still go
    // thought he entire pipeline. The only difference is that the weights
    // in a certain direction are suppressed doing it in 2Dim, will not
    // speed up the program.
    if  (twodim == 1) {
        dZ = 1000 * dZ;
    }

    // ========================================================================
    // Fix datatype issues
    nifti_image *nim_inputf = copy_nifti_as_float32(nim_inputfi);
    float *nim_inputf_data = static_cast<float*>(nim_inputf->data);

    nifti_image *nim_grad = copy_nifti_as_float32(nim_gradi);
    float *nim_grad_data = static_cast<float*>(nim_grad->data);

    nifti_image *nim_roi = copy_nifti_as_float32(nim_gradi);
    float *nim_roi_data = static_cast<float*>(nim_roi->data);

    nifti_image *nim_mask = NULL;
    float *nim_mask_data = NULL;
    // ========================================================================

    if ( do_masking == 1 ) {
        // Read input dataset, including data
        nifti_image *nim_mask_input = nifti_image_read(fmaski, 1);
        if( !nim_mask_input ) {
            fprintf(stderr,"** failed to read NIfTI from '%s'\n", fmaski);
            return 2;
        }

        if (keep_masked_voxels) {
            cout << "  Masked-out voxel will be untouched instead of zero." << endl;
        } else {
            cout << "  Masked-out voxel will be zero." << endl;
        }

        // Quickfix for issue #29
        nim_mask = copy_nifti_as_float32(nim_mask_input);
        nim_mask_data = static_cast<float*>(nim_mask->data);

        for(int t=0; t<size_t; ++t) {
            for(int z=0; z<size_z; ++z) {
                for(int y=0; y<size_y; ++y) {
                    for(int x=0; x<size_x; ++x) {
                        int voxel_i = nxyz * t + nxy * z + nx * x + y;
                        *(nim_roi_data + voxel_i) = (float)(*(nim_mask_data  + voxel_i));
                    }
                }
            }
        }
    }

    // ========================================================================
    // MAKE allocating necessary files
    // ========================================================================
    nifti_image *smoothed = nifti_copy_nim_info(nim_inputf);
    nifti_image *gausweight = nifti_copy_nim_info(nim_inputf);
    smoothed->datatype = NIFTI_TYPE_FLOAT32;
    gausweight->datatype = NIFTI_TYPE_FLOAT32;
    smoothed->nbyper = sizeof(float);
    gausweight->nbyper = sizeof(float);
    smoothed->data = calloc(smoothed->nvox, smoothed->nbyper);
    gausweight->data = calloc(gausweight->nvox, gausweight->nbyper);

    // The gaus weight is a geometry factor and only needs to be estimated
    // once (not for every time step)
    gausweight->nt = 1;
    gausweight->nvox = gausweight->nvox / size_t;
    float *smoothed_data = (float*)smoothed->data;
    float *gausweight_data = (float*)gausweight->data;

    cout << "  Time dimension of smoothed output file:  " << smoothed->nt << endl;

    int vic = max(1.,2. * FWHM_val/dX );  // ignore if voxel is too far away
    float dist_i = 0.;
    cout << "  vic: " << vic <<  endl;
    cout << "  FWHM_val: " <<  FWHM_val<<  endl;

    float temp_wight_factor = 0;  // to store temp values, so I don't make the same computations over and over again.

    // ========================================================================
    // Finding the range of gradient values
    // ========================================================================

    // Values that I need to characterize the local signals in the vicinity.
    float local_val = 0;
    int NvoxInVinc = (2*vic+1)*(2*vic+1)*(2*vic+1);
    double vec1[NvoxInVinc];
    for(int it = 0; it < NvoxInVinc; it++) vec1[it] = 0;
    float grad_stdev = 0;
    float value_dist = 0;

    // For estimation and output of program process and how much longer it will take.
    int nvoxels_to_go_across = size_z * size_x * size_y;
    int running_index = 0;
    int pref_ratio = 0;

    if (size_z * size_x * size_y > 32767) {
        cout << "  Nr. voxels is bigger than the range of int the time estimation will be wrong." << endl;
    }

    if ( do_masking == 1 ) {
        nvoxels_to_go_across = 0;
        for(int iz=0; iz<size_z; ++iz) {
            for(int iy=0; iy<size_y; ++iy) {
                for(int ix=0; ix<size_x; ++ix) {
                    if (*(nim_roi_data + nxy * iz + nx * iy + ix) > 0 ) {
                        nvoxels_to_go_across += 1;
                    }
                }
            }
        }
    }

    cout << "  The number of voxels to go across = "<< nvoxels_to_go_across << endl;

    // ========================================================================
    // Smoothing loop
    // ========================================================================

    cout << "  Big smoothing loop is being done now..." << endl;

    for(int iz=0; iz<size_z; ++iz) {
        for(int iy=0; iy<size_x; ++iy) {
            for(int ix=0; ix<size_y; ++ix) {

                if ( !( !(*(nim_roi_data + nxy*iz + nx*ix + iy ) > 0) && (do_masking == 1)  ) ) {

                    // This is to write out how many more voxels I have to go through.
                    running_index++;
                    if ((running_index * 100) / nvoxels_to_go_across != pref_ratio ) {
                        cout << "\r    "<<(running_index * 100) / nvoxels_to_go_across << "% is done." << flush;
                        pref_ratio = (running_index * 100) / nvoxels_to_go_across;
                    }

                    // I am cooking in a clean kitchen.
                    *(gausweight_data + nxy * iz + nx * ix + iy)  = 0;
                    *(smoothed_data + nxy * iz + nx * ix + iy)  = 0;
                    NvoxInVinc = 0;
                    local_val = *(nim_grad_data + nxy * iz + nx * ix + iy);

                    // Examining the environment and determining what
                    // the signal intensities are and what its distribution are
                    for (int iz_i=max(0, iz-vic); iz_i<=min(iz+vic, size_z-1); ++iz_i) {
                        for (int iy_i=max(0, iy-vic); iy_i<=min(iy+vic, size_x-1); ++iy_i) {
                            for (int ix_i=max(0, ix-vic); ix_i<=min(ix+vic, size_y-1); ++ix_i) {
                                vec1[NvoxInVinc] = (double)*(nim_grad_data + nxy * iz_i + nx * ix_i + iy_i);
                                NvoxInVinc++;
                            }
                        }
                    }

                    // The standard deviation of the signal valued in the
                    // vicinity. This is necessary to normalize how many voxels
                    // are contributing to the local smoothing.
                    // grad_stdev = (float)gsl_stats_sd(vec1, 1, NvoxInVinc);
                    grad_stdev = (float ) ren_stdev (vec1, NvoxInVinc);

                    for(int iz_i=max(0, iz-vic); iz_i<=min(iz+vic, size_z-1); ++iz_i) {
                        for(int iy_i=max(0, iy-vic); iy_i<=min(iy+vic, size_x-1); ++iy_i) {
                            for(int ix_i=max(0, ix-vic); ix_i<=min(ix+vic, size_y-1); ++ix_i) {
                                dist_i = dist((float)ix, (float)iy, (float)iz, (float)ix_i, (float)iy_i, (float)iz_i,dX,dY,dZ);
                                value_dist = fabs(local_val - *(nim_grad_data + nxy * iz_i + nx * ix_i + iy_i));

                                temp_wight_factor = gaus(dist_i,FWHM_val ) * gaus(value_dist, grad_stdev * selectivity) / gaus(0, grad_stdev * selectivity);

                                // The gaus data are important to avoid local scaling differences, when the kernel size changes. E.g. at edge of images.
                                // this is a geometric parameter and only need to be calculated for one time point.
                                // this might be avoidable, if the gaus fucnction is better normaliced.
                                *(gausweight_data + nxy * iz + nx * ix + iy) += temp_wight_factor;

                                for(int it=0; it<size_t; ++it) {  // loop across lall time steps
                                    *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) += *(nim_inputf_data + nxyz * it + nxy * iz_i + nx * ix_i + iy_i) * temp_wight_factor;
                                }
                            }
                        }
                    }

                    // Scaling the signal intensity with the overall gaus leakage
                    if (*(gausweight_data  + nxy*iz + nx*ix  + iy  ) > 0 ) {
                        for(int it=0; it<size_t; ++it) {
                            *(smoothed_data + nxyz * it + nxy * iz + nx * ix + iy) /= *(gausweight_data + nxy * iz + nx * ix + iy);
                        }
                    }
                }
            }
        }
    }
    cout << endl;

    // Fill in masked-out voxel with original input values
    if (keep_masked_voxels) {
        for(int t=0; t<size_t; ++t) {
            for(int z=0; z<size_z; ++z) {
                for(int y=0; y<size_y; ++y) {
                    for(int x=0; x<size_x; ++x) {
                        int voxel_i = nxyz * t + nxy * z + nx * x + y;
                        if (*(nim_roi_data + voxel_i) == 0) {
                            *(smoothed_data + voxel_i) = *(nim_inputf_data + voxel_i);
                        }
                    }
                }
            }
        }
    }

    if (!use_outpath) {
        fout = finfi;
    }
    save_output_nifti(fout, "smoothed", smoothed, true, use_outpath);

    return 0;
}
