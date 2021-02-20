

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_SKEW: Calculates mean, standard deviation, tSNR, skew, kurtosis, \n"
    "         and autocorrelation of timeseries. This is helpful for\n"
    "         artifact hunting (e.g. ghosting).\n"
    "\n"
    "         It also claculates image SNR (Glover and Lai 1998 and Feinberg 2013).\n"
    "           The even- and odd-numbered time points of a fMRI time series are separately averaged, \n"
    "           and the sum and difference of these two images were calculated.\n"
    "           The spatial SETEV can be used as a proxy for thermal noise. \n"
    "\n"
    "Usage:\n"
    "    LN_SKEW -input Nulled_intemp.nii \n"
    "    ../LN_SKEW -input lo_BOLD_intemp.nii \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii or nii.gz) time series.\n"
    "    -output : (Optional) Output filename, including .nii or\n"
    "              .nii.gz, and path if needed. Overwrites existing files.\n"    
    "\n"
    "Notes:\n"
    "    Applications of this program are described in this blog post: \n"
    "    <http://layerfmri.com/QA>\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false ;
    char  *fout = NULL ;
    char *fin = NULL;
    int ac;
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

    // Read input dataset
    nifti_image * nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_SKEW");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_time = nii_input->nt;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocate new nifti
    nifti_image* nii_skew = nifti_copy_nim_info(nii);
    nii_skew->nt = 1;
    nii_skew->nvox = nii->nvox / size_time;
    nii_skew->datatype = NIFTI_TYPE_FLOAT32;
    nii_skew->nbyper = sizeof(float);
    nii_skew->data = calloc(nii_skew->nvox, nii_skew->nbyper);
    float* nii_skew_data = static_cast<float*>(nii_skew->data);

    nifti_image* nii_kurt = copy_nifti_as_float32(nii_skew);
    float* nii_kurt_data = static_cast<float*>(nii_kurt->data);

    nifti_image* nii_autocorr = copy_nifti_as_float32(nii_skew);
    float* nii_autocorr_data = static_cast<float*>(nii_autocorr->data);

    nifti_image* nii_conc = copy_nifti_as_float32(nii_skew);
    float* nii_conc_data = static_cast<float*>(nii_conc->data);

    nifti_image* nii_mean = copy_nifti_as_float32(nii_skew);
    float* nii_mean_data = static_cast<float*>(nii_mean->data);

    nifti_image* nii_stdev = copy_nifti_as_float32(nii_skew);
    float* nii_stdev_data = static_cast<float*>(nii_stdev->data);

    nifti_image* nii_tSNR = copy_nifti_as_float32(nii_skew);
    float* nii_tSNR_data = static_cast<float*>(nii_tSNR->data);
    
    nifti_image* nii_NOISE = copy_nifti_as_float32(nii_skew);
    float* nii_NOISE_data = static_cast<float*>(nii_NOISE->data);
    
    nifti_image* nii_GRAD = copy_nifti_as_float32(nii_skew);
    float* nii_GRAD_data = static_cast<float*>(nii_GRAD->data);

    nifti_image* nii_NOISESTDEV = copy_nifti_as_float32(nii_skew);
    float* nii_NOISESTDEV_data = static_cast<float*>(nii_NOISESTDEV->data);

    // ========================================================================
    cout << "  Calculating skew, kurtosis, and autocorrelation..." << endl;

    double vec1[size_time];
    double vec2[size_time];
    int voxel_i = 0; 

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_time; ++it) {
                    vec1[it] =
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_skew_data + voxel_i) = ren_skew(vec1, size_time);
                *(nii_kurt_data + voxel_i) = ren_kurt(vec1, size_time);
                *(nii_autocorr_data + voxel_i) = ren_autocor(vec1, size_time);
                *(nii_mean_data + voxel_i) =  ren_average(vec1, size_time);
                *(nii_stdev_data + voxel_i) = ren_stdev(vec1, size_time);
                *(nii_tSNR_data + voxel_i) = ren_average(vec1, size_time)/ren_stdev(vec1, size_time);
            }
        }
    }

    for (int voxel_i = 0; voxel_i < nxyz ; voxel_i++) {
      if ((nii_tSNR->scl_slope) != 0)  *(nii_tSNR_data + voxel_i) /=  (nii_tSNR->scl_slope) ; 
    }

    if (!use_outpath) fout = fin;

    save_output_nifti(fout, "skew", nii_skew, true);
    save_output_nifti(fout, "kurt", nii_kurt, true);
    save_output_nifti(fout, "autocorr", nii_autocorr, true);
    save_output_nifti(fout, "mean", nii_mean, true);
    save_output_nifti(fout, "stedev", nii_stdev, true);
    save_output_nifti(fout, "tSNR", nii_tSNR, true);
    // ========================================================================
    cout << "  Calculating correlation with everything..." << endl;

    for (int it = 0; it < size_time; ++it) {
        vec1[it] = 0;
        vec2[it] = 0;
    }

    // Mean time course of everything
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_time; ++it) {
                    vec1[it] +=
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i)
                                            / nxyz);
                }
            }
        }
    }

    // Voxel-wise corelation to mean of everything
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix <size_x; ++ix) {
                voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_time; ++it)   {
                    vec2[it] =
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_conc_data + voxel_i) = ren_correl(vec1, vec2, size_time);
            }
        }
    }
    save_output_nifti(fout, "overall_correl", nii_conc, true);
    
    
    
    
    // ========================================================================
    cout << "  Calculating image SNR ..." << endl;
    //size_time = 20;
    if (size_time%2 == 1) size_time = size_time -1  ;  // make sure its and odd number of time points 
    
    for (int it = 0; it < size_time-1 ; it = it + 2 )   {
        for (int voxel_i = 0; voxel_i < nxyz ; voxel_i++) {
                    //vec2[it] =  static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                *(nii_NOISE_data + voxel_i) += static_cast<double>(*(nii_data + nxyz * it     + voxel_i)) ;
                *(nii_NOISE_data + voxel_i) -= static_cast<double>(*(nii_data + nxyz * (it+1) + voxel_i)) ;
        }
    }
  
  
  // normalicing to time course duration
    for (int voxel_i = 0; voxel_i < nxyz ; voxel_i++) {
                *(nii_NOISE_data + voxel_i) = *(nii_NOISE_data + voxel_i) / sqrt((double) (size_time)/2 ) ;
    }
    
    
    save_output_nifti(fout, "noise", nii_NOISE, true);

//-------------------------------------
    // estimating local gradient of mean  
    
    int vinc_counter = 0 ;
    int vinc_x = 0 , vinc_y = 0 , vinc_z = 0;
    int vic = 1 ; // this will result in 26 neighbors (27 voxels) and is sufficient for decent STDEV estimation
    
     
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix <size_x; ++ix) {
                vinc_counter = 0; 
                
                for (int iz_i=max(0, iz-vic); iz_i<=min(iz+vic, size_z-1); ++iz_i) {
                    for (int iy_i=max(0, iy-vic); iy_i<=min(iy+vic, size_y-1); ++iy_i) {
                        for (int ix_i=max(0, ix-vic); ix_i<=min(ix+vic, size_x-1); ++ix_i) {
                                vec1[vinc_counter] = *(nii_mean_data + nxy * iz_i + nx * iy_i + ix_i) ;  
                                vinc_counter++;
                        }
                    }
                }
                *(nii_GRAD_data + nxy * iz + nx * iy + ix) = ren_stdev(vec1, vinc_counter); 
            }
        }
    }
    save_output_nifti(fout, "local_gradient", nii_GRAD, true);

    //-------------------------------------
    // estimateing local image SNR  


    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix <size_x; ++ix) {
                vinc_counter = 0; 
                
                for (int iz_i=max(0, iz-vic); iz_i<=min(iz+vic, size_z-1); ++iz_i) {
                    for (int iy_i=max(0, iy-vic); iy_i<=min(iy+vic, size_y-1); ++iy_i) {
                        for (int ix_i=max(0, ix-vic); ix_i<=min(ix+vic, size_x-1); ++ix_i) {
                                vec1[vinc_counter] = *(nii_NOISE_data + nxy * iz_i + nx * iy_i + ix_i) ;  
                                vinc_counter++;
                        }
                    }
                }
               // *(nii_NOISESTDEV_data + nxy * iz + nx * iy + ix) =  ren_stdev(vec1, vinc_counter); 
                *(nii_NOISESTDEV_data + nxy * iz + nx * iy + ix) = *(nii_mean_data + nxy * iz + nx * iy + ix) / ren_stdev(vec1, vinc_counter) ;
            }
        }
    }
    
    for (int voxel_i = 0; voxel_i < nxyz ; voxel_i++) {
        if ((nii_NOISESTDEV->scl_slope) != 0) *(nii_NOISESTDEV_data + voxel_i) /=  (nii_NOISESTDEV->scl_slope) ; 
    }

    save_output_nifti(fout, "imageSNR", nii_NOISESTDEV, true);

    cout << "  Finished." << endl;
    return 0;
}
