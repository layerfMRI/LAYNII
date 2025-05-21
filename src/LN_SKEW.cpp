#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN_SKEW: Calculates mean, standard deviation, tSNR, skew, kurtosis, \n"
    "         autocorrelation of timeseries and more. This is helpful for\n"
    "         artifact hunting (e.g. ghosting).\n"
    "\n"
    "Usage:\n"
    "    LN_SKEW -input Nulled_intemp.nii \n"
    "    ../LN_SKEW -input lo_BOLD_intemp.nii \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii or nii.gz) time series.\n"
    "    -output : (Optional) Output filename, including .nii or .nii.gz\n"
    "              and path if needed. Overwrites existing files.\n"    
    "\n"
    "Notes:\n"
    "    - Applications of this program are described in this blog post:\n"
    "      <http://layerfmri.com/QA>\n"
    "    - It also calculates image SNR (Glover and Lai 1998 and Feinberg 2013).\n"
    "    - The even and odd time points are separately averaged. The sum and\n"
    "      difference of these two images are calculated.\n"
    "    - The spatial STDEV can be used as a proxy for thermal noise.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    bool use_outpath = false;
    char  *fout = NULL;
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
    const uint64_t size_x = nii_input->nx;
    const uint64_t size_y = nii_input->ny;
    const uint64_t size_z = nii_input->nz;
    uint64_t size_time = nii_input->nt;
    const uint64_t nx = nii_input->nx;
    const uint64_t nxy = nii_input->nx * nii_input->ny;
    const uint64_t nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

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
    cout << "    Calculating skew, kurtosis, and autocorrelation..." << endl;

    std::vector<double> vec1(size_time);
    std::vector<double> vec2(size_time);
    std::vector<double> vecl(27);  // Vector for spatial gradient (number of voxel neighbours)

    uint64_t voxel_i = 0; 
    for (uint64_t iz = 0; iz < size_z; ++iz) {
        for (uint64_t iy = 0; iy < size_y; ++iy) {
            for (uint64_t ix = 0; ix < size_x; ++ix) {
                voxel_i = nxy * iz + nx * iy + ix;
                for (uint64_t it = 0; it < size_time; ++it) {
                    vec1[it] = static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_skew_data + voxel_i) = ren_skew(vec1.data(), size_time);
                *(nii_kurt_data + voxel_i) = ren_kurt(vec1.data(), size_time);
                *(nii_autocorr_data + voxel_i) = ren_autocor(vec1.data(), size_time);
                *(nii_mean_data + voxel_i) =  ren_average(vec1.data(), size_time);
                *(nii_stdev_data + voxel_i) = ren_stdev(vec1.data(), size_time);
                *(nii_tSNR_data + voxel_i) = ren_average(vec1.data(), size_time) / ren_stdev(vec1.data(), size_time);
            }
        }
    }

    for (uint64_t voxel_i = 0; voxel_i < nxyz; voxel_i++) {
        if ( (nii_tSNR->scl_slope) != 0 ) {
            *(nii_tSNR_data + voxel_i) /=  (nii_tSNR->scl_slope); 
        }
        if ( *(nii_tSNR_data + voxel_i) != *(nii_tSNR_data + voxel_i) ) {
            *(nii_tSNR_data + voxel_i) = 0; // filtering NaNs
        }
    }

    if (!use_outpath) fout = fin;

    save_output_nifti(fout, "skew", nii_skew, true);
    save_output_nifti(fout, "kurt", nii_kurt, true);
    save_output_nifti(fout, "autocorr", nii_autocorr, true);
    save_output_nifti(fout, "mean", nii_mean, true);
    save_output_nifti(fout, "stdev", nii_stdev, true);
    save_output_nifti(fout, "tSNR", nii_tSNR, true);

    // ========================================================================
    cout << "    Calculating correlation with everything..." << endl;

    for (uint64_t it = 0; it < size_time; ++it) {
        vec1[it] = 0;
        vec2[it] = 0;
    }

    // Mean time course of everything
    for (uint64_t iz = 0; iz < size_z; ++iz) {
        for (uint64_t iy = 0; iy < size_y; ++iy) {
            for (uint64_t ix = 0; ix < size_x; ++ix) {
                voxel_i = nxy * iz + nx * iy + ix;
                for (uint64_t it = 0; it < size_time; ++it) {
                    vec1[it] += static_cast<double>(*(nii_data + nxyz * it + voxel_i) / nxyz);
                }
            }
        }
    }

    // Voxel-wise corelation to mean of everything
    for (uint64_t iz = 0; iz < size_z; ++iz) {
        for (uint64_t iy = 0; iy < size_y; ++iy) {
            for (uint64_t ix = 0; ix <size_x; ++ix) {
                voxel_i = nxy * iz + nx * iy + ix;
                for (uint64_t it = 0; it < size_time; ++it)   {
                    vec2[it] = static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_conc_data + voxel_i) = ren_correl(vec1.data(), vec2.data(), size_time);
            }
        }
    }
    save_output_nifti(fout, "overall_correl", nii_conc, true);

    // ========================================================================
    cout << "    Calculating image SNR ..." << endl;
    if (size_time%2 == 1) {
        size_time = size_time - 1;  // Make sure it has even number of time points 
    }

    for (uint64_t it = 0; it < size_time-1; it = it + 2 )   {
        for (uint64_t voxel_i = 0; voxel_i < nxyz; voxel_i++) {
            *(nii_NOISE_data + voxel_i) += static_cast<double>(*(nii_data + nxyz * it     + voxel_i));
            *(nii_NOISE_data + voxel_i) -= static_cast<double>(*(nii_data + nxyz * (it+1) + voxel_i));
        }
    }
  
    // Normalizing to time course duration
    for (uint64_t voxel_i = 0; voxel_i < nxyz; ++voxel_i) {
        *(nii_NOISE_data + voxel_i) = *(nii_NOISE_data + voxel_i) / sqrt((double) (size_time)/2 );
    }

    save_output_nifti(fout, "noise", nii_NOISE, true);

    // ------------------------------------------------------------------------
    // Estimating local gradient of mean  
    int64_t vic_counter = 0;  // Vicinity counter
    int64_t vic = 1;  // This will result in 26 neighbors (27 voxels) and is sufficient for decent STDEV estimation
    int64_t sx = static_cast<int64_t>(size_x);
    int64_t sy = static_cast<int64_t>(size_y);
    int64_t sz = static_cast<int64_t>(size_z);
    for (int64_t iz = 0; iz < sz; ++iz) {
        for (int64_t iy = 0; iy < sy; ++iy) {
            for (int64_t ix = 0; ix <sx; ++ix) {
                vic_counter = 0; 
                for (int64_t iz_i = std::max(int64_t(0), iz-vic); iz_i <= std::min(iz+vic, sz-int64_t(1)); ++iz_i) {
                    for (int64_t iy_i = std::max(int64_t(0), iy-vic); iy_i <= std::min(iy+vic, sy-int64_t(1)); ++iy_i) {
                        for (int64_t ix_i = std::max(int64_t(0), ix-vic); ix_i <= std::min(ix+vic, sx-int64_t(1)); ++ix_i) {
                            vecl[vic_counter] = *(nii_mean_data + nxy * iz_i + nx * iy_i + ix_i);  
                            vic_counter++;
                        }
                    }
                }
                *(nii_GRAD_data + nxy * iz + nx * iy + ix) = ren_stdev(vecl.data(), vic_counter); 
            }
        }
    }
    save_output_nifti(fout, "local_gradient", nii_GRAD, true);
    cout << "    Estimating local image SNR ..." << endl;

    //-------------------------------------------------------------------------
    // Estimating local image SNR
    for (int64_t iz = 0; iz < sz; ++iz) {
        for (int64_t iy = 0; iy < sy; ++iy) {
            for (int64_t ix = 0; ix <sx; ++ix) {
                vic_counter = 0; 
                for (int64_t iz_i = std::max(int64_t(0), iz-vic); iz_i <= std::min(iz+vic, sz-int64_t(1)); ++iz_i) {
                    for (int64_t iy_i = std::max(int64_t(0), iy-vic); iy_i <= std::min(iy+vic, sy-int64_t(1)); ++iy_i) {
                        for (int64_t ix_i = std::max(int64_t(0), ix-vic); ix_i <= std::min(ix+vic, sx-int64_t(1)); ++ix_i) {
                            vecl[vic_counter] = *(nii_NOISE_data + nxy * iz_i + nx * iy_i + ix_i);  
                            vic_counter++;
                        }
                    }
                }
                // *(nii_NOISESTDEV_data + nxy * iz + nx * iy + ix) =  ren_stdev(vec1.data(), vic_counter); 
                *(nii_NOISESTDEV_data + nxy * iz + nx * iy + ix) = *(nii_mean_data + nxy * iz + nx * iy + ix) / ren_stdev(vecl.data(), vic_counter);
            }
        }
    }

    for (uint64_t voxel_i = 0; voxel_i < nxyz; voxel_i++) {
        if ( (nii_NOISESTDEV->scl_slope) != 0 ) {
            *(nii_NOISESTDEV_data + voxel_i) /= (nii_NOISESTDEV->scl_slope);
        } 
        if ( *(nii_NOISESTDEV_data + voxel_i) != *(nii_NOISESTDEV_data + voxel_i) ) {
            *(nii_NOISESTDEV_data + voxel_i) = 0;
        }
    }

    save_output_nifti(fout, "imageSNR", nii_NOISESTDEV, true);

    cout << "Finished." << endl;
    return 0;
}
