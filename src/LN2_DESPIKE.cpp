#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN_DESPIKE: Voxel-wise despiking for timeseries based on robust z-scores.\n"
    "            Outliers are imputed with simple averaging of the neighboring\n"
    "            time points (t-1 t+1).\n"
    "\n"
    "Usage:\n"
    "    LN_DESPIKE -input Nulled_intemp.nii \n"
    "    ../LN_DESPIKE -input lo_BOLD_intemp.nii \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii or nii.gz) time series.\n"
    "    -sigma  : Voxels with more extreme z-scores than this parameter will be imputed.\n"
    "              Both postive and negative extremes are considered. Default is '2.58'\n"
    "              which means the imputed value is more extreme than 99%% of other values\n"
    "              in that voxel's timecourse.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "    -output : (Optional) Output filename, including .nii or .nii.gz\n"
    "              and path if needed. Overwrites existing files.\n"    
    "\n"
    "\n");
    return 0;
}


int main(int argc, char * argv[]) {
    nifti_image *nii1 = NULL;
    bool use_outpath = false, mode_debug = false;
    char *fin1 = NULL, *fout = NULL;
    float sigma = 3;
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
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-sigma")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -sigma\n");
                return 1;
            }
            sigma = std::stoi(argv[ac]);
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            use_outpath = true;
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    log_welcome("LN2_ZSCORE");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;
    const uint32_t size_time = nii1->nt;

    const uint64_t end_x = size_x - 1;
    const uint64_t end_y = size_y - 1;
    const uint64_t end_z = size_z - 1;
    const uint64_t end_time = size_time - 1;

    const uint64_t nxyz = size_z * size_y * size_x;
    const uint64_t nxyzt = size_z * size_y * size_x * size_time;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii_input = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Allocate new nifti for 4D output
    nifti_image* nii_zscore = copy_nifti_as_float32(nii_input);
    float* nii_zscore_data = static_cast<float*>(nii_zscore->data);

    nifti_image* nii_mask = copy_nifti_as_int16(nii_input);
    int16_t* nii_mask_data = static_cast<int16_t*>(nii_mask->data);

    nifti_image* nii_output = copy_nifti_as_float32(nii_input);
    float* nii_output_data = static_cast<float*>(nii_output->data);

    // Allocate new nifti for 3D outputs
    nifti_image* nii_median = nifti_copy_nim_info(nii_input);
    nii_median->nt = 1;
    nii_median->nvox = nxyz;
    nii_median->datatype = NIFTI_TYPE_FLOAT32;
    nii_median->nbyper = sizeof(float);
    nii_median->data = calloc(nii_median->nvox, nii_median->nbyper);
    float* nii_median_data = static_cast<float*>(nii_median->data);

    nifti_image* nii_mad = copy_nifti_as_float32(nii_median);
    float* nii_mad_data = static_cast<float*>(nii_mad->data);

    // Voxels of interest that will be imputed
    nifti_image* nii_counts = copy_nifti_as_int16(nii_median);
    int16_t* nii_counts_data = static_cast<int16_t*>(nii_counts->data);

    // ========================================================================
    cout << "  Computing voxel-wise medians..." << endl;
    for (uint64_t i = 0; i < nxyz; ++i) {
        vector <float> temp_vec(size_time);
        for (uint64_t t = 0; t < size_time; ++t) {
            temp_vec[t] = *(nii_input_data + i + t*nxyz);
        }

        // Calculate Median
        std::sort(temp_vec.begin(), temp_vec.end());
        if (size_time % 2 != 0) {  // Odd number of elements: return the middle one
            *(nii_median_data + i) =  temp_vec[size_time / 2];
        } else {  // Even number of elements: average the two middle ones
            *(nii_median_data + i) = ( temp_vec[size_time / 2 - 1] + temp_vec[size_time / 2] ) / 2.0;
        }
    }

    cout << "    Saving..." << endl;
    save_output_nifti(fout, "despike-median", nii_median, false);

    // ========================================================================
    cout << "  Computing median absolute deviations (MAD)..." << endl;
    for (uint64_t i = 0; i < nxyz; ++i) {
        vector <float> temp_vec(size_time);
        for (uint64_t t = 0; t < size_time; ++t) {
            temp_vec[t] = std::abs( *(nii_input_data + i + t*nxyz) - *(nii_median_data + i) );
        }

        // Calculate Median
        std::sort(temp_vec.begin(), temp_vec.end());
        if (size_time % 2 != 0) {  // Odd number of elements: return the middle one
            *(nii_mad_data + i) =  temp_vec[size_time / 2];
        } else {  // Even number of elements: average the two middle ones
            *(nii_mad_data + i) = ( temp_vec[size_time / 2 - 1] + temp_vec[size_time / 2] ) / 2.0;
        }
    }

    if (mode_debug) {
        cout << "    Saving..." << endl;
        save_output_nifti(fout, "despike-mad", nii_mad, false);
    }

    // ========================================================================
    cout << "  Computing robust z-scores..." << endl;
    for (uint64_t t = 0; t < size_time; ++t) {
        for (uint64_t i = 0; i < nxyz; ++i) {
            *(nii_zscore_data + i + t*nxyz) = *(nii_input_data + i + t*nxyz) - *(nii_median_data + i);
            *(nii_zscore_data + i + t*nxyz) *= 0.6745;
            *(nii_zscore_data + i + t*nxyz) /= *(nii_mad_data + i);
        }
    }

    if ( mode_debug ) {
        cout << "    Saving..." << endl;
        save_output_nifti(fout, "despike-robust_zscore", nii_zscore, false);
    }

    // ========================================================================
    cout << "  Finding extreme voxels..." << endl;
    float counter_total = 0;
    for (uint64_t t = 0; t < size_time; ++t) {
        for (uint64_t i = 0; i < nxyz; ++i) {
            if (*(nii_zscore_data + i + t*nxyz) >= sigma || *(nii_zscore_data + i + t*nxyz) <= -sigma) {
                *(nii_mask_data + i + t*nxyz) = 1;
                *(nii_counts_data + i) += 1;
                counter_total += 1;
            } else {
                *(nii_mask_data + i + t*nxyz) = 0;
            }
        }
    }

    cout << "    Saving..." << endl;
    save_output_nifti(fout, "despike-counts", nii_counts, false);

    // ========================================================================
    cout << "  Imputing..." << endl;
    uint64_t ix, iy, iz, it;
    for (uint64_t i = 0; i < nxyzt; ++i) {

        // Impute with average of the neighboring time points
        if (*(nii_mask_data + i) != 0) {
            std::tie(ix, iy, iz, it) = ind2sub_4D_64(i, size_x, size_y, size_z);

            if (it > 0 && it < end_time) {
                uint64_t j = sub2ind_4D_64(ix, iy, iz, it-1, size_x, size_y, size_z);
                uint64_t k = sub2ind_4D_64(ix, iy, iz, it+1, size_x, size_y, size_z);
                *(nii_output_data + i) = ( *(nii_input_data + j) + *(nii_input_data + k) ) / 2.0;
            } else if (it == end_time) {
                uint64_t j = sub2ind_4D_64(ix, iy, iz, it-1, size_x, size_y, size_z);
                *(nii_output_data + i) = *(nii_input_data + j);
            } else if (it == 0) {
                uint64_t k = sub2ind_4D_64(ix, iy, iz, it+1, size_x, size_y, size_z);
                *(nii_output_data + i) = *(nii_input_data + k);
            }
        }
    }

    cout << "    Saving..." << endl;
    save_output_nifti(fout, "despike-simple", nii_output, false);

    cout << "Finished." << endl;
    return 0;
}
