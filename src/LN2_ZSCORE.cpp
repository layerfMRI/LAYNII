#include "../dep/laynii_lib.h"


int show_help(void) {
    printf(
    "LN_ZSCORE: Standardize timeseries (z-score). It is calculated by subtracting\n"
    "           the voxel-wise mean from each voxel value and then dividing the \n"
    "           difference by the voxel-wise standard deviation.\n"
    "\n"
    "Usage:\n"
    "    LN_ZSCORE -input Nulled_intemp.nii \n"
    "    ../LN_ZSCORE -input lo_BOLD_intemp.nii \n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii or nii.gz) time series.\n"
    "    -mean   : (Optional) Output voxel-wise mean.\n"
    "    -std    : (Optional) Output voxel-wise standard deviation.\n"     
    "    -output : (Optional) Output filename, including .nii or .nii.gz\n"
    "              and path if needed. Overwrites existing files.\n"    
    "\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    nifti_image *nii1 = NULL;
    bool use_outpath = false, mode_mean=false, mode_std=false;
    char *fin1 = NULL, *fout = NULL;
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
        } else if (!strcmp(argv[ac], "-mean")) {
            mode_mean = true;
        } else if (!strcmp(argv[ac], "-std")) {
            mode_std = true;
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

    const uint32_t nxyz = size_z * size_y * size_x;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii_input = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Allocate new nifti for 4D output
    nifti_image* nii_zscore = copy_nifti_as_float32(nii_input);
    float* nii_zscore_data = static_cast<float*>(nii_zscore->data);

    // Allocate new nifti for 3D outputs
    nifti_image* nii_mean = nifti_copy_nim_info(nii_input);
    nii_mean->nt = 1;
    nii_mean->nvox = nxyz;
    nii_mean->datatype = NIFTI_TYPE_FLOAT32;
    nii_mean->nbyper = sizeof(float);
    nii_mean->data = calloc(nii_mean->nvox, nii_mean->nbyper);
    float* nii_mean_data = static_cast<float*>(nii_mean->data);

    nifti_image* nii_std = copy_nifti_as_float32(nii_mean);
    float* nii_std_data = static_cast<float*>(nii_std->data);

    // ========================================================================
    cout << "  Computing voxel-wise mean..." << endl;
    const float N = static_cast<float>(size_time);
    for (uint64_t t = 0; t < size_time; ++t) {
        for (uint64_t i = 0; i < nxyz; ++i) {
            *(nii_mean_data + i) += *(nii_input_data + i + t*nxyz) / N;
        }
    }

    if (mode_mean) {
        cout << "    Saving..." << endl;
        save_output_nifti(fout, "mean", nii_mean, false);        
    }

    // ========================================================================
    cout << "  Computing voxel-wise standard deviation..." << endl;
    for (uint64_t t = 0; t < size_time; ++t) {
        for (uint64_t i = 0; i < nxyz; ++i) {
            float temp = *(nii_mean_data + i) - *(nii_input_data + i + t*nxyz);
            *(nii_std_data + i) += (temp*temp) / N;
        }
    }

    for (uint64_t i = 0; i < nxyz; ++i) {
        *(nii_std_data + i) = sqrt(*(nii_std_data + i));
    }

    if (mode_std) {
        cout << "    Saving..." << endl;
        save_output_nifti(fout, "std", nii_std, false);
    }

    // ========================================================================
    cout << "  Z-scoring..." << endl;
    for (uint64_t t = 0; t < size_time; ++t) {
        for (uint64_t i = 0; i < nxyz; ++i) {
            *(nii_zscore_data + i + t*nxyz) = *(nii_input_data + i + t*nxyz) - *(nii_mean_data + i);
            *(nii_zscore_data + i + t*nxyz) /= *(nii_std_data + i);
        }
    }

    cout << "    Saving..." << endl;
    save_output_nifti(fout, "zscore", nii_zscore, false);

    cout << "Finished." << endl;
    return 0;
}
