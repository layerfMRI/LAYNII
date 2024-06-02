#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_REGRESS_OUT: Regress one timeseries from another, voxel-wise.\n"
    "\n"
    "Usage:\n"
    "    LN2_REGRESS_OUT -input input.nii\n"
    "    ../LN2_REGRESS_OUT -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input1 : First timeseries nifti (4D)"
    "    -input2 : Second timeseries nifti (4D)"
    "    -output : (Optional) Output basename for all outputs.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fin2 = NULL, *fout = NULL;
    int ac;
    bool mode_debug = false;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input1")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input1\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-input2")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input2\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
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

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    nii2 = nifti_image_read(fin2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }

    log_welcome("LN2_REGRESS_OUT");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;
    const uint32_t size_time = nii1->nt;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
    nifti_image* nii_input1 = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input1_data = static_cast<float*>(nii_input1->data);
    nifti_image* nii_input2 = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii2);
    float* nii_input2_data = static_cast<float*>(nii_input2->data);

    // Prepare output image
    nifti_image* nii_residual = copy_nifti_as_float32(nii_input1);
    float* nii_residual_data = static_cast<float*>(nii_residual->data);

    // Set to zero
    for (uint32_t i = 0; i != nr_voxels*size_time; ++i) {
        *(nii_residual_data + i) = 0;
    }

    nifti_image* nii_predicted = copy_nifti_as_float32(nii_residual);
    float* nii_predicted_data = static_cast<float*>(nii_predicted->data);

    // Create a 4D nifti image for point distances
    nifti_image *nii_intercept = NULL;
    float *nii_intercept_data = NULL;

    nii_intercept = nifti_copy_nim_info(nii1);
    nii_intercept->dim[0] = 4;  // For 3D nifti
    nii_intercept->dim[1] = size_x;
    nii_intercept->dim[2] = size_y;
    nii_intercept->dim[3] = size_z;
    nii_intercept->dim[4] = 1;
    nifti_update_dims_from_array(nii_intercept);
    nii_intercept->nvox = nr_voxels;
    nii_intercept->nbyper = sizeof(float);
    nii_intercept->data = calloc(nii_intercept->nvox, nii_intercept->nbyper);
    nii_intercept_data = static_cast<float*>(nii_intercept->data);

    nifti_image* nii_slope = copy_nifti_as_float32(nii_intercept);
    float* nii_slope_data = static_cast<float*>(nii_slope->data);

    nifti_image* nii_meanx = copy_nifti_as_float32(nii_intercept);
    float* nii_meanx_data = static_cast<float*>(nii_meanx->data);
    nifti_image* nii_meany = copy_nifti_as_float32(nii_intercept);
    float* nii_meany_data = static_cast<float*>(nii_meany->data);

    // ========================================================================
    cout << "  Calculating means..." << endl;
    // ========================================================================
    float n = static_cast<float>(size_time);
    for (uint32_t t = 0; t != size_time; ++t) {
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_meany_data + i) += *(nii_input1_data + i + nr_voxels*t) / n;
            *(nii_meanx_data + i) += *(nii_input2_data + i + nr_voxels*t) / n;
        }
    }

    // ========================================================================
    cout << "  Calculating slope and intercept..." << endl;
    // ========================================================================
    for (uint32_t i = 0; i != nr_voxels; ++i) {  // Loop across voxels
        float y_mean = *(nii_meany_data + i);
        float x_mean = *(nii_meanx_data + i);
        float term1 = 0, term2 = 0;

        for (uint32_t t = 0; t != size_time; ++t) {  // Loop across time points
            float y = *(nii_input1_data + i + nr_voxels*t);
            float x = *(nii_input2_data + i + nr_voxels*t);
            term1 += (x - x_mean) * (y - y_mean);
            term2 += (x - x_mean) * (x - x_mean);
        }

        *(nii_slope_data + i) = term1 / term2;
        *(nii_intercept_data + i) = y_mean - *(nii_slope_data + i) * x_mean;
    }

    save_output_nifti(fout, "slope", nii_slope, true);
    save_output_nifti(fout, "intercept", nii_intercept, true);

    // ========================================================================
    cout << "  Computing fitted timeseries..." << endl;
    // ========================================================================
    for (uint32_t i = 0; i != nr_voxels; ++i) {  // Loop across voxels
        for (uint32_t t = 0; t != size_time; ++t) {  // Loop across time points
            float slope = *(nii_slope_data + i);
            float intercept = *(nii_intercept_data + i);
            float x = *(nii_input2_data + i + nr_voxels*t);
            *(nii_predicted_data + i + nr_voxels*t) = intercept +  slope * x;
        }
    }
    save_output_nifti(fout, "fitted", nii_predicted, true);

    // ========================================================================
    cout << "  Computing residuals..." << endl;
    // ========================================================================
    for (uint32_t i = 0; i != nr_voxels; ++i) {  // Loop across voxels
        for (uint32_t t = 0; t != size_time; ++t) {  // Loop across time points
            float y_fitted = *(nii_predicted_data + i + nr_voxels*t);
            float y = *(nii_input1_data + i + nr_voxels*t);

            *(nii_residual_data + i + nr_voxels*t) = y - y_fitted;
        }
    }
    save_output_nifti(fout, "residuals", nii_residual, true);

    cout << "\n  Finished." << endl;
    return 0;
}
