#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>

int show_help(void) {
    printf(
    "LN2_UVD_LSTSQR: Least squares fit within UVD cylinders.\n"
    "\n"
    "Usage:\n"
    "    LN2_UVD_LSTSQR -values activation.nii -coord_uv uv_coord.nii -coord_d layers_equidist.nii -radius 3 -height 0.25\n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -values   : Nifti image with values that will be filtered.\n"
    "                For example an activation map or anatomical T1w images.\n"
    "    -coord_uv : A 4D nifti file that contains 2D (UV) coordinates.\n"
    "                For example LN2_MULTILATERATE output named 'UV_coords'.\n"
    "    -coord_d  : A 3D nifti file that contains cortical depth measurements or layers.\n"
    "                For example either LN2_LAYERS output named 'metric'.\n"
    "    -radius   : Radius of cylinder that will be passed over UV coordinates.\n"
    "                In units of UV coordinates, which often are in mm.\n"
    "    -height   : Height of cylinder that will be passed over D (depth)\n"
    "                coordinates. In units of normalized depth metric, which\n"
    "                are often in 0-1 range. The cylinder is centered around each voxel\n"
    "                therefore, to ensure all depth is included, this parameter should be\n"
    "                set to 2 when normalized depth metrics are being used.\n"
    "    -output   : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL;
    int ac;
    float radius = 3, height = 0.25;

    // Process user options
    if (argc < 2) return show_help();

    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-help", 5)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-values")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -values\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-coord_uv")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_uv\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-coord_d")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_d\n");
                return 1;
            }
            fin3 = argv[ac];
        } else if (!strcmp(argv[ac], "-radius")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -radius\n");
                return 1;
            }
            radius = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-height")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -height\n");
                return 1;
            }
            height = atof(argv[ac]);
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
        fprintf(stderr, "** missing option '-values'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-coords_uv'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-coords_d'\n");
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
    nii3 = nifti_image_read(fin3, 1);
    if (!nii3) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin3);
        return 2;
    }

    log_welcome("LN2_UVD_LSTSQR");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);

    // Get dimensions of input
    const int nr_voxels = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* coords_uv = copy_nifti_as_float32(nii2);
    float* coords_uv_data = static_cast<float*>(coords_uv->data);
    nifti_image* coords_d = copy_nifti_as_float32(nii3);
    float* coords_d_data = static_cast<float*>(coords_d->data);

    // ========================================================================
    // Prepare outputs
    nifti_image* nii_slope = copy_nifti_as_float32(nii_input);
    float* nii_slope_data = static_cast<float*>(nii_slope->data);
    nifti_image* nii_intercept = copy_nifti_as_float32(nii_input);
    float* nii_intercept_data = static_cast<float*>(nii_intercept->data);
    nifti_image* nii_samples = copy_nifti_as_float32(nii_input);
    float* nii_samples_data = static_cast<float*>(nii_samples->data);
    nifti_image* nii_residuals = copy_nifti_as_float32(nii_input);
    float* nii_residual_data = static_cast<float*>(nii_residuals->data);

    // Zero output niftis
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_slope_data + i) = 0;
        *(nii_intercept_data + i) = 0;
        *(nii_samples_data + i) = 0;
        *(nii_residual_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    int nr_voi = 0;  // Voxels of interest
    vector <int> vec_voi_id;
    vector <float> vec_u, vec_v, vec_d, vec_val;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(coords_uv_data + i) != 0){
            vec_voi_id.push_back(i);
            vec_u.push_back(*(coords_uv_data + nr_voxels*0 + i));
            vec_v.push_back(*(coords_uv_data + nr_voxels*1 + i));
            vec_d.push_back(*(coords_d_data + i));
            vec_val.push_back(*(nii_input_data + i));
            nr_voi += 1;
        }
    }

    // ========================================================================
    // Visit each voxel
    // ========================================================================
    cout << "  Fitting..." << endl;

    float half_height = height / 2;
    float radius_sqr = radius * radius;
    for (int i = 0; i != nr_voi; ++i) {
        // --------------------------------------------------------------------
        // Cylinder windowing in UVD space
        // --------------------------------------------------------------------
        vector <float> vec_y;
        for (int j = 0; j != nr_voi; ++j) {
            // Compute distances relative to reference UVD
            if (abs(vec_d[i] - vec_d[j]) < half_height) {  // Check height
                float dist_uv = (vec_u[i] - vec_u[j])*(vec_u[i] - vec_u[j])
                    + (vec_v[i] - vec_v[j])*(vec_v[i] - vec_v[j]);
                if (dist_uv < radius_sqr) {  // Check Euclidean distance
                    vec_y.push_back(vec_val[j]);
                }
            }
        }

        int n = vec_y.size();
        if (n > 1) {
            // ----------------------------------------------------------------
            // TODO: Sort vector by depth
            // ----------------------------------------------------------------

            // ----------------------------------------------------------------
            // Create design matrix
            // ----------------------------------------------------------------
            // TODO(Faruk): All ones flat profile for now. I need to find a way
            // to allow users to choose this.
            vector <float> vec_x(n, 1.0);

            // ----------------------------------------------------------------
            // Compute the least-squares solution to a linear matrix equation
            // ----------------------------------------------------------------
            float vec_x_sum = std::accumulate(vec_x.begin(), vec_x.end(), 0);
            float vec_y_sum = std::accumulate(vec_y.begin(), vec_y.end(), 0);
            float vec_x_avg = vec_x_sum / n;
            float vec_y_avg = vec_y_sum / n;


            float var_x = 0;
            float cov_xy = 0;
            float slope = 0, intercept_y = 0, residual = 0;

            for (int j = 0; j != n; ++j) {
                float temp = vec_x[j] - vec_x_avg;
                var_x += temp * temp;
                cov_xy += temp * (vec_y[j] - vec_y_avg);
            }

            if (var_x != 0) {
                slope = cov_xy / var_x;
                intercept_y = vec_y_avg - (slope * vec_x_avg);
            } else {
                slope = 0;  // avoid nans
                intercept_y = vec_y_avg;
            }

            // Compute residuals
            for (int j = 0; j != n; ++j) {
                float y_hat = vec_x[j] * slope + intercept_y;
                residual += (vec_y[j] - y_hat) * (vec_y[j] - y_hat);
            }
            residual /= n;

            // ----------------------------------------------------------------
            // Write median inside nifti
            // ----------------------------------------------------------------
            cout << "\r    " << i << "/" << nr_voi << " | n = " << n << " "
            << " | y_avg = " << vec_y_avg << flush;

            *(nii_slope_data + vec_voi_id[i]) = slope;
            *(nii_intercept_data + vec_voi_id[i]) = intercept_y;
            *(nii_samples_data + vec_voi_id[i]) = n;
            *(nii_residual_data + vec_voi_id[i]) = residual;
        }
    }
    cout << endl;

    save_output_nifti(fout, "UVD_lstsqr_slope", nii_slope, true);
    save_output_nifti(fout, "UVD_lstsqr_intercept", nii_intercept, true);
    save_output_nifti(fout, "UVD_lstsqr_samples", nii_samples, true);
    save_output_nifti(fout, "UVD_lstsqr_residuals", nii_residuals, true);

    cout << "\n  Finished." << endl;
    return 0;
}
