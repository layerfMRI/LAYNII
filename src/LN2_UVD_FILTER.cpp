#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>
#include <vector>
#include <algorithm>

int show_help(void) {
    printf(
    "LN2_UVD_FILTER: Median filter using flat coordinates (UV) and depth (D).\n"
    "                       Passes a cylinder through UVD coordinates.\n"
    "\n"
    "Usage:\n"
    "    LN2_UVD_FILTER -values activation.nii -coord_uv uv_coord.nii -coord_d layers_equidist.nii -radius 3 -height 0.25\n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -values    : Nifti image with values that will be filtered.\n"
    "                 For example an activation map or anatomical T1w images.\n"
    "    -coord_uv  : A 4D nifti file that contains 2D (UV) coordinates.\n"
    "                 For example LN2_MULTILATERATE output named 'UV_coords'.\n"
    "    -coord_d   : A 3D nifti file that contains cortical depth measurements or layers.\n"
    "                 For example either LN2_LAYERS output named 'metric'.\n"
    "    -radius    : Radius of cylinder that will be passed over UV coordinates.\n"
    "                 In units of UV coordinates, which often are in mm.\n"
    "    -height : height/height of cylinder that will be passed over D (depth)\n"
    "                 coordinates. In units of normalized depth metric, which are often in\n"
    "                 0-1 range.\n"
    "    -output    : (Optional) Output basename for all outputs.\n"
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

    log_welcome("LN2_UVD_FILTER");
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
    nifti_image* nii_output = copy_nifti_as_float32(nii_input);
    float* nii_output_data = static_cast<float*>(nii_output->data);

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
    // Visit each voxel to check their coordinate
    float half_height = height / 2;
    float radius_sqr = radius * radius;
    for (int i = 0; i != nr_voi; ++i) {
        cout << "\r    " << i * 100 / nr_voi << " %" << flush;
        vector <float> temp_vec;

        // --------------------------------------------------------------------
        // Cylinder windowing in UVD space
        for (int j = 0; j != nr_voi; ++j) {
            // Compute distances relative to reference UVD
            if (abs(vec_d[i] - vec_d[j]) < half_height) {  // Check height
                float dist_uv = (vec_u[i] - vec_u[j])*(vec_u[i] - vec_u[j])
                    + (vec_v[i] - vec_v[j])*(vec_v[i] - vec_v[j]);
                if (dist_uv < radius_sqr) {  // Check Euclidean distance
                    temp_vec.push_back(vec_val[j]);
                }
            }
        }

        // --------------------------------------------------------------------
        // Find median
        int n = temp_vec.size();
        float m;
        if (n % 2 == 0) {  // even
            std::nth_element(temp_vec.begin(),
                             temp_vec.begin() + n / 2,
                             temp_vec.end());

            std::nth_element(temp_vec.begin(),
                             temp_vec.begin() + (n - 1) / 2,
                             temp_vec.end());

            m = (temp_vec[n / 2] + temp_vec[(n - 1) / 2]) / 2.0;

        } else {  // odd
            std::nth_element(temp_vec.begin(),
                             temp_vec.begin() + n / 2,
                             temp_vec.end());

            m = temp_vec[n / 2];
        }
        // --------------------------------------------------------------------

        // Write median inside nifti
        *(nii_output_data + vec_voi_id[i]) = m;
    }
    cout << endl;

    save_output_nifti(fout, "UVD_median_filter", nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
