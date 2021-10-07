#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>
#include <vector>
#include <algorithm>

int show_help(void) {
    printf(
    "LN2_UVD_MEDIAN_FILTER: Median filter using flat coordinates (UV) and depth (D).\n"
    "                       Passes a cylinder through UVD coordinates.\n"
    "\n"
    "Usage:\n"
    "    LN2_UVD_MEDIAN_FILTER -values activation.nii -coord_uv uv_coord.nii -coord_d layers_equidist.nii -radius 3 -height 0.25\n"
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
    "    -height    : Height of cylinder that will be passed over D coordinates.\n"
    "                 In units of normalized depth metric, which are often in\n"
    "                 0-1 range.\n"
    "    -output    : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL;
    int ac;
    float radius = 0.2, height = 0.5;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
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

    log_welcome("LN2_UVD_MEDIAN_FILTER");
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
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(coords_uv_data + i) != 0){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int* voi_id;
    voi_id = (int*) malloc(nr_voi*sizeof(int));

    // Fill in indices to be able to remap from subset to full set of voxels
    int ii = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(coords_uv_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Visit each voxel to check their coordinate
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);

        float ref_u = *(coords_uv_data + nr_voxels*0 + i);
        float ref_v = *(coords_uv_data + nr_voxels*1 + i);
        float ref_d = *(coords_d_data + i);

        // Visit every voxel to check if they are within the cylinder
        vector <float> temp_vec;

        for (int jj = 0; jj != nr_voi; ++jj) {
            int j = *(voi_id + jj);

            float u = *(coords_uv_data + nr_voxels*0 + j);
            float v = *(coords_uv_data + nr_voxels*1 + j);
            float d = *(coords_d_data + j);
            float val = *(nii_input_data + j);

            // Compute distances relative to reference UVD
            float dist_d = abs(ref_d - d);
            if (dist_d < height) {  // first do the easy 1D check
                float dist_uv = sqrt(pow(ref_u - u, 2) + pow(ref_v - v, 2));
                if (dist_uv < radius) {  // second do 2D circle check
                    temp_vec.push_back(val);
                }
            }
        }

        // Find median
        auto m = temp_vec.begin() + temp_vec.size()/2;
        std::nth_element(temp_vec.begin(), m, temp_vec.end());
        float median = temp_vec[temp_vec.size()/2];
        // std::cout << "The median is " << median << '\n';

        // Write inside nifti
        *(nii_output_data + i) = median;
    }

    save_output_nifti(fout, "UVD_median_filter", nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
