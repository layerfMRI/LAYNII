#include "../dep/laynii_lib.h"
#include <sstream>
#include <vector>
#include <algorithm>

int show_help(void) {
    printf(
    "LN2_SPECIFICITY: Compute a voxel-wise measure of functional specificity\n"
    "                 from a 4D fMRI response matrix. This measure quantifies\n"
    "                 how specific (selective) a voxel responds when multiple tasks\n"
    "                 are evaluated (e.g., betas, percent signal change, t-statistics)\n"
    "                 This method is useful as a complement to the 'Winner Takes All'\n" 
    "                 approach, providing a measure of how specific is voxel's response\n"  
    "                 to the winning condition. \n"
    "                 Specificity is normalized between 0-1.\n"
    "\n"
    "Usage:\n"
    "    LN2_SPECIFICITY -input input.nii\n"
    "    ../LN2_SPECIFICITY -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : 4D matrix of dimensions (X, Y, Z, N) where\n"
    "             (X, Y, Z) are the spatial dimensions of the brain volume\n"  
    "             and N is the number of task conditions (e.g. fMRI task responses)\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n"
    "Citation:\n"
    "    - Pizzuti, A., Huber, L., Gulban, O.F, Benitez-Andonegui A., Peters, J., Goebel R.,\n"
    "      (2023). Imaging the columnar functional organization of\n"
    "      human area MT+ to axis-of-motion stimuli using VASO at 7 Tesla.\n"
    "      Cerebral Cortex. <https://doi.org/10.1093/cercor/bhad151>\n"
    "\n"
    "NOTES: \n" 
    "    Specificity is based on the cosine similarity between the voxel response profile\n"  
    "    and a reference axis representing a maximally specific response pattern.\n"
    "    By default, negative values are zeroed before computation.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fin2 = NULL, *fout = NULL;
    int ac;
    bool mode_debug = false;

    // Process user options
    if (argc < 2) return show_help();
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

    log_welcome("LN2_SPECIFICITY - Ale WIP");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;
    const uint32_t size_time = nii1->nt;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    // ========================================================================
    nifti_image* nii_input1 = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input1_data = static_cast<float*>(nii_input1->data);
   
    // Prepare output image - this should be 3D 
    nifti_image* nii_specificity = copy_nifti_as_float32(nii_input1);
    float* nii_specificity_data = static_cast<float*>(nii_specificity->data);

    for (int i = 0; i != nr_voxels*size_time; ++i) {
        *(nii_specificity_data + i) = 0;
    }

    // // ========================================================================
    cout << " Calculating specificity..." << endl;
    // // ========================================================================
    
    const float ONEPI = 3.14159265358979f;

    // Dynamically create reference vector `v` of length `size_time`
    vector<float> v(size_time, 0.0f);
    v[size_time - 1] = 1.0f;             // Last element is 1, others are 0
    const float norm_v = 1.0f;

    // Computing the max angle with lowest specificity between: 
    // reference vector: v and vector of ones: n (equally responding to all conditions)
    // dot prod (v, n) = 1 ; norm v = 1, norm n = sqrt(size_time)
    float max_angle = std::acos(1.0 / std::sqrt(size_time)) * 180.0 / ONEPI;

    // Process each voxel
    for (uint32_t i = 0; i != nr_voxels; ++i) {   // Loop across voxels
        vector<float> u(size_time, 0.0f);

        for (uint32_t t = 0; t != size_time; ++t) {  // Loop across vector components    
            float val = *(nii_input1_data + i + nr_voxels*t);
                if (val < 0.0) {
                    val = 0;
                }  
            u[t] = val;         
        }

        // Sort the values
        sort(u.begin(), u.end());

        // Compute cosine similarity
        float dot_product = 0.0, norm_u = 0.0;
        for (size_t i = 0; i < size_time; ++i) {
            dot_product += u[i] * v[i];
            norm_u += u[i] * u[i];
        }

        if (norm_u > 0) {
            float cosine = dot_product / (sqrt(norm_u) * sqrt(norm_v));
            cosine = min(1.0f, max(-1.0f, cosine)); // Clip to valid range

            // Convert to degrees
            float angle_degree = acos(cosine) * 180.0f / ONEPI;

            // Normalize and invert
            float vox_spec = 1.0f - (angle_degree / max_angle);

            // Store the computed L2 norm in the 3D output
            *(nii_specificity_data + i) = vox_spec;
        }
    }
    save_output_nifti(fout, "specificity", nii_specificity, true);

    cout << "\n  Finished." << endl;
    return 0;
}