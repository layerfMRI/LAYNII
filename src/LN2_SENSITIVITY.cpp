#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_SENSITIVITY: Compute a voxel-wise measure of functional sensitivity\n"
    "                 from a 4D nifti containing fMRI responses to multiple\n"
    "                 (at least 2) tasks (e.g. betas, percent signal changes...).\n"
    "                 This method provides a measure of how strongly a voxel responds\n" 
    "                 to different tasks (overall responsiveness of a voxel).\n"
    "\n"
    "Usage:\n"
    "    LN2_SENSITIVITY -input input.nii\n"
    "    ../LN2_SENSITIVITY -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : 4D matrix of dimensions (X, Y, Z, N) where (X, Y, Z) are\n"
    "              the spatial dimensions of the brain volume and N is the number\n"
    "              of task conditions (e.g. fMRI task responses)\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n"
    "Citation:\n"
    "    - Pizzuti, A., Huber, L., Gulban, O.F, Benitez-Andonegui A., Peters, J., Goebel R.,\n"
    "      (2023). Imaging the columnar functional organization of human area MT+ to \n"
    "      axis-of-motion stimuli using VASO at 7 Tesla. Cerebral Cortex.\n"
    "      <https://doi.org/10.1093/cercor/bhad151>\n"
    "\n"
    "NOTES: \n" 
    "    Sensitivity is based on the magnitude (ln2norm) of a voxel's response profile.\n"
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

    log_welcome("LN2_SENSITIVITY");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint64_t size_x = nii1->nx;
    const uint64_t size_y = nii1->ny;
    const uint64_t size_z = nii1->nz;
    const uint64_t size_time = nii1->nt;

    const uint64_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues and prepare 3D Nifti output
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Allocate new nifti for 3D outputs
    nifti_image* nii_output = nifti_copy_nim_info(nii_input);
    nii_output->nt = 1;
    nii_output->nvox = nr_voxels;
    nii_output->datatype = NIFTI_TYPE_FLOAT32;
    nii_output->nbyper = sizeof(float);
    nii_output->data = calloc(nii_output->nvox, nii_output->nbyper);
    float* nii_output_data = static_cast<float*>(nii_output->data);

    // // ========================================================================
    cout << "\n  Calculating sensitivity..." << endl;
    // // ========================================================================
    // Compute L2 norm (Euclidean norm) across the time dimension
    for (uint64_t i = 0; i != nr_voxels; ++i) {  // Loop across voxels
        float sum_sq = 0;
        for (uint64_t t = 0; t != size_time; ++t) {  // Loop across time points 
            float val = *(nii_input_data + i + nr_voxels*t);
            if (val < 0.0) {
                val = 0;
            }
            sum_sq += val * val;
        }
        // Store the computed L2 norm in the 3D output
        *(nii_output_data + i) = sqrt(sum_sq);
    }

    cout << "    Saving..." << endl;
    save_output_nifti(fout, "sensitivity", nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
