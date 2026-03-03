#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_RECIPROCAL: Compute reciprocal (1/x) of each voxel value.\n"
    "\n"
    "Usage:\n"
    "    LN2_RECIPROCAL -input input.nii.gz\n"
    "    ../LN2_RECIPROCAL -input input.nii.gz\n"
    "\n"
    "Options:\n"
    "    -help    : Show this help.\n"
    "    -input   : Nifti image that will be used to compute reciprocals.\n"
    "               This can be a 4D nifti.\n"
    "    -thr_min : Clip smaller values to this value. Default is '1.0'.\n"
    "    -thr_max : (Optional) Clip bigger values to this value. Default is '2000'.\n"
    "    -scale   : Multiply the reciprocal with this value. Default is '1000000'.\n"
    "               Set this to `1` if you dont want any scaling."
    "    -output  : (Optional) Output basename for all outputs.\n"
    "\n"
 "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fin2 = NULL, *fout = NULL;
    int ac;
    float THR_MIN = 1.0, THR_MAX = 2000, SCL = 1000000;
    bool mode_thr_max = false, mode_scl = false;

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
        } else if (!strcmp(argv[ac], "-thr_min")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -thr_min\n");
                return 1;
            }
            THR_MIN = atof(argv[ac]);    
        } else if (!strcmp(argv[ac], "-thr_max")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -thr_max\n");
                return 1;
            }
            THR_MAX = atof(argv[ac]);    
        } else if (!strcmp(argv[ac], "-scl")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -scl\n");
                return 1;
            }
            mode_scl = true;
            SCL = atof(argv[ac]);    
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

    log_welcome("LN2_RECIPROCAL");
    log_nifti_descriptives(nii1);

    std::cout << "  Minimum threshold for clipping is: " << THR_MIN << std::endl;
    if ( mode_thr_max ) {
        std::cout << "  Maximum threshold clipping is: " << THR_MAX << std::endl;
    }
    if ( mode_scl ) {
        std::cout << "  Scaling factor is: " << SCL << std::endl;
    }

    // Get dimensions of input
    const uint64_t size_x = nii1->nx;
    const uint64_t size_y = nii1->ny;
    const uint64_t size_z = nii1->nz;
    const uint64_t size_time = nii1->nt;

    const uint64_t nxyzt = size_z * size_y * size_x * size_time;

    // ========================================================================
    // Fix input datatype issues and prepare 3D Nifti output
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // ========================================================================
    cout << "\n  Clipping small values..." << endl;
    for (uint64_t i = 0; i != nxyzt; ++i) {
        if (*(nii_input_data + i) < THR_MIN) {
            *(nii_input_data + i) = THR_MIN;
        }
    }

    // ========================================================================
    cout << "\n  Computing reciprocals..." << endl;
    for (uint64_t i = 0; i != nxyzt; ++i) {
        *(nii_input_data + i) = 1.0 / *(nii_input_data + i);
    }

    // ========================================================================
    if ( mode_thr_max ) {
        cout << "\n  Clipping large values..." << endl;
        for (uint64_t i = 0; i != nxyzt; ++i) {
            if (*(nii_input_data + i) > THR_MAX) {
                *(nii_input_data + i) = THR_MAX;
            }
        }
    }

    // ========================================================================
    cout << "\n  Scaling values..." << endl;
    for (uint64_t i = 0; i != nxyzt; ++i) {
        *(nii_input_data + i) *= SCL;
    }

    // ========================================================================
    cout << "    Saving..." << endl;
    save_output_nifti(fout, "recip", nii_input, true);

    cout << "\n  Finished." << endl;
    return 0;
}
