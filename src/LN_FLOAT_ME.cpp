
#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_FLOAT_ME: Convert nifti datatype to FLOAT32. \n"
    "\n"
    "Usage:\n"
    "    LN_FLOAT_ME -input data_file.nii \n"
    "    LN_FLOAT_ME -input data_file.nii -output output_filename.nii \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Dataset that should be shorted data.\n"
    "    -output : (Optional) Output name. Overwrites existing files.\n"
    "\n");
    return 0;
}

int main(int argc, char *argv[]) {
    bool use_outpath = false;
    char *fin = NULL, *fout = NULL;
    int ac;
    if (argc < 3) return show_help();

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
            fout = fin;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
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
    nifti_image *nii = nifti_image_read(fin, 1);
    if (!nii) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_FLOAT_ME");
    log_nifti_descriptives(nii);

    // Get dimensions of input
    const int size_x = nii->nx;
    const int size_y = nii->ny;
    const int size_z = nii->nz;
    const int size_t = nii->nt;
    const int nr_voxels = size_t * size_z * size_y * size_x;

    // Cast input data to float
    nifti_image *nii_new = copy_nifti_as_float32(nii);
    float* nii_new_data = static_cast<float*>(nii_new->data);

    // Handle scaling factor effects
    float scl_slope = nii->scl_slope;
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_new_data + i) *= scl_slope;
    }
    nii_new->scl_slope = 1.;

    // Save
    if (!use_outpath) fout = fin;
    save_output_nifti(fout, "float", nii_new, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}
