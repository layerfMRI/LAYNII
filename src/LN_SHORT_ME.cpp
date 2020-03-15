

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_SHORT_ME: Convert nifti datatype to SHORT16. \n"
    "\n"
    "Usage:\n"
    "    LN_SHORT_ME -input data_file.nii \n"
    "    LN_SHORT_ME -input data_file.nii -output output_filename.nii \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Dataset that should be shorted data.\n"
    "    -output : (Optional) Output filename. If this parameter is not set \n"
    "              the original file will be overwritten.\n"
    "\n");
    return 0;
}

int main(int argc, char *argv[]) {
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

    log_welcome("LN_SHORT_ME");
    log_nifti_descriptives(nii);

    // Cast input data to float
    nifti_image *nii_new = copy_nifti_as_int16(nii);
    save_output_nifti(fout, "", nii_new, true);

    cout << "  Finished." << endl;
    return 0;
}
