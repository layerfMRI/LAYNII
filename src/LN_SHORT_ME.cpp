

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_SHORT_ME: Convert nifti datatype to SHORT. \n"
    "   In order not to loose too much depth resolution for \n"
    "   floating point values, the nii-header slope is scaled by 1000 \n"
    "\n"
    "Usage:\n"
    "    LN_SHORT_ME -input data_file.nii \n"
    "    LN_SHORT_ME -input data_file.nii -output output_filename.nii \n"
    "    ../LN_SHORT_ME -input lo_VASO_act.nii -output short.nii\n" 
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Dataset that should be shorted data.\n"
    "    -output : (Optional) Output filename, including .nii or\n"
    "              .nii.gz, and path if needed. Overwrites existing files.\n"    
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

    log_welcome("LN_SHORT_ME");
    log_nifti_descriptives(nii);

    // Cast input data to short (int16)
    nifti_image *nii_new = copy_nifti_as_float16(nii);
    save_output_nifti(fout, "short", nii_new, true, use_outpath);

    cout << "  Finished." << endl;
    return 0;
}
