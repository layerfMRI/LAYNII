

#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_DEBUGGING: Short example of Layering.\n"
    "\n"
    "    This program demonstrates how to read a NIfTI-2 dataset.\n"
    "    Set output filenames and write a NIfTI-2 dataset, all via the\n"
    "    standard NIfTI C library.\n"
    "\n"
    "Usage:\n"
    "    LN_DEBUGGING -rim rim.nii \n"
    "\n"
    "Options:\n"
    "    -help               : Show this help.\n"
    "    -disp_float_example : Show some voxel's data.\n"
    "    -rim border         : Specify input dataset.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    nifti_image *nii1 = NULL;
    char *fin = NULL, *fout = NULL;
    int ac, disp_float_eg = 0;
    if (argc < 2) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-rim")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin = argv[ac];  // no string copy, just pointer assignment
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-rim'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin, 1);

    log_welcome("LN_DEBUGGING");
    log_nifti_descriptives(nii1);

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii2 = recreate_nii_with_float_datatype(nii1);

    // ========================================================================

    // Output file name
    const char* outfilename = "nii2.nii";
    if (nifti_set_filenames(nii2, outfilename , 1, 1)) {
        return 1;
    }
    log_output(outfilename);
    nifti_image_write(nii2);

    cout << "  Finished." << endl;
    return 0;
  }
