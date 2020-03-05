

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_FLOAT_ME: Short convert dataset to FLOAT32. \n"
    "\n"
    "    This program convert any nii file in FLOATs. This is very useful \n"
    "    when you want the reduce the file size after ANTS was applied and \n"
    "    that have a small dynamic range.\n"
    "\n"
    "Usage:\n"
    "    LN_FLOAT_ME -input data_file.nii -output output_filename.nii \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Dataset that should be shorted data.\n"
    "    -output : (Optional) Output filename. If this parameter is not set \n"
    "              the original file will be overwritten.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *f_in = NULL, *f_out = NULL;
    int ac, do_outputnaming = 0;
    if (argc < 3) {  // Typing '-help' is sooo much work
        return show_help();
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            f_in = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            do_outputnaming = 1;
            cout << "  Writing output file with a different name." << endl;
            f_out = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!f_in) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii = nifti_image_read(f_in, 1);
    if (!nii) {
        fprintf(stderr, "** failed to read layer NIfTI from '%s'\n", f_in);
        return 2;
    }

    log_welcome("LN_FLOAT_ME");
    log_nifti_descriptives(nii);

    // Get dimensions of input
    const int size_z = nii->nz;
    const int size_x = nii->nx;
    const int size_y = nii->ny;
    const int size_t = nii->nt;
    const int nr_voxels = size_t * size_z * size_y * size_x;

    nifti_image* nii_new = nifti_copy_nim_info(nii);
    nii_new->datatype = NIFTI_TYPE_FLOAT32;
    nii_new->nbyper = sizeof(float);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    float* nii_new_data = static_cast<float*>(nii_new->data);

    // ========================================================================
    // Fixing potential problems with different input datatypes //
    if (nii->datatype == NIFTI_TYPE_FLOAT32) {
        float* nii_data = static_cast<float*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    }
    if (nii->datatype == NIFTI_TYPE_FLOAT64) {
        double* nii_data = static_cast<double*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    }
    if (nii->datatype == NIFTI_TYPE_INT16) {
        int16_t* nii_data = static_cast<int16_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    }
    if (nii->datatype == NIFTI_TYPE_INT32) {
        int32_t* nii_data = static_cast<int32_t*>(nii->data);
        for (int i = 0; i < nr_voxels; ++i) {
            *(nii_new_data + i) = static_cast<float>(*(nii_data + i));
        }
    }

    // ------------------------------------------------------------------------
    cout << "  Writing output... " << endl;
    if (do_outputnaming) {
        // Assign nifti_image fname/iname pair, based on output filename
        // (request to 'check' image and 'set_byte_order' here)
        if (nifti_set_filenames(nii_new, f_out, 1, 1)) {
            return 1;
        }
    }
    nifti_image_write(nii_new);
    cout << "  Finished." << endl;
    return 0;
}
