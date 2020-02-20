
#include "./common.h"
#include "./utils.h"

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
    char *input_filename = NULL, *fout = NULL;
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
                fprintf(stderr, "** missing argument for -layer_file\n");
                return 1;
            }
            input_filename = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            do_outputnaming = 1;
            cout << "  Will write output file with a different name." << endl;
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!input_filename) {
        fprintf(stderr, "** missing option '-landmarks'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image* nim_input_r = nifti_image_read(input_filename, 1);
    if (!nim_input_r) {
        fprintf(stderr, "** failed to read layer NIfTI image from '%s'\n", input_filename);
        return 2;
    }

    log_welcome("LN_FLOAT_ME");
    log_nifti_descriptives(nim_input_r);

    // Get dimsions of input
    int sizeSlice = nim_input_r->nz;
    int sizePhase = nim_input_r->nx;
    int sizeRead = nim_input_r->ny;
    int nrep = nim_input_r->nt;
    int nx = nim_input_r->nx;
    int nxy = nim_input_r->nx * nim_input_r->ny;
    int nxyz = nim_input_r->nx * nim_input_r->ny * nim_input_r->nz;

    // nim_mask->datatype = NIFTI_TYPE_FLOAT32;
    // nim_mask->nbyper = sizeof(float);
    // nim_mask->data = calloc(nim_mask->nvox, nim_mask->nbyper);

    nifti_image* nim_input = nifti_copy_nim_info(nim_input_r);
    nim_input->datatype = NIFTI_TYPE_FLOAT32;
    nim_input->nbyper = sizeof(float);
    nim_input->data = calloc(nim_input->nvox, nim_input->nbyper);
    float *nim_input_data = (float *) nim_input->data;

    //////////////////////////////////////////////////////////////
    // Fixing potential problems with different input datatypes //
    //////////////////////////////////////////////////////////////
    if (nim_input_r->datatype == NIFTI_TYPE_FLOAT32) {
        float *nim_input_r_data = (float *) nim_input_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_input_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_input_r->datatype == NIFTI_TYPE_FLOAT64) {
        double  *nim_input_r_data = (double *) nim_input_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_input_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    if (nim_input_r->datatype == NIFTI_TYPE_INT16) {
        short *nim_input_r_data = (short *) nim_input_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy < sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float)(*(nim_input_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                        //if (*(nim_input_data + nxyz * 0 + nxy * islice + nx * ix + iy) > 0) cout << * (nim_input_r_data + nxyz * 0 + nxy * islice + nx * ix + iy) << endl;
                    }
                }
            }
        }
    }
    if (nim_input_r->datatype == NIFTI_TYPE_INT32) {
        int *nim_input_r_data = (int *) nim_input_r->data;
        for (int it = 0; it < nrep; ++it) {
            for (int islice = 0; islice < sizeSlice; ++islice) {
                for (int iy = 0; iy <sizePhase; ++iy) {
                    for (int ix = 0; ix < sizeRead; ++ix) {
                        *(nim_input_data + nxyz * it + nxy * islice + nx * ix + iy) = (float) (*(nim_input_r_data + nxyz * it + nxy * islice + nx * ix + iy));
                    }
                }
            }
        }
    }
    cout << "  Writing output... " << endl;
    // Output file name
    // const char *fout_4 = "leaky_layers.nii";
    // if (nifti_set_filenames(leak_layer, fout_4, 1, 1)) {
    //     return 1;
    // }
    // nifti_image_write(leak_layer);

// const char *fout_5 = "input_file.nii";
// if (nifti_set_filenames(nim_inputf, fout_5, 1, 1)) {
//     return 1;
// }
// nifti_image_write(nim_inputf);

    if (do_outputnaming) {
        // Assign nifti_image fname/iname pair, based on output filename
        // (request to 'check' image and 'set_byte_order' here)
        if (nifti_set_filenames(nim_input, fout, 1, 1)) {
            return 1;
        }
    }
    nifti_image_write(nim_input);
    cout << "  Finished." << endl;
    return 0;
}
