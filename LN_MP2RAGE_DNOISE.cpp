

#include "./laynii_lib.h"

int show_help(void) {
    printf(
    "LN_MP2RAGE_DNOISE : Denoising MP2RAGE data.\n"
    "\n"
    "    This program removes some of the background noise in MP2RAGE,\n"
    "    UNI images to make themn look like MPRAGE images. This is done\n"
    "    without the phase information. See O’Brien K.R. et al. (2014)\n"
    "    Robust T1-Weighted Structural Brain Imaging and Morphometry\n"
    "    at 7T Using MP2RAGE. PLoS ONE 9(6): e99676.\n"
    "    <DOI:10.1371/journal.pone.0099676>\n"
    "\n"
    "Usage:\n"
    "    LN_MP2RAGE_DNOISE -INV1 INV1.nii -INV2 INV2.nii -UNI UNI.nii -beta 0.2\n"
    "\n"
    "Options\n"
    "    -help       : Show this help.\n"
    "    -INV1       : Nifti (.nii) file of the first inversion time.\n"
    "    -INV2       : Nifti (.nii) file of the second inversion time.\n"
    "    -UNI        : Nifti (.nii) of MP2RAGE UNI. Expecting SIEMENS \n"
    "                  unsigned integer 12 values between 0-4095. \n"
    "    -beta value : Regularization term. Default is 0.2.\n"
    "    -output     : (Optional) Custom output name. \n"
    "\n"
    "Note: This program supports INT16, INT32 and FLOAT32. \n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    float SIEMENS_f = 4095.0;  // uint12 range 0-4095
    char* fmaski = NULL, *fout = NULL, *finfi_1 = NULL, *finfi_2 = NULL;
    char* finfi_3 = NULL;
    int ac, custom_output = 0;
    float beta = 0.2;
    if (argc < 3) return show_help();  // Typing '-help' is sooo much work

    // Process user options: 4 are valid presently
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-beta")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -betan");
                return 1;
            }
            beta = atof(argv[ac]);
            // cout << " I will do gaussian temporal smoothing " << endl;
        } else if (!strcmp(argv[ac], "-INV1")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -INV1\n");
                return 1;
            }
            finfi_1 = argv[ac];
        } else if (!strcmp(argv[ac], "-INV2")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -INV2\n");
                return 1;
            }
            finfi_2 = argv[ac];
        } else if (!strcmp(argv[ac], "-UNI")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -UNI\n");
                return 1;
            }
            finfi_3 = argv[ac];
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            custom_output = 1;
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!finfi_1) {
        fprintf(stderr, "** missing option '-INV1'\n");
        return 1;
    }
    if (!finfi_2) {
        fprintf(stderr, "** missing option '-INV2'\n");
        return 1;
    }
    if (!finfi_3) {
        fprintf(stderr, "** missing option '-UNI '\n");
        return 1;
    }

    // Read input dataset, including data
    nifti_image* nii1 = nifti_image_read(finfi_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", finfi_1);
        return 2;
    }

    nifti_image* nii2 = nifti_image_read(finfi_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", finfi_2);
        return 2;
    }

    nifti_image* nii3 = nifti_image_read(finfi_3, 1);
    if (!nii3) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", finfi_3);
        return 2;
    }

    log_welcome("LN_MP2RAGE_DNOISE");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);

    // Get dimensions of input
    int size_x = nii1->nx;  // phase
    int size_y = nii1->ny;  // read
    int size_z = nii1->nz;  // slice
    int size_t = nii1->nt;  // time
    int nr_voxels = size_t * size_z * size_y * size_x;

    // ========================================================================
    // Fix datatype issues

    nifti_image* nii_inv1 = recreate_nii_with_float_datatype(nii1);
    float* nii_inv1_data = static_cast<float*>(nii_inv1->data);

    nifti_image* nii_inv2 = recreate_nii_with_float_datatype(nii2);
    float* nii_inv2_data = static_cast<float*>(nii_inv2->data);

    nifti_image* nii_uni = recreate_nii_with_float_datatype(nii3);
    float* nii_uni_data = static_cast<float*>(nii_uni->data);

    // Allocate output nifti files
    nifti_image* dddenoised = recreate_nii_with_float_datatype(nii1);
    float* dddenoised_data = static_cast<float*>(dddenoised->data);

    nifti_image* phaseerror = recreate_nii_with_float_datatype(nii1);
    float* phaseerror_data = static_cast<float*>(phaseerror->data);

    // ========================================================================
    // Big calculation across all voxels

    int i = 0;
    float val_uni = *(nii_uni_data + i);
    float val_inv1 = *(nii_inv1_data + i);
    float val_inv2 = *(nii_inv2_data + i);
    float new_uni1, new_uni2, val_uni_wrong;

    beta = beta * SIEMENS_f;
    for (i = 0; i < nr_voxels; ++i) {
        // Skip nan or zero voxels
        if (isnan(val_uni) || val_uni == 0 || val_uni == 0.0) {
            *(phaseerror_data + i) = 0;
            *(dddenoised_data + i) = 0;
        } else {
            // Scale UNI to range of -0.5 to 0.5 (as in O’Brien et al. [2014])
            val_uni = val_uni / SIEMENS_f - 0.5;

            if (val_uni < 0) {
                new_uni1 = val_inv2 * (1. / (2. * val_uni)
                                       + sqrt(1. / pow(2 * val_uni, 2) - 1.));
            } else {
                new_uni1 = val_inv2 * (1. / (2. * val_uni)
                                       - sqrt(1. / pow(2 * val_uni, 2) - 1.));
            }

            // Eq. 2 in O’Brien et al. [2014].
            new_uni2 = (new_uni1 * val_inv2 - beta)
                       / ((pow(new_uni1, 2) + pow(val_inv2, 2) + 2. * beta));

            // Scale back the value range
            *(dddenoised_data + i) =  (new_uni2 + 0.5) * SIEMENS_f;

            // ----------------------------------------------------------------
            // Border enhance
            val_uni_wrong = val_inv1 * val_inv2
                            / (pow(val_inv1, 2) + pow(val_inv2, 2));

            *(phaseerror_data + i) = val_uni_wrong;
            // ----------------------------------------------------------------
        }
    }

    dddenoised->scl_slope = nii_uni->scl_slope;

    if (nii_uni->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WARNING   WARNING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " ##    Why would you do such a thing?    ## " << endl;
        cout << " #####   WARNING   WARNING   WARNING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

    if (custom_output == 1) {
        string outfilename = (string) (fout);
        log_output(outfilename.c_str());
        const char* fout_1 = outfilename.c_str();
        if (nifti_set_filenames(dddenoised, fout_1, 1, 1)) {
            return 1;
        }
    } else {
        string prefix = "denoised_";
        string filename = (string) (finfi_3);
        string outfilename = prefix + filename;
        log_output(outfilename.c_str());
        if (nifti_set_filenames(dddenoised, outfilename.c_str(), 1, 1)) {
            return 1;
        }
    }
    nifti_image_write(dddenoised);

    // const char* fout_2 = "border_enhance";
    string prefix = "border_enhance_";
    string filename = (string) (finfi_3);
    string outfilename = prefix + filename;
    if (nifti_set_filenames(phaseerror, outfilename.c_str(), 1, 1)) {
        return 1;
    }
    nifti_image_write(phaseerror);

    cout << "  Finished." << endl;
    return 0;
}
