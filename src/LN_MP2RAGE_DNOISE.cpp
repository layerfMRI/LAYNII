
#include "./common.h"
#include "./utils.h"

int show_help(void) {
    printf(
    "LN_MP2RAGE_DNOISE : Denoising MP2RAGE data.\n"
    "\n"
    "    This program removes some of the background noise in MP2RAGE, \n"
    "    UNI images to make themn look like MPRAGE images. This is done \n"
    "    without needing to have the phase information. This is following \n"
    "    the paper O’Brien KR et al. (2014) Robust T1-Weighted Structural \n"
    "    Brain Imaging and Morphometry at 7T Using MP2RAGE. PLoS ONE 9(6): \n"
    "    e99676. <doi:10.1371/journal.pone.0099676> \n"
    "\n"
    "Usage:\n"
    "    LN_MP2RAGE_DNOISE -INV1 INV1.nii -INV2 INV2.nii -UNI UNI.nii -beta 0.2\n"
    "\n"
    "Options\n"
    "    -help       : Show this help.\n"
    "    -INV1       : Nifti (.nii) file of the first inversion time.\n"
    "    -INV2       : Nifti (.nii) file of the second inversion time.\n"
    "    -UNI        : Nifti (.nii) of MP2RAGE UNI. Expecting SIEMENS \n"
    "                  values between 0-4095. \n"
    "    -beta value : Regularization term. Default is 0.2.\n"
    "    -output     : (Optional) Custom output name. \n"
    "\n"
    "Note: This program supports INT16, INT32 and FLOAT32. \n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    float SIEMENS_f = 4095.0;
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
    int nx = nii1->nx;
    int nxy = nii1->nx * nii1->ny;
    int nxyz = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix datatype issues

    nifti_image* nii_inv1 = recreate_nii_with_float_datatype(nii1);
    float* nii_inv1_data = static_cast<float*>(nii_inv1->data);

    nifti_image* nii_inv2 = recreate_nii_with_float_datatype(nii2);
    float* nii_inv2_data = static_cast<float*>(nii_inv2->data);

    nifti_image* nii_uni = recreate_nii_with_float_datatype(nii3);
    float* nii_uni_data = static_cast<float*>(nii_uni->data);

    // ========================================================================

    /////////////////////////////////////
    // Make allocating necessary files //
    /////////////////////////////////////
    nifti_image* phaseerror = nifti_copy_nim_info(nii_inv1);
    nifti_image* dddenoised = nifti_copy_nim_info(nii_inv1);
    phaseerror->datatype = NIFTI_TYPE_FLOAT32;
    dddenoised->datatype = NIFTI_TYPE_FLOAT32;
    phaseerror->nbyper = sizeof(float);
    dddenoised->nbyper = sizeof(float);
    phaseerror->data = calloc(phaseerror->nvox, phaseerror->nbyper);
    dddenoised->data = calloc(dddenoised->nvox, dddenoised->nbyper);
    float* phaseerror_data = static_cast<float*>(phaseerror->data);
    float* dddenoised_data = static_cast<float*>(dddenoised->data);

    nifti_image* uni1  = nifti_copy_nim_info(nii_inv1);
    nifti_image* uni2  = nifti_copy_nim_info(nii_inv1);
    uni1->datatype = NIFTI_TYPE_FLOAT32;
    uni2->datatype = NIFTI_TYPE_FLOAT32;
    uni1->nbyper = sizeof(float);
    uni2->nbyper = sizeof(float);
    uni1->data = calloc(uni1->nvox, uni1->nbyper);
    uni2->data = calloc(uni2->nvox, uni2->nbyper);
    // float* uni1_data = static_cast<float*>(uni1->data);
    // float* uni2_data = static_cast<float*>(uni2->data);

    ///////////////////////////////////////
    // Big calculation across all voxels //
    ///////////////////////////////////////
    // float sign_ = 0;

    beta = beta * SIEMENS_f;

    float inv2val = 0;
    // float inv1val = 0;
    float unival = 0;
    float wrong_unival = 0;

    float uni1val_calc = 0;
    float uni2val_calc = 0;
    // float denoised_wrong = 0;

    /////////////////////////////////////////////////////////
    // Scaling UNI to range of -0.5 to 0.5 as in the paper //
    /////////////////////////////////////////////////////////
    FOR_EACH_VOXEL_TZYX
        unival = (*(nii_uni_data + VOXEL_ID) - SIEMENS_f * 0.5) / SIEMENS_f;
        // inv1val = *(nii_inv1_data + VOXEL_ID);
        inv2val = *(nii_inv2_data + VOXEL_ID);
        wrong_unival = *(nii_inv1_data + VOXEL_ID) * *(nii_inv2_data + VOXEL_ID) / (*(nii_inv1_data + VOXEL_ID)* *(nii_inv1_data + VOXEL_ID) + *(nii_inv2_data + VOXEL_ID)* *(nii_inv2_data + VOXEL_ID));

        // sign_ = unival;
        // *(nii_uni_data + VOXEL_ID) / *(phaseerror_data + VOXEL_ID);
        // if (sign_ <= 0) {
        //     *(nii_inv1_data + VOXEL_ID) = -1 * *(nii_inv1_data + VOXEL_ID);
        // }

        // denoised_wrong = (*(nii_inv1_data + VOXEL_ID) * *(nii_inv2_data + VOXEL_ID) - beta) / (*(nii_inv1_data + VOXEL_ID) * *(nii_inv1_data + VOXEL_ID) + *(nii_inv2_data + VOXEL_ID) * *(nii_inv2_data + VOXEL_ID) + 2. * beta);
        // denoised_wrong = (denoised_wrong +0.5) * SIEMENS_f;
        *(phaseerror_data + VOXEL_ID) = wrong_unival;

        uni1val_calc = inv2val * (1 / (2 * unival) + sqrt(1 / (4 * unival * unival) - 1));
        uni2val_calc = inv2val * (1 / (2 * unival) - sqrt(1 / (4 * unival * unival) - 1));

        if (unival > 0) {
            uni1val_calc = uni2val_calc;
        }

        // if (!(uni1val_calc > SIEMENS_f || uni1val_calc < SIEMENS_f)) uni1val_calc = inv1val;

        // *(uni1_data + VOXEL_ID) = uni1val_calc;
        // *(uni2_data + VOXEL_ID) = uni2val_calc;
        // *(phaseerror_data + VOXEL_ID) = unival;

        *(dddenoised_data + VOXEL_ID) = ((uni1val_calc * inv2val - beta) / (uni1val_calc * uni1val_calc + inv2val * inv2val + 2. * beta) + 0.5) * SIEMENS_f;
    END_FOR_EACH_VOXEL_TZYX

    dddenoised->scl_slope = nii_uni->scl_slope;

    if (nii_uni->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " ## Why would you do such a thing????    ## " << endl;
        cout << " #####   WARNING   WANRING   WANRING  ##### " << endl;
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
        string outfilename = prefix+filename;
        log_output(outfilename.c_str());
        const char* fout_1 = outfilename.c_str();
        if (nifti_set_filenames(dddenoised, fout_1, 1, 1)) {
            return 1;
        }
    }

    nifti_image_write(dddenoised);

    const char* fout_2 = "Border_enhance.nii";
    if (nifti_set_filenames(phaseerror, fout_2, 1, 1)) {
        return 1;
    }
    nifti_image_write(phaseerror);

    cout << "  Finished." << endl;
    return 0;
}
