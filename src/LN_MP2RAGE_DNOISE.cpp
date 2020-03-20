
#include "../dep/laynii_lib.h"

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
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    float SIEMENS_f = 4095.0;  // uint12 range 0-4095
    char *fout = NULL, *fin1 = NULL, *fin2 = NULL, *fin3 = NULL;
    int ac;
    int have_output = 0;
    float beta = 0.2;
    if (argc < 3) return show_help();

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-beta")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -beta");
                return 1;
            }
            beta = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-INV1")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -INV1\n");
                return 1;
            }
            fin1 = argv[ac];
        } else if (!strcmp(argv[ac], "-INV2")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -INV2\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-UNI")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -UNI\n");
                return 1;
            }
            fin3 = argv[ac];
            fout = fin3;
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            have_output = 1; 
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-INV1'\n");
        return 1;
    }
    if (!fin2) {
        fprintf(stderr, "** missing option '-INV2'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-UNI '\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }
    nifti_image* nii3 = nifti_image_read(fin3, 1);
    if (!nii3) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin3);
        return 2;
    }

    log_welcome("LN_MP2RAGE_DNOISE");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;
    const int nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii_inv1 = copy_nifti_as_float32(nii1);
    float* nii_inv1_data = static_cast<float*>(nii_inv1->data);
    nifti_image* nii_inv2 = copy_nifti_as_float32(nii2);
    float* nii_inv2_data = static_cast<float*>(nii_inv2->data);
    nifti_image* nii_uni = copy_nifti_as_float32(nii3);
    float* nii_uni_data = static_cast<float*>(nii_uni->data);

    // Allocate output nifti files
    nifti_image* nii_denoised = copy_nifti_as_float32(nii1);
    float* nii_denoised_data = static_cast<float*>(nii_denoised->data);
    nifti_image* nii_phaseerr = copy_nifti_as_float32(nii1);
    float* nii_phaseerr_data = static_cast<float*>(nii_phaseerr->data);

    // ========================================================================
    // Big calculation across all voxels
    beta = beta * SIEMENS_f;
    for (int i = 0; i != nr_voxels; ++i) {
        float val_uni = *(nii_uni_data + i);
        float val_inv1 = *(nii_inv1_data + i);
        float val_inv2 = *(nii_inv2_data + i);
        float new_uni1, new_uni2, val_uni_wrong;

        // Skip nan or zero voxels
        if (isnan(val_uni) || val_uni == 0 || val_uni == 0.0) {
            *(nii_phaseerr_data + i) = 0;
            *(nii_denoised_data + i) = 0;
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
            *(nii_denoised_data + i) =  (new_uni2 + 0.5) * SIEMENS_f;

            // ----------------------------------------------------------------
            // Border enhance
            val_uni_wrong = val_inv1 * val_inv2
                            / (pow(val_inv1, 2) + pow(val_inv2, 2));

            *(nii_phaseerr_data + i) = val_uni_wrong;
            // ----------------------------------------------------------------
        }
    }

    nii_denoised->scl_slope = nii_uni->scl_slope;

    // TODO(Faruk): This looks redundant, need to ask Renzo why it is needed.
    if (nii_uni->scl_inter != 0) {
        cout << " ########################################## " << endl;
        cout << " #####   WARNING   WARNING   WARNING  ##### " << endl;
        cout << " ## the NIFTI scale factor is asymmetric ## " << endl;
        cout << " ##    Why would you do such a thing?    ## " << endl;
        cout << " #####   WARNING   WARNING   WARNING  ##### " << endl;
        cout << " ########################################## " << endl;
    }

if (have_output == 1) {
    save_output_nifti(fout, "denoised", nii_denoised, true);
}
if (have_output == 0) {
    string prefix = "denoised_" ;
    string filename = (string) (fout) ;
    string outfilename = prefix+filename ;
    cout <<" writing "<< outfilename << endl; 
      const char  *fout_1=outfilename.c_str() ;
    if( nifti_set_filenames(nii_denoised, fout_1 , 1, 1) ) return 1;
    nifti_image_write( nii_denoised );

    //if( nifti_set_filenames(nii_denoised, outfilename , 1, 1) ) return 1;
    //nifti_image_write( nii_denoised );
}
    
    save_output_nifti(fout, "border_enhance", nii_phaseerr, true);

    cout << "  Finished." << endl;
    return 0;
}
