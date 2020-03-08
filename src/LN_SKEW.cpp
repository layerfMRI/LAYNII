

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_SKEW: Calculates skew, kurtosis, and autocorrelation of timeseries.\n"
    "         This is helpful for artifact hunting (e.g. ghosting).\n"
    "\n"
    "Usage:\n"
    "    LN_SKEW -input Nulled_intemp.nii  \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Nifti (.nii) time series.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac;
    if (argc < 2) return show_help();

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
        }
    }

    if (!fin) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset
    nifti_image * nii_input = nifti_image_read(fin, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin);
        return 2;
    }

    log_welcome("LN_SKEW");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int size_t = nii_input->nt;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // ========================================================================
    // Fix data type issues
    nifti_image* nii = copy_nifti_as_float32(nii_input);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocate new nifti
    nifti_image* nii_skew = nifti_copy_nim_info(nii);
    nii_skew->nt = 1;
    nii_skew->nvox = nii->nvox / size_t;
    nii_skew->datatype = NIFTI_TYPE_FLOAT32;
    nii_skew->nbyper = sizeof(float);
    nii_skew->data = calloc(nii_skew->nvox, nii_skew->nbyper);
    float* nii_skew_data = static_cast<float*>(nii_skew->data);

    nifti_image* nii_kurt = copy_nifti_as_float32(nii_skew);
    float* nii_kurt_data = static_cast<float*>(nii_kurt->data);

    nifti_image* nii_autocorr = copy_nifti_as_float32(nii_skew);
    float* nii_autocorr_data = static_cast<float*>(nii_autocorr->data);

    nifti_image* nii_conc = copy_nifti_as_float32(nii_skew);
    float* nii_conc_data = static_cast<float*>(nii_conc->data);

    // ========================================================================
    cout << "  Calculating skew, kurtosis, and autocorrelation..." << endl;

    double vec1[size_t];
    double vec2[size_t];

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_t; ++it) {
                    vec1[it] =
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_skew_data + voxel_i) = ren_skew(vec1, size_t);
                *(nii_kurt_data + voxel_i) = ren_kurt(vec1, size_t);
                *(nii_autocorr_data + voxel_i) = ren_autocor(vec1, size_t);
            }
        }
    }

    save_output_nifti(fin, "skew", nii_skew, true);
    save_output_nifti(fin, "kurt", nii_kurt, true);
    save_output_nifti(fin, "autocorr", nii_autocorr, true);

    // ========================================================================
    cout << "  Calculating correlation with everything..." << endl;

    for (int it = 0; it < size_t; ++it) {
        vec1[it] = 0;
        vec2[it] = 0;
    }

    // Mean time course of everything
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_t; ++it) {
                    vec1[it] +=
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i)
                                            / nxyz);
                }
            }
        }
    }

    // Voxel-wise corelation to mean of everything
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix <size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0; it < size_t; ++it)   {
                    vec2[it] =
                        static_cast<double>(*(nii_data + nxyz * it + voxel_i));
                }
                *(nii_conc_data + voxel_i) = ren_correl(vec1, vec2, size_t);
            }
        }
    }
    save_output_nifti(fin, "overall_correl", nii_conc, true);

    cout << "  Finished." << endl;
    return 0;
}
