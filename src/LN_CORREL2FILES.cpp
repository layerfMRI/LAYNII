


#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_CORREL2FILES: Estimate the voxel wise correlation of two timeseries.\n"
    "\n"
    "    This program is motivated by Eli Merriam comparing in hunting down \n"
    "    voxels that out of phase for VASO and BOLD. \n"
    "\n"
    "Usage:\n"
    "    LN_CORREL2FILES -file1 file1.nii -file2 file2.nii \n"
    "\n"
    "Options:\n"
    "    -help  : Show this help.\n"
    "    -file1 : First time series.\n"
    "    -file2 : Second time series with should have the same dimensions \n"
    "             as first time series.\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {
    // nifti_image* nim_input=NULL;
    char* fin_1 = NULL, *fin_2 = NULL;
    int ac;
    if (argc < 2) {  // Typing '-help' is sooo much work
       return show_help();
    }

    // Process user options
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-file1")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -file1\n");
                return 1;
            }
            fin_1 = argv[ac];  // Assign pointer, no string copy
        } else if (!strcmp(argv[ac], "-file2")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -file2\n");
                return 1;
            }
            fin_2 = argv[ac];  // Assign pointer, no string copy
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin_1) {
        fprintf(stderr, "** missing option '-file1'\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-file2'\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_1);
        return 2;
    }
    nifti_image*nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", fin_2);
        return 2;
    }

    log_welcome("LN_CORREL2FILES");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

    // Get dimensions of input
    int size_z = nii1->nz;
    int size_x = nii1->nx;
    int size_y = nii1->ny;
    int size_t = nii1->nt;
    int nx = nii1->nx;
    int nxy = nii1->nx * nii1->ny;
    int nxyz = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii1_temp = copy_nifti_as_float32(nii1);
    float* nii1_temp_data = static_cast<float*>(nii1_temp->data);
    nifti_image* nii2_temp = copy_nifti_as_float32(nii2);
    float* nii2_temp_data = static_cast<float*>(nii2_temp->data);

    // Allocate new nifti
    nifti_image* correl_file = copy_nifti_as_float32(nii1_temp);
    float* correl_file_data = static_cast<float*>(correl_file->data);

    // ========================================================================

    double vec_file1[size_t];
    double vec_file2[size_t];

    // Loop across voxels
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                for (int it = 0;  it < size_t; ++it) {
                    vec_file1[it] = *(nii1_temp_data + voxel_i);
                    vec_file2[it] = *(nii2_temp_data + voxel_i);
                }
                *(correl_file_data + voxel_i) =
                    ren_correl(vec_file1, vec_file2, size_t);
            }
        }
    }

    save_output_nifti(fin_1, "correlated", correl_file, true);

    cout << "  Finished." << endl;
    return 0;
}
