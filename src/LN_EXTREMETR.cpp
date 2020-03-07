

#include "../dep/laynii_lib.h"
#include <limits>

int show_help(void) {
    printf(
    "LN_EXTREMETR: To find extreme the TR time point.\n"
    "\n"
    "    This program tries to find the maximal/minimal value of a time \n"
    "    series and writes out the timepoint of that TR vor every voxel.\n"
    "\n"
    "Usage:\n"
    "    LN_EXTREMETR -input file.nii \n"
    "\n"
    "Options:\n"
    "    -help  : Show this help.\n"
    "    -input : Input time series.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin_1 = NULL;
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
            fin_1 = argv[ac];  // Assign pointer, no string copy
        }
    }
    if (!fin_1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    // Read input dataset, including data
    nifti_image* nii_in = nifti_image_read(fin_1, 1);
    if (!nii_in) {
      fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_1);
      return 2;
    }

    log_welcome("LN_EXTREMETR");
    log_nifti_descriptives(nii_in);

    // Get dimensions of input
    const int size_x = nii_in->nx;
    const int size_y = nii_in->ny;
    const int size_z = nii_in->nz;
    const int size_t = nii_in->nt;
    const int nx = size_x;
    const int nxy = size_x * size_y;
    const int nxyz = size_x * size_y * size_z;

    // ========================================================================
    // Fix datatype issues
    nifti_image* nii = copy_nifti_as_float32(nii_in);
    float* nii_data = static_cast<float*>(nii->data);

    // Allocate new nifti images
    nifti_image* nii_max = nifti_copy_nim_info(nii);
    nii_max->nt = 1;
    nii_max->datatype = NIFTI_TYPE_FLOAT32;
    nii_max->nbyper = sizeof(float);
    nii_max->nvox = nii_max->nt * nxyz;
    nii_max->data = calloc(nii_max->nvox, nii_max->nbyper);
    float* nii_max_data = static_cast<float*>(nii_max->data);

    nifti_image* nii_min = nifti_copy_nim_info(nii);
    nii_min->nt = 1;
    nii_min->datatype = NIFTI_TYPE_FLOAT32;
    nii_min->nbyper = sizeof(float);
    nii_min->nvox = nii_min->nt * nxyz;
    nii_min->data = calloc(nii_min->nvox, nii_min->nbyper);
    float* nii_min_data = static_cast<float*>(nii_min->data);

    // ========================================================================
    int TR_max = 0, TR_min = 0;

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                float max_val = 0;
                float min_val = std::numeric_limits<float>::max();
                for (int it = 0; it < size_t; ++it) {
                    int voxel_j = nxyz * it + nxy * iz + nx * iy + ix;
                    if (*(nii_data + voxel_j) > max_val) {
                        max_val = *(nii_data + voxel_j);
                        TR_max = it;
                    }
                    if (*(nii_data + voxel_j) < min_val) {
                        min_val = *(nii_data + voxel_j);
                        TR_min = it;
                    }
                }
                *(nii_min_data + voxel_i) = TR_min;
                *(nii_max_data + voxel_i) = TR_max;
           }
        }
    }
    save_output_nifti(fin_1, "MaxTR", nii_max, true);
    save_output_nifti(fin_1, "MinTR", nii_min, true);

    cout << "  Finished." << endl;
    return 0;
}
