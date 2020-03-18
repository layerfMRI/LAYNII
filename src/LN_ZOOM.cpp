

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_ZOOM: Reduces the nii file to only include the rectangular space \n"
    "         that is different from zero.\n"
    "\n"
    "Usage: \n"
    "    LN_ZOOM -mask mask.nii -input file_to_be_zoomed.nii \n"
    "\n"
    "Options:\n"
    "\n"
    "    -help  : Show this help.\n"
    "    -input : Nifti (.nii) file that should be zoomed (e.g. with \n"
    "             multiple time points).\n"
    "    -mask  : Nifti (.nii) file that determines the region of interest\n"
    "             (e.g. the layer mask with one time point).\n"
    "\n"
    "Note: Written for Insub. \n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    char* fin_1 = NULL, *fin_2 = NULL;
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
            fin_1 = argv[ac];
        } else if (!strcmp(argv[ac], "-mask")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -mask\n");
                return 1;
            }
            fin_2 = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }
    if (!fin_1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!fin_2) {
        fprintf(stderr, "** missing option '-mask'\n");
        return 1;
    }

    // Read input dataset
    nifti_image* nii1 = nifti_image_read(fin_1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_1);
        return 2;
    }
    nifti_image* nii2 = nifti_image_read(fin_2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin_2);
        return 2;
    }

    log_welcome("LN_ZOOM");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);  // mask

    // Get dimensions of input
    const int size_z = nii1->nz;
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_time = nii1->nt;
    const int nx = nii1->nx;
    const int nxy = nii1->nx * nii1->ny;
    const int nxyz = nii1->nx * nii1->ny * nii1->nz;

    // ========================================================================
    // Fix datatype issues
    nifti_image* nim_file_1 = copy_nifti_as_float32(nii1);
    float* nii1_data = static_cast<float*>(nim_file_1->data);
    nifti_image* nim_file_2 = copy_nifti_as_float32(nii2);
    float* nii2_data = static_cast<float*>(nim_file_2->data);

    // ========================================================================
    // Initialize min-max inversely to increas-decrease in the loop
    int min_z = 10000, max_z = 0;
    int min_x = 10000, max_x = 0;
    int min_y = 10000, max_y = 0;
    float mask_val = 0;
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                mask_val = *(nii2_data + nxy * iz + nx * iy + ix);
                if (mask_val > 0 && ix > max_x) max_x = ix;
                if (mask_val > 0 && iy > max_y) max_y = iy;
                if (mask_val > 0 && iz > max_z) max_z = iz;
                if (mask_val > 0 && ix < min_x) min_x = ix;
                if (mask_val > 0 && iy < min_y) min_y = iy;
                if (mask_val > 0 && iz < min_z) min_z = iz;
            }
        }
    }
    cout << "  x range is " << min_x << "-" << max_x << endl;
    cout << "  y range is " << min_y << "-" << max_y << endl;
    cout << "  z range is " << min_z << "-" << max_z << endl;
    const int new_size_z = max_z - min_z + 1;
    const int new_size_x = max_x - min_x + 1;
    const int new_size_y = max_y - min_y + 1;
    const int new_nr_voxels = new_size_x * new_size_y * new_size_z;

    // ========================================================================
    // Handle new (zoomed) nifti
    nifti_image* nii_new = nifti_copy_nim_info(nim_file_1);
    nii_new->nx = new_size_x;
    nii_new->ny = new_size_y;
    nii_new->nz = new_size_z;
    nii_new->nvox = new_nr_voxels;
    nii_new->datatype = NIFTI_TYPE_FLOAT32;
    nii_new->nbyper = sizeof(float);
    nii_new->data = calloc(nii_new->nvox, nii_new->nbyper);
    float* nii_new_data = static_cast<float*>(nii_new->data);

    // ========================================================================
    // Copy voxel from the bigger input nifti image
    const int nx_2 = nii_new->nx;
    const int nxy_2 = nii_new->nx * nii_new->ny;
    const int nxyz_2 = nii_new->nx * nii_new->ny * nii_new->nz;
    for (int it = 0; it < size_time; ++it) {
        for (int iz = min_z; iz <= max_z; ++iz) {
            for (int iy = min_y; iy <= max_y; ++iy) {
                for (int ix = min_x; ix <= max_x; ++ix) {
                    int voxel_i = nxyz * it + nxy * iz + nx * iy + ix;
                    int voxel_j = nxyz_2 * it
                                  + nxy_2 * (iz - min_z)
                                  + nx_2 * (iy - min_y)
                                  + (ix - min_x);
                    *(nii_new_data + voxel_j) = *(nii1_data + voxel_i);
                }
            }
        }
    }
    save_output_nifti(fin_1, "zoomed", nii_new, true);

    cout << "Finished!" << endl;
    return 0;
}
