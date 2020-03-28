

#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN_RAGRUG: Generating a file where every voxel has a specific value \n"
    "           that is indicates its position.\n"
    "\n"
    "    This program is inspired by Kendrick Kay's way of visualizing \n"
    "    which position at the surface comes from which voxel. See Fig. 4 of\n"
    "    <https://www.biorxiv.org/content/early/2018/06/03/337667>.\n"
    "\n"
    "Usage:\n"
    "    LN_RAGRUG -input any_file.nii \n"
    "\n"
    "Options:\n"
    "    -help  : Show this help.\n"
    "    -input : Nifti (.nii) file. This program will use the dimension of \n"
    "             this file to generate a Rag Rug file accordingly.\n"
    "\n");
    return 0;
}

int main(int argc, char * argv[]) {
    char *fin = NULL;
    int ac;
    if (argc < 3) return show_help();

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
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
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

    log_welcome("LN_RAGRUG");
    log_nifti_descriptives(nii_input);

    // Get dimensions of input
    int size_x = nii_input->nx;
    int size_y = nii_input->ny;
    int size_z = nii_input->nz;
    int nx = nii_input->nx;
    int nxy = nii_input->nx * nii_input->ny;
    int nxyz = nii_input->nx * nii_input->ny * nii_input->nz;

    // ========================================================================
    // Fixing potential problems with different input datatypes
    nifti_image* nii = copy_nifti_as_float32(nii_input);

    // Allocating new nifti images
    nifti_image* ragrug = nifti_copy_nim_info(nii);
    ragrug->nt = 1;
    ragrug->datatype = NIFTI_TYPE_INT32;
    ragrug->nbyper = sizeof(int32_t);
    ragrug->data = calloc(ragrug->nvox, ragrug->nbyper);
    ragrug->scl_slope = 1;
    int32_t* ragrug_data = static_cast<int32_t*>(ragrug->data);

    nifti_image* coord = nifti_copy_nim_info(nii);
    coord->datatype = NIFTI_TYPE_INT32;
    coord->dim[0] = 4;  // For proper 4D nifti
    coord->dim[1] = size_x;
    coord->dim[2] = size_y;
    coord->dim[3] = size_z;
    coord->dim[4] = 3;
    nifti_update_dims_from_array(coord);
    coord->nvox = nii_input->nvox * 3;
    coord->nbyper = sizeof(int32_t);
    coord->data = calloc(coord->nvox, coord->nbyper);
    coord->scl_slope = 1;
    int32_t* coord_data = static_cast<int32_t*>(coord->data);

    cout << "  Nr. voxels (input)  = " << nii_input->nvox << endl;
    cout << "  Nr. voxels (ragrug) = " << ragrug->nvox << endl;
    cout << "  Nr. voxels (coord)  = " << coord->nvox << endl;

    // ========================================================================
    cout << "  Filling nii with spatial values..." << endl;

    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int voxel_i = nxy * iz + nx * iy + ix;
                *(ragrug_data + voxel_i) = 0;

                *(coord_data + nxyz * 0 + voxel_i) = ix;
                *(coord_data + nxyz * 1 + voxel_i) = iy;
                *(coord_data + nxyz * 2 + voxel_i) = iz;

                if (ix%2 == 0) {
                    *(ragrug_data + voxel_i) += 4;
                }
                if (iy%2 == 0) {
                    *(ragrug_data + voxel_i) += 2;
                }
                if (iz%2 == 0) {
                    *(ragrug_data + voxel_i) += 1;
                }
            }
        }
    }
    save_output_nifti(fin, "ragrug", ragrug, true);
    save_output_nifti(fin, "coordinates", coord, true);

    cout << "  Finished." << endl;
    return 0;
}
