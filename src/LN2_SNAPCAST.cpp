#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_SNAPCAST: Create snapshots of isometric projections over the cardinal axes.\n"
    "              It is intended to be used as a fast way to visualize subdural\n"
    "              and leptomeningeal brain surace renders.\n"
    "\n"
    "Usage:\n"
    "    LN2_SNAPCAST -input T2star.nii.gz -mask brain_mask.nii.gz\n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -input    : A 3D nifti file that contains values to project.\n"
    "                Only non zero voxels will be used.\n"
    // "    -mask     : (Optional) 3D binary nifti file that will be used to mask the input.\n"
    "    -output   : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        // } else if (!strcmp(argv[ac], "-mask")) {
        //     if (++ac >= argc) {
        //         fprintf(stderr, "** missing argument for -mask\n");
        //         return 1;
        //     }
        //     radius = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }

    if (!fin1) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    log_welcome("LN2_SNAPCAST");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;
    const uint32_t size_time = nii1->nt;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // ========================================================================
    // Compute the new dimensions of output
    // ========================================================================
    const uint32_t size_max = std::max(size_z, std::max(size_x, size_y));
    const uint32_t size_x_out = size_max;
    const uint32_t size_y_out = size_max;
    const uint32_t size_z_out = 6;  // For 6 cube faces (cardinal data axes)
    const uint32_t size_time_out = size_time;

    const uint32_t nr_voxels_out = size_z_out * size_y_out * size_x_out;

    // ========================================================================
    // Allocating new nifti for output 3D nifti
    // ========================================================================
    // Prepare the output nifti
    nifti_image* nii_output = nifti_copy_nim_info(nii_input);
    nii_output->datatype = NIFTI_TYPE_FLOAT32;
    nii_output->dim[0] = 4;  // For proper 4D nifti
    nii_output->dim[1] = size_x_out;
    nii_output->dim[2] = size_y_out;
    nii_output->dim[3] = size_z_out;
    nii_output->dim[4] = size_time_out;
    nifti_update_dims_from_array(nii_output);
    nii_output->nvox = nr_voxels_out;
    nii_output->nbyper = sizeof(float);
    nii_output->data = calloc(nii_output->nvox, nii_output->nbyper);
    nii_output->scl_slope = 1;
    float* nii_output_data = static_cast<float*>(nii_output->data);

    log_nifti_descriptives(nii_output);

    for (int i = 0; i != nr_voxels_out; ++i) {
        *(nii_output_data + i) = 0;
    }

    // ========================================================================
    // Compute rays
    // ========================================================================
    // Loop across rays on x plane
    for (int z = 0; z != size_z; ++z) {
        for (int y = 0; y != size_y; ++y) {

            // Loop across a ray
            for (int x = 0; x != size_x; ++x) {
                int i = sub2ind_3D(x, y, z, size_x, size_y);
                if ( *(nii_input_data + i) != 0 ) {

                    int k = sub2ind_3D(y, z, 0, size_max, size_max);
                    *(nii_output_data + k) = *(nii_input_data + i);
                }
            }

            // Loop across a ray in reverse
            for (int x = size_x-1; x >= 0; --x) {
                int i = sub2ind_3D(x, y, z, size_x, size_y);
                if ( *(nii_input_data + i) != 0 ) {

                    int k = sub2ind_3D(y, z, 1, size_max, size_max);
                    *(nii_output_data + k) = *(nii_input_data + i);
                }
            }

        }
    }

    // ------------------------------------------------------------------------
    // Loop across rays on y plane
    for (int z = 0; z != size_z; ++z) {
        for (int x = 0; x != size_x; ++x) {

            // Loop across a ray
            for (int y = 0; y != size_y; ++y) {
                int i = sub2ind_3D(x, y, z, size_x, size_y);
                if ( *(nii_input_data + i) != 0 ) {

                    int k = sub2ind_3D(x, z, 2, size_max, size_max);
                    *(nii_output_data + k) = *(nii_input_data + i);
                }
            }

            // Loop across a ray in reverse
            for (int y = size_y-1; y >= 0; --y) {
                int i = sub2ind_3D(x, y, z, size_x, size_y);
                if ( *(nii_input_data + i) != 0 ) {

                    int k = sub2ind_3D(x, z, 3, size_max, size_max);
                    *(nii_output_data + k) = *(nii_input_data + i);
                }
            }

        }
    }

    // ------------------------------------------------------------------------
    // Loop across rays on z plane
    for (int y = 0; y != size_y; ++y) {
        for (int x = 0; x != size_x; ++x) {

            // Loop across a ray
            for (int z = 0; z != size_z; ++z) {
                int i = sub2ind_3D(x, y, z, size_x, size_y);
                if ( *(nii_input_data + i) != 0 ) {

                    int k = sub2ind_3D(x, y, 4, size_max, size_max);
                    *(nii_output_data + k) = *(nii_input_data + i);
                }
            }

            // Loop across a ray in reverse
            for (int z = size_z-1; z >= 0; --z) {
                int i = sub2ind_3D(x, y, z, size_x, size_y);
                if ( *(nii_input_data + i) != 0 ) {

                    int k = sub2ind_3D(x, y, 5, size_max, size_max);
                    *(nii_output_data + k) = *(nii_input_data + i);
                }
            }
        }
    }

    // ========================================================================
    // Save
    cout << "  Saving output..." << endl;
    save_output_nifti(fout, "snapcast", nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
