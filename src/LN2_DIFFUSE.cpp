
#include "../dep/laynii_lib.h"
#include <limits>

int show_help(void) {
    printf(
    "LN2_DIFFUSE: WORK IN PROGRESS... EXTREMELY EXPERIMENTAL!\n"
    "\n"
    "Usage:\n"
    "    LN2_DIFFUSE -input data.nii -nr_steps 5\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -input        : Input image.\n"
    "    -nr_steps   : Number of diffusion steps.\n"
    "    -debug        : (Optional) Save extra intermediate outputs.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    uint16_t ac, nr_steps = 5;
    bool mode_debug = false;

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
        } else if (!strcmp(argv[ac], "-nr_steps")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -nr_steps\n");
            } else {
                nr_steps = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
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

    log_welcome("LN2_DIFFUSE");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    const float dX = nii1->pixdim[1];
    const float dY = nii1->pixdim[2];
    const float dZ = nii1->pixdim[3];

    // Short diagonals
    // const float dia_xy = sqrt(dX * dX + dY * dY);
    // const float dia_xz = sqrt(dX * dX + dZ * dZ);
    // const float dia_yz = sqrt(dY * dY + dZ * dZ);
    // // Long diagonals
    // const float dia_xyz = sqrt(dX * dX + dY * dY + dZ * dZ);

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Output image
    nifti_image* nii_diffuse = copy_nifti_as_float32(nii_input);
    float* nii_diffuse_data = static_cast<float*>(nii_diffuse->data);
    // Set to zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_diffuse_data + i) = 0;
    }

    // Prepare required nifti images
    nifti_image * nii_grad = nifti_copy_nim_info(nii_diffuse);
    nii_grad->datatype = NIFTI_TYPE_FLOAT32;
    nii_grad->dim[0] = 4;  // For proper 4D nifti
    nii_grad->dim[1] = size_x;
    nii_grad->dim[2] = size_y;
    nii_grad->dim[3] = size_z;
    nii_grad->dim[4] = 26;
    nifti_update_dims_from_array(nii_grad);
    nii_grad->nvox = nii_input->nvox * 26;
    nii_grad->nbyper = sizeof(float);
    nii_grad->data = calloc(nii_grad->nvox, nii_grad->nbyper);
    nii_grad->scl_slope = 1;
    float* nii_grad_data = static_cast<float*>(nii_grad->data);

    // ========================================================================
    // Gradients
    // ========================================================================
    uint32_t ix, iy, iz, j;
    float d;

    for (uint32_t n = 0; n != nr_steps; ++n) {

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

            // ------------------------------------------------------------
            // 1-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 0 + i) = d;
            }
            if (ix < end_x) {
                j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 1 + i) = d;
            }
            if (iy > 0) {
                j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 2 + i) = d;
            }
            if (iy < end_y) {
                j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 3 + i) = d;
            }
            if (iz > 0) {
                j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 4 + i) = d;
            }
            if (iz < end_z) {
                j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 5 + i) = d;
            }

            // ------------------------------------------------------------
            // 2-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 6 + i) = d;
            }
            if (ix > 0 && iy < end_y) {
                j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 7 + i) = d;
            }
            if (ix < end_x && iy > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 8 + i) = d;
            }
            if (ix < end_x && iy < end_y) {
                j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 9 + i) = d;
            }
            if (iy > 0 && iz > 0) {
                j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 10 + i) = d;
            }
            if (iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 11 + i) = d;
            }
            if (iy < end_y && iz > 0) {
                j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 12 + i) = d;
            }
            if (iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 13 + i) = d;
            }
            if (ix > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 14 + i) = d;
            }
            if (ix < end_x && iz > 0) {
                j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 15 + i) = d;
            }
            if (ix > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 16 + i) = d;
            }
            if (ix < end_x && iz < end_z) {
                j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 17 + i) = d;
            }

            // ------------------------------------------------------------
            // 3-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 18 + i) = d;
            }
            if (ix > 0 && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 19 + i) = d;
            }
            if (ix > 0 && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 20 + i) = d;
            }
            if (ix < end_x && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 21 + i) = d;
            }
            if (ix > 0 && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 22 + i) = d;
            }
            if (ix < end_x && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 23 + i) = d;
            }
            if (ix < end_x && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 24 + i) = d;
            }
            if (ix < end_x && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                d = *(nii_input_data + i) - *(nii_input_data + j);
                *(nii_grad_data + nr_voxels * 25 + i) = d;
            }
        }

        // save_output_nifti(fout, "gradients", nii_grad, true);

        // ========================================================================
        // Diffuse
        // ========================================================================
        float d_0, d_1, d_2, d_3;

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
            d_1 = (*(nii_input_data + i) / 26.) * 1.;
            d_2 = (*(nii_input_data + i) / 26.) * (1. - std::sqrt(2)/2.);
            d_3 = (*(nii_input_data + i) / 26.) * (1. - std::sqrt(3)/2.);
            d_0 = *(nii_input_data + i) - ((d_1 * 6.) + (d_2 * 12.) + (d_3 * 8.));

            // ------------------------------------------------------------
            // 1-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_1;
            }
            if (ix < end_x) {
                j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_1;
            }
            if (iy > 0) {
                j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_1;
            }
            if (iy < end_y) {
                j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_1;
            }
            if (iz > 0) {
                j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_1;
            }
            if (iz < end_z) {
                j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_1;
            }

            // ------------------------------------------------------------
            // 2-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (ix > 0 && iy < end_y) {
                j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (ix < end_x && iy > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (ix < end_x && iy < end_y) {
                j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (iy > 0 && iz > 0) {
                j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (iy < end_y && iz > 0) {
                j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (ix > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (ix < end_x && iz > 0) {
                j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (ix > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }
            if (ix < end_x && iz < end_z) {
                j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_2;
            }

            // ------------------------------------------------------------
            // 3-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }
            if (ix > 0 && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }
            if (ix > 0 && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }
            if (ix < end_x && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }
            if (ix > 0 && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }
            if (ix < end_x && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }
            if (ix < end_x && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }
            if (ix < end_x && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += d_3;
            }

            *(nii_diffuse_data + i) += d_0;
        }

        // Replace previous image & reset diffusion step
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_input_data + i) = *(nii_diffuse_data + i);
        }
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_diffuse_data + i) = 0;
        }

    }
    save_output_nifti(fout, "diffuse", nii_input, true);

    cout << "\n  Finished." << endl;
    return 0;
}
