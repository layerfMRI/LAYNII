
#include "../dep/laynii_lib.h"
#include <limits>

int shoq_help(void) {
    printf(
    "LN2_DIFFUSE: WORK IN PROGRESS... EXTREMELY EXPERIMENTAL!\n"
    "\n"
    "Usage:\n"
    "    LN2_DIFFUSE -input data.nii -nr_steps 5\n"
    "\n"
    "Options:\n"
    "    -help     : Show this help.\n"
    "    -input    : Input image.\n"
    "    -nr_steps : Number of diffusion steps.\n"
    "    -debug    : (Optional) Save extra intermediate outputs.\n"
    "    -output   : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    uint16_t ac, nr_steps = 5;
    bool mode_debug = false;

    // Process user options
    if (argc < 2) return shoq_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return shoq_help();
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
        *(nii_diffuse_data + i) = 0.;
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
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 0 + i) = d;
            }
            if (ix < end_x) {
                j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 1 + i) = d;
            }
            if (iy > 0) {
                j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 2 + i) = d;
            }
            if (iy < end_y) {
                j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 3 + i) = d;
            }
            if (iz > 0) {
                j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 4 + i) = d;
            }
            if (iz < end_z) {
                j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 5 + i) = d;
            }

            // ------------------------------------------------------------
            // 2-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 6 + i) = d;
            }
            if (ix > 0 && iy < end_y) {
                j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 7 + i) = d;
            }
            if (ix < end_x && iy > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 8 + i) = d;
            }
            if (ix < end_x && iy < end_y) {
                j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 9 + i) = d;
            }
            if (iy > 0 && iz > 0) {
                j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 10 + i) = d;
            }
            if (iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 11 + i) = d;
            }
            if (iy < end_y && iz > 0) {
                j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 12 + i) = d;
            }
            if (iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 13 + i) = d;
            }
            if (ix > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 14 + i) = d;
            }
            if (ix < end_x && iz > 0) {
                j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 15 + i) = d;
            }
            if (ix > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 16 + i) = d;
            }
            if (ix < end_x && iz < end_z) {
                j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 17 + i) = d;
            }

            // ------------------------------------------------------------
            // 3-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 18 + i) = d;
            }
            if (ix > 0 && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 19 + i) = d;
            }
            if (ix > 0 && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 20 + i) = d;
            }
            if (ix < end_x && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 21 + i) = d;
            }
            if (ix > 0 && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 22 + i) = d;
            }
            if (ix < end_x && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 23 + i) = d;
            }
            if (ix < end_x && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 24 + i) = d;
            }
            if (ix < end_x && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                d = std::abs(*(nii_input_data + i) - *(nii_input_data + j));
                *(nii_grad_data + nr_voxels * 25 + i) = d;
            }
        }

        // save_output_nifti(fout, "gradients", nii_grad, true);

        // ========================================================================
        // Diffuse
        // ========================================================================
        float d_1, d_2, d_3;
        float q_01, q_02, q_03, q_04, q_05, q_06, q_07, q_08, q_09, q_10, q_11;
        float q_12, q_13, q_14, q_15, q_16, q_17, q_18, q_19, q_20, q_21, q_22;
        float q_23, q_24, q_25, q_26, q_00;

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
            d_1 = (*(nii_input_data + i) / 26.) * 1.;
            d_2 = (*(nii_input_data + i) / 26.) * (1. - std::sqrt(2)/2.);
            d_3 = (*(nii_input_data + i) / 26.) * (1. - std::sqrt(3)/2.);

            // Quantities that will be moved
            q_01 = *(nii_grad_data + 0 * nr_voxels + i) * d_1;
            q_02 = *(nii_grad_data + 1 * nr_voxels + i) * d_1;
            q_03 = *(nii_grad_data + 2 * nr_voxels + i) * d_1;
            q_04 = *(nii_grad_data + 3 * nr_voxels + i) * d_1;
            q_05 = *(nii_grad_data + 4 * nr_voxels + i) * d_1;
            q_06 = *(nii_grad_data + 5 * nr_voxels + i) * d_1;

            q_07 = *(nii_grad_data + 6 * nr_voxels + i) * d_2;
            q_08 = *(nii_grad_data + 7 * nr_voxels + i) * d_2;
            q_09 = *(nii_grad_data + 8 * nr_voxels + i) * d_2;
            q_10 = *(nii_grad_data + 9 * nr_voxels + i) * d_2;
            q_11 = *(nii_grad_data + 10 * nr_voxels + i) * d_2;
            q_12 = *(nii_grad_data + 11 * nr_voxels + i) * d_2;
            q_13 = *(nii_grad_data + 12 * nr_voxels + i) * d_2;
            q_14 = *(nii_grad_data + 13 * nr_voxels + i) * d_2;
            q_15 = *(nii_grad_data + 14 * nr_voxels + i) * d_2;
            q_16 = *(nii_grad_data + 15 * nr_voxels + i) * d_2;
            q_17 = *(nii_grad_data + 16 * nr_voxels + i) * d_2;
            q_18 = *(nii_grad_data + 17 * nr_voxels + i) * d_2;

            q_19 = *(nii_grad_data + 18 * nr_voxels + i) * d_3;
            q_20 = *(nii_grad_data + 19 * nr_voxels + i) * d_3;
            q_21 = *(nii_grad_data + 20 * nr_voxels + i) * d_3;
            q_22 = *(nii_grad_data + 21 * nr_voxels + i) * d_3;
            q_23 = *(nii_grad_data + 22 * nr_voxels + i) * d_3;
            q_24 = *(nii_grad_data + 23 * nr_voxels + i) * d_3;
            q_25 = *(nii_grad_data + 24 * nr_voxels + i) * d_3;
            q_26 = *(nii_grad_data + 25 * nr_voxels + i) * d_3;

            q_00 = *(nii_input_data + i) - (q_01 + q_02 + q_03 + q_04 + q_05
                + q_06 + q_07 + q_08 + q_09 + q_10 + q_11 + q_12 + q_13 + q_14
                + q_15 + q_16 + q_17 + q_18 + q_19 + q_20 + q_21 + q_22 + q_23
                + q_24 + q_25 + q_26);

            // ------------------------------------------------------------
            // 1-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_01;
            }
            if (ix < end_x) {
                j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_02;
            }
            if (iy > 0) {
                j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_03;
            }
            if (iy < end_y) {
                j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_04;
            }
            if (iz > 0) {
                j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_05;
            }
            if (iz < end_z) {
                j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_06;
            }

            // ------------------------------------------------------------
            // 2-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_07;
            }
            if (ix > 0 && iy < end_y) {
                j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_08;
            }
            if (ix < end_x && iy > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_09;
            }
            if (ix < end_x && iy < end_y) {
                j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                *(nii_diffuse_data + j) += q_10;
            }
            if (iy > 0 && iz > 0) {
                j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_11;
            }
            if (iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_12;
            }
            if (iy < end_y && iz > 0) {
                j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_13;
            }
            if (iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_14;
            }
            if (ix > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_15;
            }
            if (ix < end_x && iz > 0) {
                j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_16;
            }
            if (ix > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_17;
            }
            if (ix < end_x && iz < end_z) {
                j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_18;
            }

            // ------------------------------------------------------------
            // 3-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_19;
            }
            if (ix > 0 && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_20;
            }
            if (ix > 0 && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_21;
            }
            if (ix < end_x && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_22;
            }
            if (ix > 0 && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_23;
            }
            if (ix < end_x && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_24;
            }
            if (ix < end_x && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                *(nii_diffuse_data + j) += q_25;
            }
            if (ix < end_x && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                *(nii_diffuse_data + j) += q_26;
            }

            *(nii_diffuse_data + i) += q_00;
        }

        // Replace previous image & reset diffusion step
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_input_data + i) = *(nii_diffuse_data + i);
        }
        for (uint32_t i = 0; i != nr_voxels; ++i) {
            *(nii_diffuse_data + i) = 0.;
        }

    }
    save_output_nifti(fout, "diffuse", nii_input, true);

    cout << "\n  Finished." << endl;
    return 0;
}
