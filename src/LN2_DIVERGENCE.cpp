#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_DIVERGENCE: WIP Compute divergence on the gradients of a scalar image.\n"
    "\n"
    "Usage:\n"
    "    LN2_DIVERGENCE -input input.nii\n"
    "    ../LN2_DIVERGENCE -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help           : Show this help.\n"
    "    -input          : Nifti image that will be used to compute gradients.\n"
    "                      This can be a 4D nifti. in 4D case, 3D gradients\n"
    "                      will be computed for each volume.\n"
    "    -circular       : Gradients are computed using circular difference. Needed\n"
    "                      when the input contains e.g. phase values (0 to 2*pi).\n"
    "                      Input range is assumed to be 2*pi.\n"
    "    -circular_int13 : (Optional) Cast the input range from [-4096 4096] to [0 2*pi].\n"
    "                      This option is often needed with Siemens phase images as they\n"
    "                      commonly appear to be uint12 range with scl_slope = 2, and\n"
    "                      scl_inter = -4096 in the header. Meaning that the intended range\n"
    "                      is int13, even though the data type is uint16 and only int12 portion\n"
    "                      is used to store the phase values.\n"
    "    -output         : (Optional) Output basename for all outputs.\n"
    "    -debug          : (Optional) Save extra intermediate outputs.\n"
    "\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    bool mode_circular = false, mode_circular_int13 = false, mode_debug = false;

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
        } else if (!strcmp(argv[ac], "-circular")) {
            mode_circular = true;
        } else if (!strcmp(argv[ac], "-circular_int13")) {
            mode_circular_int13 = true;
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
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

    log_welcome("LN2_DIVERGENCE");
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
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_float32_with_scl_slope_and_scl_inter(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // Prepare output image
    nifti_image* nii_gra_x = copy_nifti_as_float32(nii_input);
    float* nii_gra_x_data = static_cast<float*>(nii_gra_x->data);

    nifti_image* nii_gra_y = copy_nifti_as_float32(nii_input);
    float* nii_gra_y_data = static_cast<float*>(nii_gra_y->data);

    nifti_image* nii_gra_z = copy_nifti_as_float32(nii_input);
    float* nii_gra_z_data = static_cast<float*>(nii_gra_z->data);


    nifti_image* nii_gra2_x = copy_nifti_as_float32(nii_input);
    float* nii_gra2_x_data = static_cast<float*>(nii_gra2_x->data);

    nifti_image* nii_gra2_y = copy_nifti_as_float32(nii_input);
    float* nii_gra2_y_data = static_cast<float*>(nii_gra2_y->data);

    nifti_image* nii_gra2_z = copy_nifti_as_float32(nii_input);
    float* nii_gra2_z_data = static_cast<float*>(nii_gra2_z->data);


    nifti_image* nii_divergence = copy_nifti_as_float32(nii_input);
    float* nii_divergence_data = static_cast<float*>(nii_divergence->data);

    // Set to zero
    for (uint32_t i = 0; i != nr_voxels*size_time; ++i) {
        *(nii_gra_x_data + i) = 0;
        *(nii_gra_y_data + i) = 0;
        *(nii_gra_z_data + i) = 0;
        *(nii_gra2_x_data + i) = 0;
        *(nii_gra2_y_data + i) = 0;
        *(nii_gra2_z_data + i) = 0;
        *(nii_divergence_data + i) = 0;
    }

    // ========================================================================
    // Convert ranges to 0 to 2*pi if opted for
    // ========================================================================
    if (mode_circular_int13) {
        cout << "  Casting [-4096 4096] to [0 2*pi] range..." << endl;
        for (uint32_t i = 0; i != nr_voxels*size_time; ++i) {
            float k = *(nii_input_data + i);
            *(nii_input_data + i) = ((k + 4096) / 8192) * 2*3.14159265358979323846;
        }
    }

    if (mode_debug) {
        save_output_nifti(fout, "float32", nii_input, false);
    }

    // ========================================================================
    // Find connected clusters
    // ========================================================================
    cout << "  Computing gradients..." << endl;

    uint32_t ix, iy, iz, it, j, k;

    if (mode_circular == false) {
        for (uint32_t t = 0; t != size_time; ++t) {
            cout << "    Volume: " << t+1 << "/" << size_time << endl;

            for (uint32_t i = 0; i != nr_voxels; ++i) {
                tie(ix, iy, iz, it) = ind2sub_4D(i, size_x, size_y, size_z);
                float gra_x = 0, gra_y = 0, gra_z = 0;
                float g21 = 0, g22 = 0, g23 = 0, g24 = 0, g25 = 0, g26 = 0;
                float g31 = 0, g32 = 0, g33 = 0, g34 = 0;
                float count_1 = 0, count_2 = 0, count_3 = 0;

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0 && ix < end_x) {
                    j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
                    gra_x += *(nii_input_data + j) - *(nii_input_data + k);
                }
                if (iy > 0 && iy < end_y) {
                    j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
                    gra_y += *(nii_input_data + j) - *(nii_input_data + k);
                }
                if (iz > 0 && iz < end_z) {
                    j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
                    gra_z += *(nii_input_data + j) - *(nii_input_data + k);
                }

                *(nii_gra_x_data + i) += gra_x;
                *(nii_gra_y_data + i) += gra_y;
                *(nii_gra_z_data + i) += gra_z;
            }
        }
        cout << "  Saving gradients..." << endl;
        if (mode_debug) {
            save_output_nifti(fout, "gra_x", nii_gra_x, true);
            save_output_nifti(fout, "gra_y", nii_gra_y, true);
            save_output_nifti(fout, "gra_z", nii_gra_z, true);
        }
    // } else {
    }
    //     cout << "  Circular difference mode (-pi to pi range) is selected..." << endl;
    //     const float ONEPI = 3.14159265358979f;
    //     const float TWOPI = 2.0f * 3.14159265358979f;

    //     for (uint32_t t = 0; t != size_time; ++t) {
    //         cout << "    Volume: " << t+1 << "/" << size_time << endl;

    //         for (uint32_t i = 0; i != nr_voxels; ++i) {
    //             tie(ix, iy, iz, it) = ind2sub_4D(i+nr_voxels*t, size_x, size_y, size_z);
    //             float gra_x, gra_y, gra_z, diff1, diff2, diff3;

    //             // ------------------------------------------------------------
    //             // 1-jump neighbours
    //             // ------------------------------------------------------------
    //             if (ix > 0 && ix < end_x) {
    //                 j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
    //                 k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
    //                 diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
    //                 diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
    //                 diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
    //                 gra_x = std::min(diff1, std::min(diff2, diff3));
    //             }
    //             if (iy > 0 && iy < end_y) {
    //                 j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
    //                 k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
    //                 diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
    //                 diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
    //                 diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
    //                 gra_y = std::min(diff1, std::min(diff2, diff3));
    //             }
    //             if (iz > 0 && iz < end_z) {
    //                 j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
    //                 k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
    //                 diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
    //                 diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
    //                 diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
    //                 gra_z = std::min(diff1, std::min(diff2, diff3));
    //             }

    //             *(nii_gra_x_data + i) += gra_x;
    //             *(nii_gra_y_data + i) += gra_y;
    //             *(nii_gra_z_data + i) += gra_z;
    //         }
    //     }

    //     if (mode_debug) {
    //     cout << "  Saving gradients..." << endl;
    //         save_output_nifti(fout, "gra_x_circular", nii_gra_x, true);
    //         save_output_nifti(fout, "gra_y_circular", nii_gra_y, true);
    //         save_output_nifti(fout, "gra_z_circular", nii_gra_z, true);
    //     }
    // }

    // ========================================================================
    // Compute Divergence on the gradient vector field
    // ========================================================================
    for (uint32_t t = 0; t != size_time; ++t) {
        cout << "    Volume: " << t+1 << "/" << size_time << endl;

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz, it) = ind2sub_4D(i+nr_voxels*t, size_x, size_y, size_z);
            float gra_x, gra_y, gra_z;

            // ------------------------------------------------------------
            // 1-jump neighbours
            // ------------------------------------------------------------
            if (ix > 0 && ix < end_x) {
                j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
                gra_x = *(nii_gra_x_data + j) - *(nii_gra_x_data + k);
            }
            if (iy > 0 && iy < end_y) {
                j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
                gra_y = *(nii_gra_y_data + j) - *(nii_gra_y_data + k);
            }
            if (iz > 0 && iz < end_z) {
                j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
                gra_z = *(nii_gra_z_data + j) - *(nii_gra_z_data + k);
            }

            // NOTE[Faruk]: Remove these for optimal RAM in next iteration
            *(nii_gra2_x_data + i) = gra_x;
            *(nii_gra2_y_data + i) = gra_y;
            *(nii_gra2_z_data + i) = gra_z;

            *(nii_divergence_data + i) = gra_x + gra_y + gra_z;
        }
    }

    if (mode_debug) {
        save_output_nifti(fout, "gra2_x", nii_gra2_x, true);
        save_output_nifti(fout, "gra2_y", nii_gra2_y, true);
        save_output_nifti(fout, "gra2_z", nii_gra2_z, true);
    }

    save_output_nifti(fout, "divergence", nii_divergence, true);

    cout << "\n  Finished." << endl;
    return 0;
}
