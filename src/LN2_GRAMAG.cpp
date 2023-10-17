#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_GRAMAG: Compute gradient magnitude image.\n"
    "\n"
    "Usage:\n"
    "    LN2_GRAMAG -input input.nii\n"
    "    ../LN2_GRAMAG -input input.nii\n"
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
    "Reference / further reading:\n"
    "    [See Figure 1 from] Gulban, O.F., Schneider, M., Marquardt, I., \n"
    "    Haast, R.A.M., De Martino, F., 2018. A scalable method to improve \n"
    "    gray matter segmentation at ultra high field MRI. PloS one 13, e0198335.\n"
    "    <https://doi.org/10.1371/journal.pone.0198335>\n"
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

    log_welcome("LN2_GRAMAG");
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
    nifti_image* nii_gramag = copy_nifti_as_float32(nii_input);
    float* nii_gramag_data = static_cast<float*>(nii_gramag->data);

    // Set to zero
    for (uint32_t i = 0; i != nr_voxels*size_time; ++i) {
        *(nii_gramag_data + i) = 0;
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

                // ----------------------------------------------------------------
                // 1-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0 && ix < end_x) {
                    j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
                    gra_x += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_1 += 1;
                }
                if (iy > 0 && iy < end_y) {
                    j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
                    gra_y += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_1 += 1;
                }
                if (iz > 0 && iz < end_z) {
                    j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
                    gra_z += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_1 += 1;
                }

                // ----------------------------------------------------------------
                // 2-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0 && iy > 0 && ix < end_x && iy < end_y) {
                    j = sub2ind_4D(ix-1, iy-1, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy+1, iz, it, size_x, size_y, size_z);
                    g21 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_2 += 1;
                }
                if (ix > 0 && iy < end_y && ix < end_x && iy > 0) {
                    j = sub2ind_4D(ix-1, iy+1, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy-1, iz, it, size_x, size_y, size_z);
                    g22 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_2 += 1;
                }
                if (iy > 0 && iz > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_4D(ix, iy-1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy+1, iz+1, it, size_x, size_y, size_z);
                    g23 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_2 += 1;
                }
                if (iy > 0 && iz < end_z && iy < end_y && iz > 0) {
                    j = sub2ind_4D(ix, iy-1, iz+1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy+1, iz-1, it, size_x, size_y, size_z);
                    g24 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_2 += 1;
                }
                if (ix > 0 && iz > 0 && ix < end_x && iz < end_z) {
                    j = sub2ind_4D(ix-1, iy, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy, iz+1, it, size_x, size_y, size_z);
                    g25 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_2 += 1;
                }
                if (ix < end_x && iz > 0 && ix > 0 && iz < end_z) {
                    j = sub2ind_4D(ix+1, iy, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix-1, iy, iz+1, it, size_x, size_y, size_z);
                    g26 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_2 += 1;
                }

                // ----------------------------------------------------------------
                // 3-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0 && ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_4D(ix-1, iy-1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy+1, iz+1, it, size_x, size_y, size_z);
                    g31 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_3 += 1;
                }
                if (ix > 0 && iy > 0 && iz < end_z && ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_4D(ix-1, iy-1, iz+1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy+1, iz-1, it, size_x, size_y, size_z);
                    g32 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_3 += 1;
                }
                if (ix > 0 && iy < end_y && iz > 0 && ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_4D(ix-1, iy+1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy-1, iz+1, it, size_x, size_y, size_z);
                    g33 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_3 += 1;
                }
                if (ix < end_x && iy > 0 && iz > 0 && ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_4D(ix+1, iy-1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix-1, iy+1, iz+1, it, size_x, size_y, size_z);
                    g34 += std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    count_3 += 1;
                }

                // ----------------------------------------------------------------
                // Compute magnitude
                // *(nii_gramag_data + i) = sqrt(gra_x*gra_x + gra_y*gra_y + gra_z*gra_z);
                // *(nii_gramag_data + i) = sqrt(g21*g21 + g22*g22 + g23*g23 + g24*g24 + g25*g25 + g26*g26);
                // *(nii_gramag_data + i) = sqrt(g31*g31 + g32*g32 + g33*g33 + g34*g34);

                // Average rate of change across spheres (shells)
                // NOTE[Faruk]: This is a bit of experimental thinking... Need to think
                // if thinking neighbors as separate spherical shells has some benefits...
                *(nii_gramag_data + i) += (gra_x + gra_y + gra_z) / count_1;
                *(nii_gramag_data + i) += (g21 + g22 + g23 + g24 + g25 + g26) / count_2;
                *(nii_gramag_data + i) += (g31 + g32 + g33 + g34) / count_3;
                *(nii_gramag_data + i) /= (count_1 / 3) + (count_2 / 6) + (count_3 / 4);
            }
        }
        cout << "  Saving output..." << endl;
        save_output_nifti(fout, "gramag", nii_gramag, true);
    } else {
        cout << "  Circular difference mode (-pi to pi range) is selected..." << endl;
        const float ONEPI = 3.14159265358979f;
        const float TWOPI = 2.0f * 3.14159265358979f;

        for (uint32_t t = 0; t != size_time; ++t) {
            cout << "    Volume: " << t+1 << "/" << size_time << endl;

            for (uint32_t i = 0; i != nr_voxels; ++i) {
                tie(ix, iy, iz, it) = ind2sub_4D(i+nr_voxels*t, size_x, size_y, size_z);
                float gra_x, gra_y, gra_z, diff1, diff2, diff3;
                float g21, g22, g23, g24, g25, g26;
                float g31, g32, g33, g34;
                float count_1 = 0, count_2 = 0, count_3 = 0;

                // ----------------------------------------------------------------
                // 1-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0 && ix < end_x) {
                    j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    gra_x = std::min(diff1, std::min(diff2, diff3));
                    count_1 += 1;
                }
                if (iy > 0 && iy < end_y) {
                    j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    gra_y = std::min(diff1, std::min(diff2, diff3));
                    count_1 += 1;
                }
                if (iz > 0 && iz < end_z) {
                    j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    gra_z = std::min(diff1, std::min(diff2, diff3));
                    count_1 += 1;
                }

                // ----------------------------------------------------------------
                // 2-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0 && iy > 0 && ix < end_x && iy < end_y) {
                    j = sub2ind_4D(ix-1, iy-1, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy+1, iz, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g21 = std::min(diff1, std::min(diff2, diff3));
                    count_2 += 1;
                }
                if (ix > 0 && iy < end_y && ix < end_x && iy > 0) {
                    j = sub2ind_4D(ix-1, iy+1, iz, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy-1, iz, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g22 = std::min(diff1, std::min(diff2, diff3));
                    count_2 += 1;
                }
                if (iy > 0 && iz > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_4D(ix, iy-1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy+1, iz+1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g23 = std::min(diff1, std::min(diff2, diff3));
                    count_2 += 1;
                }
                if (iy > 0 && iz < end_z && iy < end_y && iz > 0) {
                    j = sub2ind_4D(ix, iy-1, iz+1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix, iy+1, iz-1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g24 = std::min(diff1, std::min(diff2, diff3));
                    count_2 += 1;
                }
                if (ix > 0 && iz > 0 && ix < end_x && iz < end_z) {
                    j = sub2ind_4D(ix-1, iy, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy, iz+1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g25 = std::min(diff1, std::min(diff2, diff3));
                    count_2 += 1;
                }
                if (ix < end_x && iz > 0 && ix > 0 && iz < end_z) {
                    j = sub2ind_4D(ix+1, iy, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix-1, iy, iz+1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g26 = std::min(diff1, std::min(diff2, diff3));
                    count_2 += 1;
                }

                // ----------------------------------------------------------------
                // 3-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0 && ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_4D(ix-1, iy-1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy+1, iz+1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g31 = std::min(diff1, std::min(diff2, diff3));
                    count_3 += 1;
                }
                if (ix > 0 && iy > 0 && iz < end_z && ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_4D(ix-1, iy-1, iz+1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy+1, iz-1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g32 = std::min(diff1, std::min(diff2, diff3));
                    count_3 += 1;
                }
                if (ix > 0 && iy < end_y && iz > 0 && ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_4D(ix-1, iy+1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix+1, iy-1, iz+1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g33 = std::min(diff1, std::min(diff2, diff3));
                    count_3 += 1;
                }
                if (ix < end_x && iy > 0 && iz > 0 && ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_4D(ix+1, iy-1, iz-1, it, size_x, size_y, size_z);
                    k = sub2ind_4D(ix-1, iy+1, iz+1, it, size_x, size_y, size_z);
                    diff1 = std::abs(*(nii_input_data + j) - *(nii_input_data + k));
                    diff2 = *(nii_input_data + j) - *(nii_input_data + k) + TWOPI;
                    diff3 = *(nii_input_data + k) - *(nii_input_data + j) + TWOPI;
                    g34 = std::min(diff1, std::min(diff2, diff3));
                    count_3 += 1;
                }
            // Average differences across neighbourhood types
            *(nii_gramag_data + i+nr_voxels*t) += (gra_x + gra_y + gra_z) / count_1;
            *(nii_gramag_data + i+nr_voxels*t) += (g21 + g22 + g23 + g24 + g25 + g26) / count_2;
            *(nii_gramag_data + i+nr_voxels*t) += (g31 + g32 + g33 + g34) / count_3;
            *(nii_gramag_data + i+nr_voxels*t) /= (count_1 / 3) + (count_2 / 6) + (count_3 / 4);
            }
        }
        cout << "  Saving output..." << endl;
        save_output_nifti(fout, "gramag_circular", nii_gramag, true);
    }

    cout << "\n  Finished." << endl;
    return 0;
}
