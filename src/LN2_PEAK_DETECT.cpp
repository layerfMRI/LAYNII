#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_PEAK_DETECT: WIP... maximum filter to detect image peaks.\n"
    "\n"
    "Usage:\n"
    "    LN2_PEAK_DETECT -values activation.nii -max\n"
    "\n"
    "Options:\n"
    "    -help      : Show this help.\n"
    "    -values    : Nifti image with values that will be filtered.\n"
    "                 For example an activation map or anatomical T1w images.\n"
    "    -max       : (Default) Detect peaks with maximum filter.\n"
    "    -min       : Detect peaks with minimum filter.\n"
    "    -output    : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char* argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    bool mode_max = true, mode_min = false;

    // Process user options
    if (argc < 2) return show_help();

    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-help", 5)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-values")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -values\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-max")) {
            mode_max = true;
            mode_min = false;
        } else if (!strcmp(argv[ac], "-min")) {
            mode_max = false;
            mode_min = true;
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
        fprintf(stderr, "** missing option '-values'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    log_welcome("LN2_PEAK_DETECT");
    log_nifti_descriptives(nii1);

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;
    const uint32_t size_z = nii1->nz;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;
    const uint32_t end_z = size_z - 1;

    const uint32_t nr_voxels = size_z * size_y * size_x;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);

    // ========================================================================
    // Prepare outputs
    nifti_image* nii_output = copy_nifti_as_float32(nii_input);
    float* nii_output_data = static_cast<float*>(nii_output->data);

    // ========================================================================
    // Visit each voxel to check their 27 neighbours

    for (uint32_t i = 0; i != nr_voxels; ++i) {
        float ref = *(nii_input_data + i);
        float n_max = *(nii_input_data + i);
        float n_min = *(nii_input_data + i);

        uint32_t ix, iy, iz, j;
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

        // ------------------------------------------------------------
        // 1-jump neighbours
        // ------------------------------------------------------------
        if (ix > 0) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x) {
            j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iy > 0) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iy < end_y) {
            j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iz > 0) {
            j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iz < end_z) {
            j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        // ------------------------------------------------------------
        // 2-jump neighbours
        // ------------------------------------------------------------
        if (ix > 0 && iy > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix > 0 && iy < end_y) {
            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iy > 0) {
            j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iy < end_y) {
            j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iy > 0 && iz > 0) {
            j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iy < end_y && iz > 0) {
            j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iz > 0) {
            j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix > 0 && iz < end_z) {
            j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iz < end_z) {
            j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }

        // ------------------------------------------------------------
        // 3-jump neighbours
        // ------------------------------------------------------------
        if (ix > 0 && iy > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix > 0 && iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix > 0 && iy < end_y && iz > 0) {
            j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iy > 0 && iz > 0) {
            j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix > 0 && iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iy < end_y && iz > 0) {
            j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }
        if (ix < end_x && iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
            if (*(nii_input_data + j) > n_max) {
                n_max = *(nii_input_data + j);
            } else if (*(nii_input_data + j) < n_min) {
                n_min = *(nii_input_data + j);
            }
        }

        // Write results inside nifti
        if (mode_max) {
            if (n_max > ref || ref == 0) {
                *(nii_output_data + i) = 0;
            } else {
                *(nii_output_data + i) = 1;
            }
        } else if (mode_min) {
            if (n_min < ref || ref == 0) {
                *(nii_output_data + i) = 0;
            } else {
                *(nii_output_data + i) = 1;
            }
        }
    }

    save_output_nifti(fout, "peaks", nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
