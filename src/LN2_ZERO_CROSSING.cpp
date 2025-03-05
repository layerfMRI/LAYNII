#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_ZERO_CROSSING: Find zero crossing voxels in a scalar field.\n"
    "Usage:\n"
    "    LN2_ZERO_CROSSING -values equidist_metric.nii.gz -domain mask.nii.gz \n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -values : 3D nifti file that contains a scalar field. E.g. metric file\n"
    "              outputs from LN2_LAYERS program, or coordinate outputs from\n"
    "              LN2_MULTILATERATE program.\n"
    "    -domain : 3D nifti file that contains non-zero voxel where zero crossings\n"
    "              will be computed. In other words, a mask file.\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fin2 = NULL, *fout = NULL;
    int ac;

    // Process user options
    if (argc < 2) return show_help();
    for (ac = 1; ac < argc; ac++) {
        if (!strncmp(argv[ac], "-h", 2)) {
            return show_help();
        } else if (!strcmp(argv[ac], "-values")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -values\n");
                return 1;
            }
            fin1 = argv[ac];
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-domain")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -domain\n");
                return 1;
            }
            fin2 = argv[ac];
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
    if (!fin2) {
        fprintf(stderr, "** missing option '-domain'\n");
        return 1;
    }

    // Read input dataset, including data
    nii1 = nifti_image_read(fin1, 1);
    if (!nii1) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }
    nii2 = nifti_image_read(fin2, 1);
    if (!nii2) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin2);
        return 2;
    }


    log_welcome("LN2_ZERO_CROSSING");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

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
    nifti_image* nii_values = copy_nifti_as_float32(nii1);
    float* nii_values_data = static_cast<float*>(nii_values->data);
    nifti_image* nii_domain = copy_nifti_as_int32(nii2);
    int32_t* nii_domain_data = static_cast<int32_t*>(nii_domain->data);
    free(nii2);

    // Output nifti
    nifti_image* nii_out = copy_nifti_as_int32(nii1);
    int32_t* nii_out_data = static_cast<int32_t*>(nii_out->data);
    free(nii1);

    // Clean output array
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_out_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_domain_data + i) != 0) {
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_domain_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Find zero crossing neighboring voxels
    // ========================================================================
    cout << "\n  Finding zero crossing neighbour voxels..." << endl;

    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        uint32_t ix, iy, iz, i, j;
        i = *(voi_id + ii);
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

        // Check sign changes between neighbouring voxels
        float m = *(nii_values_data + i);
        float n;

        // --------------------------------------------------------------------
        // 1-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x) {
            j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iy > 0) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iy < end_y) {
            j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iz > 0) {
            j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iz < end_z) {
            j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        // --------------------------------------------------------------------
        // 2-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && iy > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix > 0 && iy < end_y) {
            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iy > 0) {
            j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iy < end_y) {
            j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iy > 0 && iz > 0) {
            j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iy < end_y && iz > 0) {
            j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iz > 0) {
            j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix > 0 && iz < end_z) {
            j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iz < end_z) {
            j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }

        // --------------------------------------------------------------------
        // 3-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && iy > 0 && iz > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix > 0 && iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix > 0 && iy < end_y && iz > 0) {
            j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iy > 0 && iz > 0) {
            j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix > 0 && iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iy > 0 && iz < end_z) {
            j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iy < end_y && iz > 0) {
            j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
        if (ix < end_x && iy < end_y && iz < end_z) {
            j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
            n = *(nii_values_data + j);
            if (*(nii_domain_data + j) != 0) {
                if (signbit(m) - signbit(n) != 0) {
                    *(nii_out_data + i) = 1;
                }
            }
        }
    }

    // Save output nifti
    save_output_nifti(fout, "zero_crossing", nii_out, true);
    cout << "\n  Finished." << endl;
    return 0;
}
