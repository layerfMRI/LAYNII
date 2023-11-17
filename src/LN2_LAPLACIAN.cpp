#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_LAPLACIAN: Compute image Laplacian.\n"
    "\n"
    "Usage:\n"
    "    LN2_LAPLACIAN -input input.nii\n"
    "    ../LN2_LAPLACIAN -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -input         : Nifti image that will be used to compute gradients.\n"
    "                     This can be a 4D nifti. in 4D case, 3D gradients\n"
    "                     will be computed for each volume.\n"
    "    -output        : (Optional) Output basename for all outputs.\n"
    "    -debug         : (Optional) Save extra intermediate outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
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

    log_welcome("LN2_LAPLACIAN");
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

    // Prepare outputs
    nifti_image *nii_laplacian = copy_nifti_as_float32(nii_input);
    float *nii_laplacian_data = static_cast<float*>(nii_laplacian->data);

    // Set to zero
    for (uint32_t i = 0; i != nr_voxels*size_time; ++i) {
        *(nii_laplacian_data + i) = 0;
    }

    // Prepare intermediate outputs
    nifti_image *nii_gra_x = copy_nifti_as_float32(nii_laplacian);
    float *nii_gra_x_data = static_cast<float*>(nii_gra_x->data);
    nifti_image *nii_gra_y = copy_nifti_as_float32(nii_laplacian);
    float *nii_gra_y_data = static_cast<float*>(nii_gra_y->data);
    nifti_image *nii_gra_z = copy_nifti_as_float32(nii_laplacian);
    float *nii_gra_z_data = static_cast<float*>(nii_gra_z->data);

    // ========================================================================
    // Gradients
    // ========================================================================
    cout << "  Computing first derivatives..." << endl;

    uint32_t ix, iy, iz, it, j, k;

    for (uint32_t t = 0; t != size_time; ++t) {
        cout << "\r    Volume: " << t+1 << "/" << size_time << flush;

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz, it) = ind2sub_4D(i+nr_voxels*t, size_x, size_y, size_z);
            float gra_x = 0, gra_y = 0, gra_z = 0;

            // ----------------------------------------------------------------
            // 1-jump neighbours
            // ----------------------------------------------------------------
            if (ix > 0 && ix < end_x) {
                j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
                gra_x = *(nii_input_data + j) - *(nii_input_data + k);
            }
            if (iy > 0 && iy < end_y) {
                j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
                gra_y = *(nii_input_data + j) - *(nii_input_data + k);
            }
            if (iz > 0 && iz < end_z) {
                j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
                gra_z = *(nii_input_data + j) - *(nii_input_data + k);
            }

            // ----------------------------------------------------------------
            *(nii_gra_x_data + i+nr_voxels*t) = gra_x;
            *(nii_gra_y_data + i+nr_voxels*t) = gra_y;
            *(nii_gra_z_data + i+nr_voxels*t) = gra_z;                                
        }
    }
    cout << endl;

    if (mode_debug) {
        cout << "  Saving intermediate outputs..." << endl;
        save_output_nifti(fout, "gradient_x", nii_gra_x, true);
        save_output_nifti(fout, "gradient_y", nii_gra_y, true);
        save_output_nifti(fout, "gradient_z", nii_gra_z, true);
    }

    // ========================================================================
    // Compute Laplacian
    // ========================================================================
    cout << "  Computing second derivatives for Laplacian..." << endl;

    for (uint32_t t = 0; t != size_time; ++t) {
        cout << "\r    Volume: " << t+1 << "/" << size_time << flush;

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz, it) = ind2sub_4D(i+nr_voxels*t, size_x, size_y, size_z);
            float gra_xx = 0, gra_yy = 0, gra_zz = 0;

            // ----------------------------------------------------------------
            // 1-jump neighbours
            // ----------------------------------------------------------------
            if (ix > 0 && ix < end_x) {
                j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
                gra_xx = *(nii_gra_x_data + j) - *(nii_gra_x_data + k);
            }
            if (iy > 0 && iy < end_y) {
                j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
                gra_yy = *(nii_gra_y_data + j) - *(nii_gra_y_data + k);
            }
            if (iz > 0 && iz < end_z) {
                j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
                gra_zz = *(nii_gra_z_data + j) - *(nii_gra_z_data + k);
            }

            // ----------------------------------------------------------------
            *(nii_laplacian_data + i+nr_voxels*t) = gra_xx + gra_yy + gra_zz;
        }
    }
    cout << endl;

    cout << "  Saving output..." << endl;
    save_output_nifti(fout, "laplacian", nii_laplacian, true);
    
    cout << "\n  Finished." << endl;
    return 0;
}
