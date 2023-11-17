#include "../dep/laynii_lib.h"
#include <sstream>

int show_help(void) {
    printf(
    "LN2_PHASE_GRADIENTS: Compute gradients of a phase image using circular difference.\n"
    "\n"
    "Usage:\n"
    "    LN2_PHASE_GRADIENTS -input input.nii\n"
    "    ../LN2_PHASE_GRADIENTS -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help          : Show this help.\n"
    "    -input         : Nifti image that will be used to compute gradients.\n"
    "                     This can be a 4D nifti. in 4D case, 3D gradients\n"
    "                     will be computed for each volume.\n"
    "    -int13         : (Optional) Cast the input range from [-4096 4096] to [0 2*pi].\n"
    "                     This option is often needed with Siemens phase images as they\n"
    "                     commonly appear to be uint12 range with scl_slope = 2, and\n"
    "                     scl_inter = -4096 in the header. Meaning that the intended range\n"
    "                     is int13, even though the data type is uint16 and only int12 portion\n"
    "                     is used to store the phase values.\n"
    "    -merge_outputs : (Optional) Save gradients as a 4D file. Only works\n"
    "                     with 3D images.\n"
    "    -output        : (Optional) Output basename for all outputs.\n"
    "    -debug         : (Optional) Save extra intermediate outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    bool mode_int13 = false, mode_debug = false, mode_merge_outputs = false;

    nifti_image *nii_gra = NULL;
    float *nii_gra_data = NULL;

    nifti_image *nii_gra_x = NULL, *nii_gra_y = NULL, *nii_gra_z = NULL;
    float *nii_gra_x_data = NULL, *nii_gra_y_data = NULL, *nii_gra_z_data = NULL;

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
        } else if (!strcmp(argv[ac], "-int13")) {
            mode_int13 = true;
        } else if (!strcmp(argv[ac], "-merge_outputs")) {
            mode_merge_outputs = true;
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

    log_welcome("LN2_PHASE_GRADIENTS");
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

    // Prepare output images
    if (mode_merge_outputs) {
        cout << "  Input is a 3D image." << endl;
        // Create a 4D nifti image for point distances
        nii_gra = nifti_copy_nim_info(nii_input);
        nii_gra->dim[0] = 4;  // For proper 4D nifti
        nii_gra->dim[1] = size_x;
        nii_gra->dim[2] = size_y;
        nii_gra->dim[3] = size_z;
        nii_gra->dim[4] = 3;
        nifti_update_dims_from_array(nii_gra);
        nii_gra->nvox = nr_voxels * 3;
        nii_gra->nbyper = sizeof(float);
        nii_gra->data = calloc(nii_gra->nvox, nii_gra->nbyper);
        
        nii_gra_data = static_cast<float*>(nii_gra->data);

        // Set to zero
        for (uint32_t i = 0; i != nr_voxels*3; ++i) {
            *(nii_gra_data + i) = 0;
        }
    } else {
        cout << "  Input is a 4D image (e.g. timeseries)." << endl;
        nii_gra_x = copy_nifti_as_float32(nii_input);
        nii_gra_x_data = static_cast<float*>(nii_gra_x->data);
        nii_gra_y = copy_nifti_as_float32(nii_input);
        nii_gra_y_data = static_cast<float*>(nii_gra_y->data);
        nii_gra_z = copy_nifti_as_float32(nii_input);
        nii_gra_z_data = static_cast<float*>(nii_gra_z->data);

        // Set to zero
        for (uint32_t i = 0; i != nr_voxels*size_time; ++i) {
            *(nii_gra_x_data + i) = 0;
            *(nii_gra_y_data + i) = 0;
            *(nii_gra_z_data + i) = 0;        
        }
    }

    // ========================================================================
    // Convert ranges to 0 to 2*pi if opted for
    // ========================================================================
    if (mode_int13) {
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

    const float ONEPI = 3.14159265358979f;
    const float TWOPI = 2.0f * 3.14159265358979f;

    for (uint32_t t = 0; t != size_time; ++t) {
        cout << "\r    Volume: " << t+1 << "/" << size_time << flush;

        for (uint32_t i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz, it) = ind2sub_4D(i+nr_voxels*t, size_x, size_y, size_z);
            float gra_x=0, gra_y=0, gra_z=0, diff1=0, diff2=0, diff3=0, a=0, b=0;

            // ----------------------------------------------------------------
            // 1-jump neighbours
            // ----------------------------------------------------------------
            if (ix > 0 && ix < end_x) {
                j = sub2ind_4D(ix-1, iy, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix+1, iy, iz, it, size_x, size_y, size_z);
                a = *(nii_input_data + j);
                b = *(nii_input_data + k);
                diff1 = b - a;
                diff2 = TWOPI + b - a;
                diff3 = b - (TWOPI + a);
                if (std::abs(diff2) < std::abs(diff1)) {
                    gra_x = diff2;
                } else if (std::abs(diff3) < std::abs(diff1)) {
                    gra_x = diff3;
                } else {
                    gra_x = diff1;
                }
            }
            if (iy > 0 && iy < end_y) {
                j = sub2ind_4D(ix, iy-1, iz, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy+1, iz, it, size_x, size_y, size_z);
                a = *(nii_input_data + j);
                b = *(nii_input_data + k);
                diff1 = b - a;
                diff2 = TWOPI + b - a;
                diff3 = b - (TWOPI + a);
                if (std::abs(diff2) < std::abs(diff1)) {
                    gra_y = diff2;
                } else if (std::abs(diff3) < std::abs(diff1)) {
                    gra_y = diff3;
                } else {
                    gra_y = diff1;
                }
            }
            if (iz > 0 && iz < end_z) {
                j = sub2ind_4D(ix, iy, iz-1, it, size_x, size_y, size_z);
                k = sub2ind_4D(ix, iy, iz+1, it, size_x, size_y, size_z);
                a = *(nii_input_data + j);
                b = *(nii_input_data + k);
                diff1 = b - a;
                diff2 = TWOPI + b - a;
                diff3 = b - (TWOPI + a);
                if (std::abs(diff2) < std::abs(diff1)) {
                    gra_z = diff2;
                } else if (std::abs(diff3) < std::abs(diff1)) {
                    gra_z = diff3;
                } else {
                    gra_z = diff1;
                }
            }

            if (mode_merge_outputs) {
                *(nii_gra_data + i+nr_voxels*0) = gra_x;
                *(nii_gra_data + i+nr_voxels*1) = gra_y;
                *(nii_gra_data + i+nr_voxels*2) = gra_z;                                
            } else {
                *(nii_gra_x_data + i+nr_voxels*t) = gra_x;
                *(nii_gra_y_data + i+nr_voxels*t) = gra_y;
                *(nii_gra_z_data + i+nr_voxels*t) = gra_z;                                
            }

        // Average differences across neighbourhood types
        // *(nii_gramag_data + i+nr_voxels*t) = (gra_x + gra_y + gra_z) / count_1;
        }
    }
    cout << endl;

    // ========================================================================
    cout << "  Saving output..." << endl;

    if (mode_merge_outputs) {
        save_output_nifti(fout, "gradients", nii_gra, true);
    } else {
        save_output_nifti(fout, "gradient_x", nii_gra_x, true);
        save_output_nifti(fout, "gradient_y", nii_gra_y, true);
        save_output_nifti(fout, "gradient_z", nii_gra_z, true);
    }

    cout << "\n  Finished." << endl;
    return 0;
}
