
#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>
#include <set>

int show_help(void) {
    printf(
    "LN2_PROPAGATE: EXPERIMENTAL - Propagate voxels until for a number of discrete\n"
    "               steps. Expects isotropic input.\n"
    "\n"
    "Usage:\n"
    "    LN2_PROPAGATE -input counts.nii\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -input        : Expects integer for now.\n"
    "    -max_dist     : (Optional) Maximum distance, in integers.\n"
    "    -debug        : (Optional) Save extra intermediate outputs.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}


int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    bool mode_debug = false;
    int max_dist = std::numeric_limits<int>::max();
    int iter_smooth = 0;

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
        } else if (!strcmp(argv[ac], "-max_dist")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -max_dist\n");
            } else {
                max_dist = atof(argv[ac]);
            }
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

    log_welcome("LN2_PROPAGATE");
    log_nifti_descriptives(nii1);

    if (max_dist < std::numeric_limits<float>::max()) {
        cout << "    Maximum distance is: " << max_dist << endl;
    } else
        cout << "    No maximum distance is selected." << endl;

    // Get dimensions of input
    const uint32_t size_x = nii1->nx;
    const uint32_t size_y = nii1->ny;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;

    const uint32_t nr_voxels = size_x * size_y;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_int32(nii1);
    int32_t* nii_input_data = static_cast<int32_t*>(nii_input->data);

    nifti_image* nii_output = copy_nifti_as_int32(nii1);
    int32_t* nii_output_data = static_cast<int32_t*>(nii_output->data);


    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_output_data + i) = 0;
    }
    // ========================================================================
    // WIP
    // ========================================================================
    cout << "\n  Start discrete propagation..." << endl;

    uint32_t ix, iy, iz, j;
    for (int d = 0; d != max_dist; ++d) {
        for (int i = 0; i != nr_voxels; ++i) {
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

            if (*(nii_input_data + i) == 0) {
                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    if (*(nii_input_data + j) != 0){
                        *(nii_output_data + i) = 1;
                    }
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    if (*(nii_input_data + j) != 0){
                        *(nii_output_data + i) = 1;
                    }
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    if (*(nii_input_data + j) != 0){
                        *(nii_output_data + i) = 1;
                    }
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    if (*(nii_input_data + j) != 0){
                        *(nii_output_data + i) = 1;
                    }
                }
            } else {
                *(nii_output_data + i) = 1;
            }
        }

        for (int i = 0; i != nr_voxels; ++i) {
            *(nii_input_data + i) += *(nii_output_data + i);
            *(nii_output_data + i) = 0;
        }
    }

    // Add number of points into the output tag
    save_output_nifti(fout, "floodFill", nii_input, true);

    cout << "\n  Finished." << endl;
    return 0;
}
