#include "../dep/laynii_lib.h"
#include <sstream>
#include <set>

int show_help(void) {
    printf(
    "LN2_SKELETONIZE: Compute skeleton.\n"
    "\n"
    "!!!EXPERIMENTAL!!!\n"
    "\n"
    "Usage:\n"
    "    LN2_SKELETONIZE -input input.nii\n"
    "    ../LN2_SKELETONIZE -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Binary nifti image.\n"
    "              will be computed for each volume.\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "\n"
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

    log_welcome("LN2_SKELETONIZE");
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
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_int16(nii1);
    int16_t* nii_input_data = static_cast<int16_t*>(nii_input->data);

    // Prepare output image
    nifti_image* nii_output = copy_nifti_as_int16(nii_input);
    int16_t* nii_output_data = static_cast<int16_t*>(nii_output->data);

    // Set to zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_output_data + i) = 0;
    }

    // Temporary voxels
    nifti_image* nii_temp = copy_nifti_as_int16(nii_output);
    int16_t* nii_temp_data = static_cast<int16_t*>(nii_temp->data);

    // ========================================================================
    // Skeletonize WIP
    // ========================================================================
    cout << "  Computing..." << endl;
 
    std::set<uint32_t> set_initial;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_input_data + i) == 1) {
            set_initial.insert(i);
        }
    }

    int count = 1;
    while (!set_initial.empty()) {
    cout << "  " << count << endl;

    // ========================================================================
    // Step 1: Separate fully connected voxels
    // ========================================================================
    std::set<uint32_t> set_connected, set_candidate;
    uint32_t ix, iy, iz, j, k;
    for (const auto i : set_initial) {
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
        uint8_t count = 0;
        // --------------------------------------------------------------------
        // 1-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && ix < end_x) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            k = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            if (*(nii_input_data + j) == 1) count += 1;
            if (*(nii_input_data + k) == 1) count += 1;
        }
        if (iy > 0 && iy < end_y) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            k = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            if (*(nii_input_data + j) == 1) count += 1;
            if (*(nii_input_data + k) == 1) count += 1;
        }
        // --------------------------------------------------------------------
        // 2-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && ix < end_x) {
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            k = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
            if (*(nii_input_data + j) == 1) count += 1;
            if (*(nii_input_data + k) == 1) count += 1;
        }
        if (iy > 0 && iy < end_y) {
            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            k = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            if (*(nii_input_data + j) == 1) count += 1;
            if (*(nii_input_data + k) == 1) count += 1;
        }

        if (count == 8) {
            set_connected.insert(i);
        } else {
            set_candidate.insert(i);
        }
    }

    // (Optional) Write intermediate output
    for (const auto i : set_connected) {
        *(nii_temp_data + i) = 2;
    }
    for (const auto i : set_candidate) {
        *(nii_temp_data + i) = 1;
    }

    // ========================================================================
    // Step 2: Eliminate candidates that are neighbouring fully connected
    // ========================================================================
    std::set<uint32_t> set_candidate2;
    for (const auto i : set_candidate) {
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
        bool detector = false;
        // ----------------------------------------------------------------
        // 1-jump neighbours
        // ----------------------------------------------------------------
        if (ix > 0 && ix < end_x) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            k = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            if (*(nii_temp_data + j) == 2) detector = true;
            if (*(nii_temp_data + k) == 2) detector = true;
        }
        if (iy > 0 && iy < end_y) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            k = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            if (*(nii_temp_data + j) == 2) detector = true;
            if (*(nii_temp_data + k) == 2) detector = true;
        }
        // ----------------------------------------------------------------
        // 2-jump neighbours
        // ----------------------------------------------------------------
        if (ix > 0 && ix < end_x) {
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            k = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
            if (*(nii_temp_data + j) == 2) detector = true;
            if (*(nii_temp_data + k) == 2) detector = true;
        }
        if (iy > 0 && iy < end_y) {
            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            k = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            if (*(nii_temp_data + j) == 2) detector = true;
            if (*(nii_temp_data + k) == 2) detector = true;
        }

        if (detector == false) {
            set_candidate2.insert(i);
        } 
    }

    // (Optional) Write intermediate output
    for (const auto i : set_candidate2) {
        *(nii_temp_data + i) = 3;
    }

    // ========================================================================
    // Step 3: Replace the initial set with fully connected voxels
    // ========================================================================
    // Only operate on the fully connected voxels in the next iteration
    for (const auto i : set_candidate) {
        *(nii_input_data + i) = 0;
        *(nii_temp_data + i) = 0;
    }

    // Clear temporary data
    for (const auto i : set_connected) {
        *(nii_temp_data + i) = 0;
    }

    // Write out the determined skeleton pieces
    for (const auto i : set_candidate2) {
        *(nii_output_data + i) = count;
    }

    // Refresh the sets in preparation to the next iteration
    set_initial.swap(set_connected);
    set_connected.clear();
    set_candidate.clear();
    set_candidate2.clear();

    count += 1;

    }

    // ========================================================================
    cout << "  Saving output..." << endl;
    save_output_nifti(fout, "test", nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
