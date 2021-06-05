#include "../dep/laynii_lib.h"

int show_help(void) {
    printf(
    "LN2_BORDERIZE: Reduce rim file to its borders.\n"
    "\n"
    "Usage:\n"
    "    LN2_BORDERIZE -input rim.nii\n"
    "    LN2_BORDERIZE -input rim.nii -jumps 3\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Any nifti file with integers. For instance segmentation results,\n"
    "              parcellations, or 'winner maps'.\n"
    "    -jumps  : (Optional) 1, 2 or 3 jump neighbourhood. Default is 1.\n"
    "              1 gives thinnest borders and 3 gives thicknest borders, because:\n"
    "              1 jump means voxels touching all faces will be zeroed.\n"
    "              2 jump means voxels touching all faces and edges will be zeroed.\n"
    "              3 jump means voxels touching all faces, edges, and corners will be zeroed.\n"
    "    -label  : (Optional) An integer. When given, output will only contain\n"
    "              the borders of voxels labeled with this value\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {
    bool use_outpath = false, mask_label = false;
    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac, jumps = 1, label = 0;

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
        } else if (!strcmp(argv[ac], "-jumps")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -jumps\n");
            } else {
                jumps = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-label")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -label\n");
            } else {
                mask_label = true;
                label = atof(argv[ac]);
            }
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
            use_outpath = true;
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

    log_welcome("LN2_BORDERIZE");
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
    nifti_image* nii_rim = copy_nifti_as_int32(nii1);
    int32_t* nii_rim_data = static_cast<int32_t*>(nii_rim->data);

    // Prepare output
    nifti_image* nii_borders = copy_nifti_as_int32(nii_rim);
    int32_t* nii_borders_data = static_cast<int32_t*>(nii_borders->data);

    // Setting zero
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(nii_borders_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // ------------------------------------------------------------------------
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) > 0){
            nr_voi += 1;
        }
    }

    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));
    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_rim_data + i) > 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Borders
    // ========================================================================
    cout << "\n  Finding border voxels..." << endl;

    uint32_t i, j, ix, iy, iz;
    bool switch_border = false;

    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        i = *(voi_id + ii);  // Map subset to full set
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

        if (jumps >= 1){
            // --------------------------------------------------------
            // 1-jump neighbours
            // --------------------------------------------------------
            if (ix > 0) {
                j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x) {
                j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iy > 0) {
                j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iy < end_y) {
                j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iz > 0) {
                j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iz < end_z) {
                j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
        }

        if (jumps >= 2) {
            // --------------------------------------------------------
            // 2-jump neighbours
            // --------------------------------------------------------
            if (ix > 0 && iy > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix > 0 && iy < end_y) {
                j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iy > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iy < end_y) {
                j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iy > 0 && iz > 0) {
                j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iy < end_y && iz > 0) {
                j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iz > 0) {
                j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iz < end_z) {
                j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
        }

        if (jumps >= 3) {
            // --------------------------------------------------------
            // 3-jump neighbours
            // --------------------------------------------------------
            if (ix > 0 && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix > 0 && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix > 0 && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iy > 0 && iz > 0) {
                j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix > 0 && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iy > 0 && iz < end_z) {
                j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iy < end_y && iz > 0) {
                j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
            if (ix < end_x && iy < end_y && iz < end_z) {
                j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                if (*(nii_rim_data + j) != *(nii_rim_data + i)) {
                    switch_border = true;
                }
            }
        }

        // Copy voxels that neighbor a different value than itself
        if (switch_border) {
            *(nii_borders_data + i) = *(nii_rim_data + i);
        }

        // Reset switch
        switch_border = false;
    }

    // ------------------------------------------------------------------------
    if (mask_label) {
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            i = *(voi_id + ii);  // Map subset to full set
            tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

            if (*(nii_borders_data + i) != label) {
                *(nii_borders_data + i) = 0;
            }
        }
    }
    // ------------------------------------------------------------------------

    save_output_nifti(fout, "borders", nii_borders, true, use_outpath);

    cout << "\n  Finished." << endl;
    return 0;
}
