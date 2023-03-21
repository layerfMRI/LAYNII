
#include "../dep/laynii_lib.h"
#include <sstream>
#include <algorithm>
#include <set>

int show_help(void) {
    printf(
    "LN2_NEIGHBORS: Find first order neighbors of each label.\n"
    "\n"
    "Usage:\n"
    "    LN2_NEIGHBORS -input input.nii\n"
    "    ../LN2_NEIGHBORS -input input.nii\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -input        : Integer nifti image.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;

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

    log_welcome("LN2_NEIGHBORS");
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
    nifti_image* nii_input = copy_nifti_as_int32(nii1);
    int32_t* nii_input_data = static_cast<int32_t*>(nii_input->data);

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    uint32_t nr_voi = 0;  // Voxels of interest
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_input_data + i) != 0){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int32_t* voi_id;
    voi_id = (int32_t*) malloc(nr_voi*sizeof(int32_t));

    // Fill in indices to be able to remap from subset to full set of voxels
    uint32_t ii = 0;
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        if (*(nii_input_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Find unique labels
    // ========================================================================
    set<uint32_t> set_labels;

    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        uint32_t i = *(voi_id + ii);  // Map subset to full set
        set_labels.insert(*(nii_input_data + i));
    }
    cout << "  Unique labels: [ ";
    for (uint32_t value : set_labels) {
        cout << value << " ";
    }
    cout << "]" << endl;
    cout << "  Number of unique labels: " << set_labels.size() << endl;

    // ========================================================================
    // Find connected clusters
    // ========================================================================
    cout << "  Start finding neighbors (3-jump neighborhood)..." << endl;
    uint32_t i, j, ix, iy, iz;

    // Loop through all unique labels
    for (int k : set_labels) {
        set<uint32_t> set_neighbors;
        // Loop though voxels
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            // Map subset to full set
            i = *(voi_id + ii);
            if (*(nii_input_data + i) == k) {

                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                // ----------------------------------------------------------------
                // 1-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0) {
                    j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x) {
                    j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iy > 0) {
                    j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iy < end_y) {
                    j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iz > 0) {
                    j = sub2ind_3D(ix, iy, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iz < end_z) {
                    j = sub2ind_3D(ix, iy, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                // ----------------------------------------------------------------
                // 2-jump neighbours
                // ----------------------------------------------------------------
                if (ix > 0 && iy > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix > 0 && iy < end_y) {
                    j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iy > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iy < end_y) {
                    j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix, iy-1, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix, iy-1, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix, iy+1, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix, iy+1, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iz > 0) {
                    j = sub2ind_3D(ix+1, iy, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }

                // --------------------------------------------------------
                // 3-jump neighbours
                // --------------------------------------------------------
                if (ix > 0 && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix-1, iy-1, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix > 0 && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy-1, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix > 0 && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix-1, iy+1, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iy > 0 && iz > 0) {
                    j = sub2ind_3D(ix+1, iy-1, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix > 0 && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix-1, iy+1, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iy > 0 && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy-1, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iy < end_y && iz > 0) {
                    j = sub2ind_3D(ix+1, iy+1, iz-1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
                if (ix < end_x && iy < end_y && iz < end_z) {
                    j = sub2ind_3D(ix+1, iy+1, iz+1, size_x, size_y);
                    set_neighbors.insert(*(nii_input_data + j));
                }
            }
        }

        // Remove unwanted elements
        set_neighbors.erase(0);
        set_neighbors.erase(k);

        // TODO: Make this part export to a text file
        cout << "    Label " << k << " neighbors: ";
        for (uint32_t value : set_neighbors) {
            cout << value << " ";
        }
        cout << endl;
    }

    cout << "\n  Finished." << endl;
    return 0;
}
