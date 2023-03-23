
#include "../dep/laynii_lib.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <vector>
#include <string>


int show_help(void) {
    printf(
    "LN2_NEIGHBORS: Find first order neighbors of each label, export a text or a nifti file.\n"
    "\n"
    "Usage:\n"
    "    LN2_NEIGHBORS -input input.nii\n"
    "    ../LN2_NEIGHBORS -input input.nii\n"
    "    ../LN2_NEIGHBORS -input input.nii -export_nifti\n"
    "\n"
    "Options:\n"
    "    -help         : Show this help.\n"
    "    -input        : Integer nifti image.\n"
    "    -export_nifti : (Optional) Export the neighbor information as a 4D nifti file.\n"
    "                    This is off by default because it might result in very large files.\n"
    "                    The resulting file contains the initial labels as the first volume.\n"
    "                    Other volumes contain the labels of the neighbors for each voxel.\n"
    "                    Note that different labels can have different number of neighbors.\n"
    "                    Therefore, later volumes can contains more zeros.\n"
    "    -output       : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}


int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    bool export_nifti = false;

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
        } else if (!strcmp(argv[ac], "-export_nifti")) {
            export_nifti = true;
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
    // ========================================================================
    nifti_image* nii_input = copy_nifti_as_int32(nii1);
    int32_t* nii_input_data = static_cast<int32_t*>(nii_input->data);

    // TODO[Faruk]: I cannot think of a way to avoid this right now but I think
    // this nifti can be avoided to decrease RAM load, when needed. E.g. I
    // might hold a small vector that holds the indices of the labels. But I
    // need to query that vector by the label value. Grumble grumble...
    nifti_image* idx_label = copy_nifti_as_int32(nii_input);
    int32_t* idx_label_data = static_cast<int32_t*>(idx_label->data);
    for (uint32_t i = 0; i != nr_voxels; ++i) {
        *(idx_label_data + i) = 0;
    }

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
    set<int> set_labels;

    for (uint32_t ii = 0; ii != nr_voi; ++ii) {
        uint32_t i = *(voi_id + ii);  // Map subset to full set
        set_labels.insert(*(nii_input_data + i));
    }

    cout << "  Unique labels: ";
    for (int value : set_labels) {
        cout << value << " ";
    }
    cout << "\n" << endl;
    cout << "  Number of unique labels: " << set_labels.size() << "\n" << endl;

    // Prepare a vector of vectors to hold the neighbournood information
    // NOTE[Faruk]: This is basically like an excel sheet, rows by columns.
    // But the columns are 'jagged'
    std::vector<std::vector<int>> vec_neighbors;

    // ========================================================================
    // Prepare text output
    // ========================================================================
    // Parse output path
    string path = fout;
    std::string dir, file, basename, sep, csv_path_out;
    auto pos1 = path.find_last_of('/');
    if (pos1 != string::npos) {  // For Unix
        sep = "/";
        dir = path.substr(0, pos1);
        file = path.substr(pos1 + 1);
    } else {  // For Windows
        pos1 = path.find_last_of('\\');
        if (pos1 != string::npos) {
            sep = "\\";
            dir = path.substr(0, pos1);
            file = path.substr(pos1 + 1);
        } else {  // Only the filename
            sep = "";
            dir = "";
            file = path;
        }
    }

    // Parse filename
    auto const pos2 = file.find_first_of('.');
    if (pos2 != string::npos) {
        basename = file.substr(0, pos2);
    }

    // Prepare output path
    csv_path_out = dir + sep + basename + "_neighbors" + ".csv";

    // Prepare file
    std::ofstream output_file(csv_path_out);
    if (!output_file.is_open()) {
        std::cout << "  Unable to open text file!\n";
        return 1;
    }

    // ========================================================================
    // Find first order neighbors
    // ========================================================================
    cout << "  Start finding neighbors (3-jump neighborhood)..." << endl;
    uint32_t i, j, ix, iy, iz, max_nr_neighbors = 0;

    // Loop through all unique labels
    int c = 0;
    for (int k : set_labels) {
        set<uint32_t> set_neighbors;
        // Loop though voxels
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            // Map subset to full set
            i = *(voi_id + ii);
            if (*(nii_input_data + i) == k) {
                *(idx_label_data + i) = c;

                tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);

                // ------------------------------------------------------------
                // 1-jump neighbours
                // ------------------------------------------------------------
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
                // ------------------------------------------------------------
                // 2-jump neighbours
                // ------------------------------------------------------------
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

                // ------------------------------------------------------------
                // 3-jump neighbours
                // ------------------------------------------------------------
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

        cout << "    Label " << k << " neighbors: ";
        for (uint32_t value : set_neighbors) {
            cout << value << " ";
        }
        cout << endl;

        // Insert the set as a vector into the main vector
        vector<int> vec_temp(set_neighbors.begin(), set_neighbors.end());
        vec_neighbors.push_back(vec_temp);

        // Update maximum number of neighbors (useful for preparing 4D output)
        if (max_nr_neighbors < set_neighbors.size()) {
            max_nr_neighbors = set_neighbors.size();
        }

        c += 1;
    }
    cout << endl;
    cout << "  Maximum number of neighbors:" << max_nr_neighbors << endl;

    // ====================================================================
    // Export to a CSV (comma separated value) text file
    // ====================================================================

    // Set first row as column titles
    output_file << "Label" << ",";
    for (int i = 0; i != max_nr_neighbors; ++i) {
        output_file << "Neighbor-" << i+1 << ",";
    }
    output_file << "\n";

    // Insert values in each row
    c = 0;
    for (int value : set_labels) {
        output_file << value << ",";
        for (int m = 0; m != vec_neighbors[c].size(); ++m) {
            output_file << vec_neighbors[c][m] << ",";
        }
        output_file << "\n";
        c += 1;
    }

    output_file.close();

    // ========================================================================
    // Export a 4D nifti output
    // ========================================================================
    if (export_nifti) {
        cout << "  Exporting nifti..." << endl;
        nifti_image* nii_output = nifti_copy_nim_info(nii_input);
        nii_output->dim[0] = 4;  // For proper 4D nifti
        nii_output->dim[1] = size_x;
        nii_output->dim[2] = size_y;
        nii_output->dim[3] = size_z;
        nii_output->dim[4] = max_nr_neighbors + 1;  // +1 for the initial label
        nifti_update_dims_from_array(nii_output);
        nii_output->nvox = nr_voxels * (max_nr_neighbors + 1);
        nii_output->nbyper = sizeof(int32_t);
        nii_output->data = calloc(nii_output->nvox, nii_output->nbyper);
        int32_t* nii_output_data = static_cast<int32_t*>(nii_output->data);

        // --------------------------------------------------------------------
        for (uint32_t ii = 0; ii != nr_voi; ++ii) {
            i = *(voi_id + ii);  // Map subset to full set

            // First volume is the input labels
            *(nii_output_data + i) = *(nii_input_data + i);

            // Populate the neighbors
            j = *(idx_label_data + i);
            for (int m = 0; m != vec_neighbors[j].size(); ++m) {
                *(nii_output_data + nr_voxels*(m+1) + i) = vec_neighbors[j][m];
            }
        }

        save_output_nifti(fout, "neighbors", nii_output, true);
    }

    cout << "\n  Finished." << endl;
    return 0;
}
