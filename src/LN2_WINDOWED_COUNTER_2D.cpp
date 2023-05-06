
#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>
#include <set>
#include <unordered_set>


int show_help(void) {
    printf(
    "LN2_WINDOWED_COUNTER_2D: Count uniquely labeled voxels using circular windows.\n"
    "                         Expects isotropic input.\n"
    "\n"
    "!!! EXPERIMENTAL !!!\n"
    "\n"
    "Usage:\n"
    "    LN2_WINDOWED_COUNTER_2D -input counts.nii\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Expects integer for now.\n"
    "    -radius : (Optional) Maximum distance, in integers.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}


int main(int argc, char*  argv[]) {
    nifti_image *nii_input = NULL;
    char *fin1 = NULL, *fout = NULL;
    int ac;
    bool mode_debug = false;
    int RADIUS = 45;

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
        } else if (!strcmp(argv[ac], "-radius")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -radius\n");
            } else {
                RADIUS = atof(argv[ac]);
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
    nii_input = nifti_image_read(fin1, 1);
    if (!nii_input) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin1);
        return 2;
    }

    log_welcome("LN2_WINDOWED_COUNTER_2D");
    log_nifti_descriptives(nii_input);

    std::ostringstream tag_rad;
    tag_rad << RADIUS;

    // Get dimensions of input
    const uint32_t size_x = nii_input->nx;
    const uint32_t size_y = nii_input->ny;

    const uint32_t end_x = size_x - 1;
    const uint32_t end_y = size_y - 1;

    const uint32_t nr_voxels = size_x * size_y;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii1 = copy_nifti_as_int32(nii_input);
    int32_t* nii1_data = static_cast<int32_t*>(nii1->data);

    // Prepare output nifti
    nifti_image* nii2 = copy_nifti_as_int32(nii1);
    int32_t* nii2_data = static_cast<int32_t*>(nii2->data);
    for (int i = 0; i != nr_voxels; ++i) {
        *(nii2_data + i) = 0;
    }
    nifti_image* nii3 = copy_nifti_as_int32(nii2);
    int32_t* nii3_data = static_cast<int32_t*>(nii3->data);

    // Prepare coordinates
    nifti_image* coord_x = copy_nifti_as_float32(nii_input);
    float* coord_x_data = static_cast<float*>(coord_x->data);
    nifti_image* coord_y = copy_nifti_as_float32(nii_input);
    float* coord_y_data = static_cast<float*>(coord_y->data);

    float ix, iy, iz;
    for (int i = 0; i != nr_voxels; ++i) {
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
        *(coord_x_data + i) = ix;
        *(coord_y_data + i) = iy;
    }

    if (mode_debug) {
        save_output_nifti(fout, "coord-x", coord_x, true);
        save_output_nifti(fout, "coord-y", coord_y, true);
    }

    // ========================================================================
    // Prepare voxels of interest that will be looped through
    // ========================================================================
    // NOTE[Faruk]: Hardcoding 3 step jumps for now.
    float radius_inscribed = (size_x/2) * (size_x/2);
    std::set<uint32_t> idx;
    for (int x = 1; x < size_x; x+=3) {
        for (int y = 1; y < size_y; y+=3) {
            uint32_t i = x * size_y + y;  // sub2ind_2D

            // Check if the voxel is in the inscribed circle.
            // NOTE: Assumes equal data dimensions. e.g. 1000 x 1000
            float xx = *(coord_x_data + i) - (size_x-1)/2;
            float yy = *(coord_y_data + i) - (size_y-1)/2;
            if ((xx*xx + yy*yy) < radius_inscribed) {
                idx.insert(i);
            }
        }
    }

    std::cout << "  Number of voxels : " << nr_voxels << std::endl;
    std::cout << "  Number of samples: " << idx.size() << std::endl;

    for (auto i = idx.begin(); i != idx.end(); ++i) {
        *(nii2_data + *i) = 1;
    }

    if (mode_debug) {
        save_output_nifti(fout, "samples", nii2, true);
    }

    // ========================================================================
    // Evaluate windows
    // ========================================================================
    cout << "\n  Start counting unique voxels within windows..." << endl;

    float RADSQR = RADIUS*RADIUS;

    std::unordered_set<int32_t> unique_ids;
    unique_ids.reserve(idx.size());

    int progressInterval = 100;
    int progressCounter = 0;
    int k = 1;
    float jx, jy;
    for (auto i = idx.begin(); i != idx.end(); ++i) {
        ix = *(coord_x_data + *i);
        iy = *(coord_y_data + *i);

        for (auto j = idx.begin(); j != idx.end(); ++j) {
            jx = *(coord_x_data + *j) - ix;
            jy = *(coord_y_data + *j) - iy;

            if (abs(jx) < RADIUS && abs(jy) < RADIUS) {
                if ((jx*jx + jy*jy) < RADSQR) {
                    unique_ids.insert(*(nii1_data + *j));
                }
            }
        }
        unique_ids.erase(0);
        *(nii2_data + *i) = static_cast<int32_t>(unique_ids.size());
        unique_ids.clear();

        // Output progress in intervals iterations
        progressCounter++;
        if (progressCounter == progressInterval) {
            std::cout << "\r    Processed: " << k << " out of " << idx.size() << std::flush;
            progressCounter = 0;
        }
        ++k;
    }
    std::cout << "\r    Processed: " << k-1 << " out of " << idx.size() << std::endl;

    if (mode_debug) {
        save_output_nifti(fout, "first_samples", nii2, true);
    }

    // ========================================================================
    // Fill unsampled voxels
    // ========================================================================
    cout << "\n  Start filling in unsampled voxels..." << endl;

    uint32_t j;
    for (auto i = idx.begin(); i != idx.end(); ++i) {
        tie(ix, iy, iz) = ind2sub_3D(*i, size_x, size_y);

        // --------------------------------------------------------------------
        // 1-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0) {
            j = sub2ind_3D(ix-1, iy, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }
        if (ix < end_x) {
            j = sub2ind_3D(ix+1, iy, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }
        if (iy > 0) {
            j = sub2ind_3D(ix, iy-1, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }
        if (iy < end_y) {
            j = sub2ind_3D(ix, iy+1, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }

        // --------------------------------------------------------------------
        // 2-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && iy > 0) {
            j = sub2ind_3D(ix-1, iy-1, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }
        if (ix > 0 && iy < end_y) {
            j = sub2ind_3D(ix-1, iy+1, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }
        if (ix < end_x && iy > 0) {
            j = sub2ind_3D(ix+1, iy-1, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }
        if (ix < end_x && iy < end_y) {
            j = sub2ind_3D(ix+1, iy+1, iz, size_x, size_y);
            *(nii2_data + j) = *(nii2_data + *i);
        }
    }

    // Save
    save_output_nifti(fout, "counts_rad-"+tag_rad.str(), nii2, true);

    cout << "\n  Finished." << endl;
    return 0;
}
