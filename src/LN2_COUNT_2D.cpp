
#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>
#include <set>
#include <unordered_set>


int show_help(void) {
    printf(
    "LN2_COUNT_2D: Count uniquely labeled voxels using circlular windows.\n"
    "              Expects isotropic input.\n"
    "\n"
    "!!! EXPERIMENTAL !!!\n"
    "\n"
    "Usage:\n"
    "    LN2_COUNT_2D -input counts.nii\n"
    "\n"
    "Options:\n"
    "    -help   : Show this help.\n"
    "    -input  : Expects integer for now.\n"
    "    -domain : Expect bool.\n"
    "    -radius : (Optional) Maximum distance, in integers.\n"
    "    -debug  : (Optional) Save extra intermediate outputs.\n"
    "    -output : (Optional) Output basename for all outputs.\n"
    "\n");
    return 0;
}


int main(int argc, char*  argv[]) {
    nifti_image *nii1 = NULL, *nii2 = NULL;
    char *fin1 = NULL, *fin2 = NULL, *fout = NULL;
    int ac;
    bool mode_debug = false;
    int iter_smooth = 0;
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
        } else if (!strcmp(argv[ac], "-domain")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -domain\n");
                return 1;
            }
            fin2 = argv[ac];
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

    log_welcome("LN2_COUNT_2D");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);

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
    nifti_image* nii_domain = copy_nifti_as_int32(nii2);
    int32_t* nii_domain_data = static_cast<int32_t*>(nii_domain->data);

    // Prepare coordinates
    nifti_image* coord_x = copy_nifti_as_float32(nii1);
    float* coord_x_data = static_cast<float*>(coord_x->data);
    nifti_image* coord_y = copy_nifti_as_float32(nii1);
    float* coord_y_data = static_cast<float*>(coord_y->data);

    uint32_t ix, iy, iz, jx, jy, jz;
    for (int i = 0; i != nr_voxels; ++i) {
        tie(ix, iy, iz) = ind2sub_3D(i, size_x, size_y);
        *(coord_x_data + i) = ix;
        *(coord_y_data + i) = iy;
    }

    nifti_image* nii_output = copy_nifti_as_int32(nii1);
    int32_t* nii_output_data = static_cast<int32_t*>(nii_output->data);

    for (int i = 0; i != nr_voxels; ++i) {
        *(nii_output_data + i) = 0;
    }
    // ========================================================================
    // Prepare voxels of interest that will be looped through
    // ========================================================================
    std::set<uint32_t> idx;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(nii_domain_data + i) != 0) {
            idx.insert(i);
        }
    }
    std::cout << "  Number of non-zero voxels within domain: " << idx.size() << std::endl;

    // ========================================================================
    // Evaliuate windows
    // ========================================================================
    cout << "\n  Start counting unique voxels within windows..." << endl;

    float RADSQR = RADIUS*RADIUS;

    std::unordered_set<int32_t> unique_ids;
    unique_ids.reserve(idx.size());

    int k = 1;
    for (auto i = idx.begin(); i != idx.end(); ++i) {
        ix = *(coord_x_data + *i);
        iy = *(coord_y_data + *i);
        std::cout << "\r    Iteration " << k << " of " << idx.size()  << std::flush;

        for (auto j = idx.begin(); j != idx.end(); ++j) {
            jx = *(coord_x_data + *j) - ix;
            jy = *(coord_y_data + *j) - iy;
            float norm = jx*jx + jy*jy;
            if (norm < RADSQR) {
                unique_ids.insert(*(nii_input_data + *j));
            }
        }
        *(nii_output_data + *i) = static_cast<int32_t>(unique_ids.size());
        unique_ids.clear();
        ++k;
    }

    std::cout << std::endl;

    // Add number of points into the output tag
    if (mode_debug) {
        save_output_nifti(fout, "coord-x", coord_x, true);
        save_output_nifti(fout, "coord-y", coord_y, true);
    }
    std::ostringstream tag_rad;
    tag_rad << RADIUS;
    save_output_nifti(fout, "counts_rad-"+tag_rad.str(), nii_output, true);

    cout << "\n  Finished." << endl;
    return 0;
}
