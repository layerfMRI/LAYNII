#include "../dep/laynii_lib.h"
#include <limits>
#include <sstream>

int show_help(void) {
    printf(
    "LN2_PATCH_FLATTEN: FLatten a patch of cortex using 2D flat coordinate\n"
    "                   and cortical a depth measurement.\n"
    "\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "!!! WORK IN PROGRESS... EXPERIMENTAL PROGRAM. !!!\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "\n"
    "Usage:\n"
    "    LN2_PATCH_FLATTEN -input input.nii -coord_uv uv.nii -coord_d depth.nii -out_grid 50\n"
    "\n"
    "Options:\n"
    "    -help       : Show this help.\n"
    "    -input      : TODO.\n"
    "    -domain     : TODO.\n"
    "    -coord_uv   : TODO.\n"
    "    -coord_d    : TODO.\n"
    "    -out_grid_u : TODO.\n"
    "    -out_grid_v : TODO.\n"
    "    -out_grid_d : TODO.\n"
    "    -debug      : (Optional) Save extra intermediate outputs.\n"
    "    -output     : (Optional) Output basename for all outputs.\n"
    "\n"
    "Notes:\n"
    "    - This program is written for 3D images. We might add 2D image support\n"
    "      in the future depending on requests.\n"
    "\n");
    return 0;
}

int main(int argc, char*  argv[]) {

    nifti_image *nii1 = NULL, *nii2 = NULL, *nii3 = NULL, *nii4 = NULL;
    char *fin1 = NULL, *fout = NULL, *fin2=NULL, *fin3=NULL, *fin4=NULL;
    int ac;
    int grid_u = 10, grid_v = 10, grid_d = 3;
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
        } else if (!strcmp(argv[ac], "-coord_uv")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_uv\n");
                return 1;
            }
            fin2 = argv[ac];
        } else if (!strcmp(argv[ac], "-coord_d")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coord_d\n");
                return 1;
            }
            fin3 = argv[ac];
        } else if (!strcmp(argv[ac], "-out_grid_u")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -out_grid_u\n");
                return 1;
            }
            grid_u = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-domain")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -domain\n");
                return 1;
            }
            fin4 = argv[ac];
        } else if (!strcmp(argv[ac], "-out_grid_v")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -out_grid_v\n");
                return 1;
            }
            grid_v = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-out_grid_d")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -out_grid_d\n");
                return 1;
            }
            grid_d = atof(argv[ac]);
        } else if (!strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            fout = argv[ac];
        } else if (!strcmp(argv[ac], "-debug")) {
            mode_debug = true;
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
        fprintf(stderr, "** missing option '-coords_uv'\n");
        return 1;
    }
    if (!fin3) {
        fprintf(stderr, "** missing option '-coords_d'\n");
        return 1;
    }
    if (!fin4) {
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
    nii3 = nifti_image_read(fin3, 1);
    if (!nii3) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin3);
        return 2;
    }
    nii4 = nifti_image_read(fin4, 1);
    if (!nii4) {
        fprintf(stderr, "** failed to read NIfTI from '%s'\n", fin4);
        return 2;
    }

    log_welcome("LN2_PATCH_FLATTEN");
    log_nifti_descriptives(nii1);
    log_nifti_descriptives(nii2);
    log_nifti_descriptives(nii3);
    log_nifti_descriptives(nii4);

    // Get dimensions of input
    const int size_x = nii1->nx;
    const int size_y = nii1->ny;
    const int size_z = nii1->nz;

    const int nr_voxels = size_z * size_y * size_x;
    const int nr_cells = grid_u * grid_v;

    // ========================================================================
    // Fix input datatype issues
    nifti_image* nii_input = copy_nifti_as_float32(nii1);
    float* nii_input_data = static_cast<float*>(nii_input->data);
    nifti_image* coords_uv = copy_nifti_as_float32(nii2);
    float* coords_uv_data = static_cast<float*>(coords_uv->data);
    nifti_image* coords_d = copy_nifti_as_float32(nii3);
    float* coords_d_data = static_cast<float*>(coords_d->data);
    nifti_image* domain = copy_nifti_as_int32(nii4);
    int32_t* domain_data = static_cast<int32_t*>(domain->data);

    // ------------------------------------------------------------------------
    // Prepare outputs
    nifti_image* out_cells = copy_nifti_as_int32(nii_input);
    int32_t* out_cells_data = static_cast<int32_t*>(out_cells->data);

    for (int i = 0; i != nr_voxels; ++i) {
        *(out_cells_data + i) = 0;
    }

    // ------------------------------------------------------------------------
    // Allocating new nifti for flat images
    nifti_image* flat_cells = nifti_copy_nim_info(nii1);
    flat_cells->datatype = NIFTI_TYPE_INT32;
    flat_cells->dim[0] = 4;  // For proper 4D nifti
    flat_cells->dim[1] = grid_u;
    flat_cells->dim[2] = grid_v;
    flat_cells->dim[3] = grid_d;
    flat_cells->dim[4] = 1;
    nifti_update_dims_from_array(flat_cells);
    flat_cells->nvox = nr_cells;
    flat_cells->nbyper = sizeof(int32_t);
    flat_cells->data = calloc(flat_cells->nvox, flat_cells->nbyper);
    flat_cells->scl_slope = 1;
    int32_t* flat_cells_data = static_cast<int32_t*>(flat_cells->data);

    for (int i = 0; i != nr_cells; ++i) {
        *(flat_cells_data + i) = 0;
    }

    // Flat input projection
    nifti_image* flat_input = copy_nifti_as_float32(flat_cells);
    float* flat_input_data = static_cast<float*>(flat_input->data);
    // Flat projection density map
    nifti_image* flat_density = copy_nifti_as_float32(flat_cells);
    float* flat_density_data = static_cast<float*>(flat_density->data);

    // ------------------------------------------------------------------------
    // NOTE(Faruk): This section is written to constrain the big iterative
    // flooding distance loop to the subset of voxels. Required for substantial
    // speed boost.
    // Find the subset voxels that will be used many times
    int nr_voi = 0;  // Voxels of interest
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(domain_data + i) != 0){
            nr_voi += 1;
        }
    }
    // Allocate memory to only the voxel of interest
    int* voi_id;
    voi_id = (int*) malloc(nr_voi*sizeof(int));

    // Fill in indices to be able to remap from subset to full set of voxels
    int ii = 0;
    for (int i = 0; i != nr_voxels; ++i) {
        if (*(domain_data + i) != 0){
            *(voi_id + ii) = i;
            ii += 1;
        }
    }

    // ========================================================================
    // Find coordinate ranges
    float min_u = std::numeric_limits<float>::max();
    float max_u = std::numeric_limits<float>::min();
    float min_v = std::numeric_limits<float>::max();
    float max_v = std::numeric_limits<float>::min();

    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);

        // Check U coordinate min & max
        if (*(coords_uv_data + nr_voxels*0 + i) < min_u) {
            min_u = *(coords_uv_data + nr_voxels*0 + i);
        }
        if (*(coords_uv_data + nr_voxels*0 + i) > max_u) {
            max_u = *(coords_uv_data + nr_voxels*0 + i);
        }

        // Check V coordinate min & max
        if (*(coords_uv_data + nr_voxels*1 + i) < min_v) {
            min_v = *(coords_uv_data + nr_voxels*1 + i);
        }
        if (*(coords_uv_data + nr_voxels*1 + i) > max_v) {
            max_v = *(coords_uv_data + nr_voxels*1 + i);
        }
    }
    cout << "\n  U coordinate min & max: " << min_u << " | " << max_u << endl;
    cout << "\n  V coordinate min & max: " << min_v << " | " << max_v << endl;

    // ========================================================================
    // Visit each voxel to check their coordinate
    for (int ii = 0; ii != nr_voi; ++ii) {
        int i = *(voi_id + ii);

        float u = *(coords_uv_data + nr_voxels*0 + i);
        float v = *(coords_uv_data + nr_voxels*1 + i);

        // Normalize coordinates to 0-1 range
        u = (u - min_u) / (max_u - min_u);
        v = (v - min_v) / (max_v - min_v);
        // Scale with grid size
        u *= static_cast<float>(grid_u - 1);
        v *= static_cast<float>(grid_v - 1);
        // Cast to integer (floor & cast)
        int cell_idx_u = static_cast<int>(u);
        int cell_idx_v = static_cast<int>(v);

        // Write cell index to output
        int j = grid_u * cell_idx_v + cell_idx_u;  // Cell index
        *(out_cells_data + i) = j;

        // Write visited voxel value to flat cell
        *(flat_input_data + j) += *(nii_input_data + i);
        *(flat_density_data + j) += 1;
    }

    save_output_nifti(fout, "UV_cells", out_cells, true);

    // Take the mean of each projected cell value
    for (int i = 0; i != nr_cells; ++i) {
        if (*(flat_density_data + i) > 1) {
            *(flat_input_data + i) /= *(flat_density_data + i);
        }
    }

    save_output_nifti(fout, "flat_values", flat_input, true);
    save_output_nifti(fout, "flat_density", flat_density, true);

    // ========================================================================
    // Visit each voxel to check their coordinate


    cout << "\n  Finished." << endl;
    return 0;
}
